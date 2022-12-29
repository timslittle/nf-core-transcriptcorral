# asmrna_identityclass4a.pl
# asmrna_dupfilter3c.pl identityclass() cleaned up of debug crap UPD1912

# use strict;
our(%aasize, %sizeval, %trsize, %aaqual, %oids, %oidsof, %cdsoff); # all from readSizes(); one hash? == sizeval
our(%aablast, %aablastref,$naabl); $naabl=0;
our(%aacluster,%aaclustermain); # globals for readAAcdhit
our(%better, %outrows, %validids, %bspans); ## %validids was %outids; 
our($aconsensus, %aconsensus); # AACONS

BEGIN{ warn "# required asmrna_identityclass4a.pl\n"; }
#======== 20.01.14 moved asmrna_identityclass4a.pl back to  asmrna_dupfilter4a.pl ===========================
# BEGIN{ require "asmrna_identityclass4a.pl"; } # in asmrna_dupfilter4a.pl 

=item identityclass

  from  aabugs4qual/tsaevg/cdsidentclass.sh
  input from above putspans, sorted by ^Qclen, ^Align, vTtlen:
     qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits) 

  update classes per aabugs4qual/tsaevg/classparts.pl
  .. need blastp/blast table input for this..
  .. need 3-4 final categories:  keep-main, keep-alts, trash-fragments (alt,noalt) 
     .. keep includes partial prot where good enough.
     .. trash-alts maybe (not fragment but not distinct enough)
     .. add poor cut to main,alt,noclass per below a,b qualities.
     
  cl_main.alttr and cl_alts, keep if: 
  a. homology unique or best (near best for alt),  or
  b. aasize >= 45% of main and aacomplete? and not aasame-as-main?
  c. main == fragment if aasize < MINAAFULL or (aasize < MINAAPART and partial/utrpoor)

  identityclass parameters ... see above now
    $TINYALN=$ENV{mina}||50; 
    $IS_CDSALIGN=$ENV{cdsw}||0; # default using tralign, tests better than cds
    $PHIALN=$ENV{ahi}||98; # was 65; # use mina?
    $PHI=$ENV{phi}||99; $PMID=$ENV{pmid}||90; $PLOW=$ENV{plow}||80; 
    $ALTFRAG: add isfrag pct; now 50 (0.5)
  
=cut

sub identityclass {
  my($outh, $infile, $insorted)= @_;

  use constant TEST3 => 1; # 13aug09 test fixes to alt/main classing # UPD1912: all accepted now
  use constant TEST1602 => 1; ## UPD 2016.02 problem w/ dup equal mains, equivalence after 1st see both as mains, 1 to be alt,
  use constant TEST1603 => 0; #UPD1912 off, only for now-removed eqgene measures, was 1; ## test use here eqgenes/eqexons??
  
  # our(%validids,$IS_CDSALIGN);
	my $havevalidids= scalar(%validids)?1:0; ## need this?

  ### infile is putspans() alntab: 
  ###  Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits
  ###    1     2     3   4     5     6     7     8    9
  ### sort ^Qclen, ^Align, vTtlen  : orig
  ### sort ^Qclen, ^Align, vQtlen, vTtlen, =Qid, =Tid : now
  ### maybe ^Ident before =Tid ?  after ^Align ?
  ### should it be ^Qclen, =Qid, ^Align, vTtlen :  no want top Aligns first
  ### added vQtlen before vTtlen, else get utrbad before good ***
  # my $ALNSORTORD='-k2,2nr -k7,7nr -k6,6n '; #orig
  # my $ALNSORTORD='-k2,2nr -k7,7nr -k6,6n -k1,1 -k4,4'; #add IDs to order ties

  my $ALNSORTORD='-k2,2nr -k7,7nr -k3,3n -k6,6n -k1,1 -k4,4';  #add vQtlen/3 before vTt/6 IDs to order ties
  # unless($IS_CDSALIGN) { #..not sure what is right, want cds to play role in best choice
  #   $ALNSORTORD='-k2,2nr -k7,7nr -k3,3n -k6,6n -k1,1 -k4,4';  # test same as IS_CDSALIGN  
  # }
  
  my($inh, %class,%bestmatch,%ismain, %mainof, %havepair); 
  if(ref($infile)) { $inh= $infile; } # STDIN,.. sort ??
  else {
    my $tmpdir=$ENV{TMPDIR};
    $tmpdir= './' unless($tmpdir and $tmpdir ne "/tmp"); ## { $tmpdir=`pwd`; chomp($tmpdir); }
    if($insorted) { open(IN,$infile) or die "reading $infile"; $inh=*IN; }
    else { open(IN,"sort -T $tmpdir $ALNSORTORD $infile |") or die "reading $infile"; $inh=*IN; }
  }
  
  use constant UPD1912doswapc => 0; # 0, turn this off, likely bug
  use constant UPD1912fixNOMAIN => 1; # use DEBUGnomain == $ENV{fixnomain} for testing
  use constant ANTI1912c => 1; ## antisense checks
  my $PAL_ANTIMIN= 60; #ANTI1912c opt? needs test w/ refs
  # ref_arath: antisense pAlign for same gene, diff gene == 60 is good
  #  sameg.pal	n:146 av:75 med:65 81 90;  quartiles: .25, .50, .75
  #  diffg.pal	n:333 av:66 med:44 70 91
  
  my $DEBUGnomain= $ENV{fixnomain}||0; # or UPD1912fixNOMAIN
  warn "# UPD1912doswapc => ",UPD1912doswapc," UPD1912fixNOMAIN => ",UPD1912fixNOMAIN," ANTI1912c => ",ANTI1912c," fixnomain=",$DEBUGnomain,"\n";
    
  my($lastd)=("");
  while(<$inh>) {  # input sorted align table
    next if(/^Qid|^\W/); chomp; 
    my @v= split"\t"; 
    my($qd,$qc,$qw,$td,$tc,$tw,$aln,$iden,$bits)= @v; 
    my($isfrag,$aclass,$alnmax,$pid,$pal,$samesize)= (0) x 10;
    
    my $isanti= ($IS_CDSALIGN and $bits < 0)?1:0;
    my $antiflag= ($isanti) ? "/-sense" : "/."; ## append to bestmatch ??
    $isfrag= $aclass="";
    $samesize=($tc == $qc)?1:0; # OR ? abs($tc - $qc) < $NHIALN
    my $swapqt=($tc > $qc)?1:0;
    
    if(AACONS and UPD1912doswapc) {
      if($nconsensus and $samesize) {
        #UPD1912 : Is this swapqt a problem, for dup loci, almost identicals? changes their sort order some what randomly
        # .. Could replace this aacons boost by adding .acons as decimal to alntab Qclen, Tclen so sort will favor tc.tcon > qc.nocon
        my $tdcon= $aconsensus{$td}||0; my $qdcon=$aconsensus{$qd}||0;
        $swapqt=1 if($tdcon > $qdcon);
      }
    }
    if($swapqt){ ($qd,$td)=($td,$qd); ($qc,$qw,$tc,$tw)=($tc,$tw,$qc,$qw); } # swap to greater cds?

    # eqgene classifier : drop for now
		# if(EQGENE_CHANGES_NOALN && $lastd && $lastd ne $qd) {}
				
		$lastd=$qd; # here?
    if($td eq "self" or $qd eq $td) { # putspan fix for no-align cases
      $td=$qd; # $aclass="noclass"; << should be this, but ??
      $bestmatch{$qd}="$td,100/100/self1" unless($bestmatch{$qd});
      next;  # can lose eqgene if no further td/qd ..
    }
    
    my($qsize,$tsize)= ($IS_CDSALIGN) ? ($qc,$tc) : ($qw,$tw);
    # note: alnmax is min(qsize,tsize) not max
    $alnmax= ($qsize>$tsize and $tsize>0)?$tsize:($qsize>0)?$qsize:$tsize;  
    $isfrag= ($tsize < $ALTFRAG*$qsize)?"frag":"";
    
    $pid= ($aln<1)?0: int(0.5+ 100*$iden/$aln); $pid=100 if($pid>100);
    $pal= ($alnmax<1)?0 : int(0.5+ 100*$aln/$alnmax); $pal=100 if($pal>100);
    
    my $tinyalndiff= ((($alnmax - $aln) < $NHIALN)) ? 1 : 0; # TEST3 add, use w/ $pal

    ##UPD1912 drop eqgenes parts for now
		#u $havepair{$qd}{$td}= $pal; #UPD1912:not used ?? part of EQGENE; change this to pal align val for altpar TEST1603?
		#u if($neqgene>0) {}
		
    if(ANTI1912c and $isanti) { # UPD1912 maybe .. TEST .. NOT HERE, need bestmatch for output trclass
      $antiflag= "/." if($pal < $PAL_ANTIMIN); # turn off flag
      #   next if($pal < 60); # force separate loci for rev-aa, partial cds align??
      #    # treat as NO ALIGN ? skip before bestmatch NO cant skip bestmatch
    } 

    my $pidalnval="$pid/$pal$antiflag"; # old: $pid/$pal
    my $tisaltpar=0; # treat like $skiptinyaln for now
    if(TEST1603) {
      $tisaltpar=1 if($pidalnval=~/altpar\d/); #?  and $qclass // not altparx\d
    }
    if(TEST1603) { # UPD1912:off, this may be sort of working right now..
    $bestmatch{$qd}="$td,$pidalnval" unless($bestmatch{$qd} and not($bestmatch{$qd}=~/altpar\d/));
    $bestmatch{$td}="$qd,$pidalnval" unless($bestmatch{$td} and not($bestmatch{$td}=~/altpar\d/)); 
    } else {    
    $bestmatch{$qd}="$td,$pidalnval" unless($bestmatch{$qd});
    $bestmatch{$td}="$qd,$pidalnval" unless($bestmatch{$td});  
    }
        
    my $qclass= $class{$qd}||""; my $tclassD= $class{$td}||"";
    #D warn "#dnomain1 $qd,$qc $td,$tc pid=$pid, pal=$pal, al=$aln, qcla=$qclass, tcla=$tclassD\n" if($DEBUGnomain);

    if($samesize and $qclass =~ /alt/ and $ismain{$td}) { # check/stop alt1 <> alt2 assigns
      next if($bestmatch{$qd} =~ /$td,/); # problem here? for many equal bestmatch
    }

    #..........    
    ## reclassify? althi > althi100 for 2ndary alts; althi100 = unneeded duplicates : TEST
    
    my $skiptinyaln=0;
    if(UseTINYALNBASES) {
      if($alnmax >= $MINCDS) {  
        my $minbaseover= ($alnmax >= $TINYALNCDSLEN)? $TINYALNBASES : $alnmax * 0.10; ##$MIN_PCTOVERLAP;
        $skiptinyaln=($aln < $minbaseover)?1:0;
      } else {
        $isfrag= "frag" unless($isfrag); #or "tiny" ? # dont skip assign alt/frag to tiny cds, but maybe always set isfrag ?
      }
    } else {  
      $skiptinyaln= ($pal < $TINYALN)?1:0;
    }
  
    # set alt aclass here, or not.
    if( $skiptinyaln or $tisaltpar) { } ## defer: $aclass="noalign$isfrag";  
    elsif( $tc == 0 and $qc > 0 and $pid >= $PLOW ) { $aclass="frag0aa"; } # tc==0 special case "frag0"
    elsif( $pid >= $PHI ) { 
      $aclass= ($isfrag)?"parthi":"althi"; # Primary classifier; change parthi to althipart ? althi$isfrag ?
      $aclass .= "1" if($pal >= $PHIALN or $tinyalndiff); # althi1 == hi-align + hi-ident, nearly same (identicals gone here)
      }
    elsif( $pid >= $PMID ) { $aclass="altmid$isfrag"; }
    elsif( $pid >= $PLOW ) { $aclass="altlo$isfrag"; }
    # else { $aclass=""; } # otherwise this pair are not aligned well enough to call same locus

    if(ANTI1912c and $isanti) { # UPD1912 maybe .. TEST .. NOT WORKING ??
      $aclass="" if($pal < $PAL_ANTIMIN); # force separate loci for rev-aa, partial cds align??
      $pidalnval =~ s/$antiflag/./ unless($aclass);
    } 
    
=item samesize/samebestmatch 
  problem of 2+mains w/ several essentially identicals, 
  after 2 mains created, then row linking them needs to break tie
  
  * problems losing main links with large num 910 $samesize pairs .. sort order isn't enough
  
=cut

    if(my $ocl= $class{$td}) {  
      if($ocl =~ m/^alt|^part/) { 
        $class{$td}= $aclass if($aclass eq "althi1"); 
        next;   # TEST1602  ok 
        # UPD1912, check this, NOMAIN major group now, bad
        # maybe want to drop into FINDMAIN below, add qd main if missing, but dont change old td class to aclass?
        #no if($DEBUGnomain == 2) {  ## UPD1912fixNOMAIN .. TEST 
        #no  $aclass=$ocl unless($aclass eq "althi1"); # and continue into FINDMAIN
        #no } else { next;  }  # 
      } elsif($ocl =~ m/^main/) { 
        if($qclass =~ m/^main/ or $tisaltpar) { 
          # skip to below aclass set, test5, yes, do this also?
          } 
        else { next; } # * YES, need this : test4
      } 
    }
      
    if($aclass) { 
      my $qmain=$qd; my $attr= $pidalnval; ## "$pid/$pal$antiflag";
      ## use constant FINDMAINID => 1;  
      if($class{$qd}) {   
        my $qnext= $qd; 
        my $more= ($bestmatch{$qnext})?1:0;
        
        my %qdid=( $qd => 1, $td => 1);
        while($more) { 
          $more=0;
          my($qbest,$qpal)=split",",$bestmatch{$qnext};
          
          if($qbest) {
          $qbest=0 if($ismain{$qnext} and not $ismain{$qbest}); # TEST1603 want always? maybe bad?    
                
          #n if($DEBUGnomain) { 
          #n# this way adds more /qmain tags, doesnt resolve alt1<=>alt2 lost mains, for samesize near identicals
          #n$qbest=0 if($class{$qnext} eq "main" and  $class{$qbest} ne "main"); # TEST 1912
          #n} else {
          #n $qbest=0 if($ismain{$qnext} and not $ismain{$qbest}); # TEST1603 want always? maybe bad?          
          #n}
          }
          
          if($qbest and not $qdid{$qbest}) {
            $qnext= $qbest; $qdid{$qbest}=1; 
            $more=1 if($class{$qnext} and $bestmatch{$qnext});
          }
        }
        if($qnext ne $qd) { $qmain= $qnext; 
          unless(UPD1912fixNOMAIN) { $attr="$pidalnval/$qmain"; } # $DEBUGnomain
          } 
      }
      
      $ismain{$qmain}++;   
      $class{$td}=$aclass; #? problem here when class{td} == main set before, should not reset to alt
      
      # TEST1602 this way works.  if($bestmatch{$td}) must be true from above
      if(not $bestmatch{$td} or $bestmatch{$td} =~ /^$qd/){ $bestmatch{$td}="$qd,$attr"; }
      #old if($bestmatch{$td} =~ /^$qd/) { $bestmatch{$td}="$qd,$attr"; }
      #old elsif(not $bestmatch{$td}) { $bestmatch{$td}="$qd,$attr"; }

      ## add to prevent 2+ circular alta <> altb <> altc with no main ..
      $class{$qmain}="main" unless($class{$qmain}); # should always be:  $qc >= $tc//$samesize..
      if(UPD1912fixNOMAIN) { $mainof{$td}= $qmain unless($mainof{$td}); }
      #D warn "#dnomain2 qmain:$qmain $td,$tc  qmcla=$class{$qmain}, tcla=$class{$td}\n" if($DEBUGnomain);
    }
      
  } close($inh);

  # UPD1912fixNOMAIN/DEBUGnomain .. findmain() may work, does on simple test case w/ many samesize dups
  # should it be used for other than %class keys? probably not, mainof only set for %class'ed
  sub findmain { 
    my($td,$mainof)=@_;
    my $nextd=$td; my(%mdid);
    while (my $md= $mainof->{$nextd}) { last if($mdid{$md}++); $nextd=$md; }
    return ($nextd eq $td) ? "" : $nextd;
  }
  
  # END:  print trclass table; add more fields to output: aaqual, aablast, tr,aa sizes?
  { my($q,$pal,$c,$d);
  foreach $d (sort keys %class) {
    ($q,$pal)=split",",$bestmatch{$d}; $c= $class{$d}; # alts mostly here in class list
    if(UPD1912fixNOMAIN) { my($md)= findmain($d, \%mainof); $pal .= "/$md" if($md and $md ne $q); }
    my @cla= classifytr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  foreach $d (sort grep{ not($class{$_}) } keys %ismain) { # mains, not in class list
    ($q,$pal)=split",",$bestmatch{$d}; $c= "main";    
    #d if(UPD1912fixNOMAIN) { my($md)= findmain($d, \%mainof); $pal .= "/$md" if($md and $md ne $q); }
    my @cla= classifytr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  foreach $d (sort grep{ not($class{$_} or $ismain{$_}) } keys %bestmatch) { # unclassified by shared aligns
    ($q,$pal)=split",",$bestmatch{$d}; $c= "noclass";  
    #d if(UPD1912fixNOMAIN) { my($md)= findmain($d, \%mainof); $pal .= "/$md" if($md and $md ne $q); }
    my @cla= classifytr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  }  
  
}


=item classifytr

  classifier of locus primary, alternate, fragment, redundant 
  using identityclass() collection of overlapping transcripts (CDS or full tr)  
  focused on CDS qualities, maybe should be classifyCDS() ..
  
=cut


sub classifytr {
  my($tid,$cla,$qid,$pidal)= @_;
  $qid||="0"; $pidal||="0";

  # use constant { kAATINY => 1, kAAUTRBAD => 2, kAAUTRPOOR => 4, kAADUP => 8, kAAGAPS => 16, kNONCODE => 32, };
  # use constant NOTPOORBAD => kAATINY + kAADUP + kAAGAPS; # - kAAUTRPOOR - kAAUTRBAD - kAAGAPS
  # use constant NOTPOOR => kAATINY + kAADUP + kAAUTRBAD + kAAGAPS; # - kAAUTRPOOR
  
  my($pidi,$pali)= $pidal =~ m/^(\d+).(\d+)/;
  my($isanti)= $pidal =~ m/\-sense/ ? 1 : 0;
  my $aw= $aasize{$tid} || 0;
  my $tqual= $aaqual{$tid} || ""; #? parse for "aasize,pcds,aaqual"
  my $tw= $trsize{$tid} || 1; 
  my $pcds= int(300*$aw/$tw);  
  my $butr= ($tw - 3*$aw); # okutr if butr<=300p, modify BAD_CDSUTR for small aw
  my $ispoor= ($aw < $AAMIN or ($tqual =~ m/partial/ and $aw < $AAPART))?kAATINY:0;
  my $tbits= $aablast{$tid} || "0,0";
  my($tbscore,$tbref)=split",",$tbits;
  my $aacons= (AACONS) ? ($aconsensus{$tid}||0) : 0;
  my $mapqloc= ($nineqgene) ? $gmapqual->{$tid} : 0;  
  my($mapqual,$maploc)= ($mapqloc)?(split"\t",$mapqloc):(0,0);
  # new mapqual: 85a,100i,2469l,4x,Spl:29%,KB668900 == %align,%ident,bases,exons,Spl:split

  my($tbrscore,$tbrefbest)= ($tbref) ? (split",", $aablastref{$tbref}) : (0,0); 
  my $rbits2= $aablastref{$tid} || "0,0";
  my($tbrscore2,$tbref2)= split",", $rbits2; # this is revhash tid => score,ref refbest

  my $aaclus= $aacluster{$tid} || "0,0";  
  my($aamainid,$aamainpi)= split",", $aaclus; # ="$mainid,$pi";
  if($aamainid eq $tid) {  $aamainid=""; $aamainpi=0; }
  ## FIXME: aamain can be bad/missing; fix before this, need cds/tr align info to know if aamainid is good.
  
  if($tqual =~ m/gapbad/) { $ispoor |= kAAGAPS; }  # gaps discrepancy w/ aaqual utrbad and pcds utrbad 
  if($butr >= $MINUTR) {  
    if($tqual =~ m/utrbad/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRBAD; } 
    elsif($tqual =~ m/utrpoor/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRPOOR; }
    else { $ispoor |= ($pcds < $BAD_CDSUTR)? kAAUTRBAD: ($pcds < $OK_CDSUTR)? kAAUTRPOOR : 0; }
  }

  my $tcode=0;
  if(CODEPOT1607) {
    $tcode= $sizeval{$tid}{'codepot'}||"";  
    if($tcode and $tcode =~ /^Noncode/) {  
    $tcode=0 if($tbscore>0);
    $tcode=0 if($aw > 139 or ($aw > 99 and not $ispoor)); #?? which? many utrorf/utrbad for aa>99
    $ispoor |= kNONCODE if($tcode and $tcode =~ /^Noncode/);  
    }
  }

  if($aamainpi > $AADUP_IDENT) { # TEST aacluster added info on dup prots; AADUP_IDENT=98 default
    $cla.="a2" unless($cla=~/hi[012]/);  
    $ispoor |= kAADUP if($cla =~ /althi/);  
  }
  unless( $ispoor & kAADUP ) {
    if($pidal =~ m/altmap\d+xeq/ and $cla =~ /althi1|part/) {
      $ispoor |= kAADUP;   
    }
  }
  my $aadupflag= ($ispoor & kAADUP and $aamainid) ? ",aadup:$aamainid" : "";  

  # CHECKME: adding aablast kept 40k more okay, all althi/ahia2 + 2k parthi
  # .. are these true uniq aablast or just althi aadups ?
  my $MINBLASTSCORE= 60; # aablast only? bitscore always?
  my $keepdrop="";
  unless( $tbscore == 0 or $tbits=~/^0,0/) {
    my $risbest= (($tbrefbest eq $tid) or ($tbscore >= $tbrscore))?1:0;# was ($tbrefbest eq $tid), allow 2+ same score
    my $risgood= ($tbrscore > $MINBLASTSCORE and $tbscore >= 0.90 * $tbrscore)?1:0;
    my $risok= ($tbref2 and $tbrscore2 >= $MINBLASTSCORE)?1:0;
    if($risok and $risgood) { $risgood=0 if($tbrscore2 > $tbrscore); } #?
    $ispoor = $ispoor & NOTTINY if($risgood); ## FIXME1712: clear kAATINY for $risok/risgood
    my $isaadup=($ispoor & kAADUP);
    if($risbest) { 
      $tbits.=",refbest"; $keepdrop.="okay:refbest,"; $ispoor=0 unless($isaadup); }
    elsif($risgood) { # was 3rd
      $tbits.=",refgood"; unless($isaadup){ $keepdrop.="okay:refgood,"; $ispoor=0;} } # maybe2 ok ?
    elsif($risok) {  # was 2nd ? why superceed refgood? 2nd ref match?
      $tbits.=",refok"; unless($isaadup){ $keepdrop.="okay:refok,"; $ispoor=0;} }
  }
  
  if($ispoor > kAATINY and $cla !~ /althi|parthi/) {  
    if(UPD1908) { } # maybe rescue althi also here
    if($aw >= $AAMINBAD) { $ispoor = $ispoor & NOTPOORBAD; } ##  ^ (kAAGAPS+kAAUTRPOOR+kAAUTRBAD);  ##?? ^ XOR bad
    elsif($aw >= $AAMINPOO) { $ispoor = $ispoor & NOTPOOR; } ##  ^ (0+kAAUTRPOOR);
  }

  # cdsantisense, try to guess true cases from low align, not NONCODE, aacons, other ?
  if(UPD1912 and $isanti and ($ispoor or $cla =~ /parthi/)) {
    my $antiok= 0;
    $antiok += 1 if($pali < 75);  $antiok += 2 if($pali < 55);
    $antiok += $aacons if($aacons > 1);
    if($antiok) { 
      $ispoor = $ispoor & (kAATINY + kAADUP + kAAGAPS + kNONCODE);
      if($cla =~ /parthi/ and $antiok > 2) { $cla="altmidfrag"; } # change to not part? "altmid$isfrag";  
    }
  }
  
  if(AACONS) {
    if($aacons > 0) {
      if(UPD1908 and $cla =~ /parthi|frag0aa/) { } # maybe also skip aacons rescue of poor-noclass
      else { $keepdrop .="okay:aacons$aacons,"; }
      $aadupflag .= ",aacons:$aacons";  #? remove aadup:mainid or not
    }
  }

  # FIXME1805: rescue short/part cds with good/best match to refgenes
  if($cla =~ /parthi|frag0aa/) {
    $keepdrop.= "drop";
    
  } elsif($cla =~ /althi/) {
    if(UPD1908) {  
        if($ispoor == kNONCODE) { $cla =~ s/a2//; $cla =~ s/hi1/hi/; $cla= $cla."nc"; $keepdrop.="okay"; } ## dont rescue NC alt if other poor qual
        $keepdrop.= (($ispoor & NOTPOORBAD) or ($ispoor and $pali >= 95))?"drop":"okay";
    } else {
        $keepdrop.= ($ispoor)?"drop":"okay";  # ispoor vs main size?
    }

  } elsif($cla =~ /main/) {
    if(UPD1908) { if($ispoor & kNONCODE) { $cla =~ s/a2//; $cla= $cla."nc"; $keepdrop.="okay"; } } # rescue NC main if other ispoor
    $keepdrop.= ($ispoor)?"drop":"okay";   
    
  } elsif($cla =~ /noclass/) {
    if(UPD1908) { if($ispoor == kNONCODE) { $cla= $cla."nc"; $keepdrop.="okay"; } } # dont rescue NC noclass if other ispoor
    $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop? ditto to altmid
    
  } else { # other altmid/low 
    $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop this class, ispoor maybe uniq prot: maybeok?
  }
  
  my $okay;
  if($keepdrop =~ /drop/) {
    if($keepdrop =~ /okay/) { $okay= "maybeok"; } else { $okay= "drop"; }
  } else {
    $okay= "okay";
  }
  
  $tbits.= $aadupflag if($aadupflag); # defer, AFTER refbest/good/..
  $tbits.= ",chrmap:$mapqual,$maploc" if($maploc and $mapqual);
  $tbits.= ",pflag:".int($ispoor); # DEBUG info, but keep
  $tbits="aaref:$tbits" unless($tbits=~/^0,0/); #FIXME: add tag to tbits output
  
  # if(defined $eqflag{$tid}) { # TEST1603
  #   my $eqfl= $eqflag{$tid}{$qid}||""; 
  #   if($eqfl) { $eqfl="$qid/$eqfl,"; }
  #   my @q= grep{ $_ ne $qid } sort keys %{$eqflag{$tid}}; 
  #   $eqfl .= join",",map{ "$_/".$eqflag{$tid}{$_} }@q;  
  #   $tbits.= ",feq:$eqfl";
  # }
  
  return (wantarray) ? ($tid,$okay,$cla,$qid,$pidal,$tqual,$tbits) : $okay;
}


1;
