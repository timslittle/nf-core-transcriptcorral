# asmrna_overlocusclass4a.pl

# use strict;
our(%aasize, %sizeval, %trsize, %aaqual, %oids, %oidsof, %cdsoff); # all from readSizes(); one hash? == sizeval
our(%aablast, %aablastref,$naabl); $naabl=0;
our(%aacluster,%aaclustermain); # globals for readAAcdhit
our(%better, %outrows, %validids, %bspans); ## %validids was %outids; 
our($aconsensus, %aconsensus); # AACONS

=item sub classifyFullTr
  
  further classing, after -CDSALIGN with identityclass()/classifytr(),
  using -noCDS (ie. full self-mrna-blastn table of aligns, CDS-main/alt classified)
  
  add classifyTrFull/TrFrag/ that adds tr overlaps outside(?) of CDS x CDS:
  separate CDS v UTR overlaps, extra locus classes depend on how much of which type
    .. need CDS qualities: hoscore (aablast), aasize, %CDS (coding pot)
    .. trfrag overlaps CDS of trgood > drop/relcass as partof or altof trgood.
      .. unless trfrag.UTR only over trgood.CDS and trfrag.CDS has good code
    .. trfrag overlaps only UTR of trgood, depends on trfrag coding quals: other locus or partof trgood
      .. trfrag.CDS overlaps trgood.UTR, trfrag > partof:trgood unless has strong code quals
      .. trfrag.UTR only over trgood.UTR, trfrag > other locus (code/noncode?)
    .. trfrag overlaps trfrag, depends on CDS v UTR, code quals

=item sub overlocusclass

  temp in asmrna_dupfilter3_overlocsub.pm
  now in asmrna_overlocusclass4a.pl
  
  inputs (via what? overLocusclass instead of identityclass? ) : 
    1. cds.qual with sizes cdslen,aaqual,trlen,cds-offs, and maybe cds-codepot (Code/Noncode)
    2. self-blastn btall table for all mrna x mrna, excluding already classified main x alt set,
        i.e. only diff-locus overlaps
    3. aablast btall table of ref-homol scores [option?]
  
=cut

#not this way# 
# require "asmrna_dupfilter3_overlocsub.pm";
#.........
# asmrna_dupfilter3_overlocsub.pm
# new subs for  evigene/scripts/rnaseq/asmrna_dupfilter3.pl

# use strict;
# use warnings;
#?? package main;

=item classifyFullTr
  from asmrna_dupfilter3.pl:sub classifytr()

  further classing, after -CDSALIGN with identityclass()/classifytr(),
  using -noCDS (ie. full self-mrna-blastn table of aligns, CDS-main/alt classified)
  inputs from overLocusclass
    
  add classifyTrFull/TrFrag/ that adds tr overlaps outside(?) of CDS x CDS:
  separate CDS v UTR overlaps, extra locus classes depend on how much of which type
    .. need CDS qualities: hoscore (aablast), aasize, %CDS (coding pot)
      
      readSizes() collects from aa/cds/mrna.qual: 
        aasize{id}, aaqual{}, trsize{}, oids{}/oidsof{}, cdsoff{} also now
      
    .. trfrag overlaps CDS of trgood > drop/relcass as partof or altof trgood.
      .. unless trfrag.UTR only over trgood.CDS and trfrag.CDS has good code
    .. trfrag overlaps only UTR of trgood, depends on trfrag coding quals: other locus or partof trgood
      .. trfrag.CDS overlaps trgood.UTR, trfrag > partof:trgood unless has strong code quals
      .. trfrag.UTR only over trgood.UTR, trfrag > other locus (code/noncode?)
    .. trfrag overlaps trfrag, depends on CDS v UTR, code quals

=item overlocusclass

  from asmrna_dupfilter3.pl:sub identityclass()
  
  inputs (via what? overLocusclass instead of identityclass? ) : 
    1. cds.qual with sizes cdslen,aaqual,trlen,cds-offs, and maybe cds-codepot (Code/Noncode)
    2. self-blastn btall table for all mrna x mrna, excluding already classified main x alt set,
        i.e. only diff-locus overlaps
    3. aablast btall table of ref-homol scores [option?]

=cut

sub trevidence {
  my($tid,$hiid)= @_; # ,$cla,$qid,$pal
  $hiid||=0; ## $hiid= ( $pid >= $PHI ); ## eg >= 99% ident
  
  my $aw= $aasize{$tid} || 0;
  my $tqual= $aaqual{$tid} || ""; #? parse for "aasize,pcds,aaqual"
  ## add if there:
  my $tbits= $aablast{$tid} || "0,0";
  my($tbscore,$tbref)=split",",$tbits;
  
  my $mapqloc= ($nineqgene) ? $gmapqual->{$tid} : 0;  
  my($mapqual,$maploc)= ($mapqloc)?(split"\t",$mapqloc):(0,0);

  my($tbrscore,$tbrefbest)= ($tbref) ? (split",", $aablastref{$tbref}) : (0,0); 
  my $rbits2= $aablastref{$tid} || "0,0";
  my($tbrscore2,$tbref2)= split",", $rbits2; # this is revhash tid => score,ref refbest

  my $ispoor= ($aw < $AAMIN or ($tqual =~ m/partial/ and $aw < $AAPART))?kAATINY:0;
  my $tw= $trsize{$tid} || 1; 
  my $pcds= int(300*$aw/$tw);  my $butr= ($tw - 3*$aw); # okutr if butr<=300p, modify BAD_CDSUTR for small aw
  
  if($tqual =~ m/gapbad/) { $ispoor |= kAAGAPS; }  # gaps discrepancy w/ aaqual utrbad and pcds utrbad 
  if($butr >= $MINUTR) { # ~ 300b "fixed" average utr sizes
    # should adjust %cds bad for cds size, for large cds> 10k, pcds >=90%, for small cds < 300, pcds may be 33%
    if($tqual =~ m/utrbad/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRBAD; } 
    elsif($tqual =~ m/utrpoor/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRPOOR; }
    else { $ispoor |= ($pcds < $BAD_CDSUTR)? kAAUTRBAD: ($pcds < $OK_CDSUTR)? kAAUTRPOOR:0; }
    # if($butr <= 300) {} elsif($pcds < $BAD_CDSUTR) { $ispoor |= kAAUTRBAD; } elsif($pcds < $OK_CDSUTR) { $ispoor |= kAAUTRPOOR; } 
  }

  my $tcode= 0;  
  if(CODEPOT1607) {
  $tcode= $sizeval{$tid}{'codepot'}||""; # $tcode =~ /^Noncode/ set what? kAATINY+kAAUTRBAD
  $tcode=0 if($tbscore>0);
  # $tcode=0 if($aw > $AAPART);# UPD1908 TESTING,not very accurate but use to separate too many tiny aa
  $tcode=0 if($aw > 139 or ($aw > 99 and not $ispoor)); #?? which? many utrorf/utrbad for aa>99
  $ispoor |= kNONCODE if($tcode =~ /^Noncode/); #?
  }
  
  my $MINBLASTSCORE= 60; # aablast only? bitscore always?
  my $keepdrop=""; # not used here..
  unless( $tbscore == 0 or $tbits=~/^0,0/) {
    my $risbest= (($tbrefbest eq $tid) or ($tbscore >= $tbrscore))?1:0;# was ($tbrefbest eq $tid), allow 2+ same score
    my $risgood= ($tbrscore > $MINBLASTSCORE and $tbscore >= 0.90 * $tbrscore)?1:0;
    my $risok=   ($tbrscore > $MINBLASTSCORE and $tbscore >= 0.50 * $tbrscore)?1:0;
       $risok |= ($tbref2 and $tbrscore2 >= $MINBLASTSCORE)?1:0;
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
  
  if($ispoor > kAATINY and not $hiid) { ## $cla !~ /althi|parthi/ == $hiid
    if($aw >= $AAMINBAD) { $ispoor = $ispoor & NOTPOORBAD; } ##  ^ (kAAGAPS+kAAUTRPOOR+kAAUTRBAD);  ##?? ^ XOR bad
    elsif($aw >= $AAMINPOO) { $ispoor = $ispoor & NOTPOOR; } ##  ^ (0+kAAUTRPOOR);
  }

  $tbits.= ",chrmap:$mapqual,$maploc" if($maploc and $mapqual); #? skip her?
  $tbits.= ",pflag:$ispoor"; # DEBUG info, but keep
  $tbits="aaref:$tbits" unless($tbits=~/^0,0/);

  return($ispoor,$tqual,$tbits); # .. others
}


sub classifyFullTr {
  my($tid,$cla,$qid,$pidal)= @_;
  my($aw,$ispoor,$tqual,$tbits,$mapqual,$maploc)= (0) x 9;
  my $keepdrop="";

  #   use constant { kAATINY => 1, kAAUTRBAD => 2, kAAUTRPOOR => 4, kAADUP => 8, kAAGAPS => 16, kNONCODE => 32, };
  #   use constant NOTPOORBAD => kAATINY + kAADUP; # - kAAUTRPOOR - kAAUTRBAD - kAAGAPS
  #   use constant NOTPOOR => kAATINY + kAADUP + kAAUTRBAD + kAAGAPS; # - kAAUTRPOOR
  $qid||="0"; $pidal||="0";
  my ($aclass,$fclass) =  split"/",$cla; ## ??
  $cla= $aclass; #?

if(1) {
  ($ispoor,$tqual,$tbits)= trevidence($tid,$cla =~ m/althi|parthi/);

  if($cla =~ /cull/) { $ispoor |= kAADUP; } # # add new flag kDUPEXONS ?
  if($cla =~ /althi1|part/ and $pidal =~ m/altmap\d+xeq/) { $ispoor |= kAADUP; } # kDUPEXONS?

  my $aaclus= $aacluster{$tid};  
  if($aaclus) {
    my($aamainid,$aamainpi)= split",", $aaclus; # ="$mainid,$pi";
    if($aamainid eq $tid) {  $aamainid=""; $aamainpi=0; }
    if($aamainpi > $AADUP_IDENT) { # TEST aacluster added info on dup prots; AADUP_IDENT=98 default
      $cla.="a2" unless($cla=~/hi[012]/); # hi1/a2 redundant info
      $tbits.=",aadup:$aamainid";
      $ispoor |= kAADUP if($cla =~ /althi/); # only for /althi|parthi/
    }
  } else { $aaclus ||= "0,0"; }
  
  if($tbits =~ m/aaref:/) {
    my $isaadup=($ispoor & kAADUP);
    if($tbits=~/,refbest/) { $keepdrop.="okay:refbest,"; $ispoor=0 unless($isaadup); }
    elsif($tbits=~/,refgood/ and not $isaadup) { $keepdrop.="okay:refgood,";  $ispoor=0;}
    elsif($tbits=~/,refok/ and not $isaadup) { $keepdrop.="okay:refok,";  $ispoor=0; }
  }
  
} else {  
  
  $aw= $aasize{$tid} || 0;
  $tqual= $aaqual{$tid} || ""; #? parse for "aasize,pcds,aaqual"
  $tbits= $aablast{$tid} || "0,0";
  my($tbscore,$tbref)=split",",$tbits;
  
  my $mapqloc= ($nineqgene) ? $gmapqual->{$tid} : 0;  
  ($mapqual,$maploc)= ($mapqloc)?(split"\t",$mapqloc):(0,0);

  my($tbrscore,$tbrefbest)= ($tbref) ? (split",", $aablastref{$tbref}) : (0,0); 
  my $rbits2= $aablastref{$tid} || "0,0";
  my($tbrscore2,$tbref2)= split",", $rbits2; # this is revhash tid => score,ref refbest

  my $aaclus= $aacluster{$tid} || "0,0";  
  my($aamainid,$aamainpi)= split",", $aaclus; # ="$mainid,$pi";
  if($aamainid eq $tid) {  $aamainid=""; $aamainpi=0; }
  
  my $ispoor= ($aw < $AAMIN or ($tqual =~ m/partial/ and $aw < $AAPART))?kAATINY:0;
  my $tw= $trsize{$tid} || 1; 
  my $pcds= int(300*$aw/$tw);  my $butr= ($tw - 3*$aw); # okutr if butr<=300p, modify BAD_CDSUTR for small aw

  if($tqual =~ m/gapbad/) { $ispoor |= kAAGAPS; }  # gaps discrepancy w/ aaqual utrbad and pcds utrbad 
  if($butr >= $MINUTR) { # ~ 300b "fixed" average utr sizes
    # should adjust %cds bad for cds size, for large cds> 10k, pcds >=90%, for small cds < 300, pcds may be 33%
    if($tqual =~ m/utrbad/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRBAD; } 
    elsif($tqual =~ m/utrpoor/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRPOOR; }
    else { $ispoor |= ($pcds < $BAD_CDSUTR)? kAAUTRBAD: ($pcds < $OK_CDSUTR)? kAAUTRPOOR : 0; }
    # if($butr <= 300) {} elsif($pcds < $BAD_CDSUTR){ $ispoor |= kAAUTRBAD; } elsif($pcds < $OK_CDSUTR) { $ispoor |= kAAUTRPOOR; } 
  }
  
  if($pidal =~ m/altmap\d+xeq/ and $cla =~ /althi1|part/) {
    $ispoor |= kAADUP;  # add new bitflag kDUPEXONS ?
  }
  
  if($aamainpi > $AADUP_IDENT) { # TEST aacluster added info on dup prots; AADUP_IDENT=98 default
    $cla.="a2" unless($cla=~/hi[012]/); # hi1/a2 redundant info
    $tbits.=",aadup:$aamainid";
    $ispoor |= kAADUP if($cla =~ /althi/); # only for /althi|parthi/
  }

  my $MINBLASTSCORE= 60; # aablast only? bitscore always?
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
    if($aw >= $AAMINBAD) { $ispoor = $ispoor & NOTPOORBAD; } ##  ^ (kAAGAPS+kAAUTRPOOR+kAAUTRBAD);  ##?? ^ XOR bad
    elsif($aw >= $AAMINPOO) { $ispoor = $ispoor & NOTPOOR; } ##  ^ (0+kAAUTRPOOR);
  }
}

      ## classifyFullTr() parts .. do AFTER find cross-locus alignments
  my $fclact="";
  if($fclass) {  
    
    ## use ispoor ...
    # my $thoval= $aablast{$tid} || 0;  # my($qhoscore,$qhoref)=split",",$thoval;
    # my $pcds= ($tw>0) ? 100*$tc/$tw : $OK_CDSUTR;
    # my $tcdsisbad=0;
    # if($thoval) { $tcdsisbad=0; }
    # elsif($tcode=~/^Noncode/) { $tcdsisbad=$qcode; }
    # elsif($tc < 2*$MINCDS and $pcds <= $BAD_CDSUTR) { $tcdsisbad=2; }
    # elsif($tc < $MINCDS and $pcds < $OK_CDSUTR) { $tcdsisbad=1; }
    
    ## oneloc == cdsovcds or (utrovcds and tcdsisbad)
    ## twoloc == (utrovutr or cdsovutr) and not tcdsisbad
    ## ambig  == others?
    if($fclass =~ /cdsovcds/ or ($fclass =~ /utrovcds/ and $ispoor)) { 
      if($fclass=~/part/ or $ispoor) { $fclact="drop.$fclass"; } # frag/part
      else { $fclact="alt.$fclass"; } ## large set here
    } elsif($fclass =~ /utrovutr|cdsovutr/ and not $ispoor) {
      $fclact="locus.$fclass"; # keep? use orig locus class ?
    } elsif($fclass =~ /altof|mainof/) { #  and not $ispoor
      $fclact=$fclass; # "altof"; #?  mainof > mainlocus ??
    } else {
      $fclact="locusmaybe.$fclass"; # ispoor? includes .notover1
    }
  }
  $fclact ||= "notover2"; #?
  
  #  all main+noclass/notover (no utr overlaps) that are cds-ispoort (Noncode) should be kept as ncRNA
  
  my $claret="$cla/$fclact"; # return this?
  if($cla =~ /parthi|frag0aa/ or $fclact =~ /^drop/) {
    $keepdrop.= "drop";
    # $cla=$fclact unless($cla =~ /parthi|frag/);
  
  } elsif($fclact =~ /^locus/) { # locusmaybe?
    $keepdrop.= "okay"; #? includes notover1 + ispoor == noncode?
    if($ispoor and $cla =~ /main|noclass/) { $claret="noncode/$fclact"; }
    # $cla=$fclact;
    
  } elsif($cla =~ /althi/ or $fclact =~ /^alt/) {
    $keepdrop.= ($ispoor)?"drop":"okay";  # ispoor vs main size?
    # $cla=$fclact unless($cla =~ /alt/);
    
  } elsif($cla =~ /main/) { # mainlocus? other term for best of overlocus aligns?
    if($ispoor and $fclact =~ /notover1/) { $ispoor=0; $claret="noncode/$fclact"; }
    $keepdrop.= ($ispoor)?"drop":"okay";  # ??
    
  } elsif($cla =~ /noclass/) {
    if($ispoor and $fclact =~ /notover1/) { $ispoor=0; $claret="noncode/$fclact"; }
    $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop? ditto to altmid
    
  } else { # other altmid/low ; part ?
    $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop this class, ispoor maybe uniq prot: maybeok?
  }
  
  
  my $okay;
  if($keepdrop =~ /drop/) {
    if($keepdrop =~ /okay/) { $okay= "maybeok"; } else { $okay= "drop"; }
  } else {
    $okay= "okay";
  }
  
  # now added in trevidence() .. flag to skip this?
  unless($tbits=~/chrmap:/) { $tbits.= ",chrmap:$mapqual,$maploc" if($maploc and $mapqual); }
  unless($tbits=~/pflag:/) { $tbits.= ",pflag:$ispoor"; } # DEBUG info, but keep
  unless($tbits=~/^aaref:/) { $tbits="aaref:$tbits" unless($tbits=~/^0,0/); } 
  
  ## add feq: back for overlocus cullExonEq() test??  
  if(defined $eqflag{$tid} and not($tbits =~ /feq:/)) { 
    my $eqfl= $eqflag{$tid}{$qid}||""; 
    if($eqfl) { $eqfl="$qid/$eqfl,"; }
    my @q= grep{ $_ ne $qid } sort keys %{$eqflag{$tid}}; 
    $eqfl .= join",",map{ "$_/".$eqflag{$tid}{$_} }@q;  
    $tbits.= ",feq:$eqfl";
  }
  
  ## maybe add for output cds.qual stats tid/qid (debug?)
  return (wantarray) ? ($tid,$okay,$claret,$qid,$pidal,$tqual,$tbits) : $okay;
}

=item example overlocusclass

  see asmrna_dupfilter3_overlocsub.pm

eg5:evg12aedesmpub.m4class
Aedesg12bEVm000320t19   okay    althi/alt.cdsovcds      Aedesg12bEVm000215t3    100/7/. 1171,73%,complete       aaref:aaref:2110,AAEL010606-PA,pflag:0,pflag:0
Aedesg12bEVm000320t20   okay    althi/alt.cdsovcds      Aedesg12bEVm000215t3    100/7/. 1171,72%,complete       aaref:aaref:2110,AAEL010606-PA,pflag:0,pflag:0
    100/7 << tiny cds over 7%, should it be called alt?
.. these are oid-alts, however, should use that
Aedesg12bEVm000215t1	Aedesg1EVm000222t1	Aedesg12bEVm000215	1	main	2232,90%,complete	100/100/./altmap100xeq	aaref:4532,AAEL018134-PA,refgood,chrmap:100a,100i,6699l,9x,supercont1.241:1295228-1328778:-,pflag:0,feq:Aedesg1EVm000222t2/altmapxe100.100,
Aedesg12bEVm000215t2	Aedesg1EVm000222t2	Aedesg12bEVm000215	2	althi1maybe	2232,88%,complete	100/100/./altmap100xeq	aaref:4533,AAEL018134-PA,refbest,chrmap:100a,100i,6699l,9x,supercont1.241:1295228-1328778:-,pflag:8,feq:Aedesg1EVm000222t1/altmapxe100.100,
Aedesg12bEVm000215t3	Aedesg1EVm000222t3	Aedesg12bEVm000215	3	althi	2084,85%,partial3	100/94/.	aaref:3944,AAEL018134-PA,chrmap:100a,100i,6160l,9x,Spl:5%,supercont1.490,supercont1.241:1295550-1320691:-,pflag:0
>> diff locus? m000215t3 is split-map problem case, 0320t19 may be accurate call, diff from  m000215t[12]
>> should drop any t19 problem cases, larger alts ok
Aedesg12bEVm000320t19	Aedesg1EVm000222t4	Aedesg12bEVm000320	19	cullalthi1	1171,73%,complete	100/92/./altmap52	aaref:2110,AAEL010606-PA,chrmap:99a,100i,3512l,9x,supercont1.490:478821-528287:-,pflag:0,feq:Aedesg1EVm000328t1/altmap52.36,Aedesg1EVm000222t3/altpar9.0.0,Aedesg1EVm000222t5/altmapxe99.100,Aedesg1EVm000328t10/altmap84.73,Aedesg1EVm0
Aedesg12bEVm000320t1	Aedesg1EVm000328t1	Aedesg12bEVm000320	1	main	2012,81%,complete	100/95/.	aaref:3698,AAEL010606-PA,refgood,chrmap:99a,100i,6038l,22x,Spl:2%,supercont1.1301,supercont1.490:463365-690982:-,pflag:0
Aedesg12bEVm000320t2	Aedesg1EVm000328t2	Aedesg12bEVm000320	2	althi	2008,81%,complete	100/93/.	aaref:3731,AAEL010606-PA,refgood,chrmap:97a,100i,6026l,21x,supercont1.490:463365-690982:-,pflag:0
Aedesg12bEVm000320t3	Aedesg1EVm000328t3	Aedesg12bEVm000320	3	althi	1999,81%,complete	100/95/.	aaref:3846,AAEL010606-PA,refbest,chrmap:100a,100i,5999l,22x,supercont1.490:463365-690982:-,pflag:0
Aedesg12bEVm000320t12	Aedesg1EVm000328t12	Aedesg12bEVm000320	12	althi	1251,74%,complete	100/93/.	aaref:2342,AAEL010606-PA,chrmap:99a,100i,3752l,11x,supercont1.490:478821-550496:-,pflag:0
.. cds.qual
Aedesg12bEVm000215t1	6699	0	2232,90%,complete	7387	686-7384:.	Aedesg1EVm000222t1,tidbaedes2sr6feo2ridk97Loc29099
Aedesg12bEVm000215t2	6699	0	2232,88%,complete	7602	632-7330:.	Aedesg1EVm000222t2,tidbaedes2sr4mao2ridk41Loc45077
Aedesg12bEVm000215t3	6160	92	2084,85%,partial3	7299	1047-7298:.	Aedesg1EVm000222t3,aedes2sr4ma4dsoapk25loc10166t3
Aedesg12bEVm000320t12	3752	4	1251,74%,complete	5025	31-3786:.	Aedesg1EVm000328t12,aedes2fe6nm4dsoapk27loc22685t20
Aedesg12bEVm000320t19	3512	4	1171,73%,complete	4770	16-3531:.	Aedesg1EVm000222t4,aedes2fe6nm4dsoapk27loc22685t32

  
=cut

=item defaults

$AAMIN =$ENV{aamin}||30; #was 40;   # for aacomplete, utrok

$TINYALN = 20; # %align, was 25; was $ENV{mina}||50; ignore less than this == MINAL
$TINYALNBASES= $ENV{MINALIGN}||90; # was $MINALIGN= 90; # REUSE//NOT USED NOW; use just TINYALN; change to several levels
$TINYALNCDSLEN= 3*$TINYALNBASES;
$MINCDS = $ENV{MINCDS} || 90; # or 3*MINAA  

$ALTFRAG= 0.5;

$PHIALN=$ENV{ahi} || 98; # was 99; was 65 !!;  DROP THIS; use mina?
$NHIALN=$ENV{nhi} || 9; # what? for short cds cancel few base diff < PHIALN

$PHI = $ENV{phi} ||99; 
$PMID= $ENV{pmid}||90; 
$PLOW= $ENV{plow}||80;  

$OK_CDSUTR_default= 60;
$BAD_CDSUTR_default= 30; # % CDS/trlen

=cut

sub overlocus_eqgene {
  my($qd,$td,$pal)=@_;
  my($palmap,$antiflag,$eqflag)=(0,"","");
  
  $palmap= $eqgenes->{$qd}{$td}||0; # note this is pct of td aligned to qd

  ## oids: $palmap= $eqgenes->{$qd}{$td}||0; # note this is pct of td aligned to qd

  ## problem here, when oidof(qid) sameloc oidof(td) wont have eqgene entry??
  ## cancel this, oid-alts do have eqgene entries.. but some dont overlap  
  #   if(0 and %oidsof) {
  #     # $td= $oids{$td}||$td; $qd= $oids{$qd}||$qd; 
  #     ## this isn't good enough test, paralts diff loci show up w/ same oid
  #     ## need glocation test of sameoid, add to eqgene tables?
  #     my ($qoid)= split",", $oidsof{$qd};
  #     my ($toid)= split",", $oidsof{$td};
  #     my($qog,$tog)=map{ my $g=$_; $g=~s/t\d+$//; $g; } ($qoid,$toid);
  #     if($qog and $qog eq $tog) {
  #       $antiflag .="/altoid";
  #       $palmap= $pal if($palmap<$pal); # if($palmap <= 0) ?
  #       # return( $pal,$palmap,$antiflag );
  #     }
  #   }
  
  ## if($neqgene>0) 

  ## simple way:
  # if($palmap>0) { $antiflag .="/altmap$palmap"; $eqflag="altmap$palmap";
  #	if(1 && $palmap > $pal) { $pal=$palmap; } ##? not this change EQGENE_OVERRIDES_ALN
  #}

  ## hard way:
  if($palmap>0) { 
    my $xeq= $eqexons->{$qd}{$td}||0; # only for ($neqexon>0)
    my $XPHI = 95; my $XPLO= 3;
    # xeq care about (a) xeq>= identity == not alt but redundant, 
    #   (b) xeq <= noxover, not alt but paralog maybe, (c) middle = usual alt
    
    #* set aln=0/TINYALNBASES when pal=0/TINYALN
    #*? change 'paralt' to 'altpar' to avoid other evg parse problems?
    if($neqexon>0 and $xeq >= $XPHI) { 
      #n# $antiflag .="/altmapxe$palmap.$xeq"; #?? ugh, want this way for other uses
      $antiflag .="/altmap$palmap"."xeq/$qd"; #orig..
      $eqflag="altmapxe$palmap.$xeq"; # for report? or use only antiflag?
      $pal=$palmap if($palmap >= $PHI and $pal < $PMID);  # FIXME bad annot
      }
    elsif($neqexon>0 and $xeq < $XPLO) { # problems here.. dont reset pal yet
      $antiflag .="/altparx$pal/$qd";  #later? $pal=3; $aln=13;
      $eqflag="altparx$palmap.$xeq"; # for report
      }
    else { 
      (my $qg=$qd)=~s/t\d+$//;  
      unless($td=~/^$qg/) {
        $antiflag .="/altmap$palmap/$qd"; #orig .. add .$xeq";
        $eqflag="altmap$palmap.$xeq"; # for report
        $pal=$palmap if($palmap >= $PHI and $pal < $PMID);  # FIXME bad annot
        }
      }
    }
    
  # elsif($pal>0 and exists( $eqgenes->{$td}) and exists( $eqgenes->{$qd})) { # defined or exists??
  #   ## bad here for nomap alts ** why isnt defined/exists working? for ids not in eqgene
  #   if($gmapqual->{$qd} and $gmapqual->{$td}) {
  #     my $revmap= $eqgenes->{$td}{$qd}||0; 
  #     ## unmapped gene cases are problem here.. need what? defined ($eqgenes->{td}) 
  #     $antiflag .="/altpar$pal" if($revmap<3);  # $pal=3; $aln=13; #? change later?
  #     $eqflag="altpar$pal.$palmap.$revmap"; # for report
  #   } else {
  #     #lots# warn "#DBG: bad eqgene $qd,$td\n" if($debug);
  #   }
  # }
    
  ##}

  $eqflag{$td}{$qd}=$eqflag if($eqflag);
  return( $pal,$palmap,$antiflag,$eqflag );
}


use constant minXEQ => 99; # exon-equal, high to ensure keep valid alts 
use constant minXCP => 95; # cds-align min for cullExonEq
use constant minXAL => 90; # chr-align min for cullExonEq

=item cull1ExonEq 
  from trclass2mainalt.pl, change for 1 compare

  bugs:
  FIXME: this is culling differing alts, at least by aasize,maploc, gmap-nexon, ..
  use note attrs more. chrmap:
  these are md-cullers, but not same alt form as longer alts
  need also note param of md culler 
   Anofunz4kEVm000028t7	okay	mainlocus/notover2	Anofunz4kEVm000028t1	100/68/.	2766,93%,complete	
     aaref:4249,AGAP002523-PA,refbest,chrmap:99a,99i,8297l,9x,KB669425:23470-33108:-,pflag:0
   Anofunz4kEVm000028t22	maybeok	parthi/altof	Anofunz4kEVm000028t1	100/70/.	2064,69%,complete	
     aaref:4022,AGAP002523-PA,refgood,chrmap:100a,99i,6195l,8x,KB669425:25506-33108:-,pflag:0
   
  Anofunz4kEVm000028t4	drop	althicull/altof	Anofunz4kEVm000028t1	100/59/.	3292,86%,complete	
    aaref:3773,AGAP002523-PA,refok,chrmap:97a,99i,9879l,15x,Spl:28%,KB668673,KB669425:24898-33108:-,pflag:0,feq:Anofunz4kEVm000028t22/altmapxe100.100
  Anofunz4kEVm000028t5	drop	althicull/altof	Anofunz4kEVm000028t1	100/76/.	2786,93%,complete	
   aaref:3764,AGAP002523-PA,refok,chrmap:99a,99i,8361l,10x,KB669425:23470-33108:-,pflag:0,feq:Anofunz4kEVm000028t22/altmapxe100.100,Anofunz4kEVm000028t6/altmapxe99.100,Anofunz4kEVm000028t7/altmapxe99.100
  Anofunz4kEVm000028t6	drop	althicull/altof	Anofunz4kEVm000028t1	100/68/.	2785,93%,complete	
    aaref:4238,AGAP002523-PA,refgood,chrmap:99a,99i,8354l,10x,KB669425:23470-33108:-,pflag:0,feq:Anofunz4kEVm000028t7/altmapxe99.100

  ** maybe two loci here, diff aarefs, splitmaps, diff scafs
  >> revised cull1x better, only 028t5 is drop.cull, 028t6 becomes new main (refgood),
     028t22 becomes maybeok	parthi/altof t1/t6
     t1,t2 longer aa, but smaller aaref:hoscore become alts of t6 ?? could be poor AGAP prot
     
=cut

sub cull1ExonEq {  
  my($md,$mdcl,$td,$tcl,$feq)=@_;  # ,$note
  my($culled,$galn)=("",0); 
  return "" if($mdcl=~/cull/);

  #now from caller: my($feq)= ($note=~/feq:([^;:\s]+)/)?$1:0;
  
  ## overlocus_eqgene now returns these
  ##    $antiflag .="/altmap$palmap.xeq"; #orig..
  ##    $eqflag{$td}{$qd}="altmapxe$palmap.$xeq"; # for report? or use only antiflag?

  ## my $note= $notes->{$td};
  my($tispoor,$tqual,$note) = trevidence($td,1);    
  my($mispoor,$mqual,$mnote)= trevidence($md,1);   

  return "" if( ($note=~/refbest/) or not ($tcl=~/althi|parthi/) 
        or ($mispoor>0 and $tispoor==0) );
  
  # pick out all of note parts, mnote also ??
  my($tgaln,$txn,$tspl,$mgaln,$mxn,$mspl,$cmap)=(0) x 9;
  ($cmap)= $note =~ m/chrmap:([^;:\s]+)/; # ,scaff[:]span[:]or .. not needed 
    ($tgaln)= ($cmap =~ m/(\d+)a,/)?$1:0;
    ($txn)  = ($cmap =~ m/(\d+)x,/)?$1:0;
    ($tspl) = ($cmap =~ m/Spl:(\d+)/)?$1:0;
  ($cmap)= $mnote =~ m/chrmap:([^;:\s]+)/; 
    ($mgaln)= ($cmap =~ m/(\d+)a,/)?$1:0;
    ($mxn)  = ($cmap =~ m/(\d+)x,/)?$1:0;
    ($mspl) = ($cmap =~ m/Spl:(\d+)/)?$1:0;
  
  my @xe= grep /altmapxe/, split",",$feq; # 1 only here? xd == md
  for my $xe (@xe) {
    my($xcp,$xev)= $xe =~ m/altmapxe(\d+).(\d+)/; # ugh, diff syntax now from overloc
    next unless($xev and $xev >= minXEQ and $xcp >= minXCP);
    next if($tgaln < minXAL or ($mxn and $txn > $mxn) or $tspl > 9);
    $culled= $tcl."cull"; # was "cull$tcl/$xe";
    return $culled; # or last;
  }

  return $culled;  
}

## debug add for asmrna_dupfilter3_overlocsub.pm, dont need these
sub _minb { my($x,$y)=@_; return ($y < $x) ? $y : $x; }
sub _maxb { my($x,$y)=@_; return ($y > $x) ? $y : $x; }


sub overlocusclass {
  my($outh, $infile, $insorted, $iscdsalign)= @_;
  $iscdsalign||=0; ## default NOT $IS_CDSALIGN

  ## ?? infile == blast table -outfmt 7, this sort won't work..
  ## redo alntab to include offsets $qab,$qae,$tab,$tae ?
  
  my $ALNSORTORD='-k2,2nr -k7,7nr -k3,3n -k6,6n -k1,1 -k4,4';  #add vQtlen/3 before vTt/6 IDs to order ties
  unless($iscdsalign) { # sort tlen, not clen 
    # .. not sure what is right, want cds to play role in best choice
    # .. input blast align k7 is for tlen, not clen, but want to choose best by long clen > long tlen,talign
    #t1.NO: $ALNSORTORD='-k3,3nr -k7,7nr -k2,2nr -k6,6nr -k1,1 -k4,4';  # ^Qtlen,^Align,^Qclen,vTtlen,Qid,Tid
    $ALNSORTORD='-k2,2nr -k7,7nr -k3,3n -k6,6n -k1,1 -k4,4';  # test same as iscdsalign  
  }
  
  my($inh, %class,%bestmatch,%ismain, %havepair); 
  if(ref($infile)) { $inh= $infile; } # STDIN,.. sort ??
  else {
    if($insorted) { open(IN,$infile) or die "reading $infile"; $inh=*IN; }
    else { open(IN,"sort $ALNSORTORD $infile |") or die "reading $infile"; $inh=*IN; }
  }
  
  my($lastd)=("");
  while(<$inh>) { # maybe local table, sorted, or from file
    next if(/^Qid|^\W/); chomp; 
    my @v= split"\t"; 
    my($qd,$td,$qc,$qw,$qcoff,$qcode,$qcoffb,$qcoffe, $qab,$qae,$qspan,
       $tc,$tw,$tcoff,$tcode,$tcoffb,$tcoffe, $tab,$tae,$tspan, $tor, 
       $aln,$iden,$bits,$pidn,$mis,$gap) = (0) x 30;
    my($isfrag,$aclass)=("") x 4;
    my($alnmax,$pid,$pal,$samesize)= (0) x 10;

    ## -bits == reverse align : bits not used here for scoring.. best choice other than adding column
   
if(0) {    ## no good for sort input..
    ($qd,$td,$bits,$pidn,$aln,$mis,$gap,$qab,$qae,$tab,$tae)= @v[0,1,-1,2,3,4,5, 6,7,8,9]; #blast table ?
         # 6-9 =  q. start, q. end, t. start, t. end, 

    # $iden= _maxb(0,$aln-$mis-$gap); # better? $iden = $pidn/100 * $aln;
    $iden= int($pidn/100 * $aln); # or $aln-$mis-$gap
    if($tab>$tae) { $bits=-$bits; $tor=-1; ($tab,$tae)=($tae,$tab); }
    ($qc,$qw,$qcoff,$qcode) = map{ $sizeval{$qd}{$_}||0 } qw( cdsize  trsize cdsoff  codepot);
    ($tc,$tw,$tcoff,$tcode) = map{ $sizeval{$td}{$_}||0 } qw( cdsize  trsize cdsoff  codepot);
    
} else {
    ($qd,$qc,$qw,$td,$tc,$tw,$aln,$iden,$bits,$qspan,$tspan)= @v;  #*** need blast align start,end vals
    ($qab,$qae)=split "-",$qspan; 
    ($tab,$tae)=split "-",$tspan; 
    $pidn= ($aln<1)? 0: int(0.5+ 100*$iden/$aln); $pidn=100 if($pidn>100);

    ($qcoff,$qcode) = map{ $sizeval{$qd}{$_}||0 } qw( cdsoff  codepot);
    ($tcoff,$tcode) = map{ $sizeval{$td}{$_}||0 } qw( cdsoff  codepot);
}

    # skip to next early on these
    if($td eq "self" or $qd eq $td) { # putspan fix for no-align cases
      $td=$qd; $bestmatch{$qd}="$td,100/100" unless($bestmatch{$qd});
      next;   
    }
    
    $qcoff= $cdsoff{$qd} unless($qcoff);
    $tcoff= $cdsoff{$td} unless($tcoff);

=item tcoff (only?) reading bug
    td == self bug?
  Use of uninitialized value $tcoffe in subtraction (-) at evigene/scripts/rnaseq/asmrna_dupfilter3.pl line 1224, <IN> line 2541139.
  Use of uninitialized value $tcoffb in subtraction (-) at evigene/scripts/rnaseq/asmrna_dupfilter3.pl line 1224, <IN> line 2541139.

=cut

=item special case near identicals
  .. maybe best handled here with swap to best, but need all evidence 1st.

  ## also check if class{$td}/class{$qd} ??
  ## do AFTER 1st tc > qc swap and new tcaln calcs..
  my $nearlysame=(abs($tc-$qc) < 30 and $tcaln >= 0.90*$qc)?1:0;    
  if($nearlysame) { 
    my $doswap=0;
    my($tispoor,$tqual,$tbits)= trevidence($td); # ,$tcla =~ m/althi|parthi/
    my($qispoor,$qqual,$qbits)= trevidence($qd); #,$qcla =~ m/althi|parthi/
    my($tho,$qho)=map{ my($ho)= (m/aaref:(\d+)/)?$1:0; $ho; } ($tbits,$qbits);
    if($tho > 5+$qho) { $doswap=1; } elsif($qho > 5+$tho) { $doswap=0; }
    elsif($tispoor > $qispoor) { $doswap=1; } 
    else { $doswap=0; }
    if($doswap) {
      ($qd,$td)=($td,$qd); 
      ($qc,$qw,$tc,$tw)=($tc,$tw,$qc,$qw); 
      ($qcoff,$qcode,$tcoff,$tcode)=($tcoff,$tcode,$qcoff,$qcode);
      ($qab,$qae,$tab,$tae)= ($tab,$tae,$qab,$qae);
      ($qgd,$tgd)=($tgd,$qgd); 
      ($qcoffb,$qcoffe,$tcoffb,$tcoffe)=($tcoffb,$tcoffe,$qcoffb,$qcoffe);
      ($qcaln,$qualn,$tcaln,$tualn)= ($tcaln,$tualn,$qcaln,$qualn);
    }
  }
  
eg2:evg12aedesmpub.m3class 
Aedesg12bEVm007420t1	okay	main/	Aedesg12bEVm068781t1	100/92/.	473,92%,complete	aaref:969,AAEL011510-PA,refgood,pflag:0
  ^^vv m068781t1 (3 alts) ~= m007420t1 (1 alt), same size, 781t has bettr aaref score, should replace m007420t1 
Aedesg12bEVm068781t1	okay	althi/alt.cdsovcds	Aedesg12bEVm007420t1	100/92/.	473,92%,complete	aaref:975,AAEL011510-PA,refbest,pflag:0
eg3:
Aedesg12bEVm004960t1	okay	main/	Aedesg12bEVm068777t1	99/91/.	605,91%,partial5	aaref:975,AAEL014089-PA,refgood,pflag:0
  ^^vv m068777t1 better aaref than m004960t1, same size, same locus
Aedesg12bEVm068777t1	okay	althi/alt.cdsovcds	Aedesg12bEVm004960t1	99/91/.	605,90%,partial5	aaref:989,AAEL014089-PA,refbest,pflag:0


=cut

    ## Q presumed larger than T (sort), if need swap, do before all other calcs.
    use constant SWAPBIGCDS => 1; 
    if(SWAPBIGCDS && $tc > $qc) { 
      ($qd,$td)=($td,$qd); 
      ($qc,$qw,$tc,$tw)=($tc,$tw,$qc,$qw); 
      ($qcoff,$qcode,$tcoff,$tcode)=($tcoff,$tcode,$qcoff,$qcode);
      ($qab,$qae,$tab,$tae)= ($tab,$tae,$qab,$qae);
      }
 
    my($qgd,$tgd)=map{ my $g=$_; $g=~s/t\d+$//; $g; } ($qd,$td);
    my $isaltof = ($qgd eq $tgd)?1:0; # fixme, fclass == altof
    my $ismainof= ($isaltof and $td=~/t1$/)?1:0; # want this? avoid reclassing t1mains, for now
    
use constant doALTCULL => 1;

    if(not doALTCULL and $isaltof) { 
      # dont care here for known alts??
      ## DO this after below cds-over aclass = althi 
      ## add cullExonEq() test? but need input with gmap feq:ID/altmapxe87.100 .. 
      # if(0 and $neqgene>0) {
      #   my( $palg,$palmap,$antiflagg,$eqflagg ) = overlocus_eqgene($qd,$td,$pal) ; # this sets feq altmapxe..
      #   ## $antiflag .="/altmapxe$palmap.$xeq"."xeq"; 
      #   if($eqflagg =~ m/altmapxe\d/) {
      #     ## also need aaref:val,refbest/good/ok,  $note=~/feq:([^;:\s]+)/,  in antiflag == note
      #     ## aaref from trevidence
      #     my($ispoor,$tqual,$tbits)= trevidence($td,1); # $cla =~ m/althi|parthi/
      #     my($culls)= cull1ExonEq($qd,$td,"althi","$tbits,feq:$qd/$eqflagg"); 
      #     # my $culls= ($CULLXEQ) ? cullExonEq($md,\@ad,\%alt,\%notes) : {}; #?? here
      #   }
      # }
      
      next;
    }
    
    ($qcoffb,$qcoffe)=split "-",$qcoff;
    ($tcoffb,$tcoffe)=split "-",$tcoff;
 
    my($qcaln,$qualn,$tcaln,$tualn)=(0) x 4; # split align to cds/utr portions, for each seq
    if($qae < 1) { #  or $qcoffe < 1 ## fail??
      $qcaln= _minb($aln, $qc);  $qualn= _minb($aln, $qw - $qc);
    } else {
      $qcaln= _maxb(0, _minb($qae,$qcoffe) - _maxb($qab,$qcoffb));
      $qualn= _maxb(0, $qae - $qcoffe) + _maxb(0, $qcoffb - $qab);
    }
    if($tae < 1) { #  or $tcoffe < 1
      $tcaln= _minb($aln, $tc);  $tualn= _minb($aln, $tw - $tc);
    } else {
      $tcaln= _maxb(0, _minb($tae,$tcoffe) - _maxb($tab,$tcoffb));
      $tualn= _maxb(0, $tae - $tcoffe) + _maxb(0, $tcoffb - $tab);
    }

    ## Q presumed larger than T 
    ## do AFTER 1st tc > qc swap and new tcaln calcs..
    ## also check if class{$td}/class{$qd} ??
    ## problems for aaref best but not nearlysame, for cdsover cases. want to swap anyway
    
    my $nearlysame=(abs($tc-$qc) < 30 and $tcaln >= 0.90*$qc)?1:0;  #?? tcaln >= 0.50*qc   
    my($thoval)= split",", ($aablast{$td} || "0");
    my($qhoval)= split",", ($aablast{$qd} || "0");
    $nearlysame=1 if($thoval > 5+$qhoval); # force swap check .. small diff?
    
    if($nearlysame) { 
      my $doswap=0;
      my($tispoor,$tqual,$tbits)= trevidence($td); # ,$tcla =~ m/althi|parthi/
      my($qispoor,$qqual,$qbits)= trevidence($qd); #,$qcla =~ m/althi|parthi/
      my($tho,$qho)=map{ my($ho)= (m/aaref:(\d+)/)?$1:0; $ho; } ($tbits,$qbits);
      my $hodiff=($tispoor and not $qispoor)?15:5;
      $hodiff += 10 if($isaltof and $qd=~/t1$/); # avoid swap out t1main
      if($tho > $hodiff+$qho) { $doswap=1; } # any quals but hoval to check?
      elsif($qho > 5+$tho) { $doswap=0; }
      # elsif($tispoor < $qispoor) { $doswap=1; } 
      elsif( ($qispoor & (kAATINY+kAAUTRBAD)) > ($tispoor & (kAATINY+kAAUTRBAD)) ) { $doswap=1; } 
      else { $doswap=0; }
      if($doswap) {
        ($qd,$td)=($td,$qd); 
        ($qc,$qw,$tc,$tw)=($tc,$tw,$qc,$qw); 
        ($qcoff,$qcode,$tcoff,$tcode)=($tcoff,$tcode,$qcoff,$qcode);
        ($qab,$qae,$tab,$tae)= ($tab,$tae,$qab,$qae);
        ($qgd,$tgd)=($tgd,$qgd); 
        ($qcoffb,$qcoffe,$tcoffb,$tcoffe)=($tcoffb,$tcoffe,$qcoffb,$qcoffe);
        ($qcaln,$qualn,$tcaln,$tualn)= ($tcaln,$tualn,$qcaln,$qualn);
      }
    }
    
    #?? aln as only-cds-overlap? need another val: qcaln x tcaln
    # overloc classes from q,t CDS v UTR aligns
    #   .. trfrag overlaps CDS of trgood > drop/relcass as partof or altof trgood.
    #     .. unless trfrag.UTR only over trgood.CDS and trfrag.CDS has good code
    #   .. trfrag overlaps only UTR of trgood, depends on trfrag coding quals: other locus or partof trgood
    #     .. trfrag.CDS overlaps trgood.UTR, trfrag > partof:trgood unless has strong code quals
    #     .. trfrag.UTR only over trgood.UTR, trfrag > other locus (code/noncode?)
    #   .. trfrag overlaps trfrag, depends on CDS v UTR, code quals

    $pid= int(0.5+$pidn); # replaces: ($aln<1)?0: int(0.5+ 100*$iden/$aln); $pid=100 if($pid>100);
    my $alnfull=$aln;
    
    ## use "part" not "frag" qualifier?
    ## should this use only cdslen? $ispart=($tc < 0.50 * $qc)
    # my $ispart= ($tw < $ALTFRAG*$qw)?"part":""; # see below; altfrag == 0.5
    my $ispart= ($tc < $ALTFRAG*$qc)?"part":""; # see below; altfrag == 0.5
    
    #? $samesize=($tw == $qw and $tc == $qc)?2:($tc == $qc)?1:0;    
    # my($tispoor,$tqual,$tbits)= trevidence($td,$tcla =~ m/althi|parthi/);
    # my($qispoor,$qqual,$qbits)= trevidence($qd,$qcla =~ m/althi|parthi/);
    
    #============= overlocus classing ====================
    my $fclass="";
    # my $alnsig = ($alnfull >= 0.10 * _min($qw,$tw))?1:0;
    if($isaltof) { # ($qgd eq $tgd) 
      $fclass=($ismainof)?"mainof":"altof";  # handle these other way, but separate if NO/minor overlap? == paralogs??

    # } elsif($qcaln >= $TINYALNBASES) { # == 90b ** ADD pctalnmax == qc/tc; $TINYALN == 20 now, too big?
    } elsif($qcaln >= $TINYALNBASES and ($alnfull >= 0.10 * $tw)) { # == 90b ** ADD pctalnmax == qc/tc 
      if($tcaln < $NHIALN) { # == 9b
        $fclass="utrovcds$ispart"; # this must be >= TINYALNBASES ?
      }
      else { $fclass="cdsovcds$ispart"; } # == frag
      
    } elsif($qualn >= $TINYALNBASES and ($alnfull >= 0.10 * $tw)) {
      if($tcaln < $NHIALN) { $fclass="utrovutr$ispart";  } # this must be >= TINYALNBASES ?
      else { $fclass="cdsovutr$ispart"; } # == frag
    }
    unless($fclass) {
      if( $alnfull < $TINYALNBASES or ($alnfull < 0.10 * _min($qw,$tw)) ) {
        $fclass = "notover1"; 
      }
    }
    
    ## classifyFullTr() parts .. do AFTER find cross-locus alignments
    ## oneloc == cdsovcds or (utrovcds and tcdsisbad)
    ## twoloc == (utrovutr or cdsovutr) and not tcdsisbad
    ## ambig  == others?
    #   my $fclact="";
    #   if($fclass =~ /cdsovcds/ or ($fclass =~ /utrovcds/ and $ispoor)) { 
    #     if($fclass=~/part/ or $ispoor) { $fclact="drop.$fclass"; } # frag/part
    #     else { $fclact="alt.$fclass"; }
    #   } elsif($fclass =~ /utrovutr|cdsovutr/ and not $ispoor) {
    #     $fclact="locus.$fclass"; # keep? use orig locus class ?
    #   } else {
    #     $fclact="maybelocus.$fclass";
    #   }
 
   
    #============= cdsalign classing ====================
    #above# my $alnfull=$aln;
    $aln= _maxb($qcaln,$tcaln); # ?? which min/max or qcaln?
    
    # UPD1912: caller blastn -strand plus removes all -sense aligns
    my $antiflag= ($iscdsalign and $bits < 0) ? "/-sense" : "/."; ## bestmatch="id,99/89/-sense" for antisense?
    $samesize=($tc == $qc)?1:0; # OR tinyalndiff? abs($tc - $qc) < $NHIALN
    
    ## eqgene classifier, not for overlocus ?? or yes?
		
		$havepair{$qd}{$td}++; #?? change this to pal align val ?
		$lastd=$qd;

    ## ABOVE now; 
    # if($td eq "self" or $qd eq $td) { # putspan fix for no-align cases
    #   $td=$qd; $bestmatch{$qd}="$td,100/100" unless($bestmatch{$qd});
    #   next;   
    # }

    my $qcds= ($qw>0) ? 100*$qc/$qw : $OK_CDSUTR;
    my $tcds= ($tw>0) ? 100*$tc/$tw : $OK_CDSUTR;
    my $qutrbad= ($qcds >= $OK_CDSUTR)?0:1;  
    my $tutrbad= ($tcds >= $OK_CDSUTR)?0:1;  

    my($qsize,$tsize)= ($iscdsalign) ? ($qc,$tc) : ($qw,$tw);
    $alnmax= ($qsize>$tsize and $tsize>0)?$tsize:($qsize>0)?$qsize:$tsize;  

    #deflt: $ALTFRAG= 0.5; .. change this for overlocus?
    $isfrag= ($tsize < $ALTFRAG*$qsize)?"frag":"";
    
    $pal= ($alnmax<1)?0 : int(0.5+ 100*$aln/$alnmax); $pal=100 if($pal>100);
    # my $palq= ($qsize<1)?0 : int(0.5+ 100*$aln/$qsize); $palq=100 if($palq>100);
    # my $palt= ($tsize<1)?0 : int(0.5+ 100*$aln/$tsize); $palt=100 if($palt>100);
    
    my $skiptinyaln=0;
    my $tisaltpar=0; # treat like $skiptinyaln for now
    my $tinyalndiff= ((($alnmax - $aln) < $NHIALN)) ? 1 : 0; #yes: TEST3 && 

		$havepair{$qd}{$td}= $pal; #?? change this to pal align val for altpar ?

    # for both isaltof/notaltof, add antiflag to pidval.
    my( $eqpal,$eqpalmap,$eqantiflag,$eqflagg ) 
        = ($neqgene>0)? overlocus_eqgene($qd,$td,$pal) : (0,0,0,0);

    if($isaltof) { # $neqgene>0 and 
      $antiflag .= $eqantiflag if($eqpalmap>0);  # ismainof?? /$qd ?? wrong got self:
      # Aedesg12bEVm000004t3	drop	althicull/altof	Aedesg12bEVm000004t1	100/92/./altmap98xeq/Aedesg12bEVm000004t3
    } else {
      # my( $palg,$palmap,$antiflagg,$eqflagg ) = overlocus_eqgene($qd,$td,$pal) ;
      if($eqpalmap>0) { #always?
        $pal= $eqpal; $antiflag .= $eqantiflag;
      } else {
        # cancel == paralog/alt
        $tisaltpar=1;  $antiflag .= "/noteqgene"; # $skiptinyaln= 1;
        $fclass =~ s/cdsovcds/utrovutr/; #??
        # bug: Aedesg12bEVm002491t1	drop	noclass/drop.cdsovcdspart	Aedesg12bEVm000291t1	99/77/./noteqgene
      }
    }
    
    my $pidalnval="$pid/$pal$antiflag"; # old: $pid/$pal
     ## .. change these to use palq, palt : palign per q,t size ?
    $bestmatch{$qd}="$td,$pidalnval" unless($bestmatch{$qd});
    $bestmatch{$td}="$qd,$pidalnval" unless($bestmatch{$td});  
        
    my $qclass= $class{$qd}||"";
    if($samesize and $qclass =~ /alt/ and $ismain{$td}) { # check/stop alt1 <> alt2 assigns
      next if($bestmatch{$qd} =~ /$td,/); # problem here? for many equal bestmatch
    }

    if(UseTINYALNBASES) { ## == 1
      if($alnmax >= $MINCDS) { # recall alnmax is min(qsize,tsize); never skip tiny prots, MINCDS =~ 90 bases ?
        my $minbaseover= ($alnmax >= $TINYALNCDSLEN)? $TINYALNBASES : $alnmax * 0.10; ##$MIN_PCTOVERLAP;
        $skiptinyaln=($aln < $minbaseover)?1:0;
      } else {
        $isfrag= "frag" unless($isfrag); #or "tiny" ? # dont skip assign alt/frag to tiny cds, but maybe always set isfrag ?
      }
    } else { # old, pct-TINYALN is a problem for large cds  
      $skiptinyaln= ($pal < $TINYALN)?1:0;
    }
  
    # ** change for fclass, eg. isfrag not relevant for separate loci, no cds overlap..
    unless($fclass =~ /cdsovcds|altof/) { 
      $isfrag=""; 
      if($fclass =~ /utrovcds|utrovutr|cdsovutr/) {
        $skiptinyaln=1; # yes, means dont reclassify for no/tiny cds-overlap, only for cdsovcds
        }
      }

    # set alt aclass here, or not.
    if( $skiptinyaln or $tisaltpar) { } ## defer: $aclass="noalign$isfrag";  
    elsif( $tc == 0 and $qc > 0 and $pid >= $PLOW ) { $aclass="frag0aa"; } # tc==0 special case "frag0"
    elsif( $pid >= $PHI ) { 
      $aclass= ($isfrag)?"parthi":"althi"; # Primary classifier; change parthi to althipart ? althi$isfrag ?
      $aclass .= "1" if($pal >= $PHIALN or $tinyalndiff); # althi1 == hi-align + hi-ident, nearly same (identicals gone here)
      }
    elsif( $pid >= $PMID ) { $aclass="altmid$isfrag"; } # partmid for isfrag ?
    elsif( $pid >= $PLOW ) { $aclass="altlo$isfrag"; }  # partlo for isfrag ?

    if(doALTCULL and $isaltof and not $ismainof and $aclass =~ /althi|parthi/) { 
      ## DO this after below cds-over aclass = althi 
      ## add cullExonEq() test? but need input with gmap feq:ID/altmapxe87.100 .. 
      ##  $eqpal,$eqpalmap,$eqantiflag,$eqflagg 
      if($eqflagg) { # ($neqgene>0)
        # my( $palg,$palmap,$antiflagg,$eqflagg ) = overlocus_eqgene($qd,$td,$pal) ; # this sets feq altmapxe..
        if($eqflagg =~ m/altmapxe\d/) {  
          ## need aaref:val,refbest/good/ok from trevidence,  $note=~/feq:([^;:\s]+)/,  
          ## add trevidence(qd) to get chrmap mapqual,loc for compare .. here? or cull1ExonEq
          my $qclass= $class{$qd}||""; 
          ## now in cull1x: my($ispoor,$tqual,$tbits)= trevidence($td,1); # $aclass =~ m/althi|parthi/
          my($culls)= cull1ExonEq($qd,$qclass,$td,$aclass,$eqflagg); # "$tbits,feq:$qd/$eqflagg" 
          if($culls) { $aclass= $culls; } # done above:  $antiflag .= $eqantiflag; 
        }
      }

    }

    #yes: if(TEST3)
    if(my $ocl= $class{$td} and not $aclass =~ /cull/) {  
      #yes: if(TEST1602)   
      if($ocl =~ m/^alt|^part/) { $class{$td}= $aclass if($aclass eq "althi1"); next; } # test3
      elsif($ocl =~ m/^main/) { 
        #* problem of 2+mains w/ several essentially identicals,  break tie here **
        if($qclass =~ m/^main/ or $tisaltpar) { 
          # skip to below aclass set, test5, yes, do this also?
          } 
        else { next; } # * YES, need this : test4
      } 
    }

    # add/rep fclass to aclass
    if($aclass and $fclass) { $aclass.="/$fclass"; }
    elsif($fclass) { $aclass="noclass/$fclass"; } #?? main/ noclass/ ; change part/ to other
    ## m4class most okay set
    # 10941 okay	althi/alt   << these are presumed prior mistakes, cdsover for diff loci, qualify w/ last-oid equiv?
    # 13132 okay	main/notover
    # 44160 okay	noclass/notover
    # 23274 okay	part/locus  << noclass/locus ?

    if($aclass) { 
      my $qmain=$qd; my $attr= $pidalnval; ## "$pid/$pal$antiflag";

      if($class{$qd}) {  #  =~ /alt|part/
        my $qnext= $qd; 
        my $more= ($bestmatch{$qnext})?1:0;
        my %qdid=( $qd => 1, $td => 1);
        while($more) { 
          $more=0;
          my($qbest,$qpal)=split",",$bestmatch{$qnext};
          
          #yes: if(TEST1603)
          $qbest=0 if($ismain{$qnext} and not $ismain{$qbest}); # want always? maybe bad?
          
          if($qbest and not $qdid{$qbest}) {
            $qnext= $qbest; $qdid{$qbest}=1; 
            $more=1 if($class{$qnext} and $bestmatch{$qnext});
          }
        }
        if($qnext ne $qd) { $qmain= $qnext; $attr="$pidalnval/$qmain"; } #was $attr.="/$qmain";
        ## WARNING: now attr has /-sense or /. before /qmain;  evgmrna2tsa2.pl needs to know...
      }
      
      $ismain{$qmain}++;   
      $class{$td}=$aclass; #? problem here when class{td} == main set before, should not reset to alt
      #yes: if(TEST1602)  # this way works.  if($bestmatch{$td}) must be true from above
      if($bestmatch{$td} =~ /^$qd/) { $bestmatch{$td}="$qd,$attr"; }
      elsif(not $bestmatch{$td}) { $bestmatch{$td}="$qd,$attr"; }

      #yes: if(TEST3) #add to prevent 2+ circular alta <> altb <> altc with no main ..
      $class{$qmain}="mainlocus" unless($class{$qmain}); # should always be:  $qc >= $tc//$samesize..

    }
      
  } close($inh);

  # fclass: cdsovcds,utrovcds,cdsovutr,utrovutr + part/isfrag
  ## problem here using q/pal when those are tiny cds of mostly utr overlap (for main/noclass at least)
  ## bestmatch is no longer best cds match (with alts), but best cross-locus match
  # END:  print trclass table; add more fields to output: aaqual, aablast, tr,aa sizes?
  { my($q,$pal,$c,$d);
  foreach $d (sort keys %class) {
    ($q,$pal)=split",",$bestmatch{$d}; $c= $class{$d}||"altclass"; # any missing?
    my @cla= classifyFullTr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  foreach $d (sort grep{ not($class{$_}) } keys %ismain) {
    ## FIXME: "main" here doesn't mean that, these are the better compared to overlapped loci-tr
    ## change to "oklocus" ? okovlocus?  "drop" classing from other aspects like ispoor Noncode val
    ($q,$pal)=split",",$bestmatch{$d}; $c= "mainlocus"; # was "main";    
    my @cla= classifyFullTr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  foreach $d (sort grep{ not($class{$_} or $ismain{$_}) } keys %bestmatch) {
    ($q,$pal)=split",",$bestmatch{$d}; $c= "noclass";  
    my @cla= classifyFullTr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  }  
  
}

1;

