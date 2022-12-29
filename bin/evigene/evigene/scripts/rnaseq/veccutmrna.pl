#!/usr/bin/env perl
# veccutmrna.pl
## test/work script for evigene/scripts/evgmrna2tsa.pl

use FindBin;
use lib ($FindBin::Bin,"$FindBin::Bin/.."); # assume evigene/scripts/rnaseq/thisscript.pl
my $EVIGENES="$FindBin::Bin/.."; #??

use strict;
use cdna_proteins;

my $vecsuf= $ENV{vecsuf}|| "vector.tab2";
my $outsuf= $ENV{outsuf}|| "uvcut2.mrna";
my $namesuf= $ENV{namesuf}|| "names";
my $debug= $ENV{debug}|| 0;
#drop# my $UVGAP=  $ENV{uvgap}||0;   ## NO: leave nnn gaps for cut vec, preserve codon%3 size?
my $MINSIZE= $ENV{min} || 150;
my @inmrna= @ARGV;

my(%uvcut);
my($sufaa,$sufcds,$sufuncut);
$sufaa=$outsuf; unless($sufaa=~s/\.mrna/.aa/){ $sufaa.=".aa"; }
$sufcds=$outsuf; unless($sufcds=~s/\.mrna/.cds/){ $sufcds.=".cds"; }
$sufuncut=$outsuf; unless($sufuncut=~s/cut/uncut/){ $sufuncut.=".uncut"; }

# my $APPcdnabest= findevigeneapp("cdna_bestorf.pl"); # allow ENV/path substitutions?
my $APPcdnabest="$EVIGENES/cdna_bestorf.pl";

my( %genenames,%genenamepct,%genedbxref,%namedgenes,%cddnames); # for neworf tests. from evgmrna2tsa
 
foreach my $mf (@inmrna) {
  (my $nam = $mf) =~ s/\.mrna.*//;
  my $vtab="$nam.mrna.$vecsuf";
  my $outf="$nam.$outsuf";
  my $outaa="$nam.$sufaa"; my $outcds="$nam.$sufcds"; my $outuncut="$nam.$sufuncut";
  
  my $genenamef="$nam.$namesuf"; # what? evg path is above okayset/mrna, evgmrna2tsa.pl gets it
  unless ( -f $vtab ) { warn "#uvcut missing $vtab\n"; next; }
  warn "#uvcut $mf with $vtab to $outf\n"; # if $debug;
  
  my($ngot,$nin)= parse_genenames($genenamef);
  warn "#uvcut names $genenamef n=$ngot\n"; # if $debug;
  
  %uvcut=(); 
  open(F,$vtab); while(<F>){ my($id,$b,$e,$vt)=split"\t"; $uvcut{$id} .= "$b\t$e\t$vt\n"; } close(F);
  if( $mf =~ /\.gz$/ ) { open(IN,"gunzip -c $mf |"); } else { open(IN,$mf); }

  # FIXME: add OUTUNCUT for non-uvector set, later combined
  # ? add gapclean() / trimNNNends() ? from evgmrna2tsa.pl:putseq()
  # .. needs to adjust CDSoffset, maybe CDS if partial with NNN in ENDGAP; easier to recalc prot as w/ uvcut
  
  open(OUT,'>',$outf);  open(AAOUT,'>',$outaa);  open(CDSOUT,'>',$outcds); open(OUTUNCUT,'>',$outuncut); 
  my($id,$fa,$hd,$didput,$nput,$nskip,$nbadcut,$nnochange)=(0)x9;
  while(<IN>) {
    if(/^>(\S+)/) { my $d=$1; 
      if($id) { $didput= putf($id,$hd,$fa); $nput++ if($didput>0); $nbadcut++ if($didput<0); $nskip++ if($didput==0); } 
      $hd=$_; $fa=""; $id=($uvcut{$d})?$d:0; $nnochange++ unless($id); } 
    elsif($id) { chomp; $fa.=$_; } 
    print OUTUNCUT $_ unless($id);
  } 
  if($id) { $didput= putf($id,$hd,$fa); $nput++ if($didput>0); $nbadcut++ if($didput<0); $nskip++ if($didput==0); } 

  close(IN); close(OUT); close(AAOUT); close(CDSOUT); close(OUTUNCUT);
  warn "#uvcut nochange=$nnochange, ncut=$nput, ndrop=$nskip, nbadcut=$nbadcut to $outf\n"; # if $debug;
} # end @inmrna



 
sub putf { 
  my($id,$hd,$fain)= @_;

  my $retval=0;
  my $olen=length($fain);
  my($cb,$ce)= $hd =~ m/offs=(\d+).(\d+)/;
  my($oldaaq)= $hd =~ m/aalen=(\S+)/;
  my $oldorflen = 1 + (($cb>$ce)? $cb - $ce : $ce - $cb);
  my($fac,$fav,$facgap,$lue,$lub,$ucut,$rue,$rub,$vectype);
  $facgap=$fac=$fav=""; $lue=$lub=$ucut=$rue=$rub=0; 
  
  my @vec=split"\n", $uvcut{$id}; my $nv=@vec;
  
  ## fixmed: not wrong order, but overlapped;
  #MISORDER catfishvel3ik45Loc262t5 uvec 1/2, 528-561 .. 549-567
  #MISORDER catfishvel3ik45Loc513t4 uvec 1/2, 1-27 .. 15-34
  if($nv>1) { 
    @vec= sort { $a <=> $b } @vec; # paste overlaps
    my @vec2=(); my $ncut=0;
    for(my $i=@vec - 1; $i>0; $i--) { 
       my($ub,$ue,$vty)=split"\t",$vec[$i-1];
       my($xb,$xe,$xty)=split"\t",$vec[$i]; 
       next unless($xb>0 and $xe >= $xb); #? bad data
       if($xb <= $ue) { $ub=$xb if($xb<$ub); $vec[$i-1]=join("\t",$ub,$xe,$vty); $ncut++; }  
       else { unshift @vec2, $vec[$i]; }
       }
    @vec= ($vec[0], @vec2);  $nv=@vec;
  } #? ok
  
  for(my $i=0; $i<$nv; $i++) { 
    my($ub,$ue,$vty)=split"\t",$vec[$i]; 
    my $uw=1+$ue-$ub; 
    $ucut += $uw;  ## add ucutInCDS += nb for ue - ub in ce - cb
    $vectype= $vty unless($vectype); # take 1st only?
    if($i == 0) {  $rub=$ub; if($ub>1) { my $s=substr($fain,0,$ub - 1); $fac.=$s; $facgap.=$s; } }
    else { 
      my $wnu= $ub - 1 - $lue; if($wnu>0) { my $s= substr($fain,$lue,$wnu); $fac.=$s; $facgap.=$s; }
      warn "#MISORDER $id uvec $i/$nv, $lub-$lue .. $ub-$ue\n" if($ub < $lue);
      }
    
##?? add UVGAP only if in cds-span? if($ub < $ce and $ue > $cb) 
## For ~3 cases, this is WRONG, leave out UVgap and get better new aa
## In 1 case, moderate uvec cuts cds5 makes much worse prot; 
## socatfishv1k95loc34193t1 (origaa=446,full and 100% match zfish gene; cutaa=176,utrbad;)
## vecscreen finds ~39 bp align w/ 2+ mismatch in cds5span
## w/o uvcut, full align to alpha-1,6-mannosyl-glycoprotein 2-beta-N-acetylglucosaminyltransferase [Danio rerio]
## .. use mrna.ann.txt homol-align vs vecsreen-medium to decide not to cut?
## FIXmaybe: some vector spans surrounded by NNN ; cut those?

    ## try both, use neworf to pick best : seems to work  # if($UVGAP)
    if(1) { my $nug=3 + ($uw % 3); my $uvg= substr("nnnnnn",0,$nug); $facgap .= $uvg; }
    
    if($i == $nv-1) { $rue=$ue; if($ue < $olen) { my $s= substr($fain,$ue); $fac.=$s; $facgap.=$s; } }
    $fav.="n" if($i>0 and $ub>$lue+1); $fav .= substr($fain,$ub-1,$uw); 
    $lue=$ue; $lub=$ub; 
    }
  
  ## insert here: endtrimNNN($fac); .. ($facgap)  : this for all mRNA with NNN, not just uvecs
  ## vectype="endtrim" for non-UVector but endtrim; ..
  
  $fac =~ s/[Nn]+$//; $fac =~ s/^[Nn]+//; #? yes
  $facgap =~ s/[Nn]+$//; $facgap =~ s/^[Nn]+//; #? yes
  
  
  $vectype =~ s/ match//;  my $vecNotStrong= ($vectype =~ /Strong/i)?0:1;
  
  my $clen=length($fac); my $fl=""; 
  if($rub < $ce and $rue > $cb) { if($rub <= $cb+2 and $rue<$ce) { $fl.="cds5"; }
  elsif($rue >= $ce-2 and $rub > $cb) { $fl.="cds3"; } else { $fl.="cdsin"; } }
  if($rue > $olen-9) { $fl.="end3"; } elsif($rub <= 9) { $fl.="end5"; }
  
  #FIXME2: need bestorf -partial option (always?) so that uvcut doesnt trigger further chop for complete-aa
  #FIXME: here? add protein qual check: cdna_bestorf($fac) : 
  ## if( aasize << cutsize/3? AND vectype = Moderate? and aahomol > minalign) 
  ##   then cancel cutvec; keep orig w/ note; check homol-align also?
  
  my $expect_neworflen= $oldorflen - $ucut; # not quite right for ucut in UTR
  my ($newaahdr,$newaa,$newcdshdr,$newcds,
      $newaalen,$newpcds,$newcompl,$neworflen)= getbestorf($id,$hd,$fac,$clen, $expect_neworflen);
  my $note=""; my $badcut=0;

  if($neworflen < $oldorflen) {
    my $cglen= length($facgap);
    if($cglen > $clen) {
    my ($newaahdr1,$newaa1,$newcdshdr1,$newcds1,
       $newaalen1,$newpcds1,$newcompl1,$neworflen1)= getbestorf($id,$hd,$facgap,$cglen, $expect_neworflen);
    my $nga= $newaa1 =~ tr/X/X/; # my $ngc= $newcds1 =~ tr/N/N/;   
    my $newleng= $newaalen1 - $nga;
    if( $newleng > $newaalen) { 
      $note .= "framefixlen=$newleng-$newaalen; ";
      ($newaahdr,$newaa,$newcdshdr,$newcds, $newaalen,$newpcds,$newcompl,$neworflen) = 
        ($newaahdr1,$newaa1,$newcdshdr1,$newcds1, $newaalen1,$newpcds1,$newcompl1,$neworflen1);
      ($fac,$clen)= ($facgap,$cglen);   
      } 
    }
  }
        
  warn "#getbestorf: $id oldaa=$oldaaq; new.$newaahdr\n"; # if $debug;
        
  if($neworflen > $oldorflen) { # debug check these
    my $namepct= $genenamepct{$id} || 0; # is structured: 99%,123/345,678
    my $nameref= $genedbxref{$id} || 0; #  
    my $named= $genenames{$id} || 0; #  
    my($npct,$naln)= $namepct=~m/^(\d+)\D+(\d+)/;
    my $odiff= $neworflen-$oldorflen;
    warn "#named:  $id LONGER dorf=+$odiff, named=$namepct,$nameref,$named;\n" if($npct>=50 or $nameref =~ /:/); # if $debug; 
  }

  
  if($neworflen < $oldorflen - 9) { # check named align; also debug check names for neworflen >> oldorflen ?
    #  %genenames=%genenamepct=%genedbxref
    my $namepct= $genenamepct{$id} || 0; # is structured: 99%,123/345,678
    my $nameref= $genedbxref{$id} || 0; #  
    my $named= $genenames{$id} || 0; #  
    my $odiff= $neworflen-$oldorflen;
    
    ## ?? ignore 'CDD:' as main name, it lacks related species homolog
    ## .. not sure ignore CDD: is right yet.
    # my $okname= 1; 
    my $okname= (not $named or $named =~ /^CDD:|^hypothetical|^uncharacterized|^na$/i) ? 0 : 1;
    my($npct,$naln)= $namepct=~m/^(\d+)\D+(\d+)/;
    warn "#named:  $id dorf=$odiff, named=$namepct,$nameref,$named;\n" if($npct>=50 or $nameref =~ /:/); # if $debug;
    
    my($minpct,$minaln)= ($vecNotStrong) ? (66,50) : (90,90);
    ## for vecIsStrong, want to measure size of uvcut vs size of align? really need to know if uvcut is align-supported
    ## but dont yet have ref-align spans as inputs
    
    #FIXME: adjust badcut to accept Strong when also partial-cds at cut end?)
    # .. eg: shrimpt nbadcut=26 but 20 are uv=Strong, cut ~30 cds; probably accept all Strong uvcut
    # .. also have some perfect huge vector matches w/ complete prot: NUBP1 = Ecoli vec
    if($okname and not $vecNotStrong) {
      if($oldaaq =~ /partial/ and $neworflen >= 0.95 * $expect_neworflen) { $okname=0; }
      elsif($oldaaq !~ /partial/  and $ucut > 99) { $okname=0; } # this is a long-Strong uvec, named likely vector protein
    }
     
    ## ?? need levels of npct test here, dont know unless npct ~100 if uvcut affects alignment
    if($okname and $npct >= $minpct and $naln >= $minaln and $newaalen < 0.99*$naln) { ##  ??
      $fl .= ",refalignlosscut:$odiff";
      $badcut++; $note.="loss named=$namepct,$nameref,$named; "; # d=$odiff, 
    }
  }

  ## Disallow  $vectype =~ /Strong/ here?   ; skip this if refalignlosscut ?
  ## FIXME: add note for any largish change in orflen, smaller/bigger, whether uvec or trimNNN
  if(not $badcut and ($newaalen<1 or $neworflen < 0.90 * $expect_neworflen)) {
    my $odiff= $neworflen-$oldorflen;
    $fl .= ",shortorfcut:$odiff"; # note even for Strong? yes. 
    if($vecNotStrong) { $badcut++; $note.="loss orflen=$neworflen-$oldorflen; "; }
  }
  
  my $uvhdr="uvcut=$vectype,$ucut,$rub-$rue,$fl;";
  if($badcut) {  
    # complain, maybe cancel uvcut; not badcut if $vectype=~/Strong/ always?
    warn "#BAD neworf $uvhdr $id: $note\n";
    $retval= -1;
  }
  
  chomp($hd); $hd=~s,clen=,clen=$clen/,; 
  if(length($fav)>45) { $fav=substr($fav,0,20)."..".substr($fav,-20); } ## dont stick LONG fav in header ?
  $hd=~s/ / $uvhdr uvfa=$fav; /; 
  $fac =~ s/(.{60})/$1\n/g; $fac.="\n" unless($fac=~/\n$/);
  
  ## MINSIZE should change w/ named value, as per evgmrna2tsa.pl:putseq() 
  ## unique(name)/strong homol-name should keep shorter clen; but w/ annotation to support short len
  if($clen>=$MINSIZE) { 
    print OUT $hd,"\n",$fac; 
    print AAOUT ">$id $newaahdr; $uvhdr\n$newaa\n"; 
    print CDSOUT ">$id $newcdshdr; $uvhdr\n$newcds\n"; 
    $retval= ($badcut)? -1 : 1;
  } else { 
    warn "#DROP too short clen=$clen: $hd\n"; 
    $retval= ($badcut)? -1 : 0;
  }
  return $retval;
}    



sub getbestorf {
  my($id, $hdr,$cdnain,$cdnasize, $oldorfcutlen)= @_;
  
  # my $aaseq = makename($cdnaseq,".aa");
  # my $cdsseq= makename($cdnaseq,".cds");
  # my $cmd="$APPcdnabest -nostop -cdna $cdnaseq -aaseq $aaseq -cdsseq $cdsseq"; # -minaa=$MINAA 
  
  if($cdnasize > 60) {
  ## from cdna_bestorf.pl:cdna_bestorf()
  my $fullpart= "best.fwdstrand"; # bothstrand or only fwdstrand ? this is mrna, oriented
  ## ^^ need longorf vs fullorf option
  ## -- use if uvcut-long-partial > uvcut-best but < origaa (dont want > origaa)
  
  my( $bestorf)= getBestProt2( $fullpart, $cdnain);  
  # my( $orfprot, $orient)= orfParts($bestorf, [qw(protein orient)]);
  if(ref($bestorf)) {
    my $orient= undef;
    my $asfasta=1; # ($action =~ /fasta/)
    my($aalen,$pcds,$compl,$orflen,$orfprothdr,$orfprotfa)
          = proteindoc($bestorf,$cdnasize,$orient,$asfasta);
          
    if($orflen < 0.98*$oldorfcutlen and $compl =~ /complete/) {
      my($longorf)= getBestProt2( "long.fwdstrand", $cdnain);  # ask for longest partial orf
      my($aalen1,$pcds1,$compl1,$orflen1,$orfprothdr1,$orfprotfa1)
           = proteindoc($longorf,$cdnasize,$orient,$asfasta);
     if($orflen1 > $orflen) {
      $bestorf= $longorf;
      ($aalen,$pcds,$compl,$orflen,$orfprothdr,$orfprotfa)=
        ($aalen1,$pcds1,$compl1,$orflen1,$orfprothdr1,$orfprotfa1);
     }      
    }
          
    my $cdsfa= $bestorf->{sequence};  $cdsfa  =~ s/(.{60})/$1\n/g;
    (my $cdshdr=$orfprothdr) =~ s/ / type=cds; /;
    ## NOTE: >id is not in these hdr
    return($orfprothdr,$orfprotfa,$cdshdr,$cdsfa, $aalen,$pcds,$compl,$orflen);          
  }
  }
  return("","","","",0,0,0,0);          
}



sub parse_genenames  # from evgmrna2tsa.pl
{
  my($genenames)= @_;
  my($ngot,$nin)=(0,0);
  # returns in globals: (%genenames,%genenamepct,%genedbxref,%namedgenes,%cddnames) 
  %genenames=%genenamepct=%genedbxref=%namedgenes=%cddnames=();
  
  open( my $inh, $genenames) or warn("ERR: parse_genenames reading $genenames\n");
  while(<$inh>) { 
    next unless(/^\w/ and /\t/);
    chomp; $nin++;
    my($id,$name,$pctalign,$refid,$repid)=split"\t"; # may have only id, name
    my $xtra; ($name,$xtra)=split";",$name,2; 
    $name =~ s/\s+$//;
    
#     if($pctalign =~/^\d/ and $pctalign < $MIN_NAMEIDENT) { # use MIN_IDLIKE : name-like ? got some '0%' align
#       if($pctalign >= $MIN_IDLIKE) { unless($name =~ /\blike|^Uncharacterized/) {
#         $name =~ s/, putative//; 
#         unless( $name =~ s/\s+protein$/-like protein/ ) { $name .= '-like'; } ## fixme: 'xxxx protein -like'
#         }
#       } else { next; } ## should we preserve for ann.txt table ?
#     }
    
    $refid =~ s/RefID://; 
    $genedbxref{$id} .= "$refid," if($refid);
    $genedbxref{$id} .= "$repid," if($repid and not $refid =~ /^CDD:/);
    $namedgenes{$name} .= "$id,"; #? if($pctalign >= $MIN_NAMEIDENT); # for uniq name retention
    $cddnames{$id}= $name if($name =~ /CDD:/ and not $cddnames{$id});
    unless($genenames{$id} and $name =~ /CDD:/) { # or refid =~ /CDD:/
      $pctalign ||= 0; $refid ||= 0;
      $genenames{$id}= $name;  $ngot++;
      $genenamepct{$id}= $pctalign;
    }
  } close($inh);
  
  return($ngot,$nin);
}


__END__

=item named bad new orf cases

  .. how many of these named losses can be rescued with nnn+ uvgap replacements to keep same cds-frame?
  uvcut4:   pogonus only real problem, nbadcut=501 too many to check in person

  ggrep -c 'loss named=' */log.uvcut3e
banana1all3cf/log.uvcut3e:3
catfish1all4cf/log.uvcut3e:9
litova1all3f/log.uvcut3e:2
locust1evgcf/log.uvcut3e:3
pogonus1all3cf/log.uvcut3e:917   << many need checking
shrimpt1evgf/log.uvcut3e:51
whitefly1evgcf/log.uvcut3e:73    << ditto checks
zticktr2acf/log.uvcut3e:0

  tail -n1 */log.uvcut4b  : new version veccut
==> banana1all3cf/log.uvcut4b <==
#uvcut nochange=297568, ncut=110, ndrop=0, nbadcut=2 to okayset/banana1all3.uvcut4.mrna
  banana nbadcut=2 should cancel uvcut
  
==> catfish1all4cf/log.uvcut4b <==
#uvcut nochange=117531, ncut=281, ndrop=8, nbadcut=4 to okayset/catfish1all4.uvcut4.mrna
  catfish nbadcut=4, should cancel 1, ok 3
  
==> litova1all3f/log.uvcut4b <==
#uvcut nochange=125621, ncut=39, ndrop=0, nbadcut=2 to okayset/litova1all3.uvcut4.mrna
  litova nbadcut=2 Sacsin as below, 
    .. both have frameshift fix for cut;
    .. new uvcut-cds-partial is only 21 nt shorter, in both
    .. could accept uvcut w/o much refalign loss
  
==> locust1evgcf/log.uvcut4b <==
#uvcut nochange=75275, ncut=23, ndrop=1, nbadcut=3 to okayset/locust1all5asm.uvcut4.mrna
  locust nbadcut=3, one is Strong should accept (-30 cds loss for TTP repeat prot)
    other 2 moderate, -42 cds loss, need aaref checking: locust1vel5k35Loc10590t1, locust1sop4p4k23loc19454t1
    
==> pogonus1all3cf/log.uvcut4b <==
#uvcut nochange=115190, ncut=6836, ndrop=8, nbadcut=501 to okayset/pogonus1all3.uvcut4.mrna
  pogonus only real problem, nbadcut=501 too many to check in person
      237 are uvcut=Strong, should accept
      
==> shrimpt1evgf/log.uvcut4b <==
#uvcut nochange=81862, ncut=884, ndrop=7, nbadcut=26 to okayset/allstrimpt1.uvcut4.mrna
  shrimpt nbadcut=26 but 20 are uv=Strong, cut ~30 cds; probably accept all Strong uvcut
  (adjust badcut to accept Strong when also partial-cds at cut end?)
  
==> whitefly1evgcf/log.uvcut4b <==
#uvcut nochange=177356, ncut=3125, ndrop=18, nbadcut=38 to okayset/whitefly1all8.uvcut4.mrna
        26 are uvcut=Strong, should accept (including huge Ecoli vector, 2 cases NUBP1)

==> zticktr2acf/log.uvcut4b <==
#uvcut nochange=92743, ncut=22, ndrop=0, nbadcut=0 to okayset/ztick4cvel_pt3.uvcut4.mrna

#........
banana:
  n=3 named loss 
log.uvcut4b
#BAD neworf uvcut=Moderate,41,1-41,cds5end5,refalignlosscut:-33; bananavel2k25Loc52226t1: loss named=100%,368/355,346,CDD:201778,TAIR:AT3G20790.1,,NAD(P)-binding Rossmann-fold superfamily protein; 
  >> DEFINITELY cancel uvcut bananavel2k25Loc52226t1
  .. aaref Phytophthora bits=521; Phytophthora sojae aligns 1-346; complete oxidoreductase (this is a plant, what kind?)
  .. uvcut-aa bits=513; Phyto. sojae aligns 11-337, partial
  >> has longer uvcut4-aa-partial
  #getbestorf: bananavel2k25Loc52226t1 oldaa=346,91%,complete;; new.aalen=335,91%,partial5; clen=1100; strand=+; offs=2-1009;

#BAD neworf uvcut=Moderate,45,945-989,cds3end3,refalignlosscut:-24; bananavel2k25Loc1943t3: loss named=92%,227/247,224,TAIR:AT5G47120.1,,BAX inhibitor 1; 
  >> DEFINITELY cancel uvcut bananavel2k25Loc1943t3
  #getbestorf: bananavel2k25Loc1943t3 oldaa=224,68%,complete;; new.aalen=217,68%,partial3; clen=944; strand=+; offs=292-942;
  .. aaref bits=330, 1-220  align Cucumis sativus
  .. uvcut-aaref bits=326, 1-217 align
  bananavel2k25Loc1943t3  945     989     Moderate match  UniVec_Core  21/21 ident gnl|uv|NGB00379.1:1-23 Illumina 3' RNA Adapter 
  
#getbestorf: bananavel2k25Loc3425t1 oldaa=445,86%,partial5;; new.aalen=429,85%,partial5; clen=1502; strand=+; offs=2-1291;
  >>has longer uvcut4-aa-partial than before, uvcut accepted
  
  uvcut3:   
  >> DEFINITELY cancel uvcut
#getbestorf: bananavel2k25Loc52226t1 oldaa=346,91%,complete;; new=aalen=289,79%,complete; clen=1100; strand=+; offs=140-1009;
#named:  bananavel2k25Loc52226t1 named=100%,368/355,346,CDD:201778,TAIR:AT3G20790.1,,NAD(P)-binding Rossmann-fold superfamily protein;
#BAD neworf uvtype=Moderate match, bananavel2k25Loc52226t1: loss orflen=870-1041; loss named=100%,368/355,346,CDD:201778,TAIR:AT3G20790.1,,NAD(P)-binding Rossmann-fold superfamily protein; 
--
  >> Probably cancel uvcut
#getbestorf: bananavel2k25Loc1943t3 oldaa=224,68%,complete;; new=aalen=217,68%,partial3; clen=944; strand=+; offs=292-942;
#named:  bananavel2k25Loc1943t3 named=92%,227/247,224,TAIR:AT5G47120.1,,BAX inhibitor 1;
#BAD neworf uvtype=Moderate match, bananavel2k25Loc1943t3: loss named=92%,227/247,224,TAIR:AT5G47120.1,,BAX inhibitor 1; 
--
  >> Probably cancel uvcut
  -- but check align of new shorter-complete prot
#getbestorf: bananavel2k25Loc3425t1 oldaa=445,86%,partial5;; new=aalen=377,75%,complete; clen=1502; strand=+; offs=158-1291;
#named:  bananavel2k25Loc3425t1 named=94%,428/454,445,TAIR:AT5G19180.1,,E1 C-terminal related 1;
#BAD neworf uvtype=Moderate match, bananavel2k25Loc3425t1: loss orflen=1134-1338; loss named=94%,428/454,445,TAIR:AT5G19180.1,,E1 C-terminal related 1; 
--

  
catfish:  
  n=9 named loss 
  uvcut4b says 4 are bad cuts:
#BAD neworf uvcut=Moderate,39,144-182,cds5,refalignlosscut:-810; socatfishv1k95loc34193t1: loss named=100%,450/450,446,CDD:191181,DRERI:ENSDARG00000052408,UniProt:A4IG58,,Vertebrate mannosyl (Alpha-1,6-)-glycoprotein beta-1,2-N-acetylglucosaminyltransferase (MGAT2, zgc:162268); 
  >> DEFINITELY cancel uvcut socatfishv1k95loc34193t1
#BAD neworf uvcut=Moderate,50,134-1779,cds3end3,shortorfcut:-219; socatfishv1k25loc6796t5: loss orflen=1356-1575; 
  SHOULD ok this uvcut; >> 2 uvcuts, damage cds, which also has NNN spans;
  .. but problem is partly discrepant vecscreen results; www.NCBI gives diff answer to same spans
  #getbestorf: socatfishv1k25loc6796t5 oldaa=524,88%,complete;; 
       old.aalen=524,88%,complete; clen=1779; strand=+; offs=57-1631;
       new.aalen=451,78%,complete; clen=1729; strand=+; offs=245-1600;
  socatfishv1k25loc6796t5 134     164     Moderate match  UniVec_Core 31/31 ident to gnl|uv|NGB00362.1:1-61 Illumina Paired End PCR Primer 2.0
  socatfishv1k25loc6796t5 1761    1779    Moderate match  UniVec_Core : 19/19 ident to gnl|uv|NGB00363.1:1-34 Illumina Multiplexing PCR Primer 2.0 
    >> ncbi.vecs gets same spans, calls 1st Strong 31/31 ident Illumina Paired End PCR Primer 2.0
  
#BAD neworf uvcut=Moderate,60,1-60,cds5end5,refalignlosscut:-60; catfishtrin1loc147960c0t7: loss named=78%,931/1196,919,CDD:191990,DRERI:ENSDARG00000062472,UniProt:Q1LYM8,,Mouse and  SIN3 B, transcriptional regulator () (SIN3B); 
    ^^ SHOULD ok uvcut; lost 20aa from partial5' end, but not good aa
  #getbestorf: catfishtrin1loc147960c0t7 oldaa=919,100%,partial;; new.aalen=899,100%,partial; clen=2697; strand=+; offs=1-2697;
  catfishtrin1loc147960c0t7       1       60      Moderate match  UniVec_Core 53/59 ident to gnl|uv|NGB00361.1:1-92 Illumina PCR Primer 
  blastp: best hit=tilapia; uncut bits=1513; uvcut bits=1515 << take the uvcut
  
#BAD neworf uvcut=Moderate,23,555-577,cds3end3,refalignlosscut:-24; catfishvel3ik55Loc47029t1: loss named=78%,138/177,143,CDD:29029,DRERI:ENSDARG00000060130,UniProt:B3DHB2,,Dual specificity phosphatase 3 (Vaccinia virus phosphatase VH1-related); 
   ^^ SHOULD ok uvcut; lost 8aa from partial3'end
   catfishvel3ik55Loc47029t1       555     577     Moderate match  UniVec_Core 23/23 ident to gnl|uv|NGB00361.1:1-92 Illumina PCR Primer 

  uvcut3: 
  >> maybe/maybenot cancel uvcut
#getbestorf: socatfishv1k21loc455063t1 oldaa=143,62%,partial3;; new=aalen=127,99%,partial; clen=383; strand=+; offs=1-381;
#named:  socatfishv1k21loc455063t1 named=73%,149/204,143,CDD:176994,DRERI:ENSDARG00000043453,UniProt:Q6PC80,,Ribosomal protein S5;
#BAD neworf uvtype=Strong match, socatfishv1k21loc455063t1: loss named=73%,149/204,143,CDD:176994,DRERI:ENSDARG00000043453,UniProt:Q6PC80,,Ribosomal protein S5; 
--
  >> maybe/maybenot  cancel uvcut
#getbestorf: socatfishv1k45loc0t6 oldaa=304,44%,complete-utrbad;; new=aalen=197,29%,complete-utrbad; clen=2041; strand=+; offs=443-1036;
#named:  socatfishv1k45loc0t6 named=100%,266/232,304,CDD:197473,DRERI:ENSDARG00000094077,UniProt:Q1LUQ6,,Containing trypsin domains;
#BAD neworf uvtype=Strong match, socatfishv1k45loc0t6: loss orflen=594-915; loss named=100%,266/232,304,CDD:197473,DRERI:ENSDARG00000094077,UniProt:Q1LUQ6,,Containing trypsin domains; 
--
  >> DEFINITELY cancel uvcut
#getbestorf: socatfishv1k95loc34193t1 oldaa=446,84%,complete;; new=aalen=176,34%,complete-utrbad; clen=1552; strand=+; offs=919-1449;
#named:  socatfishv1k95loc34193t1 named=100%,450/450,446,CDD:191181,DRERI:ENSDARG00000052408,UniProt:A4IG58,,Vertebrate mannosyl (Alpha-1,6-)-glycoprotein beta-1,2-N-acetylglucosaminyltransferase (MGAT2, zgc:162268);
#BAD neworf uvtype=Moderate match, socatfishv1k95loc34193t1: loss orflen=531-1341; loss named=100%,450/450,446,CDD:191181,DRERI:ENSDARG00000052408,UniProt:A4IG58,,Vertebrate mannosyl (Alpha-1,6-)-glycoprotein beta-1,2-N-acetylglucosaminyltransferase (MGAT2, zgc:162268); 

--
  >> Strong uvcut, CDD: only name loss, check this one
#getbestorf: socatfishv1k21loc431915t1 oldaa=118,99%,partial;; new=aalen=39,45%,complete-utrbad; clen=262; strand=+; offs=14-133;
#named:  socatfishv1k21loc431915t1 named=98%,65/66,118,CDD:201275,,CDD: HisKA, His Kinase A (phospho-acceptor) domain;
#BAD neworf uvtype=Strong match, socatfishv1k21loc431915t1: loss orflen=120-354; loss named=98%,65/66,118,CDD:201275,,CDD: HisKA, His Kinase A (phospho-acceptor) domain; 
--
  >> probably ok uvcut, lost only 20aa/919 from partial prot
#getbestorf: catfishtrin1loc147960c0t7 oldaa=919,100%,partial;; new=aalen=899,100%,partial; clen=2697; strand=+; offs=1-2697;
#named:  catfishtrin1loc147960c0t7 named=78%,931/1196,919,CDD:191990,DRERI:ENSDARG00000062472,UniProt:Q1LYM8,,Mouse and  SIN3 B, transcriptional regulator () (SIN3B);
#BAD neworf uvtype=Moderate match, catfishtrin1loc147960c0t7: loss named=78%,931/1196,919,CDD:191990,DRERI:ENSDARG00000062472,UniProt:Q1LYM8,,Mouse and  SIN3 B, transcriptional regulator () (SIN3B); 
--
  >>  ok Strong uvcut, partial before, 82% align, lost only 13/214 aa.
#getbestorf: catfishvel3ik45Loc1006t3 oldaa=214,99%,partial3;; new=aalen=201,99%,partial3; clen=607; strand=+; offs=3-605;
#named:  catfishvel3ik45Loc1006t3 named=82%,204/249,214,CDD:179811,DRERI:ENSDARG00000019778,UniProt:Q6DHL6,,40S ribosomal protein S6;
#BAD neworf uvtype=Strong match, catfishvel3ik45Loc1006t3: loss named=82%,204/249,214,CDD:179811,DRERI:ENSDARG00000019778,UniProt:Q6DHL6,,40S ribosomal protein S6; 
--
  >> maybe cancel uvcut
  .. new aalen < old align, has ref-species (zfish) align; is this cancel uvcut or not?
  .. add $vectype =~ Strong
#getbestorf: catfishvel3ik55Loc47029t1 oldaa=143,74%,partial3;; new=aalen=135,73%,partial3; clen=554; strand=+; offs=149-553;
#named:  catfishvel3ik55Loc47029t1 named=78%,138/177,143,CDD:29029,DRERI:ENSDARG00000060130,UniProt:B3DHB2,,Dual specificity phosphatase 3 (Vaccinia virus phosphatase VH1-related);
#BAD neworf uvtype=Moderate match, catfishvel3ik55Loc47029t1: loss named=78%,138/177,143,CDD:29029,DRERI:ENSDARG00000060130,UniProt:B3DHB2,,Dual specificity phosphatase 3 (Vaccinia virus phosphatase VH1-related); 


Locust:
  n=3 named loss  
  >> probably ok these 3 : small aa loss to partial aa;
  .. check locust1vel5k35Loc10590t1 loss named=98%, Moderate uvtype.
#getbestorf: locust1Svel1K31L7072t3 oldaa=191,68%,partial5;; new=aalen=181,67%,partial5; clen=813; strand=+; offs=2-547;
#named:  locust1Svel1K31L7072t3 named=99%,185/187,191,CDD:29151,UniRef50_F4WTM8,UniProt:F4WTM8_ACREC,,Tetratricopeptide repeat protein 36-like protein;
#BAD neworf uvtype=Strong match, locust1Svel1K31L7072t3: loss named=99%,185/187,191,CDD:29151,UniRef50_F4WTM8,UniProt:F4WTM8_ACREC,,Tetratricopeptide repeat protein 36-like protein; 
--
  >> maybe cancel uvcut .. check new align vs old
#getbestorf: locust1vel5k35Loc10590t1 oldaa=550,99%,partial;; new=aalen=536,99%,partial; clen=1611; strand=+; offs=2-1609;
#named:  locust1vel5k35Loc10590t1 named=98%,544/557,550,CDD:31971,UniRef90_J9K9Z4,UniProt:J9K9Z4_ACYPI,,Alkaline phosphatase;
#BAD neworf uvtype=Moderate match, locust1vel5k35Loc10590t1: loss named=98%,544/557,550,CDD:31971,UniRef90_J9K9Z4,UniProt:J9K9Z4_ACYPI,,Alkaline phosphatase; 
--
#getbestorf: locust1sop4p4k23loc19454t1 oldaa=175,99%,partial3;; new=aalen=161,99%,partial; clen=484; strand=+; offs=2-484;
#named:  locust1sop4p4k23loc19454t1 named=81%,172/213,175,UniRef90_Q24040,UniProt:BGB_DROME,,Big brother;
#BAD neworf uvtype=Moderate match, locust1sop4p4k23loc19454t1: loss named=81%,172/213,175,UniRef90_Q24040,UniProt:BGB_DROME,,Big brother; 

litova/whiteshrimp
  only n=2 named loss, same gene
  >> cancel uvcut, moderate only; major loss of ref-prot align
#getbestorf: sowshrimp2v2k21loc74886t1 oldaa=4266,92%,complete;; new=aalen=2391,51%,complete-utrpoor; clen=13815; strand=+; offs=5749-12924;
#named:  sowshrimp2v2k21loc74886t1 named=100%,4678/4579,4266,CDD:128987,UniRef50_Q9NZJ4,UniProt:SACS_HUMAN,,Sacsin;
#BAD neworf uvtype=Moderate match, sowshrimp2v2k21loc74886t1: loss orflen=7176-12801; loss named=100%,4678/4579,4266,CDD:128987,UniRef50_Q9NZJ4,UniProt:SACS_HUMAN,,Sacsin; 
--
#getbestorf: litovavel3k21Loc3747t3 oldaa=4238,92%,complete;; new=aalen=2363,51%,complete-utrpoor; clen=13733; strand=+; offs=5752-12843;
#named:  litovavel3k21Loc3747t3 named=100%,4678/4579,4238,CDD:128987,UniRef50_Q9NZJ4,UniProt:SACS_HUMAN,,Sacsin;
#BAD neworf uvtype=Moderate match, litovavel3k21Loc3747t3: loss orflen=7092-12717; loss named=100%,4678/4579,4238,CDD:128987,UniRef50_Q9NZJ4,UniProt:SACS_HUMAN,,Sacsin; 
  .. but vector.span is smallish, has bad aa orf effects: frameshift* 
sowshrimp2v2k21loc74886t1       5712    5737    Moderate match  UniVec_Core
litovavel3k21Loc3747t3  5715    5740    Moderate match  UniVec_Core
  
shrimpt1evgf/Tiger shrimp; 
  n=55 named loss  
   -- includes these Riboprots, all Strong uvec
  but n=12 new are longer aa, some however going complete=>partial
    -- check shrimt1vel1k41Loc2411t22 oldaa=422,58%,complete-utrpoor;; new=aalen=537,76%,partial5;
        >> 73%,423/582,422 Beta-N-acetylglucosaminidase NAG3 
    -- check shrimt1vel1k41Loc364t4 oldaa=644,86%,partial3;; new=aalen=720,99%,partial; 
        >> 76%,644/847,644 Glycogen phosphorylase, liver form; CDD: GT1_Glycogen_Phosphorylase, 
        
#BAD neworf uvtype=Strong match, shrimt1vel1k45Loc1249t2: loss named=100%,184/184,195,CDD:30440,UniRef50_P18621,UniProt:RL17_HUMAN,,60S ribosomal protein L17; 
#BAD neworf uvtype=Strong match, shrimt1vel1k55Loc236t1: loss named=93%,165/178,164,CDD:30534,UniRef50_J6EVU9,UniProt:J6EVU9_TRIAS,,40s ribosomal protein s15; 
#BAD neworf uvtype=Strong match, shrimt1vel1k41Loc2t4777: loss named=83%,124/150,135,CDD:189924,UniRef50_D3BBQ4,UniProt:D3BBQ4_POLPA,,40S ribosomal protein S26; 
#BAD neworf uvtype=Strong match, shrimt1vel1k45Loc211t1: loss named=93%,165/178,160,CDD:30534,UniRef50_J6EVU9,UniProt:J6EVU9_TRIAS,,40s ribosomal protein s15; 
 
-- 
#getbestorf: shrimt1vel1k41Loc145t1 oldaa=371,49%,complete-utrpoor;; new=aalen=367,49%,complete-utrpoor; clen=2241; strand=+; offs=115-1218;
#named:  shrimt1vel1k41Loc145t1 named=100%,405/399,371,CDD:29152,UniRef50_B7Z0T0,UniProt:B7Z0T0_DROME,,Melanization protease 1;
#BAD neworf uvtype=Moderate match, shrimt1vel1k41Loc145t1: loss named=100%,405/399,371,CDD:29152,UniRef50_B7Z0T0,UniProt:B7Z0T0_DROME,,Melanization protease 1; 
--
  >> check this new,old align
#getbestorf: shrimt1vel1k45Loc1249t2 oldaa=195,91%,partial3;; new=aalen=182,90%,partial3; clen=603; strand=+; offs=56-601;
#named:  shrimt1vel1k45Loc1249t2 named=100%,184/184,195,CDD:30440,UniRef50_P18621,UniProt:RL17_HUMAN,,60S ribosomal protein L17;
#BAD neworf uvtype=Strong match, shrimt1vel1k45Loc1249t2: loss named=100%,184/184,195,CDD:30440,UniRef50_P18621,UniProt:RL17_HUMAN,,60S ribosomal protein L17; 
--
#getbestorf: shrimt1vel1k55Loc236t1 oldaa=164,88%,partial5;; new=aalen=148,87%,partial5; clen=508; strand=+; offs=2-448;
#named:  shrimt1vel1k55Loc236t1 named=93%,165/178,164,CDD:30534,UniRef50_J6EVU9,UniProt:J6EVU9_TRIAS,,40s ribosomal protein s15;
#BAD neworf uvtype=Strong match, shrimt1vel1k55Loc236t1: loss named=93%,165/178,164,CDD:30534,UniRef50_J6EVU9,UniProt:J6EVU9_TRIAS,,40s ribosomal protein s15; 
--
  >> check this one
#getbestorf: shrimt1vel1k41Loc4606t2 oldaa=671,94%,partial3;; new=aalen=660,94%,partial3; clen=2093; strand=+; offs=112-2091;
#named:  shrimt1vel1k41Loc4606t2 named=100%,713/684,671,CDD:201371,UniRef50_Q7PWB1,UniProt:RETM_ANOGA,,Real-time;
#BAD neworf uvtype=Strong match, shrimt1vel1k41Loc4606t2: loss named=100%,713/684,671,CDD:201371,UniRef50_Q7PWB1,UniProt:RETM_ANOGA,,Real-time; 
--  
  >> check this one
#getbestorf: shrimt1vel1k45Loc175t2 oldaa=411,75%,partial5;; new=aalen=397,75%,partial5; clen=1585; strand=+; offs=2-1195;
#named:  shrimt1vel1k45Loc175t2 named=100%,406/403,411,CDD:200598,UniRef50_A8W492,UniProt:A8W492_TRICA,,Chitin deacetylase 6;
#BAD neworf uvtype=Strong match, shrimt1vel1k45Loc175t2: loss named=100%,406/403,411,CDD:200598,UniRef50_A8W492,UniProt:A8W492_TRICA,,Chitin deacetylase 6;   
--


Whitefly
  n=73 named loss

--
  >> UVector! really big uvcut, no new prot? is this Strong uv or valid gene?  
vecs: gnl|uv|FJ160466.1:1899-7260-49 BAC cloning vector pEZ BAC   Identities:1736/1736(100%)  
** Ecoli protein is perfect 100% aa match
#getbestorf: sowhiteflyv8k61loc90865t1 oldaa=391,67%,complete;; new=
#named:  sowhiteflyv8k61loc90865t1 named=92%,293/320,391,CDD:31385,UniRef50_P53384,UniProt:NUBP1_HUMAN,,Cytoso
lic Fe-S cluster assembly factor NUBP1;
#BAD neworf uvtype=Strong match, sowhiteflyv8k61loc90865t1: loss orflen=0-1176; loss named=92%,293/320,391,CDD
:31385,UniRef50_P53384,UniProt:NUBP1_HUMAN,,Cytosolic Fe-S cluster assembly factor NUBP1; 
--
  >> UVector! ditto; really big uvcut, no new prot? is this Strong uv or valid gene?
#getbestorf: whiteflyvel7ik65Loc169t3 oldaa=322,88%,partial3;; new=
#named:  whiteflyvel7ik65Loc169t3 named=92%,155/168,322,CDD:213813,UniRef50_I3L531,UniProt:I3L531_HUMAN,,Cytos
olic Fe-S cluster assembly factor NUBP1;
#BAD neworf uvtype=Strong match, whiteflyvel7ik65Loc169t3: loss orflen=0-966; loss named=92%,155/168,322,CDD:2
13813,UniRef50_I3L531,UniProt:I3L531_HUMAN,,Cytosolic Fe-S cluster assembly factor NUBP1; 
--
  >> probably cancel big uvcut Moderate match, check
#getbestorf: whitefly1vel5k35Loc7929t1 oldaa=323,48%,complete-utrpoor;; new=aalen=92,13%,complete-utrbad; clen
=1994; strand=+; offs=861-1139;
#named:  whitefly1vel5k35Loc7929t1 named=67%,280/416,323,CDD:144775,UniRef50_Q9UQV4,UniProt:LAMP3_HUMAN,,Lysos
ome-associated membrane glycoprotein 3;
#BAD neworf uvtype=Moderate match, whitefly1vel5k35Loc7929t1: loss orflen=279-972; loss named=67%,280/416,323,
CDD:144775,UniRef50_Q9UQV4,UniProt:LAMP3_HUMAN,,Lysosome-associated membrane glycoprotein 3; 
--
  >> probably accept uvcut; Known bug gene Wnt but aligns to smaller prot
  vecs: 24/24 ident to gnl|uv|NGB00325.1:1-39 Adaptor used in I.M.A.G.E. library NIH_MGC_126 and other libraries 
  -- no ref-aalign to 5'uvcut
#getbestorf: whitefly1vel5k45Loc107t9 oldaa=406,96%,partial5;; new=aalen=391,96%,partial5; clen=1220; strand=+; offs=1-1176;
#named:  whitefly1vel5k45Loc107t9 named=100%,429/386,406,CDD:201009,UniRef50_Q17NQ1,UniProt:Q17NQ1_AEDAE,,Wnt;
#BAD neworf uvtype=Strong match, whitefly1vel5k45Loc107t9: loss named=100%,429/386,406,CDD:201009,UniRef50_Q17
NQ1,UniProt:Q17NQ1_AEDAE,,Wnt; 
--
  >> NO keep uvcut: maybe cancel uvcut, check
  vecs: 30/30 ident end3: gnl|uv|NGB00362.1:1-61 Illumina Paired End PCR Primer 2.0 
  -- no ref-aalign to 3'uvcut; longest align mouse 1..132, ref size 132aa
#getbestorf: whitefly1vel5k35Loc3370t2 oldaa=143,82%,partial3;; new=aalen=133,80%,partial3; clen=493; strand=+
; offs=93-491;
#named:  whitefly1vel5k35Loc3370t2 named=96%,137/142,143,CDD:192052,UniRef50_E0VP34,UniProt:E0VP34_PEDHC,,Phos
phatidylinositol N-acetylglucosaminyltransferase subunit P, putative;
#BAD neworf uvtype=Strong match, whitefly1vel5k35Loc3370t2: loss named=96%,137/142,143,CDD:192052,UniRef50_E0V
P34,UniProt:E0VP34_PEDHC,,Phosphatidylinositol N-acetylglucosaminyltransferase subunit P, putative; 
--
  >> maybe cancel uvcut, check
  whitefly1vel5k35Loc9242t3 also
#getbestorf: whitefly1vel5k35Loc9242t2 oldaa=533,88%,partial3;; new=aalen=524,87%,partial3; clen=1789; strand=
+; offs=218-1789;
#named:  whitefly1vel5k35Loc9242t2 named=100%,530/530,533,CDD:201077,UniRef50_O75795,UniProt:UDB17_HUMAN,,UDP-
glucuronosyltransferase 2B17;
#BAD neworf uvtype=Strong match, whitefly1vel5k35Loc9242t2: loss named=100%,530/530,533,CDD:201077,UniRef50_O7
5795,UniProt:UDB17_HUMAN,,UDP-glucuronosyltransferase 2B17; 
--
  >> DEFINITELY cancel uvcut, check; Selenocysteine problem? no
  -- refaligns are full 1-440, agambie, peaaphid, same size:446,440aa
  -- refalings to uvcut-aa is -100 bits of full-aa, complete aa-aligns refute uvector
  vecs:  26/26 ident inside to gnl|uv|U29899.1:1847-4463 Cloning vector pACT2 MatchmakerII 
  whitefly1vel5k55Loc4281t3 also
#getbestorf: whitefly1vel5k55Loc4281t2 oldaa=440,64%,complete;; new=aalen=423,63%,complete; clen=2013; strand=
+; offs=154-1425;
#named:  whitefly1vel5k55Loc4281t2 named=98%,437/445,440,CDD:184450,UniRef50_D3PJW2,UniProt:D3PJW2_9MAXI,,Sele
nocysteine lyase;
#BAD neworf uvtype=Moderate match, whitefly1vel5k55Loc4281t2: loss named=98%,437/445,440,CDD:184450,UniRef50_D
3PJW2,UniProt:D3PJW2_9MAXI,,Selenocysteine lyase; 
--
  >> maybe cancel uvcut, check 
#getbestorf: whitefly1vel6k65Loc12474t1 oldaa=533,91%,partial3;; new=aalen=524,91%,partial3; clen=1723; strand
=+; offs=152-1723;
#named:  whitefly1vel6k65Loc12474t1 named=100%,530/530,533,CDD:201077,UniRef50_O75795,UniProt:UDB17_HUMAN,,UDP
-glucuronosyltransferase 2B17;
#BAD neworf uvtype=Strong match, whitefly1vel6k65Loc12474t1: loss named=100%,530/530,533,CDD:201077,UniRef50_O
75795,UniProt:UDB17_HUMAN,,UDP-glucuronosyltransferase 2B17; 
--
  >> maybe cancel uvcut, check 
  vecs: 45/46 ident to end3: gnl|uv|NGB00362.1:1-61 Illumina Paired End PCR Primer 2.0
  -- longest refalign 1-140, reflens 148 aa
#getbestorf: whiteflyvel7ik55Loc2495t1 oldaa=155,80%,partial3;; new=aalen=140,78%,partial3; clen=532; strand=+
; offs=112-531;
#named:  whiteflyvel7ik55Loc2495t1 named=99%,147/148,155,CDD:30549,UniRef50_P48160,UniProt:RL27A_DICDI,,60S ri
bosomal protein L27a;
#BAD neworf uvtype=Strong match, whiteflyvel7ik55Loc2495t1: loss named=99%,147/148,155,CDD:30549,UniRef50_P481
60,UniProt:RL27A_DICDI,,60S ribosomal protein L27a; 
--


Pogonus beetle
  n=917 named loss; 
  
  n=112 longer newaa
    pogonusvel1pk45Loc3668t4 oldaa=338,56%,complete-utrpoor;; new=aalen=458,78%,complete;
      >> 100%,323/323,338   8-oxoguanine DNA glycosylase
    pogonusvel1pk45Loc2830t4 oldaa=161,65%,complete;; new=aalen=221,99%,partial;
    pogonusvel2nk35Loc1321t71 oldaa=207,52%,complete-utrpoor;; new=aalen=235,60%,complete;
      >> 83%,204/245,207 Peptidyl-prolyl cis-trans isomerase

  named loss .. more Riboprots
    .. many are shorted by 10-20 aa for ~100% aligns to refs, often uvtype=Strong, 
    .. could be either uvcut or not; need to check how many?
    .. do we need ref-prot aligns for uvcut span?
    
#getbestorf: sobeetlepogo1pk39loc29349t1 oldaa=254,79%,complete;; new=aalen=136,44%,complete-utrbad; clen=925;
 strand=+; offs=29-439;
#named:  sobeetlepogo1pk39loc29349t1 named=100%,245/204,254,CDD:179717,UniRef50_P46782,UniProt:RS5_HUMAN,,40S 
ribosomal protein S5;
#BAD neworf uvtype=Strong match, sobeetlepogo1pk39loc29349t1: loss orflen=411-765; loss named=100%,245/204,254
,CDD:179717,UniRef50_P46782,UniProt:RS5_HUMAN,,40S ribosomal protein S5; 

#getbestorf: sobeetlepogo1pk45loc72059t1 oldaa=265,93%,complete;; new=aalen=254,93%,partial3; clen=816; strand
=+; offs=54-815;
#named:  sobeetlepogo1pk45loc72059t1 named=99%,259/262,265,CDD:211392,UniRef50_Q6C6D1,UniProt:NSA2_YARLI,,Ribo
some biogenesis protein NSA2;
#BAD neworf uvtype=Strong match, sobeetlepogo1pk45loc72059t1: loss named=99%,259/262,265,CDD:211392,UniRef50_Q
6C6D1,UniProt:NSA2_YARLI,,Ribosome biogenesis protein NSA2; 

  vecs: 37/37 ident end5 to gnl|uv|NGB00361.1:1-92 Illumina PCR Primer 
  refaln: tribol: 8-424aa, bits=581; Choristoneura: 3-423aa, bits=476
  refaln uvcut-aa: bits:543 tribol; 493 next; (351aa)
  refaln uvcut-aa-partial: bits:581 tribol;  (412aa, align 11-411)  Tetropium: align 9-412
  **  uvcut-aa-partial is best answer
  -- conflict here, uvcut chops about 13 aa at front, with some align to ref-genes but not perfect aa align  
  uvcut-aa definitely is missing leading ~80aa; try uvcut-aa-partial5
  >pogonusvel1pk45Loc1196t5 aalen=412,79%,partial5; clen=1557; strand=+; offs=1-1239; uvcut=1-37,37,cds5end5; uvfa=ATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT; type=cdna; aa c  strand=+; offs=2-1276;
  VKHQLTDDFIRPISKTTSTANSDSITSRQQDTNDDRSNKKLFRLRKYHRSLSTVPSRAVP

#getbestorf: pogonusvel1pk45Loc1196t5 oldaa=424,79%,partial5;; new=aalen=351,67%,complete; clen=1557; strand=+
; offs=184-1239;
#named:  pogonusvel1pk45Loc1196t5 named=97%,427/438,424,CDD:173833,UniRef50_Q56CY8,UniProt:Q56CY8_ANTGR,,Farne
syl diphosphate synthase;
#BAD neworf uvtype=Strong match, pogonusvel1pk45Loc1196t5: loss orflen=1056-1275; loss named=97%,427/438,424,C
DD:173833,UniRef50_Q56CY8,UniProt:Q56CY8_ANTGR,,Farnesyl diphosphate synthase; 

#getbestorf: pogonusvel2nk33Loc8140t1 oldaa=208,99%,partial;; new=aalen=194,99%,partial; clen=584; strand=+; o
ffs=2-583;
#named:  pogonusvel2nk33Loc8140t1 named=97%,196/202,208,CDD:148815,UniRef50_Q86BF9,UniProt:Q86BF9_DROME,,Odora
nt-binding protein 59a;
#BAD neworf uvtype=Moderate match, pogonusvel2nk33Loc8140t1: loss named=97%,196/202,208,CDD:148815,UniRef50_Q8
6BF9,UniProt:Q86BF9_DROME,,Odorant-binding protein 59a; 


#getbestorf: pogonusvel3xk63Loc6931t1 oldaa=504,92%,partial3;; new=aalen=493,92%,partial3; clen=1603; strand=+
; offs=125-1603;
#named:  pogonusvel3xk63Loc6931t1 named=100%,544/541,504,CDD:177721,UniRef50_Q960X4,UniProt:TIP60_DROME,,Histo
ne acetyltransferase Tip60;
#BAD neworf uvtype=Strong match, pogonusvel3xk63Loc6931t1: loss named=100%,544/541,504,CDD:177721,UniRef50_Q96
0X4,UniProt:TIP60_DROME,,Histone acetyltransferase Tip60; 


      
=cut
  
