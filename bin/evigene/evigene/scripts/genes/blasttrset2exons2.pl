#!/usr/bin/env perl
# evigene/scripts/genes/blasttrset2exons.pl

=item about

  infer exons from self-blast-align of alternate transcript set
  input is high-ident basic-local (exon) aligned transcripts as from 
    evigene tr2aacds tmpfiles/trasm-self98.blastn
  use transcript pair align spans from makeblastscore.pl -showspan
  record common align break-points of alt-trans gene set (inferred exon splice sites)
  output table of such, as 'ic:splice-site-nums,..' and exon splice points, suited to evigene overinchain2locus.pl
  similar to, replacement of  evigene/scripts/genes/inexchains_gcb.pl
  
=item usage

  a. prepare input table with exon align spans from trasm-self.blastn
  
  $evigene/scripts/makeblastscore3.pl -pIDENTMIN 99.999 -pmin 0.01 -showspan=2 -tall \
    -sizes ref_arath16apnrcd1x.cds.qual  ref_arath16apnrcd1x-self98.blastn.gz > ref_arath16apnrcd1x-self100.btall
 
  b. exon table using publicset/pubid info,  input aligned exon table sorted by biggest ref (col7=reflen, col2=refid)
  
  sort -k7,7nr -k2,2 -k6,6nr -k1,1 ref_arath16apnrcd1x-self100.btall | \
  env pubids=ref_arath16ap.pubids pfrag=0.80 minunalign=20 $evigene/scripts/genes/blasttrset2exons2.pl \
    > ref_arath16apnrcd1x.pubfrag80u20exon2a
  
  c. version exons2 includes exon chains from relevant part of evigene/scripts/genes/overinchain2locus.pl

  d. output exon table 
    pubid, trasm_oid, xn:exon_chain, ichain:locus:symbolic,  align-spans
  Atref9Evm000001t1	AT1G67120.2	xn:g1t1x1,g1t1x9589,g1t1x9599,g1t1x9610,g1t1x16203	ichain:l2551.t1:ACBCA	1-16203;1-9589,9610-16203;
  Atref9Evm000001t2	AT1G67120.1	xn:g1t1x1,g1t1x9589,g1t1x9610,g1t1x16203	icalt:l2551.t2:ACCA	1-9589,9590-16182;1-16182;
    
=cut

use strict;

my $debug = $ENV{debug}||0; 
my $pFRAG= $ENV{pfrag}||0.50; #was 0.30; # fragments test, need to check ref set

my $SKIPSELF = $ENV{skipself}||0; # skip self, any value to self align? YES, dont use skipself

=item SKIPSELF

  skip self-align (td eq rd), any value to self align? tsp == rsp == 1-length
  wc ref_arath16apnrcd1x.*u20exontab
   22556   90224 2529300 ref_arath16apnrcd1x.ssu20exontab # all singletons dropped, eg AT1G01010.1 
   40499  161996 4114257 ref_arath16apnrcd1x.u20exontab

  need to check cases where 1 alt has extended start or stop, ie 5' alt start forms, eg. perfect-frag subset
  perffrags: 
    AT1G01110.1/AT1G01110.2 << kept past cd1
    dropped in cdhit1: AT1G01580.1/AT1G01580.2 AT1G04240.1/AT1G04240.2 
      AT1G01420.1/AT1G01420.2  AT1G01390.1/AT1G01390.2 .. dropped in cdhit1

  ** Need to keep self align span to keep 5' alt starts, perfect frags
noSKIPSELF ref_arath16apnrcd1x.u20exontab
  AT1G01110.1	ic:i49,i50	xn:g6848t2x490,g6848t2x1584	1-1095;1-1095;
  AT1G01110.2	ic:i51,i49,i50	xn:g6848t2x1,g6848t2x490,g6848t2x1584	1-1584;490-1584;
SKIPSELF ref_arath16apnrcd1x.ssu20exontab
  AT1G01110.1	ic:i37,i38	xn:g3147t2x490,g3147t2x1584	1-1095;   << same splice patt w/o self end points
  AT1G01110.2	ic:i37,i38	xn:g3147t2x490,g3147t2x1584	490-1584; <<

=cut

my $BREAKPT_SLOP= 3; # bp offset to match other align breakpoints
my $MINUNALIGNBP=$ENV{minunalign}||20; #was 30; 

=item MINUNALIGNBP
  maybe want smaller, 20? or larger, 60? test
  dup ichains ref_arath16apnrcd1x.u20exontab vs ref_arath16apnrcd1x.u60exontab
  ref_arath16apnrcd1x.u60exontab =  1648
  ref_arath16apnrcd1x.u20exontab  = 1086

==> ref_arath16apnrcd1x.u20exontab <==
AT1G01010.1	ic:i01,i02	xn:g10404t1x1,g10404t1x1290	1-1290; << uniq singleton gene, self span, want that?
AT1G01020.1	ic:i03,i04,i05,i06,i07,i08,i09,i10	xn:g18736t1x1,g18736t1x142,g18736t1x144,g18736t1x234,g18736t1x293,g18736t1x353,g18736t1x512,g18736t1x738	1-738;144-738;142-738;1-512;142-234,353-512;
AT1G01020.2	ic:i03,i04,i05,i06,i11,i08,i09,i12	xn:g18736t1x1,g18736t1x142,g18736t1x144,g18736t1x234,g18736t2x293,g18736t1x353,g18736t1x512,g18736t2x576	1-512;144-512;142-512;1-576;142-234,353-576;
AT1G01020.3	ic:i13,i05,i06,i14,i08,i09,i10	xn:g18736t3x1,g18736t1x144,g18736t1x234,g18736t3x266,g18736t1x353,g18736t1x512,g18736t1x738	117-711;1-711;117-711;117-485;117-207,326-485;
AT1G01020.5	ic:i04,i05,i06,i15,i08,i09,i10	xn:g18736t1x142,g18736t1x144,g18736t1x234,g18736t5x152,g18736t1x353,g18736t1x512,g18736t1x738	1-597;3-597;1-597;1-371;1-93,212-371;
AT1G01020.6	ic:i04,i05,i06,i08,i09,i12	xn:g18736t1x142,g18736t1x144,g18736t1x234,g18736t1x353,g18736t1x512,g18736t2x576	1-91,92-251;3-91,92-251;1-91,92-251;1-91,92-315;1-315;

==> ref_arath16apnrcd1x.u60exontab <==
AT1G01010.1	ic:i01,i02	xn:g10404t1x1,g10404t1x1290	1-1290;
AT1G01020.1	ic:i03,i04,i05,i06,i07,i08,i09,i10	xn:g18736t1x1,g18736t1x142,g18736t1x144,g18736t1x234,g18736t1x293,g18736t1x353,g18736t1x512,g18736t1x738	1-738;144-738;142-738;1-512;142-234,353-512;
AT1G01020.2	ic:i03,i04,i05,i06,i11,i08,i09,i12	xn:g18736t1x1,g18736t1x142,g18736t1x144,g18736t1x234,g18736t2x293,g18736t1x353,g18736t1x512,g18736t2x576	1-512;144-512;142-512;1-576;142-234,353-576;
AT1G01020.3	ic:i13,i05,i06,i14,i08,i09,i10	xn:g18736t3x1,g18736t1x144,g18736t1x234,g18736t3x266,g18736t1x353,g18736t1x512,g18736t1x738	117-711;1-711;117-711;117-485;117-207,326-485;
AT1G01020.5	ic:i04,i05,i06,i15,i08,i09,i10	xn:g18736t1x142,g18736t1x144,g18736t1x234,g18736t5x152,g18736t1x353,g18736t1x512,g18736t1x738	1-597;3-597;1-597;1-371;1-93,212-371;
AT1G01020.6	ic:i04,i05,i06,i08,i09,i12	xn:g18736t1x142,g18736t1x144,g18736t1x234,g18736t1x353,g18736t1x512,g18736t2x576	1-91,92-251;3-91,92-251;1-91,92-251;1-91,92-315;1-315;

=cut

my $INSORTED= $ENV{sorted}||0; # later
my $pubidtab= $ENV{pubids}||"";

my $IS_CDSALIGN= $ENV{utralign} ? 0 : 1; 
my $ADD_UNIQX  = $ENV{nouniq} ? 0 : 1; # insertUniqExons() bugs? needs check
#o my $NOTALLCOMMSETS= $ENV{allcommsets}?0:1; # commonexons, default ON??
my $NOTALLCOMMSETS= $ENV{notallcommsets}?1:0; # commonexons, default ON??
#o my $BYGENE    = $ENV{bygene} ? 1 : 0; # only when have pubids, process w/ common exon info per locus
my $BYGENE    = ($pubidtab)?1:0; $BYGENE=0 if($ENV{notbygene}); # only when have pubids, process w/ common exon info per locus

use constant { MAYBEFRAG => 1, ISFRAG => 3 };
my(%GI, $GI, %gval, %xon, %gxon, %gxord, %gxun, %gxalign, %isfrag, %fragof);
my($nin,$nok,$nfrag,$skipdiffloc,$skipnotloc,$skipdidfrag,$ndupout)= (0) x 9;

sub MAIN_BEGIN {}

#* FIXME: add sort input here unless -sorted flag
# sort -k7,7nr -k2,2 -k6,6nr -k1,1 ref_arath16apnrcd1x-self98.btall
# hassle w/ input ARGV or STDIN?
unless($INSORTED) {
  # open(IN,"sort -k7,7nr -k2,2 -k6,6nr -k1,1 $infile |") or die "sorting input $infile";
  # $inh= *IN;
}


my($npubtr,$psizeh,$gscoreh,$locush,$podh,$pubqualh,$aaqualh,$altscoreh)
   = readPubidTab($pubidtab) if($pubidtab and -f $pubidtab);


while(<>) {
  next if(/^Query|^\W/); 
  my @v=split; 
  # SORTED by longest rd, so frag td are at end of rd list
  my($td,$rd,$bs,$ida,$aln,$tw,$rw,$tsp,$rsp)=@v; # makeblastscore.pl table with align spans
  my $isself=($td eq $rd);
  next if($SKIPSELF and $isself); # skip self, any value to self align? tsp == rsp == 1-length
  $nin++;
    
  # havepubids/locus will keep only those of blastself in okayset .. want that? ie ignore already dropped alts
  # ALSO replace td,rd with pubids so idparts() has locus ids
  if($npubtr) {
    unless($locush->{$td} and $locush->{$rd}) { $skipnotloc++; next; }
    if($locush->{$td} ne $locush->{$rd}) { $skipdiffloc++; next; }
  }
  
  if($isfrag{$td} or $isfrag{$rd}) { $skipdidfrag++; next; } # did already, to longest rd, no dup calls
  my $tdfrag= ($tw < $pFRAG * $rw)?MAYBEFRAG:0; # also check tspan for 1 or 2 aligns covering most of td to larger rd

  my($tgi,$tg,$ti)= idparts($td);
  my($rgi,$rg,$ri)= idparts($rd);
    
  my($ro)= ($rsp=~s/:(.)$//)?$1:0;  
  $tsp =~ s,^.*/,,; 
  if($isself) { my $oldv=$gval{$td}||""; $gval{$td}= "$tsp;".$oldv; }
  else { $gval{$td}.="$tsp;";  }
  my @tsp=split",",$tsp; my @rsp=split",",$rsp; my $lte=0;

  use constant FRAG_EXONMAX => 2; # or 2 or 3?
  if($tdfrag and @tsp <= FRAG_EXONMAX) { #was ( @tsp < 3) 
    # maybe @tsp <= 1, as @tsp <= 2 includes valid skip-exon alts : measure @rsp separation?
    my($tb,$te1)= split"-",$tsp[0];  
    my($tb2,$te)= split"-",$tsp[-1]; 
    my $fok= ($tb < $MINUNALIGNBP and ($tw - $te) < $MINUNALIGNBP)?1:0; 
    $fok=0 if($fok and @tsp>1 and ($tb2 - $te1) >= $MINUNALIGNBP);
    
    if($fok) {
      my($rb,$re1)= split"-",$rsp[0];  
      my($rb2,$re)= split"-",$rsp[-1];     
      #?? special case for alt 5'starts, tb<min if cds-align only ?? 
      #?? long diff (ref-end - te-end) would be frag rather than alt start
      # IS_CDSALIGN: allow tb < min, rb >> 1 as altstart, but frag if te << re
      # BUT there can be alt3stops, as with intron-retained or skipped near 3'end, within min unalign 
      # better to use utr-extended align
      $fok=0 if( $IS_CDSALIGN and ($rw - $re) < $MINUNALIGNBP);
      $fok=0 if( @rsp>1 and ($rb2 - $re1) >= $MINUNALIGNBP); # measure of ref exon skipped ?
    }
    if($fok) { $tdfrag=ISFRAG; $nfrag++; }
    else { } #not-frag?    
  }
  $nok++ unless($tdfrag == ISFRAG);
  
  for my $i (0..$#rsp) { 
    my($rb,$re)= split"-",$rsp[$i]; 
    my($tb,$te)= split"-",$tsp[$i];  
    my @be=($rb,$tb,$re,$te);  my $validxi=0;
    while( my @rt=splice(@be,0,2) ) { 
      my $xir=$rgi.$rt[0]; my $xit=$tgi.$rt[1]; my $txi=$rt[1];
      my $xi=$xon{$xir}||$xon{$xit}; 
      
      unless($xi){ 
        # FIXME slop in blast align break point. check more than +/- 1bp? YES, up to 3bp?
        for(my $j=1; $j<=$BREAKPT_SLOP; $j++) {
          my($rm1,$rp1)=($rt[0]-$j,$rt[0]+$j); 
          my $xj=$xon{$rgi.$rm1} || $xon{$rgi.$rp1}; # only check ref rgi?
          if($xj){ $xi= $xj; last; }
        }
        ## old
        # my($rm1,$rp1)=($rt[0]-1,$rt[0]+1); 
        # my $xj=$xon{$rgi.$rm1}||$xon{$rgi.$rp1}; $xi= $xj; #o: $xi = $xj || $xir; 
        } 

        ## do this even if have xi match to old splice? yes
      if($tdfrag == ISFRAG) { # dont create new splice, but record frag .. how?
        $fragof{$rd}{$td} .= "$txi/$xir;"; ## bug here?
        $isfrag{$td}=$rd; # more info?
        next; # no xon,gxord entries
      } 
      unless($xi){ $xi=$xir; } #  new breakpoint/splice-site
      
      $xon{$xir}=$xon{$xit}=$xi; $validxi++;
      #o: $gxon{$rd}{$xi}++; $gxon{$td}{$xi}++; # not used now?
      unless($gxord{$td}{$txi}){ $gxord{$td}{$txi}= $xi; } #?? want for isself ?
      }  
      
    unless($isself) { # ** FORGOT THIS
    # CHECK: gxun may be adding mistaken unaligned  
    # add if valid xi, gxalign{$td}{$te}=$tb; then check gxun for those ??
    #o if($validxi>1){ $gxalign{$td}{$tb}=$te; } # no, can have many tb->te1,te2,..
    #oo if($validxi>1){ $gxalign{$td}{$tb}{$te}++; } # can have many tb->te1,te2,.. FIXME: rd not td here, so unalign td is kept
    if($validxi>1){ $gxalign{$rd}{$rb}{$re}++; } # can have many tb->te1,te2,.. FIXME: rd not td here, so unalign td is kept
    
    if($lte and $tb - $lte >= $MINUNALIGNBP and not ($tdfrag == ISFRAG)) { 
      $gxun{$td}{$lte}=$tb;  #o: $gxun{$td}{$lte+3}=$tb-3; 
      #?instead? gxun{$td}{$lte}= $gxord{$td}{$lte};  gxun{$td}{$tb}= $gxord{$td}{$tb};
      }
    $lte=$te; # set 1st lte==1 to get endpt, do last te end also?
    }
    }
}
warn "# nin=$nin, nok=$nok, nfrag=$nfrag, nskipnotloc=$skipnotloc, nskipdupfrag=$skipdidfrag, nskipdiffloc=$skipdiffloc\n"; # if debug?    


sub MAIN_END {}
%xon=(); # done with, save mem?

my($XID, %xdid, %xidid, %didid);
# pulled from overinchain2locus.pl
my(%icn,%inc,%dic,%xnchain,%xchainloc,%xchainstat);  # intron/exon splice chains and parts, with tr id

my($rux)= insertUniqExons() if($ADD_UNIQX);

my($rcc)= collectExonChains(); # set xnchain, dic, inc, icn
my($rcl)= assignChainLoci(); # %xchainloc : option?

use constant COMMX_MAXTOP => 9; # 8 or 9? or what? makes diff to commx(), reduce for small nt?
$BYGENE=0 unless($npubtr);

if($BYGENE) {
 
  my( %altset);
  # my @gids= sort values %$locush;
  for my $td (sort keys %gxord) {
    my $pd= pubidof($td);
    # my $gid= $locush->{$td}||$td; #should exist
    my($gid,$ti)= evglocaltid($pd);  
    $altset{$gid}{$ti}= $td;
  }
  
  for my $gid (sort keys %altset) {
    my @ti= sort{$a <=> $b} keys %{$altset{$gid}};
    my @loctr= map{ $altset{$gid}{$_} } @ti;
    my $ntr=@loctr;

    ## handle $USE_ICLOCUS here, separate iclocus in @gid/altset, assignChainLoci info
    
    # my $xn= $xnchain{$id}||"noxchain"; # for commx from collectExonChains
    my @locxchains= map{ $xnchain{$_} || "na" } @loctr;
    my $ntoptr= ($ntr>COMMX_MAXTOP)?COMMX_MAXTOP:$ntr;
    
    my($ncomm,$commx,$commx_ar,$xncomm_hash)= commexons($gid,$ntoptr,\@locxchains); 
    # commex needs @xn exon chain info per loctr from %xchainloc, xnchain, dic, inc, icn,

    print_locus_exontab($gid,\@loctr,$ncomm,$commx,$commx_ar,$xncomm_hash);
  }
  
} else {
  print_exontab();
}

#---------- subs -----------

sub commxval {
  my($gid,$loctr,$ncomm,$commx,$commx_ar,$xncomm_hash)= @_;
  my @commxval=();
  my %commx= map{ $_ => 1 } @$commx_ar; # commx ',.,' spacers??
  my $nt= @$loctr;
  my $xncommtot=0; map{ $xncommtot++ if($xncomm_hash->{$_}); } keys %$xncomm_hash;
  ## return ("noalts") if($ncomm<1 or $nt<2); # return something? "noalts":"nocomm";
  
  for my $td (@$loctr) {
    my $xn=$xnchain{$td}||"";
    my @xn= split",",$xn; # * should be same num @xn == @sym
    my $nc=0; 
    my $hascomm= 0;  
     
    if($ncomm==0) {
      $hascomm= ($nt < 2)?"noalts":"nocomm";
    } elsif($xn =~ m/$commx/) {  #?? commx ',.,' spacers
      for(my $j=0; $j<@xn; $j++) { if($commx{$xn[$j]}){ $nc++; } else { } } # push @csym,$sym[$j];  push @csym,'.';
      $hascomm= $nc."cs"; # skip# ."c:".join("",@csym); 
    } else { 
      # check subsets, off by start/end, or inserted uniq exons
      my($nok,$nmiss)=(0,0);
      for(my $j=0; $j<@xn; $j++) {
        if($commx{$xn[$j]}){ $nc++; $nok++; } #skip: push @csym,$sym[$j]; 
        else { $nmiss++ if($nok and $nok<$ncomm); } #skip  push @csym,'.'; 
      } 
      if($nc>0) { my $pc=($nmiss>0 or $nok < $ncomm)?"ps":"cs"; $hascomm="$nc$pc"; } #skip: :.join("",@csym)
      else { $hascomm="nocomm"; }
    }
    
    my $scb=0; map{ if($xncomm_hash->{$_}>0){ $scb++; } } @xn;
    my ($pcxspan)= pctof( ($hascomm eq "noalts" or $ncomm==0)? 0 : ($hascomm eq "nocomm")? 0 : $nc/$ncomm);
    my ($pcexons)= pctof( ($xncommtot>0) ? $scb/$xncommtot: 0);
    
    ## quick fixup: span==0 for few alts, missed shorter comspan; should fix commx()
    if($pcxspan < 1 and $pcexons > 49 and $hascomm eq "nocomm") { $pcxspan=$pcexons; $hascomm=$scb."cx";} #??
    
    my $commval = "$pcxspan%cs,$pcexons%cx,$hascomm";
       ## $commval.= "$hascomm,cx:$scb/$xncommtot"; #?
    push @commxval, $commval;
  }
  return \@commxval;
}

sub pctof{ my @p; for (@_) { my $p=int(0.5 + 100*$_); push @p, ($p>100)?100:$p; } return @p; }

sub commexons  {  # from commonexons.pl
  my($gn,$nm,$xnchain_ar)=@_; # use @xn global
  
  # handle no alts / noclass quickly
  if($nm < 1) {
    return(0,"",[],{}); #return($ncomm,$xcomm,\@xcomm,\%xnscomm);
  }
  
  #UPD also count common exons as per blasttrset2exons2:altchainsStats()
  # .. all xns{xi} >= pCOM * ($nm+1)
  use constant pCOM => 0.75; # blasttrset2exons2:altchainsStats
  
  my(%xns,%xno,%xnoj,%xnpair); my $maxnx=0;
  for (my $i=0; $i<=$nm; $i++) { 
    my @xi=split",", $xnchain_ar->[$i]; 
    my $nxi=@xi; $maxnx=$nxi if($nxi>$maxnx);
    for(my $j=0; $j<$nxi; $j++) { # keep order
      my $xi=$xi[$j];
      $xns{$xi}++; 
      if($nxi == $maxnx) { $xnoj{$xi}=$j; $xno{$j}=$xi; } # will overwrite xni .. pick longest?
      $xnpair{$xi}{$xi[$j+1]}++ if($j<$nxi-1);
      }
  }
  
  my $mincomxns = 1 + int($nm * pCOM); # should be 7 or 8 for nm max = 9,10
  my %xnscomm= map{ my $v=($xns{$_} >= $mincomxns)?1:0; $_ => $v; } keys %xns; # return this hash

  my @xic=  sort{ $xns{$b}<=>$xns{$a} or $xnoj{$a} <=> $xnoj{$b} or $a cmp $b } keys %xns;
  
  my @xcomm=(); 
  my $mincomm= 1 + int( $maxnx * 0.25 ); #?? what, 1 is enough, want longest commx and max alts/commx
  my $n0test=($maxnx < 3)? 1 : ($maxnx < 6)? 3 : ($maxnx < 9)? 6 : 9; #was 3;
  
  my(@xcommset);
  for(my $i0=0, my $more=1; $i0<$n0test and $more; $i0++) {  # max may not be part of common set? try 0,1,2 ?
    my $xmax= $xic[$i0];
    next unless(exists $xnpair{$xmax});
    my $cmax= $xns{$xmax};  # what of small num: 1,2 ?
    if($cmax < 2){ $more=0; last; }
    my $cmax2= 1 + int($cmax * 0.60); # 0.80 .. 0.75 .. 0.70 .60 ? calc cut from @xic/xns val instead?
    @xcomm=($xmax);
    # ** k counter prevents inf loop here.. why is xia never 0?
    for(my $xia=$xmax, my $k=0; $xia and $k<19; ++$k) {
      my @xib= sort{ $xnpair{$xia}{$b} <=> $xnpair{$xia}{$a}} keys %{$xnpair{$xia}}; 
      $xia=0; last if(@xib == 0);
      for my $xib (@xib) { 
        if($xns{$xib} >= $cmax2) { push @xcomm,$xib; $xia=$xib if(exists $xnpair{$xib}); last; } 
      }
    }
    if(@xcomm>1) { push @xcommset,[@xcomm]; } # recheck if none are >= mincomm
    $more=0 if(@xcomm >= $mincomm and $NOTALLCOMMSETS); #? or run to end of loop, collect all xcommset spans
  }
  
  if(@xcommset > 1 and (@xcomm < $mincomm or not $NOTALLCOMMSETS) ) {
    my($xcbest)= sort{ scalar(@$b) <=> scalar(@$a) } @xcommset; # longest? or/and most alts?
    @xcomm= @$xcbest;
  }
  
  my $ncomm= scalar( @xcomm );
  my $xcomm= join",",@xcomm; #dont need?
  return($ncomm,$xcomm,\@xcomm,\%xnscomm);
}

=item insertUniqExons

  CHECK: gxun may be adding mistakes, aligned exons called as unaligned w/ fake splice in middle of exon
  find this pattern
      tref:  A-B-C-D
      tnew:  A-B-E-C  << -E- is unalined extra, ref has only B-C splice
  gxord has A,B,C align span end-points for both tref,tnew; 
  gxun supposed to collect E insert exon spans  
  BUT gxun now has offset +3/-3 from orig align ends/splices, instead should
  Measure distance of B-end to C-start for both tref, tnew, to detect E unaligned span
  
=cut
    
sub insertUniqExons { # add gxun to gxord
  my $nux=0;
  for my $td (sort keys %gxord) { 
    next unless($gxun{$td});
    my($tgi)= idparts($td);  
    my @xab= sort{$a<=>$b} keys %{$gxalign{$td}}; my $nxab=@xab;
    for my $xub (sort{$a<=>$b} keys %{$gxun{$td}}) { 
      my $xuok=1; 
      my $xue=$gxun{$td}{$xub}||0; # above: gxun{$td}{$lte}=$tb; 
      
      ## need to check all of $gxalign{$td} that may span xub-xue
      for(my $i=0; $i<$nxab and $xuok; $i++) {
        my $xab= $xab[$i]; last if($xab > $xub+$BREAKPT_SLOP);
        my @xae= sort{$b<=>$a} keys %{$gxalign{$td}{$xab}};
        for(my $j=$#xae; $j>=0; $j--) {
          my $xae=$xae[$j];
          if($xae+$BREAKPT_SLOP >= $xue and $xab-$BREAKPT_SLOP <= $xub) { $xuok=0; last; }
        }
      }
      ##OLD: add this check gxalign
      #   my $xae=$gxalign{$td}{$xub}||0; 
      #   if($xae and abs($xae-$xue)<=$BREAKPT_SLOP) { $xuok=0;  }
      #   if($xuok){
      #   for(my $j=1; $j<=$BREAKPT_SLOP; $j++) { 
      #     my($rm1,$rp1)=($xub-$j,$xub+$j); 
      #     $xae= $gxalign{$td}{$rm1}||$gxalign{$td}{$rp1}||0;
      #     if($xae and abs($xae-$xue)<=$BREAKPT_SLOP) { $xuok=0; last; }
      #     } 
      #   }
      #   if($xuok) { 
      #     for(my $j=$xub+$BREAKPT_SLOP; $j<=$xue-$BREAKPT_SLOP; $j++) { if($gxord{$td}{$j}){ $xuok=0; last; } } 
      #   }
      
      if($xuok){ my $xum=int(($xub+$xue)/2);  $gxord{$td}{$xum}=$tgi.$xum.'u'; $nux++; } #*?? 'u' tag ok?
    } 
  }
	warn "#insertUniqExons= $nux\n" if $debug;
	return($nux);
}


sub _sortichain{ return $icn{$b}<=>$icn{$a} or $a cmp $b } # sort by longest splice chain 


sub altchainsStats { 
  my($ic, $altchains)= @_;
  ## change val from 1 to stats: nx, ncommon, nshared, nuniq
  # %inc{exon}{xnchain}= trid
  # %icn{xnchain} = num-exons;  %dic{xnchain}{trid}=1

  use constant pCOM => 0.75; # guess
  my (%xitype,%sumx);
  my @ics=  sort _sortichain keys %{$altchains};
  my $nic= @ics; my $ncom= pCOM * $nic;
  for my $ic (@ics) {
    my @icxn= split",", $ic;
    my(%xtype); ##$xtype{'a'}= scalar(@icnx);
    for my $i (@icxn) {
      my @icx=  grep{ $altchains->{$_} } keys %{$inc{$i}};
      my $nicx= @icx;
      my $xtype;
      if($nicx == 1) { $xtype='u'; } # uniq
      elsif($nicx >= $ncom) { $xtype='c'; } # common
      else { $xtype='s'; } # shared
      $xtype{$xtype}++; $xtype{'a'}++;
      $xitype{$i}{$xtype}++; $xitype{$i}{'a'}++; #?
    }
    my $xtypes= join"", map{ $sumx{$_}+=$xtype{$_}; $_ . $xtype{$_} } sort keys %xtype;
    $altchains->{$ic}= $xtypes; # replaces '1' val
  }
  
  my $sumx= join"", map{ $_ . $sumx{$_} } sort keys %sumx; #? dup counts per ic with same exon
  my @ix= sort keys %xitype; # diff way, count per exon, its type(s)
  my %xtype; for my $ix (@ix){ my @t = sort keys %{$xitype{$ix}}; map{ $xtype{$_}++ } @t; }
  my $xtypes= join"", map{ $_ . $xtype{$_} } sort keys %xtype;
  return($nic,$xtypes,$sumx);
}

sub altchainsOfChain { 
  # collect all tr that contain shared exons,  from %inc{exon}{xnchain}= trid
  my($ic, $altchains)= @_;
  $altchains->{$ic}= 1;  # ic == exon/intron chain of transcript: x1,x2,x3,..,xn
  for my $i (split",", $ic) {
    my @oc= sort _sortichain grep{ not $altchains->{$_} } keys %{$inc{$i}};  
    map{ $altchains->{$_}= 1; } @oc; 
    map{ altchainsOfChain($_, $altchains); } @oc; # recurse
    }
}

sub collectExonChains {
  for my $td (sort keys %gxord) { 
    my %xdid=(); 
    my @xi= sort{ $a <=> $b } keys %{$gxord{$td}};  # per transcript.splice
    my @xn= grep/\w/, map{ my $x=$gxord{$td}{$_}||""; ($xdid{$x}++)?"":$x; } @xi; # to common splice id
    my $xnchain=join",",@xn;  # == ic
    $xnchain{$td}=$xnchain; # for output
    $icn{$xnchain}=scalar(@xn); $dic{$xnchain}{$td}++;
    map{ $inc{$_}{$xnchain}=$td; } @xn; # NOTE: may be many td/xnchain, use @td = keys dic{xnchain}
  }
  # sets globals: (%icn,%dic,%inc);
  my($ntd, $nchain)=(scalar(keys %xnchain), scalar(keys %icn));
	warn "#collectExonChains= $nchain of $ntd ids\n" if $debug;
}

=item about assignChainLoci

  after collectExonChains()
  from overinchain2locus.pl:print_ichains()
  
=cut

sub assignChainLoci { 

  # modify to record, not print, ichain info, from print_ichains()
  my($iloc,$ialt,%locid,%otds, %sumc, %suma); # globals?
  my @ic= sort _sortichain keys %icn; 
  for my $ic (@ic) { 
    my ($tdc,@tdd)= grep{ not $locid{$_}} sort keys %{$dic{$ic}}; 
    next unless($tdc); 
    $iloc++; $ialt=1; 
    my $lid="l$iloc.t$ialt";
    map{ $locid{$_}=$lid } ($tdc,@tdd);  $suma{$ialt}++;
    my %otds=();   
    
    my $locchains="/$ialt:$ic/";
    my %ocall= ();
    altchainsOfChain($ic, \%ocall);
    my($nic,$xtypes,$sumx)= altchainsStats($ic, \%ocall); # UPD
    my @oc= grep{$_ ne $ic} sort _sortichain keys %ocall;  

    my %oty=();
    $oty{'alt'}=1; # $tdc ichain/main ; 
    if(@tdd){ $oty{'dup'} = scalar(@tdd); } # add @tdd icdup count in %oty
    
    for my $oc (@oc) { 
      # upd here, also test if this oc has : see altchainsStats()
      #  a. uniq exons, 
      #  b. uniq chain, ie not just skipped exons of longer chain
    
      my $isub= index($locchains,$oc);
      my $ot= ($isub>0) ? "sub" : "alt";
      if($ot eq "alt") { 
        ++$ialt; $locchains .="$ialt:$oc/"; $lid="l$iloc.t$ialt"; $suma{$ialt}++;
      } else { 
        my $salt=$ialt; # subchain of which alt?
        my $j= rindex($locchains,":",$isub); my $k= rindex($locchains,"/",$j); 
        $salt= substr($locchains,$k+1,$j-$k-1) if($k>=0 and $j>$k); 
        $lid="l$iloc.t$salt"; 
      } 
      
      my @tds=grep{ not $locid{$_}} sort keys %{$dic{$oc}}; 
      if(@tds){ my $tds=join",",@tds; $otds{$oc}="$ot\t$lid\t$tds"; $oty{$ot}++; map{ $locid{$_}= $lid; } @tds; }
    }      
    
    my($strh,$locsum)=  chainsStringsOf($ic,\%otds);

    # * print ichain, ialt, isub set for this ic chain *
    #__ print above before alts 
    my $ns= 1 + $ic =~ tr/,/,/;
    my $salt= join",", map{ my $n=$oty{$_}; "ic$_:$n" } sort keys %oty;
 
    # UPD: add stats in %xchainlocstat{id}, locsum and per td: imax, salt, exon counts: common, alt(shared), uniq      
    
    #O print "\n#icloc\t$locid{$tdc}\t$salt; imax:$ns\n";# too long: iruns:$locsum;
    my $tdlist=$tdc; 
    my $si= $strh->{$ic}||".";
    my $icstat= $ocall{$ic}||"ns"; #UPD
    if(@tdd) { $tdlist=join",",$tdc,@tdd; }
    #O print join("\t","ichain",$locid{$tdc},$tdlist,$ns."i,".$ic,$si)."\n"; 
    #a my $xlocval=  join("\t","ichain",$locid{$tdc},$tdlist,$ns."i,".$ic,$si); #? leave out tds, oc ??
    my $ictype="ichain";
    my $statval= "$icstat;sum:$xtypes/$sumx;imax:$ns;$salt;iruns:$locsum"; #?? all the summary for main tdc
    for my $id ($tdc,@tdd) {
      my $xlocval=  join(":",$ictype,$locid{$id},$si); #? leave out tds, oc ??
      $xchainloc{$id}= $xlocval;
      $xchainstat{$id}= $statval;
      $sumc{$ictype}++;
      $ictype="icdup";
      $statval=$icstat; # for dups 
    }
    
    #__ print in alts loop
    for my $oc (sort _sortichain keys %otds) { 
      my($ot,$lid,$tds)= split"\t",$otds{$oc}; # tds can be list, all with same locid{}
      my $ns= 1 + $oc =~ tr/,/,/;  
      my $si=$strh->{$oc}||".";
      my $icstat= $ocall{$oc}||"ns"; #UPD
      ## tds == list of td ids, some of @tdd
      #O print join("\t","ic$ot",$lid,$tds,$ns."i,".$oc,$si)."\n";  
      #a my $xlocval=  join("\t","ic$ot",$lid,$tds,$ns."i,".$oc,$si); #? leave out tds, oc ??
      for my $id (split",",$tds){
        my $xlocval=  join(":","ic$ot",$lid,$si); #? leave out tds, oc ??
        $xchainloc{$id}= $xlocval;
        $xchainstat{$id}= "$icstat"; #??
        $sumc{"ic$ot"}++; 
        $ot="dup";
        }
      } 
  }
  
  ## summary: class sum sumc,  alt hist? t1..t20 ?
  if($debug){ print STDERR "#assignChainLoci\n"; # warn
    print STDERR "#n_class:"; for my $c (qw(ichain icalt icsub icdup)){ print STDERR" $c=",$sumc{$c}||0; } print STDERR "\n";
    print STDERR "#n_alts :"; for my $i (1..20){ print STDERR " t$i=",$suma{$i}||0; } print STDERR "\n";
    }
}

sub print_locus_exontab {
  my($gid,$loctr_ar,$ncomm,$commx,$commx_ar,$xncomm_hash)= @_;
  
  my $ntr= scalar(@$loctr_ar);
  my $commxval_ar=[];
  if($ncomm) {
   ($commxval_ar)= commxval($gid,$loctr_ar,$ncomm,$commx,$commx_ar,$xncomm_hash);
  }
  
  # for my $id (@$loctr_ar) 
  for(my $i=0; $i<$ntr; $i++) {
    my $id= $loctr_ar->[$i];
    my $xn= $xnchain{$id}||"noxchain"; #? or/and xchainloc{$id}
    my $xloc= $xchainloc{$id}||"noxloc";
    my $xstat= $xchainstat{$id}||"noxstat";
    my $commxval= ($ncomm)? $commxval_ar->[$i] : "nocommx"; # noalts missed
    
    my $isfrag= $isfrag{$id}||0; # any chance of this here?
    my $dval=$gval{$id}||"na"; 
    my $pd= pubidout($id);  
    if($isfrag){ $dval="fragof1:$isfrag/$dval"; } # check for, shouldnt see
    
    if($didid{$id}) { $ndupout++; }
    else { 
      print join("\t","$pd$id","xn:$xn",$xloc,$commxval,$xstat,$dval)."\n";
      # print "$pd$id\txn:$xn\t$xloc\t$xstat\t$dval\n"; 
      $didid{$id}++; }
    
    if(exists $fragof{$id}) {
      my @frid= sort keys %{$fragof{$id}};
      for my $fd (@frid) {
        $dval=$gval{$fd}||"na"; 
        my $xnf= $fragof{$id}{$fd}; # txi/xir; ..
        my $pfd= pubidout($fd);  
        my $xloc= $xchainloc{$fd}||"noxloc";  #??
        my $xstat= $xchainstat{$fd}||"noxstat";
        my $fcommxval= "nocommxf";# frag no 
        if($didid{$fd}) { $ndupout++; }
        else { 
          print join("\t","$pfd$fd","xnf:$xnf",$xloc,$fcommxval,$xstat,"fragof:$id/$dval")."\n";
          # print "$pfd$fd\txnf:$xnf\t$xloc\t$xstat\tfragof:$id/$dval\n"; 
          $didid{$fd}++; }
      }
    }
  } 
}


sub print_exontab {
  
  my @ids=sort keys %gxord; # sort by pubids for output?
  if($npubtr) { 
    #NOTE: this sort doesnt include frags, output below at 1st/longest fragof
    my @sid= sort{ 
      my($ga,$ta)=evglocaltid( pubidof($a));  
      my($gb,$tb)=evglocaltid( pubidof($b)); 
      ($ga cmp $gb or $ta <=> $tb) } @ids;
    @ids=@sid;
    }

if(1) {
  for my $id (@ids) { 
    my($gid)= ($npubtr) ? evglocaltid( pubidof($id)) : $id; 
    print_locus_exontab($gid,[$id],0);
  }
} else {  
  for my $id (@ids) { 
    
    ## this is collectExonChains():
    ## .. $icn{$xnchain}=scalar(@xn); $dic{$xnchain}{$td}++;
    ## ?? need $xnchain{id}=xnchain for here?
    my $xn= $xnchain{$id}||"noxchain"; #? or/and xchainloc{$id}
    my $xloc= $xchainloc{$id}||"noxloc";
    my $xstat= $xchainstat{$id}||"noxstat";
    
    # my %xdid=(); 
    # my @xi= sort{ $a <=> $b } keys %{$gxord{$id}};  
    # my @xn= grep/\w/, map{ my $x=$gxord{$id}{$_}||"na"; ($xdid{$x}++)?"":$x; } @xi;
    # my $xn=join",",@xn; 
    
    #?? drop $icx, with use of assignChainLoci() need only $xn == ichain, @xn, @xi
    # my @ix= map{ my $i= $xidid{$_}||++$XID; $xidid{$_}=$i; sprintf "i%02d",$i; } @x;
    # my $icx=join",",@ix; 
    
    my $isfrag= $isfrag{$id}||0; # any chance of this here?
    my $dval=$gval{$id}||"na"; 
    my $pd= pubidout($id); # ($npubtr)? $podh->{$id}."\t" :"";
    if($isfrag){ $dval="fragof1:$isfrag/$dval"; } # check for, shouldnt see
    if($didid{$id}) { $ndupout++; }
    else { print "$pd$id\txn:$xn\t$xloc\t$xstat\t$dval\n"; $didid{$id}++; }
    #o: else { print "$pd$id\tic:$icx\txn:$xn\t$dval\n"; $didid{$id}++; }
    # ic:i1,i2,i9,i3,.. is format wanted by other soft
    # xn:g1t2x213,g1t2x439,.. is gene-alt-splicept syntax; dval is input align spans
    
    if(exists $fragof{$id}) {
      my @frid= sort keys %{$fragof{$id}};
      # BUG where: @fragof{id} includes not-frags, including longest..
      for my $fd (@frid) {
        $dval=$gval{$fd}||"na"; 
        my $xnf= $fragof{$id}{$fd}; # txi/xir; ..
        my $icf="na";#??
        my $pfd= pubidout($fd); #($npubtr)? $podh->{$fd}."\t" :"";
        my $xloc= $xchainloc{$fd}||"noxloc";  #??
        my $xstat= $xchainstat{$fd}||"noxstat";
        if($didid{$fd}) { $ndupout++; }
        else { print "$pfd$fd\txnf:$xnf\t$xloc\t$xstat\tfragof:$id/$dval\n"; $didid{$fd}++; }
        #o else { print "$pfd$fd\ticf:$icf\txnf:$xnf\tfragof:$id/$dval\n"; $didid{$fd}++; }
      }
    }
  } 
} # else(1) off
}


sub idparts {
  my($id)=@_;
  $id= pubidof($id); #?? ok here
  my($g,$i)= $id =~ m/^(\w+)[t\.](\d+)$/; # FIXME presumes altnum is tail number, t1,t2 or \.1, \.2 ..
  my $gi=$GI{$g}||++$GI; $GI{$g}=$gi; 
  return("g".$gi."t".$i."x", $g, $i);  # simplified gene/alt local id
}

sub evglocaltid{ my($id)=@_; return ($id=~m/(\w+)[t\.](\d+)$/)?($1,$2):($id,99999); }
sub pubidof { my($id)=@_; if($npubtr){ return $podh->{$id}||$id; } else { $id; } }
sub pubidout { my($id)=@_; if($npubtr){ my $pd=$podh->{$id}||"na"; return "$pd\t"; } else { ""; } }


sub readPubidTab { # from genes/inexchains_gcb.pl
	my($pubids)= @_;
	my(%alen,%gscore,%pod,%oda,%aaqual,%locus);  
	our(%altscore,@locids);
	my($ntr,$lpg)=(0,0);

  sub raltscore {
    our(%altscore,@locids);
    my $topas=1; my %lord; 
    for my $i (0..$#locids) { $lord{$locids[$i]}=$i; }
    @locids = sort{ $altscore{$b}<=>$altscore{$a} or $lord{$a}<=>$lord{$b} } @locids;
    for my $i (0..$#locids) { 
      my $lid=$locids[$i]; my $as= $altscore{$lid}||1; 
      if($i==0) { $topas=$as; $altscore{$lid}=1.0; }
      else { $altscore{$lid}= $as/$topas; }  # altscore rel to top score
    }
    @locids=();
  }
  
  return 0 unless(open(F,$pubids));
	while(<F>) {
		next unless(/^\w/);  
		chomp; my @v=split"\t"; 
		##  Public_mRNA_ID originalID PublicGeneID AltNum Class  AAqual pIdAln Notes; Notes=aaref:xxx,chrmap:xxx,pflag:nnn,feq:xxx,...
		my($pd,$od,$cl,$aaq,$pia,$aaref)=@v[0,1,4,5,6,7];  # same as readTrclass
	  $ntr++;
	  
    my($pg,$ti)=  evglocaltid($pd);  # per locus alt qual scores
	  if($lpg and $lpg ne $pg) { raltscore(); @locids=(); }
	  
    my $pflag="";
		#? keep/use pflag qual? 
    # aaref:4625,Sobic.001G542400.3.p,refbest,..
    if($aaref =~ m/aaref:([^;\s]+)/) {  $aaref=$1; }
		$aaref =~ s/^0,0[,]?//; 
		$aaref =~ s/(,pflag:\d+).*/$1/; $pflag=$1;
		$aaref =~ s/,aadup:[^,;\s]*//; # ,aadup:idnnn,refbest, ; ok here? keep ref(best|good|ok) flag
		$aaref ||= "0";
		
		#* add locus/alttr qual scores, for Ig loc reclass? i.e. per loc, each alt is nn% qual of 1st/best alt,
		# so dont reclass fragments/junk as new loci; assume input.pubids is locus/alt sorted
		# abs qual= gscore + aw; maybe + aacomplete?, pflag>0, 
		
		my($aw)= $aaq=~m/^(\d+)/;
		my($gscore)= ($aaref=~m/(\d+)/)?$1:0; # other score wt sum?
		my $aqual= ($aaq=~/partial|bad/)? 0.80: 1.0;
		my $tqual= $gscore + $aw * $aqual;
		
		my $val=join"\t",$aaq,$pia,$cl,$aaref; # add?
		$pod{$pd}=$od; $pod{$od}=$pd;
		$locus{$pd}=$locus{$od}=$pg; # UPD19 added
		$gscore{$pd}= $gscore{$od}=$gscore;
		$alen{$pd}= $alen{$od}= $aw;
		$oda{$pd}= $oda{$od}= $val;
	  $aaqual{$pd}= $aaqual{$od}= $aaq; #? skip readAaQual, but that adds gaps not here
	  
		$altscore{$pd}= $tqual; # global  
		push @locids, $pd; $lpg=$pg;
	  
		} close(F);
		
  if(@locids) { raltscore(); @locids=(); }
	warn "#readPubidTab($pubids)= $ntr\n";
	return ($ntr,\%alen,\%gscore,\%locus,\%pod,\%oda,\%aaqual,\%altscore); #? 
}

sub chainsStringsOf { # from overinchain2locus.pl
  my($mainic,$otds)= @_;
  my(%istr, %lij, %ljc, %lic, %cin, %irun);
  
  my $isrev=($mainic =~ m/\d(r|rk)\b/)?1:0;
  my @oc= sort _sortichain keys %$otds;
  for my $ic ($mainic,@oc) {
    for my $i (split",", $ic) { my($ib)= ($i=~m/i(\d+)/)?$1:$i; 
      $ljc{$ib}++; $lij{$ib}=$i; $lic{$ib}{$ic}++; }
    }
    
  my @in= sort{ $a<=>$b } keys %lij;
  my $nin0= $#in;
  my @krun = (0..$nin0);
  my @lrun = (1) x (1+$nin0);
  my (@inc); my $krun=0;
  for(my $k=0; $k<=$nin0; $k++, $krun++) { 
    my $ib= $in[$k];    
    my $fi= $ljc{$ib};
    my $lic= join",", sort keys %{$lic{$ib}}; # need ichains of in to say if same run..
    $krun[$k]= $krun; $lrun[$krun]=1; # $frun[$k]= $fi;
    for(my $k2= $k+1; $k2 <= $nin0; $k2++) {
      my $jb=$in[$k2];
      my $ljc= join",", sort keys %{$lic{$jb}}; # need ichains of in to say if same run..
      if($ljc{$jb} == $fi and $ljc eq $lic){ $k=$k2; $krun[$k]= $krun; $lrun[$krun]++; } else { last; } 
      }
  }

  my $locusfreq=""; # not useful.. but need %cin set
  my (@frun);
  #UPD: offset c-char by min(@krun) so always have ABCD as 1st set
  my $kcmin=-1; 
  for(my $k=0; $k<=$nin0; $k++) { my $kc= $krun[$k]; $kcmin=$kc if($kcmin<0 or $kc<$kcmin); }
  
  for(my $k=0; $k<=$nin0; $k++) {     
    my $kc= $krun[$k]; my $kcc= $kc - $kcmin;  
    my $c= ($kcc>51)? '_' : ($kcc>=26) ? chr(97+($kcc-26)) : chr(65+$kcc);
    my $ib= $in[$k];
    my $in= $lij{ $ib };
    my $fi= $ljc{ $ib }; my $lrun=$lrun[$kc];
    $frun[$kc]= "$c$kc"."f$fi"."l$lrun";  
    $cin{$in}= $c; # need this
  }
  $locusfreq= (($isrev)?"rev_":"fwd_"). join",",@frun; # for each ichain? or just locus sum?
  
  for my $ic ($mainic,@oc) {
    my @cic;
    for my $in (split",", $ic) { 
      my $c= $cin{$in}; $c||="-"; 
      push @cic, $c;
      }
    $istr{$ic}=  join"",@cic;
  }
  return(\%istr, $locusfreq);
}

__END__

=item FIXME pubids/evglocus

  add here use of publicset/xxx.pubids, or equivalent locus class table
  to filter out spurious alt aligns.
  
=item FIXME fragments

  frag trs align to subset of full alt, esp. long gene, and need classing as frags.
  this is like perfect-frag cdhit -c1 classing, 
  differs in that frag-align can have non-identical ends, or a few miss/indels internally.
 
 * UPD pfrag=0.30 maybe way to small, frags can be largish but subset span in main tr,
    with "small" end diffs .. need to check that those end diffs are small enough to call frag
    try pfrag=0.75 or 0.80 and end-diff < 20 bp, or < 10bp for largish frags?
     
  ** BUT this also will capture the alt-5'start exons that differ outside of CDS unless
  input blast align includes those 5'UTR extended cds aligns.
  
  eg case008t.pubids, nalt=66 (36 drop for short/utrbad), has 1 refat best match, of 4 refat alts;
  AtEvg3g000008t1	tidbarab6roo3dno1sridk77Loc3527	main	3600,91%,complete	100/100/.	aacons:3,pflag:0,tscore:3600,	tidbarab6roo3dno1sridk77Loc3527
  align bspans of tidbarab6roo3dno1sridk77Loc3527, self=1-10803;	
  .. most are single span align to shorter tr, a few w/ multi-exon align may be valid alts
  .. modify here to skip/flag subset aligns, that span most of frag tr
  251-10803;337-10803;531-10803;730-10803;979-10803;1492-10803;1-9113;1-8805;251-9113;251-8805;531-9113;
  2341-10803;531-8805;730-9113;730-8805;979-9113;979-8805;1492-9113;2011-9568;1492-8805;
  2854-8805,8997-9114,9115-9298;1-5802;251-5802;337-5802;730-5802;895-5802;5821-10803;979-5802;2548-7239;
  1-3987;1-3899;3910-7453;7443-7501,7519-10803;7456-10803;1-2934;531-2934;730-2934;531-2572;979-2934;
  7246-8763,8853-8988,9113-9298;8545-9873;8872-10803;4708-6297;1-1442;1492-2934;337-1740;2518-3888;
  10122-10803;9649-10803;1-1078;5203-6287;8872-9873;1-897;8596-9453;337-1078;1-707;2281-2934;337-963;
  531-1078;9649-10215;251-707;9358-9873;7513-7883,7891-7996;531-963;9934-10410;7513-7884,7914-7996;2518-2934;
  730-1078;339-680;2011-2234;1808-2064;1115-1360;10174-10410;9649-9873;7963-8178;1-182;339-475;339-493;
  
  eg frags to g008 in btall, 100aa to 300aa, to one 600aa x this 3600aa main tr
  .. all are almost fully aligned frags to main, with some slop at ends
 grep ^tidbarab6roo3dno1sridk77Loc3527 tmpfiles/bevgc_evg3weed1rdnrcd1x-self100.btall | grep arab6roo3dntatrinLocDN40350c3
tidbarab6roo3dno1sridk77Loc3527	arab6roo3dntatrinLocDN40350c3g2t1	2555	1329	1329	10803	1938	8545-9873/8545-9873	610-1938:+
tidbarab6roo3dno1sridk77Loc3527	arab6roo3dntatrinLocDN40350c3g2t5	1927	1002	1002	10803	1002	8872-9873/8872-9873	1-1002:+
tidbarab6roo3dno1sridk77Loc3527	arab6roo3dntatrinLocDN40350c3g1t3	1258	654	654	10803	687	2281-2934/2281-2934	1-654:+
tidbarab6roo3dno1sridk77Loc3527	arab6roo3dntatrinLocDN40350c3g3t4	1206	627	627	10803	639	337-963/337-963	1-627:+
tidbarab6roo3dno1sridk77Loc3527	arab6roo3dntatrinLocDN40350c3g2t3	992	516	516	10803	516	9358-9873/9358-9873	1-516:+
tidbarab6roo3dno1sridk77Loc3527	arab6roo3dntatrinLocDN40350c3g1t1	802	417	417	10803	450	2518-2934/2518-2934	1-417:+
tidbarab6roo3dno1sridk77Loc3527	arab6roo3dntatrinLocDN40350c3g3t8	473	246	246	10803	246	1115-1360/1115-1360	1-246:+
tidbarab6roo3dno1sridk77Loc3527	arab6roo3dntatrinLocDN40350c3g2t6	433	225	225	10803	225	9649-9873/9649-9873	1-225:+
tidbarab6roo3dno1sridk77Loc3527	arab6roo3dntatrinLocDN40350c3g1t5	431	224	224	10803	306	2011-2234/2011-2234	1-224:+
  
=item FIXME paralogs

  trasm-self98.blastn includes paralog aligns;
  need either to reformulate this pipe to separate out paralogs, or
  reblast trasm -self100 to get only perfect aligned exons, and 
  add criteria ? for min-align to count (exon size), same as MINUNALIGNBP ?
  
  * use makeblastscore.pl -pidentmin 99.999 : cuts out some paralts, maybe all feasible
    ref_arath16apnrcd1x.i100u20exontab.icloc  -pident 100
      #sum nloc=26889, nuni=19396, nalt=6694, npar=799
    ref_arath16apnrcd1x.u20exontab.icloc
      #sum nloc=26085, nuni=18394, nalt=6505, npar=1186

  * add to makeblastscore output of mismatches/indels per span?
  .. need 1+ exon span of 100% ident to call as alt vs paralog/hetero/read-err
  
  eg. 'AT5G28400|AT3G61780|AT5G28320' ref_arath16ap.pubids 
  Note: 99% ident is not 100%, all true alts should have 100% align common exons  
Atref9Evm000717t1    AT3G61780.1     main    1121,86%,complete       99/23/. pflag:0,tscore:1121,
Atref9Evm000717t2    AT5G28400.1     alt     973,100%,complete       99/18/. pflag:0,tscore:973,
Atref9Evm000717t3    AT5G28320.1     alt     927,100%,complete       99/23/. pflag:0,tscore:927,
  zgrep '^AT3G61780.1' ref_arath16apnrcd1x-self98.blastn.gz
AT3G61780.1     AT3G61780.1     100.000 3366    0       0       1       3366    1       3366    0.0     6472
AT3G61780.1     AT5G28400.1     98.692  535     7       0       2271    2805    1941    2475    0.0     988
AT3G61780.1     AT5G28320.1     98.872  532     6       0       2260    2791    1744    2275    0.0     988
AT3G61780.1     AT5G28320.1     98.113  106     2       0       615     720     495     600     1.39e-47        192
  ref_arath16apnrcd1x.u20exontab.icloc
ichain  l6.t1   AT5G28320.1     FFGHBHIBJJKLLLDDDMDN
icalt   l6.t2   AT5G28400.1     FFOHHOJJOLLLDODDDN
icalt   l6.t3   AT3G61780.1     ABBCDDDDE  <<DD^^ shared exons, BUT has mismatches
  
=item usage

$evigene/scripts/makeblastscore3.pl -pmin 0.01 -CDSSPAN  -showspan=2 -tall \
  -sizes ref_arath16apnrcd1x.cds.qual  ref_arath16apnrcd1x-self98.blastn.gz > ref_arath16apnrcd1x-self98.btall
 ADD opt for self100.btall: -pIDENTMIN 99.999 or 100
 
#* FIXME: sort in this prog unless -sorted flag
sort -k7,7nr -k2,2 -k6,6nr -k1,1 ref_arath16apnrcd1x-self98.btall |\
  $evigene/scripts/genes/blasttrset2exons.pl  > ref_arath16apnrcd1x.exontab

cut -f1,2 ref_arath16apnrcd1x.exontab | $evigene/scripts/genes/overinchain2locus.pl \
  > ref_arath16apnrcd1x.exontab.icloc

head caseG77080a.selfbtall 
Query   Source  Bits    Ident   Align   Qlen    Slen    Qspan   Sspan
AT1G77080.2     AT1G77080.8     1029    535     535     579     672     1-579/1-185,232-291,292-579     1-185,244-305,385-672:+
AT1G77080.2     AT1G77080.4     1025    533     533     579     591     1-579/1-185,232-579     1-185,244-591:+
AT1G77080.2     AT1G77080.5     829     431     431     579     522     1-477/1-185,232-477     1-185,244-489:+
AT1G77080.2     AT1G77080.9     673     350     350     579     510     232-579/232-291,292-579 82-143,223-510:+
AT1G77080.2     AT1G77080.6     477     248     248     579     450     232-477/232-291,292-477 82-143,223-408:+
AT1G77080.2     AT1G77080.7     473     246     246     579     369     232-477/232-477 82-327:+
AT1G77080.2     AT1G77080.3     398     207     207     579     426     232-436/232-291,292-436 82-143,223-367:+
AT1G77080.3     AT1G77080.9     706     367     367     426     510     1-367/1-367     1-367:+
AT1G77080.3     AT1G77080.6     706     367     367     426     450     1-367/1-367     1-367:+

grep -v 'AT1G77080.1.  ' caseG77080a.xonchainsu | cut -f1,2,3
AT1G77080.2	ic:i10,i02,i15,i03,i04,i05,i06,i07,i08,i09	xn:g1t8x1,g1t8x185,g1t2x208,g1t8x244,g1t8x305,g1t8x385,g1t8x529,g1t8x534,g1t8x570,g1t8x672
AT1G77080.3	ic:i01,i02,i03,i16,i17,i06,i14	xn:g1t9x1,g1t8x185,g1t8x244,g1t4x303,g1t3x183,g1t8x529,g1t12x437
AT1G77080.4	ic:i10,i02,i18,i03,i04,i05,i06,i07,i08,i09	xn:g1t8x1,g1t8x185,g1t4x214,g1t8x244,g1t8x305,g1t8x385,g1t8x529,g1t8x534,g1t8x570,g1t8x672
AT1G77080.5	ic:i10,i02,i19,i03,i04,i05,i06,i07,i08	xn:g1t8x1,g1t8x185,g1t5x214,g1t8x244,g1t8x305,g1t8x385,g1t8x529,g1t8x534,g1t8x570
AT1G77080.6	ic:i01,i02,i03,i16,i20,i06,i07,i08,i12	xn:g1t9x1,g1t8x185,g1t8x244,g1t4x303,g1t6x183,g1t8x529,g1t8x534,g1t8x570,g1t11x531
AT1G77080.7	ic:i01,i02,i03,i04,i05,i06,i07,i08,i12	xn:g1t9x1,g1t8x185,g1t8x244,g1t8x305,g1t8x385,g1t8x529,g1t8x534,g1t8x570,g1t11x531
AT1G77080.8	ic:i10,i02,i21,i03,i04,i22,i05,i06,i07,i08,i09	xn:g1t8x1,g1t8x185,g1t8x214,g1t8x244,g1t8x305,g1t8x345,g1t8x385,g1t8x529,g1t8x534,g1t8x570,g1t8x672
AT1G77080.9	ic:i01,i02,i03,i16,i23,i06,i07,i08,i09	xn:g1t9x1,g1t8x185,g1t8x244,g1t4x303,g1t9x183,g1t8x529,g1t8x534,g1t8x570,g1t8x672


=item blasttrset2exons2 : merge ichains parts

exons2:  xn: chain should be same for old,new; added ichain col4, removed ic: col3
  ichain col is squish of prior ichains table:
    AT1G67120.2	ichain:l2551.t1:ACBCA
    AT1G67120.1	icalt:l2551.t2:ACCA   <<  subchain of t1, skips B exon: need note
    
    AT3G02260.1	ichain:l57.t1:ACDACCDACCC  : 2 uniq exons
    AT3G02260.4	icalt:l57.t2:ACBACCACCC  : 1 uniq exon, skips 1 of t1
    AT3G02260.3	icalt:l57.t3:ACACCACCC  << t3 same chain, subchain of t2, skips B,D exons
    AT3G02260.2	icalt:l57.t3:ACACCACCC  << t3 .. call  this redundant *?
    
    AT4G17140.3	icalt:l19.t2:HCMANJHPDCLF
    AT4G17140.2	icalt:l19.t3:HCNHPECKLOF
    AT4G17140.1	ichain:l19.t1:HCGANIHPCKLBF

  longest cds <> longest ichain, 
  common xpat here: DGMAC + extras at 3'end; check some parsing? DGMACE.KB > DGMACEa, DGMACE.KI > DGMACEb,   
    AT3G63340.2	icalt:l350.t6:DGJMAN	1-3228;1-811,909-3228;1-811,909-1537;1-811,909-1537;1-811,909-1537;1-811,909-1537;
    AT3G63340.7	icalt:l350.t7:DGMAN	1-809,810-3129;1-3129;1-1438;1-1438;1-1438;1-1438;
    AT3G63340.1	icalt:l350.t2:DGMACEKB	1-811,812-1438;1-1438;1-3126;1-3030;1-3007;1-2544;
    AT3G63340.12	icalt:l350.t3:DGMAHFO	fragof:AT3G63340.1/1-811,812-1438,1515-1587;1-1438,1515-1587;1-1515,1516-1587;
    AT3G63340.8	ichain:l350.t1:DGMACEKI	1-811,812-1438;1-1438;1-3030;1-3042;1-3007;1-2544;
    AT3G63340.3	icalt:l350.t4:DGMACEL	1-811,812-1438;1-1438;1-3007;1-3007;1-3015;1-2544;
    AT3G63340.11	icalt:l350.t5:DGMACP	1-811,812-1438;1-1438;1-2544;1-2544;1-2544;1-2571;
     
old blast2exons, recheck
sort -k7,7nr -k2,2 -k6,6nr -k1,1 ref_arath16apnrcd1x-self100.btall | env pubids=ref_arath16ap.pubids pfrag=0.80 minunalign=20 $evigene/scripts/genes/blasttrset2exons.pl > ref_arath16apnrcd1x.pubfrag80u20exons08
#readPubidTab(ref_arath16ap.pubids)= 39866
# nin=82758, nok=75866, nfrag=697, nskipnotloc=2596, nskipdupfrag=3180, nskipdiffloc=419

head ref_arath16apnrcd1x.pubfrag80u20exons07 | cut -f2,4
AT1G67120.2	xn:g1t1x1,g1t1x9589,g1t1x9599,g1t1x9610,g1t1x16203
AT1G67120.1	xn:g1t1x1,g1t1x9589,g1t1x9610,g1t1x16203
AT3G02260.1	xn:g2t1x1,g2t1x292,g2t1x305,g2t1x319,g2t1x593,g2t1x7235,g2t1x7245,g2t1x7255,g2t1x7667,g2t1x7680,g2t1x15297
AT3G02260.4	xn:g2t1x1,g2t1x292,g2t2x305,g2t1x319,g2t1x593,g2t1x7235,g2t1x7255,g2t1x7667,g2t1x7680,g2t1x15297
AT3G02260.3	xn:g2t1x1,g2t1x292,g2t1x319,g2t1x593,g2t1x7235,g2t1x7255,g2t1x7667,g2t1x7680,g2t1x15297
AT3G02260.2	xn:g2t1x1,g2t1x292,g2t1x319,g2t1x593,g2t1x7235,g2t1x7255,g2t1x7667,g2t1x7680,g2t1x15297
AT5G23110.1	xn:g3t1x1,g3t1x14121
AT4G17140.3	xn:g4t1x1,g4t1x2638,g4t1x2661,g4t1x2684,g4t1x2738,g4t1x3759,g4t1x4780,g4t1x9633,g4t1x9653,g4t1x9674,g4t1x9750,g4t1x12660
AT4G17140.2	xn:g4t1x1,g4t1x2638,g4t1x2738,g4t1x4780,g4t1x9633,g4t2x9545,g4t1x9674,g4t2x9640,g4t1x9750,g4t2x9695,g4t1x12660
AT4G17140.1	xn:g4t1x1,g4t1x2638,g4t3x2650,g4t1x2684,g4t1x2738,g4t3x2730,g4t1x4780,g4t1x9633,g4t1x9674,g4t2x9640,g4t1x9750,g4t3x9687,g4t1x12660

new blast2exons2,
sort -k7,7nr -k2,2 -k6,6nr -k1,1 ref_arath16apnrcd1x-self100.btall | env pubids=ref_arath16ap.pubids pfrag=0.80 minunalign=20 $evigene/scripts/genes/blasttrset2exons2.pl > ref_arath16apnrcd1x.pubfrag80u20exon2a 
#readPubidTab(ref_arath16ap.pubids)= 39866
# nin=82758, nok=75866, nfrag=697, nskipnotloc=2596, nskipdupfrag=3180, nskipdiffloc=419

head ref_arath16apnrcd1x.pubfrag80u20exon2a | cut -f2,3,4
AT1G67120.2	xn:g1t1x1,g1t1x9589,g1t1x9599,g1t1x9610,g1t1x16203	ichain:l2551.t1:ACBCA
AT1G67120.1	xn:g1t1x1,g1t1x9589,g1t1x9610,g1t1x16203	icalt:l2551.t2:ACCA
AT3G02260.1	xn:g2t1x1,g2t1x292,g2t1x305,g2t1x319,g2t1x593,g2t1x7235,g2t1x7245,g2t1x7255,g2t1x7667,g2t1x7680,g2t1x15297	ichain:l57.t1:ACDACCDACCC
AT3G02260.4	xn:g2t1x1,g2t1x292,g2t2x305,g2t1x319,g2t1x593,g2t1x7235,g2t1x7255,g2t1x7667,g2t1x7680,g2t1x15297	icalt:l57.t2:ACBACCACCC
AT3G02260.3	xn:g2t1x1,g2t1x292,g2t1x319,g2t1x593,g2t1x7235,g2t1x7255,g2t1x7667,g2t1x7680,g2t1x15297	icalt:l57.t3:ACACCACCC
AT3G02260.2	xn:g2t1x1,g2t1x292,g2t1x319,g2t1x593,g2t1x7235,g2t1x7255,g2t1x7667,g2t1x7680,g2t1x15297	icalt:l57.t3:ACACCACCC
AT5G23110.1	xn:g3t1x1,g3t1x14121	ichain:l25276.t1:AA
AT4G17140.3	xn:g4t1x1,g4t1x2638,g4t1x2661,g4t1x2684,g4t1x2738,g4t1x3759,g4t1x4780,g4t1x9633,g4t1x9653,g4t1x9674,g4t1x9750,g4t1x12660	icalt:l19.t2:HCMANJHPDCLF
AT4G17140.2	xn:g4t1x1,g4t1x2638,g4t1x2738,g4t1x4780,g4t1x9633,g4t2x9545,g4t1x9674,g4t2x9640,g4t1x9750,g4t2x9695,g4t1x12660	icalt:l19.t3:HCNHPECKLOF
AT4G17140.1	xn:g4t1x1,g4t1x2638,g4t3x2650,g4t1x2684,g4t1x2738,g4t3x2730,g4t1x4780,g4t1x9633,g4t1x9674,g4t2x9640,g4t1x9750,g4t3x9687,g4t1x12660	ichain:l19.t1:HCGANIHPCKLBF

=item tests

  arabidopsis/tr2aacds_test1908f/try9bestevgc/
  sort -k7,7nr -k2,2 -k6,6nr -k1,1 bevgc_evg3weed1rdnrcd1x-self100.btall | env pubids=../publicset/bevgc_evg3weed1rd.pubids pfrag=0.33 minunalign=20 $evigene/scripts/genes/blasttrset2exons.pl > bevgc_evg3weed1rdnrcd1x.pubfrg20exons
  #readPubidTab(../publicset/bevgc_evg3weed1rd.pubids)= 401404
  # nin=10220502, nok=5097769, nfrag=204776, nskipnotloc=4781955, nskipdiffloc=136002
  .. upd same stats, fixed frag pubids
  # nin=10220502, nok=5097769, nfrag=204776, nskipnotloc=4781955, nskipdiffloc=136002
  .. upd skip dup frag ids
  # nin=10220502, nok=4459571, nfrag=72573, nskipnotloc=4781955, nskipdiffloc=136002
  532756 bevgc_evg3weed1rdnrcd1x.pubfrg20exons01
  424293 bevgc_evg3weed1rdnrcd1x.pubfrg20exons02 : frag dup bug fixed ?? 
  401404 bevgc_evg3weed1rdnrcd1x.pubfrg20exons   : upd3 dedup

  .. upd3 dups
  # nin=10220502, nok=4459571, nfrag=72573, nskipnotloc=4781955, nskipdupfrag=770401, nskipdiffloc=136002
  
  * Bug got dup pubids, from frags 
grep '^AtEvg3g000001t' bevgc_evg3weed1rdnrcd1x.pubfrg20exons | cut -f1 | sort | uniq -c | grep -v ' 1 ' | head
  11 AtEvg3g000001t70
  19 AtEvg3g000001t81
   5 AtEvg3g000001t84
    .. still have some dup pubids : AtEvg3g000001t81 fragof1
    .. must be first entry for not-frag? then fragof?
AtEvg3g000001t81	arab6roo3dntr6akvelvk43Loc903t27	ic:i364872	xn:g1t1x1033	fragof1:arab6roo3dntr6akvelvk77Loc1301t6/1-82,83-2922;
AtEvg3g000001t81	arab6roo3dntr6akvelvk43Loc903t27	icf:na	xnf:82/g1t1x1115;83/g1t1x1211;2922/g1t1x4050;	fragof:arab6roo3dntr6akvelvk77Loc1301t6/1-82,83-2922;
     
  ..
  arabidopsis/tr2aacds_test1908f/try9refat/
  sort -k7,7nr -k2,2 -k6,6nr -k1,1 ref_arath16apnrcd1x-self100.btall | \
   env pfrag=0.33 minunalign=20 $evigene/scripts/genes/blasttrset2exons.pl > ref_arath16apnrcd1x.defragu20exontab
  nin=82758, nok=82579, nfrag=179, nskipnotloc=0, nskipdiffloc=0 
   env pubids=ref_arath16ap.pubids pfrag=0.33 minunalign=20
  nin=82758, nok=79715, nfrag=28, nskipnotloc=2596, nskipdiffloc=419

cat ref_arath16apnrcd1x.pubfragu20exons | cut -f1,2,4 | sort -k1,1 | head -20
Atref9Evm000001t1	AT1G67120.2	xn:g1t2x1,g1t2x9589,g1t2x9599,g1t2x9610,g1t2x16203
Atref9Evm000001t2	AT1G67120.1	xn:g1t2x1,g1t2x9589,g1t2x9610,g1t2x16203

Atref9Evm000002t1	AT3G02260.1	xn:g2t1x1,g2t1x292,g2t1x305,g2t1x319,g2t1x593,g2t1x7235,g2t1x7245,g2t1x7255,g2t1x7667,g2t1x7680,g2t1x15297
Atref9Evm000002t2	AT3G02260.4	xn:g2t1x1,g2t1x292,g2t4x305,g2t1x319,g2t1x593,g2t1x7235,g2t1x7255,g2t1x7667,g2t1x7680,g2t1x15297
Atref9Evm000002t3	AT3G02260.3	xn:g2t1x1,g2t1x292,g2t1x319,g2t1x593,g2t1x7235,g2t1x7255,g2t1x7667,g2t1x7680,g2t1x15297
Atref9Evm000002t4	AT3G02260.2	xn:g2t1x1,g2t1x292,g2t1x319,g2t1x593,g2t1x7235,g2t1x7255,g2t1x7667,g2t1x7680,g2t1x15297

Atref9Evm000003t1	AT5G23110.1	xn:g3t1x1,g3t1x14121
Atref9Evm000004t1	AT4G17140.3	xn:g4t3x1,g4t3x2638,g4t3x2661,g4t3x2684,g4t3x2738,g4t3x3759,g4t3x4780,g4t3x9633,g4t3x9653,g4t3x9674,g4t3x9750,g4t3x12660
Atref9Evm000004t2	AT4G17140.2	xn:g4t3x1,g4t3x2638,g4t3x2738,g4t3x4780,g4t3x9633,g4t2x9545,g4t3x9674,g4t2x9640,g4t3x9750,g4t2x9695,g4t3x12660
Atref9Evm000004t3	AT4G17140.1	xn:g4t3x1,g4t3x2638,g4t1x2650,g4t3x2684,g4t3x2738,g4t1x2730,g4t3x4780,g4t3x9633,g4t3x9674,g4t2x9640,g4t3x9750,g4t1x9687,g4t3x12660

Atref9Evm000005t1	AT1G48090.2	xn:g5t2x1,g5t2x73,g5t2x676,g5t2x709,g5t2x742,g5t2x6548,g5t2x6588,g5t2x6628,g5t2x12522
Atref9Evm000005t2	AT1G48090.3	xn:g5t3x1,g5t2x73,g5t2x676,g5t3x652,g5t2x742,g5t2x6548,g5t3x6531,g5t2x6628,g5t2x12522
Atref9Evm000005t3	AT1G48090.1	xn:g5t2x1,g5t2x73,g5t2x676,g5t1x709,g5t2x742,g5t1x6546,g5t2x6548,g5t2x6628,g5t2x12522
Atref9Evm000005t4	AT1G48090.6	xn:g5t3x1,g5t2x73,g5t2x676,g5t2x742,g5t1x6546,g5t6x6465,g5t2x12522
Atref9Evm000005t5	AT1G48090.5	xn:g5t3x1,g5t2x73,g5t2x676,g5t2x742,g5t2x6548,g5t2x6628,g5t1x6546,g5t2x12522

Atref9Evm000006t1	AT2G17930.1	xn:g6t1x1,g6t1x11577
Atref9Evm000007t1	AT4G36080.1	xn:g7t1x1,g7t1x4511,g7t1x4557,g7t1x4602,g7t1x4629,g7t1x11505
Atref9Evm000007t2	AT4G36080.3	xn:g7t1x1,g7t1x4511,g7t3x4532,g7t1x4557,g7t1x4629,g7t1x11505
Atref9Evm000007t3	AT4G36080.2	xn:g7t1x1,g7t1x4511,g7t1x4602,g7t2x4525,g7t1x4629,g7t1x11505

=cut
