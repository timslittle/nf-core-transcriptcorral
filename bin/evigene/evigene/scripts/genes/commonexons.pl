#!/usr/bin/perl
# commonexons.pl < geneset.exontab

=item about 

  find common/constituative exons shared by longest alternates
  from geneset.exontab of  evigene/scripts/genes/blasttrset2exons2.pl
  
  exontab is ordered by evgene locus, longest alts first

=item FIXME paralogs

  some of evg classed loci contain paralogs, with diff ichain loci ids
  using evg gene id will call those diff icloci as 'nocomm'
  use ic locus ids to rectify ??    eg.
  >> paralog locus A, ichain:l917
Atref9Evm000688t1       AT4G00930.1     ichain:l917.t1:BBEDFC   3c:BBE...       a6c3s2u1;sum:a9c3s2u4/a23c14s5u4;imax:6;icalt:4,icsub:1;iruns:fwd_A0f1l1,B1f5l2,C2f1l1,D3f3l1,E4f4l1,F5f2l1,G6f1l1,H7f1l1
Atref9Evm000688t2       AT4G00930.4     icalt:l917.t3:BBEG      3c:BBE. a4c3u1
Atref9Evm000688t3       AT4G00930.2     icsub:l917.t1:BBEDF     3c:BBE..        a5c3s2
Atref9Evm000688t4       AT4G00930.3     icalt:l917.t2:BBEDH     3c:BBE..        a5c3s1u1
Atref9Evm000688t5       AT4G00930.5     icalt:l917.t4:BBA       2p:BB.  a3c2u1
  >> paralog locus B, ichain:l3616
Atref9Evm000688t6       AT5G37190.1     ichain:l3616.t1:EACD    nocomm  a4c2s1u1;sum:a5c2s1u2/a10c6s2u2;imax:4;icalt:2,icsub:1;iruns:fwd_A0f3l1,B1f1l1,C2f2l1,D3f1l1,E4f3l1
Atref9Evm000688t7       AT5G37190.2     icsub:l3616.t1:EAC      nocomm  a3c2s1
Atref9Evm000688t8       AT5G37190.3     icalt:l3616.t2:EAB      nocomm  a3c2u1
    
    separate by ic locus in putg(), call putgsub(icA), putgsub(icB)
=cut  

use strict;

my $USE_ICLOCUS= $ENV{iclocus}||0; # default on ? seems to work on at ref
my $ONLY_ALTS=  $ENV{onlyalts}||0;
my $DROP_FRAG= ($ENV{keepfrag})?0:1;
my $NOTALLCOMMSETS= $ENV{allcommsets}?0:1; # test
my $debug=  $ENV{debug}||0;

my($lgd,@tval,@xn);
while(<>) {
  next if(/^\W/);
  my @v=split;
  my($pd,$td,$xn,$icv,$icstat,$aspans)=@v;
  $xn=~s/^xn://; $v[2]= $xn; # drop tag
  ## FIXME fragof $aspans tag should be skipped for commx .. 
  (my $gd=$pd) =~ s/t\d+$//; 
  if($lgd ne $gd){ putg($lgd) if($lgd); @tval=(); @xn=(); }  
  push @tval, \@v;  
  push @xn, $xn;  $lgd=$gd;
}
putg($lgd);

#------ subs -----------

sub commx  {
  my($gn,$nm)=@_; # use @xn global
  
  # handle no alts / noclass quickly
  if($nm < 1) {
    return(0,"",[],[],{}); #return($ncomm,$xcomm,\@xcomm,\@xcommset,\%xnscomm);
  }
  
  #UPD also count common exons as per blasttrset2exons2:altchainsStats()
  # .. all xns{xi} >= pCOM * ($nm+1)
  use constant pCOM => 0.75; # blasttrset2exons2:altchainsStats
  
  my(%xns,%xno,%xnoj,%xnpair); my $maxnx=0;
  for (my $i=0; $i<=$nm; $i++) { 
    my @xi=split",", $xn[$i]; 
    my $nxi=@xi; $maxnx=$nxi if($nxi>$maxnx);
    for(my $j=0; $j<$nxi; $j++) { # keep order
      my $xi=$xi[$j];
      $xns{$xi}++; 
      #o $xno{$xi}{$j}++; $xnoj{$j}{$xi}++; # not used now
      if($nxi == $maxnx) { $xnoj{$xi}=$j; $xno{$j}=$xi; } # will overwrite xni .. pick longest?
      #? or should this be xnpair{$xi}{$lxi} ? lxi = last xi
      $xnpair{$xi}{$xi[$j+1]}++ if($j<$nxi-1);
      #? $xnpairb{$xi}{$xi[$j-1]}++ if($j>0);
      }
  }
  
  my $mincomxns = 1 + int($nm * pCOM); # should be 7 or 8 for nm max = 9,10
  my %xnscomm= map{ my $v=($xns{$_} >= $mincomxns)?1:0; $_ => $v; } keys %xns; # return this hash

  my @xic=  sort{ $xns{$b}<=>$xns{$a} or $xnoj{$a} <=> $xnoj{$b} or $a cmp $b } keys %xns;
  
  my @xcomm=(); 
  my $mincomm= 1 + int( $maxnx * 0.25 ); #?? what, 1 is enough, want longest commx and max alts/commx
  ##not all xic: int( scalar(@xic) * 0.25 ); #?? what
  my $n0test=($maxnx < 3)? 1 : ($maxnx < 6)? 3 : ($maxnx < 9)? 6 : 9; #was 3;
  #o my $n0test=($maxnx < 3)? 1 : 3;
  #x my $n0test=($maxnx < 3)? 1 : ($maxnx > 30) ? 14 : int($maxnx/2);
  ## case: Hstra9Evm000248t hangs here at n0test > 3, with maxnx >= 25
  ##  ** n0test matters for wacky trasm, 9 is better than 3, gets 3'end common missed w/ few tests
  ## but larger nums cause probs w/ some .. how/where?
  
  # my $xicb = 1 + int(scalar(@xic) * 0.40); # index top quarter?
  # my $cmax3= 1 + $xns{ $xic[$xicb] }; # replace cmax2 ? # not useful
  
  my(@xcommset);
  for(my $i0=0, my $more=1; $i0<$n0test and $more; $i0++) {  # max may not be part of common set? try 0,1,2 ?
    my $xmax= $xic[$i0];
    next unless(exists $xnpair{$xmax});
    my $cmax= $xns{$xmax};  # what of small num: 1,2 ?
    if($cmax < 2){ $more=0; last; }
    my $cmax2= 1 + int($cmax * 0.60); # 0.80 .. 0.75 .. 0.70 .60 ? calc cut from @xic/xns val instead?
    ## loop here for @xcomm chain
    @xcomm=($xmax);
    # $backup= ($xnoj1{$xmax} > $xnoj1{$xic[$i0+1]})?1:0; #?? build common backwards, should test both ends from max
    #?? go thru all of xnpair span? may have comm at two end points, need some xcomm spacer, '.'?
    #.. maybe only sure way is to call commx() 2nd time, w/ 1st commx cut from chains.
    # ** k counter prevents inf loop here.. why is xia never 0?
    for(my $xia=$xmax, my $k=0; $xia and $k<19; ++$k) {
      my @xib= sort{ $xnpair{$xia}{$b} <=> $xnpair{$xia}{$a}} keys %{$xnpair{$xia}}; 
      $xia=0; last if(@xib == 0);
      #n my $xib= $xib[0] || 0;
      #n if($xns{$xib} >= $cmax2) { push @xcomm,$xib; $xia=$xib if(exists $xnpair{$xib}); } #?? only 1 xib
      #SS else { push @xcomm,'.'; $xia=$xib if($xnpair{$xib}); } # spacer to find 2nd comm span ?? problematic
      ## use this way --vv-- sometimes 1st xib poor choice
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
  
  #SS while($xcomm[-1] eq '.'){ pop(@xcomm); } # end spacers
  #SS my $ncomm= scalar( grep{ $_ ne '.' } @xcomm );
  ## maybe want \@xcommset
  my $ncomm= scalar( @xcomm );
  my $xcomm= join",",@xcomm; #dont need?
  return($ncomm,$xcomm,\@xcomm,\@xcommset,\%xnscomm);
}

use constant MAXTOP => 8; # 8 or 9? or what? makes diff to commx(), reduce for small nt?

sub putg {
  my($gn)=@_;
  my $nt= @tval - 1; 
  
    ## FIXME fragof $aspans tag should be skipped for commx .. change to @xn, $nm needed from @tval
    ## DROP_FRAG maybe want in output table, w/o commx testing ?? need to pass @ok index list to commx()
  my @fragval=();
  if($DROP_FRAG) {
    my (@ok,@fragi); my($ndrop)= (0);
    for (my $i=0; $i<=$nt; $i++) {
      my($pd,$td,$xn,$icv,$icstat,$aspans)= @{ $tval[$i] };       
      if($aspans =~ /^frag/) { push @fragi, $i; } else { push @ok, $i; }
    }
    if(@ok < @tval) {
      $ndrop += @tval - @ok;
      @fragval= @tval[@fragi]; # preserve output
      @tval= @tval[@ok]; @xn= @xn[@ok];
      $nt= @tval - 1; 
      warn "#drop_frag $gn n=$ndrop\n" if $debug; # seems ok
    }
  }  
  
  if($USE_ICLOCUS and $nt > 0) { 
    my(%igene);
    for (my $i=0; $i<=$nt; $i++) { # all tval here, not just nm top set
      my($pd,$td,$xn,$icv,$icstat,$aspans)= @{ $tval[$i] };
      my($ictype,$icloc,$icsym)=split":",$icv; 
      (my $icg=$icloc) =~ s/\.t\d+//;
      push @{$igene{$icg}}, $i;
    }
    my @icg= sort keys %igene;
    # separate data in @tval, @xn
    if(@icg>1) {
      my @tval_all= @tval; my @xn_all= @xn;
      my @reta;
      for my $icg (@icg) {
        my @i= @{$igene{$icg}};
        @tval= @tval_all[@i]; @xn= @xn_all[@i];
        my @ret= putoneg($gn.'_'.$icg, \@fragval); 
        push @reta, @ret; @fragval=(); # once only, no icloc parse?
      }
      return @reta;
    }
  }
  
  return putoneg($gn, \@fragval);
}

sub putoneg {
  my($gn, $fragset)=@_;
  my $nt= @tval - 1; 
  my $nm= ($nt > MAXTOP) ? MAXTOP : $nt; 
  return 0 if($nt < 1 and $ONLY_ALTS);
    ## FIXME fragof $aspans tag should be skipped for commx .. change to @xn, $nm needed from @tval

  my($ncomm,$commx,$commx_ar,$commsets_arar,$xncomm_hash)= commx($gn, $nm); # @commx = split",",$commx
  my %commx= map{ $_ => 1 } @$commx_ar; # commx ',.,' spacers??
  ## commsets_arar == array of array of @commx, may be 2+ common spans
  my $xncommtot=0; map{ $xncommtot++ if($xncomm_hash->{$_}); } keys %$xncomm_hash;
   
  for (my $i=0; $i<=$nt; $i++) { # all tval here, not just nm top set
  
    my($pd,$td,$xn,$icv,$icstat,$aspans)= @{ $tval[$i] };

    my($ictype,$icloc,$icsym)=split":",$icv; 
    # my($sa,$sc,$ss,$su)=map{ my($c)= $icstat=~m/$_(\d+)/?$1:0; $c; } qw( a c s u);

    ## convert $commx/@commx to icsym patt
    my @xn= split",",$xn; # * should be same num @xn == @sym
    my @sym= split"",$icsym;
    my @csym=(); my $nc=0; 

    my $hascomm= 0; # ($ncomm==0)? -1 : ($xn =~ m/$commx/)? $ncomm."c" : 0;
    if($ncomm==0) {
      $hascomm= ($nt < 1)?"noalts":"nocomm";
    } elsif($xn =~ m/$commx/) {  #?? commx ',.,' spacers
      for(my $j=0; $j<@xn; $j++) { if($commx{$xn[$j]}){ $nc++; push @csym,$sym[$j]; } else { push @csym,'.';} } 
      $hascomm= $nc."c:".join("",@csym); 
    } else { 
      # check subsets, off by start/end, or inserted uniq exons
      my($nok,$nmiss)=(0,0);
      for(my $j=0; $j<@xn; $j++) {
        #SS if($commx[$j] eq '.') {  push @csym,'.'; $nok++; } 
        #SS elsif
        if($commx{$xn[$j]}){ $nc++; push @csym,$sym[$j]; $nok++; } 
        else { push @csym,'.'; $nmiss++ if($nok and $nok<$ncomm); } 
      } 
      if($nc>0) { my $pc=($nmiss>0 or $nok < $ncomm)?"p":"c"; $hascomm="$nc$pc:".join("",@csym); }
      else { $hascomm="nocomm"; }
    }
    
    # another xcomm count, as per  icstat $sc
    my $scb=0; map{ if($xncomm_hash->{$_}>0){ $scb++; } } @xn;
    
    #?? add $nc/$ncomm pct of comm as qual stat
    #?? output "xn:$xn" preserve input tag?
    
    my ($pcxspan)= pctof( ($hascomm eq "noalts" or $ncomm==0)? 0 : ($hascomm eq "nocomm")? 0 : $nc/$ncomm);
    my ($pcexons)= pctof( ($xncommtot>0) ? $scb/$xncommtot: 0);
    
    #o my $commval= $hascomm; # ($hascomm <= 0)?"nocomm":"$commx.$hascomm";
    #p my $commval="$hascomm,cx:$scb/$xncommtot"; #?
    #q my $commval = ($pcexons > $pcxspan)? "$pcexons%cx,$pcxspan%cs," : "$pcxspan%cs,$pcexons%cx,";
    my $commval = "$pcxspan%cs,$pcexons%cx,";
       $commval.= "$hascomm,cx:$scb/$xncommtot"; #?
    print join("\t",$pd,$td,$xn,$icv,$commval,$icstat,$gn)."\n"; # ,$aspans ??** want in output? should out == in but added commx cols?
  }
  
  # ONLY_ALTS above will skip frags:   return 0 if($nt < 1 and $ONLY_ALTS);
  if(ref($fragset) and @$fragset > 0) {
    #o my $commval= "nocommfrag"; # modify? nocommfrag; and nocommnoalt above?
    my $fragt="fragof:"; # fixme
    my $nft= scalar(@$fragset) - 1;
    for (my $i=0; $i<=$nft; $i++) { # all tval here, not just nm top set
      my($pd,$td,$xn,$icv,$icstat,$aspans)= @{ $fragset->[$i] };
      # ($fragt)= $aspans =~ m/(frag\w+:\w+)/;
      my @xnf= split",",$xn; my $scb=0; map{ if($xncomm_hash->{$_}>0){ $scb++; } } @xnf;
      my($pcexons)= pctof( ($xncommtot>0) ? $scb/$xncommtot: 0);
      my $commval= "$pcexons%cx,0%cs,nocommfrag";
      print join("\t",$pd,$td,$xn,$icv,$commval,$icstat,$fragt.$gn)."\n"; # ,$aspans ??** want in output? should out == in but added commx cols?
    }
  }
  return($nt); # + nfragt?
}


sub pctof{ my @p; for (@_) { my $p=int(0.5 + 100*$_); push @p, ($p>100)?100:$p; } return @p; }
 
=item input exontab

 egrep 'Hsref9Evm000001t[2345] ' try9refhum/tmpfiles/ref_human18nc.pf50u20b.exon9j
 
Hsref9Evm000001t2	# pubid
  humang07273t8	  # trid/oid
  # exon chain, local-id string
  xn:g1t1x1,g1t1x1395,g1t1x1666,g1t1x1801,g1t1x10303,g1t1x11254,g1t1x11311,g1t1x14091,g1t1x14371,g1t2x26621u,g1t1x31847,g1t1x32095,g1t1x32886,g1t1x33173,g1t1x33739,g1t1x33996,g1t1x34299,g1t1x34532,g1t1x34610,g1t1x34787,g1t1x35472,g1t1x35546,g1t2x34762u,g1t1x35875,g1t1x35958,g1t2x35927u,g1t1x37281,g1t1x37457,g1t2x36603u,g1t1x39041,g1t1x39123,g1t1x39902,g1t1x40294,g1t1x40474,g1t1x40876,g1t1x40926,g1t1x107976	
  # chain class:locus.altnum:exonchars
  icalt:l2.t2:RfpiTG_fxD_et_bJwnB_r_luZSQLWI_YFE_Ny	
  a37c13s20u4	# exon type counts: all,common,shared,uniq
  # gene-alt seq align break points, 1st is self-span
  1-106869;1-1398,1399-10305,10306-34524,34525-34850,34851-106869;1-1398,1399-33351,36174-36350,37934-38795,39187-106869;1-33048,33662-33839,38016-106869;1-33047,33584-33840,38016-106869;1-10363,13143-33351,36174-106869;1-1398,1399-10363,13143-33351,36174-106869;1-32225,32791-33048,33662-33839,38016-106869;1-1398,1399-10363,13143-33351,36174-36350,37934-38795,39187-106869;1-31147,31938-32225,32791-33048,33662-33840,39367-39769,39819-106869;1-30899,32791-33048,33662-33840,39367-106869;1-30898,31936-32225,32791-33048,33662-33840,39367-106869;1-1398,1399-1669,1804-10306,10307-13423;1-1398,1399-1669,1804-10306,10307-13423;1-13423;1-13423,39819-106869;1-1398,1399-1669,1804-13423;36176-106869;

Hsref9Evm000001t2	xn:g1t1x1,g1t1x1395,g1t1x1666,g1t1x1801,g1t1x10303,g1t1x11254,g1t1x11311,g1t1x14091,g1t1x14371,g1t2x26621u,g1t1x31847,g1t1x32095,g1t1x32886,g1t1x33173,g1t1x33739,g1t1x33996,g1t1x34299,g1t1x34532,g1t1x34610,g1t1x34787,g1t1x35472,g1t1x35546,g1t2x34762u,g1t1x35875,g1t1x35958,g1t2x35927u,g1t1x37281,g1t1x37457,g1t2x36603u,g1t1x39041,g1t1x39123,g1t1x39902,g1t1x40294,g1t1x40474,g1t1x40876,g1t1x40926,g1t1x107976
  icalt:l2.t2:RfpiTG_fxD_et_bJwnB_r_luZSQLWI_YFE_Ny	a37c13s20u4
Hsref9Evm000001t3	xn:g1t1x1,g1t1x1395,g1t1x1666,g1t1x1801,g1t1x10303,g1t1x11254,g1t1x11311,g1t1x14091,g1t1x14371,g1t1x31847,g1t1x32095,g1t1x32886,g1t1x33173,g1t1x33739,g1t1x33996,g1t1x34299,g1t1x34532,g1t1x34610,g1t1x34787,g1t3x34174,g1t3x34259,g1t1x37281,g1t3x34361,g1t1x37457,g1t1x39041,g1t3x34526,g1t1x39123,g1t1x39902,g1t1x40294,g1t1x40474,g1t1x40876,g1t1x40926,g1t1x107976
  icalt:l2.t3:RfpiTG_fx_et_bJwnB_UvQcLIo_YFE_Ny	a33c13s20
Hsref9Evm000001t4	xn:g1t1x1,g1t1x1395,g1t1x1666,g1t1x1801,g1t1x10303,g1t1x11254,g1t1x11311,g1t1x14091,g1t1x14371,g1t4x24395u,g1t1x31847,g1t1x32095,g1t1x32886,g1t4x32225,g1t4x32791,g1t1x33996,g1t1x34610,g1t3x34174,g1t1x34787,g1t5x33290,g1t3x34361,g1t1x39123,g1t1x39902,g1t1x40294,g1t4x34915,g1t1x40876,g1t1x40926,g1t1x107976
  icalt:l2.t7:RfpiTG_fxm_etVjJBU_Ac_YFP_Ny	a28c12s15u1
Hsref9Evm000001t5	xn:g1t1x1,g1t1x1395,g1t1x1666,g1t1x1801,g1t1x10303,g1t1x11254,g1t1x11311,g1t1x14091,g1t1x14371,g1t5x24393u,g1t1x31847,g1t1x32095,g1t1x32886,g1t1x33173,g1t1x33739,g1t1x33996,g1t1x34532,g1t1x34787,g1t3x34259,g1t3x34361,g1t5x33521u,g1t1x39123,g1t4x33485,g1t5x33979u,g1t1x39902,g1t1x40294,g1t1x40474,g1t1x40876,g1t1x40926,g1t1x107976
  icalt:l2.t6:RfpiTG_fxk_et_bJn_vcq_hzYFE_Ny	a30c13s13u4

=cut
