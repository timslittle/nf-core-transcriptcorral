#!/usr/bin/env perl
# evigene/genoasm/blgapend2tab.pl  

=item usage

perl blgapend2tab.pl \
 dmag20sk6tmaca83dc.max10kgap.agp \
 dmag20sk6tmaca83dc.max10kgapends.sk6tfin.blastn \
 > dmag20sk6tmaca83dc.max10kgapends.sk6tfin.degaptab

degaptab:
GapID	gapW	fillW	okfill	fillID	fillB	fillE	fillO	gapAGP
gap1_dmag20sk6tmaca83v_scf29	513	488	ok	dmag20sk6tmafin_ct025199	102990	103478	f	dmag20sk6tmaca83v_scf29:1889:2401:N:513
gap2_dmag20sk6tmaca83v_scf29	3811	4111	ok	dmag20sk6tmafin_ct025199	104703	108814	f	dmag20sk6tmaca83v_scf29:3633:7443:N:3811
gap3_dmag20sk6tmaca83v_scf29	531	3205	poor	dmag20sk6tmafin_ct023615	129381	132586	f	dmag20sk6tmaca83v_scf29:9159:9689:N:531

=item FIXME block dupl fill spans 

  e.g. 2 gaps in hi dupx TE span align to same large fill span

egrep  '(_scf46647|_scf46612)     ' evdgap_dmag20sk6tmaca83dc.degapalltab
gap1_dmag20sk6tmaca83v_scf46612	1911	158716	poor	dmag21ni7ma2call_sc84	196385	355101	f	dmag20sk6tmaca83v_scf46612:388:2298:N:1911
gap1_dmag20sk6tmaca83v_scf46647	1064	196028	poor	dmag21ni7ma2call_sc84	159007	355035	r	dmag20sk6tmaca83v_scf46647:406:1469:N:1064
gap2_dmag20sk6tmaca83v_scf46647	1003	193170	poor	dmag21ni7ma2call_sc84	161185	354355	r	dmag20sk6tmaca83v_scf46647:1838:2840:N:1003
                                ^^^ -- ^^^^ hi ovdup gap expanded by 100x-200x, somewhat bad as 2 gaps align to same fill span.
                                            should replace both gap1,2,interseq w/ 1 fill span
gap2_dmag20sk6tmaca83v_scf46647	1003	-280	bad	dmag20sk6tmafin_ct000519	655	375	r	dmag20sk6tmac 

  * FIXME self-overlap gap fill spans are bad also
gap4_dmag20sk6tmaca83v_scf244057	1032	2143551	poor	c5_dmag20sk6tmaca83v_scf244057	201	2143752	r	dmag20sk6tmaca83v_scf244057:1076245:1077276:N:1032
gap3_dmag20sk6tmaca83v_scf234029	1205	1945742	poor	c3_dmag20sk6tmaca83v_scf234029	201	1945943	r	dmag20sk6tmaca83v_scf234029:4369765:4370969:N:1205
  
=cut

use strict;

my $CTORD=$ENV{ctord} || $ENV{ctorder} || ""; # ctga,ctgb,ctgc ..
my @cord= split/[,\s]/,$CTORD;
my($id,$ofid,$ofb,$ofe,$clen,$gapw,$gid,$didhdr,$lido)=(0) x 9;
my(@bl,%fillfrom,%agp);

sub _sortgv { 
 my $ak= ($$a[3] eq 'bad')?3:($$a[3] eq 'poor')?2:1;
 my $bk= ($$b[3] eq 'bad')?3:($$b[3] eq 'poor')?2:1;
 my $act= $$a[4]; my $bct= $$b[4];
 my $ao=9; if($CTORD){ for my $i (0..$#cord) { my $co=$cord[$i]; if($act =~ m/$co/) { $ao=$i; last; } } }
 my $bo=9; if($CTORD){ for my $i (0..$#cord) { my $co=$cord[$i]; if($bct =~ m/$co/) { $bo=$i; last; } } }
 return ($ak <=> $bk or $ao <=> $bo or $act cmp $bct); 
}
      
# $dofillsort= $ENV{fillagp} or $ENV{tabsort};
if($ENV{fillagp} or $ENV{tabsort}) {

# cat dmag20sk6tmaca83dc.max10kgapends.{selfctg,sk6tfin,ni7ma}.degaptab | egrep -v '#|^==|^GapID' | \
#  sed 's/^gap//; s/_/	/;' | sort -k2,2 -k1,1n -k6,6 | sed 's/^/gap/; s/	/_/; ' \
#  > dmag20sk6tmaca83dc.max10kgapends.degap3tab
  my($nctg,$ngap,$nfill,$tabhdr)= (0) x 9; 
  my (@tabhdr,%arow,%grow);
  while(<>) {
    my @v=split;
    if(/^GapID/) { $tabhdr=$_; @tabhdr=@v; next; }
    elsif(/^\W/){ next; }
    # in-agp: 
    my($gsc,$gi,$i);
    # W: dmag20sk6tmaca83v_scf10096	5519	5636	3	W	c2_dmag20sk6tmaca83v_scf10096	1	118	+
    if($v[4] eq "W"){ ($gsc,$i)= @v[0,3]; $gi= 1+int($i/2); $arow{$gsc}{$gi}{W}=[@v]; $nctg++; } # ($gid=$v[5]) =~ s/^c/gap/; 
    # N: dmag20sk6tmaca83v_scf10096	1328	5518	2	N	4191	scaffold	yes	paired-ends
    elsif($v[4] eq "N") { ($gsc,$i)=@v[0,3]; $gi=int($i/2); $arow{$gsc}{$gi}{N}=[@v]; $ngap++; } # $agv=join(":",@v[0,1,2,4,5]); $agp{$gid}=$agv;  }  
    elsif($v[0] =~ /^gap/) { # or $v[3] =~ m/ok|poor|bad|^[a-z]/
    my($gid,$gapw,$fillw,$okfill,$fillid,$fillb,$fille,$fillo,$gapagp)=@v;
    ($gi,$gsc)=split"_",$gid,2; $gi=~s/^gap//;
    push @{$grow{$gsc}{$gi}}, [@v]; $nfill++;
    }
  }
  
  # $ENV{fillagp}:
  if($ngap and $nctg) { 
  print "# Fill.agp ctorder=$CTORD \n";
  for my $gsc (sort keys %arow) {
    for my $gi (sort keys %{$arow{$gsc}} ) {
      my $wagp= $arow{$gsc}{$gi}{W};
      print join("\t",@$wagp)."\n"; 
      my $nagp= $arow{$gsc}{$gi}{N}; # miss?
      #bad: if($nagp and my @gv= @{ $grow{$gsc}{$gi} } ) 
      if($nagp and $grow{$gsc}{$gi}) {
        my @gv= @{ $grow{$gsc}{$gi} };
        my($svt)= sort _sortgv @gv;
        # F: dmag20sk6tmaca83v_scf10096	1328	5518	2	F	fillctg_scf10096	1	118	+
        my @fagp= @$nagp; $fagp[4]= "F"; # Fill not W
        my($fok,$fid,$fb,$fe,$fo)= @{$svt}[3,4,5,6,7]; 
        if($fok ne 'bad') {
          if($fo eq "r" or $fo eq '-') { $fo='-'; } # do we rev fb,fe?
          elsif($fo eq "f") { $fo='+'; }
          @fagp[5,6,7,8]= ($fid,$fb,$fe,$fo); 
          print join("\t",@fagp)."\n"; 
          print "#o:".join("\t",@$nagp)."\n"; $nagp=0;
          }
        }
      if($nagp) {  print join("\t",@$nagp)."\n"; }
    }
  }
  }
  
  if( $ENV{tabsort} ) {
  print "# Sorted.degaptab \n$tabhdr";
  for my $gsc (sort keys %grow) {
    for my $gi (sort keys %{$grow{$gsc}}) {
      my @gv= @{ $grow{$gsc}{$gi} };
      for my $sv (sort _sortgv @gv) { print join("\t",@$sv)."\n"; } # add order num??
    }
  }
  }
  
  exit;
}

while(<>) {
  if(/^\W/){ 
    if(/# Query: (\S+)/){ my $d=$1; 
      putbl($id,$gapw,@bl) if($id and @bl);
      # fixme: change hdr, no gapw=, use agp N gap len
# Query: gap1_dmag20sk6tmaca83v_scf10096 len=4509; orig=dmag20sk6tmaca83v_scf10096:1127-5636;
# Query: gap2_dmag20sk6tmaca83v_scf10096 len=1981; orig=dmag20sk6tmaca83v_scf10096:5518-7499;
      
      @bl=();  $id=$d; ($ofid,$ofb,$ofe,$clen,$gapw)= (0) x 9;
      if( my($ofs)=m/ofs=([^;\s]+)/) { ($ofid,$ofb,$ofe)=split":",$ofs; }
      ($clen)=m/len=(\d+)/; if(($gapw)=m/gapw=(\d+)/){ $gapw--; }
      } 
    next; 
  }

  my @v=split; 
  if(@v<10){ # in=agp
    # W: dmag20sk6tmaca83v_scf10096	5519	5636	3	W	c2_dmag20sk6tmaca83v_scf10096	1	118	+
    if($v[4] eq "W"){ ($gid=$v[5]) =~ s/^c/gap/; my @cv=@v; }  
    # N: dmag20sk6tmaca83v_scf10096	1328	5518	2	N	4191	scaffold	yes	paired-ends
    elsif($v[4] eq "N") { my $agv=join(":",@v[0,1,2,4,5]); $agp{$gid}=$agv;  }  
    next; 
  } else {
    push @bl, [@v];
    # ($ct,$sc,$pi,$al,$mi,$nd,$cb,$ce,$sb,$se,$ev,$bs,$sw)=@v; 
    # $lct=$ct; $lsc=$sc; 
  }
}

putbl($id,$gapw,@bl) if($id and @bl); 

sub max{ $_[0]>$_[1] ? $_[0] : $_[1]; } 
sub min{ $_[0]<$_[1] ? $_[0] : $_[1] ;}

=item okfill

** add qual to degaptab: ok1, ok2 (2nd fill source ~same as ok1), 
   bad (eg. -117 fillw): bad/poor for -fill or tiny +fill, or huge +fill 
   bad huge is rel to gapw, eg fillw/gapw >= 2? >= 5 .. depends on gapw size, ie 900f/300g ok, but 90000f/30000g not ok?
   .. gap classes: tiny < 100, small 100..300?, medium: 300..3k/5k, large 5k-15k, huge > 15k
=cut

use constant rFillTooHuge => 50; #? 100; # fillw/gapw,  lower?
sub okfill { # decide if align is valid fill
  my($gapw,$fillw,)= @_;
  my $ok="ok"; 
  if($fillw <= 0 or $gapw <= 0) { 
    $ok="bad"; 
  } else {
    my $pfg= $fillw/$gapw;  # my $dfg= $fillw - $gapw;
    #^^ change to bad when fillw >>> gapw, eg. >= 100x
    if($pfg >= rFillTooHuge) { 
      $ok="bad";  # for all gap sizes?
    } elsif($gapw < 100) { # tiny
      if($pfg >= 4 or $pfg < 0.10){ $ok= "poor"; }
    } elsif($gapw < 300) { # small
      if($pfg >= 5 or $pfg < 0.10){ $ok= "poor"; }
    } elsif($gapw < 5000) { # medium
      if($pfg >= 5 or $pfg < 0.30){ $ok= "poor"; }
    } elsif($gapw > 15_000) { # huge, gapw likely very uncertain
      if($pfg >= 5 or $pfg < 0.20){ $ok= "poor"; }
    } else { # large, 5k..15k
      if($pfg >= 5 or $pfg < 0.30){ $ok= "poor"; }
    }
  }
  return $ok;
}

sub putbl {
  my($id,$gapw,@bl)=@_;  
  my %sc; my @sc=grep/\w/, map{ my $s=$_->[1]; ($sc{$s}++)?"":$s; } @bl; 
  my @sbl= sort{ $$a[1] cmp $$b[1] or $$a[6]<=>$$b[6] } @bl;
  for my $s (@sc) { 
    my @tbl=grep { $_->[1] eq $s } @sbl; my $i=0; my $j=$i+1; 
    my $ok=(@tbl>1 and $tbl[$j]->[6] > 100+$tbl[$i]->[7]); 
    if($ok){ 
      my $si=$tbl[$i]; my $sj=$tbl[$j];
      my @sil= @{$si}[1,8,9,6,7]; my @sjl= @{$sj}[1,8,9,6,7]; my $rev=($sil[2] > $sjl[1])?"r":"";
      my $ciw=1 + $sil[4] - $sil[3]; my $cjw=1 + $sjl[4] - $sjl[3];
      (my $ido=$id)=~s/gap\d+.//;  # ido should == ofid
      my $ctgap= ($ofid and $ofb) ? join":",$ofid,$ofb+$sil[4] + 1,$ofe - $cjw - 1 : "";  
      if($rev){ my @t=@sil; @sil=@sjl; @sjl=@t;} 
      my $gb= max(@sil[1,2]); my $ge= min(@sjl[1,2]); $rev||="f";
      my $gsc=$sil[0];  $gb++; $ge--; my $sgapw=$ge-$gb; 
      my $agp=$agp{$id}||""; 
      if($agp){ my @agp=split":",$agp; $gapw ||= $agp[4]; }
      my $okf= okfill($gapw, $sgapw);
      # if($gsc eq $ido) { $okf .= "sameid"; }  # bad? skip? for self blast
      if($gsc =~ m/$ido/) { $okf .= "sameid"; $okf =~ s/poor/bad/; }  # bad? skip? for self blast
      # bad: gap3_dmag20sk6tmaca83v_scf234029 over c3_dmag20sk6tmaca83v_scf234029
      # check dup gsc:gb:ge spans
      if($okf =~ /bad/) { 
        #ignore fillfrom
      } elsif($fillfrom{$gsc}) {
        my @fgb= sort{$a<=>$b} keys %{$fillfrom{$gsc}}; 
        my ($gbo,$geo)=(0,0);
        for my $fb (@fgb) { #? last if($fb > $ge); 
          my $fe= $fillfrom{$gsc}{$fb}; 
          if($fe > $gb) { ($gbo,$geo)=($fb,$fe); last; } 
        }
        if($geo > $gb and $gbo < $ge) { $okf="badover"; }
        if($gb ne $gbo){ $fillfrom{$gsc}{$gb}= $ge; }
      } else {
        $fillfrom{$gsc}{$gb}=$ge;
      }
      my $hdr= join("\t",qw(GapID gapW fillW okfill fillID fillB fillE fillO gapAGP));
      my $val= join("\t",$id, $gapw, $sgapw, $okf, $gsc, $gb, $ge, $rev); # only this + agp?
      $val.= "\t$agp" if($agp);
      $val.= "\t$ctgap" if(0 and $ctgap);
      print "# $ido ...\n" if($lido ne $ido); $lido=$ido;
      print $hdr,"\n" if(1 > $didhdr++); 
      print $val,"\n"; 
      return $val unless($okf =~ /bad/); #? both
      } 
    } 
  return 0; 
} 

=item outputs

dmag20sk6tmaca83dc.max10kgapends.sk6tfin.degaptab
# dmag20sk6tmaca83v_scf29 ...
GapID	gapW	fillW	okfill	fillID	fillB	fillE	fillO	gapAGP
gap1_dmag20sk6tmaca83v_scf29	513	488	ok	dmag20sk6tmafin_ct025199	102990	103478	f	dmag20sk6tmaca83v_scf29:1889:2401:N:513
gap2_dmag20sk6tmaca83v_scf29	3811	4111	ok	dmag20sk6tmafin_ct025199	104703	108814	f	dmag20sk6tmaca83v_scf29:3633:7443:N:3811
gap3_dmag20sk6tmaca83v_scf29	531	3205	poor	dmag20sk6tmafin_ct023615	129381	132586	f	dmag20sk6tmaca83v_scf29:9159:9689:N:531
gap4_dmag20sk6tmaca83v_scf29	153	2370	poor	dmag20sk6tmafin_ct023615	130674	133044	f	dmag20sk6tmaca83v_scf29:12022:12174:N:153
gap5_dmag20sk6tmaca83v_scf29	938	2304	ok	dmag20sk6tmafin_ct023615	134555	136859	f	dmag20sk6tmaca83v_scf29:13689:14626:N:938
gap6_dmag20sk6tmaca83v_scf29	8176	8929	ok	dmag20sk6tmafin_ct025199	119570	128499	f	dmag20sk6tmaca83v_scf29:17402:25577:N:8176
# dmag20sk6tmaca83v_scf48 ...
gap1_dmag20sk6tmaca83v_scf48	9386	9285	ok	dmag20sk6tmafin_ct023245	449	9734	f	dmag20sk6tmaca83v_scf48:1231:10616:N:9386

dmag20sk6tmaca83dc.max10kgapends.degap3tab
gap1_dmag20sk6tmaca83v_scf1004	3869	5895	ok	c1_dmag20sk6tmaca83v_ctg57992371	6696	12591	f	dmag20sk6tmaca83v_scf1004:1500:5368:N:3869
gap1_dmag20sk6tmaca83v_scf1004	3869	5895	ok	dmag20sk6tmafin_ct022421	20838	26733	f	dmag20sk6tmaca83v_scf1004:1500:5368:N:3869
gap1_dmag20sk6tmaca83v_scf1004	3869	5907	ok	dmag21ni7ma2call_sc93	704172	710079	f	dmag20sk6tmaca83v_scf1004:1500:5368:N:3869


=cut