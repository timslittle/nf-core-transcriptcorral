#!/usr/bin/env perl
# evigene/scripts/prot/aanstat.pl
# see evigene/scripts/prot/aastat.sh
# add count aa>=100 (or opt)

use strict;

my $top=$ENV{top}||1000; 
my $main=$ENV{main}||0;
my $smallaa=$ENV{small}||100; 
my $NOTALL=$ENV{notall}||0;
my $SHOWHEAD=$ENV{header}||0;
my $STDN50=$ENV{topn50}?0:$NOTALL; # change for $NOTALL ?
my $OUTH=*STDOUT;

my $GTOP=($main eq "t" or $main =~ /long/)?1:0;
my(%gdid);

sub notmain {
  my($id)=@_;
  if($GTOP) { $id=~s/utrorf//; $id=~s/[t\.]\d+$//; return($gdid{$id}++);  } 
  elsif($main) { return($id=~/$main$/)?0:1; } 
  else { return 0; }
}

my $xh=""; 
$xh.=" and all" unless($NOTALL);
$xh.=" for main-isoform of genes" if($main);
print $OUTH  "# AA-quality summary for $top longest-aa$xh, nbig>=$smallaa.aa\n";
for my $aasize (@ARGV) { 

  my($name,$istemp)=(0) x 9; 
  if($aasize =~ /^stdin|^-$/i) {
    my $aasizetmp="/tmp/aastat$$.in"; # fixme better tempfile() loc
    open(AA,'>',$aasizetmp); 
    while(<STDIN>){ print AA $_; } close(AA);
    $aasize= $aasizetmp; $name=$ENV{name}||"stdin"; $istemp=1;
  } else {
    $name=`basename $aasize .aa.qual`; chomp($name); $name =~ s/\.aa//; 
  }
  unless(-s $aasize) { warn "error: empty aasize=$aasize\n"; next;}
  
  our(@aw,@nn,@aq,@trlen); @aw=@nn=@aq=@trlen=(); 
  my($nok,$n,$sw,$sn,$nt,$nc,$nbig)=(0) x 9; 
  open(AA,"sort -k2,2nr $aasize |"); # FIXME: sort TMPDIR=localdir
  while(<AA>) { 
    next unless(/^\w/); 
    my @v= split; $nt++;
    my($id,$aw,$nn,$aqual,$trlen,$cdoff,$oid)= @v;
    # if($idclass{$id} =~ /okay/) {}
    my $ok=($main and notmain($id)) ? 0 : 1;
    if($ok) { 
      $nok++; push @aw,$aw; push @nn,$nn; push @aq,$aqual; push @trlen,$trlen;
      $nbig++ if($aw >= $smallaa);
    } 
  } 
  close(AA); unlink($aasize) if($istemp);
  
  # stat top1k and allokay
  sub aastat { 
    my ($name,$n,$ntot,$nbig)= @_; $n=@aw if($n > @aw); 
    my($sw,$sn,$nc,$msum,$m50)=(0) x 9; 
    my $n1=$n-1;
    for my $i (0..$n1) { 
      $sw+=$aw[$i]; $sn+=$nn[$i]; $nc++ if($aq[$i]=~/complete/);
      # $msum+=$tlens[$i]; 
      }
    my $aw=int($sw/$n); my $an=int(10*$sn/$n)/10;  my($mx,$md,$mi)= @aw[0,int($n/2),$n1];

    # trlen.n50 here, calc on *all* @trlen by default (ie std calc), but opt to use same n as aatop1k
    ## @trlen not sorted same as @aw aalen
    my $trn1=($STDN50)? @trlen - 1 : $n-1;
    my @tlens= sort{ $b <=> $a } @trlen[0..$trn1];
    my $m50stat="";
    if( $tlens[0] > 0 ) { # have trlen data?
      for my $i (0..$trn1) { $msum+= $tlens[$i]; }
      my $mhalf=$msum/2;  my $mrunt=0; 
      for my $i (0..$trn1) { my $tl= $tlens[$i]; $mrunt+=$tl; if($mrunt >= $mhalf){ $m50=$tl; last; } } # $m50="$tl,i$i"; #debug
      $m50stat=" trn50=$m50;";
    }
    my $ntv=($NOTALL) ? "nt=$ntot;" : "n=$n;";
    print $OUTH  "$name\t $ntv average=$aw; median=$md; min,max=$mi,$mx; nfull=$nc; nbig=$nbig;$m50stat gaps=$sn,$an\n"; #sum=$sw;     
    # as per aastat.sh, add nfull=$nc
    # add nsmall=nnn, nbig=mmm : both?  
  }
  
  print $OUTH  "# AA-quality summary of $aasize: longest $top and all\n" if($SHOWHEAD);
  aastat("$name.top",($nok<$top)?$nok:$top, $nok,$nbig);
  aastat("$name.all",$nok,$nok,$nbig) unless($NOTALL or $nok<$top);

}

__END__

=item usage

# old w/o n50
perl -ne 'BEGIN{ @I= map{ $_ - 1 } (2,3,6,7,9,10); } 
if(/^#tsa.(\w+)/){ $tc=$1; if($tab=$tab{$tc}) { ($nt,$ave,$nf)= m/(?:nt|average|nfull)=(\d+)/g; 
($amax)=m/max=\d+,(\d+)/; print "$tab\t$ave.a1k,$amax.am\t$nt\n";  } } 
else { if(/^\s+code/) { s/^\s+/item\t/; }  @v=split; ($tc)=@vv=@v[@I]; $tab{$tc}=join"\t",@vv; 
if($tc eq "code") { print $tab{$tc},"\taaqual\tntr\n"; } }' \
  $bg/evigenes/transraterr/transrate_trqual_gr16s2t.txt tsaorf_G*.aastat \
  > $bg/evigenes/transraterr/transrate_trqual_gr16_aastat.tab

# new tabulate n50 w/ aanstat.pl, for nt=155 items

/bio/bio-grid/evigenes/trqual_transratef/
env notall=2 $evigene/scripts/prot/aanstat.pl wwqual/*.aa.qual \
   > tsaorf_trtab_nu.aastat

perl -ne 'BEGIN{ @I= map{ $_ - 1 } (2,3,6,7,9,10); } 
if(/^tsa.(\w+)/){ $tc=$1; if($tab=$tab{$tc}) { ($nt,$ave,$nf,$n50)= m/(?:nt|average|nfull|trn50)=(\d+)/g; 
($amax)=m/max=\d+,(\d+)/; print "$tab\t$ave.a1k,$amax.am\t$n50\t$nt\n";  } } 
else { if(/^\s+code/) { s/^\s+/item\t/; }  @v=split; ($tc)=@vv=@v[@I]; $tab{$tc}=join"\t",@vv; 
if($tc eq "code") { print $tab{$tc},"\taaqual\tn50\tntr\n"; } }' \
  $bg/evigenes/transraterr/transrate_trqual_gr16s2t.txt \
  tsaorf_trtab_nu.aastat \
  > transrate_trqual_gr16_aastatnu.tab

=cut
