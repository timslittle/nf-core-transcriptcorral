#!/usr/bin/env perl
# gnodes_sumasmgcn2rtab.pl

=item usage

opts:  env awt=1 gwt=1 missr=1 neg=1  .. better

head -25 *gnode*/{*test8f_SRR*_sum.txt,*test8f7sumgcn_genesum.txt} | \
 env awt=1 gwt=1 missr=1 neg=1  $gw/daphnia/gnodes_genecopyf/gnodes_sumasmgcn2rtab.pl \
   > gnodes_asmgcn_miss3e.rtab

.. env rplot=1 ..  gnodes_sumasmgcn2rtab.pl > gnodes_asmgcn_miss3e.R

=item update
   
   21SEP27: upd for covtab8j, sumgenecov upd w/ lCopy cols

=cut 

#? default on: neg=1 missr=1 awt=1 gwt=1

unless($ENV{nodef}) {
$DONEG= $DOMISR= $DO_AWEIGHT= $DO_GWEIGHT= 1;
} else {
$DONEG=$ENV{neg}||0; $DOMISR=$ENV{missr}||0;
$DO_AWEIGHT= $ENV{awt}||$ENV{aweight}||0;
$DO_GWEIGHT= $ENV{gwt}||$ENV{gweight}||0;
}
$DOGL=$ENV{glev}||0; 
$DROPCLASS=$ENV{dropclass}||"CAll"; # not for R plots.

$DO_SHOWS = ($ENV{nosrc})?0:1; $DO_SHOWD = $ENV{showd}||0;
$DO_RSCRIPT=$ENV{rplot}||0;
$YMIN=$ENV{ymin}||-40; #  was -60 for dapmag min ?? maybe calc from all input data range?
$YMIN=-$YMIN if($YMIN>0);
my($pxlo_min,$pxlo_max)= (9999,-9999);

my @asum= qw(allasm uniqasm dupasm cdsasm); # opt?
# allasm uniqasm dupasm cdsasm CDSann CDSbus TEann RPTann NOann
if($ENV{aann}){
  @asum= qw(allasm uniqasm dupasm CDSann TEann RPTann); # opt?
}
$asumpatt=join"|",@asum;


my $outh= *STDOUT; # change?
if($DO_RSCRIPT) {
  rplot($outh, "", 1 > $nRplot++); # header 1st in case DO_SHOWD
}

sub fnsource {
  # local( $_ )= @_; 
  my($fn)= @_;
  $fnsource= $fn; 
  @fn=split"/",$fn; if(@fn>2){ ($fg)=grep/20gnode/,@fn; $fb=$fn[-1]; $fn="$fg/$fb" if($fg); } # fixme: dmag20gnodes/testxxx/dataset
  unless( $fn =~ s,20gnode./,:, ) { @fn=split("/",$fn); $fn=$fn[-1]; }
  map { 
    s/:arath18tair_chr/:at18chr/; s/:arath/:at/; 
    s/at18chr_(he|bwa)/at21${1}_/; s/aweed:at/arath:/; 
    s/:drosmel6/:/; s/:dropse/:/; 
    s/:daphpulex_pa42v2/:19ml/; s/:daphplx_gasm16ml/:16ml/; s/:dplx/:/;  
    s/:(dapcar|daphcari)/:/; s/:dmag(.....).*/:$1/; s/dmag:20sk4/dmag:20maca/; 
    s/^d(mag|car|plx)/da$1/;  s/maca.*/maca/; 
  } ($fn); 
  $fn=~s/_.*//; 
  ($fs,$fa)=split":",$fn; 
  unless($fa) { $fs=~s/dro/dr/; $fs=substr($fs,0,10); }
  elsif($fs=~/^dro/){ $fs=~s/dro/dr/; $fs=substr($fs,0,4); $fa=substr($fa,0,6);  }
  else { $fs=substr($fs,0,3); $fa=substr($fa,0,7); }
  $fsa=$fs.$fa;
  $ngall=$datend=0; $nsrc++;
  if($row{$fsa}{'Aall'} and $row{$fsa}{'C01'}) { # fixme : need new fsa name
    $fsa .= ++$fsadup;
  }
  $src{$fsa} .= $fnsource.", ";
  print "#src:$fsa = $fnsource,$nsrc\n" if($DO_SHOWS or $DO_SHOWD);
}

$datend=1;
while(<>) {

  if(/^==> (\S+)/) { fnsource($1);
#     $fnsource= $fn=$1; 
#     unless( $fn =~ s,20gnode./,:, ) { @fn=split("/",$fn); $fn=$fn[-1]; }
#     map{ 
#     s/:arath18tair_chr/:at18chr/; s/:arath/:at/; 
#     s/at18chr_(he|bwa)/at21${1}_/; s/aweed:at/arath:/; 
#     s/:drosmel6/:/; s/:dropse/:/; 
#     s/:daphpulex_pa42v2/:19ml/; s/:daphplx_gasm16ml/:16ml/; s/:dplx/:/;  
#     s/:(dapcar|daphcari)/:/; s/:dmag(.....).*/:$1/; s/dmag:20sk4/dmag:20maca/; 
#     s/^d(mag|car|plx)/da$1/;  s/maca.*/maca/; } ($fn); $fn=~s/_.*//; 
#     ($fs,$fa)=split":",$fn; 
#     unless($fa) { $fs=~s/dro/dr/; $fs=substr($fs,0,10); }
#     elsif($fs=~/^dro/){ $fs=~s/dro/dr/; $fs=substr($fs,0,4); $fa=substr($fa,0,6);  }
#     else { $fs=substr($fs,0,3); $fa=substr($fa,0,7); }
#     $fsa=$fs.$fa;
#     $ngall=$datend=0; $nsrc++;
#     if($row{$fsa}{'Aall'} and $row{$fsa}{'C01'}) { # fixme : need new fsa name
#       $fsa .= ++$fsadup;
#     }
#     $src{$fsa} .= $fnsource.", ";
#     print "#src:$fnsource,$nsrc\n" if($DO_SHOWS or $DO_SHOWD);
    
  } elsif(/^ (All|[12])/) { # CLevel  genecnsum/sumgcn_genesum
    next if($datend);
    my($cl,$ng,$gc,$xc,$ac,$lcp,$prmis,$xlohi)= (0) x 9;
    my @v=split; 
    # FIXME: covtab8j, sumgenecov upd w/ lCopy cols
    # CLevel	nGene	gCopy	xCopy	aCopy	lCopy	aMiss	x0,xLo,xEq,xHi
    if(@v > 7){ 
      ($cl,$ng,$gc,$xc,$ac,$lcp,$prmis,$xlohi)= @v;
      map{ s/,.*//; } ($ac,$lcp); # may be tuple: ave,mdn
      # ($ac,$lcp)= ($lcp,$ac); #? yes or no : NO, not yet, lCopy is unreliable for complex dups
    } else {
      ($cl,$ng,$gc,$xc,$ac,$prmis,$xlohi)= @v; $lcp= $ac; # ?
    }

    print "#gcn:$_" if($DO_SHOWD);
    @xlh=split",",$xlohi; shift @xlh if(@xlh>3); # x0 first always now?
    
    if($DO_GWEIGHT) {  # wt xlh by *all ng* ?? maybe need to recalc using expect gcn - obs gcn from sumgcn.sumgenetab
    $ngall=$ng if($cl eq "All"); # *should* be 1st row
    ($pxlo,$pxeq,$pxhi)= map{  ($ngall < 1)? 0 : sprintf"%.2f",100*$_/$ngall; } @xlh;
    
    } else {
    ($pxlo,$pxeq,$pxhi)= map{  sprintf"%.2f",100*$_/$ng; } @xlh;
    }
    ($prmis,$nrmis)=split",",$prmis; $prmis=~s/%//;# $prmis=~s/%,.*//; 
    $nrtot= ($prmis < 0.1 ) ? 0 : int($nrmis / ($prmis/100));
    $prhit= sprintf"%.2f", (100 - $prmis);
    $prmis= sprintf"%.2f",$prmis;   
    
    $cl="0$cl" if($cl<10 and $cl ne "All"); 
    $cgn="C$cl"; $cgv="G$cl"; $cgmis=$cgn."mis";
    if($DONEG){ $prmis= -$prmis; $pxlo= -$pxlo; } 
    $row{$fsa}{$cgn} = [$pxlo,$pxeq,$pxhi,$prmis, ];
    if($DOMISR) { $row{$fsa}{$cgmis} = [$prmis,$prhit,$nrmis,$nrtot, ]; }
    if($DOGL){ $row{$fsa}{$cgv} = [$ng,$gc,$ac ]; }
    $pxlo_min= $pxlo if($pxlo<$pxlo_min); $pxlo_max= $pxlo if($pxlo>$pxlo_max);
    # print join("\t",qw(SppAsm GcnLevel pXlo pXeq pXhi pRMiss))."\n" if(1>$hd++); 
    # print join("\t",$fsa,"C$cl",$pxlo,$pxeq,$pxhi,$prmis)."\n"; 

    
  } elsif(/^Field\s/){ 
    $datend=2;    
   
  } elsif(/^\s*Genes Copynum Summary/) { # genecnsum.txt top
    #  Genes Copynum Summary by Copy level, for nGene=13918
    if($datend>0){ fnsource($ARGV); $datend=0; }
  
  } elsif(/^\s*Source=/ and /KUlow=/) { # chrasm_sum.txt top  
    # Source=Dropse20uc, KUlow=88.7, KUhigh=94, FlowcytSize=161-180 Mb Formula_LN/C=166.8-176.7 Mb (dropse20chrs)
    if($datend>0){ fnsource($ARGV); $datend=0; }
  } 
  
  # option: allasm|uniqasm|dupasm|cdsasm .. try TEann CDSann RPTann NOann
  #o elsif(/^(allasm|uniqasm|dupasm|cdsasm)/) # chrasm_sum.txt
  elsif(/^($asumpatt)/)  #? allasm uniqasm dupasm CDSann TEann RPTann
  { 
    next if($datend);
    print "#asm:$_" if($DO_SHOWD);
    ($casm,$omb)= @v=split"\t"; 
    ($emb,$xcopy)= (@v<7) ? @v[2,3] : @v[4,5];
    
    #o $casm=~s/asm//; $casm="A$casm";
    $casm=~s/asm|ann//; $casm="A".lc($casm);
    
    $pxeq= ($xcopy<0.01)? 0 : sprintf"%.2f", 100/$xcopy;    
    $pxlo=  sprintf"%.2f",100 - $pxeq; $pxlo= -$pxlo if($DONEG);
    $row{$fsa}{$casm}=[$pxlo, $pxeq, $omb, $emb]; #?? ,
        
  } elsif(/^\s*Total span=/) {
    # Total span=163.2 Mb, covspan=163.2, gaps=0 Mb for assembly dropse20chrs
  
  } elsif(/^\s*Genome Size Est=/) {
    # Genome Size Est= 168.4-178.5 Mb (Nread), 166.8-176.7 Mb (Maprd), for readset SRR11460802ab,
    #  for Size=LN/C, Cov=88.7,94, N_reads=105560784, N_maprd=104503221,99.0%, L_readlen=150
  
  } elsif(/^\s*Uniq Conserved Gene Cover/) {
    # Uniq Conserved Gene Cover: median=88.7, ave=89.1, sem=3.62, n=612
    $datend=1;
  }
  
}

# print $outh join("\t",qw(SppAsm Item pXlo pXeq pXhi pMiss ))."\n" if(1>$hd++); 
# print join("\t",qw(SppAsm Item Ngene gCopy aCopy ))."\n" if(1>$hd++); 
# print join("\t",qw(SppAsm Item pXlo ObsM EstM ))."\n" if(1>$hd++); 
for $fsa (sort keys %row) {
  putoneasm($outh, $fsa);
  # for $cl (sort keys %{$row{$fsa}}) { $rv= $row{$fsa}{$cl};  print $outh join("\t",$fsa,$cl,@$rv)."\n"; }
}


sub putoneasm {
  my($outh, $fsa)= @_;
  my @cl= sort keys %{$row{$fsa}};
  @hdr=qw(SppAsm Item pXlo pXeq pXhi pMiss );
  
  if($DO_AWEIGHT){
    # $row{$fsa}{$casm}=[$pxlo, $pxeq, $omb, $emb]; #?? ,
    if($allmb = $row{$fsa}{'Aall'}->[3] ) {
      @arv= map{ $row{$fsa}{$_} } grep /^A/, @cl;
      for $arv (@arv) { ($pxlo,$pxeq,$omb,$emb)= @$arv; 
        $plomb= sprintf"%.2f",100*($emb - $omb) / $allmb; $plomb= -$plomb if($DONEG);
        $arv->[0]= $plomb;
      }
    }
  }
  
  if($DOMISR) {  # sum(nrmis)/sum(nrtot), not all cgmis = C1mis,C2mis,..
    # $row{$fsa}{$cgmis} = [$prmis,$prhit,$nrmis,$nrtot, ]; 
    @arv= map{ $row{$fsa}{$_} } grep { m/^C[1-9].+mis$/ } @cl; # ** SKIP Call,C0 for sum, not valid addition
    my($smis,$stot)=(0,0);
    for $arv (@arv) { ($prmis,$prhit,$nrmis,$nrtot, )= @$arv; $smis+=$nrmis; $stot+=$nrtot; }
    $prmis= ($stot<1)?0: sprintf"%.2f", 100*$smis/$stot; $prhit= 100 - $prmis;
    if($DONEG){ $prmis= -$prmis; }
    $row{$fsa}{'Cmiss'}= [$prmis,$prhit,$smis,$nrtot, ];
    @cla= grep { not m/^C.+mis$/ } @cl;
    push @cla,"Cmiss"; @cl=@cla;
  }
  
  if($DO_GWEIGHT){
  
  }
  
  if($DROPCLASS) { @cl= grep{ not m/$DROPCLASS/ } @cl; }
  # @cl reorder: Aall  Acds  Adup  Auniq > Aall Auniq Adup Acds
  map{ s/Auniq/Ab/; s/Adup/Ac/; s/Acds/Ae/; } @cl; # s/Ate/At/; s/Arpt/Ar/; 
  @cl= sort @cl;
  map{ s/Ab/Auniq/; s/Ac/Adup/; s/Ae/Acds/; } @cl;
  
  if($DO_RSCRIPT) {
    my $rval = join(" ",@hdr)."\n";
    for $cl (@cl) { $rv= $row{$fsa}{$cl}; $rval .= join(" ",$fsa,$cl,@$rv)."\n"; }
    rplot($outh, $rval, 1 > $nRplot++);
  } else {
    print $outh join("\t",@hdr)."\n" if(1>$hdr++); 
    for $cl (@cl) { $rv= $row{$fsa}{$cl}; print $outh join("\t",$fsa,$cl,@$rv)."\n"; }
  }
  
}

=item rplot

orig 2v.. change
----
tab2v <- function(r, tabd) { 
 v= tabd[r,c("Item","pXlo","pMiss")]; v[,1]= as.character(v[,1]); mm= 10 * min(v[,3]); 
 vv=rbind( v[c(1,4,3,2),], c("Cmiss",mm,0), v[5:7,]); vv[,2]= as.numeric(vv[,2]); vv; 
}

for (sp in 1:4) {  
 spn=levels(tabd[,1])[sp]; r=tabd$SppAsm == spn & tabd$Item != "CAll";  v=tab2v(r,tabd); 
 pdf(paste("asmgcxlo1a_",spn,".pdf",sep=""),6,4);  
 barplot( v[,2], ylab="% Deficit", ylim=c(-90,0),
   main=paste("Chrasm/GeneCN for",spn), names.arg=v[,1], cex.names=0.85, space=0.55); 
 mtext("Chr Asm", side=1, line=0, adj=0.25); mtext("Gene CN", side=1,line=0, adj=0.75); 
 dev.off(); 
}

=item rplot upd

vd=read.table(header=T, text="
SppAsm  Item    pXlo    pXeq    pXhi    pMiss
ara20ma Aall    -30.41  69.44   119     171
ara20ma Auniq   -7.02   90.09   98      110
ara20ma Adup    -23.98  33.67   21      62
ara20ma Acds    -11.7   72.46   52      72
ara20ma C01     0       94.64   0.75    -0.9
ara20ma C02-9   -1.57   1.25    0.09    -2.9
ara20ma C10-99  -0.2    0.04    0.01    -4.6
ara20ma Cmiss   -1.68   98.32   598445  18724705
");

acdeficitplot(vd);

=cut

sub rplot {
  my($outh, $datatable, $first)=@_;
  
# plname? gnodes_ chr_gene_deficits_ ?? asmgcn_miss_ ?? test: asmgcxlo3d_
  if($first) {
  
    # plotsize: 6,4 < should be nrow, max(4, abs(YMIN)/10) ?
    my $bbs='\\\\';
    my $rtop=<<"EOS";
#!/usr/bin/env Rscript

acdeficitplot <- function(vd, plname="asmgcn_miss_",ymin=$YMIN, gscale=4) {
  pnames= vd[,2]; rcn= grep("^C",pnames); 
  vd[rcn,3]= gscale*vd[rcn,3]; 
  vcol=rep("gray",nrow(vd)); vcol[rcn]="green"; 
  rm= grep("^Cmiss",pnames); if(length(rm)>0){ vcol[rm]="red2"; }
  pnames= gsub("^A(${bbs}w+)","${bbs}U${bbs}1",pnames,perl=T)
  spn=as.character(vd[1,1]);
  plw= 0.80 * max(7, nrow(vd)); plh= 0.10 * max(40,abs(ymin));
  pdf(paste(plname,spn,".pdf",sep=""),plw,plh);  # 6,4 inch min
  barplot( vd[,3], ylab="% Deficit", ylim=c(ymin,0), col=vcol,
     main=paste("Chrasm/GeneCN for",spn), names.arg=pnames, cex.names=0.85, space=0.55); 
  gticks=seq(0,ymin,-10); axis(4, at=gticks,  labels=gticks/gscale);
  mtext("Chr Asm", side=1, line=0, adj=0.25); mtext("Gene CN", side=1,line=0, adj=0.75); 
  dev.off()
}

EOS
    print $outh $rtop;
  }
  return unless($datatable); # can print top before datatable
  
  my $rplot=<<"EOS";
vd=read.table(header=T, text=\"
$datatable\");
acdeficitplot(vd);
  
EOS
    print $outh $rplot;
}


__END__

=item try2

opts:  env awt=1 gwt=1 missr=1 neg=1  .. better
barplots: need 2 y-axes: Chrpart ylim=(0,-50?) , GeneCN ylim=c(0,-15?)

head -25 dmag20gnodes*/{*test8f_SRR*_sum.txt,*test8f7sumgcn_genesum.txt} | \
 env awt=1 gwt=1 missr=1 neg=1  $gw/daphnia/gnodes_genecopyf/gnodes_sumasmgcn2rtab.pl \
   > gnodes_asmgcn_miss3e.rtab

.. env rplot=1 ..  gnodes_sumasmgcn2rtab.pl > gnodes_asmgcn_miss3e.R

SppAsm  Item    pXlo    pXeq    pXhi    pMiss
dam19sk Aall    -44.44  55.56   115     207
dam19sk Auniq   -10.14  81.30   87      108
dam19sk Adup    -34.78  27.70   28      100
dam19sk Acds    -35.27  41.67   52      125
dam19sk C01     0       74.28   1.94    -0.55
dam19sk C02-9   -11.09  2.53    0.03    -2.5
dam19sk C10-99  -2.79   0.00    0.00    -1.6
dam19sk CAll    -14.17  76.80   8.14    -6.2
dam19sk Cmiss   -3.92   96.08   1160977 15864596

dam20ma Aall    -14.72  85.47   197     231
dam20ma Auniq   0.43    100.00  98      97
dam20ma Adup    -14.72  74.07   99      133
dam20ma Acds    -9.09   77.52   76      97
dam20ma C01     0       67.44   8.98    -0.03
dam20ma C02-9   -7.77   5.57    0.45    -0.01
dam20ma C10-99  -2.72   0.12    0.00    -0.01
dam20ma CAll    -10.84  73.13   15.95   -0.02
dam20ma Cmiss   0       100     6005    0

aweed20gnodes
SppAsm  Item    pXlo    pXeq    pXhi    pMiss
ara18ch Aall    -32.77  67.11   119     177
ara18ch Auniq   -14.12  81.30   108     133
ara18ch Adup    -18.76  25.77   11.6    44.8
ara18ch Acds    -12.99  68.49   50      73
ara18ch C01     0       95.43   0.37    0
ara18ch C02-9   -1.17   1.67    0.11    0
ara18ch C10-99  -0.06   0.15    0.04    -0.01
ara18ch CAll    -1.27   97.26   1.47    0
ara18ch Cmiss   0       100     1259    0

ara20ma Aall    -30.41  69.44   119     171
ara20ma Auniq   -7.02   90.09   98      110
ara20ma Adup    -23.98  33.67   21      62
ara20ma Acds    -11.7   72.46   52      72
ara20ma C01     0       94.64   0.75    -0.9
ara20ma C02-9   -1.57   1.25    0.09    -2.9
ara20ma C10-99  -0.2    0.04    0.01    -4.6
ara20ma CAll    -1.81   95.93   1.78    -1.7
ara20ma Cmiss   -1.68   98.32   598445  18724705

dros
SppAsm  Item    pXlo    pXeq    pXhi    pMiss
drmeref Aall    -16.87  83.33   138     166
drmeref Auniq   -7.83   90.09   115     128
drmeref Adup    -9.04   60.24   23      38
drmeref Acds    -5.42   74.07   28      37
drmeref C01     0       96.47   0.83    0
drmeref C02-9   -0.39   1.00    0.16    0
drmeref C10-99  -0.83   0.13    0.00    0
drmeref CAll    -1.22   97.60   1.18    0
drmeref Cmiss   0       100     126     0

drps20c Aall    -6.32   94.34   163     174
drps20c Auniq   -4.6    94.34   129     137
drps20c Adup    -1.72   92.59   34      37
drps20c Acds    -1.72   90.09   33      36
drps20c C01     0       92.30   0.47    0
drps20c C02-9   -0.68   4.63    0.24    0
drps20c C10-99  -0.57   0.94    0.01    0
drps20c CAll    -1.31   97.87   0.82    0
drps20c Cmiss   0       100     257     0

=cut
  
 ls */{*test8f_SRR*_sum.txt,*test8f7sumgcn_genesum.txt}
aweed20gnodes/arath18tair_chr_bwatest8f_SRR10178325_sum.txt
aweed20gnodes/arath18tair_chr_hetest8f_SRR3703081_sum.txt
aweed20gnodes/arath18tair_chr_test8f_SRR10178325_sum.txt
aweed20gnodes/arath20max_chr_test8f_SRR10178325_sum.txt
aweed20gnodes/at18chr_bwatest8f7sumgcn_genesum.txt
aweed20gnodes/at18chr_hetest8f7sumgcn_genesum.txt
aweed20gnodes/at18chr_test8f7sumgcn_genesum.txt
aweed20gnodes/at20max_test8f7sumgcn_genesum.txt
dcar20gnodes/dapcar20maca1i_fin_test8f7sumgcn_genesum.txt
dcar20gnodes/dapcar20maca1i_fin_test8f_SRR10389283a_sum.txt
dcar20gnodes/daphcari20chr_test8f7sumgcn_genesum.txt
dcar20gnodes/daphcari20chr_test8f_SRR10389283a_sum.txt
dmag20gnodes/dmag14bgi2vtop5k_test8f7sumgcn_genesum.txt
dmag20gnodes/dmag14bgi2vtop5k_test8f_SRR7825549b_sum.txt
dmag20gnodes/dmag15nwb2asm_test8f7sumgcn_genesum.txt
dmag20gnodes/dmag15nwb2asm_test8f_SRR7825549b_sum.txt
dmag20gnodes/dmag19skasm_test8f7sumgcn_genesum.txt
dmag20gnodes/dmag19skasm_test8f_SRR7825549b_sum.txt
dmag20gnodes/dmag20sk4maca20ok_test8f7sumgcn_genesum.txt
dmag20gnodes/dmag20sk4maca20ok_test8f_SRR7825549b_sum.txt
dplx20gnodeq/daphplx_gasm16ml_test8f7sumgcn_genesum.txt@
dplx20gnodeq/daphplx_gasm16ml_test8f_SRR13333791_sum.txt@
dplx20gnodeq/daphpulex_pa42v2_test8f7sumgcn_genesum.txt@
dplx20gnodeq/daphpulex_pa42v2_test8f_SRR13333791_sum.txt@
dplx20gnodeq/dplx20maca4pkr_dc_test8f7sumgcn_genesum.txt@
dplx20gnodeq/dplx20maca4pkr_dc_test8f_SRR13333791_sum.txt@
dromel20gnodes/drosmel6ref_chr_test8f7sumgcn_genesum.txt
dromel20gnodes/drosmel6ref_chr_test8f_SRR11460802ab_sum.txt
dropse20gnodes/dropse20chrs_test8f7sumgcn_genesum.txt
dropse20gnodes/dropse20chrs_test8f_SRR11813283a_sum.txt

 

dgmm:gnodes20f:% head -15  dropse20gnodes/{*_test8f_SRR*_sum.txt,*_test8f7sumgcn_genesum.txt}

==> dropse20gnodes/dropse20chrs_test8f_SRR11813283a_sum.txt <==
Source=Dropse20uc, KUlow=88.7, KUhigh=94, FlowcytSize=161-180 Mb Formula_LN/C=166.8-176.7 Mb (dropse20chrs)
_____   ______  Low     Low     High    High    Total    
Item_   Obs.Mb  Est.Mb  xCopy   Est.Mb  xCopy   tCopy   Description
---------------------------------------------------------------------------
allasm  163     164     1.00    174     1.06    12.18   measured assembly
_total  163     164     .       174     .       .       total assembly 
uniqasm 129     129     1.00    137     1.06    1.06    asm with unique gDNA
dupasm  34      35      1.02    37      1.08    54.04   asm with multimap gDNA
cdsasm  33      34      1.05    36      1.11    2.45    asm with CDS-mapped gDNA
CDSann  27      29      1.05    30      1.12    1.49    asm with CDS annotations
CDSbus  1.7     1.7     1.00    1.8     1.06    1.07    asm with unique BUSCO orlog CDS
TEann   0.2     0.1     0.95    0.2     1.00    40.41   asm with Transposons
RPTann  18.3    18.5    1.01    19.6    1.07    4.74    asm with simple Repeats
NOann   119     118     0.99    125     1.05    15.60   asm without annotations
UCGmed  0.1     0.1     1.01    0.1     1.07    1.07    
---------------------------------------------------------------------------
  xCopy = excess read copy depth; tCopy = assembly multi-map copy depth
  Total span=163.2 Mb, covspan=163.2, gaps=0 Mb for assembly dropse20chrs
  Genome Size Est= 168.4-178.5 Mb (Nread), 166.8-176.7 Mb (Maprd), for readset SRR11460802ab,
   for Size=LN/C, Cov=88.7,94, N_reads=105560784, N_maprd=104503221,99.0%, L_readlen=150
  Uniq Conserved Gene Cover: median=88.7, ave=89.1, sem=3.62, n=612

==> dropse20gnodes/dropse20chrs_test8f7sumgcn_genesum.txt <==
# Chromosome Assembly Gene-Copynum Counts by Copy level
# N Genes valid=14329, zero-cover=1, all=14330 (zero:gClass=zero or gNread<1)
 ---------------------------------------
 Genes Copynum Summary by Copy level, for nGene=14329
 aCopy=Chrasm, gCopy=Gene set, xCopy=aCopy/gCopy, xLo,xEq,xHi= xCopy<1,=1,>1
 aMiss=Missed gene reads on chrasm, % of gene reads
 CLevel nGene   gCopy   xCopy   aCopy   aMiss   xLo,xEq,xHi
 All    14329   1.6     1.0     1.3     0.00063%,130    188,14024,117
 0      14      0.6     2.0     1.1     0%,0    0,0,14
 1      13294   1.0     1.0     1.0     0.00048%,72     0,13226,68
 2-9    795     3.0     1.0     2.9     0.0031%,51      98,663,34
 10-99  217     19.9    0.8     13.2    0.00014%,4      81,135,1
 99-499 9       277.1   0.0     1.7     0.00021%,3      9,0,0
 ---------------------------------------

>> above 2rtab gnodes20sumasmgcn1a.rtab
SppAsm  Item    pXlo    pXeq    pXhi    pMiss
drps20c Aall    -5.66   94.34   0       0
drps20c Acds    -9.91   90.09   0       0
drps20c Adup    -7.41   92.59   0       0
drps20c Auniq   -5.60   94.40   0       0
drps20c C01     0       99.49   0.51    0
drps20c C02-9   -12.33  83.40   4.28    0
drps20c C10-99  -37.33  62.21   0.46    0
drps20c CAll    -1.31   97.87   0.82    0
..
drps20c C02-9   -12.33  83.40   4.28    0
drps20c C02-9   -5.74  xxx.   xx.    0  : % gene total wt = 100 * (1/xc) * ngc / ngt = 100/0.97 * 795/14329
  795/14329 = 5% of total in C02, xCopy deficit = 1/0.97 = 1.030 .. not right, need (1/xc - 1/1)?
  ie, for 5% genes in C2 level, there is an 0.74% deficit in C2 copies, relative to total of genes?
  C2.xplo=0.684% for drps, vs -12.33% old; = 100*xlo/ngt may be best answer
  C2.xplo=1.17% for arath, vs -39.58% old pxlo
  C2.xplo=11.1% for dmag19sk, vs -80.44% old
  C2.xplo=7.77% for dmag20ma, vs -56.33 old
  
   or pxlo-gwt= (2.9/3.0)[ac/gc] * 98,663,34 / 14329[gtot] ??

perl -e '$clab="C2.dpse"; @x=(98,663,34); $xc=2.9/3.0; $ngc=795; $ngt=14329; $difngc=(1/$xc)*$ngc; $pxc=sprintf"%.3g", 100*$difngc /$ngt; @p=map{ sprintf"%.3g", 100*$_/$ngt } @x; $xc=sprintf"%.2f",$xc; print "$clab: ngc=$ngc,ngt=$ngt,xc=$xc, pxc=$pxc, x=@x, xp=@p; \n";'
C2.dpse: ngc=795,ngt=14329,xc=0.97, pxc=5.74, x=98 663 34, xp=0.684 4.63 0.237; 

---
aweed20gnodes/at18chr_test8f7sumgcn_genesum.txt
 CLevel nGene   gCopy   xCopy   aCopy   aMiss   xLo,xEq,xHi
 All    27347   1.3     1.0     1.1     0.0043%,825     347,26597,403
 1      26200   1.0     1.0     1.0     0.00066%,88     0,26098,102
 2-9    806     3.1     0.8     2.3     0.0043%,32      319,458,29
 10-99  68      21.0    1.0     17.5    0.0099%,314     16,41,11

perl -e '$clab="C2.at18"; @x=(319,458,29); $xc=2.3/3.1; $ngc=806; $ngt=27347; $difngc=(1/$xc)*$ngc; $pxc=sprintf"%.3g", 100*$difngc /$ngt; @p=map{ sprintf"%.3g", 100*$_/$ngt } @x; $xc=sprintf"%.2f",$xc; print "$clab: ngc=$ngc,ngt=$ngt,xc=$xc, pxc=$pxc, x=@x, xp=@p; \n";'
C2.at18: ngc=806,ngt=27347,xc=0.74, pxc=3.97, x=319 458 29, xp=1.17 1.67 0.106; 

==> dmag20gnodes/dmag19skasm_test8f7sumgcn_genesum.txt <==
 CLevel	nGene	gCopy	xCopy	aCopy	aMiss	x0,xLo,xEq,xHi
 All	25322	3.0	1.0	1.1	6.2%,983605	224,3589,19448,2061
 1	19352	1.0	1.0	1.0	0.55%,37472	52,0,18808,492
 2-9	3492	3.4	0.5	1.3	2.5%,80098	36,2809,640,7
 10-99	718	28.9	0.1	2.0	1.6%,59802	11,707,0,0

==> dmag20gnodes/dmag20sk4maca20ok_test8f7sumgcn_genesum.txt <==
 CLevel	nGene	gCopy	xCopy	aCopy	aMiss	x0,xLo,xEq,xHi
 All	25322	3.0	1.2	1.4	0.02%,3356	20,2744,18518,4040
 1	19352	1.0	1.2	1.2	0.028%,1902	0,0,17078,2274
 2-9	3492	3.4	0.7	2.1	0.0089%,292	0,1967,1410,115
 10-99	718	28.9	0.2	4.7	0.012%,455	0,688,30,0

perl -e '$clab="C2.dam19"; @x=(2809,640,7); $xc=1.3/3.4; $ngc=3492; $ngt=25322; $difngc=(1/$xc)*$ngc; $pxc=sprintf"%.3g", 100*$difngc /$ngt; @p=map{ sprintf"%.3g", 100*$_/$ngt } @x; $xc=sprintf"%.2f",$xc; print "$clab: ngc=$ngc,ngt=$ngt,xc=$xc, pxc=$pxc, x=@x, xp=@p; \n";'
C2.dam19: ngc=3492,ngt=25322,xc=0.38, pxc=36.1, x=2809 640 7, xp=11.1 2.53 0.0276; 

perl -e '$clab="C2.dam20"; @x=(1967,1410,115); $xc=2.1/3.4; $ngc=3492; $ngt=25322; $difngc=(1/$xc)*$ngc; $pxc=sprintf"%.3g", 100*$difngc /$ngt; @p=map{ sprintf"%.3g", 100*$_/$ngt } @x; $xc=sprintf"%.2f",$xc; print "$clab: ngc=$ngc,ngt=$ngt,xc=$xc, pxc=$pxc, x=@x, xp=@p; \n";'
C2.dam20: ngc=3492,ngt=25322,xc=0.62, pxc=22.3, x=1967 1410 115, xp=7.77 5.57 0.454; 


#-----
 
  # Moments for all genes in sumgcn_genesum.txt .. want these?
  Field   Median  Mean    N       StDev   Skew    Range   Correlations
  gLen    1182    1606    14329   2335    0.545   96,56544
  gNread  829     1219    14329   3536    0.331   52,218862
  gCopy   1       1.57    14329   7.74    0.220   0.400,320       r.gLen:0.121,r.gNread:0.838
  aCopy   1       1.31    14329   2.19    0.428   1,33.1  r.gCopy:0.354,r.gLen:0.378,r.gNread:0.270
  aCopyRd 1       1.31    14329   2.17    0.428   1,32.1  r.gCopy:0.354,r.gLen:0.381,r.gNread:0.271
  gCM     88      107     14329   491     0.116   30,27672
  aCovT   75      93.3    14329   228     0.241   3,10316 r.gCM:0.809,r.gLen:0.274,r.gNread:0.809
  aReadT  894     4389    14329   47274   0.222   52,2158093      r.gCM:0.220,r.gLen:0.0979,r.gNread:0.310
  aReadZ  0       0.0091  14329   0.150   0.182   0,8

=item try1
 head -25  aweed*gnodes/{*test8f_SRR*_sum.txt,*test8f7sumgcn_genesum.txt} |  $gw/daphnia/gnodes_genecopyf/gnodes_sumasmgcn2rtab.pl
 
SppAsm  Item    pXlo    pXeq    pXhi    pMiss
ara18ch Aall    67.11   119     177
ara18ch Acds    68.49   50      73
ara18ch Adup    25.77   11.6    44.8
ara18ch Auniq   81.30   108     133
ara18ch C01     0.00    99.61   0.39    0.00
ara18ch C02-9   39.58   56.82   3.60    0.00
ara18ch C10-99  23.53   60.29   16.18   0.01
ara18ch CAll    1.27    97.26   1.47    0.00
ara18ch G01     26200   1.0     1.0
ara18ch G02-9   806     3.1     2.3
ara18ch G10-99  68      21.0    17.5
ara18ch GAll    27347   1.3     1.1

ara20ma Aall    69.44   119     171
ara20ma Acds    72.46   52      72
ara20ma Adup    33.67   21      62
ara20ma Auniq   90.09   98      110
ara20ma C01     0.00    98.78   0.78    0.90
ara20ma C02-9   53.23   42.56   2.98    2.90
ara20ma C10-99  79.41   14.71   5.88    4.60
ara20ma CAll    1.81    95.93   1.78    1.70
ara20ma G01     26200   1.0     1.0
ara20ma G02-9   806     3.1     1.9
ara20ma G10-99  68      21.0    13.6
ara20ma GAll    27347   1.3     1.2

ara21bw	Aall	68.03	119	176
ara21bw	Acds	68.97	51	74
ara21bw	Adup	25.77	11.4	44.4
ara21bw	Auniq	81.97	108	132
ara21bw	C01	0.00	99.65	0.35	0.00
ara21bw	C02-9	50.96	45.31	3.73	0.00
ara21bw	C10-99	35.00	56.67	8.33	0.01
ara21bw	CAll	1.97	96.66	1.37	0.01
ara21bw	G01	26040	1.0	1.0
ara21bw	G02-9	991	3.0	2.1
ara21bw	G10-99	60	20.1	12.0
ara21bw	GAll	27347	1.3	1.1

ara21he	Aall	66.67	119	179
ara21he	Acds	63.29	49	78
ara21he	Adup	31.35	15.2	48.5
ara21he	Auniq	80.00	104	130
ara21he	C01	0.00	99.11	0.89	0.00
ara21he	C02-9	55.60	36.53	7.87	0.01
ara21he	C10-99	31.11	57.78	11.11	0.00
ara21he	CAll	2.05	93.33	4.62	0.00
ara21he	G01	24488	1.0	1.0
ara21he	G02-9	928	2.8	2.9
ara21he	G10-99	45	21.9	16.6
ara21he	GAll	26398	1.3	1.1


=cut
