#!/usr/bin/env perl
# gnodes_std2genecov.pl

=item about gnodes_std2genecov

  this is gnodes_sam2genecov.pl  for samtools depth data

  samtools depth cds.bam > geneid gbase covdepth table (aka pileup table)
  using st depth -g 0x100 for multimap aCovT, st depth -G 0x100 for onemap aCovM,
  so input table is (id base covT covM)

  accessory input from samtools flagstats, for readlen, nreads, nmapreads
    
=item UPD220222

=cut

use strict;
use Getopt::Long;  

use constant UPD220222 => 1; # == UPD22FEB
use constant UPD21JUN => 1; #  
use constant UPD21DEC => 1; #  test gene family read-share matrix

my $debug= $ENV{debug}||1;
my $pMAXDUP= 0.05; #?? need tests
my $RDLEN= $ENV{rdlen}||150; #?? guess, need here
my $MIN_IDENT = 0.40; # UPD7e was 0.65; # lo is best? << keep same as gnodes_sam2covtab.pl, or lower for gene x DNA align deficit
my $MIN_DUPIDENT = 0.98; #UPD7e was 0.98; # was .99/1; hi is best? or 1.0; # == ident equal to 1st/top align, lower if desired
   ## gnodes_sam2covtab7e now has  MIN_DUPIDENT=0.999; lower that?, raise this? or what?
my $TrimZEROEnds = $ENV{trimzero}||0; # UPD220227 TEST option to reduce partial read cover at ends, introns of CDS align
    ## tests say TrimZEROEnds is not useful as is, minor improvments, test another way?
    
my $UCGenePatt='UCG|busco';
my $ContamGenePatt='contam'; #UPD220222, for plants, other?
my($optFullAlignFilter,$opt1000filter,$optUerrFilter)=(1, 1, 0); # default fullalignfilter=on, others?

my ($crclassf,$rdstats,$outtab,$intab,$dosumtab)=("") x 9; 
my @intab; # opt to allow many -bam inputs, -bam one.bam two.bam three.bam ... or .sam ;  

my $optok= GetOptions( 
  'input=s', \@intab, # bam or sam or STDIN ? for samtools view call, add opt many @inbam
  'stats|readstats|flagstats=s',\$rdstats, 
  'output|covtab=s',\$outtab, 
  'crclassf|idclassf=s',\$crclassf, # alt table chr => class for savereadids
  # 'maxerr=s',\$MINE, #?? reuse for NM:i:(err) filter # replace w/ minident ?
  # 'minident=s',\$MIN_IDENT, #? only used to set MINALN ?
  'maxdup=s', \$pMAXDUP,  #? reuse MIN_DUPIDENT for dup reads/total reads filter |mindupident
  # 'mindupident=s', \$MIN_DUPIDENT,  # as per sam2covtab; pMAXDUP = 1 - $MIN_DUPIDENT ??
  # 'minrdlen|minreadlen=i',\$MINRDLEN, # add .. not same as MINALN ??
  'rdlen|readlen=i',\$RDLEN, # NEED;
  'UCGeneAnnot=s', \$UCGenePatt, # 
  'filterFullAlign!', \$optFullAlignFilter, # -nofilterFull to turn off
  'filter1000bases!', \$opt1000filter, # -nofilterFull to turn off
  # 'filterUerr!', \$optUerrFilter, 
  # 'heterozygosity!', \$DOHETZ, #? drop this
  # 'HOMIN=s', \$HOMIN, # for DOHETZ
  #'summary|sumtab!', \$dosumtab, 
  'debug!', \$debug, 
  );

if(@ARGV) { push @intab, grep(/\.(tab)$/, @ARGV); }

unless($optok) {
  die "usage: gnodes_std2genecov.pl -in genes_readmap.stdepth.tab -output genes_readmap.covtab
    or: gnodes_std2genecov.pl < genes_readmap.stdepth.tab > genes_readmap.covtab
    opts: -idclass genes.idclass : gene id table w/ UCG|busco class flag
    -maxdup=$pMAXDUP portion for uniq class, 
  ";
}

my($ngene,  $nmapid, $stw, $scm, $srdlen, $srdlen_nomap, $lrid)= (0) x 9; # $nreadid, $nr
my($n_readid, $nrdlen_tot, $n_nomap, $n_mapok, $n_mapbad, $n_dupbad, $n_mapnotucg, $n_rdtooshort, # other: $n_dupbad, $n_partb,  
   $n_intron, $n_insert, $n_delete, $n_softclip, $n_mismatch)= (0) x 19; # globals?
my(%tw, %nrd, %nrdlen, %loc, %locdup, %cmd);
my(%rdshare,%genetab); #UPD21DEC

sub MAIN_stub {}

# if($dosumtab){
#   ucgcov_sumtab(); 
#   exit;
# }  

my $outh= *STDOUT; if($outtab) { rename($outtab,"$outtab.old") if(-s $outtab); open($outh, '>', $outtab) or die "writing $outtab"; }

my $topline= "#gnodes_sam2genecov options: minident=$MIN_IDENT, mindupid=$MIN_DUPIDENT, maxdup=$pMAXDUP,  trimends=$TrimZEROEnds, idclass=$crclassf\n"; # maxerr=$MINE,
print $outh $topline; warn $topline if($debug and $outtab);

# * idclassh only used for isUCGene()
my($nidclass,$idclassh,$idclasslist)= ($crclassf) ? read_idclass($crclassf,1) : (0); # cds,te class by id, $idclass->{id} = class

my($nreads, $nmap)= readmapstats($rdstats); # == $nrdlen_tot,$n_mapok globals

use constant UCGENE_RDFILT => 0; #UPD21Jun04:off, 1 == skip read counts for not-UCGene(id)
my $testUCG= ($nidclass>0 and $UCGenePatt)?1:0;

# $intab= shift @intab if(@intab and not $intab);
# are intabs separate read sets, << assume this, for now
# or same reads, diff data slices? covm, covt, covcds, covte ?? need opts

my($nid,$lid,@covt,@covm)=(0,0);
my($inok,$inh,$intab);
if(@intab) { ($inok,$inh,$intab)= openRead(shift @intab); }
else { ($inok,$inh,$intab)= (1,*STDIN,'STDIN'); }

while($inok) {
  # warn "# read($intab)=$inok\n" if($debug); 
  
  while(<$inh>){
    next if(/^\W/);
    my($id,$ib,@cov)=split; my $ncov=@cov;
    
    if($id ne $lid) {
      $nid++ if(putgened($lid,\@covm,\@covt));
      @covm=@covt=();
    }
    
    my ($ct,$cm,$cu,$ccds)=@cov;
    if($ncov==1) { $cm=$ct; }
    push @covm, $cm; push @covt,$ct;
    $lid=$id;
  }
  $nid++ if(putgened($lid,\@covm,\@covt));
  close($inh); 
  warn "# read($intab)=$nid\n" if($debug); 
  ($inok,$inh,$intab)= openRead(shift @intab);
}


my ($topinfo,$ucmwMedn,$ucmwAve,$ucmwSE)=("",0,0,0);
use constant { kUCGfilter => 1, kIsUniqfilter => 2, k1000filter => 4, k2000filter => 8,
              kFullAlignFilter => 16, kUerrFilter => 32, };

if(UPD21JUN) {
  #>> only if have idclass w/ busco/UCG annots == my $testUCG= ($nidclass>0 and $UCGenePatt)?1:0;

  my $cucg_filt= kIsUniqfilter ; # ?opt for k2000filter
  $cucg_filt += k1000filter if($opt1000filter); #def on
  $cucg_filt += kFullAlignFilter if($optFullAlignFilter); #def on
  $cucg_filt += kUerrFilter if($optUerrFilter);
  my($btopinfo,$bucmwMedn,$bucmwAve,$bucmwSE)= topsummary_UCG("Measured Unique Genes",$cucg_filt); #  >= 1k long

  if($testUCG) {
    $cucg_filt=  kUCGfilter;
    $cucg_filt += k1000filter if($opt1000filter); #def on
    $cucg_filt += kFullAlignFilter if($optFullAlignFilter);
    $cucg_filt += kUerrFilter if($optUerrFilter);
    my($atopinfo,$aucmwMedn,$aucmwAve,$aucmwSE)= topsummary_UCG("UCG Class Genes", $cucg_filt); #  >= 1k long: now label all Filters:

    my $adiff= abs( $aucmwMedn - $aucmwAve);
    my $bdiff= abs( $bucmwMedn - $bucmwAve);
    if($adiff <= $bdiff) {
      $topinfo= $atopinfo . $btopinfo;
      ($ucmwMedn,$ucmwAve,$ucmwSE)= ($aucmwMedn,$aucmwAve,$aucmwSE); 
    } else {
      $topinfo= $btopinfo . $atopinfo;
      ($ucmwMedn,$ucmwAve,$ucmwSE)= ($bucmwMedn,$bucmwAve,$bucmwSE); #<< FIXME need choice from stats
    }
  } else {
    ($topinfo,$ucmwMedn,$ucmwAve,$ucmwSE)= ($btopinfo,$bucmwMedn,$bucmwAve,$bucmwSE);
  }
  print $outh $topinfo; warn $topinfo if($debug and $outtab);
} 


genestable($outh, $ucmwMedn,$ucmwAve,$ucmwSE);

# if(UPD21DEC) { 
#   my $gfouth= *STDOUT; 
#   if($outtab) {
#     my $gfout= $outtab; $gfout =~ s/\.\w+$//; $gfout.=".genefamtab"; # replaces .genexcopy ?
#     rename($gfout,"$gfout.old") if(-s $gfout); 
#     open($outh, '>', $gfout) or die "writing $gfout"; 
#   }
#   print $gfouth "# genefam_table \n";  
#   genefam_table($gfouth, $ucmwMedn,$ucmwAve,$ucmwSE);
# }

my $botinfo= join", ", grep /\w/, map{ my $bn=`basename $_ .bam`; chomp($bn); ($bn)?$bn:""; } (@intab);
print $outh "# intab: ",$botinfo,"\n" if($botinfo);  

#======================

use constant C_MEDIAN => 1; # use median not ave, tested both, ave influenced by extremes
use constant USE_CNZMEDIAN => 1; # test, then switch, looks better than cmw .. need for UCG C_MEDIAN also
use constant { 
  sZERO => 9, # for sums
  cZERO => 2, # for C vals
  ZERO => 0, # for cov depth
  NoZEROS => 1, # skip depth == 0 as gaps
  spanZERO => 50, # percent of gene span covered, maybe higher: 66%? 75% most appear >90%, poorqual < 75%
  XCOPYisDUP => 1.66,
  cSKEW => 1.85, # C.nz.ave >> C.nz.median
  };

sub _max{ return($_[0] < $_[1])?$_[1]:$_[0]; }
sub _min{ return($_[0] > $_[1])?$_[1]:$_[0]; }
sub _minNot0{ return($_[0] == 0 or $_[0] > $_[1]) ? $_[1] : $_[0]; }

=item isUCGene

  isUCGene if idclass flag 'UCG|busco' says it is unique-conserved-gene
  
  UPD2202: add contam class genes, ie from MT/Chloroplast, likely highcopy num spurious for some uses
  11 contam,MT  AT2G07827t1
   6 contam,MT  AT2G07835t1
  19 contam,Pltd        AT2G07565t1
   4 contam,Pltd        AT2G07732t1
  
=cut

sub isUCGene { # expand w/ idclass patt
  my($id)= @_;
  if($testUCG) { # allow also for all-are-ucg option
    my $idclass= $idclassh->{$id} || "";
    return($idclass =~ m/$UCGenePatt/i)?1 : ($idclass =~ m/$ContamGenePatt/i)? -1 : 0; # UPD2202
    # return($idclass =~ m/$UCGenePatt/i)?1:0; # 
  } else {
    return 1; # no choice
  }  
}

#      putgened($lid,\@covm,\@covt) if($lid);
sub putgened {
  my($td,$covm,$covt,$covu,$ccds,$cte,$cunk)=@_;
  return unless($td and ref($covm) =~ /ARRAY/);

  use constant USE_TCOV => 1; # FIXME: this way?

  my($tw,$noz,$md,$av,$sd,$skw,$pcov,$sum)= nmastat(@$covm);
  $tw{$td}= $tw;
  $locdup{$td}{'N'}= $cmd{$td}= $sum;
  $locdup{$td}{'M'}= $sum; # new val
  
  ## reuse locdup
  # my $sln= $locdup{$td}{'N'}||0;  
  # my $sld= $locdup{$td}{'D'}||0;
 
  my($neq,$sumu)=(0,0);
  my($tnv,$tnoz,$tmd,$tav,$tsd,$tskw,$tpcov,$tsum)=
    ($tw,$noz,$md,$av,$sd,$skw,$pcov,$sum); #??
  if($covt and @$covt) {
    ($tnv,$tnoz,$tmd,$tav,$tsd,$tskw,$tpcov,$tsum)= nmastat(@$covt);
    $locdup{$td}{'D'}= $tsum - $sum;
    if(USE_TCOV) {
    $locdup{$td}{'N'}= $cmd{$td}= $tsum; #?? FIXME is this proper N sum for genexcopy tab?
    }
      
    for(my $i=0; $i<$tnv; $i++) {
      my $ct=$$covt[$i];      
      if($ct>0 and $ct == $$covm[$i]) { $neq++; $sumu += $ct; }
      # else? sumu += what ?? fraction of 
    }
    $locdup{$td}{'U'}= $sumu;
    $locdup{$td}{'Un'}= $neq; 
  } else {
    ($neq,$sumu)=($noz,$sum);
  }
  
  my $copyn= ($sum<1)?0:$tsum/$sum;  
  my $uok= ($sum<sZERO or 100*$pcov < spanZERO)?"zero"
      :(abs($skw) > cSKEW)?"skew"
      # :($copyNumEst >= XCOPYisDUP)?"dupx"
      :($copyn >= XCOPYisDUP)?"dupx"
      :( 1 - $sumu/$sum > $pMAXDUP)?"dups" #? or 1 - $neq/$noz; 
      :"uniq"; 
  
  if(USE_TCOV) { #?? include both tmd.., md vals here?
  $genetab{$td}=[$tw,0,$copyn, $tmd, $uok, $tav,$tmd,$tnoz,100*$tpcov,$tskw]; #? reus this
  } else {
  $genetab{$td}=[$tw,0,$copyn, $md, $uok, $av,$md,$noz,100*$pcov,$skw]; #? reus this
  }
  #my $cnz= sprintf"%.2f,%d,%d,%.0f,%.1f",$av,$md,$noz,100*$pcov,$skw;
  # $genetab{$td}=[$tw,0,$copyn, $md, $uok, $cnz]; #? reus this
  # $genetab{$td}=[$tw,$rc,$copyNumEst, $CMint, $uok,"$cnz,$cnzmed,$nnz,$cspan"]; #? reus this
  return (wantarray) ? ($tw,$md,$uok) : $tw;
}

sub nmastat { 
  my @v=@_; 
  my($nv,$noz,$md,$av,$sd,$skw,$pcov,$s,$ss,$nh,$i,$v,$ns) = (0) x 19; 
  $nv=@v; $noz= scalar( grep{ $_ > ZERO } @v );  
  $nh= (NoZEROS)?$noz:$nv; 
  return($nv, (0) x 7) if($nh == 0);
  my @sv=sort{$b <=> $a} @v;  
  my $md=$sv[ int($nh/2)]; 
  if($TrimZEROEnds) {
    my $ezero= $TrimZEROEnds; # min cov here > 0, like 1,2,3 ?
    for($i=0; $i<$nv; $i++){ 
      my $j=($i>0)?$i-1:$i; my $k=($i<$nv-1)?$i+1:$i; 
      my($lv,$v,$fv)=@v[$j,$i,$k]; 
      if($v>$ezero and $lv>$ezero and $fv>$ezero){ $s+=$v; $ss+=$v*$v; $ns++; } 
    }
    $noz= $ns; #? also this, since sum s has only ns 
  } else {
  $ns=$nh;
  for($i=0; $i<$nh; $i++){ $v=$sv[$i]; $s+=$v; $ss+=$v*$v; } 
  }
  if($ns>0){ $av=$s/$ns; $sd=($ns<2)?0:($ss - $av*$av)/($ns-1); $sd=sqrt($sd); } 
  $skw=$av - $md; $pcov=($noz/$nv); 
  return ($nv,$noz,$md,$av,$sd,$skw,$pcov,$s); 
}
  
sub genestable {
  my($outh,$ucmwMedn,$ucmwAve,$ucmwSE)= @_;
  my $didhdr=0;
  my $Cmap= (C_MEDIAN) ? $ucmwMedn : $ucmwAve; # was: $mcmw : $ucmt;
  
  # main gene data globals: %tw, %nrd %cmd  %loc{$td}[0..tlen]{DUN}
  for my $td (sort keys %tw) { 
    my $tw=$tw{$td}; my $cmd=$cmd{$td}||0; 
    my $rc=$nrd{$td}||0; my $rclen=$nrdlen{$td}||0;
    my $cmw= $cmd/$tw; my $crl= $rclen/$tw; # drop int(); was int($rc*$RDLEN/$tw); 
  
    my @lns=(); # median of cov N val for gene span bases, maybe should replace $cmd{td}/tw, cmp to cnz
    my($sld,$slu,$sln,$slerr,$slun,$neq,$nnz, $shet, $shetho)=(0) x 19;

    my($ttw,$trc,$tcn,$tcm,$tcla,$tcav,$tcmm,$tnz,$tspan,$tskw)= @{$genetab{$td}};
    #?? also add covM vals: $mcn,$mcm?,$mnz,$mspan,$mskw
    
#  my $cnz= sprintf"%.2f,%d,%d,%.0f,%.1f",$av,$md,$noz,100*$pcov,$skw;
#  $genetab{$td}=[$tw,0,0, $md, 0, $cnz]; #? reus this
#  $genetab{$td}=[$tw,$rc,$copyNumEst, $CMint, $uok,"$cnz,$cnzmed,$nnz,$cspan"]; #? reus this

    $slerr= $locdup{$td}{'Uerr'}||0;
    $sln= $locdup{$td}{'N'}||0; 
    $sld= $locdup{$td}{'D'}||0; 
    $slu= $slun= $locdup{$td}{'U'}||0; 
    $neq= $locdup{$td}{'Un'}||0; 
    $nnz=$tnz;
    
    #>> change this per gene loc{td}[i] to locdup{td}{(D U N)} ?? need neq, nnz ?
    # for(my $i=0; $i<$tw; $i++) { 
    #   my($ld,$lu,$ln)= map{ $loc{$td}[$i]{$_}||0 } qw(D U N); 
    #   $sld+=$ld; $slu+=$lu; 
    #   if($ln>0){ 
    #     $sln+=$ln; $nnz++; push @lns,$ln; 
    #     if($lu == $ln){ $slun+=$lu; $neq++; }
    #   }
    # }

    my $ceq=($neq<sZERO)?0:sprintf"%.1f", $slun/$neq;
    my $cnz=($nnz<sZERO)?0:sprintf"%.1f", $sln/$nnz;
    my $cnzmed=$tcm; 
    # if($nnz>=sZERO){ @lns=sort{$b<=>$a}@lns; $cnzmed= sprintf"%.0f", $lns[int($nnz/2)]; }
    my $cspan= sprintf"%.1f",$tspan; # ($tw<sZERO)?0:sprintf"%.1f", 100*$nnz/$tw;
    
    my $merr=($slu<sZERO)?0:sprintf"%.1f",  100*$slerr/$slu;
    # my $phet=($shetho<sZERO)?0:sprintf"%.0f,%d/%d",  100*$shet/$shetho,$shet,$shetho;

    my($CMint,$copyNumEst)=(0,0);
    if(USE_CNZMEDIAN) {
      $CMint= $cnzmed; $copyNumEst= sprintf"%.1f", $cnzmed/$Cmap; #est from cmw or cnz = $sln/$nnz ?
    } else {
      $CMint=  sprintf"%.0f", $cmw; # not int($cmw);
      $copyNumEst= sprintf"%.1f", $cmw/$Cmap; #est from cmw or cnz = $sln/$nnz ?
    }

    my $cskewd= ($cnz > cSKEW * $cnzmed)?1:0; # what cSKEW cut off?
    #?? $cskewd= (abs($tskw) > cSKEW)?1:0;
    
    # ^^ modify uok class: sld/sln>maxdup == dupls, xcopy >= 1.66 == duplx
    my $uok= ($sln<sZERO or $cspan < spanZERO)?"zero"
      :($cskewd)?"skew" # this is too many? displacing dupx
      :($copyNumEst >= XCOPYisDUP)?"dupx"
      :($sld/$sln > $pMAXDUP)?"dups" #? keep or not
      :"uniq"; 

    my $ucgf="";  # add isUCGene() flag ?? which col? S.D,U,N,E?
    if($testUCG) { my $gcl= isUCGene($td); $ucgf= ($gcl==1)?",UCG":($gcl==-1)?",contam":""; }
    # if($testUCG) { $ucgf= (isUCGene($td))?",UCG":""; }
    
    unless($didhdr++) {
      my @hcols= qw(Gene_ID Glen Nread Rdlen tCopy C.M Uniq Merr C.nz S.Dup,Unq,Nt,Equ,UCG);
      print $outh join("\t",@hcols)."\n";
    } 
    my @cols=($td,$tw,$rc,$rclen, $copyNumEst, $CMint, $uok, $merr, "$cnz,$cnzmed,$nnz,$cspan", "$sld,$slu,$sln,$slun$ucgf");  
    # push @cols,$phet if($DOHETZ);
    print $outh join("\t",@cols)."\n"; # off: ,"$ceq,$neq"

    if(UPD21DEC) { $genetab{$td}=[$tw,$rc,$copyNumEst, $CMint, $uok,"$cnz,$cnzmed,$nnz,$cspan"]; } # add $ucgf == UCG|contam|0
  }
}  


sub topsummary_UCG {
  my($sumlabel, $whichUniqGeneFilter)=@_;
  $whichUniqGeneFilter= kUCGfilter unless($whichUniqGeneFilter);
  $sumlabel ||= "UCG Class Genes";
  my $skipALLsum= (UPD21JUN)?1:0; 
  
  my $ngene= scalar(keys %tw); # tw now includes not-isUCGene
  my $stw=0; for my $td (sort keys %tw){ $stw += $tw{$td}; }

  my $cspan_min = 90; #? option? dont want too low,
  
  use constant pUERR_MAX =>  0.33; #? test opts? 0.25? 0.50?
  my $uerr_max= 999999;
  if($whichUniqGeneFilter  & kUerrFilter) { # was unless($ALN_SUB_NMI)  # drop this, no good (?)
    my @uerr= sort{ $b <=> $a } map{ $locdup{$_}{'Uerr'}||0 } keys %locdup; #? limit ids to UCG here?
    my $nerr= @uerr;
    $uerr_max= $uerr[ int( pUERR_MAX * $nerr) ]; # drop worst 1/3, 1/4 ?
  }
   
  # totals including dup genes
  my $rct=($stw<sZERO)?0:sprintf "%.1f", $srdlen / $stw; #was $n_mapok * $RDLEN / $stw ;
  my $cmt=($stw<sZERO)?0:sprintf "%.1f", $scm / $stw;  
  my $ardlen= ($n_mapok<sZERO)?0:sprintf "%.0f", $srdlen / $n_mapok; 
  
  my ($n_contam, $slen_contam,$srdlen_contam)= (0,0,0); #UPD22mar
  my $srdlen_tot= $srdlen+$srdlen_nomap; 
  my $nrdlen_tot= $n_mapok + $n_nomap + $n_mapnotucg + $n_mapbad; #  nrdlen_tot should == n_readid **

  my $labelFILT="Filters:";
  $labelFILT .= "UCG," if($whichUniqGeneFilter & kUCGfilter);
  $labelFILT .= "size>=2k," if($whichUniqGeneFilter & k2000filter);
  $labelFILT .= "size>=1k," if($whichUniqGeneFilter & k1000filter);
  $labelFILT .= "covspan>=$cspan_min," if($whichUniqGeneFilter & kFullAlignFilter);
  $labelFILT .= "mErr<$uerr_max," if($whichUniqGeneFilter & kUerrFilter);
  
  my($urdlen,$utw,$ucm,$ucmw,$ssucmw)=(0) x 9; 
  my(@crl,@cmw);
  for my $td (sort keys %tw) { 
    my $tw= $tw{$td} or next; 
    my $gok=1;
    
    my $gclass= isUCGene($td); # UPD2202: -1 = contam, 1 = UCG, 0 = other
    # (UPD21JUN) 
    if($whichUniqGeneFilter & kUCGfilter) { $gok=0 unless( $gclass > 0); }
    elsif($whichUniqGeneFilter & kIsUniqfilter) { } # == ($sln < sZERO or $sld/$sln > $pMAXDUP) below
    if($whichUniqGeneFilter & k2000filter) { $gok= 0 if($tw < 2000); }
    elsif($whichUniqGeneFilter & k1000filter) { $gok= 0 if($tw < 1000); }
    else { $gok=0 if($tw < 200); } #? < 200 lower limit, eg arabid has a few 32 bp cds, other ref prots < 100aa
 
    my($ttw,$trc,$tcn,$tcm,$tcla,$tcav,$tcmm,$tnz,$tspan,$tskw)= @{$genetab{$td}};
    # $slerr= $locdup{$td}{'Uerr'}||0;
    # $sln= $locdup{$td}{'N'}||0; 
    # $sld= $locdup{$td}{'D'}||0; 
    # $slun= $locdup{$td}{'U'}||0; 
    # $neq= $locdup{$td}{'Un'}||0; 
    # $nnz=$tnz;
   
    my $cspan= 100*$tspan; # -1; # always calc, add to sums? No use here, is in genestable()
    if($gok and $whichUniqGeneFilter &  kFullAlignFilter) {
       $gok=0 if($cspan < $cspan_min);
    }
    # if($gok and $whichUniqGeneFilter &  kFullAlignFilter) {
    #   my $nnz=0;
    #   for(my $i=0; $i<$tw; $i++) { 
    #     my $ln= $loc{$td}[$i]{'N'}||0; 
    #     $nnz++ if($ln>0);
    #   }
    #   $cspan= 100*$nnz/$tw;
    #   $gok=0 if($cspan < $cspan_min); # kFullAlignFilter
    # }
    
    ##UPD21NOV: ALN_SUB_NMI here: nmi only used for $locdup{$td}{'Uerr'} += $nmi
    if($gok and $whichUniqGeneFilter &  kUerrFilter) {
      my $uerr= $locdup{$td}{'Uerr'}||0; # == slerr below
      $gok=0 if($uerr > $uerr_max); # kUerrFilter
    }
    
    # UPD2202: measure contam for -Contam
    if($gclass == -1) {
      $n_contam++; $slen_contam += $tw;
      $srdlen_contam += $locdup{$td}{'M'}||0; # == sum bases align to contam  cmd and N includes dup cov, M is covM, 
    } 
    
    next unless($gok);
    
    my $sln= $locdup{$td}{'N'}||0;  
    my $sld= $locdup{$td}{'D'}||0;
    # check td span cover for < spanZERO as bad data, notuniq
    next if($sln < sZERO or $sld/$sln > $pMAXDUP); # $isuniq=0 
    
    my $cmd= $cmd{$td}||0; 
    my $rclen=$nrdlen{$td}||0; # any zero/missing?  
    my $crl= $rclen/$tw;  # dont int() these, lost precision..
    my $cmw= $cmd/$tw; # primary C.m val == read aligns sum/gene from addCigar; 
                       # /tw is problem for partialAlign genes, thus kFullAlignFilter;
                       # alternate fix: use $nnz, ie covered span, but that affects meaning of C.m, can use C.nz tuple vals to correct
    $urdlen += $rclen; $ucm += $cmd; $utw += $tw; 
    $ucmw += $cmw; $ssucmw += $cmw*$cmw; # use ucmw for ave, not ucm/utw
    push @crl,$crl; push @cmw,$cmw;  
  }
  
  @crl= sort{ $b <=> $a } @crl; @cmw= sort{ $b <=> $a } @cmw;
  my $nc= @cmw; my $nch= int($nc/2);
  my $mclr= $crl[$nch]; my $mcmw= $cmw[$nch];
  
  my($ucmwMedn,$ucmwAve,$ucmwSD,$ucmwSE)=($cmw[$nch],0,0,0);
  if($nc > sZERO) {
    $ucmwAve= $ucmw/$nc; 
    $ucmwSD= sqrt(($ssucmw - $ucmwAve*$ucmwAve)/($nc-1)); 
    $ucmwSE= $ucmwSD/sqrt($nc); #? change to report ucmwSD , nc ?
  }
  
  # problem of extremes in  ucm/utw, use instead meanof(@cmw) as better stat
  my $urct=($utw<sZERO)?0:sprintf "%.1f", $urdlen / $utw; #was $n_mapok * $RDLEN / $stw ;
  my $ucmt= $ucmwAve; # old: ($utw<sZERO)?0:sprintf "%.1f", $ucm / $utw;  
  
  my $nmapcdsid= $nmapid+$n_mapnotucg;
  my $pmapucg= ($n_readid<sZERO)?0:sprintf"%.1f", 100*$nmapid/$n_readid; # pmr FIXME should reduce nmapid by filters, ?? here?
  my $pmapcds= ($n_readid<sZERO)?0:sprintf"%.1f", 100*$nmapcdsid/$n_readid;
  
  my $nberr = $n_mismatch + $n_insert + $n_delete; # mid or sid; see also new locdup{td}{Uerr}
  my $pberr = ($scm<sZERO)?0:sprintf "%.2f", 100*$nberr / $scm;  
  
  my $Cmap= (C_MEDIAN) ? $ucmwMedn : $ucmwAve; # was: $mcmw : $ucmt;
  my $Genosize_mb= ($Cmap<cZERO) ? 0 :  sprintf"%.1f",  $ardlen * $n_readid / $Cmap / 1_000_000; 
  my $CDS_mb= ($Cmap<cZERO) ? 0 :  sprintf"%.1f",  $ardlen * $nmapcdsid / $Cmap / 1_000_000; 
    
  my $GSEeqn= sprintf "L*N/C= %d * %d / %.1f",$ardlen,$n_readid,$Cmap;
  
  my $addGSEalt="";
  my $srdlen_totEST= ($n_mapok>=sZERO) ? $srdlen * $n_readid / $n_mapok : 0; 

  my $pDIFFLEN= 0.999;
  my $srdlen_totc= $srdlen_tot - $srdlen_contam;
  my $dogsalt= ( $n_mapok>=sZERO and ratio($srdlen_totc, $srdlen_totEST) >= $pDIFFLEN ) ? 0 : 1;
  if($dogsalt){
    my $Genosize_alt= ($Cmap<cZERO) ? 0 :  sprintf"%.1f", $srdlen_totc / $Cmap / 1_000_000; 
    my $ardlen_alt  = ($nrdlen_tot<1)? 0 : sprintf"%.1f", $srdlen_totc / $nrdlen_tot;
    my $GSEeqn_alt  =  sprintf "LN/C= %d / %.1f for Lr: $ardlen_alt, Nr: %d = %d+ %d+ %d+ %d [ok,nomap,noucg,bad]", 
      $srdlen_totc, $Cmap,
      $nrdlen_tot, $n_mapok, $n_nomap, $n_mapnotucg, $n_mapbad; # srdlen_tot sum of all readlens counted, *should* == ave(rdlen) * n_readid
    if($srdlen_contam>0){
      my $addalt= sprintf "LN: %d, L: %d, N: %d contam", $srdlen_contam, $slen_contam, $n_contam;
      $GSEeqn_alt .= "\n# Gsize_alt removed Contam= $addalt";
    }
    $addGSEalt=  "# Gsize_alt  = $Genosize_alt Mb of $GSEeqn_alt\n";
  }
  #UPD2202: change precision .1f to .2f for the crucial Cucg vals
  my $CmapAveSE= sprintf "%.2f, %.2f +/-%.2f (mdn,ave,sem)",$ucmwMedn,$ucmwAve,$ucmwSE; # replace C.Map/W=$mcmw, $ucmt
  $mclr=int($mclr);
  
  my $alldepth=($skipALLsum)?"":"# All  Gene Cov Depth n=$ngene, C.Map/W=$cmt ave, for W.genes=$stw, LN.reads=$srdlen\n"; # C.LN/W=$rct, 
  my $showucg=($n_mapnotucg > 0) ? " Nr.ucgcds=$nmapid ($pmapucg%),":"";
  
  my $topinfo= "# $sumlabel, $nc of $ngene total, $labelFILT --------------\n" . $alldepth . 
  "# Uniq Gene Cov Depth n=$nc, C.Map/W=$CmapAveSE for W.genes=$utw, LN.reads=$urdlen\n" .   # C.LN/W=$mclr,  
  "# Genome_size= $Genosize_mb Mb of $GSEeqn, CDS_size= $CDS_mb Mb,\n" .
  $addGSEalt .
  "#   for $showucg Nr.allcds=$nmapcdsid ($pmapcds%), Nr.total=$n_readid, Lrdlen=$ardlen, MapErr=$nberr ($pberr%), Nshort=$n_rdtooshort\n"; 
  #UPD18mar: add new LN.reads=srdlen_tot, others above
  #my $srdlen_tot= $srdlen+$srdlen_nomap; 
  #my $nrdlen_tot= $n_mapok + $n_nomap + $n_mapnotucg + $n_mapbad; #  nrdlen_tot should == n_readid **

  return($topinfo,$ucmwMedn,$ucmwAve,$ucmwSE);
}

=item samtools flagstats

117824819 + 0 in total (QC-passed reads + QC-failed reads)
80316418 + 0 primary   << n_reads
36954853 + 0 secondary  << n_dupmap? yes
553548 + 0 supplementary  << frag map from primary, dont count
0 + 0 duplicates
0 + 0 primary duplicates
66227238 + 0 mapped (56.21% : N/A)   << n_dupmap = mapped - primap - supmap = 66227238 - 28718837 - 553548 = 36954853
28718837 + 0 primary mapped (35.76% : N/A)  << n_mapok
80316418 + 0 paired in sequencing
40158209 + 0 read1
40158209 + 0 read2
20529252 + 0 properly paired (25.56% : N/A)
21757004 + 0 with itself and mate mapped
6961833 + 0 singletons (8.67% : N/A)
886502 + 0 with mate mapped to a different chr
600097 + 0 with mate mapped to a different chr (mapQ>=5)

=cut

sub readmapstats {
  my($inf)=@_;
  my($ok,$inh)= openRead($inf);
  return(0) unless($ok);
  my($nrd,$nmap,$ndupmap)=(0,0,0);
  while(<$inh>){
    my($n,@v) = split;
    if(m/ primary$/ and not $nrd){ $nrd=$n; }
    if(m/ secondary$/ and not $ndupmap){ $ndupmap=$n; }
    if(m/ primary mapped / and not $nmap){ $nmap=$n; }
    # if(m/xxx/){ $RDLEN=$n; }
  } close($inh);
  #?? sum if several mapstats?
  ($nrdlen_tot,$n_mapok)=($nrd,$nmap);
  $n_readid= $nrd;
  #note: n_uniqmap is NOT n_map - n_dupmap, dupmap has all dups of each primap w/ any dup
  $n_nomap= $nrdlen_tot - $n_mapok;
  $nmapid= $n_mapok; # should be same as $n_mapok, 2nd above, be sure
  
  $srdlen= $n_mapok * $RDLEN;
  $srdlen_nomap= $n_nomap * $RDLEN;
  return($nrdlen_tot,$n_mapok);
}

sub ratio {
  my($va,$vb)=@_; 
  my($r,$inv)=(0,0);
  if($va > $vb){ ($va,$vb,$inv)=($vb,$va,1); }
  if($vb < 0.01){ $r= 0; } else { $r= $va/$vb;  }
  return (wantarray)? ($r,$inv) : $r;
}


sub read_idclass {
  my($crclassf,$allclasses,$crlenh)=@_; 
  my($nid,$ncl,$ok,$inh)=(0,0,0);
  my %crclass=(); my @crclass=();
  $allclasses||=0;

  use constant kMAXIDCLASS => 99999; # idclass limit 9, using 3-4 now : ONLY FOR covtab cols .. not for genexcopy
  my $CRTPAT=''; # no default, see sam2covtab
  
  if($crclassf and ($ok,$inh)= openRead($crclassf) and $ok) { 
    while(<$inh>){ next if(/^\W/); 
      my($cr,$crclass,@clx)=split;  # may have more columns .. keep all ie CDS,BUSCO 
      if($allclasses and @clx){ $crclass=join(" ",$crclass,@clx); }
      $crclass{$cr}= $crclass || 0;  
    } close($inh); 
    $nid= scalar(keys %crclass);
  }
  
  unless($nid>0) { @crclass=('UNK'); }
  else {
    my %crcl=(); # make list of values = classes
    map{ $crcl{$crclass{$_}}++; } keys %crclass;  # or values %crclass > not uniq
    @crclass= sort keys %crcl; # sort by count?
  }
  
  $ncl= @crclass; 
  if($ncl > kMAXIDCLASS) { # problem, cancel..
    warn "#ERR too many idclasses n=$ncl from nid=$nid of $crclassf  \n";# or sam ids x CRTPAT='$CRTPAT'
    $nid=0; %crclass=(); @crclass=();
  }
  warn "# read nid=$nid, nclass=$ncl from $crclassf\n" if($debug);
  return($nid,\%crclass,\@crclass);
}  


sub openRead {
  my($fn)=@_; my($ok,$inh)=(0); 
  return(0,undef,$fn) unless($fn and -f $fn);
  if($fn =~ /\.gz/) { $ok= open($inh,"gunzip -c $fn |"); } else { $ok= open($inh,$fn); }
  return($ok,$inh,$fn);
}

__END__



