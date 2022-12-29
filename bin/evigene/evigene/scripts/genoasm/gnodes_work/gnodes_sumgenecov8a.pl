#!/usr/bin/perl
# gnodes_sumgenecov.pl

=item usage

  gnodes_sumgenecov.pl  -chrgenetab xxx_chr.genetab -genexcopy xxx_cds.genexcopy 
     -cdscov?? xxx_cds.covtab -anngenes ?? xxx_ann.genecopyn
     -out xxx_genecopy.sumtab?
    -idclass xxx.idclass -anntab xxx.anntab 

  from gnodes_genescov.pl, replacing inputs, outputs
  
  dropse20gnodes/
  $evigene/scripts/genoasm/gnodes_genescov.pl -debug  \
    -title dropse20cdschr -asmid dropse20chrs -metad dropse20chrs.metad \
    -idclass dropse20cdste.idclass -ann dropse20chrs.fa.anntab \
    -cdscov dropse20t1cds_SRR11813283_1_bwa.cdschr7b.covtab -chrcov dropse20chrs_SRR11813283_1_bwa.cdschr7b.covtab
    
  $evigene/scripts/genoasm/gnodes_sumgenecov.pl -title at18chr_test8f2sumgcn  \
   -genexcopy arath18tair1cds_SRR10178325_b2_mim.geneycopymd -chrgenetab arath18tair_chr_SRR10178325_test8f.genetab \
   -asmid arath18tair_chr -meta arath20asm.metad -idclass arath18tair1cds.idclass 
 
  env transform=log notabout=1  $evigene/scripts/genoasm/gnodes_sumgenecov.pl \
   -title at20max_test8f3lsumgcn  \
   -genexcopy arath18tair1cds_SRR10178325_b2_mim.geneycopymd -chrgenetab arath20max_chr_SRR10178325_test8f.genetab \
   -asmid arath20max_chr -meta arath20asm.metad -idclass arath18tair1cds.idclass 

=item UPD: add Deficit table/plot summary

  from gnodes_sumasmgcn2rtab.pl
  inputs: this.genesum.txt and gnodes_covsum chr asm xxx_sum.txt
   
  output: R format table, or R plot script, of chr-asm deficits, both on chr-spans and gene copynums
  SppAsm        Item    pXlo    pXeq    pXhi    pMiss
  ara18ch       Aall    -32.77  67.11   119     177
  ara18ch       Auniq   -14.12  81.30   108     133
  ara18ch       Adup    -18.76  25.77   11.6    44.8
  ara18ch       Acds    -12.99  68.49   50      73
  ara18ch       C01     0       95.43   0.37    0
  ara18ch       C02-9   -1.17   1.67    0.11    0
  ara18ch       C10-99  -0.06   0.15    0.04    -0.01
  ara18ch       Cmiss   0       100     1259    0
  --- or --
  # Rscript
  acdeficitplot <- function(vd, plname="asmgcn_miss_",ymin=-60, gscale=4) { ... barplot() }
  vd=read.table(header=T, text="
  SppAsm Item pXlo pXeq pXhi pMiss
  ara18ch Aall -32.77 67.11 119 177
  ara18ch Auniq -14.12 81.30 108 133
  ara18ch Adup -18.76 25.77 11.6 44.8
  ara18ch Acds -12.99 68.49 50 73
  ara18ch C01 0 95.43 0.37 0
  ara18ch C02-9 -1.17 1.67 0.11 0
  ara18ch C10-99 -0.06 0.15 0.04 -0.01
  ara18ch Cmiss 0 100 314 0
  ");
  acdeficitplot(vd);
   
=cut    

use strict;
use Getopt::Long;  

use constant UPD21AUG => 1; # fixes, add aLen, assembly length vs gLen as other copy num
use constant UPD21SEP => 1;
use constant LLOCI_CALC => 2; # UPD21SEP16, LLOCI_CALC==2 now best match to gCopy, swap lloci_calc2 for calc1 (both out for now)

my $debug=$ENV{debug}||1;
my @IVAR=(5,6,7); # covtab cols ==  rd_total, rd/mmap, rd_uniq for chrasm x rd map 
my $DOMED= 0; # $ENV{median}||0; * SUCKS UP ALL covtab vals, mem pig
my $DOHIST=1; # replace median from @allvals w/ histo mode(s)/peak, est median, counts for depth vals depth[0..999]?
my $ZEROS= 0; # $ENV{zero}||0;
my $showSUM=1; # $ENV{sum}||0;
my $NTAVE= $ENV{ntave} || 0;
my $MINRD= 1; #? skip all cov rows w/ acovt < MINRD, want to measure zeros? gaps?
my $CBIN=100; # calc from rows ib dist
my $TOTALID="total"; #? "allgenes";
my $KU_OPT=0; # read from metad, or calc from genemeans, unless -CU|kucg=val
my $XHICUT= $ENV{xhicut}|| 1.7; my $XLOCUT= $ENV{xlocut}|| 0.55;# xCopy opts, change def? xhi=1.60 .. 1.66? xlo=0.50..0.54?
my $TRIMAVE= $ENV{trimave}||0; # for ave calc, trim extremes?

my $STOPTS=$ENV{'statopt'} || ""; # opt for sumgene_sumout : nozero, trim
my $TRANS= $ENV{'transform'}||""; # transform for sumgene_sumout log, sqrt, ..
my $NOTABO= $ENV{notabout} || 0;
my $Genetab8i=$ENV{'genetab8i'} || 0; # UPD21AUG gnodes_sam2covtab8i.pl

use constant kSAMPLE => 1000;  # CBIN size median sample @cbin, should be all same = default $CBIN

use constant UCG_KU_SPECIAL => 1; # special KU cov-depth calc for uniq-conserved-genes: median of median depth
my($nucg,@ucg_med_cov,@ucg_med_ids,@ucg_one_cov,%ucg_id_cov)=(0); # UCG_KU_SPECIAL globals

my $intitle= $ENV{title}||$ENV{name}||""; # was "asource"; 
my $asmid= $intitle; #? intitle for output, asmid for data
my ($sampledata,$anntable,$outmeans,$outsum,$ncpu,$crclassf,$readsetid)=("") x 9;
my ($cdscovtab,$chrcovtab,$genexcopy); # input data
my (@lvar, @genextop, @CBIN);


my $optok= GetOptions( 
  'cdscovtab=s',\$cdscovtab, 
  'chrgenetab=s',\$chrcovtab, # was chrcovtab=
  'genexcopy=s',\$genexcopy, # from sam2genecov
  'output|outsummary=s', \$outsum,
  'means=s', \$outmeans,
  'title|name=s',\$intitle,
  'asmid=s', \$asmid, # superceeds intitle for finding/using asm info
  'readsetid=s', \$readsetid, # upd21f20
  'genomedata|sumdata|sampledata|metadata=s', \$sampledata, # optional input  
  'anntable|annotations=s', \$anntable,
  'crclassf|idclassf=s',\$crclassf, # table of chr/gene id => class ? want for BUSCO, TEfam, : No, rely on anntable having these
  'minread=i', \$MINRD, # 5 default, test, see if it affects over-assembly (Xcopy < 1)
  'CU|kucg=s', \$KU_OPT,
  'XHICUT=s', \$XHICUT, 'XLOCUT=s', \$XLOCUT, # xCopy hi/lo cutoff from 1.0
  'statopt=s', \$STOPTS, 'transform=s', \$TRANS,
  'NOTABOUT!', \$NOTABO, 
  'ncpu=i', \$ncpu, # not used here
  'debug!', \$debug, 
  );

die "usage: gnodes_sumgenecov.pl -chrgenetab xxx_chr.genetab  -cdscovtab xxx_cds.covtab -anntab xxx.anntab  -name xxxcdschr 
 opts: -metadata xxx.metad -asmid xxx_asm -idclass xxx.idclass -CU=0 -minread=$MINRD -debug 
" unless($optok and ($cdscovtab or $chrcovtab or $outmeans));

if($asmid and not $intitle) { $intitle=$asmid; } # want one of these
elsif($intitle and not $asmid){ $asmid=$intitle; }
unless($readsetid) {
  my $rdset=""; # FIXME leads to messy fname if no SRRnnnn in covtabs
  for my $incov ($cdscovtab,$chrcovtab,@ARGV) {
    if($incov and $incov =~ m/(\w+).(\w+)\.covtab/) { 
      my($asmrd,$vers)=($1,$2); 
      $asmrd =~ s/($asmid|$intitle)//g;
      if($asmrd =~ m/([A-Z]RR\d+)/) { my $na=$1;  $rdset .= $na unless($rdset=~m/$na/); last; }
      elsif($asmrd=~/([A-Za-z0-9]+)/){ my $na=$1;  $rdset .= $na unless($rdset=~m/$na/); last; }
    }
  }
  if($rdset =~ m/\w\w/){ $rdset=substr($rdset,0,19) if(length($rdset>19)); $readsetid=$rdset; }
}

my $title= $intitle;
unless($outsum) { 
  ($outsum=$intitle) =~ s/\W/_/g; # title always ok?
  $outsum.="_".$readsetid if($readsetid and $outsum !~ m/$readsetid/);
  $outsum.="_genesum.txt";
}
unless($outmeans){ ($outmeans=$outsum) =~ s/.genesum.*//; $outmeans.=".genemeans"; }
#add here or putGeneCovsum: outxcopy=outname.genexcopy

warn "#genecov output to sum=$outsum\n" if $debug; #  means=$outmeans

sub MAINstub {}

my($nsamp,$masmid,%samvals); # masmid == asmid nor not?
($nsamp,$masmid)= readMetad($sampledata);

my($nann,%annvals);
($nann)= ($anntable) ? readAnnots($anntable) : (0); # format from gnodes_annotate.pl, same crID, crLoc, an, ids as covtab

my($nidclass,$idclassh,$idclasslist)= ($crclassf) ? read_idclass($crclassf, 1) : (0); 

use constant { GENEVEC => 0, GENEHASH =>1, };
use constant { STORECOV => 1, MINCOVT => 1, }; # readGenetab MINCOVT 1..9 ? skip 0, tiny covs as spurious
# read_genexcopy global gene copynum table
my($gucg_n, $gucg_med, $gucg_ave, $gucg_err, $ngenecov, $genecovh) =  
      ($genexcopy)? read_genexcopy($genexcopy): (0,0,0,0,0); 

my($ntcds,$cdslen,$cdsrdmap)= ($cdscovtab) ? readChrtab($cdscovtab, 1) : 0;  # drop, only for cdslen in genecov{id}{Glen}

my($nccds)= ($cdscovtab) ? readCovtab($cdscovtab, 1) : 0; #? useful w/ genexcopy ?

my($ncchr)= ($chrcovtab) ? readGenetab($chrcovtab, 0) : 0; # primary data + genexcopy

#...... new outputs .......
my $outnam=$outsum; $outnam=~s/.genesum.*//;

#UPD: special case genecovh.gClass $isuniq eq "zero", store? or put aside? here or output?
my @gids= sort keys %$genecovh;

# add gLen < $MINLEN here?
#   use constant { gMINLEN => 50, gMINREAD => 3 }; #?? genexcopy drop these filts? or move to outputs?

my @gidnoz= grep{ not ($genecovh->{$_}->{gClass} =~ m/zero|^0/ or $genecovh->{$_}->{gNread} < $MINRD ) } @gids;
my @gidzero= grep{ ($genecovh->{$_}->{gClass} =~ m/zero|^0/ or $genecovh->{$_}->{gNread} < $MINRD ) } @gids;
my $ngzero= @gidzero; 

#UPD21AUG20, aLen/aBins output: aLen= aBins * $CBIN ** This is over-estimate from CBIN slop, any correction?
for my $gid (@gids) { my $nb= $genecovh->{$gid}->{'aBins'}||0; $genecovh->{$gid}->{'aLen'}= $nb*$CBIN; }

#UPD21SEP08 ** Output cols add lCopy, lLoci .. lCopy2 testing SEP14
# 21SEP16: swap lCopy1/2, calc2 is best now .. old is lCopy2 now, drop

#UPD21SEP22: ? add new est lCopy * lCov/CU , should be closest to gCopy if measures are accurate
# .. also swap out aCopy of lCopy, lxCopy=lCopy*lCov/CU?, for xCopy in Copynum Summary by Copy level ?
use constant USE_lCopy_xCopy => 1; #UPD21SEP22: add lxCopy, calc xCopy from lCopy/gCopy, and? lxCopy/gCopy     
my $CUest= $KU_OPT || $gucg_med;
if(USE_lCopy_xCopy  and $CUest >= 1) {
  for my $gid (@gids) { 
    my $lcopy= $genecovh->{$gid}->{'lCopy'}||0; 
    my $lcov= $genecovh->{$gid}->{'lCov'}||0;  # FIXME: lcov is tuple: 407.9,31,100%, use median=2
    my($lcova,$lcovm)=split",",$lcov; $lcovm||=$lcova;
    my $lxcopy= sprintf "%.1f", $lcopy * $lcovm / $CUest; # dang, not sure this is * or +
    $genecovh->{$gid}->{'lxCopy'} = $lxcopy; # ugh, output in sumgene_tabout(), count in sumgene_counts()
    }
}

# aCovZ == 0 for all? check
# my @FLDOold= qw( gLen gNread gCopy gCM gCnz gClass  aCopy aCovT aCovM aCovZ  mCopy mCovT aCopyRd aReadT aReadM aReadZ);
my @FLDO= qw( gLen gNread gCopy gCM gCnz gClass aCopy aCovT aCovM aCopyRd aReadT aReadM aReadZ aLen lCopy lCov  lLoci );
push @FLDO, 'lxCopy' if(USE_lCopy_xCopy);
if(exists $genecovh->{$gids[1]}->{'lCopy2'}) { push @FLDO, qw(lCopy2 lCov2 lLoci2); } #? lCov/lCovM ?

unless($NOTABO) {
my($ntabo,$tabo)= sumgene_tabout(\@gids, \@FLDO,  $outnam.".sumgenetab");
  warn "#sumgene_tabout: n=$ntabo (nzero=$ngzero) $tabo\n" if $debug;
}

#my $STOPTS=$ENV{'statopt'} || ""; # opt for sumgene_sumout : nozero, trim
#my $TRANS= $ENV{'transform'}||""; # transform for sumgene_sumout log, sqrt, ..

# my @FLDSUMold= qw( gLen gNread gCopy aCopy mCopy aCopyRd  gCM aCovT mCovT aReadT); #?  aCovZ or aReadZ?
my @FLDSUM= qw( gLen gNread gCopy aCopy aCopyRd  gCM aCovT aReadT aReadZ aLen lCopy lCov); #?  aCovZ or aReadZ?
my %FLDCOR= ( gCopy => "gLen,gNread,lCopy", #? lCopy: yes
    aCopy => "gCopy,gLen,gNread,lCopy", aCopyRd => "gCopy,gLen,gNread,lCopy", 
    aCovT => "gCM,gLen,gNread",  aReadT => "gCM,gLen,gNread",
    ); #  mCopy => "gCopy,gLen,gNread", mCovT => "gCM,gLen,gNread",

my($stitle,$nsum)=(0,0);
# $outsum="$outnam.genesum.txt";
my($ook,$outsumh)= openOut($outsum);

# FIXME2: add gClass counts table, in sumgene_counts()? or end w/ @genextop? ie: uniq:50%,12000 dupx:25%,6000 ..
# FIXMEd: gCopy cuts for counts should ~ match gCopy cuts for moments
use constant { C0hi=>0.66, C1hi=>1.55, C2hi=>9.99, C3hi=>99.9, C4hi=>499 }; # xCopy = aCopy/gCopy cuts

## insert here: Counts for genes cover depth on chromosome assembly
$stitle="# Chromosome Assembly Gene-Copynum Counts by Copy level, Chrasm=$asmid";
my($nco,$cnbrief)= sumgene_counts(  $outsumh, $stitle, \@gids, \@gidnoz, \@gidzero);
$nsum+=1;  warn "#$cnbrief\n" if($debug);

$stitle="# Moments for all genes";
$nsum += sumgene_sumout( \@gidnoz, \@FLDSUM, $outsumh, $stitle, \%FLDCOR); # call several times, one output file

my $gidsel= sumgene_select( \@gidnoz, { gCopy => [C1hi,C3hi] }); # 1.75,99
# $stitle= "# Moments for gCopy => [1.75,99] genes";
$stitle= join" ","# Moments for gCopy => [", C1hi,C3hi,"] genes";
$nsum += sumgene_sumout( $gidsel, \@FLDSUM, $outsumh, $stitle, \%FLDCOR); # call several times, one output file

my $gidsel= sumgene_select( \@gidnoz, { gCopy => [C0hi,C1hi] }); # 0.75,1.50
$stitle= join" ","# Moments for gCopy => [", C0hi,C1hi,"] genes";
$nsum += sumgene_sumout( $gidsel, \@FLDSUM, $outsumh, $stitle, \%FLDCOR); # call several times, one output file

my $gidsel= sumgene_select( \@gidnoz, { gCopy => [0.0,C0hi] }); # 0.0,0.50 include gClass =~ zero here? or no
$stitle= join" ","# Moments for gCopy => [", 0.0,C0hi,"] genes";
$nsum += sumgene_sumout( $gidsel, \@FLDSUM, $outsumh, $stitle); # call several times, one output file

if(@genextop){
  print $outsumh "\n# Uniq Gene Coverage summary of $genexcopy \n"; 
  map{ print $outsumh $_,"\n" } @genextop;
}
close( $outsumh);
warn "#sumgene_sumout: n=$nsum $outsum\n" if $debug;

# if($nccds > 0 or $ncchr > 0) {
#   putGeneCovStats($outmeans); # revise, use genecovh table of data
#   %ucg_id_cov=(); # empty hash of all data, use outmeans
# }
# 
# # step2 usage: read outmeans, write outsum, dont need covtabs
# if(-s $outmeans) {
#   putGeneCovSum($outsum,$outmeans);  # drop or rewrite, ? want summary separate from genecov.tab
# }

#=======================================================



=item sumgenecov stats

genecov hash fields:
  UCG subset: gCM gCnz gCopy=1 gClass=uniq nGene
  gCopy gCM gCnz gLen gNread gClass  : from cds_gdna.genexcopy .. add UCG flag, to gClass?
  aCopy aCovT aCovM aCovZ      : ave from chr_cds_gdna.genetab
  mCopy mCovT mCovM            : mdn from  ""
  aCopyRd aReadT aReadM aReadZ : readgene count from chr_cds_gdna.genetab

output genecovsum table?
  gID gLen gNread gCopy gCM gCnz gClass,  aCopy aCovT aCovM aCovZ, 
    mCopy mCovT?, aCopyRd aReadT aReadM aReadZ ?
  
output genecovsum summary:
  ave, medn for a. all, b. gCopy =~ 1, c. gCopy >=2 .. <10 or 100, d. gCopy < 0.5 ~ zero, e. gCopy>99 or what
  log(ave) << option to switch all output to transformed vals? ie gene.tab + gene.sum after transform
  skew stats?
  cor(gCopy,aCopy) .. ,mCopy,aCopyRd .. gLen, gNread, [amr]CovT? CovM?

 -----
  # global($gucg_n, $gucg_med, $gucg_ave, $gucg_err, $ngenecov, $genecovh)= read_genexcopy(infile)
  $gcovtab{$gid}= { gCopy => $tcopy, gLen => $glen, gNread => $nread, gClass => $isuniq,
                         gCM => $cmed, gCnz => $cnz, };
  $genecovh->{$gid}->
 {'aCovT'} = $ct; {'aCovM'} = $cm; {'aCovZ'} = $cno; {'aCopy'} = $pca;
 {'mCovT'} = $mt; {'mCovM'} = $mm; {'mCopy'} = $pcm;
 {'aReadT'} = $mt; {'aReadM'} = $mm; {'aReadZ'} = $mno; #? is it valid
 {'aCopyRd'} = $pcm;

=cut

=item sumgene_select

  %selfield= ( gCopy => [1.75,99], gLen => [300,3000] .. )
  \@gidsel= sumgene_select(\@gids, \%selfield)
  
=cut 

sub sumgene_select {
  my($gids, $selfields)= @_;

  my($nsel)=(0); my @selid= ();
  my @fld= sort keys %$selfields;
  # mainly for simple case: @fld == 1
  if(@fld == 1) {
    my $f=$fld[0];
    #a my($slo,$shi)= $selfields->{$f}->[0,1]; # bad
    my($slo,$shi)= @{$selfields->{$f}}[0,1]; # ok
    @selid= grep{  $genecovh->{$_}->{$f} >= $slo and  $genecovh->{$_}->{$f} <= $shi  } @$gids;
  } else {
    my %fsel;
    for my $f (@fld) {
      # my @srange= @{ $selfields->{$f} };  $fsel{$f}= \@srange;
      #bad my($slo,$shi)= $selfields->{$f}->[0,1];
      my($slo,$shi)= @{$selfields->{$f}}[0,1]; #?? ok
      $fsel{$f}[0]= $slo;  $fsel{$f}[1]= $shi; 
    }
    
    for my $gid (@$gids){ 
      my($isok,$isbad)=(0,0);
      for my $f (@fld) {
        my $v= $genecovh->{$gid}->{$f} || 0;
        if($v >= $fsel{$f}[0] and $v <= $fsel{$f}[1]) { $isok++; } else { $isbad++; }
      }
      push @selid, $gid if($isok and not $isbad);  
    } 
  }
  $nsel= @selid;
  return($nsel,\@selid);  
}

=item counts: gClass added

at18chr_test8f8sumgcn_genesum.txt
 Gene Classes of
   uniq=1c, dupx=2+c, dups=partial-dup, skew=uneven coverage, zero=below minimal cover
 Subset Ngene   uniq            dupx            dups            skew    zero
 valid  27347   92.3%,25246     1.3%,350        6.3%,1713       0.1%,38 0.0%,0
 all    27444   92.0%,25246     1.3%,350        6.3%,1718       0.1%,38 0.3%,92
 ---------------------------------------

daphpulex_pa42v2_test8f8sumgcn_genesum.txt
 Gene Classes of daphplx17evgt1m_cds_SRR13333791_b8_mim.geneycopy
  uniq=1c, dupx=2+c, dups=partial-dup, skew=uneven coverage, zero=below minimal cover
 Subset Ngene   uniq            dupx            dups            skew            zero
 valid  32510   55.2%,17943     15.6%,5059      26.2%,8507      3.1%,1001       0.0%,0
 all    33037   54.3%,17943     15.3%,5062      25.9%,8555      3.0%,1004       1.4%,473
 ---------------------------------------

=item sumgene_counts
  
  $ncnt= sumgene_counts(  $outsumh, $stitle, \@gids, \@gidnoz, \@gidzero);

 FIXME2: add gClass counts table, in sumgene_counts()?  ie: uniq:50%,12000 dupx:25%,6000 ..
  dapplx19ml:  8555 dups 5062 dupx 1004 skew 17943 uniq 473 zero

  counts/gCopy as from gnodes_annotate 
  
  ==> arath18tair_chr-arath18tair1cds.genecopyn <==
  # Genes Copynum Summary by Copy level, for nGene=27359 found, 3 missed
  # tCopy=Chrasm, gCopy=Gene set, xCopy=tCopy/gCopy, x0,xLo,xEq,xHi= xCopy 0,<1,=1,>1
  #CnGene nGene   tCopy   xCopy   gCopy   x0,xLo,xEq,xHi
  #0      713     1.0     2.4     0.5     0,0,85,628
  #1      26212   1.0     1.0     1.0     1,206,25573,432
  #2-9    401     1.6     0.6     3.2     1,302,90,8
  #10-99  23      1.8     0.1     22.7    1,22,0,0
  #99-499 6       1.3     0.0     160.9   0,6,0,0
  #500+   7       4.2     0.0     831.5   0,7,0,0
  #---------------------------------------

# report these slices
my @gidnoz= grep{ not ($genecovh->{$_}->{gClass} =~ m/zero|^0/ or $genecovh->{$_}->{gNread} < $MINRD ) } @gids;
my @gidzero= grep{ ($genecovh->{$_}->{gClass} =~ m/zero|^0/ or $genecovh->{$_}->{gNread} < $MINRD ) } @gids;

# my @FLDOold= qw( gLen gNread gCopy gCM gCnz gClass  aCopy aCovT aCovM aCovZ  mCopy mCovT aCopyRd aReadT aReadM aReadZ);
my @FLDO= qw( gLen gNread gCopy gCM gCnz gClass aCopy aCovT aCovM aCopyRd aReadT aReadM aReadZ);
  
=cut

sub sumgene_counts {
  my($ctabh, $title, $gidall, $gidnoz, $gidzero, )= @_;

  use constant { C0hi=>0.66, C1hi=>1.55, C2hi=>9.99, C3hi=>99.9, C4hi=>499 }; # xCopy = aCopy/gCopy cuts
  # FIXME: gCopy cuts for counts should ~ match gCopy cuts for moments
  use constant CTOOLOW => 0.1;
  my $nall= @$gidall; my $noz= @$gidnoz; my $nzero= @$gidzero;
  print $ctabh $title,"\n";
  print $ctabh "# N Genes valid=$noz, zero-cover=$nzero, all=$nall (zero:gClass=zero or gNread<$MINRD)\n";
  
  # putcn() and putcn_missing()
  my($ncnsum, $ncmiss, %cnsum, %gclass);
  for my $id (@$gidall) {
    my($glen,$gnread,$gcopy,$gcm,$gclass)= map{ $genecovh->{$id}->{$_}||0 } qw( gLen gNread gCopy gCM gClass); 
    $gclass{$gclass}{'gidall'}++; 
  }

  #UPD21SEP09: add lCopy, drop aLen for counts
  my $xcopyis= (USE_lCopy_xCopy) ? "xCopy=lCopy/gCopy, ": "xCopy=aCopy/gCopy, ";
  for my $id (@$gidnoz) {
    my($glen,$gnread,$gcopy,$gcm,$gclass,$acopy,$acovt,$acopyrd,$areadt,$areadm,$areadz,$alen,$lcopy,$lcov,$lxcopy)= 
      map{ $genecovh->{$id}->{$_}||0 } 
      qw( gLen gNread gCopy gCM gClass aCopy aCovT aCopyRd aReadT aReadM aReadZ aLen lCopy lCov lxCopy);  

    my($xcopy,$xcopy2);
    if(USE_lCopy_xCopy) { # also add lxCopy
    $xcopy = ($gcopy<0.01)? 0: $lcopy/$gcopy;  
    $xcopy2= ($gcopy<0.01)? 0: $lxcopy/$gcopy; #?? use this also nor not
    } else {
    $xcopy= ($gcopy<0.01)? 0: $acopy/$gcopy; $xcopy2= $xcopy; 
    }
    
    #UPD21AUG20, add aLen output ?? in CLevel table: gCopy xCopy aCopy xAlen?/agLen ?
    my $acopylen= ($glen<1) ? 0 : $alen / $glen; 
    
    #? count gclass here? not per gcla as that == uniq/dupx gclass
    $gclass{$gclass}{'gidnoz'}++; 
    
    my $gcl= ($gcopy < C0hi)?"0" :(C0hi<=$gcopy and $gcopy<=C1hi)?"1" 
      :(C1hi<$gcopy and $gcopy<=C2hi)?"2-9" :(C2hi<$gcopy and $gcopy<=C3hi)?"10-99"
      :(C3hi<$gcopy and $gcopy<=C4hi)?"100-499":"500+";

    # $ceq= "xzero"; when aCovT or aReadT below MINRD?
    my $ceq= ($xcopy<C0hi)? "xlo" : ($xcopy>C1hi)?"xhi": "xeq";  
    if($areadt < $MINRD or $acopyrd < CTOOLOW) {  # MINRD = 1 default, want higher min here?
      $ncmiss++; $ceq="xzero";
    }
    
    #? add gcl="All" ?
    $ncnsum++;
    for my $gcla ($gcl,"-1") { # -1 == All for sort
      $cnsum{$gcla}{ngene}++; $cnsum{$gcla}{$ceq} ++; 
      $cnsum{$gcla}{gcopy}+= $gcopy; $cnsum{$gcla}{xcopy}+= $xcopy; 
      $cnsum{$gcla}{acopy}+= $acopy; $cnsum{$gcla}{lcopy}+= $lcopy; $cnsum{$gcla}{lxcopy}+= $lxcopy; 
      $cnsum{$gcla}{aglen}+= $acopylen;
      # add areadz count, if > minz
      $cnsum{$gcla}{missrd}+= $areadz;  $cnsum{$gcla}{areadm}+= $areadm; 
      $cnsum{$gcla}{miss5p}++ if( $areadz > 0.049*$areadm);
      #? $cnsum{$gcla}{acopyrd}+= $acopyrd; # .. areadt ?
    }
  }
  
  # putcn_sum
  my $cnbrief="";
  use constant AlwaysX0 => 1; # ncmiss only show of x0, is too confusing.
  if(1) {
    my @gcl=sort{ $a<=>$b } keys %cnsum;  
    my $ngt=", for nGene=$ncnsum, Chrasm_id=$asmid"; 
    #o: $ngt .= " found, $ncmiss missed" if($ncmiss>0); # BUT ncmiss in ncnsum, dont need ncmiss == x0 val
    my @ceql= qw(xlo xeq xhi); 
    unshift @ceql, "xzero" if(AlwaysX0 or $ncmiss>0); # can confuse to have some w/ some w/o xzero: always add?
    my $labxloeqhi= (AlwaysX0 or $ncmiss>0) ? "x0,xLo,xEq,xHi= xCopy 0,<1,=1,>1"
      : "xLo,xEq,xHi= xCopy<1,=1,>1";
    (my $labxleh= $labxloeqhi) =~ s/=.*//;
    #D if(UPD21AUG){ $labxloeqhi .= ", agLen= aLen/gLen"; }
    if(UPD21AUG){ $labxloeqhi .= ", lCopy=Chr loci"; } #UPD21SEP09
    
    $cnbrief="Genes Copynum Summary$ngt, xCopy= ";
    print $ctabh " ---------------------------------------\n";
    print $ctabh " Genes Copynum Summary by Copy level$ngt\n"; #== title
    print $ctabh " aCopy=Chrasm, gCopy=Gene set, $xcopyis $labxloeqhi\n";
    print $ctabh " aMiss=Missed gene reads on chrasm, % of gene reads\n";
    print $ctabh " ".join("\t",qw(CLevel nGene gCopy xCopy aCopy lCopy aMiss), $labxleh)."\n";
    for my $gcl (@gcl) {
      my $ng= $cnsum{$gcl}{ngene} or next;  # ng includes xzero for COPYNUMzeros, x/t reduced by zeros
      my($xc,$tc,$gc, $misp, $misr, $lcopy, $lxcopy, $aglen)= map{ sprintf"%.1f",$cnsum{$gcl}{$_}/$ng } 
        qw( xcopy acopy gcopy miss5p missrd lcopy lxcopy aglen) ;
      my($ceqn)= join",", map{ $cnsum{$gcl}{$_}||0 } @ceql; # qw(xlo xeq xhi); #? add xzero

      my $mis="";
      if(1) { # pMissRD calc 
        my($stot,$smis)= map{ $cnsum{$gcl}{$_}||0 } qw(areadm missrd);
        my $prmis= ($stot<1)? 0: sprintf"%.2g",100*$smis/$stot; 
        $mis="$prmis%,$smis";
      } else {
        $misp=100*$misp; # pct n-genes w/ > 0.05 missed cds-reads
        $mis="$misp%,$misr"; # or just $misp;
      }
      
      my $gclo= ($gcl eq "-1")?"All" : ($gcl eq "100-499")? "99-499" : $gcl;
      my $lcopyo= $lcopy; 
      if(USE_lCopy_xCopy and $lxcopy and $lxcopy ne $lcopy) { #? lxcopy here also, only if better est of $gc ?
        if( abs($gc - $lxcopy) < abs($gc - $lcopy) ) { $lcopyo= "$lcopy,lx:$lxcopy"; }
      }
      
      print $ctabh " ".join("\t",$gclo,$ng,$gc,$xc,$tc,$lcopyo,$mis,$ceqn)."\n"; 
      $cnbrief .= $gclo."c:$xc/$ng, " unless($gcl eq "0" or $gcl > 100);
    }
    print $ctabh " ---------------------------------------\n\n";
    
    print $ctabh " Gene Classes of $genexcopy\n";
    print $ctabh "  uniq=1c, dupx=2+c, dups=partial-dup, skew=uneven coverage, zero=below minimal cover\n";
    # dapplx19ml:  8555 dups 5062 dupx 1004 skew 17943 uniq 473 zero
    #  $gclass{$gclass}{'gidnoz'} and 'gidall'
    my @gcl= qw( uniq dupx dups skew zero ); my $gchdr=0;
    for my $gnz (qw(gidnoz gidall)) {
      my @gcv=(); my $ng=($gnz eq 'gidnoz')?$noz:$nall;
      for my $gcl (@gcl) { 
        my $c=$gclass{$gcl}{$gnz}||0; my $p=sprintf"%.1f",100*$c/$ng;
        push @gcv, "$p%,$c";
      }
      my $gnzl=($gnz eq 'gidnoz')?"valid":"all";
      print $ctabh " ".join("\t","Subset","Ngene",@gcl)."\n" if(1>$gchdr++);
      print $ctabh " ".join("\t",$gnzl,$ng,@gcv)."\n";
    }
    print $ctabh " ---------------------------------------\n\n";
  }
  
  # warn "#$cnbrief table to copytab\n" if($debug);
  return($ncnsum, $cnbrief);
}

=item pMissRD calc

env st=sumgene perl -ne 'BEGIN{  $NOCB=$ENV{nocb}||0;
$ST=$ENV{st}||"readgene";
$gtot=$gzero=$gmis=$gmis5p=$stot=$smd=$smis=0; } if(/^.GeneID/){
putv() if($gtot); ($fn=$ARGV) =~ s/.genetab//; } @v=split; next
if(/^\W|^total/ or $v[1] ne $ST); splice(@v,2,0,1) if($NOCB);
($id,$cr,$cb,$ct,$cm,$cmis)=@v; $gtot++; $gzero++ if($ct<2); $gmis++
if($cmis>0); $gmis5p++ if($cmis > 0.05*$cm); $stot+=$ct; $smd+=$cm;
$smis+=$cmis; END{ putv(); } sub putv{ ($pnmis,$p5mis)=map{
sprintf"%.2g%%",100*$_/$gtot; } ($gmis,$gmis5p);  $prmis=($smd<1)? 0:
sprintf"%.2g%%",100*$smis/$smd; print join("\t",qw(Ngene ZeroNg MissNg
Miss05% pMiss05 Rtot Rmtot Rmiss pRmiss Source ))."\n" if(1>$hdo++);
print
join("\t",$gtot,$gzero,$gmis,$gmis5p,$p5mis,$stot,$smd,$smis,$prmis,$
fn)."\n"; $gtot=$gzero=$gmis=$gmis5p=$stot=$smd=$smis=0; } ' \
 a*SRR10178325_test8f.genetab

Ngene   ZeroNg  MissNg  Miss05% pMiss05 Rtot            Rmtot           Rmiss   pRmiss  Source
26204   186     2162    506     1.9%    64191360        17022822        41092   0.24%   dmag14bgi2vtop5k_SRR7825549b_test8f
26203   197     3452    881     3.4%    30886818        16957033        119167  0.70%   dmag15nwb2asm_SRR7825549b_test8f
26204   376     2744    1183    4.5%    26356054        15962570        1111262 7%      dmag19skasm_SRR7825549b_test8f
26204   112     802     143     0.55%   58983806        17058858        3623    0.021%  dmag20sk4maca20ok_SRR7825549b_test8f
--
Ngene   ZeroNg  MissNg  Miss05% pMiss05 Rtot            Rmtot           Rmiss   pRmiss  Source
32839   40      6766    626     1.9%    127557198       26280569        402232  1.5%    daphplx_gasm16ml_SRR13333791_test8d
32839   83      7741    739     2.3%    349520723       26271355        250092  0.95%   daphpulex_pa42v2_SRR13333791_test8d
32839   67      7145    791     2.4%    268200134       26173426        360974  1.4%    dplx20maca4pkr_dc_SRR13333791_test8d
--
Ngene   ZeroNg  MissNg  Miss05% pMiss05 Rtot            Rmtot           Rmiss   pRmiss  Source
27380           47      2       0.0073% 98750890        19193878        825     0.0043% arath18tair_chr_SRR10178325_test8f
27380   133     3313    928     3.4%    748127607       18816037        361008  1.9%    arath20max_chr_SRR10178325_test8f
--
Ngene   ZeroNg  MissNg  Miss05% pMiss05 Rtot            Rmtot           Rmiss   pRmiss  Source
13928   7       33      0       0%      41060338        9112263         63      0.00069% drosmel6ref_chr_SRR11460802ab_test8f
14329   0       210     0       0%      87416601        29252079        681     0.0023%  dropse20chrs_SRR11813283ab_test8d

=cut

=item sumgene_sumout

at18chr_test8f2sumgcn.genesum.txt
# ALL genes ----------------------
Field   Median,Mean,N,StDev,Range________
gLen    1005    1187    27380   1483    0,16203
gNread  273     415     27380   7208    0,846276
gCopy   1       1.28    27380   9       0,683
aCopy   1       1.10    27380   1.83    0,167
mCopy   1       1.16    27380   2.99    0,167
aCopyRd 1       1.09    27380   1.74    1,166
gCM     33      41.9    27380   478     0,46641
aCovT   24      27.8    27380   136     0,10026
mCovT   26      26.6    27380   41.5    0,3414
aReadT  432     3607    27380   273060  2,30409470
aCovZ   0       0       27380   0       0,0

# gCopy => [1.75,99] genes ----------------------
Field   Median,Mean,N,StDev,Range________
gLen    582     813     705     1086    78,5730
gNread  288     2112    705     30954   20,640809
gCopy   2.60    5.15    705     8.84    1.80,92.8
aCopy   2       4.02    705     7.10    1,38
mCopy   2       6.07    705     17.4    1,167
aCopyRd 2       3.69    705     6.26    1,24.1
gCM     48      165     705     1629    12,39122
aCovT   46      104     705     484     7,10026
mCovT   45      74.4    705     112     1,498
aReadT  1716    85816   705     1476474 87,30409470
aCovZ   0       0       705     0       0,0

-- add corof( aCopy, gCopy) .. select pairwise cors
-- add parts of this? esp x0,xLo,xEq,xHi ?
==> arath20max_chr-arath18tair1cds.genecopyn <==
# Genes Copynum Summary by Copy level, for nGene=27175 found, 176 missed
# tCopy=Chrasm, gCopy=Gene set, xCopy=tCopy/gCopy, x0,xLo,xEq,xHi= xCopy 0,<1,=1,>1
#CnGene nGene   tCopy   xCopy   gCopy   x0,xLo,xEq,xHi
#0      702     1.0     2.4     0.5     0,0,95,607
#1      26212   1.0     1.0     1.0     164,293,25273,482
#2-9    401     1.7     0.6     3.2     11,277,97,16
#10-99  23      2.2     0.1     22.7    1,22,0,0
#99-499 6       2.7     0.0     160.9   0,6,0,0
#500+   7       8.9     0.0     831.5   0,7,0,0
#---------------------------------------

=cut

sub sumgene_sumout {
  my($gids, $fields, $outh, $title, $corto)= @_;

  my($no)=(0); 
  my @fld= @$fields;
  $outh= *STDOUT unless(defined($outh));
  
  # my $ofn= $outna.".sumgenesum"; #?? name suf?
  # my($ook,$outh)= openOut($ofn); return unless($ook);
  # outh param for repeat calls
  my $coradd= ($corto and ref($corto));
  my $corhd= ($coradd) ? "\tCorrelations" : "";
  
  print $outh "$title\n" if($title);
  print $outh "StatOpts=$STOPTS, Transform=$TRANS, \n" if($STOPTS or $TRANS);
  # print $outh "Field\tMedian,Mean,N,StDev,Skew,Range$corhd\n";
  print $outh join("\t",qw(Field Median Mean N StDev Skew Range))."$corhd\n";
  for my $fld (@fld) {
    my @val= map{ $genecovh->{$_}->{$fld} || 0 } @$gids;
    my @stat= meanvar_p( $STOPTS, $TRANS, \@val); # $med,$ave,$nit,$se,$range
    
    if($coradd and my $fldto= $corto->{$fld}) {
      my @fldb=split",", $fldto;  # maybe 3-5 fldb ?
      my @cors=();
      for my $fldb (@fldb) {
        my @valb= map{ $genecovh->{$_}->{$fldb} || 0 } @$gids;
        my @corab= covar_p( $STOPTS, $TRANS, \@val, \@valb); # $cor,covar : pairwise? matrix for all cor vars?
        # corab = ($nit,$cor,$covar,$ave,$aveb,$sd,$sdb)
        push @cors, "r.$fldb:$corab[1]"; 
        if($debug > 1) { push @cors, "rd.$fldb:".join(",",@corab); } #?
      }
      push @stat, join(",",@cors) if(@cors);
    }
    
    print $outh join("\t",$fld,@stat)."\n"; $no++;
  }
  print $outh "\n";
  return($no);
}


sub sumgene_tabout {
  my($gids, $fields, $ofn, )= @_;

  my($no)=(0); 
  my @fld= @$fields;
  my($ook,$outh)= openOut($ofn); return unless($ook);
  print $outh join("\t","GeneID",@fld)."\n";
  for my $gid (@$gids){ 
    my @val= map{ $genecovh->{$gid}->{$_} || 0 } @fld;
    print $outh join("\t",$gid,@val)."\n"; $no++;
  } close($outh);
  return($no,$ofn);  
}


# my($ntcds,$cdssizes)= ($cdscovtab) ? readChrtab($cdscovtab, 1) : 0; # 
# .cdschr7b.chrtab
# #ChrID	chrlen	nmap	nuniq	nmult	nmread	nomap
# dropse20uc:g117183162t1	489	722	214	508	343	0
# dropse20uc:g117183163t1	3459	2816	1884	932	2088	0

sub readChrtab { # read cds.chrtab if avail for sizes, also chr.chrtab?
  my($covtab, $iscdscov)=@_; #? for both cds,chr covtabs?
  my($nt,%len,%rdmap)=(0);
  (my $chrtab=$covtab) =~ s/covtab/chrtab/;
  open(my $inh,$covtab) or return 0;
  while(<$inh>) {
    next if(/^\W/);
    my @v=split; my($id,$clen,$nrmap,$nuniq,$nmult,$nmread)=@v;
    $len{$id}=$clen; $nt++;
    $rdmap{$id}=[$nrmap,$nuniq,$nmread]; # or all?
  } close($inh);  
  warn "#readChrtab($chrtab) nok=$nt\n" if $debug;
  return($nt,\%len,\%rdmap);
}


sub read_genexcopy {
  my($genexcopy)=@_;
  my($nucg,$cmd,$cav,$cerr,$ndup,$nuni,$ngc,$nzero,$topsum)= (0) x 9;
  my(%gcovtab);

  use constant { gMINLEN => 50, gMINREAD => 3 }; #?? genexcopy drop these filts? or move to outputs?
  
  my($ok,$inh)= openRead($genexcopy); return(0) unless($ok);
  $topsum=1;
  while(<$inh>) {
    if(/^\W/){
      # Uniq Gene Cov Depth n=954, .. C.Map/W=38.4, 38.0 +/-1.26 (mdn,ave,sem) 
      #UPD: save all of top sum stats for output
      if(/Uniq Gene Cov Depth n=(\d+)/ and not $cmd) {  
        $nucg=$1;
        # if( m/C.Map.W=([\d\.]+), ([\d\.]+)....([\d\.]+)/ ) { ($cmd,$cav,$cerr)= ($1,$2,$3);}
        if( m/C.Map.W=([\d\.]+), ([\d\.]+)\W+([\d\.]+)/ ) { ($cmd,$cav,$cerr)= ($1,$2,$3);}
        elsif( m/C.Map.W=([\d\.]+)/) { $cmd=$1; }      
        # last if($cmd);
      }
      if($topsum){ 
        s/^#//; chomp(); push @genextop,$_;
      }
    
    } else { # if(/^\w/) # Gene_ID or gene row, collect all gene rows
      # cols: Gene_ID Glen Nread Rdlen tCopy C.M Uniq Merr C.nz S.Dup,Unq,Nt,Equ,UCG
      my @v=split; $topsum=0;
      my($gid,$glen,$nread,$rdlen,$tcopy,$cmed,$isuniq,$merr,$cnz,$sdune)=@v;
      next if(/^Gene_ID/);
      #UPD: special case gClass $isuniq eq "zero", store? or put aside? here or output?
      my $GOK= 1; #was ($glen >= gMINLEN and $nread >= gMINREAD);
      
      if($GOK) { #  need to record nread == 0 cases ***
        if(GENEVEC) {
        $gcovtab{$gid}= [$tcopy,$cmed,$cnz,$glen,$nread,$isuniq];
        } else {
        # switch to named col vals, extend w/ other gene data
        $gcovtab{$gid}= { gCopy => $tcopy, gLen => $glen, gNread => $nread, gClass => $isuniq,
                         gCM => $cmed, gCnz => $cnz, };
        }
        $ngc++; 
        if($tcopy>=1.75){ $ndup++; } elsif($tcopy>=0.50){ $nuni++; } else { $nzero++; }
      } else {
        $nzero++;
      }
    }
    
  } close($inh);
  
  # maybe store in genecovh->{UCG}->{gCM} = $cmd, {gCopy}= 1, {gClass}='uniq', {gCnz}= "$cav,$cmd,span,pcov", {nGenes}= $nucg
  warn "#read_genescov: UCG med=$cmd, sd=$cerr, ave=$cav, n=$ngc\n" if $debug;
  warn "#read_genescov: ngene=$ngc, nzero=$nzero ($nuni 1-copy, $ndup 2+copy) from $genexcopy\n" if $debug;  
  # global($gucg_n, $gucg_med, $gucg_ave, $gucg_err, $ngenecov, $genecovh)= read_genexcopy(infile)
  return($nucg,$cmd,$cav,$cerr, $ngc, \%gcovtab); 
}

=item read_genexcopy

UPD 21jun11: now two C ests, first is best

# UCG Class Genes >= 1k long, 581 of 33037 total ----------------------
# Uniq Gene Cov Depth n=581, C.Map/W=39.6, 39.5 +/-1.67 (mdn,ave,sem) for W.genes=1073412, LN.reads=53677934
# Measured Unique Genes >= 1k long, 8686 of 33037 total ----------------------
# Uniq Gene Cov Depth n=8686, C.Map/W=39.1, 40.9 +/-0.53 (mdn,ave,sem) for W.genes=19138140, LN.reads=974367782

UPD21JUN from sam2genecov

Only want this info for pickKuni:
  # Uniq Gene Cov Depth n=954, .. C.Map/W=38.4, 38.0 +/-1.26 (mdn,ave,sem) 
  
 gnode8d/daphplx17evgt1m_cds_SRR13333791_b8_mim.genexcopy
#gnodes_sam2genecov options: minident=0.4, mindupid=0.98, maxdup=0.05,  idclass=daphplx17evgt1m_cds.idclass
# All  Gene Cov Depth n=33037, C.LN/W=63.4, C.Map/W=50.4 ave, for W.genes=39654210, LN.reads=2515206136
# Uniq Gene Cov Depth n=954, C.LN/W=49,  C.Map/W=38.4, 38.0 +/-1.26 (mdn,ave,sem) for W.genes=1337262, LN.reads=66509100
# Genome_size= 227.0 Mb of L*N/C= 150 * 58112590 / 38.4, CDS_size= 65.5 Mb,
# Gsize_alt  = 227.1 Mb of LN/C= 8718015246 / 38.4 for Lr: 150.0, Nr: 58112590 = 16765820+ 37506403+ 0+ 3840367 [ok,nomap,noucg,bad]
#   for Nr.ucgcds=16765820 (28.9%), Nr.allcds=16765820 (28.9%), Nr.total=58112590, Lrdlen=150, MapErr=61433951 (3.07%), Nshort=0
Gene_ID	Glen	Nread	Rdlen	tCopy	C.M	Uniq	Merr	C.nz	S.Dup,Unq,Nt,Equ,UCG
Daplx7pEVm000001t1	41844	12614	1892352	1.1	43	dups	0.9	47.9,41815,99.9	219299,1782013,2001416,1520188
Daplx7pEVm000002t1	36285	12090	1813741	1.2	46	uniq	1.8	46.8,36285,100.0	36898,1660739,1699097,1513236
Daplx7pEVm000003t1	26235	9394	1409241	1.2	48	uniq	1.1	48.5,26043,99.3	14841,1248173,1263456,1133070

=cut

=item readGenetab 
  
  1. read genetab of gnodes_sam2covtab8e.pl putGenetab8e()
     .. similar topy readCovtab tables, adds 1st col: GeneID and sumgene/readgene rows per GeneID
     .. want maybe ave/median cov depth over chr-positions, also sumgene,readgene
     .. copynum est from CovT/CovM and/or from CovT/C.UCG of sam2genecov genexcopy table
     .. also gene-chr-span vs gene-span of interest.
     .. also skew calc for ave - median diff, ie high-copy sub-spans
     
    "# cols for GeneID, ChrID, Pos: CovT CovM zero = read map counts/bin as per chr.covtab\n";
    "# cols for GeneID, sumgene,  1: allCovT allCovM noCov\n"; 
    "# cols for GeneID, readgene, 1: nrdmaps, nreads, nomap\n";

  2. read genexcopy of gnodes_sam2genecov.pl genestable()
     Gene_ID Glen Nread Rdlen tCopy C.M Uniq Merr C.nz S.Dup,Unq,Nt,Equ,UCG
     esp: C.nz = $cnzave,$cnzmed,$nnz,$cspan
          tCopy = cnzmed/Cnzmed.UCG
  
  3. read genecopyn of annotate ?? this is depreciated, not consistent w/ cds-read map cover calcs
      .. est. gene copies from cds x chr blastn, btall tabulation .. poor job w/ tandem dupls (same chr)
    GeneID____________      tCopy   Nc,s,t  Aln     Span    Glen
    Daplx7pEVm000023t1      1.0     1,1,0   14018   14215   14268
    Daplx7pEVm000024t1      1.0     1,1,0   12367   12336   12317
    Daplx7pEVm000055t1      1.0     1,1,0   12711   12708   12708

=item readGenetab test calcs

genetab8c.info
pt=dplx20maca4pkr_dc_SRR13333792_test8c; 
grep readgene $pt.genetab | perl -ne '@v=split; splice(@v,2,0,1)
if(@v<6); ($id,$pt,$nb,$ct,$cm,$cz)=@v; if($lid ne $id){ putv();
@sv=@v; } else { for $i (3,4,5){ $sv[$i]+=$v[$i]; } } $lid=$id; END{
putv(); } sub putv{ return unless(@sv); my($ct,$cm,$cno)=@sv[3,4,5];
$pno=($cm<1 or $cno == 0)?0: sprintf "%.2f%%",100 * $cno/$cm; 
$pc=($cm<1)?0:sprintf"%.1f",$ct/$cm; print join("\t",qw(GeneID Stat
One aCovT aCovM nMiss Cm pMiss))."\n" if(1>$hdo++); print
join("\t",@sv,$pc,$pno)."\n"; } ' \
 > $pt.rdgenetab

aweed20gnodes/genetab8e.info
pt=arath18tair_chr_SRR10178325_test8e

perl -ne 'next if(/^\W/); @v=split; splice(@v,2,0,1) if(@v<6);
($id,$pt,$nb,$ct,$cm,$cz)=@v; $pok=($pt =~ /sumgene|readgene/)?0:1;
next unless($pok); if($lid ne $id){ putv(); @sv=@v; @mv=([@v[3,4,5]]);
 $nv=1; } else { $lid=$id; next if($v[3]<1); $nv++; push @mv,
[@v[3,4,5]]; for $i (3,4,5){ $sv[$i]+=$v[$i]; } } $lid=$id; END{
putv(); } sub putv{ return unless(@sv); @sva= map{ sprintf"%.0f",
$_/$nv } @sv[3,4,5]; splice(@sv,3,3,@sva); my($ct,$cm,$cno)=@sva;
$pno=($cm<1 or $cno == 0)?0: sprintf "%.2f%%",100 * $cno/$cm; 
$pc=($cm<1)?$ct:sprintf"%.1f",$ct/$cm; $hn=int($nv/2); @mva=(); 
@sm=sort{$$b[0] <=> $$a[0]} @mv; $hm=$sm[$hn]; @mv=@sv; 
splice(@mv,3,3,@$hm); ($mt,$mm,$mno)=@$hm; $pcm=
($mm<1)?$mt:sprintf"%.1f",$mt/$mm;  print join("\t",qw(GeneID Stat One
aCovT aCovM nMiss Cm pMiss))."\n" if(1>$hdo++);
splice(@sv,1,2,"ave",$nv); print join("\t",@sv,$pc,$pno)."\n";
$pnom||=0; splice(@mv,1,2,"med",$nv); print
join("\t",@mv,$pcm,$pnom)."\n"; } ' \
  $pt.genetab > $pt.medgenetab

  
=cut

=item UPD 8i genetab

  gnodes_sam2covtab8i.pl replaces ChrID, Pos w/ GPos column
arath18tair_chr_SRR10178325ab_test8ibwa.genetab
#GeneID	GPos	aCovT	aCovM	aCovU	aCovZ
AT1G01010t1	readgene	373	373	373	0
AT1G01010t1	sumgene	565	565	565	0
AT1G01010t1	100	79	79	79	0
AT1G01010t1	200	82	82	82	0
AT1G01010t1	300	50	50	50	0
AT1G01010t1	400	57	57	57	0
AT1G01010t1	500	102	102	102	0
AT1G01010t1	600	84	84	84	0
AT1G01010t1	700	63	63	63	0

=cut

sub readGenetab {
  my($ingenetab,$outname)=@_;
  # ingenetab == output of sam2covtab: putGenetab8e() == GeneID x Chr x CPos, aCovT/aCovM/aMiss

  my($keyLC1,$keyLL1,$keyLV1)= (LLOCI_CALC == 2) ? ('lCopy2','lLoci2','lCov2') : ('lCopy','lLoci','lCov');
  my($keyLC2,$keyLL2,$keyLV2)= (LLOCI_CALC == 2) ? ('lCopy','lLoci','lCov') : ('lCopy2','lLoci2','lCov2');
  
  my @IC=(3,4,5); # cov value cols
  # my $gloc8i=0; # #GeneID	GPos	aCovT	aCovM	aCovU	aCovZ from putGenetab8i
  if($Genetab8i){ @IC=(2,3,4,5); }
  our($nv,$nzero,$hdo,$outh,@sv,@mv,@rdv,@suv,@hdr)=(0,0);
  our(%cloc,%lcov,%lcovgb,%lcdsgb,%lspan,@lcovall); # Genetab8i cr: loci cov

  my $outf= $outname . '.sumgenetab'; #???
  my($lid,$nout)=(0,0);

  #UPD21AUG20: add CBIN sample here, only in readCovtab() now, unused.
  # ** Also need to add aNbin or aBinN => aLen est corresponds to gLen, aLen/gLen =~ asm copy num found
  # ** need to adjust aCopy w/ aLen/gLen, now aCopy = aCovT/aCovM, 
  # ..  in clean/simple case aCovT/aCovM == aLen/gLen  
  
  #UPD21SEP05: Genetab8i changes: TRY_GENELOC == 2 adds geneloc x crloc rows, skip here or handle separately
  #    col3 key == cr:ID:LOC
  #?? dont output table here, collect in genecovh{id}{fields} ?
  # .. then print gene rows of select fields, stats by gCopy levels, etc
  use constant CR_MINCOV2 => 1; # cr_mincov vs lcdsgb cov filter spurious dup loci
  
  sub putv { if(STORECOV) { storev(); } else { printv(); } }
  sub printv { 
    our($nv,$nzero,$hdo,$outh,@sv,@mv,@rdv,@suv,@hdr);
    my $no=0;
    return $no unless($nv>0 and @sv); 
    my $gid= $sv[0];
    my @sva= map{ sprintf"%.0f", $_/$nv } @sv[@IC]; 
    @sv[@IC]= @sva; # splice(@sv,3,$#sva,@sva); 
    my($ct,$cm,$cu,$cno)=(0) x 4;
    my($mt,$mm,$mu,$mno)=(0) x 4;
    if($Genetab8i) { ($ct,$cm,$cu,$cno)=@sva; } else { ($ct,$cm,$cno)=@sva; }
    
    my $pca=($cm<1)?$ct:sprintf"%.1f",$ct/$cm; 
    my $pno=($cm<1 or $cno == 0)?0: sprintf "%.2f%%",100 * $cno/$cm; 
    if(1>$hdo++) {
      my @hdo= (@hdr >=6) ? @hdr : qw(GeneID Stat One aCovT aCovM nMiss); # nMiss => aMiss ?
      @hdo[1]= 'Stat'; push @hdo,qw( aCopy pMiss); # Cm == Cn == tCopy ?
      print $outh join("\t",@hdo)."\n";
      }
    splice(@sv,1,2,"ave",$nv); 
    print $outh join("\t",@sv,$pca,$pno)."\n"; $no++;

    if(@mv > 0) { #? median always, or option?
      my @sm=sort{$$b[0] <=> $$a[0]} @mv; 
      my $hn=int($nv/2); my $hm=$sm[$hn]; 
      my @mvo=@sv; splice(@mvo,3,3,@$hm); my($mt,$mm,$mno)=@$hm; 
      my $pcm= ($mm<1)?$mt:sprintf"%.1f",$mt/$mm;  
      my $pnom=$pno; #not for median: ($mm<1 or $mno == 0)?0: sprintf "%.2f%%",100 * $mno/$mm; 
      splice(@mvo,1,2,"med",$nv); 
      print $outh join("\t",@mvo,$pcm,$pnom)."\n";   
    }
    
    if(@rdv > 0) { # readgene counts
      if(@IC>3){  ($mt,$mm,$mu,$mno)= @rdv[@IC]; } else { ($mt,$mm,$mno)= @rdv[@IC]; }
      my $pcm= ($mm<1)?$mt:sprintf"%.1f",$mt/$mm;  
      splice(@rdv,1,1,"reads"); #,$nv
      print $outh join("\t",@rdv,$pcm,0)."\n";   
    }
    
    if(@suv > 0 and not @rdv) { # sumgene counts, skip if @rdv
      if(@IC>3){  ($mt,$mm,$mu,$mno)= @suv[@IC]; } else { ($mt,$mm,$mno)= @suv[@IC]; }
      my $pcm= ($mm<1)?$mt:sprintf"%.1f",$mt/$mm;  
      splice(@suv,1,1,"sumc"); #,$nv
      print $outh join("\t",@suv,$pcm,0)."\n";   
    }
    return $no;
  }

  sub lloci_calc2 { 
    my($gid)=@_;
    our($nv, %cloc,%lcov,%lcovgb,%lcdsgb,%lspan,@lcovall); # Genetab8i cr: loci cov
    my ($lmax,$lsum)=(0) x 9;
    my @gb= sort{$a <=> $b} keys %lcovgb;
    my $ngb= @gb; 
    my(@loc2list,%gcloc,%glcov,@covm);
    my($scovm,$ncovm)=(0,0);
    for my $gb (@gb) {
    
      #FIXME: lcovgb{gb} has 2 cov vals for cr: $crval = cds-chr cov and covm = chr-total cov,
      # want both, crval to sort best locus crcb, covm as truest crcb cov val
      # add lcdsgb{gb}{crcb} = crval == genecdscov
      my @crcb= sort{ $lcdsgb{$gb}{$b} <=> $lcdsgb{$gb}{$a} or $a cmp $b } keys %{$lcovgb{$gb}};

      # use cdscov < 0.25*lcovcds test for spurious crcb
      # dont need here a test adjoining crcb, as dup loci measured at same gb bins
      
      my(%glcovgb); my $lcovcds=0;
      for my $crcb (@crcb) { # join near bins as per calc1, but for same gbin
        my($cc,$cr,$cb)=split":",$crcb;
        my $covm=$lcovgb{$gb}{$crcb};
        
        # see below cr_mincov: try this way, not that
        if(CR_MINCOV2) {
        my $covcds= $lcdsgb{$gb}{$crcb} ||0;
        if($covcds > $lcovcds){ $lcovcds= $covcds; } 
        elsif($covcds < 0.25* $lcovcds) { next; } # skip tiny gene-cov crcb
        }
        
        $scovm += $covm; $ncovm++; push @covm,$covm;    
        use constant G2BIN => 1000; # fixme
        my $cbb=int($cb/G2BIN); 
        my $ckey="$cr:".($cbb+0); my $cloc=0; 
        for my $i (0,1,-1,2,-2) { my $ckt="$cr:".($cbb+$i); last if($cloc=$cloc{$ckt}); } 
        $cloc="$cr:".($cbb*G2BIN) unless($cloc); 
        $gcloc{$ckey}=$cloc; $glcovgb{$cloc}+=$covm; # insert $gb here?
      }
      my @loc2listgb= sort{  $glcovgb{$b} <=>  $glcovgb{$a} } keys %glcovgb;
      my $nl= @loc2listgb; $lmax= $nl if($nl>$lmax); $lsum+=$nl;
      map{ $glcov{$_} += $glcovgb{$_} } @loc2listgb; #??
      # ?? guesss loci from @crcb, joining across @gb
    }
    
    # lcopy2 DID over-estimates copy (> all others) for tab8i, due likely to overlapped gbins, tab8j may fix
    # lcopy2 this way is maybe best match to gcopy, better than lcopy, maybe better than acopy
    my $lcopy2= ($ngb<1)? 0 : $lsum / $ngb; #? maybe, or $lmax? or calc1 way? maybe test all for agreement?
    $genecovh->{$gid}->{$keyLC2} = sprintf "%.1f",$lcopy2; # drop ,lmax sprintf "%.1f,%d",$lcopy2,$lmax; # drop ,lmax
 
    @loc2list= sort{  $glcov{$b} <=>  $glcov{$a} } keys %glcov;
    #? use this lLoci1 format?  push @loclist, "$ac,$mdcov,$pspan,$l"; 
    if(@loc2list > $lcopy2){ @loc2list=splice(@loc2list,0,int(0.5+$lcopy2)); }
    $genecovh->{$gid}->{$keyLL2} = join("; ",@loc2list);
    
    @covm= sort{$b<=>$a} @covm; my $covmed= $covm[ int($ncovm/2) ];
    my $covav= ($ncovm<1)?0: $scovm/$ncovm; 
    my $pspant=0; # fixme
       $pspant = ($nv<1 or $lcopy2 < 0.1 )? 0 : 100 * ($lsum / $lcopy2) / $nv; # nv = total gene bins, lsum/lcopy2 = ave gene bins over n-loci ~= ngb
    $genecovh->{$gid}->{$keyLV2} = sprintf"%.1f,%d,%.0f%", $covav, $covmed, $pspant;

  }

  sub lloci_calc1 { # if($Genetab8i)  #UPD21SEP08 ** Output cols add lCopy, lLoc1,2,3,..9
    my($gid)=@_;
    our(%cloc,%lcov,%lcovgb,%lcdsgb,%lspan,@lcovall); # Genetab8i cr: loci cov
    my @cloc=sort{ $lcov{$b} <=> $lcov{$a} } keys %lcov; 
    my $nloc= @cloc;
    @lcovall= sort{$b <=> $a} @lcovall; my $nall= @lcovall;
    my $mdcov= $lcovall[int($nall/2)];
    $genecovh->{$gid}->{$keyLC1} = $nloc;
    #^^ change to lCopy = sum(pspan/100), ie 100p + 100p + 33p = 2.33 lCopy
    ## change {'lLoc'.(1+$i)} to one list: loc1,v,v; loc2,v,v; loc3,v,v; .. ?
    my $lcopy=0; my($ssp,$ssc)=(0,0); my @loclist=();
    for(my $i=0; $i<$nloc; $i++) {
      my $l= $cloc[$i]; 
      my $sc=$lcov{$l}; my $sp=$lspan{$l}||0; 
      my $ac=($sp<1)?0:sprintf "%.0f",$sc/$sp; 
      $ssc+= $sc; $ssp+= $sp;
      my $pspan= int(0.5 + 100*$sp/$nv); $pspan=100 if($pspan>100); #?? is this right?
      # $genecovh->{$gid}->{'lLoc'.(1+$i)} =  "$ac,$mdcov,$pspan,$l"; # ac ~~ aCovm, want median also?
      $lcopy += $pspan/100; 
      if($i < 10) { push @loclist, "$ac,$mdcov,$pspan,$l"; }
      }
    $genecovh->{$gid}->{$keyLC1} = sprintf"%.1f",$lcopy;
    $genecovh->{$gid}->{$keyLL1} = join("; ",@loclist);
    #? add col for avecov,mdcov ?over all loci
    my $pspant= ($nloc<1 or $nv<1)?0: 100*$ssp/($nloc*$nv);
    my $lcova= ($ssp<1)?0:sprintf"%.1f,%d,%.0f%", $ssc/$ssp, $mdcov, $pspant;
    $genecovh->{$gid}->{$keyLV1} = $lcova;
  }
  
  sub storev { 
    our($nv,$nzero,$hdo,$outh,@sv,@mv,@rdv,@suv,@hdr);
    our(%cloc,%lcov,%lcovgb,%lcdsgb,%lspan,@lcovall); # Genetab8i cr: loci cov
    my $no=0;
    return $no unless($nv>0 and @sv); 
    my $gid= $sv[0];
    my @sva= map{ sprintf"%.0f", $_/$nv } @sv[@IC]; 
    @sv[@IC]= @sva; # splice(@sv,3,$#sva,@sva); 
    my($ct,$cm,$cu,$cno)=(0) x 4;
    my($mt,$mm,$mu,$mno)=(0) x 4;
    if($Genetab8i) { ($ct,$cm,$cu,$cno)=@sva; } else { ($ct,$cm,$cno)=@sva; }

    my $pca=($cm<1)?$ct:sprintf"%.1f",$ct/$cm; 
    my $pno=($cm<1 or $cno == 0)?0: sprintf "%.2f%%",100 * $cno/$cm; 
      # my @hdo= (@hdr >=6) ? @hdr : qw(GeneID Stat One aCovT aCovM nMiss); # nMiss => aMiss ?
      # @hdo[1]= 'Stat'; push @hdo,qw( aCopy pMiss); # Cm == Cn == tCopy ?
     
    splice(@sv,1,2,"ave",$nv); 
    $no++;
    if(0 and $debug and $nout < 10) {
      warn "#readgenetab($gid): asv=@sv \n"; ## BAD data: gid == ''; fixed
      # oknow
      #readgenetab(AT1G01010t1): asv=AT1G01010t1 ave 13 85 85 0 
      #readgenetab(AT1G01020t1): asv=AT1G01020t1 ave 6 60 60 0 
      #readgenetab(AT1G01030t1): asv=AT1G01030t1 ave 11 98 98 0 
    }
    
    if(GENEHASH) {
      $genecovh->{$gid}->{'aCovT'} = $ct;
      $genecovh->{$gid}->{'aCovM'} = $cm;
      $genecovh->{$gid}->{'aCovZ'} = $cno;
      $genecovh->{$gid}->{'aCopy'} = $pca;
      $genecovh->{$gid}->{'aBins'} = $nv; #UPD21AUG20, output: aLen= aBins * $CBIN
      # ^^ nv bin count is over-count for aLen, many bins contain partial (ie few %) align to 100b, eg intron, end edges
      # .. any way to correct? using ct/cm deviation from median?
      # eg: Dapsim1EVm000991t1 uniq 1copy, gLen=3510 aLen=5600, medn count=34, 7 bins have 1..9 count, 8 have 10..19
      #  41 have ~full read count, ie 4100 bp  .. some extras =~ alt exons
      
      # $gcovtab{$gid}= { gCopy => $tcopy, gLen => $glen, gNread => $nread, gClass => $isuniq, gCM => $cmed, gCnz => $cnz, };
    }

    if($Genetab8i) { #UPD21SEP08 ** Output cols add lCopy, lLoc1,2,3,..9
      lloci_calc1($gid);      
      lloci_calc2($gid); # LLOCI_CALC => 2 test, sets lLoci2, FIXED: set to lCopy[1]; overest copy bias, due to overlaps?
     }
     
    if(@mv > 0) { #? median always, or option?
      my @sm=sort{$$b[0] <=> $$a[0]} @mv; 
      my $hn=int($nv/2); my $hm=$sm[$hn]; 
      my @mvo=@sv; 
      @mvo[@IC]=@$hm;  # splice(@mvo,3,3,@$hm);       
      if(@IC>3){  ($mt,$mm,$mu,$mno)= @$hm[@IC]; } else { ($mt,$mm,$mno)= @$hm[@IC]; }
      
      my $pcm= ($mm<1)?$mt:sprintf"%.1f",$mt/$mm;  
      my $pnom=$pno; #not for median: ($mm<1 or $mno == 0)?0: sprintf "%.2f%%",100 * $mno/$mm; 
      splice(@mvo,1,2,"med",$nv); 
      if(GENEHASH) {
        $genecovh->{$gid}->{'mCovT'} = $mt;
        $genecovh->{$gid}->{'mCovM'} = $mm;
        # $genecovh->{$gid}->{'mCovZ'} = $cno;
        $genecovh->{$gid}->{'mCopy'} = $pcm;
      }
    }
    
    if(@rdv > 0) { # readgene counts
      if(@IC>3){  ($mt,$mm,$mu,$mno)= @rdv[@IC]; } else { ($mt,$mm,$mno)= @rdv[@IC]; }
      my $pcm= ($mm<1)?$mt:sprintf"%.1f",$mt/$mm;  
      splice(@rdv,1,1,"reads"); #,$nv

      if(GENEHASH) {
      $genecovh->{$gid}->{'aReadT'} = $mt;
      $genecovh->{$gid}->{'aReadM'} = $mm;
      $genecovh->{$gid}->{'aReadU'} = $mu;
      $genecovh->{$gid}->{'aReadZ'} = $mno; #? is it valid
      $genecovh->{$gid}->{'aCopyRd'} = $pcm;
      }      
    }
    
    if(@suv > 0 and not @rdv) { # sumgene counts, skip if @rdv
      if(@IC>3){  ($mt,$mm,$mu,$mno)= @suv[@IC]; } else { ($mt,$mm,$mno)= @suv[@IC]; }
      my $pcm= ($mm<1)?$mt:sprintf"%.1f",$mt/$mm;  
      splice(@suv,1,1,"sumc"); #,$nv
      if(GENEHASH) {
      $genecovh->{$gid}->{'aReadT'} = $mt;
      $genecovh->{$gid}->{'aReadM'} = $mm;
      $genecovh->{$gid}->{'aReadU'} = $mu;
      $genecovh->{$gid}->{'aReadZ'} = $mno; #? is it valid
      $genecovh->{$gid}->{'aCopyRd'} = $pcm;
      }      
    }
    return $no;
  }
  
  
  my($ok,$inh)= openRead($ingenetab);
  if(STORECOV) { 
    $outf="";  #  no output STORECOV, genecovh then output table
  } elsif($ok) {
    # rename($outf,"$outf.old") if(-f $outf); open($outh,'>',$outf);
    ($ok,$outh)= openOut($outf);    
  }
  return unless($ok);

  my($lcr,$lcb)=(0,0,0);
  while(<$inh>){
  
    if(/^\W/) {
      if(/^#GeneID/i){ 
        # GeneID ChrID Pos aCovT aCovM noCov : putGenetab8f; expect, any others?
        # GeneID, GPos: aCovT aCovM aCovU aCovZero: putGenetab8i
        s/^#//; @hdr=split;  
        if($hdr[1] =~ m/GPos/i) { $Genetab8i=1; @IC=(2,3,4,5); } # putGenetab8i
        elsif($hdr[1] =~ m/ChrID/i) { $Genetab8i=0; @IC=(3,4,5); }
        }
      next;
    } 

    #?? splice(@v,2,0,1) if(@v<6); #< no, check header
    # cols for GeneID, ChrID, Pos: CovT CovM zero = read map counts/bin as per chr.covtab
    # cols for GeneID, sumgene,  1: allCovT allCovM noCov 
    # cols for GeneID, readgene, 1: nrdmaps, nreads, nomap
    
    my @v=split; 
    my($id,$pt,$nb,$ct,$cm,$cz)=@v; 
    if($lid ne $id){
      $nout+= putv() if($nv); 
      $nv=0; @suv=@rdv=@mv= @sv=();  ($lcr,$lcb)=(0,0);
      @sv=@v; for my $i (@IC){ $sv[$i]=0; } 
      %cloc=%lcov=%lcovgb=%lcdsgb=%lspan=(); @lcovall=();
    }
    
    if($pt =~ /sum|read|gene/) { 
      # collect sum,read nums for CN, aves
      if($pt =~ /readgene/) { @rdv= @v; }
      elsif($pt =~ /sumgene/) { @suv= @v; }

    } elsif($Genetab8i and ($nb =~ /^cr:/ or $cz =~ /^cr:/)) { # or $cz =~ /^cr:/ ?
      my($gid,$gb,$crcb,$covt,$covm,$covu,$crval)=(0) x 9;
      
      # note: crval (cov for cds == chr reads) != covm (all reads from chr.tab)
      # use crval to filter out low-cov spurious cr: locs, commonly 0,1,2 vals, esp drop crval==0 cases
      if(@v < 6) { ($gid,$gb,$crcb,$covm)=@v; $crval= $covm; } 
      elsif($cz =~ /^cr:/) { ($gid,$gb,$covt,$covm,$covu,$crcb,$crval)=@v; }   
      else { ($gid,$gb,$crcb,$covt,$covm,$covu,$crval)= @v; } #? prefer this way
 
      my $cr_mincov=0;
      # cr_mincov filter here may be problem for uneven/part duplications, drop? add into calc2() per gbin?
      unless(CR_MINCOV2) {      
        $cr_mincov= ($nv<1) ? 1 :  0.15 * $sv[ $IC[1] ]/$nv; # ugh, change for putGenetab8j ?
        $cr_mincov=1 if($cr_mincov < 1); #? or MINRD
      }       
      
      if($crval > $cr_mincov) { # was covm; need to filter out tiny cov cases, what mincov ? use 0.20?*$sv[1]/$nv as min-expected
      my($cc,$cr,$cb)=split":",$crcb;
      
      $lcovgb{$gb}{$crcb}=$covm; # for LLOCI_CALC==2, save till putv($gid) ?
      $lcdsgb{$gb}{$crcb}=$crval||1; # $covm; # crval is best gene-locus cov, while covm is best chr cov, use both
 
      use constant GBIN => 1000; # was 2000; # fixme, match putGenetab8i .. 
      # but 8j drops GENEBIN=2000 use, no good for huge gnomes, crcb bins are BN=100 
      # .. this is for calc1, calc2 doesnt use
      my $cbb=int($cb/GBIN); 
      my $ckey="$cr:".($cbb+0); my $cloc=0; 
      # need to merge some crcb overlaps, here? 
      for my $i (0,1,-1,2,-2) { my $ckt="$cr:".($cbb+$i); last if($cloc=$cloc{$ckt}); } 
      $cloc="$cr:".($cbb*GBIN) unless($cloc); 
 
      $cloc{$ckey}=$cloc; $lcov{$cloc}+=$covm; $lspan{$cloc}++; # insert $gb here?
      push @lcovall,$covm; 
      }  
      
#     } elsif($Genetab8i and $cz =~ /^cr:/) { # merge w/ above, only col order differs
#       my($gid,$gb,$covt,$covm,$covu,$crcb,$crval)=@v;  
#       my $cr_mincov= ($nv<1) ? 1 :  0.20 * $sv[ $IC[1] ]/$nv; # ugh, change for putGenetab8j ?
# 
#       if($crval >= $cr_mincov) {
#       my($cc,$cr,$cb)=split":",$crcb;
#       
#       $lcovgb{$gb}{$crcb}= $covm; # or "$covm\t$covt"; #?? want both or just covm ?
# 
#       my $cbb=int($cb/GBIN); 
#       my $ckey="$cr:".($cbb+0); my $cloc=0; 
#       for my $i (0,1,-1,2,-2) { my $ckt="$cr:".($cbb+$i); last if($cloc=$cloc{$ckt}); } 
#       $cloc="$cr:".($cbb*GBIN) unless($cloc); 
#       $cloc{$ckey}=$cloc; $lcov{$cloc}+=$covm; $lspan{$cloc}++; # insert $gb here?
#       push @lcovall,$covm; 
#       }
      
=item new,new
      # FIXME: this may be new,new format, mixed w/ next form, which is best?
      # ($gid,$gb,$crcb,$covt,$covm,$covu,$crval)=@v;  
   AT1G01010t1	100	cr:NC_003070.9:3900	8 : orig cr: form
   AT1G01010t1	100	cr:NC_003070.9:3900	44	44	44	8  : new2, chrcov merge
   AT1G01010t1	100	44	44	44	cr:NC_003070.9:3900	8  : new1, chrcov merge

#GeneID	GPos	aCovT	aCovM	aCovU	aCovZ
AT1G01010t1	readgene	372	372	372	0
AT1G01010t1	100	cr:NC_003070.9:3900	44	44	44	8
AT1G01010t1	100	43	43	43	0
AT1G01010t1	100	cr:NC_003070.9:3800	32	32	32	32
AT1G01010t1	200	19	19	19	0
AT1G01010t1	200	cr:NC_003070.9:4100	24	24	24	3
AT1G01010t1	200	cr:NC_003070.9:4000	40	40	40	21

=item test new tab8j cr: cols

AT1G01810t1  readgene  125     125     125     0
AT1G01810t1     100     36      36      36      0
AT1G01810t1     200     39      39      39      0
AT1G01810t1     300     37      37      37      0
AT1G01810t1     400     38      38      38      0
                    v-- chr.covtab covt,m,u --v at cr:cb loc
AT1G01810t1     100     33      33      33      cr:NC_003070.9:295300   27
AT1G01810t1     200     36      36      36      cr:NC_003070.9:295400   37
AT1G01810t1     300     37      37      37      cr:NC_003070.9:295500   35
AT1G01810t1     400     37      37      37      cr:NC_003070.9:295600   29

#FIXME 21SEP14: this GBIN test is no good for large gnomes like pig w/ longer introns
# .. ?? change to counting n cr: per gene-bin as better locus count? join across gene-bins the nearest/same-cr set?
# gnodes_sam2covtab8j.pl putGenetab8j() changes gene-bins to fixed set, removes problem of overlapping bins?
# lLoci est. then is ave(n-crbins) over gene-bins? collect all gb/crcb per gid to calc

=cut


    } else { # pt == ChrID
      my($cr,$cb)= ($Genetab8i) ? (0,$pt) : ($pt,$nb); # proper var names

      if($ct < $MINRD)  #? should be ct < 1 ?
      { 
        $nzero++; # ** ?? need to record zeros here ?? this blocks aCovZ ?
      } else { 
        $nv++; 
        push @mv, [@v[@IC]]; 
        for my $i (@IC){ $sv[$i]+=$v[$i]; } 
        #here? $cr_mincov= ($nv<1) ? 2 :  0.15 * $sv[ $IC[1] ]/$nv; # Genetab8j
      }   
      
      # copy from readCovtab; Genetab8i: cr/lcr == 0, cb == gene-bin
      if($lcr eq $cr and $lcb>0 and $cb > $lcb) {  
        if(@CBIN < kSAMPLE) { 
          my $bspan= $cb - $lcb; push @CBIN, $bspan; 
          my $nbin=@CBIN; $CBIN= $CBIN[int($nbin/2)] if($nbin>100);
          } 
      }
      ($lcr,$lcb)= ($cr,$cb);
    }
    
    $lid=$id; 
  } close($inh);
  
  $nout+= putv(); unless(STORECOV) { close($outh); }
  warn "#readGenetab: n=$nout, cbin=$CBIN, $outf\n" if $debug;

  return($nout, $outf);  
}


sub readCovtab { 
  my($covtab, $iscdscov)=@_; #? for both cds,chr covtabs?
  open(my $inh,$covtab) or die "reading $covtab";
  # sumgenecov: only iscdscov here ? store to genecovh->{gid}->{field} w/  constant  STORECOV => 1, as readGenetab
  
  my($nin,$nok,$lcr,$lcb,$lida)= (0) x 9;
  # my(@cbin); => @CBIN
  
  while(<$inh>) {
    my @v=split; 
    my $hasivar= (@v >= 8); my $hasann= (@v >= 9);
    if(/^\W/){ 
      if($hasivar and not m/[,=]/ and /Cov/i) { @lvar=@v[@IVAR] unless(@lvar); }  
      next;
    }
    next unless($hasivar); # error? 
    my($cr,$cb,@av)=@v; $nin++;
    
    my($cCDS,$cTE,$cUNK)= splice(@av,0,3); #($ct,$cm,$cu) NOW == CDSrd, TErd, UNKrd
    my($act,$acm,$acu)  = splice(@av,0,3); # == @IVAR, same  rd_total, rd/mmap, rd_uniq for chrasm x rd map 
    my($an,$ids)= (@av>0) ? splice(@av,0,2) : ("",""); # may not exist

    if($nann>0 and my $anidx= $annvals{$cr}{$cb}){
      my($anx,$idx)=split" ",$anidx;
      if($anx){ $an = ($an)? "$an,$anx" : $anx; }
      if($idx){ $ids= ($ids)? "$ids,$idx" : $idx; }
    } 
    # #? use both annvals, idclass? was elsif($nidclass)
    if($nidclass) { #  added from idclassf; use chrid or $ids      
      $ids=$cr unless($ids); #trick for cds.cov; see below ids and busco
      for my $id (split(",",$ids)) {
        if($id and my $anx= $idclassh->{$id} ) { unless($an =~ m/$anx/){ $an = ($an) ? "$an,$anx" : $anx; last; } }
      }
    }

    my $ok= ( $act >= $MINRD )?1:0; # skip? if too low
    next unless($ok);
    $nok++;
    
    # ucg_id_cov => gene_cov table : genecov{id}{cdsm,cdst,chrm,chrt}= @cov_vals (act,acm,acu?)
    # assume 20_000 to 50_000 gene ids, ave cds-length of 2_000/100 bins
    #   %genecov becomes hash-array of 50k x 4 x 20 vals =~ 4 mil vals? 10s of megabytes of mem?
    
    if($iscdscov) {
      # cr == geneid, all together
      my $ida= $cr;
      #? if($ida ne $lida) { putgene(xxx); }  
      #? gene_sum( $cr, $act,$acm,$acu);
      #x push @{$ucg_id_cov{$ida}}, $acm; # keep all? $act,$acm,$acu; rename ucg_id_cov > gene_id_cov
      push @{$ucg_id_cov{$ida}{cdsm}}, $acm; # keep all? $act,$acm,$acu; rename ucg_id_cov > gene_id_cov
      push @{$ucg_id_cov{$ida}{cdst}}, $act; # keep all? $act,$acm,$acu; rename ucg_id_cov > gene_id_cov
      $ucg_id_cov{$ida}{ncds}++;
      $lida= $ida;
      
    } else {
      if($ids and $an =~ /CDS/i) { # require CDS annot? 
        my($ida,@idc); 
        #UPD21apr26: ?? revert to count 1 id/bin locus, more skews gene sums : No, not that either
        # .. both ways are wrong: count all ids over-counts some lowqual id-hits, count 1 under-counts some paralogs
        # .. likely need better cds.anntab, separate best/2nd cds id aligns, count all best/top per bin-locus, not 2nds
        my @ids= split",",$ids; #? add ALL CDS ids ..uniq? ?
        if($nidclass) { @idc= grep{ $idclassh->{$_} =~ m/CDS/i } @ids; }
        else { @idc= grep{ $ucg_id_cov{$_} } @ids; } # this way only? assumes have cdscov, 1st read
        #old2: for $ida (@idc)
        #x ($ida) = shift @idc; # this way gives many more over-asm; 
        #x   .. arath18ccc_SRR10178325 casm<CU	1390 vs casm<CU	57 arath18cab_genesum
        #x ($ida) = sort @idc; # no change from unsort
        # if($ida) # UPD21apr26??
        for $ida (@idc)
        {
          push @{$ucg_id_cov{$ida}{chrm}}, $acm;  
          push @{$ucg_id_cov{$ida}{chrt}}, $act;  
          $ucg_id_cov{$ida}{nchr}++;
          $lida= $ida; # gene ids unordered on chrcov
        }
                
      }
    }

    if($lcr eq $cr and $lcb>0 and $cb > $lcb) { #?
      if(@CBIN < kSAMPLE) { 
        my $bspan= $cb - $lcb; push @CBIN, $bspan; 
        my $nbin=@CBIN; $CBIN= $CBIN[int($nbin/2)] if($nbin>100);
        } 
    }
    
    ($lcr,$lcb)= ($cr,$cb);
  } close($inh);

  warn "#readCovtab($covtab,cds=$iscdscov) nok=$nok, nin=$nin\n" if $debug;
  return($nok);
}

=item putGeneCovSum read genecov means table, write summary of classes
  
  sumclass ok, xcopy levels (<=0.5,1,>=2,>=3,..), locov, hicov, tooshort, ..

genemeans:
GeneID	CD_med	CD_ave	CD_nit	CD_sem	CD_rng	Ch_med	Ch_ave	Ch_nit	Ch_sem	Ch_rng	CDt_med	CDt_ave	Cht_med	Cht_ave
dropse20uc:g6903256t1	89	88.42	566	3.75	28,121	92	91.30	621	3.69	62,120	89	88.49	92	91.37
dropse20uc:g6900405t1	98	97.69	458	4.62	52,198	97	97.72	480	4.51	61,183	143	216.65	121	209.09
dropse20uc:g6902868t1	92	91.26	447	4.36	10,134	92	92.72	475	4.29	44,133	92	91.26	92	92.72
...
dropse20uc:g6898477t1	0	0	0	0	0	102	107.00	3	62.02	101,118	0	0	102	107.00
dropse20uc:g117184805t1	0	0	0	0	0	118	115.50	2	0.00	113,118	0	0	118	115.50
allgenes	88	109.06	14327	4.30	22,28894	95	98.03	14326	1.58	36,13623
  ^total now

=cut

sub diffratio {
  my($va,$vb,$mindif)=@_; $mindif||=0.66;
  my($cl,$r)=(0,0);
  return ($cl,$r) if($va<0.01 and $vb<0.01);

  if($vb < 0.01) { $r=($va>=1) ? 9 : 1; }
  else { $r= $va/$vb;  }
  $cl= ($r < 0.57)? -1 : ($r > 1.75) ? 1 : 0;
  # 1/1.75 =~ 0.57; 1/0.55 =~ 1.75
    
  return ($cl,$r); 
}

# sub putXcopy { 
#   my($Oh, $id, @ctn)=@_;
#   # @ctn= ($xcds,$txcds,$ncds, $xchr,$txchr,$ncr), only 2 cols here CDS,Chr
#  
#   sub xfmt{ my($v)=@_; my $d=($v>99)?0:($v>9)?1:2; return sprintf "%.${d}f",$v; }
#   # my $XCUT=$ENV{xcut}|| 1.7; my $LCUT= $ENV{lcut}||0.55;
#   my($hic,$hiv)=("","");   
#   print $Oh $id; 
#   while( my @c= splice(@ctn,0,3) ) {
#     $hic .= ($c[1]>$XHICUT)?"2":($c[1]<$XLOCUT)?"0":"1"; 
#     $hiv .= ($c[0]>$XHICUT)?"x":($c[0]<$XLOCUT)?"z":"n";  
#     print $Oh "\t", join(",",xfmt($c[0]),xfmt($c[1]),$c[2]);  
#     }
#   # add counts of CApatt to outsum, hash it
#   my $capatt="c$hic$hiv";
#   print $Oh "\t",$capatt,"\n"; 
#   return ($capatt);
# }

# sub putGeneCovSum {  # REPLACE/upd
#   my($outsum, $outmeans)= @_;
# 
#   use constant { kIT => 0, kTOOSHORT => 1, kDIFFBN => 4, kDIFFLEN => 2, 
#         kDIFFKU => 8, kDIFFKUn => 64, kDIFFKUC => 16, kDIFFKUCn => 32,
#         kTOTKU =>128, kTOTKUn=>256 };
#   my @klasses= (kTOOSHORT,kDIFFKU,kDIFFKUn,kDIFFKUC,kDIFFKUCn, kTOTKU,kTOTKUn, kDIFFBN,kDIFFLEN);
#   my %klabels= ( 0 => 'item', 1 => 'short',4 => 'dLENbn', 2 => 'dLENcds', 
#           8 => 'casm>CU', 64 => 'casm<CU', 16 => 'cds>CU', 32 => 'cds<CU',
#           128 => 'ctot>CU', 256 => 'ctot<CU' );
#   #? %klabels constant vals are not inserted in hash keys !!
#   #? my %klabels= ( kIT => 'item', kTOOSHORT => 'tooshort', kDIFFBN => 'diffbn', kDIFFKU => 'diffku', kDIFFLEN => 'difflen' );
#   
#   my $MINBN=2; # min bins, ie 2 x 100b
#   my $pDIFFBN=0.66; # cdsbins <> chrbins
#   my $pDIFFKU=0.66; # cov/KU 
# 
#   my $KU= $KU_OPT || $samvals{$asmid}{kucg} || $samvals{$masmid}{kucg}; # should be global -option
#   unless($KU) { $KU=  _minnot0($samvals{$asmid}{kulo}, $samvals{$asmid}{kuhi} );   }      # or _max ?
#   unless($KU){ 
#     open(F,"tail $outmeans |"); while(<F>){ my @v= split; 
#     if($v[0] eq $TOTALID){ $KU= _max($v[1], $v[6]);  } } close(F);
#   }
#   
#   # add outxcovtab => name.genemeans => name.genexcopy * or .genexcovtab ?
#   (my $outxcopy=$outmeans) =~ s/.genemeans//; $outxcopy.=".genexcopy";
#   warn "#putGeneCovSum($outmeans => $outsum, $outxcopy), asmid=$asmid, KU=$KU\n" if $debug;
# 
#   open(F,$outmeans) or return;
#   # rename($outxcopy,"$outxcopy.old") if(-f $outxcopy); 
#   # open(my $xouth,'>',$outxcopy) or die "writing $outxcopy";  
#   my($ook,$xouth)= openOut($outxcopy);
# 
#   # for putXcopy(name.genexcopy)
#   my @XCOLS=("GeneID","CDS_xC,tC,nB","Asm_xC,tC,nB","CApatt");
#   #above# my $XHICUT=$ENV{xcut}|| 1.7; my $XLOCUT= $ENV{lcut}||0.55;
# 
#   print $xouth join("\t",@XCOLS)."\n";
# 
#   my($nit);
#   my(%class,%sums,%capatt,%capats,@hd,@tcds,@tchr);
#   
#   while(<F>) {
#     next if(/^\W/);
#     my @v=split;
#     my $id=$v[0];
#     if($id =~ /^GeneID/){ @hd=@v; next; }
#     elsif($id eq $TOTALID) { 
#      @tcds=@v[1..5]; @tchr= @v[6..10]; 
#      unless($KU) { $KU= _max($v[1], $v[6]); } # median of medians, chr val most accurate
#      next; } # last in?
#     
#     $nit++;
#     my($cdsm,$cdsa,$cdsn,$cdse,$cdsr)= @v[1..5];
#     my($crm,$cra,$crn,$cre,$crr)= @v[6..10]; # may be missing?
#     my($cdsrt,$crrt)= @v[11,13];
#     # for xcopy: id, cdsm/KU,cdsrt/KU,cdsn,  crm/KU, crrt/KU, crn, cpatt: c11nn/c22xn/...
#     # FIXME: cpatt==c00zz must ignore zz = xchr,xcds
#     
#     my $cdsw= ($ntcds) ? $cdslen->{$id} : 0;
#     
#     my $cl=0; my $tooshort=0;
#     if($cdsn < $MINBN and $crn < $MINBN) { $cl |= kTOOSHORT; $tooshort=1; }
#     if($tooshort) { } # skip other classes?
# 
#     my($clwc,$rwc,$clwa,$rwa)= (0) x 9;
#     if($cdsw and not $tooshort) { my $cdswb= $cdsw/$CBIN; 
#       ($clwc,$rwc)= diffratio(  $cdsn, $cdswb, $pDIFFBN); 
#       ($clwa,$rwa)= diffratio(  $crn, $cdswb, $pDIFFBN); 
#       if($clwa or $clwc) { $cl |= kDIFFLEN; }
#       push @{$sums{'wchr'}}, $rwa;
#       # compare cds-length to readmap len, both cds, chr? measures dups and underasm ?
#       }
# 
#     # for KU set, count <1 and >1
#     my($cldif,$pdifbin)= ($tooshort)?(0,1):diffratio($cdsn, $crn, $pDIFFBN);
#     if($cldif) { $cl |= kDIFFBN; } # cov span differs, pdifbin < 1 means crn larger, has dups, pdifbin > 1 means crn missing cds?
# 
#     # tCopy vals: add total cov counts, both cdsrt, crrt? or just one?
#     # do these before xchr,xcds, if txcds & txchr < LOCUT, skip xchr,xcds test
#     my($cltcds,$txcds)=  diffratio($cdsrt, $KU, $pDIFFKU);
#     my($cltchr,$txchr)=  diffratio($crrt, $KU, $pDIFFKU);
#     if($cltchr<0) { $cl |= kTOTKUn; } # xcov != 1, for chrasm , 
#     elsif($cltchr>0) { $cl |= kTOTKU; } # xcov != 1, for chrasm
# 
#     ## xCopy vals
#     my($noxcopy)= ($tooshort or ($cltchr<0 and $cltcds<0) )?1:0;
#     my($clxc,$xcds)= ($noxcopy)?(0,1) : diffratio($cdsm, $KU, $pDIFFKU);
#     my($clxa,$xchr)= ($noxcopy)?(0,1) : diffratio($crm, $KU, $pDIFFKU);
#     if($clxc<0) { $cl |= kDIFFKUCn; } elsif($clxc>0) { $cl |= kDIFFKUC; } # xcov != 1, for cds
#     if($clxa<0) { $cl |= kDIFFKUn; }  elsif($clxa>0) { $cl |= kDIFFKU; } # xcov != 1, for chrasm
# 
#     # cds <> chr xcov, ignore this?
#     my($clxb,$xcds_chr)=  diffratio($xcds, $xchr, $pDIFFKU); 
#     
#     # want counts of xcopy per gene: on cds and on chr, count too many, too few, too short
#     # sum classes here? ave,medn for xchr,xcds, diff(cds,chr)
#     $class{0}++;
#     for my $ic (@klasses) { if($cl & $ic) { $class{$ic}++; }  }
#     push @{$sums{'difbin'}}, $pdifbin; # was difbin
#     push @{$sums{'xcds'}}, $xcds;
#     push @{$sums{'xchr'}}, $xchr;
#     push @{$sums{'ycds_cr'}}, $xcds_chr;
#     push @{$sums{'tchr'}}, $txchr;
#     
#     #add for xcopy: id, xcds=cdsm/KU,txcds=cdsrt/KU,cdsn,  xchr=crm/KU, txchr=crrt/KU, crn, cpatt: c11nn/c22xn/...
#     my($capatt)= putXcopy( $xouth, $id, $xcds, $txcds, $cdsn, $xchr, $txchr, $crn);
#     
#     $capatt{$capatt}++;
#     # for only 4 patt sym, add totals c.2../c.1../c.0.., c...x/c...n/c...z ? 
#     # should == casm>CU,casm<CU, ctot>CU,ctot<CU
#     # FIXME: c00zz must ignore zz
#     my($cap); # capats separate counter from capatt
#     ($cap=$capatt) =~ s/c...(.)/c...$1/; $capats{$cap}++; #casm = z/n/x
#     ($cap=$capatt) =~ s/c.(.)../c.$1../; $capats{$cap}++; #ctot = 0/1/2
#   }
#   
#   # rename($outsum,"$outsum.old") if(-f $outsum); 
#   # open(my $outh,'>',$outsum) or die "writing $outsum";  
#   my($ook,$outh)= openOut($outsum);
# 
#   print $outh "SUMMARY of $outmeans, asmid=$asmid \n";
#   print $outh "Params: CU=$KU, xCopy(Cov/CU)<>1.0: hi=$XHICUT, lo=$XLOCUT\n";
#  
#   #? make 2 cols, add capatt right of Class_counts? max=10 rows
#   print $outh "Class_Counts__________\n";
#   
#   my $it=0;
#   for my $ic (kIT,@klasses) {
#     my $c= $class{$ic}||0; 
#     my $lb= $klabels{$ic}||"na$ic"; # why bad here? hash{kIT} == {'kIT'} not {0}
#     print $outh "$lb\t$c"; print $outh ($it % 2 == 0) ? "\t" : "\n"; $it++; 
#   } print $outh "\n" if($it % 2 == 1); 
#   
#   my %clabs=(
#     'c.1..' => "\tasm-one-copy",
#     'c.2..' => "\tasm-two+copy",
#     'c.0..' => "\tasm-miss-copy",
#     'c...x' => "\tasm-miss-dup",
#     'c...z' => "\tasm-extra-dup",
#     'c...n' => "\tasm-norm-dup",
#   );
#   my @caprows=(); 
#   my @cap= sort{ $capatt{$b} <=> $capatt{$a} or $a cmp $b} keys %capatt;
#   for my $i (0 .. $#cap) { my $c=$cap[$i]; push @caprows, "$c\t".$capatt{$c}; }
#   @cap= sort keys %capats; # put last? first? always only 6?
#   for my $i (0 .. $#cap) { my $c=$cap[$i]; my $cl=$clabs{$c}||""; push @caprows, "$c\t".$capats{$c}.$cl; }
#   
#   my $nct=@caprows;  
#   print $outh "\nxCopypat_Counts cpatt:tCDS,tAsm,xCDS,xAsm\n";  
#   my($ir,$ic)= ($nct>36)?(9,5) :($nct>27)?(9,4) :($nct>18)?(9,3) : ($nct>9)?(9,2) : ($nct,1);
#   for my $j (1 .. $ir) { 
#     my @row=();
#     for my $k (1 .. $ic) {
#       my $kk= $ir * ($k-1) + ($j-1);
#       push @row, $caprows[ $kk ] unless($kk>=$nct); 
#     }
#     print $outh join("\t",@row)."\n";
#   }
#   
#   
#   print $outh "\nMedian,Mean,N,StDev,Range________\n";
#   for my $sk (sort keys %sums){
#     my @v=  @{$sums{$sk}};
#     my @stat= meanvar(@v); # $med,$ave,$nit,$se,$range
#     # my @cors= corcovar(\@lstat,\@lv,\@stat,\@v); # $cor,covar : pairwise? matrix for all cor vars?
#     print $outh join("\t",$sk,@stat)."\n";
#   }
#   
#   close($outh);
# }


# sub putGeneCovStats {
#   my($outmeans)= @_;
#   
#   # rename($outmeans,"$outmeans.old") if(-f $outmeans); 
#   # open(my $outmeanh,'>',$outmeans) or die "writing $outmeans";  
#   my($ook,$outmeanh)= openOut($outmeans);
#   
#   # note: CBIN * n = bases, add Mb vals?
#   my($nput,@cdsmd,@chrmd)=(0);
#   my @ids= sort{ $ucg_id_cov{$b}{ncds} <=> $ucg_id_cov{$a}{ncds} or $a cmp $b } keys %ucg_id_cov; # by cov size
#   my ($didhdr);
#   for my $ida (@ids) {
#     my $ncds= $ucg_id_cov{$ida}{ncds};
#     my $nchr= $ucg_id_cov{$ida}{nchr};
#     my @cdsm= ($ncds) ? @{$ucg_id_cov{$ida}{cdsm}} : ();
#     my @cdsstat= meanvar(@cdsm); # $med,$ave,$nit,$se,$range
#     my @cdsts  = meanvar(($ncds) ? @{$ucg_id_cov{$ida}{cdst}} : ());  
#     
#     my @chrm= ($nchr) ? @{$ucg_id_cov{$ida}{chrm}} : ();
#     my @chrstat= meanvar(@chrm); # $med,$ave,$nit,$se,$range
#     my @chrts= meanvar(($nchr) ? @{$ucg_id_cov{$ida}{chrt}} : ());  
# 
#     push @cdsmd, $cdsstat[0] if($ncds);
#     push @chrmd, $chrstat[0] if($nchr);
#     
#     unless($didhdr++) {
#     my @stath= qw(med ave nit std rng); # std was sem
#     my @gstat= map{ "CD_".$_ } @stath;
#     my @cstat= map{ "Ch_".$_ } @stath;
#     my @gstt = map{ "CDt_".$_ } @stath[0,1];
#     my @cstt = map{ "Cht_".$_ } @stath[0,1];
#     print $outmeanh join("\t","GeneID",@gstat,@cstat,@gstt,@cstt)."\n" ;
#     }
#     print $outmeanh join("\t",$ida,@cdsstat,@chrstat,@cdsts[0,1],@chrts[0,1])."\n"; $nput++;
#     
#     # calc xCopy? need chr,cds sum val?
#     # how many xcopy in @cds, how many in @chr, diff?
#   }
#   
#   my @cdstt= meanvar(@cdsmd);
#   my @chrtt= meanvar(@chrmd);
#   print $outmeanh join("\t",$TOTALID,@cdstt,@chrtt)."\n";
#   
#   close($outmeanh);
#   return($nput);
# }


=item eg Arath copynum tables
  
  -- note arath20max is ~10 Mb larger than arath18ta, and likely has some more gene copies
  .. but those stats may be obscured by uncertain filtering of what is a copy
  -- maybe bugs counting tandupl copies: tCopy from Align/Glen, but tandup addition missing?
     .. and tandup copy > 1 missing now from btall table (only 1 tandup tabulated) affects large chrs most      
       
==> arath18tair_chr-arath18tair1cds.genecopyn <==
GeneID_____     tCopy   Nc,s,t  Aln     Span    Glen   xCopy   gCopy
AT2G05510t1     0       0,0,0   0       0       315     0       0.9
AT2G07623t1     0       0,0,0   0       0       246     0       6.4
AT5G03710t1     0       0,0,0   0       0       246     0       21.6
#---------------------------------------
# Genes Copynum Summary by Copy level, for nGene=27359 found, 3 missed
# tCopy=Chrasm, gCopy=Gene set, xCopy=tCopy/gCopy, x0,xLo,xEq,xHi= xCopy 0,<1,=1,>1
#CnGene nGene   tCopy   xCopy   gCopy   x0,xLo,xEq,xHi
#0      713     1.0     2.3     0.5     0,0,93,620
#1      26212   1.0     1.0     1.0     1,227,25828,156
#2-9    401     1.4     0.5     3.2     1,350,47,3
#10-99  23      1.4     0.1     22.7    1,22,0,0
#99-499 6       1.2     0.0     160.9   0,6,0,0
#500+   7       2.6     0.0     831.5   0,7,0,0
#---------------------------------------

==> arath20max_chr-arath18tair1cds.genecopyn <==
AT5G59616t1     0       0,0,0   0       0       462     0       1.1
AT5G59650t1     0       0,0,0   0       0       2655    0       0.9
AT5G61710t1     0       0,0,0   0       0       468     0       1.8
#---------------------------------------
# Genes Copynum Summary by Copy level, for nGene=27175 found, 176 missed
# tCopy=Chrasm, gCopy=Gene set, xCopy=tCopy/gCopy, x0,xLo,xEq,xHi= xCopy 0,<1,=1,>1
#CnGene nGene   tCopy   xCopy   gCopy   x0,xLo,xEq,xHi
#0      702     1.0     2.3     0.5     0,0,101,601
#1      26212   1.0     1.0     1.0     164,317,25556,175
#2-9    401     1.3     0.5     3.2     11,329,54,7
#10-99  23      1.5     0.1     22.7    1,22,0,0
#99-499 6       1.8     0.0     160.9   0,6,0,0
#500+   7       4.5     0.0     831.5   0,7,0,0
#---------------------------------------
    
=cut
  
=item putGeneCN UPD21Jun04 copynum table
 
 # parts from gnodes_annotate
 
=cut

sub putGeneCN_OLD { 
  my($intabh, $copytab)= @_;
  
  use constant COPYNUM => 1;  # UPD21Jun04
  use constant COPYNUMzeros => 1;  # count chr-align missed genes in ave t/xCopy?
  use constant { C0hi=>0.66, C1hi=>1.55, C2hi=>9.99, C3hi=>99.9, C4hi=>499 };
  my @cn=(); our($ctabh,%cnsum,$ncnsum,%cndid,$ncmiss);
  my $ngenecov=1; my $lasttd;
  # my $copytab=$oname . ".genecopyn"; # used .genexcopy, keep that?
  
  if(COPYNUM) {
    open( $ctabh,">$copytab"); # save old?
    my @chd=qw(GeneID tCopy Nc,s,t Align Span Glen); my $xlab="";
    if($ngenecov>0) { push @chd, "xCopy","gCopy"; $xlab=" xCopy=tC/gC; gCopy=gene x DNA copy num"; }
    
    print $ctabh "#cols: tCopy = gene-align-spans / gene-length; Nc,s,t = Num copies (scaf,tandem); \n";
    print $ctabh "#cols:  Align, Span= sum align,spans; Glen = gene-length $xlab\n";
    print $ctabh "#".join("\t", @chd)."\n"; 
  }
  
  sub putcn { our($ctabh,%cnsum,$ncnsum,%cndid,$ncmiss);
    my($id,$tlen,$tcopy,$saln,$sspan,$ncp,$nscaf,$ntand)=@_;
    return 0 unless($saln > 0 and $ncp > 0);
    my @ecols=($saln,$sspan,$tlen);
    # if($ngenecov>0 and exists $genecovh->{$id}) 
    if($ngenecov>0) { 
      my $gcopy=0;
      if(GENEVEC) {
      $gcopy= (exists $genecovh->{$id}) ? $genecovh->{$id}->[0] : 0; 
      } else {
      $gcopy= (exists $genecovh->{$id}) ? $genecovh->{$id}->{'gCopy'} : 0;
      }
      
      if($gcopy<0.01) {
        push @ecols, 0,0;
      } else {  
        my $xcopy= $tcopy/$gcopy;
        push @ecols, sprintf("%.1f",$xcopy), $gcopy; 
      
        my $gcl= ($gcopy < C0hi)?"0" :(C0hi<=$gcopy and $gcopy<=C1hi)?"1" 
          :(C1hi<$gcopy and $gcopy<=C2hi)?"2-9" :(C2hi<$gcopy and $gcopy<=C3hi)?"10-99"
          :(C3hi<$gcopy and $gcopy<=C4hi)?"100-499":"500+";
        my $ceq= ($xcopy<C0hi)? "xlo" : ($xcopy>C1hi)?"xhi": "xeq"; # add xzero?
        $cnsum{$gcl}{ngene}++; $cnsum{$gcl}{$ceq} ++; $ncnsum++;
        $cnsum{$gcl}{xcopy}+= $xcopy; $cnsum{$gcl}{tcopy}+= $tcopy; $cnsum{$gcl}{gcopy}+= $gcopy;
      }
    }
    $cndid{$id}=$ncp; # or tcopy   
    print $ctabh join("\t",$id,sprintf("%.1f",$tcopy),"$ncp,$nscaf,$ntand",@ecols)."\n";
  }

  sub putcn_missing { our($ctabh,%cnsum,$ncnsum,%cndid,$ncmiss);
    $ncmiss=0;
    if($ngenecov>0) { 
      for my $id (sort keys %{$genecovh}) {
        next if($cndid{$id});
        my($gcopy,$cmed,$cnz,$tlen,$nread,$isuniq);
        if(GENEVEC) {
        ($gcopy,$cmed,$cnz,$tlen,$nread,$isuniq)= @{$genecovh->{$id}};
        } else {
        # switch to named col vals, extend w/ other gene data
        # $gcovtab{$gid}= { gCopy => $tcopy, gLen => $glen, gNread => $nread, gClass => $isuniq, gCM => $cmed, gCnz => $cnz, };
        ($gcopy,$cmed,$cnz,$tlen,$nread,$isuniq)= 
          map{ $genecovh->{$id}->{$_}||0 } qw(gCopy gCM gCnz gLen gNread gClass);
        }
        
        next if($gcopy < C0hi); # == zero in genecov also
        my @ecols=(0, 0, $tlen, 0, $gcopy); # ($saln,$sspan,$tlen, $xcopy, $gcopy);
        print $ctabh join("\t",$id, 0,"0,0,0",@ecols)."\n"; $ncmiss++;

        my $gcl= ($gcopy < C0hi)?"0" :(C0hi<=$gcopy and $gcopy<=C1hi)?"1" 
          :(C1hi<$gcopy and $gcopy<=C2hi)?"2-9" :(C2hi<$gcopy and $gcopy<=C3hi)?"10-99"
          :(C3hi<$gcopy and $gcopy<=C4hi)?"100-499":"500+";
        my $ceq= "xzero";  # add xzero here
        $cnsum{$gcl}{$ceq} ++; #NO: $ncnsum++;
        if(COPYNUMzeros) { $cnsum{$gcl}{ngene}++; $cnsum{$gcl}{gcopy}+= $gcopy; } #<< skip this sum for ave-copyn of nonzero gene set?
        # my $xcopy=0; my $tcopy=0;
        # zero: $cnsum{$gcl}{xcopy}+= $xcopy; $cnsum{$gcl}{tcopy}+= $tcopy; 
      }
    }
    return $ncmiss;    
  }  
  
  sub putcn_sum { our($ctabh,%cnsum,$ncnsum,%cndid,$ncmiss);
    return "Genes Copynum Summary=0" unless($ncnsum>0);
    my @gcl=sort{ $a<=>$b } keys %cnsum;  
    my $ngt=", for nGene=$ncnsum"; $ngt .= " found, $ncmiss missed" if($ncmiss>0);
    my @ceql= qw(xlo xeq xhi); unshift @ceql, "xzero" if($ncmiss>0);
    my $labxloeqhi= ($ncmiss>0) ? "x0,xLo,xEq,xHi= xCopy 0,<1,=1,>1"
      : "xLo,xEq,xHi= xCopy<1,=1,>1";
    (my $labxleh= $labxloeqhi) =~ s/=.*//;
    
    my $cnbrief="Genes Copynum Summary$ngt, xCopy= ";
    print $ctabh "#---------------------------------------\n";
    print $ctabh "# Genes Copynum Summary by Copy level$ngt\n";
    print $ctabh "# tCopy=Chrasm, gCopy=Gene set, xCopy=tCopy/gCopy, $labxloeqhi\n";
    print $ctabh "#".join("\t",qw(CnGene nGene tCopy xCopy gCopy), $labxleh)."\n";
    for my $gcl (@gcl) {
      my $ng= $cnsum{$gcl}{ngene} or next;  # ng includes xzero for COPYNUMzeros, x/t reduced by zeros
      my($xc,$tc,$gc)= map{ sprintf"%.1f",$cnsum{$gcl}{$_}/$ng } qw( xcopy tcopy gcopy) ;
      my($ceqn)= join",", map{ $cnsum{$gcl}{$_}||0 } @ceql; # qw(xlo xeq xhi); #? add xzero
      $gcl="99-499" if($gcl eq "100-499");
      print $ctabh "#".join("\t",$gcl,$ng,$tc,$xc,$gc,$ceqn)."\n";
      $cnbrief .= $gcl."c:$xc/$ng, " unless($gcl eq "0" or $gcl > 100);
    }
    #daphpulex_pa42v2_SRR13333791_b8_mim cnbrief=  
    #Genes Copynum Summary, for nGene=32408 found, 185 missed, xCopy= 1c:1.2/21151, 2-9c:0.8/3711, 10-99c:0.3/142, 99-499c:0.0/18
    print $ctabh "#---------------------------------------\n";
    return($cnbrief);
  }  
  
  while(<$intabh>) {
    next if(/^\W/);
    #o: my($td,$sc,$bs,$ida,$al,$tw,$sw,$tsp,$ssp)=split; 
    my($td,$sc,$bs,$ida,$al,$tw,$sw,$tsp,$ssp,$tdup,$sdup)=split; # FIXME annot btall tab not here
    
    my $nextgn= ($td ne $lasttd)?1:0;
    if(COPYNUM and $nextgn) {
      putcn(@cn);
      @cn=($td, $tw, (0) x 6); # id, Glen, tCopy, Aln, Span, Ncp, Nscaf, Ntand
    }
    $lasttd= $td;

    if(COPYNUM and ($nextgn)) { # or $bs >= $bmaxcut   .. add to @cn
      my($tsb,$tse)= $tsp=~m/^(\d+).(\d+)/;
      my $tspan= ($tse>$tsb)? 1+$tse-$tsb : $al;
      #x $cn[2] += ($tw<1) ? 1 : $tspan/$tw; #tcopy << should be align/tw instead of tspan/tw, tspan is oversized measure
      my $tcopy=  ($tw<1) ? 1 : $al/$tw;
      $cn[2] += $tcopy; # * FIXME tandup needs tcopy addition      
      $cn[3] += $al; 
      $cn[4] += $tspan; 
      $cn[6] ++; $cn[5] ++; # new scanf + all copy
      #?? only one tandup slot?
      if($tdup and $tdup ne "0d") { 
        $cn[7] ++; $cn[5] ++; # tand,all copy
        my($ttsb,$ttse)= $tdup=~m/^(\d+).(\d+)/; 
        my $ttspan= ($ttse>$ttsb)? 1+$ttse-$ttsb : $al;
        $cn[2] +=  ($tw<1) ? $tcopy : _min($ttspan/$tw,$tcopy); # * FIXME tandup needs tcopy addition 
      }
    }
      
  } # close($intabh);

  if(COPYNUM) { 
    putcn(@cn); 
    putcn_missing(); 
    my $cnsum= putcn_sum(); 
    close($ctabh); 
    warn "#$cnsum table to $copytab\n" if($debug);
  }
  
  
}

#-----

sub _max{ return($_[0] < $_[1])?$_[1]:$_[0]; }
sub _min{ return($_[0] > $_[1])?$_[1]:$_[0]; }
sub _minnot0{ return ($_[0] == 0) ? $_[1] : ($_[1] == 0) ? $_[0] : _min(@_); }

=item meanvar

 -- see evigene/scripts/genoasm/cds_meanvar.pl
 -- ** add covar/cor calc for 2+ vars
 -- maybe add transforms: log(x) ; exp(x); sqrt(x) ..
 -- see also cds_meanvar:sub trim()
 
=cut

sub transform_p {
  my($trf, $vals)= @_;
  my @tval;
  if($trf eq "log"){
    @tval= map{ ($_ > 0) ? log($_) : undef } @$vals;
  } elsif($trf eq "sqrt"){
    @tval= map{ ($_ > 0) ? sqrt($_) : undef } @$vals;
  } else {
    @tval= @$vals;
  }
  return @tval;
}

sub backtrans_p {
  my($trf, $vals)= @_;
  my @tval;
  if($trf eq "log"){
    @tval= map{ exp($_) } @$vals;
  } elsif($trf eq "sqrt"){
    @tval= map{ $_ * $_  } @$vals;
  } else {
    @tval= @$vals;
  }
  return @tval;
}

sub skew_p { my($ave,$med,$sd)=@_; my $skew= ($sd<0.001)? 0: 3 * ($ave - $med) / $sd; return $skew; }

=item skew stats

Pearson mode skewness uses the above facts to help you find out if you
have positive or negative skewness. If you have a distribution and you
know the mean, mode, and standard deviation (σ), then the Pearson mode
skewness formula is:

  (mean-mode)/stdev 

if you don’t know the mode, you won’t be able to use Pearson mode
skewness. However, the direction of skewness can be also figured out
by finding where the mean and the median are. According to Business
Statistics, this leads to a second, equivalent formula:

Pearson's second skewness coefficient (median skewness):
  3(Mean – Median) / stdev

Both equations give you results in standard deviations, which are dimensionless units of measurement from the mean.
https://www.statisticshowto.com/pearson-mode-skewness/
  
=cut

sub format_v { local $_= $_[0];
 my $d=($_ >= 100)?0:($_ >= 10)?1:($_ >= 1)?2:($_ >= 0.1)?3:4; $_= sprintf"%.${d}f",$_; $_ =~ s/\.0+$//; return $_; } 
    

sub meanvar_p {
  my($stops, $transform, $vals)= @_;

  use constant { loTRIM => 0.10, hiTRIM => 0.90, };# option? 0.05/0.95 ?
  
  my($med,$nit,$ave,$var,$sd,$se,$s,$range, $skew, $dotrim, $ncut)= (0) x 9;
  $nit= @$vals; 
  return($med,$ave,$nit,$se,$range) if($nit<1);
  
  my @sval= sort { $a <=> $b } @$vals;
  for my $i (0..$#sval) { my $v= $sval[$i]; if( $v =~ s/,.*// ) { $sval[$i]=$v; $ncut++; } }
  if($stops) {
    if($stops =~ /nozero|nonzer/)  { @sval= grep{ $_ > 0 or $_ < 0 } @sval; }
    if($stops =~ /trim/)  { $nit= @sval;
      my $lo= $nit * loTRIM; my $hi= $nit * hiTRIM;
      @sval= splice( @sval, $lo, $hi-$lo) if($nit > 9); 
    }
  }
  
  if($transform) {
    @sval= transform_p($transform,\@sval);
  }
  
  $nit= @sval;
  $med= $sval[ int($nit/2) ]; 
  
  for my $i (0..$#sval) { my $v= $sval[$i]; $s += $v; $var += $v*$v; }
  $ave=$s/$nit; if($nit>2) { 
    $var=($var - $ave*$ave)/($nit-1); $sd=sqrt($var); $se=$sd/sqrt($nit); 
    $skew= ($sd<0.001) ? 0 : 3 * ($ave - $med) / $sd; # DIV 0 err w/ $sd ??? how come
    }

  my($rlo,$rhi)= @sval[0,-1]; # or use 1st/3rd quart?
  # $range="$sval[0],$sval[-1]"; # or 1/3 quartile?  $nit*.25, $nit*.75 ?
  
  if($transform) {
    ($med, $ave, $rlo, $rhi)= backtrans_p($transform, [ $med, $ave, $rlo, $rhi ]);
  }
 
  # map{ my $d=($_ > 99)?0:($_ > 9)?1:2; $_= sprintf"%.${d}f",$_; $_ =~ s/\.0+$//; } ($med,$ave,$sd,$rlo,$rhi,$skew);
  map{ $_= format_v($_) } ($med,$ave,$sd,$rlo,$rhi,$skew);
  
  return($med,$ave,$nit,$sd,$skew,"$rlo,$rhi"); # se > $sd ; changed; added skew
}

sub covar_p {
  my($stops, $transform, $vals, $valb)= @_;

  use constant { loTRIM => 0.10, hiTRIM => 0.90, };# option? 0.05/0.95 ?
  
  my($med,$nit,$ave,$var,$sd,$se,$s,$range, $dotrim, $ncut)= (0) x 9;
  my($nitb,$cor,$covar,$aveb,$sdb,$sb,$varb)= (0) x 9;
  $nit= @$vals; $nitb= @$valb;
  return($nit,$cor,$covar,$ave,$aveb,$sd,$sdb) if($nit<3 or $nitb<3);
  
  my @sval= @$vals; my @svalb= @$valb;
  for my $i (0..$#sval) { my $v= $sval[$i]; if( $v =~ s/,.*// ) { $sval[$i]=$v; $ncut++; } }
  for my $i (0..$#svalb) { my $v= $svalb[$i]; if( $v =~ s/,.*// ) { $svalb[$i]=$v; $ncut++; } }
  # ** need pairwise sort, filters .. assume vala,valb ordered by same gid
  # my @sval= sort {$a <=> $b } @$vals;  my @svalb= sort {$a <=> $b } @$valb;
  if($stops) {
    if($stops =~ /nozero|nonzer/)  { 
      # @sval= grep{ $_ > 0 or $_ < 0 } @sval;
    }
    if($stops =~ /trim/)  { 
      # $nit= @sval;
      # my $lo= $nit * loTRIM; my $hi= $nit * hiTRIM;
      # @sval= splice( @sval, $lo, $hi-$lo) if($nit > 9); 
    }
  }
  
  if($transform) {
    @sval= transform_p($transform,\@sval);
    @svalb= transform_p($transform,\@svalb);
  }
  
  $nit= @sval; $nitb= @svalb;
  if($nit ne $nitb) {
    $nitb= $nit = _min($nit,$nitb); # check if not same, err?
  }
  # $med= $sval[ int($nit/2) ]; 
  
  for my $i (0..$#sval) {
     my $v= $sval[$i]; $s += $v; $var += $v*$v; 
     my $vb= $svalb[$i]; $sb += $vb; $varb += $vb*$vb; 
     $covar += $v * $vb;
  }
  
  $ave=$s/$nit; $aveb= $sb/$nitb;
  if($nit>2){ 
    $covar= ($covar - $ave*$aveb)/($nit - 1 ); #?? check
    $var=($var - $ave*$ave)/($nit-1); $sd=sqrt($var); # $se=$sd/sqrt($nit); 
    $varb=($varb - $aveb*$aveb)/($nitb-1); $sdb=sqrt($varb); 
    my $cden=$sd*$sdb;
    $cor= ($cden < 0.00001) ? 0: $covar/$cden; # DIV0 ERR; is that sqrt( v * vb) ? or ( sd * sdb )?
  }

  # my($rlo,$rhi)= @sval[0,-1]; # or use 1st/3rd quart?
  # $range="$sval[0],$sval[-1]"; # or 1/3 quartile?  $nit*.25, $nit*.75 ?
  
  if($transform) {
    ($ave, $aveb)= backtrans_p($transform, [ $ave, $aveb ]);
  }
 
  #map{ my $d=($_ >= 100)?0:($_ >= 10)?1:($_ >= 1)?2:($_ >= 0.1)?3:4; $_= sprintf"%.${d}f",$_; $_ =~ s/\.0+$//; } 
  #  ($cor,$covar,$ave,$aveb,$sd,$sdb);
  map{ $_= format_v($_) } ($cor,$covar,$ave,$aveb,$sd,$sdb);
  
  return($nit,$cor,$covar,$ave,$aveb,$sd,$sdb); 
}


# sub meanvar {
#   my(@vals)= @_;
# 
#   use constant { loTRIM => 0.10, hiTRIM => 0.90, };# option? 0.05/0.95 ?
#   my($med,$nit,$ave,$var,$sd,$se,$s,$range)= (0) x 9;
#   $nit= @vals; 
#   return($med,$ave,$nit,$se,$range) if($nit<1);
#   
#   my @sval= sort {$a <=> $b} @vals;
#   $med= $sval[ int($nit/2) ]; 
#   
#   my $nitall=$nit; #??
#   if($TRIMAVE and $nit>9) {
#     my $lo= $nit * loTRIM; my $hi= $nit * hiTRIM;
#     @sval= splice( @sval, $lo, $hi-$lo);
#     $nit= @sval;
#   }
#   for my $i (0..$#sval) { my $v= $sval[$i]; $s += $v; $var += $v*$v; }
#   $ave=$s/$nit; if($nit>2){ $var=($var - $ave*$ave)/($nit-1); $sd=sqrt($var); $se=$sd/sqrt($nit); }
# 
#   my($rlo,$rhi)= @sval[0,-1]; # or use 1st/3rd quart?
#   # $range="$sval[0],$sval[-1]"; # or 1/3 quartile?  $nit*.25, $nit*.75 ?
#   map{ my $d=($_ > 99)?0:($_ > 9)?1:2; $_= sprintf"%.${d}f",$_; $_ =~ s/\.0+$//; } ($med,$ave,$sd,$rlo,$rhi);
#   
#   return($med,$ave,$nit,$sd,"$rlo,$rhi"); # se > $sd ; changed
# }

sub readAnnots {
  my($anntable)= @_;
  my($nann,$lpt)=(0,0); 
  warn "#read anntable=$anntable\n" if($debug);
  open(F,$anntable); 
  while(<F>){
    next if(/^\W/); my @v= split;
    my($cid,$cib,$an,$ids)= @v; # any more than an,ids ?
    next unless($an or $ids);
    $annvals{$cid}{$cib}="$an\t$ids"; $nann++;
  } close(F);
  return($nann);
}


sub read_idclass { # for case of no anntable, but idclass, eg CDS.covtab
  my($crclassf,$allclasses,$crlenh)=@_; 
  my($nid,$ncl,$ok,$inh)=(0,0,0);
  my %crclass=();   my @crclass=();
  $allclasses||=0;

  use constant kMAXIDCLASS => 9; # idclass limit, using 3-4 now
  my $CRTPAT=''; # no default, see sam2covtab
  
  if($crclassf and ($ok,$inh)= openRead($crclassf) and $ok) { 
    while(<$inh>){ next if(/^\W/); 
      my($cr,$crclass,@clx)=split;  # may have more columns .. keep all ie CDS,BUSCO 
      if($allclasses and @clx){ $crclass=join(" ",$crclass,@clx); }
      $crclass{$cr}= $crclass || 0;  
    } close($inh); 
    $nid= scalar(keys %crclass);
  }
  #.. insert other way from sam hdr ids and CRTPAT?
  if($nid==0 and $CRTPAT and ref($crlenh)){
    for my $cr (sort keys %$crlenh) {
      my($crt)= ($cr =~ m/($CRTPAT)/)? $1 : 'UNK'; 
      $crclass{$cr}= $crt;
    }
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
    warn "#ERR too many idclasses n=$ncl from nid=$nid of $crclassf or sam ids x CRTPAT='$CRTPAT' \n";
    $nid=0; %crclass=(); @crclass=();
  }
  warn "# read nid=$nid, nclass=$ncl from $crclassf\n" if($debug);
  return($nid,\%crclass,\@crclass);
}


sub readMetad {
  my($sampledata)= @_;
  my($nsam,$aid)=(0,0); #? only 1 aid, should have default
  if($sampledata and -f $sampledata) {
    warn "#read sampledata=$sampledata\n" if($debug);
    open(F,$sampledata); # 
    while(<F>){
      next if(/^\W/);
      my($key,$val)= (m/^(\w+)\s*=\s*(.+)$/) ? ($1,$2):(0,0);
      next unless($key);
      
      if($key eq 'pt' or $key eq 'asmid') { $aid=$val; $nsam++;
        map{ $samvals{$aid}{$_}= 0; } qw(flowcyto cytomb asmtotal atotalmb asmname );
    
      } elsif($key eq 'flowcyto') { # /^flowcyto=(.+)$/)  
        # my $flowcyto= $val; # flowcyto=234-391 Mb
        $samvals{$aid}{$key}= $val;
        $samvals{$aid}{cytomb}  = ($val=~m/(\d+)/)?$1:0; # can be range 160-180
    
      } elsif($key eq 'glncformula') { # /^glncformula=(.+)$/ 
        $samvals{$aid}{$key}= $val; # glncformula=nnn Mb
        
      } elsif($key eq 'asmtotal') { # /^asmtotal=(.+)$/)   
        my $asmtotal= $val; 
        $samvals{$aid}{$key} = $asmtotal;
        $samvals{$aid}{atotalmb} = ($asmtotal=~m/(\d+)/)?$1:0; # can be range 160-180
    
      } elsif($key eq 'kulo' or $key eq 'kuhi' or $key eq 'asmname'){ 
        $val=~s/\s.*//; $samvals{$aid}{$key}= $val;  # /^kulo=(\S+)/)  # test hi<lo
        
      } else { # save all valid key = val
        $val=~s/\s.*//; $samvals{$aid}{$key} = $val;
      }
      
    } close(F);
  }
  
  # allow ENV{key} vals ?
  for my $key (qw(asmid asmname asmtotal flowcyto kulo kuhi kucg)) {
    if(my $ev= $ENV{$key}) { 
      if($key eq "asmid"){ $aid=$ev; $nsam++; } 
      $samvals{$aid}{$key}= $ev if($aid); 
      }
  }
  
  return($nsam,$aid); # return ($nsam, \%metavals); #?
}

sub openOut {
  my($fn)=@_; my($ok,$fh)=(0); 
  return($ok,undef) unless($fn and $fn=~/\w/);
  rename($fn,"$fn.old") if( -f $fn );
  $ok= open($fh,'>',$fn); unless($ok){ warn "# FAIL writing $fn\n"; }
  return($ok,$fh);
} 

sub openRead {
  my($fn)=@_; my($ok,$inh)=(0); 
  return(0,undef) unless($fn and -f $fn);
  if($fn =~ /\.gz/) { $ok= open($inh,"gunzip -c $fn |"); } else { $ok= open($inh,$fn); }
  # die "cant read $fn" unless ($ok); # warn? leave to caller
  unless($ok){ warn "# FAIL reading $fn\n"; }
  return($ok,$inh);
}

__END__



=item try2out

$evigene/scripts/genoasm/gnodes_sumgenecov.pl -title at18chr_test8f2sumgcn  -genexcopy arath18tair1cds_SRR10178325_b2_mim.geneycopymd -chrgenetab arath18tair_chr_SRR10178325_test8f.genetab  -asmid arath18tair_chr -meta arath20asm.metad -idclass arath18tair1cds.idclass 

#genecov output to sum=at18chr_test8f2sumgcn_genesum.txt means=at18chr_test8f2sumgcn.genemeans
#read sampledata=arath20asm.metad
# read nid=27444, nclass=2 from arath18tair1cds.idclass
#read_genescov: UCG med=32.7, sd=1.06, ave=32.9, n=27368
#read_genescov: ngene=27368, nzero=77 (26615 1-copy, 716 2+copy) from arath18tair1cds_SRR10178325_b2_mim.geneycopymd
#readGenetab: n=27380, 
#sumgene_tabout: n=27380 at18chr_test8f2sumgcn.sumgenetab
#sumgene_sumout: n=44 at18chr_test8f2sumgcn.genesum.txt
dgmm:aweed20gnodes:% head at18chr_test8f2sumgcn.sumgenetab

head at18chr_test8f2sumgcn.sumgenetab
GeneID          gLen    gNread  gCopy   gCM     gCnz                  gClass  aCopy   aCovT   aCovM   aCovZ   mCopy   mCovT   aCopyRd aReadT  aReadM  aReadZ
AT1G01010t1     1290    364     1.0     32,1.0  32.0,33,1290,100.0      uniq    1.0     22      22      0       1.0     24      1.0     555     555     0
AT1G01020t1     576     180     0.9     28,0.9  33.0,31,485,84.2        uniq    1.0     18      18      0       1.0     20      1.0     279     279     0
AT1G01030t1     1008    303     1.3     40,1.2  40.5,41,1008,100.0      uniq    1.0     32      32      0       1.0     34      1.0     461     461     0
AT1G01040t1     5733    1555    1.0     34,1.0  34.2,33,5733,100.0      uniq    1.0     30      30      0       1.0     30      1.0     2357    2357    0
AT1G01050t1     639     197     0.9     27,0.8  30.5,30,572,89.5        uniq    1.0     22      22      0       1.0     25      1.0     303     303     0
AT1G01060t1     1938    543     1.1     36,1.1  36.4,37,1901,98.1       uniq    1.0     29      29      0       1.0     29      1.0     819     819     0
AT1G01070t1     1098    366     1.1     38,1.2  38.3,37,1098,100.0      uniq    1.0     26      26      0       1.0     28      1.0     556     556     0
AT1G01080t1     885     190     0.8     26,0.8  28.0,27,828,93.6        uniq    1.0     21      21      0       1.0     23      1.0     292     292     0
AT1G01090t1     1287    284     0.9     29,0.9  28.8,29,1287,100.0      uniq    1.0     25      25      0       1.0     28      1.0     431     431     0

at18chr_test8f2sumgcn.genesum.txt
# ALL genes ----------------------
Field   Median,Mean,N,StDev,Range________
gLen    1005    1187    27380   1483    0,16203
gNread  273     415     27380   7208    0,846276
gCopy   1       1.28    27380   9       0,683
aCopy   1       1.10    27380   1.83    0,167
mCopy   1       1.16    27380   2.99    0,167
aCopyRd 1       1.09    27380   1.74    1,166
gCM     33      41.9    27380   478     0,46641
aCovT   24      27.8    27380   136     0,10026
mCovT   26      26.6    27380   41.5    0,3414
aReadT  432     3607    27380   273060  2,30409470
aCovZ   0       0       27380   0       0,0

# gCopy => [1.75,99] genes ----------------------
Field   Median,Mean,N,StDev,Range________
gLen    582     813     705     1086    78,5730
gNread  288     2112    705     30954   20,640809
gCopy   2.60    5.15    705     8.84    1.80,92.8
aCopy   2       4.02    705     7.10    1,38
mCopy   2       6.07    705     17.4    1,167
aCopyRd 2       3.69    705     6.26    1,24.1
gCM     48      165     705     1629    12,39122
aCovT   46      104     705     484     7,10026
mCovT   45      74.4    705     112     1,498
aReadT  1716    85816   705     1476474 87,30409470
aCovZ   0       0       705     0       0,0

# gCopy => [0.75,1.50] genes ----------------------
Field   Median,Mean,N,StDev,Range________
gLen    1044    1225    25654   1516    66,16203
gNread  281     334     25654   541     10,50654
gCopy   1       1.02    25654   1.03    0.80,1.50
aCopy   1       1.01    25654   1.02    1,7.50
mCopy   1       1.03    25654   1.11    0,38
aCopyRd 1       1.01    25654   1.02    1,7
gCM     33      33      25654   55.3    4,6773
aCovT   24      24.4    25654   46.2    3,5771
mCovT   26      25.5    25654   26.4    0,56
aReadT  434     534     25654   1760    20,242402
aCovZ   0       0       25654   0       0,0

# gCopy => [0.0,0.50] genes ----------------------
Field   Median,Mean,N,StDev,Range________
gLen    225     237     84      282     0,636
gNread  15      21.1    84      30      0,94
gCopy   0.40    0.35    84      0.39    0,0.50
aCopy   1       1.10    84      1.29    0,6
mCopy   1       0.88    84      0.95    0,2
aCopyRd 1       1.04    84      1.05    1,2.40
gCM     7       7.94    84      9.9     0,20
aCovT   8       7.71    84      8.60    0,15
mCovT   5       5.71    84      7.23    0,17
aReadT  38      46.9    84      59.4    2,171
aCovZ   0       0       84      0       0,0

=cut


