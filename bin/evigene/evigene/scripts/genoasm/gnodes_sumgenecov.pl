#!/usr/bin/perl
# gnodes_sumgenecov8j.pl

use strict;
use Getopt::Long;  

use constant UPD21AUG => 1; # fixes, add aLen, assembly length vs gLen as other copy num
use constant UPD21SEP => 1;
use constant UPD21OCT => 1; # UPD21OCT test opts: 
#  env gidclean=1 : cleaner gene-cds subset by glen, grds, no skew/zero : ** MAKE THIS DEFAULT **
#  env xcutwide=1 : XClo, XChi=0.52,1.90 opt, vs narrower 0.66,1.55, ie more eq vals: any change to chrasm order?
#  SWAP_AC_COPY=1 : cCopy better than aCopy now, USE_mCopy4aCopy=0: median val not clear than ave, drop
#  CCOV_CDSFILT=1 : filter cr: covtab cCopy data by top 50% of cds crvals, to drop out spurious dups, more reliable cxCopy

use constant LLOCI_CALC => 2; # UPD21SEP16, LLOCI_CALC==2 now best match to gCopy, swap lloci_calc2 for calc1 (both out for now)
use constant kSAMPLE => 1000;  # CBIN size median sample @cbin, should be all same = default $CBIN

my $debug=$ENV{debug}||1;
my @IVAR=(5,6,7); # covtab cols ==  rd_total, rd/mmap, rd_uniq for chrasm x rd map 
my $DOMED= 0; # $ENV{median}||0; * SUCKS UP ALL covtab vals, mem pig
my $DOHIST=1; # replace median from @allvals w/ histo mode(s)/peak, est median, counts for depth vals depth[0..999]?
my $ZEROS= 0; # $ENV{zero}||0;
my $showSUM=1; # $ENV{sum}||0;
my $MINRD= 1; #? skip all cov rows w/ acovt < MINRD, want to measure zeros? gaps?
my $CBIN=100; # calc from rows ib dist
my $TOTALID="total"; #? "allgenes";
my $KU_OPT=0; # read from metad, or calc from genemeans, unless -CU|kucg=val
my $XHICUT= $ENV{xhicut}|| 1.7; my $XLOCUT= $ENV{xlocut}|| 0.55;# xCopy opts, change def? xhi=1.60 .. 1.66? xlo=0.50..0.54?
my $TRIMAVE= $ENV{trimave}||0; # for ave calc, trim extremes?
my $STOPTS=$ENV{'statopt'} || ""; # opt for sumgene_sumout : nozero, trim
my $TRANS= $ENV{'transform'}||""; # transform for sumgene_sumout log, sqrt, ..
my $NOTABO= $ENV{notabout} || 0;
my $Genetab8i=2; # $ENV{'genetab8i'} || 0; # UPD21AUG gnodes_sam2covtab8i.pl; 2==8j
use constant DO_RPLOT => UPD21OCT;    
my($DONEG,$DOMISR,$DO_AWEIGHT,$DO_GWEIGHT) = ( 1 ) x 9; # DO_RPLOT defaults
  $DOMISR= 2; #update method
my $USE_GIDCLEAN= 1; # $ENV{gidclean}||0; # make default w/ -nogidclean opt
my $MINGLEN= 200; #FIXME $ENV{minglen} 
my $Lrdlen=$ENV{readlen}||0; # global?
my $HAS_CDSCOV= 0; my $cdscovh; # ($HAS_CDSCOV,$cdscovh)= readCovtab($cdscovtab); 

my $intitle= $ENV{title}||$ENV{name}||""; # was "asource"; 
my $asmid= $intitle; #? intitle for output, asmid for data
my ($sampledata,$anntable,$outmeans,$outsum,$ncpu,$crclassf,$readsetid)=("") x 9;
my ($cdscovtab,$chrcovtab,$genexcopy,$chrsumtext); # input data
my (@lvar, @genextop, @CBIN);

my $optok= GetOptions( 
  'cdscovtab=s',\$cdscovtab, # was UNUSED
  'chrgenetab=s',\$chrcovtab, # was chrcovtab=
  'genexcopy=s',\$genexcopy, # from sam2genecov
  'chrsummary=s',\$chrsumtext, # name_sum.txt, DO_RPLOT eg arath18tair_chr_bwatest8f_SRR10178325_sum.txt

  'output|outsummary=s', \$outsum,
  'means=s', \$outmeans,
  'title|name=s',\$intitle,
  'asmid=s', \$asmid, # UNUSED
  'readsetid=s', \$readsetid,   
  'genomedata|sumdata|sampledata|metadata=s', \$sampledata,  # UNUSED
  'anntable|annotations=s', \$anntable,  # UNUSED
  'crclassf|idclassf=s',\$crclassf, # UNUSED table of chr/gene id => class ? want for BUSCO, TEfam, : No, rely on anntable having these
  'minread=i', \$MINRD, # 5 default, test, see if it affects over-assembly (Xcopy < 1)
  'CU|kucg=s', \$KU_OPT,
  'XHICUT=s', \$XHICUT, 'XLOCUT=s', \$XLOCUT, # xCopy hi/lo cutoff from 1.0
  'statopt=s', \$STOPTS, 'transform=s', \$TRANS,
  'GIDCLEAN!', \$USE_GIDCLEAN,  'MINGLEN=i', \$MINGLEN, # for gidclean
  'noTABOUT!', \$NOTABO, 
  'ncpu=i', \$ncpu, # not used here
  'debug!', \$debug, 
  );

die "usage: gnodes_sumgenecov.pl -chrgenetab xxx_chr.genetab  -genexcopy xxx_cds.genexcopy 
 opts: -name xxxcdschr  -CU=$KU_OPT -minread=$MINRD -debug 
" unless($optok and ($cdscovtab or $chrcovtab or $outmeans));
# -idclass xxx.idclass  -cdscovtab xxx_cds.covtab  -anntab xxx.anntab   -metadata xxx.metad -asmid xxx_asm 

if($asmid and not $intitle) { $intitle=$asmid; } # want one of these
elsif($intitle and not $asmid){ $asmid=$intitle; }
unless($readsetid) {
  my $rdset=""; # FIXME leads to messy fname if no SRRnnnn in covtabs
  for my $incov ($cdscovtab,$chrcovtab,@ARGV) {
    if($incov and $incov =~ m/([A-Za-z0-9]+).([A-Za-z0-9]+)\.(covtab|genetab)/) { 
      my($asmrd,$vers)=($1,$2); 
      $asmrd =~ s/($asmid|$intitle)//g;
      if($asmrd =~ m/([A-Z]RR\d+)/) { my $na=$1;  $rdset .= $na unless($rdset=~m/$na/); last; }
      # ^^ drops SRRnnn{ab} suffix
      elsif($asmrd=~/([A-Za-z0-9]+)/){ my $na=$1;  $rdset .= $na unless($rdset=~m/$na/); last; }
    }
  }
  if($rdset =~ m/\w\w/){ $rdset=substr($rdset,0,19) if(length($rdset>19)); $readsetid=$rdset; }
} elsif($readsetid =~ /^(no|none|0)/) {
  $readsetid="";
}

# my $title= $intitle;
unless($outsum) { 
  ($outsum=$intitle) =~ s/\W/_/g; # title always ok?
  $outsum.="_".$readsetid if($readsetid and $outsum !~ m/$readsetid/);
  $outsum.="_genesum.txt";
}

unless($outmeans){ ($outmeans=$outsum) =~ s/.genesum.*//; $outmeans.=".genemeans"; }

unless($chrsumtext) {
  # arath18tair_chr_SRR10178325ab_test8jgbwa.genetab => arath18tair_chr_test8jgbwa_SRR10178325ab_sum.txt
  # dang name form change: asmid,srrid,title > asmid,title,srrid
  (my $cst=$chrcovtab) =~ s/.genetab/_sum.txt/;   
  if( -f $cst) { $chrsumtext=$cst; }
  else {
    if($readsetid and $cst =~ m/$readsetid/) { # ugh..
      if($cst =~ s/([_]?${readsetid}[a-z]*)//){ my $rsd=$1; $cst =~ s/_sum.txt/${rsd}_sum.txt/;
      if( -f $cst) { $chrsumtext=$cst; } 
      }
    }
  }
}  

sub MAINstub {}
warn "#genecov output to sum=$outsum\n" if $debug; #  means=$outmeans

# read_genexcopy global gene copynum table
my($gucg_n, $gucg_med, $gucg_ave, $gucg_err, $ngenecov, $genecovh) =  
      ($genexcopy)? read_genexcopy($genexcopy): (0,0,0,0,0); 

($HAS_CDSCOV,$cdscovh)= ($cdscovtab) ? readCovtab($cdscovtab) : (0); 
 
my($ncchr)= ($chrcovtab) ? readGenetab($chrcovtab, 0) : 0; # primary data + genexcopy

my($ncrsum,$crsumh)= ($chrsumtext)? read_crsumtext($chrsumtext) : 0; # DO_RPLOT only ?

#...... new outputs .......
my $outnam=$outsum; $outnam=~s/.genesum.*//;

my @gids= sort keys %$genecovh;
my @gidnoz= grep{ not ($genecovh->{$_}->{gClass} =~ m/zero|^0/ or $genecovh->{$_}->{gNread} < $MINRD ) } @gids;
my @gidzero= grep{ ($genecovh->{$_}->{gClass} =~ m/zero|^0/ or $genecovh->{$_}->{gNread} < $MINRD ) } @gids;

## add @gidclean subset? for better stats, remove gLen < 200|MINGLEN, zero, skew?, gNread < MINRD
#up my $USE_GIDCLEAN= $ENV{gidclean}||0; # make default w/ -nogidclean opt
#above my $MINGLEN= $ENV{minglen}||200; #FIXME
my $min2read= ($MINRD < 10) ? 10 : $MINRD;
my @gidclean= grep{ not ($genecovh->{$_}->{gClass} =~ m/skew/ 
  or $genecovh->{$_}->{gLen} < $MINGLEN or $genecovh->{$_}->{gNread} < $min2read ) } @gidnoz;

my $gidnoz= ($USE_GIDCLEAN) ? \@gidclean : \@gidnoz; # this affects all stats below
my $ngzero= scalar(@gids) - scalar(@$gidnoz); #was  @gidzero;

#UPD21AUG20, aLen/aBins output: aLen= aBins * $CBIN ** This is over-estimate from CBIN slop, any correction?
for my $gid (@gids) { my $nb= $genecovh->{$gid}->{'aBins'}||0; $genecovh->{$gid}->{'aLen'}= $nb*$CBIN; }

#UPD21SEP08 ** Output cols add lCopy, lLoci .. lCopy2 testing SEP14
# 21SEP16: swap lCopy1/2, calc2 is best now .. old is lCopy2 now, drop
#UPD21SEP22: ? add new est lCopy * lCov/CU , should be closest to gCopy if measures are accurate
# .. also swap out aCopy of lCopy, lxCopy=lCopy*lCov/CU?, for xCopy in Copynum Summary by Copy level ?

use constant USE_lCopy_xCopy => 2; # UPD21SEP22: add lxCopy, calc xCopy from lCopy/gCopy, and? lxCopy/gCopy
use constant USE_mCopy4aCopy => 0; # was =1; # UPD21OCT .. mCopy not helpful
use constant SWAP_AC_COPY => UPD21OCT; # UPD21OCT : cCopy is better than aCopy now ??
use constant PLOT_GENEBASES => 0; #Not UPD21OCT; # upd21oct31, add/replace rplot of % ngene w/ % gene mbases as better view
    #^^ doesnt look good, excess weight to uniq set w/ over-estimate aCopy/Copy by a bit.
    # .. drop PLOT_GENEBASES unless resolve problems.
    
# add similar to aCopy? axCopy = aCopy * aCovM / CUest ?     
## sumgene_counts cutoffs for copynum levels
use constant { C0hi=>0.66, C1hi=>1.55, C2hi=>9.99, C3hi=>99.9, C4hi=>499 }; # sumgene_counts xCopy = aCopy/gCopy cuts
use constant { XClo => 0.66, XChi => 1.55 }; # now same as copynum C0hi , C1hi , but for xCopy=iCopy/gCopy lo/eq/hi class
  # more permissive: { XClo => 0.52, XChi => 1.90 }; .. change to option
my($XClo,$XChi)= (XClo, XChi); # .. change to option, what?
if(my $xcut= $ENV{xcutwide}) { ($XClo,$XChi)= (0.52,1.90); }

## UPD21OCT22: test use aReadM to replace cCopy/cCovM, bad est for uniq crmap != cdsmap, eg human
## .. drop this as aReadM over-estimates CovT, see eCopy now works w/ human & daphpulex ??
#See human20gnodes_test8j_genecopy.info and dplx20gnodes_test8j_genecopy.info
# my $USE_RDMCOV = $ENV{rdmcov}||0;
# if($USE_RDMCOV and $Lrdlen > 0) {  # must precede SWAP_AC_COPY and USE_lCopy_xCopy xxCopy
#   warn "#USING rdmcov, readgene CovT est = readM * $Lrdlen / genelen\n" if($debug);
#   for my $gid (@gids) { 
#     my($glen, $rdt,$rdm, $origcc,$ct,$cm)= map{ $genecovh->{$gid}->{$_}||0 } qw(gLen aReadT aReadM cCopy cCovT cCovM);
#     # global $Lrdlen= # read len, where? read_genexcopy has it
#     # below: $origcc= ($cm<1)?0:sprintf"%.1f",$ct/$cm; 
#     if($glen > 0 and $rdm > 0 and $cm > 0) {
#       my $estCovT= $rdm * $Lrdlen / $glen;
#       my $estCC= sprintf"%.1f",$estCovT/$cm;
#       $genecovh->{$gid}->{'cCopyOld'}= $origcc;
#       $genecovh->{$gid}->{'cCopy'}= $estCC;
#       $genecovh->{$gid}->{'cCovTOld'}= $ct;
#       $genecovh->{$gid}->{'cCovT'}= sprintf"%.1f", $estCovT;
#     }
#   }
# }


        
our $CUest= $KU_OPT || $gucg_med;

# sub make_axcopy{ our($CUest);
#   my($axtag,$gid,$acopy,$acovm)=@_; 
#   map{ s/,.*//; } ($acopy,$acovm); $acopy ||= 0;
#   my $axcopy= sprintf "%.1f", $acopy * $acovm / $CUest;  
#   $genecovh->{$gid}->{axtag} = $axcopy; # ugh, output in sumgene_tabout(), count in sumgene_counts()
# }

if(USE_lCopy_xCopy  and $CUest >= 1) {
  for my $gid (@gids) { 
    my $lcopy= $genecovh->{$gid}->{'lCopy'}||0; # now tuple ? av,md,pct
    my($lcopya,$lcopym)=split",",$lcopy; $lcopym ||= $lcopya;
    my $lcov= $genecovh->{$gid}->{'lCov'}||0;  # FIXME: lcov is tuple: 407.9,31,100%, use median=2
    my($lcova,$lcovm)=split",",$lcov; $lcovm||=$lcova;
    my $lxcopy= sprintf "%.1f", $lcopya * $lcovm / $CUest; # dang, not sure this is * or +
    $genecovh->{$gid}->{'lxCopy'} = $lxcopy; # ugh, output in sumgene_tabout(), count in sumgene_counts()
#     make_axcopy('lxCopy', $gid, $lcopya, $lcovm);
    
    if(UPD21OCT) { # see new cCopy, cCovM for aCovM
    my $acopy= (USE_mCopy4aCopy) ? $genecovh->{$gid}->{'mCopy'} : $genecovh->{$gid}->{'aCopy'};   
    my $acovm= $genecovh->{$gid}->{'aCovM'}||0;   # Dang, this from sam2covtab8j is not equiv to gCov/UCGcov, is over-est
#     make_axcopy('axCopy', $gid, $acopy, $acovm);
    map{ s/,.*//; } ($acopy,$acovm); $acopy ||= 0;
    my $axcopy= sprintf "%.1f", $acopy * $acovm / $CUest;  
    $genecovh->{$gid}->{'axCopy'} = $axcopy; # ugh, output in sumgene_tabout(), count in sumgene_counts()
    }
    
    if(UPD21OCT) {
    my $acopy= $genecovh->{$gid}->{'cCopy'}||0;   
    my $acovm= $genecovh->{$gid}->{'cCovM'}||0;   # Dang, this from sam2covtab8j is not equiv to gCov/UCGcov, is over-est
#     make_axcopy('cxCopy', $gid, $acopy, $acovm);
    map{ s/,.*//; } ($acopy,$acovm); $acopy ||= 0;
    my $axcopy= sprintf "%.1f", $acopy * $acovm / $CUest;  
    $genecovh->{$gid}->{'cxCopy'} = $axcopy; # ugh, output in sumgene_tabout(), count in sumgene_counts()
    }
    if($HAS_CDSCOV){ # eCopy
    my $acopy= $genecovh->{$gid}->{'eCopy'}||0;   
    my $acovm= $genecovh->{$gid}->{'eCovM'}||0;   # Dang, this from sam2covtab8j is not equiv to gCov/UCGcov, is over-est
#     make_axcopy('exCopy', $gid, $acopy, $acovm);
    map{ s/,.*//; } ($acopy,$acovm); $acopy ||= 0;
    my $axcopy= sprintf "%.1f", $acopy * $acovm / $CUest;  
    $genecovh->{$gid}->{'exCopy'} = $axcopy; # ugh, output in sumgene_tabout(), count in sumgene_counts()
    }
    }
}

if(SWAP_AC_COPY) { 
#SWAP_AC_COPY: replace aCopy w/ cCopy as best est; ** too confusing, do this only for output sumtab cols?
# .. for public use remove poor stats, eCopy now best, keep aC? cC? 
my @ackeys= qw(aCopy aCovT aCovM axCopy cCopy cCovT cCovM cxCopy);
# UPD21OCT31: $ENV{acopyis} == 'e' now default best case, allow acopyis=a to override?
my $ack= $ENV{acopyis} || 'e'; if($ack) { map{ s/^c/$ack/ } @ackeys; } # ie, cCopy => eCopy
for my $gid (@gids) { 
  my($ac,$at,$am,$ax, $cc,$ct,$cm,$cx)= map{ $genecovh->{$gid}->{$_}||0 } @ackeys;
  @{$genecovh->{$gid}}{@ackeys} = ($cc,$ct,$cm,$cx, $ac,$at,$am,$ax); #?? perl ok
}
}

#UPD21OCT: drop aLen as useless?
my @FLDO= qw( gLen gNread gCopy gCM gCnz gClass aCopy aCovT aCovM aCopyRd aReadT aReadM aReadZ lCopy lCov  lLoci );
# if(exists $genecovh->{$gids[0]}->{'aLen'}) { push @FLDO, 'aLen'; }
push @FLDO, 'lxCopy' if(USE_lCopy_xCopy);
push @FLDO, 'axCopy' if(UPD21OCT and USE_lCopy_xCopy);
push @FLDO, qw(cCopy cCovT cCovM cxCopy)  if(UPD21OCT); #debug 
push @FLDO, qw(eCopy eCovT eCovM exCopy)  if($HAS_CDSCOV); #ddd 
unless(UPD21OCT) {
if(exists $genecovh->{$gids[1]}->{'lCopy2'}) { push @FLDO, qw(lCopy2 lCov2 lLoci2); } #? lCov/lCovM ?
}


if(UPD21OCT) { 
# add mCopy median to aCopy for tabout? USE_mCopy4aCopy?
# add eqToGcopy/eqCopy col as per debug tests: a:lo,r:lo,l:lo,xl:eq,4:lo, for gCopy x aCopy,aCopyRd,lCopy,lxCopy + others
my @cpfld= grep { m/Copy/ and not m/^(gCopy|aCopyRd)/ } @FLDO;
my @cplab= map { (my $l=$_) =~ s/Copy//; $l; } @cpfld;
push @FLDO, 'eqToGcopy';
for my $gid (@gids) { # some are empty, no data == 0
  my @cpv= map{ $genecovh->{$gid}->{$_}||0 } @cpfld;
  my $gc= $genecovh->{$gid}->{'gCopy'};
  my $eqgc= eqcopy($gc,\@cpv,\@cplab); # eqToGcopy col ~~  a:lo,r:lo,l:lo,xl:eq,4:lo
  $genecovh->{$gid}->{'eqToGcopy'}= $eqgc;
  #--
  if(0) { # USE_mCopy4aCopy ? turn off, dont want both now
  my($ac,$mc)= map{ $genecovh->{$gid}->{$_}||0 } qw( aCopy mCopy); 
  $genecovh->{$gid}->{'aCopy'} = "$ac,$mc" unless($ac == 0 and $mc == 0); 
  }
  }
}

unless($NOTABO) {
  my($ntabo,$tabo)= sumgene_tabout(\@gids, \@FLDO,  $outnam.".sumgenetab");
  warn "#sumgene_tabout: n=$ntabo (nzero=$ngzero) $tabo\n" if $debug;
}

my($stitle,$nsum)=(0,0);
my($ook,$outsumh)= openOut($outsum);

# FIXME2: add gClass counts table, in sumgene_counts()? or end w/ @genextop? ie: uniq:50%,12000 dupx:25%,6000 ..
# FIXMEd: gCopy cuts for counts should ~ match gCopy cuts for moments
#above: use constant { C0hi=>0.66, C1hi=>1.55, C2hi=>9.99, C3hi=>99.9, C4hi=>499 }; # xCopy = aCopy/gCopy cuts


## insert here: Counts for genes cover depth on chromosome assembly
$stitle="# Chromosome Assembly Gene-Copynum Counts by Copy level, Chrasm=$asmid";
$stitle .= " (Clean gene subset)" if ($USE_GIDCLEAN);
my($nco,$cnbrief,$ngrow,$growh)= sumgene_counts(  $outsumh, $stitle, \@gids, $gidnoz, \@gidzero);
$nsum+=1;  warn "#$cnbrief\n" if($debug);

# DO_RPLOT here? in sumgene_count? uses sumgene_counts and $ncrsum,$crsumh
my $outrplot="";
if(DO_RPLOT and ($ncrsum or $ngrow)) {
  my $plname=$intitle || $asmid; # want short name here
  ($outrplot=$outsum) =~ s/.genesum.*//; $outrplot.=".sumgene.rplot";  # ? .R or .Rplot or .Rscript ? 
  # plot.pdf == asmgcn_miss_$title.pdf ? should be outrplot.pdf
  my($rok)= sumgene_rplot($outrplot, $plname, $crsumh, $growh);
  warn "#sumgene_rplot $outrplot\n" if $debug; #  means=$outmeans
  unless($rok){ warn "#sumgene_rplot error\n"; $outrplot=""; }
}

my @FLDSUM= qw( gLen gNread gCopy aCopy aCopyRd  gCM aCovT aCovM aReadT aReadZ lCopy lCov); #? drop:aLen; aCovZ or aReadZ?
push @FLDSUM, qw(cCopy cCovT cCovM axCopy lxCopy cxCopy)  if(UPD21OCT); #add, debug 
push @FLDSUM, qw(eCopy eCovT eCovM exCopy) if($HAS_CDSCOV);

my %FLDCOR= ( gCopy => "gLen,gNread,lCopy", #? lCopy: yes
    aCopy => "gCopy,gLen,gNread,lCopy", aCopyRd => "gCopy,gLen,gNread,lCopy", 
    aCovT => "gCM,gLen,gNread",  aReadT => "gCM,gLen,gNread",
    ); #  mCopy => "gCopy,gLen,gNread", mCovT => "gCM,gLen,gNread",

$stitle="# Moments for all genes";
$nsum += sumgene_sumout( \@gidnoz, \@FLDSUM, $outsumh, $stitle, \%FLDCOR); # call several times, one output file

my $gidsel= sumgene_select( \@gidnoz, { gCopy => [C1hi,C3hi] }); # 1.75,99
$stitle= join" ","# Moments for gCopy => [", C1hi,C3hi,"] genes";
$nsum += sumgene_sumout( $gidsel, \@FLDSUM, $outsumh, $stitle, \%FLDCOR); # call several times, one output file

my $gidsel= sumgene_select( \@gidnoz, { gCopy => [C0hi,C1hi] }); # 0.75,1.50
$stitle= join" ","# Moments for gCopy => [", C0hi,C1hi,"] genes";
$nsum += sumgene_sumout( $gidsel, \@FLDSUM, $outsumh, $stitle, \%FLDCOR); # call several times, one output file

my $gidsel= sumgene_select( \@gidnoz, { gCopy => [0.0,C0hi] }); # 0.0,0.50 include gClass =~ zero here? or no
$stitle= join" ","# Moments for gCopy => [", 0.0,C0hi,"] genes";
$nsum += sumgene_sumout( $gidsel, \@FLDSUM, $outsumh, $stitle); # call several times, one output file

# if(@genextop){
#   print $outsumh "\n# Uniq Gene Coverage summary of $genexcopy \n"; 
#   map{ print $outsumh $_,"\n" } @genextop;
# }

close( $outsumh);
warn "#sumgene_sumout: n=$nsum $outsum\n" if $debug;

if(DO_RPLOT and $outrplot) { # run Rscript at end
  unless( bad_exec("Rscript")) {
  my $err= system("./$outrplot"); # exec here 
  warn "#gnodes_sumcovplot Rscript $outrplot, err= $err\n" if $debug;
  }
}
#=======================================================

## UNUSED input parts
# my($nsamp,$masmid,%samvals); # masmid == asmid nor not?
# ($nsamp,$masmid)= readMetad($sampledata);
#
# my($nann,%annvals);
# ($nann)= ($anntable) ? readAnnots($anntable) : (0); # format from gnodes_annotate.pl, same crID, crLoc, an, ids as covtab
# 
# my($nidclass,$idclassh,$idclasslist)= ($crclassf) ? read_idclass($crclassf, 1) : (0); 
#
# my($ntcds,$cdslen,$cdsrdmap)= ($cdscovtab) ? readChrtab($cdscovtab, 1) : 0;  # drop, only for cdslen in genecov{id}{Glen}
# 
# my($nccds)= ($cdscovtab) ? readCovtab($cdscovtab, 1) : 0; #? useful w/ genexcopy ?
#-------------------------------------------

sub readCovtab { # == readCovtab
  my($incdstab)=@_;
  my($ok,$inh)= openRead($incdstab);
  return unless($ok);
  my(%cdscov);
  my($nin,$lcr,$lcb)=(0,0,0);
  while(<$inh>){
  
    if(/^\W/) {
      # if(/^#ChrID/i){ }
      next;
    } 
    # cols for covtab: ChrID, Pos: covt covm covu aCovT aCovM aCovU 
    
    my @v=split; 
    my($id,$pos,$ct,$cm,$cu, $act,$acm,$acu)=@v; 
    if($act >= $MINRD) { # MINRD
    $cdscov{$id}{$pos}{'aCovT'}= $act;
    $cdscov{$id}{$pos}{'aCovM'}= $acm;
    $cdscov{$id}{$pos}{'aCovU'}= $acu;
    $nin++;  $lcr=$id; $lcb=$pos;
    }
  } close($inh);
  warn "#readCovtab($incdstab) n=$nin\n" if($debug);
  return($nin, \%cdscov);
}

sub readGenetab {
  my($ingenetab,$outname)=@_;
  # ingenetab == output of sam2covtab: putGenetab8e() == GeneID x Chr x CPos, aCovT/aCovM/aMiss

  my($keyLC1,$keyLL1,$keyLV1)= ('lCopy2','lLoci2','lCov2') ; #UPD21OCT : drop old loci_calc1 set?
  my($keyLC2,$keyLL2,$keyLV2)= ('lCopy','lLoci','lCov');
  
  my @IC=(3,4,5); # cov value cols
  # my $gloc8i=0; # #GeneID	GPos	aCovT	aCovM	aCovU	aCovZ from putGenetab8i
  if($Genetab8i){ @IC=(2,3,4,5); }
  our($nv,$nzero,$hdo,$outh,@sv,@mv,@rdv,@suv,@hdr)=(0,0);
  our(%cloc,%lcov,%lcovgb,%lcdsgb,%lspan,@crcov,@lcovall); # Genetab8i cr: loci cov

  my $outf= $outname . '.sumgenetab'; #???
  my($lid,$nout)=(0,0);

  use constant CR_MINCOVp =>0.50; # 0.25; # tried 0.75,0.50; >0.25 is WRONG for human gene set ?? ***
      # cr_mincov vs lcdsgb cutoff pct from top cdscov, was 0.25 - 0.50
  use constant CR_MINCOV2 => 1; # cr_mincov vs lcdsgb cov filter spurious dup loci
  use constant LCOPYMED => 0; # want lcopy av,med,span tuple or just ave?
  
  sub lloci_calc2 { 
    my($gid)=@_;
    our($nv, %cloc,%lcov,%lcovgb,%lcdsgb,%lspan,@crcov,@lcovall); # Genetab8i cr: loci cov
    my ($lmax,$lsum)=(0) x 9;
    # also: push @crcov, [$gb,$crcb,$covt,$covm,$covu,$crval]; # replace lcovgb,lcdsgb w/ this one list?

    my @gb= sort{$a <=> $b} keys %lcovgb;
    my $ngb= @gb; 
    my(@loc2list,@mloc,%gcloc,%glcov,@covm);
    my(%gbcrcb_ok); #upd.21oct24 return valid crcb/gb 
    my($scovm,$ncovm)=(0,0);
    for my $gb (@gb) {
    
      #FIXME: lcovgb{gb} has 2 cov vals for cr: $crval = cds-chr cov and covm = chr-total cov,
      # want both, crval to sort best locus crcb, covm as truest crcb cov val
      # add lcdsgb{gb}{crcb} = crval == genecdscov
      my @crcb= sort{ $lcdsgb{$gb}{$b} <=> $lcdsgb{$gb}{$a} or $a cmp $b } keys %{$lcovgb{$gb}};

      # CR_MINCOVp, use cdscov < 0.25*lcovcds test for spurious crcb
      # dont need here a test adjoining crcb, as dup loci measured at same gb bins
      
      my(%glcovgb); my $lcovcds=0;
      for my $crcb (@crcb) { # join near bins as per calc1, but for same gbin
        my $crcbOk= 1; my($cc,$cr,$cb)=split":",$crcb;
        my $covm=$lcovgb{$gb}{$crcb};
        
        # see below cr_mincov: try this way, not that
        # PROBLEM: covcds == 0 for some w/ valid, but low cover, due to parts cover = int(0.49) == 0 problem
        # ?? add here cr:cb = 1+lastcb overlap-bin filter, same as cb/G2BIN test?
        if(CR_MINCOV2) {
          my $covcds= $lcdsgb{$gb}{$crcb} ||0;
          if($covcds > $lcovcds){ $lcovcds= $covcds; } 
          elsif($covcds < CR_MINCOVp * $lcovcds) { $crcbOk=0; next; } # was 0.25, skip tiny gene-cov crcb
          #^^^ this filter bad for human dupx genes, for bad reasons getting low covcds vals for 2nd true copies,
          #     but have high/normal chrcov vals
        }
        
        # $gbcrcb_ok{$gb}{$crcb}=$crcbOk; # $covm;
        
        $scovm += $covm; $ncovm++; push @covm,$covm;    
        #ORIG: use constant G2BIN => 1000; # fixme
        use constant G2BIN => 100; #  TEST dplx
        my $cbb=int($cb/G2BIN); 
        my $ckey="$cr:".($cbb+0); my $cloc=0; 
        for my $i (0,1,-1,2,-2) { my $ckt="$cr:".($cbb+$i); last if($cloc=$cloc{$ckt}); } 
        unless($cloc){ $cloc="$cr:".($cbb*G2BIN);  } # not here:  $crcbOk=0;
        $crcbOk=0 if($glcovgb{$cloc}); # have overlap bin ??
        $gcloc{$ckey}=$cloc; $glcovgb{$cloc}+=$covm; # insert $gb here?
        $gbcrcb_ok{$gb}{$crcb}=$crcbOk; 
      }
      my @loc2listgb= sort{  $glcovgb{$b} <=>  $glcovgb{$a} } keys %glcovgb;
      my $nl= @loc2listgb;  $lsum+=$nl; push @mloc, $nl; $lmax= $nl if($nl>$lmax);
      map{ $glcov{$_} += $glcovgb{$_} } @loc2listgb; #??
      # ?? guesss loci from @crcb, joining across @gb
    }
    
    # lcopy2 DID over-estimates copy (> all others) for tab8i, due likely to overlapped gbins, tab8j may fix
    # lcopy2 this way is maybe best match to gcopy, better than lcopy, maybe better than acopy
    my $lcopy2= ($ngb<1)? 0 : $lsum / $ngb; #? maybe, or $lmax? or calc1 way? maybe test all for agreement?
    if(LCOPYMED) {
      @mloc= sort{ $b <=> $a } @mloc; 
      my $lcopymed= $mloc[ int($ngb/2) ]; #UPD add median here? doesnt seem informative == int (lcopy2)
      my $lcspan= ($nv<1)?0:int(0.5 + 100*$ngb/$nv); #? not useful
      $genecovh->{$gid}->{$keyLC2} = sprintf "%.1f,%d,%d",$lcopy2,$lcopymed,$lcspan;  
    } else {
      $genecovh->{$gid}->{$keyLC2} = sprintf "%.1f",$lcopy2; 
    }
    
    @loc2list= sort{  $glcov{$b} <=>  $glcov{$a} } keys %glcov;
    #? use this lLoci1 format?  push @loclist, "$ac,$mdcov,$pspan,$l"; 
    if(@loc2list > $lcopy2){ @loc2list=splice(@loc2list,0,int(0.5+$lcopy2)); }
    $genecovh->{$gid}->{$keyLL2} = join("; ",@loc2list);
    
    @covm= sort{$b<=>$a} @covm; my $covmed= $covm[ int($ncovm/2) ];
    my $covav= ($ncovm<1)?0: $scovm/$ncovm; 
    my $pspant=0; # fixme
       $pspant = ($nv<1 or $lcopy2 < 0.1 )? 0 : 100 * ($lsum / $lcopy2) / $nv; # nv = total gene bins, lsum/lcopy2 = ave gene bins over n-loci ~= ngb
    $genecovh->{$gid}->{$keyLV2} = sprintf"%.1f,%d,%.0f%", $covav, $covmed, $pspant;
    return($ncovm,\%gbcrcb_ok);
  }

  sub lloci_calc1 { # #UPD21SEP08 ** Output cols add lCopy, lLoc1,2,3,..9
    my($gid)=@_;
    our(%cloc,%lcov,%lcovgb,%lcdsgb,%lspan,@crcov,@lcovall); # Genetab8i cr: loci cov
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
  
  # sub putv { if(STORECOV) { storev(); } else { printv(); } }
  sub storev { 
    our($nv,$nzero,$hdo,$outh,@sv,@mv,@rdv,@suv,@hdr);
    our(%cloc,%lcov,%lcovgb,%lcdsgb,%lspan,@crcov,@lcovall); # Genetab8i cr: loci cov
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
    splice(@sv,1,2,"ave",$nv); 
    $no++;

    $genecovh->{$gid}->{'aCovT'} = $ct;
    $genecovh->{$gid}->{'aCovM'} = $cm;
    $genecovh->{$gid}->{'aCovZ'} = $cno;
    $genecovh->{$gid}->{'aCopy'} = $pca;
    $genecovh->{$gid}->{'aBins'} = $nv; #UPD21AUG20, output: aLen= aBins * $CBIN
    # ^^ nv bin count is over-count for aLen, many bins contain partial (ie few %) align to 100b, eg intron, end edges
    # .. any way to correct? using ct/cm deviation from median?
    # eg: Dapsim1EVm000991t1 uniq 1copy, gLen=3510 aLen=5600, medn count=34, 7 bins have 1..9 count, 8 have 10..19
    #  41 have ~full read count, ie 4100 bp  .. some extras =~ alt exons
     
    #UPD21OCT : maybe save [covt,covm,covu?,gb] for alt calc of aCovT,M .. sam2cov8j vals are not equal to UCG cov vals
    # .. but cr: row vals should be
    #  push @crcov, [$gb,$crcb,$covt,$covm,$covu,$crval]; # replace lcovgb,lcdsgb w/ this one list?
    # Here: add as  $genecovh->{$gid}->{'cCovT'} cCovM, cCopy? and median cmCopy? where cCovM should be UCG cov for uniqs
    # Not quite accurate enough, esp. cxCopy = cCopy * cCovM/CU
    # .. try filter by $crval >= 0.50 * $topcrval, or similar to weed out spurious cds-chr aligns

    #UPD21OCT24: use gbcrcb_ok to filter @crcov
    my($ncovm,$gbcrcb_ok)= (0,{});
    if($Genetab8i) { #UPD21SEP08 ** Output cols add lCopy, lLoc1,2,3,..9
      lloci_calc1($gid);      
      ($ncovm,$gbcrcb_ok)= lloci_calc2($gid); # LLOCI_CALC => 2 test, sets lLoci2, FIXED: set to lCopy[1]; overest copy bias, due to overlaps?
     }

=item UPD21OCT24 HAS_CDSCOV

    # UPD21OCT24: maybe CDSFILT should be more restrictive: crval >= 0.75? * topval? CR_MINCOVp
    # uniqs are getting spurious dupl due to smallish cdscov, true dups have near equal cdscov ?
    # still have problem when chr map is more uniq than cds map, so covt/covm is 1.0 for cds-dups

    HAS_CDSCOV bug here from DBG21OCT_covtsort with extreme covt outliers picked 1st
    .. maybe should use median not ave for ctc,cmc,cuc ? 
    ** YES single extreme outlier can throw off ctgb => ctc_nu, but need sums for ctgb < cdst test
    .. ie dromel Mito 100x cov val for 1 gene bin will mess this up
    .. use some median<>ave calc to exclude extremes in chr.covt, 
    #? sort gbin>covt crcov, exclude extremes (>> and <<) befor this ctgb += covt ?
    #? use aCovT/M as guide for what is extreme?
    .. Or use median covt,m over gbin as guide for extemes, ie assume few gbins have outliers
=cut

    use constant CCOV_CDSFILT => 1;
    if(UPD21OCT and @crcov>0) {
      my ($cutval,$nkeep)=(0,0);
      my($ct,$cm,$cu,$ncr)=(0) x 4;
      my($ctc,$cmc,$cuc,$ncrc,$cdst)=(0) x 5; # HAS_CDSCOV vals
      my %ctgb=(); my %gbcrv=(); 
      
      use constant DBG21OCT_covtsort => 0;
      use constant DBG21OCT_covtmedian => 1;
      
      if(DBG21OCT_covtmedian) { # remove extreme outliers in covt
        # this looks messy but test: works for drosmel outliers, retest daphplx, aweed, human genesums 
        # .. may need to drop some gbins w/ many extremes to get median ctc_nu
        @crcov= sort{ $$a[0] <=> $$b[0] or $$b[5] <=> $$a[5] or $$b[3] <=> $$a[3] } @crcov; # new: gbin>crval>covm
        my($lgb, @gb, @gbcm, @gbct, %gboutl, @cm)=(0);
        my $crend= [-1,0,0,0,0,0];
        for my $crc (@crcov, $crend) {
          my($gb,$crcb,$covt,$covm,$covu,$crval)= @$crc; 
          if($gb ne $lgb and @cm){
            push @gb, $lgb;
            my $n=@cm; my $nh=int($n/2);
            @cm= sort{$b <=> $a} @cm;
            if($n>2) {
            my($cmt,$cmm,$cmz)= @cm[0,$nh,-1];
            push @gbct, $cmt; push @gbcm, $cmm; 
            } else {
            my($cmt,$cmm)= @cm[0,-1];
            push @gbct, $cmt; push @gbcm, $cmm; 
            }
            @cm=();
          }
          push @cm,$covm;
          $lgb= $gb;
        } 
        my $ng= @gbct; my $nh= int($ng/2);
        my @gbcs = sort{$b <=> $a} @gbct;
        my $gbctm= $gbcs[$nh];
        my $outm= 4 * $gbctm;
        for my $i (0..$ng) {
          # if($gbct[$i] > $outm and $gbcm[$i] < $outm)  # ?? and $gbcm[$i] < $outm : drop this one?
          if($gbct[$i] > $outm)  # ?
          {
            my $gb= $gb[$i]; $gboutl{$gb}=1; # drop outliers from @crcov: change $gbcrcb_ok->{$gb}{$crcb}
          }
        }
        my $nout=0; my $ncrcov=@crcov; # count, cant drop all/most
        for my $crc (@crcov) {
          my($gb,$crcb,$covt,$covm,$covu,$crval)= @$crc; 
          if($gboutl{$gb} and $covm > $outm) { $nout++; }          
        }
        if($nout < 0.66 * $ncrcov) {
        for my $crc (@crcov) {
          my($gb,$crcb,$covt,$covm,$covu,$crval)= @$crc; 
          if($gboutl{$gb} and $covm > $outm) {
            $gbcrcb_ok->{$gb}{$crcb}= 0;
          }
        }
        }
                      
      } elsif(DBG21OCT_covtsort) { # TEST dplx HAS_CDSCOV; BAD for extreme outliers in covt
        @crcov= sort{ $$a[0] <=> $$b[0]  or  $$b[2] <=> $$a[2] } @crcov; # gbin>covt
      } else { # CCOV_CDSFILT
        #old: @crcov= sort{ $$b[5] <=> $$a[5] or $$b[3] <=> $$a[3] or $$a[0] <=> $$b[0] } @crcov; # old: crval>covm>gbin
        @crcov= sort{ $$a[0] <=> $$b[0] or $$b[5] <=> $$a[5] or $$b[3] <=> $$a[3] } @crcov; # new: gbin>crval>covm
      }
      
      if(CCOV_CDSFILT) {
        #O: @crcov= sort{ $$b[5] <=> $$a[5] or $$b[3] <=> $$a[3] or $$a[0] <=> $$b[0] } @crcov; # crval>covt>gbin
        my $topval= $crcov[0]->[5]; 
        $cutval= CR_MINCOVp * $topval; # << maybe this should depend on which gb
        my $nc= @crcov; $nkeep=  2; # was ($nc < 5) ? $nc: int($nc/2);
      }
      
      for my $crc (@crcov) {
        my($gb,$crcb,$covt,$covm,$covu,$crval)= @$crc;
        if($ncovm) {
          next unless($gbcrcb_ok->{$gb}{$crcb}); #UPD21OCT24, skip spurious crcb/gb
        }
        if($HAS_CDSCOV) {
          # HAS_CDSCOV bug here from DBG21OCT_covtsort with extreme covt outliers picked 1st
          
          $cdst= $cdscovh->{$gid}{$gb}{'aCovT'}||0; # this is *trick* here, sum chr-covt up to CDS-cov-total
          if($ctgb{$gb} < $cdst){ $ctgb{$gb}+=$covt; $ctc+=$covt; $cmc+=$covm; $cuc+=$covu; $ncrc++; }
        }
        if(CCOV_CDSFILT) { if($gbcrv{$gb}){ next if( $crval < CR_MINCOVp * $gbcrv{$gb}); } else { $gbcrv{$gb}=$crval; } }
        #O: if(CCOV_CDSFILT) { next if($ncr >= $nkeep and $crval < $cutval); } #? last if() for sorted @crcov; keep 1/2 of @crcov ?
        $ct+=$covt; $cm+=$covm; $cu+=$covu; $ncr++;
      }

      $pca= ($cm<1)?0:sprintf"%.1f",$ct/$cm; 
      if($ncr>1){ ($ct,$cm,$cu)= map{ sprintf"%.0f", $_/$ncr } ($ct,$cm,$cu); }
      $genecovh->{$gid}->{'cCovT'} = $ct;
      $genecovh->{$gid}->{'cCovM'} = $cm;
      $genecovh->{$gid}->{'cCopy'} = $pca;
      if($HAS_CDSCOV) {
        my $ngb= scalar(keys %ctgb); # nrcr < ngb ??
        # FIXME: eCovT:ctc is not div by ncrc, but added up to cds covt, while cmc/ncrc is right here !???
        my $ctc_nu=0; my @ctgb=(); # try median vs av;
        for my $gb (sort keys %ctgb){ my $ctgb=$ctgb{$gb}; $ctc_nu += $ctgb; push @ctgb,$ctgb; } 
        if($ncrc < 1) { ($ctc,$cmc,$cuc,$ncrc)=($ct,$cm,$cu,$ncr); } # fixup for missing
        if($ncrc > 1) { ($ctc,$cmc,$cuc)= map{ sprintf"%.0f", $_/$ncrc } ($ctc,$cmc,$cuc); }
        $ctc_nu= ($ngb < 1)? 0 : sprintf"%.0f", $ctc_nu/$ngb;
        if(1) { @ctgb= sort{ $b<=>$a } @ctgb;  $ctc_nu=$ctgb[int($ngb/2)]; } #debug test median
        $pca= ($cmc<1)?0:sprintf"%.1f",$ctc_nu/$cmc; 
        $genecovh->{$gid}->{'eCovT'} = $ctc_nu;
        $genecovh->{$gid}->{'eCovM'} = $cmc;
        $genecovh->{$gid}->{'eCopy'} = $pca;
      
      }
    }
    
#     my($ncovm,$gbcrcb_ok)= (0,{});
#     if($Genetab8i) { #UPD21SEP08 ** Output cols add lCopy, lLoc1,2,3,..9
#       lloci_calc1($gid);      
#       ($ncovm,$gbcrcb_ok)= lloci_calc2($gid); # LLOCI_CALC => 2 test, sets lLoci2, FIXED: set to lCopy[1]; overest copy bias, due to overlaps?
#      }
     
    if(@mv > 0) { #? median always, or option?
      my @sm=sort{$$b[0] <=> $$a[0]} @mv; 
      my $hn=int($nv/2); my $hm=$sm[$hn]; 
      my @mvo=@sv; 
      @mvo[@IC]=@$hm;  # splice(@mvo,3,3,@$hm);       
      if(@IC>3){  ($mt,$mm,$mu,$mno)= @$hm[@IC]; } else { ($mt,$mm,$mno)= @$hm[@IC]; }
      
      my $pcm= ($mm<1)?$mt:sprintf"%.1f",$mt/$mm;  
      my $pnom=$pno; #not for median: ($mm<1 or $mno == 0)?0: sprintf "%.2f%%",100 * $mno/$mm; 
      splice(@mvo,1,2,"med",$nv); 
      $genecovh->{$gid}->{'mCovT'} = $mt;
      $genecovh->{$gid}->{'mCovM'} = $mm;
      $genecovh->{$gid}->{'mCopy'} = $pcm;
    }
    
    if(@rdv > 0) { # readgene counts
      if(@IC>3){  ($mt,$mm,$mu,$mno)= @rdv[@IC]; } else { ($mt,$mm,$mno)= @rdv[@IC]; }
      my $pcm= ($mm<1)?$mt:sprintf"%.1f",$mt/$mm;  
      splice(@rdv,1,1,"reads"); #,$nv

      $genecovh->{$gid}->{'aReadT'} = $mt;# readgene CovT read count
      $genecovh->{$gid}->{'aReadM'} = $mm;# readgene CovM read count
      $genecovh->{$gid}->{'aReadU'} = $mu;
      $genecovh->{$gid}->{'aReadZ'} = $mno; #? is it valid
      $genecovh->{$gid}->{'aCopyRd'} = $pcm;
    }
    
    # if(@suv > 0 and not @rdv) { # sumgene counts, skip if @rdv
    #   if(@IC>3){  ($mt,$mm,$mu,$mno)= @suv[@IC]; } else { ($mt,$mm,$mno)= @suv[@IC]; }
    #   my $pcm= ($mm<1)?$mt:sprintf"%.1f",$mt/$mm;  
    #   splice(@suv,1,1,"sumc"); #,$nv
    #   $genecovh->{$gid}->{'aReadT'} = $mt;
    #   $genecovh->{$gid}->{'aReadM'} = $mm;
    #   $genecovh->{$gid}->{'aReadU'} = $mu;
    #   $genecovh->{$gid}->{'aReadZ'} = $mno; #? is it valid
    #   $genecovh->{$gid}->{'aCopyRd'} = $pcm;
    # }
    return $no;
  }
  
  
  my($ok,$inh)= openRead($ingenetab);
  $outf="";  #  no output STORECOV, genecovh then output table
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
    # cols for GeneID, ChrID, Pos: CovT CovM zero = read map counts/bin as per chr.covtab
    # cols for GeneID, sumgene,  1: allCovT allCovM noCov 
    # cols for GeneID, readgene, 1: nrdmaps, nreads, nomap
    
    my @v=split; 
    my($id,$pt,$nb,$ct,$cm,$cz)=@v; 
    if($lid ne $id){
      $nout+= storev() if($nv); 
      $nv=0; @suv=@rdv=@mv= @sv=();  ($lcr,$lcb)=(0,0);
      @sv=@v; for my $i (@IC){ $sv[$i]=0; } 
      %cloc=%lcov=%lcovgb=%lcdsgb=%lspan=(); @crcov= @lcovall=();
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
      
      if(CR_MINCOV2 or $crval >= $cr_mincov) { # was covm; need to filter out tiny cov cases, what mincov ? use 0.20?*$sv[1]/$nv as min-expected
      my($cc,$cr,$cb)=split":",$crcb;
      
      $lcovgb{$gb}{$crcb}=$covm; # for LLOCI_CALC==2, save till storev($gid) ?
      $lcdsgb{$gb}{$crcb}=$crval||1; # $covm; # crval is best gene-locus cov, while covm is best chr cov, use both
 
      #UPD21OCT : maybe save [covt,covm,covu?,gb] for alt calc of aCovT,M .. sam2cov8j vals are not equal to UCG cov vals
      # .. but cr: row vals should be
      push @crcov, [$gb,$crcb,$covt,$covm,$covu,$crval]; # replace lcovgb,lcdsgb w/ this one list?
      
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
      

    } else { # pt == ChrID
      my($cr,$cb)= (0,$pt); # ($Genetab8i) ?  : ($pt,$nb); # proper var names

      if($ct < $MINRD)  #? should be ct < 1 ?
      { 
        $nzero++; # ** ?? need to record zeros here ?? this blocks aCovZ ?
      } else { 
        $nv++; 
        push @mv, [@v]; # was @v[@IC] 
        for my $i (@IC){ $sv[$i]+=$v[$i]; } 
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
  
  $nout+= storev(); 
  warn "#readGenetab: n=$nout, cbin=$CBIN, $outf\n" if $debug;

  return($nout, $outf);  
}


sub eqcopy { 
  my($gcopy,$ccar,$cclab)= @_; 
  my $ncc=@$ccar; my @ev;
  return "na" unless($gcopy >= 0.5 and $ncc > 0); # which cutoff?
  for(my $i=0; $i<$ncc; $i++) {
    my $icopy= $ccar->[$i]; my $cl= $cclab->[$i] || (1+$i);
    my $xcopy= $icopy/$gcopy;  
    my $ceq= ($xcopy < $XClo)? 'lo' : ($xcopy > $XChi)?'hi': 'eq';  # same as for sum_couts
    push @ev,"$cl:$ceq";
  }
  return join",",@ev;
}

sub sumgene_counts {
  my($ctabh, $title, $gidall, $gidnoz, $gidzero, )= @_;

  # FIXME: gCopy cuts for counts should ~ match gCopy cuts for moments
  use constant CTOOLOW => 0.1;
  my %gidnoz= map{ $_ => 1 } @$gidnoz;
  my @gidzero= grep{ not $gidnoz{$_} } @$gidall; # recalc, gidnoz now option w/ gidclean subset
  $gidzero= \@gidzero;
  my $zerodef= ($USE_GIDCLEAN) ? "gClass=zero|skew or gNread<10 or gLen<200" : "gClass=zero or gNread<$MINRD";
  
  my $nall= @$gidall; my $noz= @$gidnoz; my $nzero= @$gidzero;
  print $ctabh $title,"\n";
  print $ctabh "# N Genes valid=$noz, poor-class=$nzero, all=$nall (poor:$zerodef)\n";
  
  # putcn() and putcn_missing()
  my($ncnsum, $ncmiss, %cnsum, %cnbases, %gclass,%eqToGcopy)= (0,0);
  my($nrrow, %rrow)= (0); # DO_RPLOT
  
  for my $id (@$gidall) {
    my($glen,$gnread,$gcopy,$gcm,$gclass)= map{ $genecovh->{$id}->{$_}||0 } qw( gLen gNread gCopy gCM gClass); 
    $gclass{$gclass}{'gidall'}++; 
  }
  
  #UPD21SEP09: add lCopy, drop aLen for counts
  my $xcopyis= (USE_lCopy_xCopy == 1) ? "xCopy=lCopy/gCopy, ": "xCopy=aCopy/gCopy, ";
  if(SWAP_AC_COPY){ 
    $xcopyis.= "aCopy.upd from Chrasm covtab, "; # SWAP_AC_COPY
    if(my $ack= $ENV{acopyis}) { $xcopyis =~ s/upd/upd=${ack}Copy/; } # ie, cCopy => eCopy
    }
  # UPD21OCT: mCopy = median for acopy
  # USE_mCopy4aCopy means swap mCopy for aCopy for xCopy num?
  for my $id (@$gidnoz) {
    my($glen,$gnread,$gcopy,$gcm,$gclass,$acopy,$mcopy,$acovt,$acopyrd,$areadt,$areadm,$areadz,$alen,
       $lcopy,$lcov,$lxcopy, $ccopy, $axcopy, $ecopy,$excopy)= 
      map{ $genecovh->{$id}->{$_}||0 } 
      qw( gLen gNread gCopy gCM gClass aCopy mCopy aCovT aCopyRd aReadT aReadM aReadZ aLen lCopy lCov lxCopy 
          cCopy axCopy eCopy exCopy);  
      # ^ FIXME: need to set above new cols?
      
    if($acopy =~ m/,/){ ($acopy)=split",",$acopy; } #UPD21OCT tacked on mcopy for tabout
    my($lcopya,$lcopym,$lcopys)=($lcopy,0,0);
    if($lcopy =~ m/,/) { #UPD21OCT, or not
      ($lcopya,$lcopym,$lcopys)=split",",$lcopy;
      $lcopy= $lcopya; # $lcopym; #median now ?? or not, median lcopy is not really diff from ave, but for no fraction
    }
    
    my($xcopy,$xcopy2);
    # USE_lCopy_xCopy == 2 , revert to xcopy = acopy/gcopy as more reliable : 
    #  21OCT20: lcopy best for human,pig .. getting ~uniq chr dna map counts (covm == covt) for dupx cds genes,
    #   ie. diff multimap results for chr vs cds, so aCovT/aCovM == 1 where lCopy = 2..9
    #  However, daphplx doesnt have this result, but aCovT/aCovM more reliable
    # .. maybe due to mix of gDNA, CDS, chrasm seq sets, diff map rates for diff strains?  
    # .. effects of lower multi-map qual of one cds dupx vs uniq chr map?
    # ?? account for w/ cds.bam dup rate vs chr.bam dup rate? or cds.covtab x chr.covtab > genetab vals?
    if(USE_lCopy_xCopy == 1) { # also add lxCopy
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
    my $ceq= ($xcopy < $XClo)? "xlo" : ($xcopy > $XChi)?"xhi": "xeq";  
    if($areadt < $MINRD or $acopyrd < CTOOLOW) {  # MINRD = 1 default, want higher min here?
      $ncmiss++; $ceq="xzero";
    }
    
    ##UPD21OCT add counts/class of this extended xCopy counts, debug only?
    #   $genecovh->{$id}->{eqToGcopy} == eqcopy($gc,\@cpv,\@cplab); # eqToGcopy col ~~  a:lo,r:lo,l:lo,xl:eq,4:lo
    if(my $eqtogc= $genecovh->{$id}->{'eqToGcopy'}) {
      my @eq=split",",$eqtogc;
      for my $gcla ($gcl,'-1') { map{ my($k,$v)=split":",$_;  $eqToGcopy{$k}{$gcla}{$v}++ if($v); } @eq; }
    }
    
    #? add gcl="All" ?
    $ncnsum++;
    for my $gcla ($gcl,'-1') { # -1 == All for sort
      $cnsum{$gcla}{ngene}++; $cnsum{$gcla}{$ceq} ++; 
      $cnsum{$gcla}{gcopy}+= $gcopy; $cnsum{$gcla}{xcopy}+= $xcopy; 
      $cnsum{$gcla}{acopy}+= $acopy; $cnsum{$gcla}{mcopy}+= $mcopy; $cnsum{$gcla}{axcopy}+= $axcopy; 
      $cnsum{$gcla}{ccopy}+= $ccopy; ## add cCopy ..
      $cnsum{$gcla}{ecopy}+= $ecopy; $cnsum{$gcla}{excopy}+= $excopy; 
      $cnsum{$gcla}{lcopy}+= $lcopy; $cnsum{$gcla}{lxcopy}+= $lxcopy; 
      $cnsum{$gcla}{aglen}+= $acopylen;
      # add areadz count, if > minz
      $cnsum{$gcla}{missrd}+= $areadz;  $cnsum{$gcla}{areadm}+= $areadm; 
      $cnsum{$gcla}{miss5p}++ if( $areadz > 0.049*$areadm);
      #? $cnsum{$gcla}{acopyrd}+= $acopyrd; # .. areadt ?
      if(PLOT_GENEBASES ) { #UPD21OCT add rplot of gene copynum x genelen = tot bases per class, gCopy vs aCopy only?
        unless($gcl >= 100) { # dont sum huge copynum set, few genes, skews view
        $cnbases{$gcla}{ngene}++; $cnbases{$gcla}{glen} += $glen;
        $cnbases{$gcla}{gcopy} += $glen * $gcopy;
        $cnbases{$gcla}{acopy} += $glen * $acopy; # == ecopy of acopy.upd ??
        }
      }      
    }
  }
  
  # putcn_sum
  use constant AlwaysX0 => 1; # ncmiss only show of x0, is too confusing.
  my @eqToGkeys= sort keys %eqToGcopy;
  my @gcl=sort{ $a<=>$b } keys %cnsum;  
  my $ngt=", for nGene=$ncnsum, Chrasm_id=$asmid"; 
  my @ceql= qw(xlo xeq xhi); 
  unshift @ceql, "xzero" if(AlwaysX0 or $ncmiss>0); # can confuse to have some w/ some w/o xzero: always add?
  my $labxloeqhi= (AlwaysX0 or $ncmiss>0) ? "x0,xLo,xEq,xHi= xCopy 0,<1,=1,>1"
    : "xLo,xEq,xHi= xCopy<1,=1,>1";
  (my $labxleh= $labxloeqhi) =~ s/=.*//;
  #D if(UPD21AUG){ $labxloeqhi .= ", agLen= aLen/gLen"; }
  if(UPD21AUG){ $labxloeqhi .= ", lCopy=Chr loci"; } #UPD21SEP09
  if(UPD21OCT and @eqToGkeys){ $labxleh.="\teqGcopy:lo,eq,hi"; }
  # if(SWAP_AC_COPY){ $labxloeqhi .= ", aCopy of cr:covtab, cCopy of gene ave."; } #  xcopyis has this
  if(SWAP_AC_COPY){ 
    my $acl="cCopy"; if(my $ack= $ENV{acopyis}) { $acl= "${ack}Copy"; } # ie, cCopy => eCopy
    print $ctabh "# aCopy of gene-chr cov, $acl of gene-ave cov, lCopy of cds-chr overlap\n"; 
  } else { print $ctabh "# aCopy of gene-ave cov, cCopy of gene-chr cov, lCopy of cds-chr overlap\n"; }

  my $cnbrief="Genes Copynum Summary$ngt, xCopy= ";
  print $ctabh " ---------------------------------------\n";
  ## UPD21OCT: moved to table footnote, dupl of # title
  #d print $ctabh " Genes Copynum Summary by Copy level$ngt\n"; #== title
  print $ctabh " ".join("\t",qw(CLevel nGene gCopy xCopy aCopy lCopy cCopy aMiss), $labxleh)."\n";

  my $ngall=$cnsum{'-1'}{ngene} || $noz; # == $noz gidnoz, must be >0
  for my $gcl (@gcl) {
    my $ng= $cnsum{$gcl}{ngene} or next;  # ng includes xzero for COPYNUMzeros, x/t reduced by zeros
    my($xc,$tc, $mcopy, $gc, $misp, $misr, $lcopy, $lxcopy, $aglen, $ccopy, $axcopy, $ecopy, $excopy)= 
      map{ sprintf"%.1f",$cnsum{$gcl}{$_}/$ng } 
      qw( xcopy acopy mcopy gcopy miss5p missrd lcopy lxcopy aglen ccopy axcopy ecopy excopy) ;
    my($ceqn)= join",", map{ $cnsum{$gcl}{$_}||0 } @ceql; # qw(xlo xeq xhi); #? add xzero

    #UPD21OCT need keys for each iCopy here: sort keys %eqToGcopy ?
    my $ceqxt=""; 
    for my $ek (@eqToGkeys) {
      my $cev= join",", map{ $eqToGcopy{$ek}{$gcl}{$_}||0 } qw(lo eq hi); #? no zero
      $ceqxt.="$ek:$cev;";
    }
    
    my($stot,$smis)= map{ $cnsum{$gcl}{$_}||0 } qw(areadm missrd);
    my $prmis= ($stot<1)? 0: sprintf"%.2g",100*$smis/$stot; 
    my $mis="$prmis%,$smis";
    
    my $gclo= ($gcl eq '-1')?"All" : ($gcl eq "100-499")? "99-499" : $gcl;

    my $acopyo= $tc; # (UPD21OCT)? $tc : "$tc,$mcopy";
    my $lcopyo= $lcopy; #>> drop this lx: maybe add tc==acopy,ax:$axcopy ?
    if(UPD21OCT and USE_lCopy_xCopy > 1) {
      if($axcopy > 0.01 and $axcopy ne $tc and abs($gc - $axcopy) < abs($gc - $tc) ) { $acopyo= "$tc,ax:$axcopy"; } 
    } elsif(USE_lCopy_xCopy) {
      if($lxcopy > 0.01 and $lxcopy ne $lcopy and  abs($gc - $lxcopy) < abs($gc - $lcopy) ) { $lcopyo= "$lcopy,lx:$lxcopy"; }
    }
    
    print $ctabh " ".join("\t",$gclo,$ng,$gc,$xc,$acopyo,$lcopyo,$ccopy,$mis,$ceqn,$ceqxt)."\n"; 
    $cnbrief .= $gclo."c:$xc/$ng, " unless($gcl eq "0" or $gcl > 100);

    if(DO_RPLOT) { #UPD21OCT, from gnodes_sumasmgcn2rtab.pl
      # FIXME: DOMISR Cmiss bad for rplot see dmag19sk too low, change: use Call gcl == -1 sum
      if($DOMISR > 1 and $gcl eq '-1') {
        my $prhit = sprintf"%.2f", (100 - $prmis);
        my $prmisf= sprintf"%.2f", $prmis;  $prmisf= -$prmisf if($DONEG);
        $rrow{'Cmiss'} = [$prmisf,$prhit,$smis,$stot, ];
      }
      unless($gcl <= 0 or $gcl >= 100) { # save 1,2-9,10-99 .. skip All,
        my $cgn= ($gcl eq '-1')? "CAll" : ($gcl<10 and $gcl>=0) ? "C0$gcl" : "C$gcl";   

        my($pxlo,$pxeq,$pxhi)=map{ my $c=$cnsum{$gcl}{$_}||0;  sprintf"%.2f",100*$c/$ngall; } qw(xlo xeq xhi); 
        my $prhit = sprintf"%.2f", (100 - $prmis);
        my $prmisf= sprintf"%.2f", $prmis;   

        if($DONEG){ $prmisf= -$prmisf; $pxlo= -$pxlo; } 
        $rrow{$cgn} = [$pxlo,$pxeq,$pxhi,$prmisf, ]; $nrrow++;
        if($DOMISR == 1) { my $cgmis= $cgn."mis"; $rrow{$cgmis} = [$prmisf,$prhit,$smis,$stot, ]; }

        #?? insert col in rrow for gene Mbases, only % aCopy - gCopy ? also both gcopy,acopy Mb?
        #?? or add as new rrow set: $bgn= B$gcl ?
        if(PLOT_GENEBASES ) { #UPD21OCT add rplot of gene copynum x genelen = tot Mbases per class
          # $cnbases{$gcla}{ngene}++; $cnbases{$gcla}{glen} += $glen;
          # $cnbases{$gcla}{gcopy} += $glen * $gcopy;
          # $cnbases{$gcla}{acopy} += $glen * $acopy; # == ecopy of acopy.upd ??
          my $mbot= $cnbases{'-1'}{gcopy}; #? is this right
          my $mbo= $cnbases{$gcl}{gcopy};
          my $mbe= $cnbases{$gcl}{acopy};
          my $peo= ($mbot<1)? 0 : int( 1000 * ($mbe - $mbo)/$mbot )/10;
          ($mbe,$mbo,$mbot)= map{ sprintf "%.1f",$_ / 1_000_000 } ($mbe,$mbo,$mbot); # change to megabases
          # fixme: keep Cmiss, diff gscale=4 ??? or replace w/ Bmiss in %bases ??
          (my $bgn=$cgn) =~ s/^C/B/;
          $rrow{$bgn} = [$peo,$mbe,$mbo,$mbot, ]; $nrrow++;
          #NOT HERE if($DOMISR) { my $bgmis= $bgn."mis"; $rrow{$bgmis} = [$prmisf,$prhit,$smis,$stot, ]; }
        }      

      }
    }
    
  }
  print $ctabh " ---------------------------------------\n";
  ## UPD21OCT: moved to table footnote, dupl of # title
  #dup# print $ctabh " Genes Copynum Summary by Copy level$ngt\n"; #== title
  print $ctabh " aCopy=Chrasm, gCopy=Gene set, $xcopyis $labxloeqhi\n";
  print $ctabh " aMiss=Missed gene reads on chrasm, % of gene reads\n";
  
  print $ctabh "\n Gene Classes of $genexcopy\n";
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

if(DO_RPLOT and $nrrow) {
  # do_rplot(\%rrow); or return \%rrow as plot data
}
  
  # warn "#$cnbrief table to copytab\n" if($debug);
  return($ncnsum, $cnbrief, $nrrow, \%rrow);
}

=item sumgene_rplot

try1:
SppAsm Item pXlo pXeq pXhi pMiss
at18test8jn Aall -32.39 68.03 119 176
at18test8jn Auniq -13.64 81.97 108 132
at18test8jn Adup -18.75 25.77 11.4 44.4
at18test8jn Acds -13.07 68.97 51 74
at18test8jn C-1 -2 96.75 1.21 -0.01  << C-1 == CAll, skip
at18test8jn C01 0 0.00 0.00 0
at18test8jn C02-9 0 0.00 0.00 0
at18test8jn C10-99 -0.07 0.14 0.00 -0.02
at18test8jn Cmiss -0.02 99.98 593 2542383

# ($ret)= sumgene_rplot($outh, $asmname, \%rowcrsum, \%rowgsum);
=cut

sub sumgene_rplot {
  my($outf, $asmname, $rowc, $rowg)= @_;
  
  # my($outpdf=$outf) =~ s/\.\w+$/.pdf/; 
  # but rscript : plname=(plname,spn,".pdf",sep="") and spn=asmname == asmgcn_at18test8jn.pdf
  # should this be plname=(asmname,plname,.pdf) == at18test8jn_asmgcn_plot.pdf
  # my $fsa=$asmname; # == 1st row in R datatable, used as label, pdf name
  
  my($nRplot,$hdr)=(0,0);
  my %row=(); 
  map{ $row{$_} = $rowc->{$_} } keys %$rowc;
  map{ $row{$_} = $rowg->{$_} } keys %$rowg;

  my @cl= sort keys %row;
  my @hdr=qw(SppAsm Item pXlo pXeq pXhi pMiss );
  my($ook,$outh)= openOut($outf); return unless($ook);
  
  if($DO_AWEIGHT){
    if(my $allmb = $row{'Aall'}->[3] ) {
      my @arv= map{ $row{$_} } grep /^A/, @cl;
      for my $arv (@arv) { my($pxlo,$pxeq,$omb,$emb)= @$arv; 
        my $plomb= sprintf"%.2f",100*($emb - $omb) / $allmb; $plomb= -$plomb if($DONEG);
        $arv->[0]= $plomb;
      }
    }
  }
  
  if($DOMISR) {  # sum(nrmis)/sum(nrtot), not all cgmis = C1mis,C2mis,..
    # $row{$cgmis} = [$prmis,$prhit,$nrmis,$nrtot, ]; 
    # FIXME: bad CMiss, replace w/
    #  $rrow{"Cmiss"} = [$prmisf,$prhit,$smis,$stot, ];
    if($row{'Cmiss'}) {
      #ok: my($prmis,$prhit,$nrmis,$nrtot, )= $row{'Cmiss'};
      if(PLOT_GENEBASES) { $row{'Bmiss'}= $row{'Cmiss'}; push @cl,"Bmiss"; } # or change R script?
    } else {
      my @arv= map{ $row{$_} } grep { m/^C[0-9].+mis$/ } @cl; # ** SKIP Call,C0 for sum, not valid addition
      my($smis,$stot)=(0,0);
      for my $arv (@arv) { my($prmis,$prhit,$nrmis,$nrtot, )= @$arv; $smis+=$nrmis; $stot+=$nrtot; }
      my $prmis= ($stot<1)?0: sprintf"%.2f", 100*$smis/$stot; my $prhit= 100 - $prmis;
      if($DONEG){ $prmis= -$prmis; }
      $row{'Cmiss'}= [$prmis,$prhit,$smis,$stot, ];
      my @cla= grep { not m/^C.+mis$/ } @cl;
      push @cla,"Cmiss";
      if(PLOT_GENEBASES) { $row{'Bmiss'}= [$prmis,$prhit,$smis,$stot, ]; push @cla,"Bmiss"; } # or change R script?
      @cl=@cla;
    }
  }
  
  if($DO_GWEIGHT){
  
  }
  
  #x if($DROPCLASS) { @cl= grep{ not m/$DROPCLASS/ } @cl; }
  # @cl reorder: Aall  Acds  Adup  Auniq > Aall Auniq Adup Acds
  map{ s/Auniq/Ab/; s/Adup/Ac/; s/Acds/Ae/; } @cl; # s/Ate/At/; s/Arpt/Ar/; 
  @cl= sort @cl;
  map{ s/Ab/Auniq/; s/Ac/Adup/; s/Ae/Acds/; } @cl;

  my $DO_RSCRIPT=DO_RPLOT; #? always
  if($DO_RSCRIPT) {  
    my $rval = join(" ",@hdr)."\n";
    for my $cl (@cl) { my $rv= $row{$cl}; $rval .= join(" ",$asmname,$cl,@$rv)."\n"; }
    rplot_crgenesum($outh, $rval, 1 > $nRplot++);
    
  } else {
    print $outh join("\t",@hdr)."\n" if(1>$hdr++); 
    for my $cl (@cl) { my $rv= $row{$cl}; print $outh join("\t",$asmname,$cl,@$rv)."\n"; }
  }
  
  close($outh); system("chmod +x $outf") if($DO_RSCRIPT);
  return(1);
}

sub rplot_crgenesum { # UPD for PLOT_GENEBASES, 21oct31
  my($outh,  $datatable, $first)=@_;

  # FIXME rplot opt, should calc from data
  my $YMIN=$ENV{ymin}||-40; #  was -60 for dapmag min ?? maybe calc from all input data range?
  $YMIN=-$YMIN if($YMIN>0);
  
# plname? gnodes_ chr_gene_deficits_ ?? asmgcn_miss_ ?? test: asmgcxlo3d_
  if($first) {
  
    #PLOT_GENEBASES: BnotC = flag to switch to deficit/excess of %Bases
    # plotsize: 6,4 < should be nrow, max(4, abs(YMIN)/10) ?
    my $bbs='\\\\';
    my $rtop=<<"EOS";
#!/usr/bin/env Rscript

acdeficitplot <- function(vd, plname="asmgcn_plot",ymin=$YMIN, gscale=4, BnotC=F) {
  if(BnotC) {
    rcn= grep("^C",vd[,2]); if(length(rcn)>0){ vd= vd[ -rcn,]; }
    pnames= vd[,2]; rcn= grep("^B",pnames);  cmiss="Bmiss"; # gscale=2; # fixme
  } else {
    rcn= grep("^B",vd[,2]); if(length(rcn)>0){ vd= vd[ -rcn,]; }
    pnames= vd[,2]; rcn= grep("^C",pnames);  cmiss="Cmiss";
  }
  vd[rcn,3]= gscale*vd[rcn,3]; 
  vcol=rep("gray",nrow(vd)); vcol[rcn]="green"; 
  rm= grep(cmiss,pnames); if(length(rm)>0){ vcol[rm]="red2"; }
  pnames= gsub("^A(${bbs}w+)","${bbs}U${bbs}1",pnames,perl=T)
  spn=as.character(vd[1,1]); # title or asmid
  plname=paste(spn,"_",plname,".pdf",sep="");
  plw= 0.80 * max(8, nrow(vd)); plh= 0.10 * max(50,10+abs(ymin));
  pdf(plname,plw,plh);  # 6,5 inch min
  barplot( vd[,3], ylab="% Deficit / Excess", ylim=c(ymin,10), col=vcol,
     main=paste("Chrasm/GeneCN for",spn), names.arg=pnames, cex.names=0.85, space=0.55); 
  gticks=seq(10,ymin,-10); axis(4, at=gticks,  labels=gticks/gscale);
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

# sub rplot_crgenesum_OLD {
#   my($outh,  $datatable, $first)=@_;
# 
#   # FIXME rplot opt, should calc from data
#   my $YMIN=$ENV{ymin}||-40; #  was -60 for dapmag min ?? maybe calc from all input data range?
#   $YMIN=-$YMIN if($YMIN>0);
#   
# # plname? gnodes_ chr_gene_deficits_ ?? asmgcn_miss_ ?? test: asmgcxlo3d_
#   if($first) {
#   
#     # plotsize: 6,4 < should be nrow, max(4, abs(YMIN)/10) ?
#     my $bbs='\\\\';
#     my $rtop=<<"EOS";
# #!/usr/bin/env Rscript
# 
# acdeficitplot <- function(vd, plname="asmgcn_plot",ymin=$YMIN, gscale=4) {
#   pnames= vd[,2]; rcn= grep("^C",pnames); 
#   vd[rcn,3]= gscale*vd[rcn,3]; 
#   vcol=rep("gray",nrow(vd)); vcol[rcn]="green"; 
#   rm= grep("^Cmiss",pnames); if(length(rm)>0){ vcol[rm]="red2"; }
#   pnames= gsub("^A(${bbs}w+)","${bbs}U${bbs}1",pnames,perl=T)
#   spn=as.character(vd[1,1]); # title or asmid
#   plname=paste(spn,"_",plname,".pdf",sep="");
#   plw= 0.80 * max(7, nrow(vd)); plh= 0.10 * max(40,abs(ymin));
#   pdf(plname,plw,plh);  # 6,4 inch min
#   barplot( vd[,3], ylab="% Deficit", ylim=c(ymin,0), col=vcol,
#      main=paste("Chrasm/GeneCN for",spn), names.arg=pnames, cex.names=0.85, space=0.55); 
#   gticks=seq(0,ymin,-10); axis(4, at=gticks,  labels=gticks/gscale);
#   mtext("Chr Asm", side=1, line=0, adj=0.25); mtext("Gene CN", side=1,line=0, adj=0.75); 
#   dev.off()
# }
# 
# EOS
#     print $outh $rtop;
#   }
#   return unless($datatable); # can print top before datatable
# 
#   my $rplot=<<"EOS";
# vd=read.table(header=T, text=\"
# $datatable\");
# acdeficitplot(vd);
#   
# EOS
#     print $outh $rplot;
# }


sub sumgene_tabout {
  my($gids, $fields, $ofn, )= @_;

  my($no)=(0); 
  my @fld= @$fields;
  my($ook,$outh)= openOut($ofn); return unless($ook);

  if(SWAP_AC_COPY) {
  my $acl="cCopy"; if(my $ack= $ENV{acopyis}) { $acl= "${ack}Copy"; } # ie, cCopy => eCopy
  print $outh "#acopyis: aCopy.col=$acl, swapped with $acl.col=orig aCopy \n" unless($acl eq "cCopy"); # FIXME nasty message
  }
  print $outh join("\t","GeneID",@fld)."\n";
  for my $gid (@$gids){ 
    my @val= map{ $genecovh->{$gid}->{$_} || 0 } @fld;
    print $outh join("\t",$gid,@val)."\n"; $no++;
  } close($outh);
  return($no,$ofn);  
}

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

sub sumgene_sumout {
  my($gids, $fields, $outh, $title, $corto)= @_;

  my($no)=(0); 
  my @fld= @$fields;
  $outh= *STDOUT unless(defined($outh));
  
  my $coradd= ($corto and ref($corto));
  my $corhd= ($coradd) ? "\tCorrelations" : "";
  
  print $outh "$title\n" if($title);
  print $outh "StatOpts=$STOPTS, Transform=$TRANS, \n" if($STOPTS or $TRANS);
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


sub _max{ return($_[0] < $_[1])?$_[1]:$_[0]; }
sub _min{ return($_[0] > $_[1])?$_[1]:$_[0]; }
sub _minnot0{ return ($_[0] == 0) ? $_[1] : ($_[1] == 0) ? $_[0] : _min(@_); }

sub skew_p { my($ave,$med,$sd)=@_; my $skew= ($sd<0.001)? 0 : 3 * ($ave - $med) / $sd; return $skew; }

sub format_v { local $_= $_[0];
  my $d=($_ >= 100)?0:($_ >= 10)?1:($_ >= 1)?2:($_ >= 0.1)?3:4; $_= sprintf"%.${d}f",$_; $_ =~ s/\.0+$//; return $_; 
} 
    

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
    #BADbug: $var=($var - $ave*$ave)/($nit-1); $sd=sqrt($var); $se=$sd/sqrt($nit); 
    $var=($var - $nit*$ave*$ave)/($nit-1); $sd=sqrt($var); $se=$sd/sqrt($nit);  #UPD2202
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
    #BB $covar= ($covar - $ave*$aveb)/($nit - 1 ); #?? check
    #BADbug $var=($var - $ave*$ave)/($nit-1); $sd=sqrt($var); # $se=$sd/sqrt($nit); 
    #BADbug $varb=($varb - $aveb*$aveb)/($nitb-1); $sdb=sqrt($varb); 
    $covar= ($covar - $nit*$ave*$aveb)/($nit - 1 ); #?? check
    $var=($var - $nit*$ave*$ave)/($nit-1); $sd=sqrt($var); # UPD2202
    $varb=($varb - $nitb*$aveb*$aveb)/($nitb-1); $sdb=sqrt($varb); 
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


sub read_crsumtext {
  my($intext)= @_;
  
  my @asum= qw(allasm uniqasm dupasm cdsasm); # opt?
  # allasm uniqasm dupasm cdsasm CDSann CDSbus TEann RPTann NOann
  if($ENV{aann}){ # is this default?  CDSann better than cdsasm
    @asum= qw(allasm uniqasm dupasm CDSann TEann RPTann); # opt?
  }
  
  my $asumpatt=join"|",@asum;
  my($nsum,$datend)=(0) x 9;
  my(%rrow);
  
  my($ok,$inh)= openRead($intext); 
  # return($nsum,\%rrow) unless $ok;
  while(<$inh>){
    
    if(/^\s*Source=/ and /KUlow=/) { # chrasm_sum.txt top  
      # Source=Dropse20uc, KUlow=88.7, KUhigh=94, FlowcytSize=161-180 Mb Formula_LN/C=166.8-176.7 Mb (dropse20chrs)
      if($datend>0){ $datend=0; }# fnsource($ARGV); 
      
    } elsif(/^($asumpatt)/)  { #? allasm uniqasm dupasm CDSann TEann RPTann
      next if($datend);
      # print "#asm:$_" if($DO_SHOWD);
      my @v=split"\t"; 
      my($casm,$omb)= @v;
      my($emb,$xcopy)= (@v<7) ? @v[2,3] : @v[4,5];
      
      #o $casm=~s/asm//; $casm="A$casm";
      $casm=~s/asm|ann//; $casm="A".lc($casm);
      
      my $pxeq= ($xcopy<0.01)? 0 : sprintf"%.2f", 100/$xcopy;    
      my $pxlo=  sprintf"%.2f",100 - $pxeq; $pxlo= -$pxlo if($DONEG);
      $rrow{$casm}=[$pxlo, $pxeq, $omb, $emb]; #?? ,
      $nsum++;    
    } elsif(/^\s*Total span=/) {
      # Total span=163.2 Mb, covspan=163.2, gaps=0 Mb for assembly dropse20chrs
    
    } elsif(/^\s*Genome Size Est=/) {
      # Genome Size Est= 168.4-178.5 Mb (Nread), 166.8-176.7 Mb (Maprd), for readset SRR11460802ab,
      #  for Size=LN/C, Cov=88.7,94, N_reads=105560784, N_maprd=104503221,99.0%, L_readlen=150
    
    } elsif(/^\s*Uniq Conserved Gene Cover/) {
      # Uniq Conserved Gene Cover: median=88.7, ave=89.1, sem=3.62, n=612
      $datend=1;
    }
  
  } close($inh);
  
  warn "#read_crsumtext: n=$nsum of $intext\n" if $debug;
  return($nsum,\%rrow);
}

sub read_genexcopy {
  my($genexcopy)=@_;
  my($nucg,$cmd,$cav,$cerr,$ndup,$nuni,$ngc,$nzero,$topsum)= (0) x 9;
  my(%gcovtab);

  # use constant { gMINLEN => 50, gMINREAD => 3 }; #?? genexcopy drop these filts? or move to outputs?
  
  #UPD21OCT : parse for Lrdlen, in L*N/C or next
  # Genome_size= 3458.6 Mb of L*N/C= 126 * 526870326 / 19.2, CDS_size= 216.8 Mb,
  #   for Nr.ucgcds=33019050 (6.3%), Nr.allcds=33019050 (6.3%), Nr.total=526870326, Lrdlen=126, MapErr=237308701 (7.82%), Nshort=0
  # $Lrdlen=0; # global?
  
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
      } elsif($Lrdlen < 9 and /Genome_size=/ and m,L.N/C=\s*(\d+),) { $Lrdlen=$1; 
      } elsif($Lrdlen < 9 and m/Lrdlen=(\d+)/) { $Lrdlen=$1; }
      
      if($topsum){ 
        s/^#//; chomp(); push @genextop,$_;
      }
    
    } else { # if(/^\w/) # Gene_ID or gene row, collect all gene rows
      # cols: Gene_ID Glen Nread Rdlen tCopy C.M Uniq Merr C.nz S.Dup,Unq,Nt,Equ,UCG
      my @v=split; $topsum=0;
      my($gid,$glen,$nread,$rdlen,$tcopy,$cmed,$isuniq,$merr,$cnz,$sdune)=@v;
      next if(/^Gene_ID/);
      
      #UPD: special case gClass $isuniq eq "zero", store? or put aside? here or output?
      # my $GOK= 1; #was ($glen >= gMINLEN and $nread >= gMINREAD);
      if(1) { #$GOK  need to record nread == 0 cases ***
        # switch to named col vals, extend w/ other gene data
        $gcovtab{$gid}= { gCopy => $tcopy, gLen => $glen, gNread => $nread, gClass => $isuniq,
                         gCM => $cmed, gCnz => $cnz, };
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

sub bad_exec {
  my @exfail=(); # my %exfail;
  for my $e (@_) {
    my $err= system("$e --help >/dev/null 2>&1");
    push @exfail, $e if($err);
  }
  return @exfail;
}

__END__

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

=item UPD 8j genetab adds cr:crid:crb acovt,m,u cdscov; drop sumgene row

#GeneID	GPos	aCovT	aCovM	aCovU	aCovZ
AT1G01010t1	readgene	372	372	372	0
AT1G01010t1	100	cr:NC_003070.9:3900	44	44	44	8
AT1G01010t1	100	43	43	43	0
AT1G01010t1	100	cr:NC_003070.9:3800	32	32	32	32
AT1G01010t1	200	19	19	19	0
AT1G01010t1	200	cr:NC_003070.9:4100	24	24	24	3
AT1G01010t1	200	cr:NC_003070.9:4000	40	40	40	21

=cut

      
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

=item meanvar

 -- see evigene/scripts/genoasm/cds_meanvar.pl
 -- ** add covar/cor calc for 2+ vars
 -- maybe add transforms: log(x) ; exp(x); sqrt(x) ..
 -- see also cds_meanvar:sub trim()
 
=cut


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

