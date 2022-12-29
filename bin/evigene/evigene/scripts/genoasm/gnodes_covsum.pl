#!/usr/bin/perl
# gnodes_covsum.pl for evigene/scripts/genoasm/ 
# from daphnia3covfilt.pl, daphnia3covtables.pl

=item usage gnodes3_covsum

  gnodes_covsum.pl -covtab $pt.covtab [-out_covmean $pt.covmean] [-out_summary $pt.covsum.txt]
    opts: -genomedata|metadata|?? daphnia3species.genomes.data

=item gnodes_covsum trial

  input dmag15nwb2asm_SRR7825549cl_1a_bwa.cdschr7b.covtab from
    evigene/scripts/genoasm/gnodes_sam2covtab.pl (now sam2covtab7b.pl)
    
  evigene/scripts/genoasm/gnodes_covsum.pl -debug -title dmag15nwb2asm_cov7b  \
   dmag15nwb2asm_SRR7825549cl_1a_bwa.cdschr7b.covtab
  
  # output to sum=dmag15nwb2asm_cov7b_sum.txt means=dmag15nwb2asm_cov7b.means
  #title allasm/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 
  #covtab span=126.1.mb, covspan=107.3.mb info: n_readid=4450506, n_mapok=79875595, n_mapbad=130272508, n_dupbad=236433975, 
  #meanvar nt=1073743, ti=allasm/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 
  #meanvar nt=857529, ti=uniqasm/dmag15nwb2asm_cov7b uniqrd=1, cdsrd=0, an=0, terd=0, unkrd=0 
   ...
  
  ==> dmag15nwb2asm_cov7b_sum.txt <==
  Source=dmag15nwb2asm_cov7b, KUlow=45, KUhigh=48.85, FlowcytSize=0 Formula_LN/C=0 (dmag15nwb2asm_cov7b)
  _____   ______  Low     Low     High    High    Total    
  Item_   Obs.Mb  Est.Mb  xCopy   Est.Mb  xCopy   tCopy   Description
  ---------------------------------------------------------------------------
  allasm  107     161     1.50    175     1.63    4.00    measured assembly
  uniqasm 86      86      1.00    93      1.09    1.09    asm with unique gDNA
  dupasm  21      75      3.51    81      3.81    15.56   asm with multimap gDNA
  cdsasm  46      62      1.36    68      1.47    2.24    asm with CDS-mapped gDNA
  teasm   8.7     21.5    2.49    23.4    2.70    6.07    asm with TE-mapped gDNA
  unkasm  10.7    57.5    5.39    62.5    5.85    27.39   asm with UnclassRepeats-mapped gDNA
  NOann   107     161     1.50    175     1.63    4.00    asm without annotations
  CDSbus  0.0     0.0     0.00    0.0     0.00    0.00    asm with unique BUSCO orlog CDS
  ---------------------------------------------------------------------------
    xCopy = excess read copy depth; tCopy = assembly multi-map copy depth (-xCopy)
  
  NOTE: NOann == allasm in this trial w/o annotations
  Rename 'UnclassRepeats' > ?? UnclassRepmodl-mapped  UnknownRM? Unclass/Unknown Repeats/Duplicates?


=item input covtab from sam2covtab7b.pl
  -- same table of chrid, cpos, covcds_3 covall_3
  .. counting dna read hits per 100 bp
  .. but covcds_3 cols change meaning: CovT/CovM/CovU now are, for chrasm, reads matched to CDS/TE/UNK seqs
  
  #cdsxchr_covtab par: minident=0.65, mindupid=0.98, softclip=0, BIN=100, part=1/8 
  #dmag19skasm_SRR7825549cl_1a_bwa.cdschr7a.covtab.pt1 n_readid 40419628,   
  #  n_mapok 7997661, n_mapbad 4321056, n_dupbad 16887553, 
  #ChrID	Pos	CovT	CovM	CovU	aCovT	aCovM	aCovU
  dapmag19sk:LG2	100	0	0	0	1	1	1
  dapmag19sk:LG2	200	0	0	0	4	4	4
  dapmag19sk:LG2	300	0	0	0	9	9	9
  dapmag19sk:LG2	400	0	0	0	12	12	12
  
  -- also have chrtab, sums per chrid (or geneid..), useful or not?
  -- total reads, maplocs are found at end as '# total_XXX' lines, per part
  dmag19skasm_SRR7825549cl_1a_bwa.cdschr7a.chrtab.pt1
  # nmap,nmult,nuniq = locations, nmread,nuniq,nomap = read counts
  # n_notcr = 0
  #ChrID        	chrlen  	nmap   	nuniq  	nmult  	nmread	nomap
  dapmag19sk:LG2	16359456	3126127	479851	2646276	650830	0
  dapmag19sk:LG1	14067088	1322635	382735	939900	455601	0
  dapmag19sk:LG3	11088946	1167344	298312	869032	353424	0
  
=item real data usage

  aweed20gnodes/
  pt=arath18tair_chr; 
  pt=arath20max_chr; 
  env kucg=33 $evigene/scripts/genoasm/gnodes_covsum2p.pl -plotchr -debug \
   -genexcopy arath18tair1cds_SRR10178325_b2_mim.geneycopy -asmid $pt  -title${pt}_testplota -anntab $pt.anntab  \
   -crclass arath18tair1cds.idclass -sumdata arath20asm.metad   ${pt}_SRR10178325_test8f.covtab
  
  env kucg=0 $evigene/scripts/genoasm/gnodes_covsum2p.pl -plotchr -debug  \
   -genexcopy arath18tair1cds_SRR3703081_b2_mim.geneycopy  -asmid $pt  -title ${pt}_testplhetz -anntab $pt.anntab  \
   -crclass arath18tair1cds.idclass -sumdata arath20asm.metad   ${pt}_SRR3703081_hetest8f.covtab
  
  dromel20gnodes/
  pt=drosmel6ref_chr;
  pt=drosmel20pi_chr;  
  env kucg=0 $evigene/scripts/genoasm/gnodes_covsum2p.pl -plotchr -debug \
    -genexcopy dromel6relt1cds_SRR11460802_b2_mim.genexcopy -asmid $pt  -title ${pt}_testplota \
   -anntab $pt.anntab  -crclass dromel6relt1cds.idclass -sumdata drosmelchr.metad ${pt}_SRR11460802_b2_mim.cc8a.covtab

  env kucg=0 $evigene/scripts/genoasm/gnodes_covsum2p.pl -plotchr -debug \
    -genexcopy dromel6relt1cds_SRR10512945_1odd.genexcopy   -asmid $pt  -title ${pt}_testplhetz\
   -anntab $pt.anntab  -crclass dromel6relt1cds.idclass -sumdata drosmelchr.metad ${pt}_SRR10512945_1odd_bwa.cc8a.covtab
      
=item usage prior daphnia3covfilt.pl, daphnia3covtables.pl 
  
  mv  all8cov3h_sum.tables  all8cov3h_sum.tables.old
  touch  all8cov3h_sum.tables
  
  for ctab in *cov3h.tab; do { 
    cnam=`basename $ctab -generdtwo_cov3h.tab`;
    onam=$cnam.cov3h.meansum; touch $onam; echo "out $onam";
    for cfilt in allasm uniqasm dupasm cdsasm CDSann CDSbus CDSuni CDSdup TEann RPTann NOann; do { 
      env filt=$cfilt name=$cnam ../daphnia3covfilt.pl $ctab >> $onam
    } done
    
    ../daphnia3covtables.pl ../daphnia3covsum.data $onam >> all8cov3h_sum.tables
  } done

=cut

use FindBin;
use lib ("$FindBin::Bin/../"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
our $EVIGENES="$FindBin::Bin/..";  

use strict;
use Getopt::Long;  

use constant UPD7f => 1; #UPD 21apr30
use constant UPD21JUN => 1; # update annots, add sam2genecov UCG CU/KU val to pickKUni as new "best" CU
use constant UPD21AUG => 1; # CONTAMann or other filter of contaminant parts/scaffolds/...; 
#UPD21AUG  add plot tables, large 50k-100k bins/chr, same parts as _sum.txt (all,uniq/dup/cds asm, CDS/TE/RPT/NO ann)
#          and Rplot function to draw cov depth xCopy graphs over chr spans, all & select parts
#   sam2genecov table, via gnodes_pipe8a == xxx.genexcopy #?? what suffix
use constant UPD22FEB => 1; # pickKUni test update, LN/C est - CONTAM

my $debug=$ENV{debug}||1;
my @IVAR=(5,6,7); # covtab cols ==  rd_total, rd/mmap, rd_uniq for chrasm x rd map 
my $DOMED= 0; # $ENV{median}||0; * SUCKS UP ALL covtab vals, mem pig
my $DOHIST=1; # replace median from @allvals w/ histo mode(s)/peak, est median, counts for depth vals depth[0..999]?
my $ZEROS= 0; # $ENV{zero}||0;
my $showSUM=1; # $ENV{sum}||0;
my $NTAVE= $ENV{ntave} || 0;
my $MINRD= 5; #? option: skip all cov rows w/ acovt < MINRD, want to measure zeros? gaps?
my $MINDUNIQ= 2; # was 1, for test of uniqspan == covt <= covu + MINDUNIQ
my $CBIN=100; # calc from rows ib dist
use constant kSAMPLE => 1000;  # CBIN size median sample @cbin, should be all same = default $CBIN
my $PLOTCHR=0; #UPD21AUG
my $PLOTBIN= $ENV{plotbin} || 50_000; #? 100_000; #UPD21AUG option? #
my $PLOT_COVT= $ENV{plotcovt} || 0; # PLOT_COVT USE_COVT
my $gnodes_rplots="$EVIGENES/genoasm/gnodes_chrplot.R";

my $KUCG= $ENV{kucg}||0; #UPD21AUG

# UPD21JUN modify UCG_KU_SPECIAL, add/replace w/ genecov KU vals if available
use constant UCG_KU_SPECIAL => 1; # special KU cov-depth calc for uniq-conserved-genes: median of median depth
  # seems to add reliability, valid non-par statistic (med-of-med) w/ error est, 
  # theory: median cover (act? or acm) per UCgene is best depth est per gene, 
  #   then median over all UCG medians (~1000) is best est for entire genome
use constant UCG_KU_MEANS => 'UCGmed'; # filt_sum()/meanvar() for UCG also, @parts == UCGmed ?
use constant UCG_FIXID => 1; # use %ucg_id_cov  : UPD21JUN change this

my($nucg,@ucg_med_cov,@ucg_med_ids,@ucg_one_cov,%ucg_id_cov)=(0); # UCG_KU_SPECIAL globals

#  #UPD21AUG UPD21aug06: add here? filter out annot == contam for bacter.contam scaffolds .. use anntab? or separate table?
#  # maybe want some stats on dropset? use filter parts w/ anntab contam flags?

my %filts= (
  allasm =>  "urd=0 crd=0 an=0",
  uniqasm => "urd=1 crd=0 an=0",
  dupasm =>  "urd=2 crd=0 an=0",
  cdsasm =>  "urd=0 crd=1 an=0",  # cmp CDSm
  #UPD add teasm, unkasm parts ? may want mix of crd=1 unkrd=1 or/and trd=1 unkrd=1
  teasm =>  "urd=0 crd=0 terd=1 an=0",   
  unkasm =>  "urd=0 crd=2 terd=2 unkrd=1 an=0", # unclassified from rptmodlr, may be CDS or TE or other;
    # only measure when crd/terd == 0, ie =2 flag here
  #o: unkasm =>  "urd=0 crd=0 unkrd=1 an=0", # special case of unclassified from repetmodlr, may be CDS or TE or other
  CDSann =>  "urd=0 crd=0 an=CDS",
    CDSbus =>  "urd=0 crd=0 an=busco",
    CDSuni =>  "urd=0 crd=0 an=CDSuni",
    CDSdup =>  "urd=0 crd=0 an=CDSdup",
  TEann =>  "urd=0 crd=0 an=TE",
  RPTann => "urd=0 crd=0 an=repeat",
  GAPann => "urd=0 crd=0 an=gap", # added, check .. may need to adjust many small gaps vs big gaps
  # leave out GAPann in @DEFparts for now .. most gaps == 0 read cov, already excluded
  NOann =>  "urd=0 crd=0 an=no(CDS|TE|repeat|gap)",
  CONTAM => "an=contam", #UPD21AUG: if use this, need others including allasm to skip an=~/contam/
);

my %filtdesc= (
  allasm =>  "measured assembly",
  uniqasm => "asm with unique gDNA",
  dupasm =>  "asm with multimap gDNA",
  cdsasm =>  "asm with CDS-mapped gDNA",   
  teasm  =>  "asm with TE-mapped gDNA",   
  unkasm  =>  "asm with UnclassRepeats-mapped gDNA",   
  CDSann =>  "asm with CDS annotations",
    CDSuni =>  "asm with unique ortholog CDS", #?? CDSuni ~= CDSbus
    CDSbus =>  "asm with unique BUSCO orlog CDS",
    CDSdup =>  "asm with paralog CDS",
  TEann =>  "asm with Transposons",
  RPTann => "asm with simple Repeats",
  GAPann =>  "asm with gaps",
  NOann =>  "asm without annotations",
  CONTAM =>  "asm contamination", #UPD21AUG
);

my @DEFparts= qw(allasm uniqasm dupasm cdsasm teasm unkasm CDSann CDSbus TEann RPTann NOann ); #  GAPann CDSuni CDSdup   cdsunan 
push @DEFparts, 'CONTAM' if(UPD21AUG);
# GAP annots need more tests, but MINRD effectively blocks all/most-gap bins;
# UNK is new an but ignore for now, seq annot equiv to unkasm


my @parts= @DEFparts; # opts change?
my $RTABLE=$ENV{rtab}||0;  my $didRhd=0; # opt R format sum.txt
my $keeps='rdacovt|rdacovm|aCovT|aCovM'; #? rdacovu ? rdcovt,m,u??

my $intitle= $ENV{title}||$ENV{name}||""; # was "asource"; 
my $asmid= $intitle; #? intitle for output, asmid for data
my $filtn= $ENV{filt}||"";
my ($sampledata,$anntable,$outmeans,$outsum,$ncpu,$crclassf,$readsetid,$genexcopy)=("") x 9;

my $optok= GetOptions( 
  #? 'covtab=s',\$covtab, 'chrtab=s',\$chrtab,
  'title|name=s',\$intitle,
  'asmid=s', \$asmid, # superceeds intitle for finding/using asm info
  'readsetid=s', \$readsetid, # upd21f20
  'filter=s',\$filtn, #<< replace w/ @parts loop over filters
  'metadata|sumdata|sampledata=s', \$sampledata, # optional input  
  'anntable|annotations=s', \$anntable,
  'genexcopy|genecovtab=s', \$genexcopy, # UPD21JUN sam2genecov covtab == $genecdsname.genexcopy; 
  'crclassf|idclassf=s',\$crclassf, # table of chr/gene id => class ? want for BUSCO, TEfam, : No, rely on anntable having these
  'output|outsummary=s', \$outsum,
  'meanvar|outmeans=s', \$outmeans,
  'minread=i', \$MINRD, # 5 default, test, see if it affects over-assembly (Xcopy < 1)
  'minduniq=i', \$MINDUNIQ, # 2 now, ?boost for crappy longreads
  'plotchr!',\$PLOTCHR, #UPD21AUG
  'KUCG=s',\$KUCG, # opt, mainly for plotchr? for pickKuni also? supplant KUlo,KUhi
  'Rformat!', \$RTABLE, 
  'median!', \$DOMED, # -nomedian to save mem  >> change to depth_histogram, count hits of D[0..999],Dk[(1000..10000)/10] ? Dx[>10k]
  'histogram|peak|mode!', \$DOHIST, # -nomedian to save mem  >> change to depth_histogram, count hits of D[0..999],Dk[(1000..10000)/10] ? Dx[>10k]
  'ncpu=i', \$ncpu, # not used here
  'debug!', \$debug, 
  );

die "usage: gnodes_covsum.pl chrasm_readset.covtab : table from gnodes_sam2covtab 
  opts: -title=label : output prefix
 -genex genecds.genexcopy : gene copy table from gnodes_genescov
 -ann asm.anntab : assembly annotations table from gnodes_annotate
 -asmid chrid -metad species.metad : assembly id and metadata about assemblies
 -plotchr : output chr-span coverage table for plots, instead of assembly coverage summary
  ... other opts ..
\n" unless($optok);

=item FIXME upd 21feb20 
  
  ** Need readset id in title and/or outsum,outmeans; 
     assume not provided in title,outsum,outmeans
  .. covtab inputs are from @ARGV or STDIN, perl <>
  .. cur covtab names are "asmid_readset_rdmap.vers.covtab"
  dmag19skasm_SRR7825548_1_bwa.cc7d.covtab  
  dmag19skasm_SRR7825549cl_1_bwa.cc7d.covtab
  .. cur gnodes_pipe call is
  $EVIGENES/genoasm/gnodes_covsum.pl  -asmid dmag19skasm -title dmag19skr48a -anntab dmag19skasm.anntab \
     -crclass dmag7fincds.idclass -sumdata dmag20asm.metad -debug dmag19skasm_SRR7825548_1_bwa.cc7d.covtab
  $EVIGENES/genoasm/gnodes_covsum.pl  -asmid dmag7fincds -title dmag7fincds .. dmag7fincds_SRR7825548_1_bwa.cc7d.covtab
   
=cut


if($asmid and not $intitle) { $intitle=$asmid; } # want one of these
elsif($intitle and not $asmid){ $asmid=$intitle; }

unless($readsetid) {
  my $rdset="";
  for my $incov (@ARGV) {
    if($incov =~ m/(\w+).(\w+)\.covtab/) { 
      my($asmrd,$vers)=($1,$2); 
      $asmrd =~ s/($asmid|$intitle)//g;
      if($asmrd=~/([A-Za-z0-9]+)/){ my $na=$1;  $rdset .= $na unless($rdset=~m/$na/); }
    }
  }
  if($rdset =~ m/\w\w/){ $rdset=substr($rdset,0,19) if(length($rdset>19)); $readsetid=$rdset; }
  #? else { $readsetid="na"; }
}


my $title= $intitle;
# [-out_covmean $pt.covmean] [-out_summary $pt.covsum.txt]
unless($outsum) { 
  ($outsum=$intitle) =~ s/\W/_/g; # title always ok?
  $outsum.="_".$readsetid if($readsetid and $outsum !~ m/$readsetid/);
  $outsum.="_sum.txt";
}
unless($outmeans){ ($outmeans=$outsum) =~ s/.sum.*//; $outmeans.=".means"; }

if($PLOTCHR){ 
  ($outmeans=$outsum) =~ s/.sum.*//; $outmeans .= ".plottab"; #??
  $outsum =~ s/.sum.*//;  $outsum .= ".plotchr.sh"; # Rscript
  @parts= grep{ not m/CDSbus/ } @parts; # FIXME: option
  warn "#covsum plotchr to plottab=$outmeans, Rscript=$outsum\n" if $debug;
} else {
  warn "#covsum output to sum=$outsum means=$outmeans\n" if $debug;
}

my($nsamp,%samvals)= (0);
($nsamp)= readMetad($sampledata) if($sampledata);
#? $samvals{$asmid}{readset}= $readsetid if($readsetid);

my($nann,$ANhasGTOP,$ANhasCONTAM,%annvals,%anntypes)= (0,0,0);
if($anntable) {
($nann)= readAnnots($anntable) ; # format from gnodes_annotate.pl, same crID, crLoc, an, ids as covtab
$ANhasGTOP=1 if($anntypes{'Gtop'});
$ANhasCONTAM=1 if($anntypes{'contam'});
}

# my($ngenex,$genexh)= ($genexcopy)?readGeneCopy($genexcopy):(0);  # UPD21JUN only for gene C.M/KUgenes val?
my($gucg_n, $gucg_med, $gucg_ave, $gucg_err)=(0) x 4; #now global for plotchr, pickKuni
if($genexcopy) {
  ($gucg_n, $gucg_med, $gucg_ave, $gucg_err)= read_genexcopy($genexcopy); 
##   my $CUniq= $KUCG || $gucg_med; $CUniq=1 if($CUniq<1);
  my $KUO=($KUCG)?"KUCG=$KUCG supercedes ":"";
  warn "#pickKuni Genes$KUO gucg_med=$gucg_med, sd=$gucg_err, ave=$gucg_ave, n=$gucg_n\n" if $debug;
}

# added for special case CDS.covtab > sum, annots in idclass tab
my($nidclass,$idclassh,$idclasslist)= ($crclassf) ? read_idclass($crclassf, 1) : (0); # cds,te class by id, $idclass->{id} = class

my($URD,$CRD,$TERD,$UNKRD,$FAN)=(0) x 9; # filter flag globals
my @FAN=(); 
setFiltVals($filtn||"allasm"); # this sets title, want

sub setFiltVals {
  my($filtn)=@_; # == @DEFparts
  # UCG_KU_MEANS special here ?
  if(UCG_KU_MEANS and $filtn eq UCG_KU_MEANS) {
    $title="$filtn/$intitle"; ($URD,$CRD,$TERD,$UNKRD,$FAN)=(0) x 9;
    $title.=" uniqrd=1, cdsrd=1"; return;
  }
  if($filtn and my $fv=$filts{$filtn}) {
    $title="$filtn/$intitle";
    ($URD,$CRD,$TERD,$UNKRD,$FAN)=(0) x 9; #* zero these globals
    for my $fk (split" ",$fv) {
      my($k,$v)=split"=",$fk; 
      $URD=$v if($k eq "urd");   $CRD=$v if($k eq "crd"); 
      $TERD=$v if($k eq "terd"); $UNKRD=$v if($k eq "unkrd");
      $FAN=$v if($k eq "an");
    }
    if($FAN and $FAN=~/\w\w/) { @FAN= grep( /\w/, split(/[,\s]/,$FAN)); } else { @FAN=(); }
    $title.=" uniqrd=$URD, cdsrd=$CRD, an=$FAN, terd=$TERD, unkrd=$UNKRD ";  
  }
}


my($nt,@sum,@nit,@svar,@sval,@dhist); # meanvar sum globals
my(%fnt,%fsum,%fnit,%fsvar,%fsval,%fdhist); # globals for filt_sum
my(%sum,@lvar); #meansToSum global, doesnt need to be

my($outhdr,$outmeanh,$outsumh)=(0);
warn "#title $title\n";  

use constant FLOOP => 1;
# FIXMEd: input all data for each filter in @DEFparts, write filt meanvar to $outmeanh 
# for cfilt in @filtset; do { env filt=$cfilt name=$cnam ../daphnia3covfilt.pl $ctab  .. }
# either a. reopen/read input covtab several times, or b. hash save each filter's table vals = FLOOP
# for b. need these/filt: my($nt,@sum,@nit,@lvar,@svar,@sval); # meanvar sum globals

my($n_readid,$n_nomap,$readlen,$n_mapok,$n_mapbad,$n_dupbad, $ninrow,$nincov,$nbiggap)= (0) x 9; # from covtab comments
## n_mapok,mapbad,dupbad are map stats ie w/ dup maps; n_readid, n_nomap are read counts, no_map were rejected by mapper,
## .. unavail is n_readok <= n_mapok, = only for no dup map .. see chrtab sums for n_readok

=item covtab comments, parse?
  
  header
#ChrID	Pos	CovT	CovM	CovU	aCovT	aCovM	aCovU

  at end, v7b : n_readid is for part now, NOT total same for all parts
    useful, sum over parts: n_mapok n_mapbad n_dupbad
#pt0.dmag19skasm_SRR7825549cl_1a_bwa.cdschr7b.covtab.pt0 n_readid 5398297, n_nomap 330240, n_notcr 0, n_mapok 7641320, n_mapbad 10749912, n_dupbad 10970001, n_mismatch 18223374, n_insert 4274267, n_delete 3226332, n_softclip 0, n_intron 0
#pt1.dmag19skasm_SRR7825549cl_1a_bwa.cdschr7b.covtab.pt1 n_readid 5398297, n_nomap 328059, n_notcr 0, n_mapok 7644057, n_mapbad 10752982, n_dupbad 10971034, n_mismatch 18226650, n_insert 4277759, n_delete 3232077, n_softclip 0, n_intron 0

=item chrtab comments, parse?
  * also want some stats from chrtab *

dmag15nwb2asm_SRR7825549cl_1a_bwa.cdschr7b.chrtab
#ChrID	chrlen	nmap	nuniq	nmult	nmread	nomap
dmag24nwb7c_contig00412	2025	622	258	364	315	0
dmag24nwb7c_contig00564	409	70	32	38	47	0
...
  * total_reads are for part subset
  * cds_reads is total of cdste.readids ?, same all pts, use to measure chrasm missed reads
    .. should cds_reads total == mapt/in? >> mapt/in here is bad count, cdsrdclass_in, will become part count of cdsrd-hit+cdsrd-miss,
       so cdsrd-miss = in - hit, or change table to list cdsrd-miss
    
#pt5. total_locs	124683373	84197125	2744540	81452585	4012331	186729
#pt5. total_reads	4637235	mapt:4012331,86.52%	uniq:2744540,59.18%	mult:1267791,27.34%	nomap:186729,4.03%
#pt5. cds_reads	24205429	mapt/in:2808223/27793511	cl123/in:1089008,425261,1293954/10515522,4276338,13001651	nomap:1601467
  
#pt6. total_locs	124668747	83891387	2743785	81147602	4011790	186655
#pt6. total_reads	4637161	mapt:4011790,86.51%	uniq:2743785,59.17%	mult:1268005,27.34%	nomap:186655,4.03%
#pt6. cds_reads	24205429	mapt/in:2807662/27793511	cl123/in:1089928,424030,1293704/10515522,4276338,13001651	nomap:1601971

=cut
  
sub MAIN_loop {}
my($lcr,$lcb,$lucg, $lcrbin, $estcrbin, $sumcov, $nsumcov)=(0) x 9;
my(@cbin,@rdlen);

##  <> INPUT is name.covtab from sam2covtab, col format
#    my($cid,$cbinloc,covcds,covte,covx, acovt, acovm, acovu, anns, ids)=split 
#    acov[tmu] are read cov counts, and covcds/te/x are cds/te/rpt?-readid matches, may be missing
#    anns, ids are likely not here but from readAnnots(anntable)

rename($outmeans,"$outmeans.old") if(-f $outmeans); 
open($outmeanh,'>',$outmeans) or die "writing $outmeans";  

# my $inh= *STDIN << use perl <> for @ARGV or STDIN 
while(<>) {
  my @v=split; 
  my $hasivar= (@v >= 6); ##  was >= 8, fix for bad data  BUG 21Aug01, missing 0 last col aCovU so @v == 7
  my $hasann= (@v >= 9);

  if(/^\W/){ 
  
    if($hasivar and not m/[,=]/ and /Cov/i) { @lvar=@v[@IVAR] unless(@lvar); }  
     # @IVAR= (5,6,7) == aCovT	aCovM	aCovU

=item BUG for unmerged inputs @lvar == global col labels
   >> 1st #cdsxchr line parsed for #ChrID ... columns
  #cdsxchr_covtab par: minident=0.65, mindupid=0.98, softclip=0, readlen=150, minaln=97, BIN=100, part=0/1 
  #dropse20cdsbusco_SRR11813283_2_bwa.cc7d.covtab n_readid 1108778, n_nomap 0, n_notcr 0, n_mapok 765004, n_mapbad 349631, n_dupbad 275, n_mismatch 319658, n_insert 1055, n_delete 1593, n_softclip 0, n_intron 0
  #ChrID	Pos	CovT	CovM	CovU	aCovT	aCovM	aCovU
  dropse20uc:g4817792t1	100	0	0	0	67	67	67
  dropse20uc:g4817792t1	200	0	0	0	89	89	89
=cut

    # #pt9.cdsxchr_covtab par: minident=0.65, mindupid=0.98, softclip=0, readlen=150, minaln=97, BIN=100, part=9/24 
    if(/^#/ and /readlen=(\d+)/) { 
      my $prdlen=$1; # part readlen may vary
      push @rdlen,$prdlen; $readlen=$prdlen; #global? store in samvals
    }  
    if(/^#/ and /n_readid/ and /n_mapok/) { # end comm
      my $ln=$_; my($ptn)= $ln =~ m/^#(\S+)/;
      my($p_readid,$p_nomap,$p_mapok,$p_mapbad,$p_dupbad)= map{ my($v)= $ln =~ (m/$_.(\d+)/)?$1:0; $v; } 
        qw(n_readid n_nomap n_mapok n_mapbad n_dupbad);
      $n_readid+= $p_readid; $n_nomap+= $p_nomap; #read stats; o $n_readid= $p_readid if($p_readid); #old way
      $n_mapok+= $p_mapok; $n_mapbad+= $p_mapbad; $n_dupbad+= $p_dupbad; # mapping stats
    }
    next; 
  }
   
  unless($hasivar) { next; }
  my($cr,$cb,@av)=@v; # split; # error if(@av < 6)

  #old: my($an,$ids)=@av[-2,-1]; my($ct,$cm,$cu,$act,$acm,$acu)=@av;  
  my($ct,$cm,$cu)= splice(@av,0,3); # NOW == CDSrd, TErd, UNKrd
  my($cCDS,$cTE,$cUNK)= ($ct,$cm,$cu); # RENAME these vars
  my($act,$acm,$acu)= splice(@av,0,3); # == @IVAR, same  rd_total, rd/mmap, rd_uniq for chrasm x rd map 
  my($an,$ids)= (@av>0) ? splice(@av,0,2) : ("",""); # may not exist

   # any more cols?
  if($nann>0 and my $anidx= $annvals{$cr}{$cb}) {
    my($anx,$idx)=split" ",$anidx;
    if($anx){ $an = ($an)? "$an,$anx" : $anx; }
    if($idx){ $ids= ($ids)? "$ids,$idx" : $idx; }
  } elsif($nidclass) { #  added from idclassf; use chrid or $ids
    $ids=$cr unless($ids); #trick for cds.cov; see below ids and busco
    for my $id (split(",",$ids)) {
      if($id and my $anx= $idclassh->{$id} ) { $an = ($an) ? "$an,$anx" : $anx; last; }
    }
  }
  
  if($an =~ m/RPT/){ $an =~ s/\bRPT\b/repeat/; } # fixup does this occur?
  #t: if($an =~ m/gap/){ $an =~ s/\bgap\b/biggap/; }  # test this for any gap? ie 10/100 bp
  
  # my $MINRD= 5; #? option? skip all cov rows w/ acovt < MINRD, want to measure zeros? gaps?
  # my $MINDUNIQ= 2; # for test of uniqspan == covt <= covu + MINDUNIQ, boost for crappy longread aligns?
  my $hascov= ($act>=$MINRD)?1:0;
  my $isuniq= ($hascov and $act <= $acu + $MINDUNIQ)?1:0;
  $ninrow++;
  $nincov++  if($hascov);
  $nbiggap++ if($an =~ /gap/); # was biggap, instead count any gap spans here; 
  # any act < MINRD but not gap?
  
  #UPD21Jun: new/more annots from gnodes_annotate: 
  #  CDS, GTop | partial < require GTop (best of several gene align), discard partial aln
  #  and 1st id1,id2 id1 is always best align, so change UCG_FIXID
  
if(UCG_KU_SPECIAL) {
  if($ids and $an =~ /busco|UCG/i) { # allow other pats? UCG ?

    if($isuniq) { # isuniq == $act <= $acu + 2
      # add nuniq == ntot filter for UCG: skip if( $nu < $nm-1) .. tests say it helps

if(UPD21JUN) { # updated annots: Gtop|partial, 1st-id
      my($ida)= split",",$ids; # * 1st id == UCG* or use idclass{id} =~ /busco|UCG/ ** now is valid assump
      $ida=0 if($an =~ /partial|lowqual/ or ($ANhasGTOP and not $an =~ /Gtop/));
      if($ida) {
        push @{$ucg_id_cov{$ida}}, $acm; # check act,acu also ? separate @ucg_t,u ?
        $lucg= $ida;
      }
} elsif(UCG_FIXID) {      

# FIXME: use idclass (if avail) and ids to decide, need id hash of @ucg_one_cov
# this doesnt change relative result for 1 test, dplx19ml x dplx20maca, diff in K.UCG; reduces n to n bu gene
# note: that dplx diff may be due to high dupl rate for busco genes: real or false, both dplx19,dplx20 have it
# ..    test refined UCG gene spans, revise annot to reject busco/ucg w/ >1 map loc, or split map, or partial < 90%
# .. spp w/ high tCopy/xCopy  for CDSbus: dplx20mac (2x), dplx19ml (2x), dcari(only 1/2x,2asm), cucumber (3-4x,2asm), 
# .. spp with CDSbus t/xCopy == 1: dropse, apis, aweed (2asm), dplx16ml

      my($ida); my @ids= split",",$ids;
      if($nidclass) { ($ida)= grep{ $idclassh->{$_} =~ m/busco|UCG/i } sort @ids; } 
      else { ($ida)= $ids[0]; } # otherwise old way
      if($ida) {
        push @{$ucg_id_cov{$ida}}, $acm; # check act,acu also ? separate @ucg_t,u ?
        $lucg= $ida;
      }
      
} else {      
      my($ida)=split",",$ids; # *assume 1st id == UCG* or use idclass{id} =~ /busco|UCG/ ** BAD assump
      if($ida ne $lucg) { 
        sum_ucg_cov($lucg,\@ucg_one_cov) if($lucg); 
        @ucg_one_cov=();
      }
      push @ucg_one_cov, $acm; # check act,acu also ? separate @ucg_t,u ?
      $lucg= $ida;
}
    }
  }
}

  if($PLOTCHR and ($lcr ne $cr or $cb > $lcrbin)) {
    put_plotrow($outmeanh,$lcr,$lcrbin,$sumcov,$nsumcov) if($lcr); # print row of part vals, means of filt_sum, then empty sums
    
    if($lcr ne $cr) { $lcrbin=$PLOTBIN; $estcrbin=0; } else { $lcrbin=$PLOTBIN + $cb; }
    $sumcov=$nsumcov=0;
  }
  
  if($lcr eq $cr and $lcb>0 and $cb > $lcb) { 
    if(@cbin < kSAMPLE) { 
      my $bspan= $cb - $lcb; push @cbin, $bspan; 
      my $nbin=@cbin; $CBIN= $cbin[int($nbin/2)] if($nbin>100);
      } 
  }
  ($lcr,$lcb)= ($cr,$cb);
  $sumcov+=$acm; $nsumcov++; # for plotrow, should be == allasm val
  
  ## >> add filter loop here, for @parts == @DEFparts
  ## >> sum() needs filter name, save/filt:  ($nt,@sum,@nit,@lvar,@svar,@sval)
if(FLOOP) {
  for my $filtp (@parts) {
    setFiltVals($filtp); # skip?

    # FIXME: gap annots need to be parsed now : skip all bins w/ gap> 33%? == 33/100 bp; count/total gap spans
    # .. handle gaps in filt_sum() or here ? need annot: if($an =~ /biggap/) { $ok=0; $gapsum{big}+=bspan }
    # FIXME: ann key mixup: RPT == (simple) repeat, UNK = unclassed repeat/te/cds_dupl

    my $ok= $hascov; # == ( $act >= $MINRD )?1:0; # ? but count act < minrd ?
    if($filtp eq 'GAPann') { $ok=1; } # ignore MINRD and biggap here only
    else { $ok=0 if($an =~ /biggap/); } # BUT this cancels GAPann from @FAN below ..

    #UPD21aug06: add here? filter out annot == contam for bacter.contam scaffolds .. use anntab? or separate table?
    if(UPD21AUG) {
    if($filtp eq 'CONTAM') { $ok=($an =~ /contam/)?1:0; }
    else {  $ok=($an =~ /contam/)?0:1; }
    }
    
    $ok= ($URD==1) ? ($ok and $isuniq) : ($URD==2) ? ($ok and !$isuniq) : $ok;  
    $ok= ($CRD==1) ? ($ok and $cCDS > 1) : ($CRD==2) ? ($ok and $cCDS < 2) : $ok;  
    $ok= ($TERD==1) ? ($ok and $cTE > 1) : ($TERD==2) ? ($ok and $cTE < 2) : $ok;  
    
    $ok= ($UNKRD==1) ? ($ok and $cUNK > 1) : ($UNKRD==2) ? ($ok and $cUNK < 2) : $ok;  
    #count UNKRD only if cCDS < 2 and cTE < 2 : this now doen in %filts:
    #    unkasm =>  "urd=0 crd=2 terd=2 unkrd=1 an=0",
    #not this way ok= ($UNKRD==1) ? ($ok and $cUNK > 1 and $cCDS < 2 and $cTE < 2) : ($UNKRD==2) ? ($ok and $cUNK < 2) : $ok;  
    
    if(@FAN>0) { #  and $hasann??? require $an, or @v > 8
      my @pat=@FAN; for my $pat (@pat){ my $pnot=0; if($pat =~ s/^(no|\-)//){ $pnot=1; }  
       if($pnot) { $ok=0 if($an =~ m/$pat/); } else { $ok=0 unless($an =~ m/$pat/); }  } 
    } 
       
    if($ok){ filt_sum( $filtp, $act,$acm,$acu); } # clearer
  } # filtp loop
  
} else {
  my $ok= $hascov; # == ($act>=$MINRD)?1:0; # ? but count act < minrd ?
  
  $ok= ($URD==1) ? ($ok and $isuniq) : ($URD==2) ? ($ok and !$isuniq) : $ok;  
  $ok= ($CRD==1) ? ($ok and $ct > 1) : ($CRD==2) ? ($ok and $ct < 2) : $ok;  
  #* add ($TERD==1) ($UNKRD==1) cases
  if(@FAN>0 and $hasann) { # also require $an, or @v > 8
    my @pat=@FAN; for my $pat (@pat){ my $pnot=0; if($pat =~ s/^(no|\-)//){ $pnot=1; }  
     if($pnot) { $ok=0 if($an =~ m/$pat/); } else { $ok=0 unless($an =~ m/$pat/); }  } 
  } 
     
  if($ok){ #  and $act>4
    if(1){ sum($act,$acm,$acu); } # clearer
  }
}

} # end loop while(<>) covtab input

if(UCG_FIXID or UPD21JUN) {     
  for my $iducg (sort keys %ucg_id_cov) {
    my $covlist= $ucg_id_cov{$iducg}; 
    sum_ucg_cov($iducg,$covlist); delete $ucg_id_cov{$iducg};
  } 
}  

if(@rdlen>1){  my $nl=@rdlen; my $sl=0; map{ $sl+=$_ } @rdlen; $readlen= int($sl/$nl); } 
$samvals{$asmid}{readlen}= $readlen if($readlen);
$samvals{$asmid}{n_readid}= $n_readid; $samvals{$asmid}{n_nomap}= $n_nomap;  

#above: my $nbin=@cbin; $CBIN= $cbin[int($nbin/2)];
my($mbrow,$mbcov,$mbgap)= map{ int( $CBIN*$_/100_000)/10 } ($ninrow,$nincov,$nbiggap);
warn "#covtab span=$mbrow.mb, covspan=$mbcov.mb, gaps=$mbgap.mb info: n_readid=$n_readid, n_nomap=$n_nomap, n_mapok=$n_mapok, n_mapbad=$n_mapbad, n_dupbad=$n_dupbad, \n";

if($PLOTCHR) { #  and ($lcr ne $cr or $cb > $lcrbin)
  put_plotrow($outmeanh,$lcr,$lcrbin,$sumcov,$nsumcov) if($lcr); # print row of part vals, means of filt_sum, then empty sums
  close($outmeanh);

  # ($outmeans=$outsum) =~ s/.sum.*//; $outmeans .= ".plottab"; #??
  # $outsum =~ s/.sum.*//;  $outsum .= ".plotchr.sh"; # Rscript
  
  my $otitle= $outsum; $otitle =~ s/.plotchr.sh//;
  rename($outsum,"$outsum.old") if(-f $outsum); 
  open($outsumh,'>',$outsum) or die "writing $outsum";
  
  # plotChrRScript( $outsumh, $outmeans, $otitle);
  my $rplot=<<"EOS";
#!/usr/bin/env Rscript
source("$gnodes_rplots");
gnodes_chrplot( "$otitle", "$outmeans");
EOS

# source("/XXXXXXX/evigene/scripts/genoasm/../genoasm/gnodes_chrplot.R");
# gnodes_chrplot( "dmag15nwb2asm_pl8f_SRR7825549b", "dmag15nwb2asm_pl8f_SRR7825549b.plottab"); #?? shorten otitle?

  print $outsumh $rplot;
  close($outsumh);
  system("chmod +x $outsum");
  unless( bad_exec("Rscript")) {
  my $err= system("./$outsum"); # exec here 
  warn "#gnodes_chrplot Rscript $outsum, err= $err\n" if $debug;
  }
  exit(0); # done?? none of below means,sum tables
}    

# FIXME: -meanvar outmeans exists, reuse instead of rewrite ? for sumtable
#above: rename($outmeans,"$outmeans.old") if(-f $outmeans); 
#above: open($outmeanh,'>',$outmeans) or die "writing $outmeans";  

if(FLOOP) {
  my @oparts= @parts;
  if(UCG_KU_MEANS) { push @oparts, UCG_KU_MEANS; }
  for my $filtp (@oparts) {
    setFiltVals($filtp); # skip?
    
  # my($nt,@sum,@nit,@svar,@sval,@dhist); # meanvar sum globals: arrays are for 3 vars: act,acm,acu
  # my(%fnt,%fsum,%fnit,%fsvar,%fsval,%fdhist); # globals for filt_sum
    $nt= $fnt{$filtp} || 0;
    warn "#meanvar nt=$nt, ti=$title\n" if($debug > 1);
    next if($nt < 1); #?
    @sum= @{ $fsum{$filtp} }; delete $fsum{$filtp};
    @nit= @{ $fnit{$filtp} }; delete $fnit{$filtp};
    @svar= @{ $fsvar{$filtp} }; delete $fsvar{$filtp};
    @sval= ($DOMED) ? @{ $fsval{$filtp} } : (); if($DOMED) { delete $fsval{$filtp}; }
    @dhist= ($DOHIST) ? @{ $fdhist{$filtp} } : (); if($DOHIST) { delete $fdhist{$filtp}; }

    meanvar($outmeanh);  
  }
  #... ^^ LOOP here to meanvar(out) for each filter before closing outmeans
} else {

  meanvar($outmeanh);  
}

close($outmeanh);

# Next: read $outmeanh, write $outsumh
rename($outsum,"$outsum.old") if(-f $outsum); 
open($outsumh,'>',$outsum) or die "writing $outsum";
meansToSum( $outsumh, $outmeans, $asmid);

# move to meansToSum()
# my($mbrow,$mbcov,$mbgap)= map{ int( $CBIN*$_/100_000)/10 } ($ninrow,$nincov,$nbiggap);
# print $outsumh "# total span=$mbrow.mb, covspan=$mbcov.mb, gaps=$mbgap.mb with n_readid=$n_readid for $asmid\n";

close($outsumh);

#------------------

use constant { DH1k => 1000, DH10k => 10000, DH2k => 2000 }; # Depth histogram cuts

sub sum { 
  my @s=@_; 
  my $nz=0; map{ $nz++ if($_); } @s; 
  return $nz unless($nz); # or $ZEROS ??
  for my $i (0..$#s) {
    my $s=$s[$i]; if($s>0 or $ZEROS) {
    $sum[$i]+=$s; $nit[$i]++; $svar[$i]+=$s*$s; 
    if($DOMED){ push @{$sval[$i]},$s; }
    if($DOHIST){ 
      my $k=($s < DH1k) ? $s # 0..999
          : ($s < DH10k) ? DH1k + int($s/10) # 1000 + (0..999) == 1000..1999
          : DH2k ;  # 2000+ what?
      $dhist[$i][$k]++; # FIXMEd need [$i] == covT/M/U col
      }
    }
  } 
  $nt++; return $nz; 
} 

sub filt_sum { # same calc as sum() but for each $flt
  my ($flt,@s)=@_; 
  my $nz=0; map{ $nz++ if($_); } @s; 
  return $nz unless($nz); # or $ZEROS ??
  for my $i (0..$#s) {
    my $s=$s[$i]; if($s>0 or $ZEROS) {
    $fsum{$flt}[$i]+=$s; $fnit{$flt}[$i]++; $fsvar{$flt}[$i]+=$s*$s; 
    if($DOMED){ push @{$fsval{$flt}[$i]},$s; }
    if($DOHIST){ 
      my $k= ($s < DH1k) ? $s : ($s < DH10k) ? DH1k + int($s/10) : DH2k ; 
      $fdhist{$flt}[$i][$k]++; # FIXME need [$i] == covT/M/U col
      }
    }
  } 
  $fnt{$flt}++; return $nz; 
} 


sub put_plotrow {
  my($outh,$lcr,$lcrbin,$sumcov,$nsumcov)=@_; # print row of part vals, means of filt_sum, then empty sums
  # also sumcov,nsumcov
  
  my $ICOV=($PLOT_COVT)?0:1; # aCovT vs aCovM
  ## see generdcov2pltab.pl putc() .. obs span == input bins (CBIN here), est span == bins * aCovM/CU (== xCopy or rcov)
  # my $bspan= $nbin * 100;  # 100 == inputBINSIZE
  # my $cspan= int($bspan * $rcov);
  my $CUniq= $KUCG || $gucg_med; $CUniq=1 if($CUniq<1);
  my $rcov= ($nsumcov<1) ? 0 : $sumcov/($CUniq * $nsumcov);
  my $obsSpan= $nsumcov * $CBIN;
  my $estSpan= int($obsSpan * $rcov);
  my $estEnd= $estcrbin + $estSpan;
  
  my @pcols= qw(ChrID ObsLoc EstLoc Xtotal);
  my @vcols= ($lcr,$lcrbin,$estEnd, sprintf("%.1f",$rcov));
  #? add Est ChrPos as above, need CU.M/KU.ucg/other uniq cov depth val 
  my @zero=(0,0,0);
  
  for my $filtp (@parts) {
    # not for output: setFiltVals($filtp); # skip?
    
    my $nt= $fnt{$filtp} || 0; delete $fnt{$filtp};
    # warn "#meanvar nt=$nt, ti=$title\n" if($debug > 1);
    # next if($nt < 1); #?
    my @sum= ($nt < 1)? @zero : @{ $fsum{$filtp} }; delete $fsum{$filtp};
    my @nit= ($nt < 1)? @zero : @{ $fnit{$filtp} }; delete $fnit{$filtp};
    #my @svar= @{ $fsvar{$filtp} }; 
    delete $fsvar{$filtp};
    #my @sval= ($DOMED) ? @{ $fsval{$filtp} } : (); 
    if($DOMED) { delete $fsval{$filtp}; }
    #my @dhist= ($DOHIST) ? @{ $fdhist{$filtp} } : (); 
    if($DOHIST) { delete $fdhist{$filtp}; }
    # meanvar($outmeanh);  
    ## PLOT cols want only @nit[1], @sum[1] == aCovM val ** OPTION to use i=0 == aCovT
    # my $i= $ICOV;
    my $nit=$nit[$ICOV]||0; 
    my $s= $sum[$ICOV]||0; 
    my $aveXCopy= ($nit<1)? 0: $s/($nit*$CUniq); # want xCopy? ave= $s/($nit*$gucg_med)
    # if($nit>3){ $ave=$s/$nit; $var=($svar[$i] - $ave*$ave)/($nit-1); $sd=sqrt($var); $se=$sd/sqrt($nit); }
    
    push @pcols, 'X'.$filtp, 'N'.$filtp;
    push @vcols, sprintf("%.1f",$aveXCopy),$nit; 
  }
  
  print $outh join("\t",@pcols)."\n" unless($outhdr++);
  print $outh join("\t",@vcols)."\n";
  $estcrbin= $estEnd;
  #?? $bat += $bspan; $eat += $cspan; 
}

# sub plotChrRScript {
#   my( $outh, $chrcovtable, $otitle)= @_;
#   #? read R subs from $evigene/scripts/genoasm/gnodes_chrplot.R
#   # rather than stuff in here
# 
#   my $rplot=<<"EOS";
# #!/usr/bin/env Rscript
# 
# source("$gnodes_rplots");
# gnodes_chrplot( $otitle, $chrcovtable);
#   
# EOS
#     print $outh $rplot;
# }

# update for DOHIST: meanvar() find peak/mode, est med of D[0..999],Dk[1k..10k],Dx[>10k]

sub meanvar { 
  my($outh)= @_;
  my($med,$nit,$ave,$var,$sd,$se,$s);
  print $outh "$title stats for nt=$nt\n"; 
  
  #   if(0 and $NTAVE) {  # not used here
  #     my $lvar=($showSUM)?"Sum":"Var"; # add Range?
  #     print $outh join("\t",qw(Item__ Mean SEM StDev),$lvar)."\n"; 
  #     for my $i (0..$#sum){ $lvar=$lvar[$i]||"it$i"; $s=$sum[$i]; 
  #      $ave=$s/$nt; $var=($svar[$i] - $ave*$ave)/($nt-1); $sd=sqrt($var); $se=$sd/sqrt($nt); 
  #     $var=$s if($showSUM);
  #     printf $outh "%-6s\t%.2f\t%.2f\t%.2f\t%.2f\n",$lvar,$ave,$se,$sd,$var; } 
  #   }
  
  my @hdr=qw(Mean SEM Nitem StDev); 
  #? new cols: Median Mean SEM Nitem StDev Mode Sum|Var|Range
  if($DOHIST) { unshift(@hdr,"Median");  push(@hdr,"Mode"); } 
  elsif($DOMED) { unshift(@hdr,"Median"); }
  push(@hdr,($showSUM)?"Sum":"Var"); # add Range?
  
  my $fmt="%.2f\t%.2f\t%d\t%.2f";
  if($DOHIST) { $fmt="% 4d\t".$fmt."\t%4d"; } 
  elsif($DOMED) { $fmt="% 4d\t".$fmt; }
  $fmt.=($showSUM)?"\t%d":"\t%.1f"; 
   
  print $outh join("\t","Item  ",@hdr)."\n"; 
  for my $i (0..$#sum) { 
    my $litem= $lvar[$i]||"it$i"; $s=$sum[$i]; $nit=$nit[$i]; ($ave,$var,$sd,$se)=(0) x 4;
    #UPD2202: Uhgg bad old err: var=(ssq - av*av) SHOULD BE (ssq - N*av*av)
    #OLD: if($nit>3){ $ave=$s/$nit; $var=($svar[$i] - $ave*$ave)/($nit-1); $sd=sqrt($var); $se=$sd/sqrt($nit); }
    if($nit>0){ $ave=$s/$nit; if($nit>1){ $var=($svar[$i] - $nit*$ave*$ave)/($nit-1); $sd=sqrt($var); $se=$sd/sqrt($nit); } }
    my @mv=($ave,$se,$nit,$sd); 
    if($DOMED) { my @ss=sort{$a <=> $b} @{$sval[$i]}; $med= $ss[ int($nit/2) ];  unshift(@mv,$med); }
    if($DOHIST) {
      my($kmed, $kpeak, $npeak, $ntot, $nhalf, $nat)=(0) x 9;
      my $kz= ($ZEROS) ? 0 : $MINRD;  # count zeros?
      #  DH1k => 1000, DH10k => 10000, DH2k => 2000 
      for my $k ($kz .. DH2k) { $ntot += $dhist[$i][$k]||0; } $nhalf= int($ntot/2);
      for my $k ($kz .. DH2k) {
        my $kr= $k;
        my $nk= $dhist[$i][$k]||0;
        if($k >= DH2k) { $kr= DH10k; } # ceiling for this
        elsif($k >= DH1k) { $kr= ($k - DH1k) * 10; } # ugh, was k*10, 
          # kr[DH1k] should be k*10 - DH1k; from DH1k + int($s/10) ; ie 2000 => 1200 ; (1200 - 1000)*10 = 2000
        if($nk > $npeak and $k < DH1k){ $npeak= $nk; $kpeak= $kr; } # peak DONT count >= DH2k or >= DH1k?
        if($nat < $nhalf){ $kmed= $kr; }
        $nat += $nk;
      }
      unless($DOMED) { unshift(@mv,$kmed); }
      push(@mv, $kpeak);
    }
    push @mv,(($showSUM)?$s:$var);
    
    printf $outh "%-6s\t$fmt\n",$litem,@mv;
    } 
  print $outh "\n"; #------ ?
}

# my($nann,%annvals)= (0);
# ($nann)= readAnnots($anntable) if($anntable); # format from gnodes_annotate.pl, same crID, crLoc, an, ids as covtab
sub readAnnots {
  my($anntable)= @_;
  my($nann,$lpt)=(0,0); 

#if(UPD21AUG and $ANhasCONTAM) : A_CONTAM => 'contam', new annot
#UPD21JUN: 
# use constant { # gnodes_annotate.pl current annot flags : export?
#     A_CDS => 'CDS', A_TOP => 'Gtop', A_PART=>'partial', A_LOWQUAL=>'lowqual',
#     A_BUSCO => 'busco', A_UCG => 'UCG',
#     A_GAP => 'gap',
#     A_TE => 'TE', # gff.type=transposon
#     A_REPEAT => 'repeat', # gff.type=Simple_repeat or repeat
#     A_UNK => 'UNK', # unknown/unclassified type (from repmask only?)
# }; 

  warn "#read anntable=$anntable\n" if($debug);
  open(F,$anntable) or do { warn "#err: missing anntable=$anntable\n"; return 0; };
  while(<F>){
    next if(/^\W/); my @v= split;
    my($cid,$cib,$an,$ids)= @v; # any more than an,ids ?
    next unless($an or $ids);
    $annvals{$cid}{$cib}="$an\t$ids"; $nann++;
    map{ $anntypes{$_}++; } split",",$an;
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

#update to pipe.pl vers
sub readMetad {
  my($sampledata)= @_;
  my($nsam,$aid)=(0,0); 
  warn "#read sampledata=$sampledata\n" if($debug);
  open(F,$sampledata) or do { warn "#err: missing sampledata=$sampledata\n"; return 0;  };
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
  
  return($nsam); # return ($nsam, \%metavals); #?
}

 
sub _max{ return($_[0] < $_[1])?$_[1]:$_[0]; }
sub _min{ return($_[0] > $_[1])?$_[1]:$_[0]; }
sub _minnot0{ return ($_[0] == 0) ? $_[1] : ($_[1] == 0) ? $_[0] : _min(@_); }
#bad: sub _minnot0{ return($_[0] > $_[1] and $_[1] > 0)?$_[1]:$_[0]; }

sub sum_ucg_cov { # UCG_KU_SPECIAL
  my($uid, $ucovs)= @_;
  my $ncov= @$ucovs;
  if($ncov>2) { # MIN_UCG_BIN => 2 
    my @cu= sort{ $b<=>$a } @$ucovs;
    my $cmed= $cu[ int($ncov/2) ];
    push @ucg_med_cov, $cmed; # global median depth for all UCG
    push @ucg_med_ids, $uid; # ? not used, or check for dups?
    $nucg=@ucg_med_cov;
if(UCG_KU_MEANS) { # also do filt_sum + meanvar ?? only for 1 val cmed per gene
    filt_sum( UCG_KU_MEANS, $cmed,$cmed,$cmed); # ?? $act,$acm,$acu
}    
  }
}

sub KU_ucg_cov { # UCG_KU_SPECIAL
  $nucg=@ucg_med_cov; # require nucg > 100? 500?
  my $MINUCG=$ENV{minucg}||90; # was 3, hit that w/ pacbio no-uniq long err aligns, bad result
  if($nucg >= $MINUCG) {
    my @cu= sort{ $b<=>$a } @ucg_med_cov;
    my $cmed= $cu[ int($nucg/2) ];
    my $cmer= 0; # what is error stat for med-of-med?
    # some stat web sez: se(med)= sqrt(var), var= ssq - ss/n
    my ($cm,$css)=(0,0); for(my $i=0; $i<$nucg; $i++)  { my $cu=$cu[$i]; $cm+=$cu; $css += $cu*$cu; }
    $cm=$cm/$nucg;  $cmer= sqrt( ($css - $cm*$cm)/$nucg );
    return($cmed, $cmer, $nucg, $cm);
  }
  return(0);
}  

#UPD22apr: large range valid 0.65,1.8.. 
sub rangeKuni {
  my($kutype,$kulo,$kuhi)=@_;
  my($kumin,$kumax)=(0,999999); # use when ku is set in metadata
my $pLOr=0.50; # was 0.80;
my $pHIr=2.00; # was 1.20 .. 1.10
  if($kutype =~ /^set/)   #o:if($kutype eq "set") 
  {  # test, replace w/ kbu_med if it seems accurate: it *might be* bogus
    $kumin= ($kulo == 0) ? $kuhi * $pLOr : $kulo * $pLOr;
    $kumax= ($kuhi == 0) ? $kulo * $pHIr : $kuhi * $pHIr;
    #o $kumax= ($kuhi == 0) ? $kulo*1.20 : ($kuhi < $kumin*1.5) ? $kuhi : $kuhi * 1.10;
  }
  return($kumin,$kumax);
}

sub read_genexcopy {
  my($genexcopy)=@_;
  my($nucg,$cmd,$cav,$cerr)= (0) x 9;
  open(F,"head -50 $genexcopy |") or return 0;
  while(<F>) {
    # Uniq Gene Cov Depth n=954, .. C.Map/W=38.4, 38.0 +/-1.26 (mdn,ave,sem) 
    if(/Uniq Gene Cov Depth n=(\d+)/) {  
      $nucg=$1;
      # if( m/C.Map.W=([\d\.]+), ([\d\.]+)....([\d\.]+)/ ) { ($cmd,$cav,$cerr)= ($1,$2,$3);}
      if( m/C.Map.W=([\d\.]+), ([\d\.]+)\W+([\d\.]+)/ ) { ($cmd,$cav,$cerr)= ($1,$2,$3);}
      elsif( m/C.Map.W=([\d\.]+)/) { $cmd=$1; }      
      last if($cmd);
    }
  } close(F);
  return($nucg,$cmd,$cav,$cerr); 
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

=item pickKuni

  UPD220227:
  * problematic try to pick best C standard depth from several stats
  * now includes general asm parts: uniqasm, NOann, allasm which may have consistent, but non-standard vals
  * test limit pickKuni to gene-sets: CMgenes (from genexcopy), KUucg that seems best est on chrasm, and CDSbus
     unless not available, then use uniqasm
   ie. add valid-rank to %kuept = [ cmd, cav, nval, edam:cmd-cav, sort-rank:3,2,1,0 ]
   OR use weight for edam, instead of sort-rank?
     eg. CDSbus,wt=1.5 vs KUucg,wt=1.25 vs CMgenes,wt=1
     best = min( edam*wt) ?
   
=cut

sub pickKuni {
  my($asmid)= @_;
  
  my %kuept=(); #was %ept=();
  
  use constant { rank_other => 0, # UPD22FEB
    rank_KUset => 3, rank_CMgenes => 3, 
    rank_KUucg => 2, rank_CDSbus => 2, rank_uniqasm => 1, 
    }; 
  
  my($kucg_med, $kucg_err, $kucg_n, $kucg_ave)= (0) x 4; # UCG_KU_SPECIAL, also genexcopy?
  #now global: my($gucg_n, $gucg_med, $gucg_ave, $gucg_err)= (0) x 4; # genexcopy
  my($kumin,$kumax)=(0,999999); # use when ku is set in metadata

  # add user opt: $KUCG here? replace kulo/kuhi?
  my $kulo= $samvals{$asmid}{kulo}|| $KUCG; # int or real val
  my $kuhi= $samvals{$asmid}{kuhi}|| 0;
  if($kulo and $kuhi) { ($kulo,$kuhi)= ($kuhi,$kulo) if($kulo > $kuhi); }
  my $kutype=($kulo or $kuhi)?"set":"";
  if($kulo or $kuhi){ $kuept{'KUset'} = [$kulo,$kuhi,1,0,rank_KUset]; } #???
  
  if(UPD21JUN and $genexcopy) {
    unless($gucg_n) {
    ($gucg_n, $gucg_med, $gucg_ave, $gucg_err)= read_genexcopy($genexcopy); 
    warn "#pickKuni Genes gucg_med=$gucg_med, sd(gucg_med)=$gucg_err, gucg_ave=$gucg_ave, kucg_n=$gucg_n\n" if $debug;
    }
    if($gucg_n and $gucg_med) {
      $kuept{'CMgenes'} = [$gucg_med,$gucg_ave,$gucg_n,0,rank_CMgenes]; 
      if(1) {
        if($kulo == 0) { $kulo= $gucg_med; }
        elsif($kulo < $gucg_med) { $kuhi= _max($kuhi,$gucg_med); }
        else { $kulo= _minnot0($kulo,$gucg_med); }
        if($kulo and $kuhi) { ($kulo,$kuhi)= ($kuhi,$kulo) if($kulo > $kuhi); }
        $kutype= "set_UCG";
      }
      $samvals{$asmid}{kucg_med}= $gucg_med;
      $samvals{$asmid}{kucg_ave}= $gucg_ave;
      $samvals{$asmid}{kucg_sem}= $gucg_err;
      $samvals{$asmid}{kucg_n}= $gucg_n;
      }
  }
  
  ($kumin,$kumax)= rangeKuni($kutype, $kulo, $kuhi);
  warn "#pickKuni samvals kulo=$kulo, kuhi=$kuhi\n" if $debug;
 
  # fixme DOHIST: mode also
  # if DOMED == 0  uav->[1] == ave not med, [2] == SEM
  my($kuasm_med,$kuasm_ave,$kuasm_n)=(0,0,0);
  if(my $pv= $sum{'uniqasm'}{'rdacovm'}) { # << FIXME rdacovm OR aCovM 
    #  $kuasm_med= $pv->[1]||0;  $kuasm_ave= ($DOMED||$DOHIST) ? $pv->[2]||0 : 0;
    ($kuasm_med,$kuasm_ave,$kuasm_n)= @{$pv}[1,2,4];
    warn "#pickKuni uniqasm kulo=$kuasm_med, kuhi=$kuasm_ave\n" if $debug;
  }
  
  my($kbu_med,$kbu_ave,$kbu_n)=(0,0,0);
  if(my $pv= $sum{'CDSbus'}{'rdacovm'}) {  # was bad: rdacovt, use rdacovm
    # $kbu_med= $pv->[1]||0; $kbu_ave= ($DOMED||$DOHIST) ? $pv->[2]||0 : 0;
    ($kbu_med,$kbu_ave,$kbu_n)= @{$pv}[1,2,4];

unless(UPD21JUN) { # dont set kulo/hi yet    
    if($kutype eq "set") {  # test, replace w/ kbu_med if it seems accurate: it *might be* bogus
      ($kumin,$kumax)= rangeKuni($kutype, $kulo, $kuhi);
      if($kbu_med >= $kumin and $kbu_med <= $kumax) {
        $kutype= "set_CDSbus";
        if($kulo == 0) { $kulo= $kbu_med; }
        else { 
          if(abs($kbu_med - $kulo) < abs($kbu_med - $kuhi)) { $kulo= $kbu_med; } 
          else { $kuhi= $kbu_med; } 
        }
      }
    }
}    
    warn "#pickKuni CDSbus kulo=$kbu_med, kuhi=$kbu_ave\n" if $debug;
  } 

  if(UCG_KU_SPECIAL and UPD21JUN) { # add KU_ucg_cov to kuept
    ($kucg_med, $kucg_err, $kucg_n, $kucg_ave)= KU_ucg_cov(); #UPD7f add ave? 
    if($kucg_n >= 90 and $kucg_med > 0) {
    $kuept{'KUucg'} = [$kucg_med,sprintf("%.2f",$kucg_ave),$kucg_n,0,rank_KUucg];

    # UPD21JUN move ahead of kc_bestpart instead of after, rely more on kc_bestpart method
    # UPD21JUN: skip this set, rely on bestpart method; one of kulo/kuhi == 0, if KUucg is way off, get wrong est
#     ($kumin,$kumax)= rangeKuni($kutype, $kulo, $kuhi);
#     if($kucg_med >= $kumin and $kucg_med <= $kumax) { # drop: $kutype =~ /set/ and 
#       if($kulo == 0) { $kulo= $kucg_med; }
#       elsif($kuhi == 0) { $kuhi= $kucg_med; }
#       else { #? leave this to kc_bestpart
#         $kulo= _minnot0($kulo,$kucg_med);  
#         $kuhi= _max($kuhi,$kucg_med);   
#       }
#       if($kulo and $kuhi) { ($kulo,$kuhi)= ($kuhi,$kulo) if($kulo > $kuhi); }
#       $kutype= "set_UCG";
    
#     } elsif($kbu_med) { # tends to be higher than kucg_med, miscalc?
#       $kulo= _minnot0($kbu_med,$kucg_med); #?? 
#       $kuhi= _max($kbu_med,$kucg_med);  
#       $kutype= "set_UCG";
#     }
    
    #? add kucg_med +/- sem as range?
    my $kusem= ($kucg_n < 4) ? 0 : sprintf"%.3f",$kucg_err/sqrt($kucg_n);
    ($kucg_err,$kucg_ave)= map{ sprintf"%.3f",$_ } ($kucg_err,$kucg_ave);
    warn "#pickKuni UCG ku_med=$kucg_med, sd(ku_med)=$kucg_err, ku_ave=$kucg_ave, sem=$kusem, kucg_n=$kucg_n\n" if $debug;
    
    unless(UPD21JUN and $gucg_med) {
      $samvals{$asmid}{kucg_med}= $kucg_med;
      $samvals{$asmid}{kucg_ave}= $kucg_ave;
      $samvals{$asmid}{kucg_sem}= $kusem;
      $samvals{$asmid}{kucg_n}= $kucg_n;
      }
    }
  }  
  
  #?? add another ku method? check other sum parts for median ~= mean, prefer over median << mean or >> mean
  # test median == mean for sumparts: uniqasm CDSbus CDSann NOann allasm 
  my($kc_med,$kc_ave,$kc_n,$kc_bestpart)=(0,0,0);
  my($kc2med,$kc2ave,$kc2n,$kc2bestpart)=($kuasm_med,$kuasm_ave,$kuasm_n,'uniqasm'); # UPD22FEB : replace below uniqasm kulo=$kuasm_med w/ kc2med
      # ^^ default to uniqasm ??
  if(UPD7f or UPD21JUN) {
    my $MINN=$ENV{minucg}||90; # see MINUCG above
    my @ptadd=qw(CDSann NOann allasm);
    $kuept{'CDSbus'} = [$kbu_med,$kbu_ave,$kbu_n,0,rank_CDSbus] if($kbu_n >= $MINN);
    $kuept{'uniqasm'}= [$kuasm_med,$kuasm_ave,$kuasm_n,0,rank_uniqasm] if($kuasm_n >= $MINN);
              
    for my $pt (@ptadd) {
      if(my $pv= $sum{$pt}{'rdacovm'}) { 
        my($pmed,$pave,$pn)= @{$pv}[1,2,4];
        $kuept{$pt}= [$pmed,$pave,$pn,0,rank_other] if($pn >= $MINN);
      }
    }
    
    my @pt= sort keys %kuept;
    for my $pt (@pt){
      my($md,$av,$n)= @{$kuept{$pt}}; 
      #o my $edam= sprintf"%.4f",abs($av - $md)/sqrt($n); # error term
      # err / sqrt(n) bad when n = num genes vs num genome bins !* not comparable, drop n?
      my $edam= sprintf"%.4f",abs($av - $md); # error term
      $kuept{$pt}->[3]= $edam;
    }
    if(UPD22FEB) {
      @pt= sort{ $kuept{$b}->[4] <=> $kuept{$a}->[4] # pick-rank, b-high wins, all zero default 
        or $kuept{$a}->[3] <=> $kuept{$b}->[3] # edma error, low wins
        or $a cmp $b }  keys %kuept;
   
    } else {
      @pt= sort{ $kuept{$a}->[3] <=> $kuept{$b}->[3] or $a cmp $b }  keys %kuept;
    }
    $kc_bestpart= $pt[0]; ($kc_med,$kc_ave,$kc_n)= @{$kuept{$kc_bestpart}};

    ## UPD21JUN: maybe pick also 2nd best if 1st == prior kulo/hi, ie CMgenes or KUucg
    # ($kc2med,$kc2ave,$kc2n,$kc2bestpart)= (@pt>1) ? (@{$kuept{$pt[1]}},$pt[1]) : (0,0,0);
    if(@pt > 1) {
      ($kc2med,$kc2ave,$kc2n,$kc2bestpart)= (@{$kuept{$pt[1]}},$pt[1]); 
      # if kc2n < MINN, try pt2 ?
    }
    if($kc2n < $MINN and $kuasm_n >= $MINN) { # UPD22FEB $kuasm_med,$kuasm_ave,$kuasm_n,0,rank_uniqasm
      ($kc2med,$kc2ave,$kc2n,$kc2bestpart)= ($kuasm_med,$kuasm_ave,$kuasm_n,"uniqasm");
    }
  
   
    ## UPD21JUN add this to sum.txt output:  $samvals{$asmid}{'pickKUsum'}
    my $pkusum="A Best C Depth Est. from part: $kc_bestpart\n";
    if(UPD22FEB) {
    $pkusum .= "Part\tCmed\tCave\tCn\tCerr\tCrank\n";
    } else {
    $pkusum .= "Part\tCmed\tCave\tCn\tCerr\n";
    }
    for my $pt (@pt){ $pkusum .= join("\t",$pt,@{$kuept{$pt}})."\n"; }
    $samvals{$asmid}{'pickKUsum'}= $pkusum;
    warn $pkusum if($debug);
    
    ($kumin,$kumax)= rangeKuni($kutype, $kulo, $kuhi);
    if($kc_med >= $kumin and $kc_med <= $kumax) {
      $kutype= "set_".$kc_bestpart;
      if($kulo == 0) { $kulo= $kc_med; }
      elsif($kuhi == 0) { $kuhi= $kc_med; }
      elsif(abs($kc_med - $kulo) < abs($kc_med - $kuhi)) { $kulo= $kc_med; } 
      else { $kuhi= $kc_med; } 
      if($kulo and $kuhi) { ($kulo,$kuhi)= ($kuhi,$kulo) if($kulo > $kuhi); }
      
      if($kulo == $kuhi and $kc2med >= $kumin and $kc2med <= $kumax) { # UPD21JUN
        if($kc2med > $kulo){ $kuhi= $kc2med; }
        else { $kulo= $kc2med; }
      }
    }
    
  }
  
  if(UCG_KU_SPECIAL and not UPD21JUN) { #moved before kc_bestpart
    ($kucg_med, $kucg_err, $kucg_n, $kucg_ave)= KU_ucg_cov(); #UPD7f add ave?
    # these should take precedence over CDSbus, kbu_med and kuasm_med
  
    ($kumin,$kumax)= rangeKuni($kutype, $kulo, $kuhi);
    
    if($kucg_med == 0) {
    
    } elsif($kutype =~ /set/ and $kucg_med >= $kumin and $kucg_med <= $kumax) {
      if($kulo == 0) { $kulo= $kucg_med; }
      elsif($kuhi == 0) { $kuhi= $kucg_med; }
      else { 
        $kulo= _minnot0($kulo,$kucg_med);  
        $kuhi= _max($kuhi,$kucg_med);   
        }
      $kutype= "set_UCG";
#     } elsif($kbu_med) { # tends to be higher than kucg_med, miscalc?
#       $kulo= _minnot0($kbu_med,$kucg_med); #?? 
#       $kuhi= _max($kbu_med,$kucg_med);  
#       $kutype= "set_UCG";
    }
    
    #? add kucg_med +/- sem as range?
    my $kusem= ($kucg_n < 4) ? 0 : sprintf"%.3f",$kucg_err/sqrt($kucg_n);
    ($kucg_err,$kucg_ave)= map{ sprintf"%.3f",$_ } ($kucg_err,$kucg_ave);
    warn "#pickKuni UCG ku_med=$kucg_med, sd(ku_med)=$kucg_err, ku_ave=$kucg_ave, sem=$kusem, kucg_n=$kucg_n\n" if $debug;
    
    unless(UPD21JUN and $gucg_med) {
    $samvals{$asmid}{kucg_med}= $kucg_med;
    $samvals{$asmid}{kucg_ave}= $kucg_ave;
    $samvals{$asmid}{kucg_sem}= $kusem;
    $samvals{$asmid}{kucg_n}= $kucg_n;
    }
  }

  if($kulo == 0) {
    $kulo= _minnot0( _minnot0($kbu_med,$kbu_ave),$kc2med); # UPD22FEB kuasm_med>kc2med check kbu_ave also ?
  }
  if($kuhi == 0) {
    $kuhi= _max($kbu_med,$kc2med);   # ? but not kbu_ave may be hi due to dupls
    $kuhi= $kulo unless($kuhi);
  }
  $kulo= $kuhi if($kuhi>0 and $kulo == 0);
  
  if($kulo and $kuhi) { ($kulo,$kuhi)= ($kuhi,$kulo) if($kulo > $kuhi); }
  return($kulo,$kuhi,$kutype); # and other? name tag?
}

sub dumpsum {
  my($outh, $asmid,$np,$nr)= @_;
  
  # %sum{ @parts } = rows of [@statv] cols = @stath
  #  $sum{$stn}{$sti}= [@statv];  # this way?

# add from input covtab stats:
# my($mbrow,$mbcov,$mbgap)= map{ int( $CBIN*$_/100_000)/10 } ($ninrow,$nincov,$nbiggap);
  my $mbCBIN= $CBIN / 1_000_000;

# warn "#covtab span=$mbrow.mb, covspan=$mbcov.mb, gaps=$mbgap.mb info: n_readid=$n_readid, n_mapok=$n_mapok, n_mapbad=$n_mapbad, n_dupbad=$n_dupbad, \n";
# total span=119.6.mb, covspan=115.1.mb, gaps=0.1.mb with n_readid=53925998 for arath18tair_chr
# FIXME: add calc Formula_LN/C= (rdlen * n_readid) / (KU * 1_000_000) = (150*53925998)/(52*1000000) = 155.55 Mb
# add range w/ KU sem=1.3, LN/C=50..54 = 150 Mb .. 162 Mb

  my($kulo,$kuhi,$kutype)= pickKuni($asmid); #NOTE kulo gives HIGHER xCopy and estMb than kuhi, 
  my $twoku= ($kulo ne $kuhi)?1:0; my $oneku= not $twoku;
  
  my $flowcyto= $samvals{$asmid}{flowcyto}|| 0; 
  my $glnc= $samvals{$asmid}{glncformula}|| 0; 
  my $asmtotal= $samvals{$asmid}{asmtotal}|| 0; 
  my $atotalmb= $samvals{$asmid}{atotalmb}|| 0; 
  my $asmname = $samvals{$asmid}{asmname}|| $asmid; 

  #  Formula_LN/C= (rdlen * n_readid) / (KU * 1_000_000)  = (150*53925998)/(52*1000000) = 155.55 Mb
  my $readlen=  $samvals{$asmid}{readlen} || $readlen; # zero? # use global n_readid 
  # est_lnc > gse_lnc? 'genome-size-est'
  # my ($est_lnc,$low_lnc)= (0,0);
  my($gsen, $gsem, $gse_lnc, $gsf_lnc,$gse_lmc,$gsf_lmc)=(0) x 9; # gse = size est/hi, gsf=size floor/lo, lnc=len * Nread/C, lmc=len * Nmaprd/C
  my $n_maptrd= $n_readid - $n_nomap;
  
  # UPD22FEB: remove CONTAM from Formula_LN/C=$glnc also
  my ($contamb,$contamCovM,$contamHi,$contamLo)= (0,0,0,0);
  if(UPD21AUG and defined $sum{'CONTAM'}) {
    my $mv= $sum{'CONTAM'}{'rdacovm'} || [];
    my($cnta,$cntn)= ($mv->[2], $mv->[4]);
    $contamb= sprintf "%.0f",  $mbCBIN * $cntn||0;
    $contamCovM= $mbCBIN * $cntn * $cnta; # esthi,lo = ccovm/($kulo,$kuhi)
    ($contamHi,$contamLo)= map{ ($_ > 0 ) ? sprintf("%.1f", $contamCovM/$_) : 0 } ($kulo,$kuhi);
  }
  
  if($readlen and $kuhi) {
    my $NLmb= $readlen * $n_readid / 1_000_000;
    my $MLmb= $readlen * $n_maptrd / 1_000_000;
    my $CONTMlab="";
    if(UPD22FEB and $contamb){
    ($gse_lnc,$gsf_lnc)= map{ ($_ > 0 ) ? sprintf("%.1f",($NLmb - $contamCovM)/$_) : 0 } ($kulo,$kuhi);
    ($gse_lmc,$gsf_lmc)= map{ ($_ > 0 ) ? sprintf("%.1f",($MLmb - $contamCovM)/$_) : 0 } ($kulo,$kuhi);
    $CONTMlab=" -Contam";
    } else {
    ($gse_lnc,$gsf_lnc)= map{ ($_ > 0 ) ? sprintf("%.1f",$NLmb/$_) : 0 } ($kulo,$kuhi);
    ($gse_lmc,$gsf_lmc)= map{ ($_ > 0 ) ? sprintf("%.1f",$MLmb/$_) : 0 } ($kulo,$kuhi);
    }
    $gsen= ($twoku) ? "$gsf_lnc-$gse_lnc" : $gse_lnc; 
    $gsem= ($twoku) ? "$gsf_lmc-$gse_lmc" : $gse_lmc; # FIXMEd: gse,gsf swapt for Nmaprd
    map{ $_ .= " Mb$CONTMlab" } ($gsen,$gsem);
    $glnc= $gsem; # $gsen; # use nread or nmapt ? replace samval ??
  }
  
  unless($RTABLE){
  print $outh "Source=$asmname, KUlow=$kulo, KUhigh=$kuhi, FlowcytSize=$flowcyto Formula_LN/C=$glnc ($asmid)\n";
  #upd.opt? replace KUlo/hi w/ C_UCG=nn +/- se,eg: C_UCG=50 +/-1.3, 
  # my($CUmed,$CUsem,$CUn)= map{ $samvals{$asmid}{$_} } qw(kucg_med kucg_sem kucg_n);
  # print $outh "Source=$asmname, C_UCG=$CUmed +/-$CUsem, FlowcytSize=$flowcyto Formula_LN/C=$glnc ($asmid)\n";
  }
  if($kulo == 0 or $kuhi == 0) { warn "ERR: kulow, kuhi\n"; return -1; }

  my $dv=join"", ('-') x 75; # w/ Descript
  my $didhd=0;
  
  if(@parts < 6) { @parts=@DEFparts; }
  for my $sp (sort keys %sum) { 
    unless(grep{ $sp eq $_ } @parts) { push @parts, $sp; }
  }
  

  for my $sp (@parts) { 
    next  unless($sum{$sp});
    my $mv= $sum{$sp}{'rdacovm'} || [];
    my $tv= $sum{$sp}{'rdacovt'} || [];
    # Item  	Median	Mean	SEM	Nitem	StDev	Sum
    my @MCOL= ($DOMED||$DOHIST) ? (1,2,4,5) : (1,1,3,4);
    my( $med,$ave,$nit,$sd)= map{ $mv->[$_]||0 } @MCOL;
    my( $tmed,$tave,$tnit,$tsd)= map{ $tv->[$_]||0 } @MCOL;
    
    my $spdesc= $filtdesc{$sp}||"";
    
    # covn,covm = xCopy ratio of excess copy num
    my($xchi,$xclo) = map{ $ave/$_ } ( $kulo, $kuhi); 
    
    my $covt= $tave/$kulo; # and kuhi?
    my $covte= $covt * $xclo; # skip: est Total Copy = obs * xCopy 
    
    my $mb= $nit * $mbCBIN; # == binsize=100/MB=1_000_000; # == 1/10_000 == 0.0001
    #^^ FIXME use CBIN: see above mb = int( $CBIN*$_/100_000)/10 
    #o: if($sp eq 'allasm' and $atotalmb>$mb){ $mb=$atotalmb; } # show both, obs mb = nit, and atotalmb
    
    
    my($xxchi,$xxclo) = ($xchi,$xclo);
    
    my $mbf=($mb > 20)? "%.0f" : "%.1f";
    my($esthi,$estlo) = map{ sprintf $mbf, $mb * $_; } ($xchi,$xclo);
    ($mb)= sprintf $mbf, $mb;
    ($xchi,$xclo,$covt,$covte)= map{ sprintf "%.2f", $_ }  ($xchi,$xclo,$covt,$covte);

    if($RTABLE){
      unless($didRhd++) {
      my @h= qw( asmname asmpart obsmb estmb estlo esthi xcopy tcopy);  
      print $outh join("\t", @h)."\n";
      }
      my $estav= sprintf $mbf, ($esthi+$estlo)/2;
      my @v= ( $asmname, $sp, $mb, $estav, $estlo, $esthi, $xchi, $covt);  
      print $outh join("\t", @v)."\n";
    } else {    
      unless($didhd++) {
        my @hda= qw( _____ ______ Low Low     High High   Total);
        my @hdb= qw( Item_ Obs.Mb Est.Mb xCopy Est.Mb xCopy tCopy);
        if($oneku){ splice(@hda,2,2); map{ s/High/    / } @hda; splice(@hdb,2,2); }
        if($spdesc) { push @hdb,"Description"; push @hda," "; }
        print $outh join("\t",@hda)."\n";
        print $outh join("\t",@hdb)."\n";
        print $outh "$dv\n";
      }
      my @v= ( $sp, $mb, $estlo,  $xclo, $esthi, $xchi, $covt, $spdesc);  
      if($oneku){ splice(@v,2,2); }
      print $outh join("\t", @v)."\n";
    }
    
    if($sp eq 'allasm' and $atotalmb>=$mb) { # add new row here
      ## dont like this way, instead use $dif = $atotalmb - $mb, add to est..
      #my($esthi,$estlo) = map{ sprintf "%.1f", $atotalmb * $_; } ($xxchi,$xxclo);
      #($mb)= map{ sprintf "%.1f", $_ }  ($atotalmb);

      my $dif= $atotalmb - $mb; # cleaner this way
      my $tlabel="total assembly ";
      if(UPD21AUG and $contamb > 0) {
        $atotalmb -= $contamb; $dif= $atotalmb - $mb;
        $tlabel .= " - contam";
      }
      
      if($RTABLE){
        my $estav= sprintf $mbf, ($esthi+$estlo)/2;
        my @v= ( $asmname, "asmtotal", $atotalmb, $estav+$dif, $estlo+$dif,  $esthi+$dif, $xchi, $covt);  
        print $outh join("\t", @v)."\n";
        #Radd: FlowcytSize=$flowcyto and?  $glnc=  glncformula 
  
        # my($flo,$fhi)= split/\W/,$flowcyto; # check have lonn-hinn or just nn
        my($flo,$fhi)= map{ $_ } $flowcyto =~ m/(\d+)/g;  $fhi=$flo unless($fhi);
        my $estav= sprintf "%.0f", ($fhi+$flo)/2;
        @v= ( $asmname, "flowcyto", $estav, $estav, $flo,  $fhi, 1, 1);  
        print $outh join("\t", @v)."\n";
  
        if($gse_lnc){ #? $glnc formula_lnc >> "$gsf_lnc-$gse_lnc" : $gse_lnc;  
          my $gsa_lnc= sprintf "%.2f", ($gsf_lnc + $gse_lnc)/2;
          @v= ( $asmname, "GLNC_est", $gsa_lnc, $gsa_lnc, $gsf_lnc, $gse_lnc, 1, 1);  
          print $outh join("\t", @v)."\n";
        }
        
      } else {
        my @v= ( "_total", $atotalmb, $estlo+$dif,  ".", $esthi+$dif, ".", ".",$tlabel);  
        if($oneku){ splice(@v,2,2); }
        print $outh join("\t", @v)."\n";
      }
    }
    
  }
  unless($RTABLE){
  print $outh "$dv\n";
  print $outh "  xCopy = excess read copy depth; tCopy = assembly multi-map copy depth\n"; # was (-xCopy)
  # FIXME: check, is this accurate? "tCopy = assembly multi-map copy depth (-xCopy)"?
  #     maybe covt contains xCopy.hi: $covt= $tave/$kulo 

  my($mbrow,$mbcov,$mbgap)= map{ int( $CBIN*$_/100_000)/10 } ($ninrow,$nincov,$nbiggap);#< globals??
  print $outh "  Total span=$mbrow Mb, covspan=$mbcov, gaps=$mbgap Mb for assembly $asmid\n";

  # FIXME: give gse range for both total N and N-nomap reads
  my $pmap= ($n_readid<1)?0:sprintf "%.1f%%",100*$n_maptrd/$n_readid;
  my $rdset=($readsetid)? " for readset $readsetid," : "";
  print $outh "  Genome Size Est= $gsen (Nread), $gsem (Maprd),$rdset\n"; #  formula LN/C,
  print $outh "   for Size=LN/C, Cov=$kulo,$kuhi, N_reads=$n_readid, N_maprd=$n_maptrd,$pmap, L_readlen=$readlen\n";
  if(UCG_KU_SPECIAL) {  
    my($med,$sem,$n,$ave)= map{ $samvals{$asmid}{$_} } qw(kucg_med kucg_sem kucg_n kucg_ave);
    print $outh "  Uniq Conserved Gene Cover: median=$med, ave=$ave, sem=$sem, n=$n\n" if($n);
  }
  
  # pickKUsum == same as debug.log table
  if(my $pkusum= $samvals{$asmid}{'pickKUsum'}) { print $outh "\n$pkusum\n"; }
  
  print $outh "\n";
  }
}

sub meansToSum {
  my($outh,$meanstab, $asmid)=@_;
  my $first=1; # collect @parts from dtable, dont reuse ?
  my @origparts=@parts; @parts=();
  my $dtable= $meanstab;
  # was for my $dtable (@ARGV)  
  my($lasmid,$sn,$stn,$np,$nr) = (0) x 9; %sum=(); 
  warn "#read datatable=$dtable\n" if($debug);
  open(F,$dtable); 
  while(<F>){
    if(/^\W/) { next; }
    my $lin=$_;
    
    if(/uniqrd=/) {
      ($sn)= m/^\s*(\S+)/; 
      ($stn)= ($sn=~s,^(\w+)/,,) ? $1 : ""; #? default
      # NOTE: $sn should == asm ID, match samvals{sn}, but may not? see intitle <> asmid
      my $ataid= $sn;
      
      unless( $samvals{$ataid}) { $ataid= $asmid if($samvals{$asmid}); }
      unless( $samvals{$ataid}) { 
        warn "# missing samvals for entry=$ataid\n";
        map{ $samvals{$ataid}{$_}= 0; } qw(flowcyto cytomb asmtotal atotalmb asmname );
        $samvals{$ataid}{asmname}=$ataid;
      }

      if($lasmid ne $ataid) { #? could get same asmid here from below; yes bad now
        dumpsum($outh, $lasmid,$np,$nr) if($lasmid and $np and %sum);  
        %sum=(); $np=$nr=0; 
        $lasmid=$ataid; 
      }

      #   unless( $samvals{$ataid}) { $ataid= $asmid if($samvals{$asmid}); }
      #   unless( $samvals{$ataid}) { 
      #     warn "# missing samvals for entry=$ataid\n";
      #     map{ $samvals{$ataid}{$_}= 0; } qw(flowcyto cytomb asmtotal atotalmb asmname );
      #     $samvals{$ataid}{asmname}=$ataid;
      #   }
      
      my($urd,$crd,$ann,$nt,$terd,$unkrd)= map{ my($v)= $lin =~ m/$_=(\w+)/ ? $1 : "0"; $v; } 
        qw( uniqrd cdsrc an nt terd unkrd);
      
      unless($stn) {  
        my $st="other"; #? no longer need this?

        # recheck %filts for these, may be wrong .. dont need unless $stn is miss/wrong
        $st="allasm" if($urd eq "0" and $crd eq "0" and not $ann);
        $st="uniqasm" if($urd eq "1" and $crd eq "0" and not $ann);
        $st="dupasm" if($urd eq "2" and $crd eq "0" and not $ann);
        
        $st="cdsasm" if($urd eq "0" and $crd eq "1" and not $ann);
        $st="teasm" if($urd eq "0" and $terd eq "1" and not $ann);
        $st="unkasm" if($urd eq "0" and $unkrd eq "1" and not $ann);

        # FIXME: change ann eq XXX to ann =~ m/XXX/
        $st="CDSann" if($urd eq "0" and $crd eq "0" and $ann eq "CDS"); 
        $st="TEann" if($urd eq "0" and $crd eq "0" and $ann eq "TE");
        $st="RPTann" if($urd eq "0" and $crd eq "0" and $ann eq "repeat");
        $st="NOann" if($urd eq "0" and $crd eq "0" and $ann eq 'no(CDS|TE|repeat)');
        $st="CONTAM" if($ann eq "contam");# UPD21AUG
        $stn=$st unless($stn);
        
        # if($stn and $st eq "other"){ $st=$stn; }
        # if($stn ne $st) { } #what ?
      }
      
      push @parts, $stn if($first==1); # was BAD.. already set this global
      $np++; $nr=0;
      
    } elsif(/^Item/) { 
      # @stath= split; 
       
    } elsif(/^\s*(rd|aCov|Cov)/) {  # this is bad now, Item names are "^\s*(aCovT|aCovM|aCovU)"
      ## need new $keeps= 'rdacovt|rdacovm|aCovT|aCovM'; 
      my @statv= split;
      my($sti)= $statv[0];
      if($sti =~ m/^aCov|^Cov/){ $sti='rd'. lc($sti); $statv[0]=$sti; } # aCovM => rdacovm
      if($sti =~ m/^($keeps)$/) {
        $sum{$stn}{$sti}= [@statv];  # this way?
        $nr++;
      }
    }
    
  } close(F);
  
  # warn "#done datatable=$dtable\n" if($debug);
  dumpsum($outh, $lasmid,$np,$nr) if($lasmid and $np and %sum); 

# move to meansToSum()
# my($mbrow,$mbcov,$mbgap)= map{ int( $CBIN*$_/100_000)/10 } ($ninrow,$nincov,$nbiggap);
# print $outsumh "# total span=$mbrow.mb, covspan=$mbcov.mb, gaps=$mbgap.mb with n_readid=$n_readid for $asmid\n";

  @parts= @origparts;
}

sub openRead {
  my($fn)=@_; my($ok,$inh)=(0); 
  return(0,undef) unless($fn and -f $fn);
  if($fn =~ /\.gz/) { $ok= open($inh,"gunzip -c $fn |"); } else { $ok= open($inh,$fn); }
  # die "cant read $fn" unless ($ok); # warn? leave to caller
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

=item UPD21JUN try2

 $evigene/scripts/genoasm/gnodes_covsum.pl  -asmid daphplx_gasm16ml -title daphplx_gasm16ml8dtest -anntab daphplx_gasm16ml.anntab -crclass ../daphplx17evgt1m_cds.idclass -sumdata ../daphplx20chrs.metad -debug -genexcopy daphplx17evgt1m_cds_SRR13333791_b8_mim.genexcopy  daphplx_gasm16ml_SRR13333791_b8_mim.cc8a.covtab
#covsum output to sum=daphplx_gasm16ml8dtest_SRR13333791_sum.txt means=daphplx_gasm16ml8dtest_SRR13333791.means
#read sampledata=../daphplx20chrs.metad
#read anntable=daphplx_gasm16ml.anntab
# read nid=33037, nclass=3 from ../daphplx17evgt1m_cds.idclass
#title allasm/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 
#covtab span=156.4.mb, covspan=140.1.mb, gaps=13.6.mb info: n_readid=58112591, n_nomap=7586800, n_mapok=64759246, n_mapbad=224995419, n_dupbad=123328715, 
#meanvar nt=1401617, ti=allasm/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=1082835, ti=uniqasm/daphplx_gasm16ml8dtest uniqrd=1, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=318782, ti=dupasm/daphplx_gasm16ml8dtest uniqrd=2, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=597212, ti=cdsasm/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=1, an=0, terd=0, unkrd=0 
#meanvar nt=158, ti=teasm/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=0, terd=1, unkrd=0 
#meanvar nt=0, ti=unkasm/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=2, an=0, terd=2, unkrd=1 
#meanvar nt=531305, ti=CDSann/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=CDS, terd=0, unkrd=0 
#meanvar nt=18943, ti=CDSbus/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=busco, terd=0, unkrd=0 
#meanvar nt=64655, ti=TEann/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=TE, terd=0, unkrd=0 
#meanvar nt=62619, ti=RPTann/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=repeat, terd=0, unkrd=0 
#meanvar nt=790670, ti=NOann/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=no(CDS|TE|repeat|gap), terd=0, unkrd=0 
#meanvar nt=859, ti=UCGmed/daphplx_gasm16ml8dtest uniqrd=1, cdsrd=1
#read datatable=daphplx_gasm16ml8dtest_SRR13333791.means
#pickKuni Genes gucg_med=38.4, sd(gucg_med)=1.26, gucg_ave=38.0, kucg_n=954
#pickKuni samvals kulo=38.4, kuhi=38.4
#pickKuni uniqasm kulo=41, kuhi=42.98
#pickKuni CDSbus kulo=44, kuhi=44.24
#pickKuni bestpart: CDSbus
#Item   Cmed    Cave    Cn      Cerr
# CDSbus        44      44.24   18941   0.0017
# uniqasm       41      42.98   1082835 0.0019
# allasm        40      48.77   1400237 0.0074
# NOann         39      45.89   789972  0.0078
# CDSann        42      50.72   530835  0.0120
# CMgenes       38.4    38.0    954     0.0130
#pickKuni UCG ku_med=44, sd(ku_med)=44.868, ku_ave=43.993, sem=1.531, kucg_n=859

>>> WRONG KUlo/hi=44 when CMgenes/UCG exists.. change one (both?)
head -33 daphplx_gasm16ml8dtest_SRR13333791_sum.txt
Source=Daphplx16ml, KUlow=44, KUhigh=44, FlowcytSize=215-264 Mb Formula_LN/C=172.2 Mb (daphplx_gasm16ml)
_____   ______                  Total    
Item_   Obs.Mb  Est.Mb  xCopy   tCopy   Description
---------------------------------------------------------------------------
allasm  140     155     1.11    1.42    measured assembly
_total  156     171     .       .       total assembly 
uniqasm 108     106     0.98    0.98    asm with unique gDNA
dupasm  32      49      1.56    2.93    asm with multimap gDNA
cdsasm  60      76      1.28    1.87    asm with CDS-mapped gDNA
teasm   0.0     0.0     0.77    0.86    asm with TE-mapped gDNA
CDSann  53      61      1.15    1.69    asm with CDS annotations
CDSbus  1.9     1.9     1.01    1.03    asm with unique BUSCO orlog CDS
TEann   6.5     9.3     1.44    1.86    asm with Transposons
RPTann  6.3     10.4    1.65    2.91    asm with simple Repeats
NOann   79      82      1.04    1.19    asm without annotations
UCGmed  0.1     0.1     1.00    1.00    
---------------------------------------------------------------------------
  xCopy = excess read copy depth; tCopy = assembly multi-map copy depth
  Total span=156.4 Mb, covspan=140.1, gaps=13.6 Mb for assembly daphplx_gasm16ml
  Genome Size Est= 198.1 Mb (Nread), 172.2 Mb (Maprd), for readset SRR13333791,
   for Size=LN/C, Cov=44,44, N_reads=58112591, N_maprd=50525791,86.9%, L_readlen=150
  Uniq Conserved Gene Cover: median=38.4, ave=38.0, sem=1.26, n=954
  
  
=item UPD21JUN try1

dgmm:gnode8d:% $evigene/scripts/genoasm/gnodes_covsum.pl  -asmid daphplx_gasm16ml -title daphplx_gasm16ml8dtest -anntab daphplx_gasm16ml.anntab -crclass ../daphplx17evgt1m_cds.idclass -sumdata ../daphplx20chrs.metad -debug -genexcopy daphplx17evgt1m_cds_SRR13333791_b8_mim.genexcopy  daphplx_gasm16ml_SRR13333791_b8_mim.cc8a.covtab
#covsum output to sum=daphplx_gasm16ml8dtest_SRR13333791_sum.txt means=daphplx_gasm16ml8dtest_SRR13333791.means
#read sampledata=../daphplx20chrs.metad
#read anntable=daphplx_gasm16ml.anntab << WARN if missing
# read nid=33037, nclass=3 from ../daphplx17evgt1m_cds.idclass
#title allasm/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 
#covtab span=156.4.mb, covspan=140.1.mb, gaps=0.mb info: n_readid=58112591, n_nomap=7586800, n_mapok=64759246, n_mapbad=224995419, n_dupbad=123328715, 

#meanvar nt=1401617, ti=allasm/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=1082835, ti=uniqasm/daphplx_gasm16ml8dtest uniqrd=1, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=318782, ti=dupasm/daphplx_gasm16ml8dtest uniqrd=2, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=597212, ti=cdsasm/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=1, an=0, terd=0, unkrd=0 
#meanvar nt=158, ti=teasm/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=0, terd=1, unkrd=0 

     vvv missing annot data ???  anntab.gz not used
#meanvar nt=0, ti=unkasm/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=2, an=0, terd=2, unkrd=1 
#meanvar nt=0, ti=CDSann/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=CDS, terd=0, unkrd=0 
#meanvar nt=0, ti=CDSbus/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=busco, terd=0, unkrd=0 
#meanvar nt=0, ti=TEann/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=TE, terd=0, unkrd=0 
#meanvar nt=0, ti=RPTann/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=repeat, terd=0, unkrd=0 

#meanvar nt=1401617, ti=NOann/daphplx_gasm16ml8dtest uniqrd=0, cdsrd=0, an=no(CDS|TE|repeat|gap), terd=0, unkrd=0 
#meanvar nt=0, ti=UCGmed/daphplx_gasm16ml8dtest uniqrd=1, cdsrd=1

#read datatable=daphplx_gasm16ml8dtest_SRR13333791.means
#pickKuni Genes gucg_med=38.4, sd(gucg_med)=1.26, gucg_ave=38.0, kucg_n=954
#pickKuni samvals kulo=38.4, kuhi=38.4
#pickKuni uniqasm kulo=41, kuhi=42.98
#pickKuni bestpart: uniqasm
#Item           Cmed    Cave    Cn      Cerr
# uniqasm       41      42.98   1082835 0.0019
# NOann         40      48.77   1400237 0.0074
# allasm        40      48.77   1400237 0.0074
# CMgenes       38.4    38.0    954     0.0130  << FIXME Cn not same here, Cerr too high 
  vvv missing annot data
#pickKuni UCG ku_med=0, sd(ku_med)=0.000, ku_ave=0.000, sem=0, kucg_n=

=cut

=item data tables

  allasm/daphplx_gasm16ml uniqrd=0, cdsrd=0, an=0,  stats for nt=1382381
  Item  	Median	Mean	SEM	Nitem	StDev	Sum
  rdacovt	  16	33.27	0.86	1382381	1016.36	45986935
  rdacovm	  15	18.78	0.17	1382330	205.71	25964243
  rdacovu	  14	16.41	0.12	1259948	136.44	20678633
  
  CDSbus/daphplx_gasm16ml uniqrd=0, cdsrd=0, an=busco,  stats for nt=16796
  Item  	Median	Mean	SEM	Nitem	StDev	Sum
  rdacovt	  17	18.09	0.15	16796	19.99	303786
  rdacovm	  17	17.29	0.14	16796	18.51	290451
  rdacovu	  17	17.14	0.14	16319	18.28	279781

=cut

=item try1b

../daphnia3covtables.pl ../daphnia3covsum.data dmag19skasm.cov3h.meansum
#read sampledata=../daphnia3covsum.data
#read datatable=dmag19skasm.cov3h.meansum
#done datatable=dmag19skasm.cov3h.meansum

Source=Daphmag19sk, KUlow=40, KUhigh=42.13, FlowcytSize=234-391 Mb  (dmag19skasm)
_____   ______  ____Low____     ____High____    Total
Item_   Obs.Mb  Est.Mb  xCopy   Est.Mb  xCopy   tCopy
-------------------------------------------------------
allasm  113.6   193.6   1.70    203.9   1.80    7.72
totalasm        123     203     1.70    213.3   1.80    7.72
uniqasm 82.9    82.2    0.99    86.5    1.04    1.05
dupasm  30.7    111.4   3.63    117.3   3.82    25.75
cdsasm  57.9    139.9   2.41    147.3   2.54    13.93
CDSann  34.6    39.2    1.13    41.3    1.19    1.46
CDSbus  1.9     1.9     1.00    2.0     1.05    1.07
CDSuni  25.0    25.6    1.02    27.0    1.08    1.14
CDSdup  4.1     5.3     1.31    5.6     1.38    2.04
TEann   1.3     2.1     1.58    2.2     1.66    3.23
RPTann  5.8     5.9     1.03    6.2     1.08    1.70
NOann   73.5    148.1   2.02    156.0   2.12    11.10
-------------------------------------------------------

=cut

=item metadata

-- use this, w/ updates
  metadata=$gw/daphnia/daphcov3gsum/cov6gset/daphnia3cov6sum.data

pt=dmag15nwb2fullasm
flowcyto=234-391 Mb
glncformula=210 Mb  
asmtotal=134 Mb
asmname=Daphmag10nb
kulo=55   # modal peak of all asm.SRR7825549cl_1a_cc6e.covtab 
kuhi=47   # cds.busco rdepth
.............................

pt=dmag14bgi2vtop5k
flowcyto=234-391 Mb
glncformula=210 Mb  
asmtotal=180 Mb
asmname=Daphmag14il
kulo=55   # modal peak of all asm.SRR7825549cl_1a_cc6e.covtab 
kuhi=47   # cds.busco rdepth
.............................

pt=dmag19skasm
flowcyto=234-391 Mb
glncformula=210 Mb  
asmtotal=123 Mb
asmname=Daphmag19sk
kulo=55   # modal peak of all asm.SRR7825549cl_1a_cc6e.covtab 
kuhi=47   # cds.busco rdepth
.............................

pt=dropse20chrs
flowcyto=161-180 Mb
glncformula=174 Mb
asmtotal=163 Mb
asmname=Dropse20uc
# data dependent Ku vals, from various results..
kulo=91 # dropse20t1cds_SRR11813283_1 busco.covtab ave Covr
kuhi=95 # dropse20chrs  CDSbusco rdacovm ave
............................................


=item gnodes_covsum daphmag15asm trial

$gw/daphnia/gnodes_covsum.pl -debug -title dmag15nwb2asm_cov7b  \
 dmag15nwb2asm_SRR7825549cl_1a_bwa.cdschr7b.covtab

# output to sum=dmag15nwb2asm_cov7b_sum.txt means=dmag15nwb2asm_cov7b.means
#title allasm/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 
#covtab span=126.1.mb, covspan=107.3.mb info: n_readid=4450506, n_mapok=79875595, n_mapbad=130272508, n_dupbad=236433975, 
#meanvar nt=1073743, ti=allasm/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=857529, ti=uniqasm/dmag15nwb2asm_cov7b uniqrd=1, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=216214, ti=dupasm/dmag15nwb2asm_cov7b uniqrd=2, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=459933, ti=cdsasm/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=1, an=0, terd=0, unkrd=0 
#meanvar nt=86982, ti=teasm/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=0, an=0, terd=1, unkrd=0 
#meanvar nt=108431, ti=unkasm/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=1 
#meanvar nt=0, ti=CDSann/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=0, an=CDS, terd=0, unkrd=0 
#meanvar nt=0, ti=CDSbus/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=0, an=busco, terd=0, unkrd=0 
#meanvar nt=0, ti=TEann/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=0, an=TE, terd=0, unkrd=0 
#meanvar nt=0, ti=RPTann/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=0, an=repeat, terd=0, unkrd=0 
#meanvar nt=1073743, ti=NOann/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=0, an=no(CDS|TE|repeat), terd=0, unkrd=0 
#read datatable=dmag15nwb2asm_cov7b.means
# missing samvals for entry=dmag15nwb2asm_cov7b

==> dmag15nwb2asm_cov7b_sum.txt <==
Source=dmag15nwb2asm_cov7b, KUlow=45, KUhigh=48.85, FlowcytSize=0 Formula_LN/C=0 (dmag15nwb2asm_cov7b)
_____   ______  Low     Low     High    High    Total    
Item_   Obs.Mb  Est.Mb  xCopy   Est.Mb  xCopy   tCopy   Description
---------------------------------------------------------------------------
allasm  107     161     1.50    175     1.63    4.00    measured assembly
uniqasm 86      86      1.00    93      1.09    1.09    asm with unique gDNA
dupasm  21      75      3.51    81      3.81    15.56   asm with multimap gDNA
cdsasm  46      62      1.36    68      1.47    2.24    asm with CDS-mapped gDNA
teasm   8.7     21.5    2.49    23.4    2.70    6.07    asm with TE-mapped gDNA
unkasm  10.7    57.5    5.39    62.5    5.85    27.39   asm with UnclassRepeats-mapped gDNA
NOann   107     161     1.50    175     1.63    4.00    asm without annotations
CDSbus  0.0     0.0     0.00    0.0     0.00    0.00    asm with unique BUSCO orlog CDS
---------------------------------------------------------------------------
  xCopy = excess read copy depth; tCopy = assembly multi-map copy depth (-xCopy)

==> dmag15nwb2asm_cov7b.means <==
allasm/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0  stats for nt=1073743
Item  	Median	Mean	SEM	Nitem	StDev	Mode	Sum
aCovT 	  46	180.03	3.59	1073743	3721.34	  45	193302323
aCovM 	  45	73.30	0.64	1071430	664.71	  45	78538839
aCovU 	  44	51.97	0.23	1037800	233.33	  45	53934516

uniqasm/dmag15nwb2asm_cov7b uniqrd=1, cdsrd=0, an=0, terd=0, unkrd=0  stats for nt=857529
Item  	Median	Mean	SEM	Nitem	StDev	Mode	Sum
aCovT 	  45	48.87	0.16	857529	149.65	  45	41910791
aCovM 	  45	48.85	0.16	857529	149.65	  45	41892932
aCovU 	  45	48.84	0.16	857529	149.64	  45	41882182

dupasm/dmag15nwb2asm_cov7b uniqrd=2, cdsrd=0, an=0, terd=0, unkrd=0  stats for nt=216214
Item  	Median	Mean	SEM	Nitem	StDev	Mode	Sum
aCovT 	  76	700.19	17.82	216214	8287.59	  48	151391532
aCovM 	  52	171.32	3.15	213901	1457.19	  48	36645907
aCovU 	  40	66.86	1.07	180271	454.86	  44	12052334

cdsasm/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=1, an=0, terd=0, unkrd=0  stats for nt=459933
Item  	Median	Mean	SEM	Nitem	StDev	Mode	Sum
aCovT 	  47	100.61	1.30	459933	880.64	  46	46272791
aCovM 	  47	66.36	0.59	459460	396.96	  46	30491472
aCovU 	  46	55.14	0.44	446618	291.62	  46	24628486

teasm/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=0, an=0, terd=1, unkrd=0  stats for nt=86982
Item  	Median	Mean	SEM	Nitem	StDev	Mode	Sum
aCovT 	  71	273.14	5.19	86982	1529.73	  46	23758559
aCovM 	  53	121.44	1.83	86548	536.99	  45	10510473
aCovU 	  44	75.66	1.35	79098	378.96	  46	5984298

unkasm/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=1  stats for nt=108431
Item  	Median	Mean	SEM	Nitem	StDev	Mode	Sum
aCovT 	  74	1232.44	35.49	108431	11688.08	  47	133634645
aCovM 	  50	263.33	6.20	106755	2025.86	  45	28111510
aCovU 	  41	81.83	1.75	95963	541.05	  43	7852715

NOann/dmag15nwb2asm_cov7b uniqrd=0, cdsrd=0, an=no(CDS|TE|repeat), terd=0, unkrd=0  stats for nt=1073743
Item  	Median	Mean	SEM	Nitem	StDev	Mode	Sum
aCovT 	  46	180.03	3.59	1073743	3721.34	  45	193302323
aCovM 	  45	73.30	0.64	1071430	664.71	  45	78538839
aCovU 	  44	51.97	0.23	1037800	233.33	  45	53934516


=item gnodes_covsum trials

  
-- cds, te seq, limited info

$gw/daphnia/gnodes_covsum.pl -debug -title dmag7fincds7b  dmag7fincds_SRR7825549cl_1a_bwa.cdschr7b.covtab
# output to sum=dmag7fincds7b_sum.txt means=dmag7fincds7b.means
#title allasm/dmag7fincds7b uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 
#covtab span=35.3.mb, covspan=32.6.mb info: n_readid=5198220, n_mapok=15668803, n_mapbad=65960367, n_dupbad=4072502, 
#meanvar nt=326382, ti=allasm/dmag7fincds7b uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=229734, ti=uniqasm/dmag7fincds7b uniqrd=1, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=96648, ti=dupasm/dmag7fincds7b uniqrd=2, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=0, ti=cdsasm/dmag7fincds7b uniqrd=0, cdsrd=1, an=0, terd=0, unkrd=0 
#meanvar nt=0, ti=teasm/dmag7fincds7b uniqrd=0, cdsrd=0, an=0, terd=1, unkrd=0 
#meanvar nt=0, ti=unkasm/dmag7fincds7b uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=1 
#meanvar nt=0, ti=CDSann/dmag7fincds7b uniqrd=0, cdsrd=0, an=CDS, terd=0, unkrd=0 
#meanvar nt=0, ti=CDSbus/dmag7fincds7b uniqrd=0, cdsrd=0, an=busco, terd=0, unkrd=0 
#meanvar nt=0, ti=TEann/dmag7fincds7b uniqrd=0, cdsrd=0, an=TE, terd=0, unkrd=0 
#meanvar nt=0, ti=RPTann/dmag7fincds7b uniqrd=0, cdsrd=0, an=repeat, terd=0, unkrd=0 
#meanvar nt=326382, ti=NOann/dmag7fincds7b uniqrd=0, cdsrd=0, an=no(CDS|TE|repeat), terd=0, unkrd=0 
#read datatable=dmag7fincds7b.means
# missing samvals for entry=dmag7fincds7b

dmag20sk4ma_tefam7b_sum.txt
* note HUGE KU vals, should use samp data w/ KU set from CDS vals for better TE Est.
Source=dmag20sk4ma_tefam7b, KUlow=326, KUhigh=768.60, FlowcytSize=0 Formula_LN/C=0 (dmag20sk4ma_tefam7b)
_____   ______  Low     Low     High    High    Total    
Item_   Obs.Mb  Est.Mb  xCopy   Est.Mb  xCopy   tCopy   Description
---------------------------------------------------------------------------
allasm  1.4     3.6     2.51    8.4     5.91    6.99    measured assembly
uniqasm 1.0     1.0     1.00    2.4     2.36    2.36    asm with unique gDNA
dupasm  0.4     2.5     6.41    6.0     15.11   18.98   asm with multimap gDNA
CDSbus  0.0     0.0     0.00    0.0     0.00    0.00    asm with unique BUSCO orlog CDS
NOann   1.4     3.6     2.51    8.4     5.91    6.99    asm without annotations
---------------------------------------------------------------------------
  xCopy = excess read copy depth; tCopy = assembly multi-map copy depth (-xCopy)

dmag7fincds7b_sum.txt
Source=dmag7fincds7b, KUlow=45, KUhigh=57.38, FlowcytSize=0 Formula_LN/C=0 (dmag7fincds7b)
_____   ______  Low     Low     High    High    Total    
Item_   Obs.Mb  Est.Mb  xCopy   Est.Mb  xCopy   tCopy   Description
---------------------------------------------------------------------------
allasm  33      41      1.27    53      1.62    2.18    measured assembly
uniqasm 23      23      1.00    29      1.28    1.28    asm with unique gDNA
dupasm  9.6     18.3    1.91    23.3    2.43    4.34    asm with multimap gDNA
CDSbus  0.0     0.0     0.00    0.0     0.00    0.00    asm with unique BUSCO orlog CDS
NOann   33      41      1.27    53      1.62    2.18    asm without annotations
---------------------------------------------------------------------------
  xCopy = excess read copy depth; tCopy = assembly multi-map copy depth (-xCopy)

head -24 dmag7fincds7b.means  
allasm/dmag7fincds7b uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0  stats for nt=326382
Item  	Median	Mean	SEM	Nitem	StDev	Mode	Sum
aCovT 	  47	98.28	1.71	326382	976.98	  49	32076692
aCovM 	  42	72.69	1.52	325644	869.74	  46	23671352
aCovU 	  45	70.39	1.72	268535	891.05	  46	18901235

uniqasm/dmag7fincds7b uniqrd=1, cdsrd=0, an=0, terd=0, unkrd=0  stats for nt=229734
Item  	Median	Mean	SEM	Nitem	StDev	Mode	Sum
aCovT 	  45	57.39	0.96	229734	458.19	  46	13185216
aCovM 	  45	57.38	0.96	229734	458.19	  46	13181553
aCovU 	  45	57.36	0.96	229734	458.19	  46	13178390

dupasm/dmag7fincds7b uniqrd=2, cdsrd=0, an=0, terd=0, unkrd=0  stats for nt=96648
Item  	Median	Mean	SEM	Nitem	StDev	Mode	Sum
aCovT 	  52	195.47	5.31	96648	1650.55	  49	18891476
aCovM 	  27	109.37	4.64	95910	1437.19	  24	10489799
aCovU 	  36	147.49	10.47	38801	2062.04	   5	5722845

NOann/dmag7fincds7b uniqrd=0, cdsrd=0, an=no(CDS|TE|repeat), terd=0, unkrd=0  stats for nt=326382
Item  	Median	Mean	SEM	Nitem	StDev	Mode	Sum
aCovT 	  47	98.28	1.71	326382	976.98	  49	32076692
aCovM 	  42	72.69	1.52	325644	869.74	  46	23671352
aCovU 	  45	70.39	1.72	268535	891.05	  46	18901235

=item try4 fixme

pt=dmag19skasm; 
$evigene/scripts/genoasm/gnodes_covsum.pl -debug -title ${pt}_cov7ban -anntable ${pt}.anntab -sumdata ../daphnia3cov6sum.data  ${pt}_SRR7825549cl_1a_bwa.cdschr7b.covtab

BUG: -title dmag19skasm_cov7ban is used as asm ID to find sumdata, but -sumdata uses other, pt=dmag19skasm
...  change -title to same asm ID? or add new opt -asmid $pt, use title for output names

dmag19skasm_cov7ban.means
allasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0  stats for nt=1149351
Item  	Median	Mean	SEM	Nitem	StDev	Mode	Sum
aCovT 	  55	129.31	8.70	1149351	9324.29	  55	148620870
aCovM 	  53	81.52	3.98	1149182	4268.80	  55	93685022
aCovU 	  52	61.94	1.62	1124683	1716.88	  55	69657402

uniqasm/dmag19skasm_cov7ban uniqrd=1, cdsrd=0, an=0, terd=0, unkrd=0  stats for nt=883110
Item  	Median	Mean	SEM	Nitem	StDev	Mode	Sum
aCovT 	  54	56.50	0.16	883110	153.78	  55	49892053
aCovM 	  54	56.46	0.16	883110	153.77	  55	49864595
dgmm:gdnareads:% head dmag19skasm_cov7ban_sum.txt
Source=dmag19skasm_cov7ban, KUlow=54, KUhigh=56.46, FlowcytSize=0 Formula_LN/C=0 (dmag19skasm_cov7ban)

=item try5 fixed asmid, anntab

** bad allasm total mb
pt=dmag19skasm; 
$evigene/scripts/genoasm/gnodes_covsum.pl -debug -asmid $pt  -title ${pt}_cov7ban \
 -anntable ${pt}.anntab -sumdata ../daphnia3cov6sum.data  ${pt}_SRR7825549cl_1a_bwa.cdschr7b.covtab

dmag19skasm_cov7ban.means
allasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0  stats for nt=614030 << was nt=1149351
Item  	Median	Mean	SEM	Nitem	StDev	Mode	Sum
aCovT 	  53	64.49	0.67	614030	522.21	  53	39601293
aCovM 	  52	55.53	0.30	613953	233.17	  53	34091461
aCovU 	  51	51.54	0.18	602286	141.67	  54	31042082

# output to sum=dmag19skasm_cov7ban_sum.txt means=dmag19skasm_cov7ban.means
#read sampledata=../daphnia3cov6sum.data
#read anntable=dmag19skasm.anntab
#title allasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 
#covtab span=123.1.mb, covspan=114.9.mb info: n_readid=5398296, n_mapok=61130093, n_mapbad=85997706, n_dupbad=87731998, 
                                                 ^^^^^^^^^^^^ 1 pt only not total; need to add like others
#meanvar nt=614030, ti=allasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 
            ^^^^^^ should == covspan 1149000
#meanvar nt=496959, ti=uniqasm/dmag19skasm_cov7ban uniqrd=1, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=117071, ti=dupasm/dmag19skasm_cov7ban uniqrd=2, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=53675, ti=cdsasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=1, an=0, terd=0, unkrd=0 
#meanvar nt=7379, ti=teasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=0, terd=1, unkrd=0 
#meanvar nt=38777, ti=unkasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=1 
#meanvar nt=516283, ti=CDSann/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=CDS, terd=0, unkrd=0 
#meanvar nt=18773, ti=CDSbus/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=busco, terd=0, unkrd=0 
#meanvar nt=100404, ti=TEann/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=TE, terd=0, unkrd=0 
#meanvar nt=3468, ti=RPTann/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=repeat, terd=0, unkrd=0 
#meanvar nt=614029, ti=NOann/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=no(CDS|TE|repeat), terd=0, unkrd=0 

#read datatable=dmag19skasm_cov7ban.means

pt=dmag19skasm
flowcyto=234-391 Mb
asmtotal=123 Mb
kulo=55   # modal peak of all asm.SRR7825549cl_1a_cc6e.covtab 
kuhi=47   # cds.busco rdepth

CDSbus/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=busco, terd=0, unkrd=0  stats for nt=18773
Item    Median  Mean    SEM     Nitem   StDev   Mode    Sum
aCovT     55    56.73   0.45    18773   61.11     54    1065008
aCovM     55    55.34   0.43    18773   59.12     54    1038867
aCovU     55    55.00   0.43    18458   58.81     54    1015237

>> poor value KUhigh=56.73 ? should be CDSbus median or/and asm.data vals
>> wrong asmid here maybe?
Source=dmag19skasm_cov7ban, KUlow=53, KUhigh=56.73, FlowcytSize=0 Formula_LN/C=0 (dmag19skasm_cov7ban)
_____   ______  Low     Low     High    High    Total    
Item_   Obs.Mb  Est.Mb  xCopy   Est.Mb  xCopy   tCopy   Description
---------------------------------------------------------------------------
allasm  61      60      0.98    64      1.05    1.22    measured assembly << ?? way too low
uniqasm 50      47      0.95    50      1.01    1.01    asm with unique gDNA
dupasm  11.7    13.1    1.12    14.0    1.19    2.08    asm with multimap gDNA
cdsasm  5.4     5.9     1.10    6.3     1.18    1.44    asm with CDS-mapped gDNA
teasm   0.7     1.2     1.57    1.2     1.68    2.28    asm with TE-mapped gDNA
unkasm  3.9     7.0     1.81    7.5     1.94    3.59    asm with UnclassRepeats-mapped gDNA
CDSann  52      99      1.93    106     2.06    3.79    asm with CDS annotations
CDSbus  1.9     1.8     0.98    2.0     1.04    1.07    asm with unique BUSCO orlog CDS
TEann   10.0    55.3    5.51    59.2    5.90    13.92   asm with Transposons
RPTann  0.3     2.8     8.03    3.0     8.60    15.44   asm with simple Repeats
NOann   61      60      0.98    64      1.05    1.22    asm without annotations
---------------------------------------------------------------------------
  xCopy = excess read copy depth; tCopy = assembly multi-map copy depth (-xCopy)

=item try6 fixed fixes
  .. n_readid is 1 pt, need to sum parts
  
  .. sum is bad now, 1 for each part:
       131  dmag19skasm_cov7ban_sum.txt

Source=Daphmag19sk, KUlow=47, KUhigh=55, FlowcytSize=234-391 Mb Formula_LN/C=210 Mb   (dmag19skasm)
_____	______	Low	Low	High	High	Total	 
Item_	Obs.Mb	Est.Mb	xCopy	Est.Mb	xCopy	tCopy	Description
---------------------------------------------------------------------------
allasm	115	170	1.48	199	1.73	2.75	measured assembly
_total	123	178	.	207	.	.	total assembly 
uniqasm	0.0	0.0	0.00	0.0	0.00	0.00	asm with unique gDNA
CDSbus	0.0	0.0	0.00	0.0	0.00	0.00	asm with unique BUSCO orlog CDS
---------------------------------------------------------------------------

Source=Daphmag19sk, KUlow=47, KUhigh=55, FlowcytSize=234-391 Mb Formula_LN/C=210 Mb   (dmag19skasm)
_____	______	Low	Low	High	High	Total	 
Item_	Obs.Mb	Est.Mb	xCopy	Est.Mb	xCopy	tCopy	Description
---------------------------------------------------------------------------
uniqasm	88	91	1.03	106	1.20	1.20	asm with unique gDNA
CDSbus	0.0	0.0	0.00	0.0	0.00	0.00	asm with unique BUSCO orlog CDS
uniqasm	88	91	1.03	106	1.20	1.20	asm with unique gDNA

  ... etc for each pt

pt=dmag19skasm; $evigene/scripts/genoasm/gnodes_covsum.pl -debug -asmid $pt  -title ${pt}_cov7ban -anntable ${pt}.anntab -sumdata ../daphnia3cov6sum.data  ${pt}_SRR7825549cl_1a_bwa.cdschr7b.covtab
# output to sum=dmag19skasm_cov7ban_sum.txt means=dmag19skasm_cov7ban.means
#read sampledata=../daphnia3cov6sum.data
#read anntable=dmag19skasm.anntab
#title allasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 
#covtab span=123.1.mb, covspan=114.9.mb info: n_readid=5398296, n_mapok=61130093, n_mapbad=85997706, n_dupbad=87731998, 
#meanvar nt=1149351, ti=allasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=883110, ti=uniqasm/dmag19skasm_cov7ban uniqrd=1, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=266241, ti=dupasm/dmag19skasm_cov7ban uniqrd=2, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=494633, ti=cdsasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=1, an=0, terd=0, unkrd=0 
#meanvar nt=98026, ti=teasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=0, terd=1, unkrd=0 
#meanvar nt=107156, ti=unkasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=1 
#meanvar nt=516283, ti=CDSann/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=CDS, terd=0, unkrd=0 
#meanvar nt=18773, ti=CDSbus/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=busco, terd=0, unkrd=0 
#meanvar nt=100404, ti=TEann/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=TE, terd=0, unkrd=0 
#meanvar nt=3468, ti=RPTann/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=repeat, terd=0, unkrd=0 
#meanvar nt=614029, ti=NOann/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=no(CDS|TE|repeat), terd=0, unkrd=0 
#read datatable=dmag19skasm_cov7ban.means

=item try7

 pt=dmag19skasm;
 $evigene/scripts/genoasm/gnodes_covsum.pl -debug -asmid $pt  -title ${pt}_cov7ban \
  -anntable ${pt}.anntab -sumdata ../daphnia3cov6sum.data  ${pt}_SRR7825549cl_1a_bwa.cdschr7b.covtab
  
# output to sum=dmag19skasm_cov7ban_sum.txt means=dmag19skasm_cov7ban.means
#read sampledata=../daphnia3cov6sum.data
#read anntable=dmag19skasm.anntab
#title allasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 

#covtab span=123.1.mb, covspan=114.9.mb info: n_readid=43186373, n_mapok=61130093, n_mapbad=85997706, n_dupbad=87731998, 

#meanvar nt=1149351, ti=allasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=883110, ti=uniqasm/dmag19skasm_cov7ban uniqrd=1, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=266241, ti=dupasm/dmag19skasm_cov7ban uniqrd=2, cdsrd=0, an=0, terd=0, unkrd=0 
#meanvar nt=494633, ti=cdsasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=1, an=0, terd=0, unkrd=0 
#meanvar nt=98026, ti=teasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=0, terd=1, unkrd=0 
#meanvar nt=107156, ti=unkasm/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=0, terd=0, unkrd=1 
#meanvar nt=516283, ti=CDSann/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=CDS, terd=0, unkrd=0 
#meanvar nt=18773, ti=CDSbus/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=busco, terd=0, unkrd=0 
#meanvar nt=100404, ti=TEann/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=TE, terd=0, unkrd=0 
#meanvar nt=3468, ti=RPTann/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=repeat, terd=0, unkrd=0 
#meanvar nt=614029, ti=NOann/dmag19skasm_cov7ban uniqrd=0, cdsrd=0, an=no(CDS|TE|repeat), terd=0, unkrd=0 
#read datatable=dmag19skasm_cov7ban.means

#pickKuni samvals kulo=47, kuhi=55    << 47-48 is still val for busco cds x gdna
#pickKuni uniqasm kulo=54, kuhi=56.46
#pickKuni CDSbus kulo=55, kuhi=56.73  << always use klo == CDSbus median as one KU value?

dmag19skasm_cov7ban_sum.txt
Source=Daphmag19sk, KUlow=47, KUhigh=55, FlowcytSize=234-391 Mb Formula_LN/C=210 Mb   (dmag19skasm)
_____   ______  Low     Low     High    High    Total    
Item_   Obs.Mb  Est.Mb  xCopy   Est.Mb  xCopy   tCopy   Description
---------------------------------------------------------------------------
allasm  115     170     1.48    199     1.73    2.75    measured assembly
_total  123     178     .       207     .       .       total assembly 
uniqasm 88      91      1.03    106     1.20    1.20    asm with unique gDNA
dupasm  27      80      2.99    93      3.50    7.89    asm with multimap gDNA
cdsasm  49      66      1.34    78      1.57    2.24    asm with CDS-mapped gDNA
teasm   9.8     25.3    2.58    29.6    3.02    6.10    asm with TE-mapped gDNA
unkasm  10.7    59.4    5.54    69.5    6.49    15.09   asm with UnclassRepeats-mapped gDNA
CDSann  52      103     1.99    120     2.32    4.28    asm with CDS annotations
CDSbus  1.9     1.9     1.01    2.2     1.18    1.21    asm with unique BUSCO orlog CDS
TEann   10.0    57.0    5.68    66.7    6.65    15.69   asm with Transposons
RPTann  0.3     2.9     8.29    3.4     9.70    17.41   asm with simple Repeats
NOann   61      62      1.01    73      1.18    1.37    asm without annotations
---------------------------------------------------------------------------
  xCopy = excess read copy depth; tCopy = assembly multi-map copy depth (-xCopy)

  ^^^ === 7b changes versus cov6g set === vvv
  .. most are ~same values, obs and est
  .. cdsasm obs smaller, due to TE annot of CDS changed,
      -- and cdsasm est xCopy lower, still highish (40% of tot instead of 60%)
  .. unkasm added : major component, adjust filt parts for CDS +/- unkasm,
      ie. call cdsasm + unkasm > cds, unkasm only as unk
  .. teasm added
  .. TEann grew much, from teasm seq align?
  .. RPTann obs dropped much .. is this annot mistake w/ TEann?
      -- but RPTann xCopy high now .. calc or annot bug?
  .. NOann dups are moved out (to unkasm probably)
  ** report list of disjuct-subset unions: 
    uniq+dup asm = allasm, cds+te+unk = allasm - noann (70/115mb), CDSann+TEann+RPTann = allasm - NOann
    total = allasm + gaps + what?
  >> filter parts probably need update for new classes, ie TE/RPT ann vs te/unk asm
  
6b.Source=Daphmag19sk, KUlow=47, KUhigh=55, FlowcytSize=234-391 Mb Formula_LN/C=210 Mb   (dmag19skasm)
_____   ______  Low     Low     High    High    Total    
Item_   Obs.Mb  Est.Mb  xCopy   Est.Mb  xCopy   tCopy   Description
---------------------------------------------------------------------------
allasm  115     170     1.48    199     1.73    2.75    measured assembly
_total  123     178     .       207     .       .       total assembly 
uniqasm 88      91      1.03    106     1.20    1.20    asm with unique gDNA
dupasm  27      80      2.99    93      3.50    7.89    asm with multimap gDNA
cdsasm  59      112     1.90    131     2.22    3.94    asm with CDS-mapped gDNA
CDSann  43      52      1.21    61      1.42    1.67    asm with CDS annotations
CDSbus  2.0     2.0     1.00    2.4     1.18    1.20    asm with unique BUSCO orlog CDS
TEann   1.3     2.1     1.61    2.5     1.89    3.34    asm with Transposons
RPTann  5.8     6.0     1.04    7.0     1.21    1.54    asm with simple Repeats
NOann   67      113     1.68    132     1.97    3.51    asm without annotations
---------------------------------------------------------------------------
  
=item try 8

  -- unkasm reduced by 1/2 due to -cds,te filt
    .. unkasm mostly overlaps cdsasm (& CDSann)
  ** NOTE: cdsasm vs CDSann and teasm vs TEann: Obs about same, but Est higher for ann of both
     .. explain that .. which is more likely true? unkasm may account for some/all? if it can be defined
   
     
dmag19skasm_cov7bap_sum.txt
Source=Daphmag19sk, KUlow=47, KUhigh=55, FlowcytSize=234-391 Mb Formula_LN/C=210 Mb   (dmag19skasm)
_____   ______  Low     Low     High    High    Total    
Item_   Obs.Mb  Est.Mb  xCopy   Est.Mb  xCopy   tCopy   Description
---------------------------------------------------------------------------
allasm  115     170     1.48    199     1.73    2.75    measured assembly
_total  123     178     .       207     .       .       total assembly 
uniqasm 88      91      1.03    106     1.20    1.20    asm with unique gDNA
dupasm  27      80      2.99    93      3.50    7.89    asm with multimap gDNA
cdsasm  49      66      1.34    78      1.57    2.24    asm with CDS-mapped gDNA
teasm   9.8     25.3    2.58    29.6    3.02    6.10    asm with TE-mapped gDNA
unkasm  4.6     38.4    8.30    45.0    9.72    23.97   asm with UnclassRepeats-mapped gDNA
CDSann  52      103     1.99    120     2.32    4.28    asm with CDS annotations
CDSbus  1.9     1.9     1.01    2.2     1.18    1.21    asm with unique BUSCO orlog CDS
TEann   10.0    57.0    5.68    66.7    6.65    15.69   asm with Transposons
RPTann  0.3     2.9     8.29    3.4     9.70    17.41   asm with simple Repeats
NOann   60      62      1.02    72      1.20    1.39    asm without annotations
---------------------------------------------------------------------------
  xCopy = excess read copy depth; tCopy = assembly multi-map copy depth (-xCopy)

# total span=123.1.mb, covspan=114.9.mb, gaps=0.mb with n_readid=43186373 for dmag19skasm

>> 7baq has gap = biggap filter .. reduces Obs/Est (and xCopy?) for most
==> dmag19skasm_cov7baq_sum.txt <==
Source=Daphmag19sk, KUlow=47, KUhigh=55, FlowcytSize=234-391 Mb Formula_LN/C=210 Mb   (dmag19skasm)
_____   ______  Low     Low     High    High    Total    
Item_   Obs.Mb  Est.Mb  xCopy   Est.Mb  xCopy   tCopy   Description
---------------------------------------------------------------------------
allasm  113     166     1.47    194     1.72    2.70    measured assembly
_total  123     176     .       204     .       .       total assembly 
uniqasm 87      90      1.03    106     1.21    1.21    asm with unique gDNA
dupasm  26      76      2.91    89      3.41    7.69    asm with multimap gDNA
cdsasm  49      66      1.35    77      1.57    2.25    asm with CDS-mapped gDNA
teasm   9.7     25.1    2.60    29.4    3.04    6.15    asm with TE-mapped gDNA
unkasm  4.5     34.8    7.72    40.8    9.03    22.56   asm with UnclassRepeats-mapped gDNA
CDSann  51      99      1.93    116     2.26    4.13    asm with CDS annotations
CDSbus  1.9     1.9     1.01    2.2     1.18    1.21    asm with unique BUSCO orlog CDS
TEann   9.9     53.3    5.40    62.4    6.32    15.00   asm with Transposons
RPTann  0.3     2.9     8.43    3.3     9.86    17.69   asm with simple Repeats
NOann   60      62      1.02    72      1.20    1.39    asm without annotations
---------------------------------------------------------------------------
  xCopy = excess read copy depth; tCopy = assembly multi-map copy depth (-xCopy)

# total span=123.1.mb, covspan=114.9.mb, gaps=9.5.mb with n_readid=43186373 for dmag19skasm

>> 7bar : GAPann, note Obs/Est for messy case of bins w/ reads + gaps, no-read gap spans more
dmag19skasm_cov7bar_sum.txt
Source=Daphmag19sk, KUlow=47, KUhigh=55, FlowcytSize=234-391 Mb Formula_LN/C=210 Mb   (dmag19skasm)
_____   ______  Low     Low     High    High    Total    
Item_   Obs.Mb  Est.Mb  xCopy   Est.Mb  xCopy   tCopy   Description
---------------------------------------------------------------------------
allasm  115     170     1.48    199     1.73    2.75    measured assembly
_total  123     178     .       207     .       .       total assembly 
uniqasm 88      91      1.03    106     1.20    1.20    asm with unique gDNA
dupasm  27      80      2.99    93      3.50    7.89    asm with multimap gDNA
cdsasm  49      66      1.34    78      1.57    2.24    asm with CDS-mapped gDNA
teasm   9.8     25.3    2.58    29.6    3.02    6.10    asm with TE-mapped gDNA
unkasm  4.6     38.4    8.30    45.0    9.72    23.97   asm with UnclassRepeats-mapped gDNA
CDSann  52      103     1.99    120     2.32    4.28    asm with CDS annotations
CDSbus  1.9     1.9     1.01    2.2     1.18    1.21    asm with unique BUSCO orlog CDS
TEann   10.0    57.0    5.68    66.7    6.65    15.69   asm with Transposons
RPTann  0.3     2.9     8.29    3.4     9.70    17.41   asm with simple Repeats
GAPann  1.9     4.2     2.21    4.9     2.58    5.11    asm with gaps
        ^^^ gaps+reads vs 9.5mb below, for gaps w/o reads; want this? slightly useful
NOann   60      62      1.02    72      1.20    1.39    asm without annotations
---------------------------------------------------------------------------
  xCopy = excess read copy depth; tCopy = assembly multi-map copy depth (-xCopy)

# total span=123.1.mb, covspan=114.9.mb, gaps=9.5.mb with n_readid=43186373 for dmag19skasm

=cut
