#!/usr/bin/perl
# gnodes_pipe.pl for evigene/scripts/genoasm/ 

=item usage 

  gnodes_pipe.pl  -chr chr.fasta -cds cds.fasta -te te.fasta 
    opts: -asmid daphplx_gasm16ml -metadata  daphnia3genomes.metadata
    opts: -idclasses cdste.idclass : table of class per cds,te ID: CDS,TE,UNK ?RPT (simple repeat) 
            and modifiers: uniq,busco,duplicate, ..
  generates script run_gnodes.sh
  
=item run steps

  from dplx20gnodes.qsub.hist

  env genolist="genome/dropse20chrs.fa.gz genome/dplx20maca1refdeg.fa.gz genome/dcar20chr.fa.gz" datad=`pwd` qsub -q batch run_repmod2b.sh
  
  env asm=daphplx17evgt1m_cds.fa reads=readsf/SRR3090572_1.fastq.gz datad=`pwd` qsub -q debug run_gnodes_dnamap.sh
  env asm=dplx20maca1r_tefams.fa reads=readsf/SRR3090572_1.fastq.gz datad=`pwd` qsub -q debug run_gnodes_dnamap.sh
  env asm=daphplx_gasm16ml.fa reads=readsf/SRR3090572_1.fastq.gz datad=`pwd` qsub -q debug run_gnodes_dnamap.sh
  
  pt=daphplx17evgt1m_cds; env opts="-savereadids -crclassf=dplx20cdste.idclass -minident=0.45"  bam=${pt}_SRR3090572_1_bwa.bam  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh
  pt=dplx20maca1r_tefams; env opts="-savereadids -crclassf=dplx20cdste.idclass -minident=0.45"  bam=${pt}_SRR3090572_1_bwa.bam  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh
  env merge_readids="daphplx17evgt1m_cds_SRR3090572_1_bwa.readids dplx20maca1r_tefams_SRR3090572_1_bwa.readids"  cdstab=dplx20cdste_SRR3090572_1_bwa.readids  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh
  
  env opts="-crclassf=dplx20cdste.idclass" cdstab=dplx20cdste_SRR3090572_1_bwa.readids  bam=daphplx_gasm16ml_SRR3090572_1_bwa.bam  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh
  
  env idclass=dplx20cdste.idclass cds=daphplx17evgt1m_cds.fa te=dplx20maca1r_tefams.fa chr=daphplx_gasm16ml.fa datad=`pwd`   qsub -q debug run_gnodes_ann.sh
  
  env asmid=daphplx_gasm16ml covtab=daphplx_gasm16ml_SRR3090572_1_bwa.cdschr7b.covtab anntab=daphplx_gasm16ml.anntab asmdata=daphplx20chrs.metad  datad=`pwd` qsub -q debug run_gnodes_sum.sh
  
  #--- 2nd read set, want to merge covtab,chrtab for new sum table

=cut

use FindBin;
use lib ("$FindBin::Bin/../"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
our $EVIGENES="$FindBin::Bin/..";  
use strict;
use Getopt::Long;  

$EVIGENES =~ s,/\w+/\.\.$,,;
# export EVIGENES=/Very/long/system/path/to/bio/evigene/scripts/genoasm/..

use constant UPD21AUG => 1; # chrplot add, other fixes
use constant UPD21JUL => 1; # gnodes_sumgenecov.pl, finally tested out
use constant UPD21JUN => 1;
# UPD21JUN: replace ucgcov w/ sam2genescov, revised annote, add RDMAPPER = minimap (options?)
# ** need still for genescov: mrna input > cdhit-est -c 0.98? reduce duplicates, then readmap, separate from cds-rdmap
# .. for improved genes copynum table
use constant UPD21OCT => 1; # UPD21OCT to gnodes_sam2covtab8j,  gnodes_sumgenecov8j, new pipe calls

# UPD21OCT public release, renaming scripts "8j" to base name: gnodes_pipe,_sam2covtab,_sumgenecov
my $VERS="8j"; #UPD21OCT; should upd to 8b,c,d,e? 
my $debug=$ENV{debug}||1;
my $ncpu=$ENV{NCPU}||16; my $maxmem=$ENV{maxmem}||"64gb"; # cluster opts
   my $walltime= $ENV{walltime} || "1:00:00";
my($doREPMASK,$doBUSCO)= (1,1); # default write scripts, turn off
my $buscodb= $ENV{buscodb} || ""; # "eukaryota"; # or not?
my $RDMAPPER=$ENV{rdmapper} ||'minimap2'; # or 'bwa-mem' or 'bowtie2' or what? global opt or setup opt ?
my $DOGENETAB=(UPD21JUN)?1:0;
   
my $gnodes1_dnamap=""; # local proc .. "$EVIGENES/genoasm/gnodes_dnamap.pl";
my $gnodes1_sam2cov="\$EVIGENES/genoasm/gnodes_sam2covtab.pl"; # 8j; sam2covtab8e, gnodes_sam2covtab8a for UPD21JUN
  # UPD21OCT: gnodes_sam2covtab8j.pl major algo revision since 21JUN, drops readids huge mem pig for pair read of cds.bam x chr.bam
my $gnodes2_ann="\$EVIGENES/genoasm/gnodes_annotate.pl";
my $gnodes3_covsum="\$EVIGENES/genoasm/gnodes_covsum.pl";

# FIXME: rename/link gnodes_sam2covtab8j > gnodes_sam2covtab, gnodes_sumgenecov8j > gnodes_sumgenecov
my $gnodes3_sam2genecov="\$EVIGENES/genoasm/gnodes_sam2genecov.pl"; #UPD21JUN replace ucgcov,genescov .. step 2a before annot now
my $gnodes4_sumgenecov="\$EVIGENES/genoasm/gnodes_sumgenecov.pl"; #8j; UPD21JUL : replaces gnodes_genescov 
  # UPD21OCT: gnodes_sumgenecov8j.pl
  # gnodes_sumgenecov collates, summarizes sam2covtab xx.genetab, sam2genecov xx.genexcopy, annotate ?genecopyn tables

my $gnodes3_ucgcov="\$EVIGENES/genoasm/gnodes_sam2ucgcov.pl"; #obsolete w/ UPD21JUN
my $gnodes4_genescov="\$EVIGENES/genoasm/gnodes_genescov.pl"; #obsolete UPD21JUN, upd21apr20
# add maybe: "\$EVIGENES/cdna_bestorf.pl" for sh_mrna2cds()

my $intitle= $ENV{title}||$ENV{name}||""; 
my $asmid; #no default $intitle; # intitle for output, asmid for data
my ($chrasm,$teseq,$cdsseq,$mrnaseq,$crclassf,$sampledata,$anntable,$datad,$gnodes_setup,$mergecov,$species)=("") x 19;
my($gidprefix,$ridprefix)=("",""); #UPD21OCT; get from data
my @reads;

my $optok= GetOptions( 
  'chrasm|assembly|genome=s',\$chrasm,
  'teseq=s',\$teseq,
  'cdsseq=s',\$cdsseq,  
  'mrnaseq=s',\$mrnaseq, # UPD21JUN, maybe replace cdsseq w/ this, mainly for precision in rd depth: mrna x rdmap > sam2genes copynum ; needs test
  'reads=s',\@reads, # many, esp SRRnnn_[12].fq pairs
  'RDMAPPER=s', \$RDMAPPER, # UPD21JUL: limit to bwa and minimap2 for now
  'outputdir|datad=s',\$datad,  
  'merge:s',\$mergecov,  
  'title|name=s',\$intitle,
  'asmid=s', \$asmid, # superceeds intitle for finding/using asm info
  'genomedata|sumdata|sampledata|metadata=s', \$sampledata, # optional input  << change to gnodes_pipe.cfg/config ?
  'anntable|annotations=s', \$anntable, # not an option?
  'idclass=s',\$crclassf, # table of chr/gene id => class ? want for BUSCO, TEfam, : No, rely on anntable having these
  'species=s', \$species, # repeatmasker opt only?
  'ncpu=i', \$ncpu, 'maxmem=s', \$maxmem, 'walltime=s', \$walltime,
  'REPMASK!', \$doREPMASK, # norepmask turns off
  'BUSCOdb=s', \$buscodb, 
  'debug!', \$debug, 
  );

## use metad/sampledata for opts: buscodb=  species=  repmaskdb==species=? and other params: title? ..
# check @ARGV for .fastq/.fq reads ?
if(@ARGV and my @rds= grep( /\.(fastq|fq)/, @ARGV)) { push @reads, @rds; }

unless($optok and (($chrasm and @reads) or $mergecov) ) {
die "usage: $0 -chrasm mychr.fa -cds mygenecds.fa -reads SRRnnn.fastq .. 
options: -title outputname -asmid chrasm_id -metadata chrasm.meta 
         -idclass genes.idclass -anntable chrasm.anntab  : in/out info tables
         -teseq mytransposons.fa  : use these instead of RepeatMasker( mychr.fa)
         -ncpu $ncpu -maxmem $maxmem -walltime $walltime : compute cluster opts
         -merge x.covtab y.covtab z.covtab : merge several part results
This gnodes_pipe makes run_scripts.sh to compute several gnodes_ steps, given input files. 
Initialize with local gnodes_setup.sh or ENV{gnodes_setup} file.
It does no computations (safe to run on 1 cpu)
Version: $VERS  
\n";
}

sub basename { my($f,$suf)=@_; $f=`basename $f`; chomp($f); $f=~s/$suf// if($suf); return($f); }

unless($asmid){ $asmid= basename($chrasm,'\.fa.*'); } #? or not
# $asmid= $intitle if($intitle and not $asmid); # no, use chrasm name

my($nmeta,$metavals,$asmidlist)= readMetad($sampledata);
# metavals->{asmid}{keys}= vals
if($nmeta and not $species){ $species = $metavals->{$asmid}{species}||""; }
if($nmeta and not $intitle){ $intitle = $metavals->{$asmid}{asmname}||""; }
$intitle=$asmid unless($intitle);

# need datad path:
my $mydir=`pwd`; chomp($mydir); # user option? make new work dir?
if($datad) {
  unless($datad =~ m,^/, or $mydir =~ m,$datad,) { $datad="$mydir/$datad"; }
} else {
  $datad=$mydir;
}  

# add new main method: -merge_covtabs -out daphpulex_pa42v2_merge3r daphpulex_pa42v2_SRR*bwa.cc7c.covtab
# also if @reads > 1, do both: gen script from gnodes_pipe_MAIN for each @reads, then gnodes_mergecov_MAIN

if(@reads > 1) {
  #UPD21apr20: change sh_dnamap(..,@reads) to map all read sets, samtools cat -o all.bam parts.bam,
  #  then use all.bam for rest of calcs
if(1) {
  my ($reads0)= shift @reads;
  my($runsh)= gnodes_pipe_MAIN( $reads0, 0, @reads );
  warn "$runsh\n";

} else {  
  my @covtabs;
  for my $i (0..$#reads) {
    my($partsh,$chrcovtab,$sumout)= gnodes_pipe_MAIN( $reads[$i], 1+$i );
    push @covtabs, $chrcovtab;
    warn "$partsh\n";
  }
  my($runsh)= gnodes_mergecov_MAIN( @covtabs);
  warn "$runsh\n";
}
  
} elsif($mergecov) {
  my @covtabs= grep(/covtab$/, @ARGV);
  my($runsh)= gnodes_mergecov_MAIN( @covtabs );
  warn "$runsh\n";
  
} else {
  my($runsh)= gnodes_pipe_MAIN( $reads[0] );
  warn "$runsh\n";
}


#=======================================================


sub gnodes_pipe_MAIN {  
  my($reads,$it,@morereads)=@_;
  # my $reads= $reads[0]; # fixme
  
  return "# NO reads to map" unless($reads);
  if($it){ $it = "_r".$it; } else { $it=""; }
  my $runsh="run_gnodes_$intitle$it.sh"; # need read $i num
  
  my $OUTSH; open($OUTSH,'>',$runsh);
  
  ($gnodes_setup)= sh_setup($OUTSH,$runsh);

  #? insert here call to sh_repeatmask() sh_buscoscan ?
  #  want results before sh_annot(), sh_dnamap(te) are run
  my($rmout,$buscout)=("",""); # options
  my $teann=$teseq;
  if($doREPMASK and not $teseq) {  
    my $rmaskdb= $metavals->{$asmid}{rmaskdb} || $species;
    ($rmout)= sh_repeatmask($OUTSH,$chrasm,$rmaskdb);
    $teann=$rmout; 
    # NOT $teseq= $rmout; # okay? << BUG at sh_dnamap($OUTSH,$teseq,$reads)
    ##  asmidx="bwax/ficari1asm_bwx"
    ##  bwa-mem2 index -p $asmidx ./rmout/ficari1asm.fa.out <<<
  }
  
  
  if($mrnaseq and not $cdsseq) { #UPD21AUG17: make cds from mrna if needed
    ($cdsseq)= sh_mrna2cds($OUTSH,$mrnaseq);
  }
  
  if($doBUSCO and $cdsseq) {  
    $buscodb ||= $metavals->{$asmid}{buscodb} || $ENV{buscodb} || "";
    ($buscout)= sh_buscoscan($OUTSH,$cdsseq,$buscodb);
  }

  
  # this can be run independently of dnamap/sam2cov .. call before those, generates $crclassf
  # UPD run this AFTER sh_sam2genecov(), ie can use that genexcopy table
  my($anout,$idclassa)=("","");
  unless(UPD21JUN) {
    ($anout,$idclassa)= sh_annot($OUTSH,$chrasm,$cdsseq,$teann,$crclassf,$buscout);
    $anntable= $anout; # is anntable option or not?
    $crclassf= $idclassa unless($crclassf);
  }
    # ? pre-make crclassf if doesnt exist and have cdsseq,teseq .. see also gnodes_annotate.pl make_idclass()
  
  #UPD21apr20: sh_dnamap_many : samtools cat -o all.bam @parts.bam, or for 1 reads does orig way
  #UPD21apr20: FIXME this way causes readids _1/_2 to be same, spurious duplicates, one bug is that
  #   cdsasm count now ~2x of before, and 2x of CDSann when should be nearly same. Other bugs not yet found.
  #   solution? add .1/.2 suffix to read ids? maybe ok but ID number is used as integer. Change IDprefix?
  #  add sam flag for pair/mate?  that means more recoding w/ sam parsers
  # .. cant easily fix here w/o edit read input files.. fix in sam2covtab ?

  # UPD21JUN sh_dnamap( @morereads==() ) is ok now; drop sh_dnamap_many

  my($gnbam,$tebam,$crbam);
  ($gnbam)= ($cdsseq)? sh_dnamap($OUTSH,$cdsseq,$reads,@morereads) : 0;
  ($tebam)= ($teseq )? sh_dnamap($OUTSH,$teseq,$reads,@morereads) : 0;
  ($crbam)= ($chrasm)? sh_dnamap($OUTSH,$chrasm,$reads,@morereads) : 0;
  
  #UPD replace sh_ucgcov()
  my($genescovtab,$ucgcovtab)=("","");
  if( UPD21JUN ) {
    # problem: sh_annot returns crclassf file name ..
    #   sh_sam2genecov covtab == $genecdsname.genexcopy; #?? what suffix
    
    #UPD21Aug01: missing crclassf means no UCG/busco calc here.  try to find it before sh_annot()
    unless($crclassf){ my $idc=$cdsseq; $idc =~ s/\..*/.idclass/;  $crclassf= $idc if(-s $idc); }
    
    if($mrnaseq and $mrnaseq ne $cdsseq) { # mrnaseq can/should be cdhit -c0.9[78] reduced transcript set
      my($mrnabam)= sh_dnamap($OUTSH,$mrnaseq,$reads,@morereads);
      ($genescovtab)= sh_sam2genecov($OUTSH, $mrnabam, "",  $crclassf, "norecalc") ;
    } elsif($gnbam) {
      ($genescovtab)= sh_sam2genecov($OUTSH, $gnbam, "",  $crclassf, "norecalc") ;
    }
  
    ($anout,$idclassa)= sh_annot($OUTSH,$chrasm,$cdsseq,$teann,$crclassf,$buscout, $genescovtab);
    $anntable= $anout; # is anntable option or not?
    $crclassf= $idclassa unless($crclassf);

  } else { # old
    ($ucgcovtab)= ($gnbam) ? sh_ucgcov($OUTSH, $gnbam, "",  $crclassf, "norecalc") : ("");
    #^ maybe insert before temap,chrmap
  }

=item UPD21AUG: BIG-PIG  gnodes pipe algo changes

  samsplit2cov for pig genome
  see genoasm/run_gdebug_samsplit2cov8g.sh
  .. note may want to estimate data size to decide if need this algo change?
  .. if reads.fastq[.gz] exist, check size to decide if big-pig split algo needed?
     or use for all/small genomes?
  -- two main changes, split-bams & add-readids-to-chr-bam, are separate,
      split-bams only for big-pig data? but add cds.ids to chr.bam for all?
      
  a. split cds,te,chr bams to ncpu parts, after sh_dnamap (includes cat -o all.bam parts*.bam)
     gnodes_samaddrdid.pl  -split -nparts=NN -bam cds|te|chr.bam
     .. need to track/use partnames here in pipeline? partnames = name.pt00.typesuf
     
  b. sh_sam2covtab( cds|te) on parts.bam, each w/ parts.readids,covtab outputs
     .. merge parts.covtab for sh_covsum(cdstab)
     
  c0. sh_merge_readids() is obsolete? now merges only cds+te subsets: $cdsreadids,$tereadids
      gnodes_sam2cov -merge readids -savereadids=$mergedids @readid_parts
      .. this creates class columns: readid  nrd  CDS TE UNK geneids 
      maybe not obsolete, but want it to work w/ parts.readids, input to samaddrdid
  
  c. gnodes_samaddrdid.pl -readids cds.parts.readids -bam chr.parts.bam -outbam chr.parts.rdid.bam
     .. replace chr.parts.bam w/ parts.rdid.bam
     .. ensure option for 2+ readid files per part: -readids  CDS.ptN.rids -readids TE.ptN.rids -readids UNK.ptN.rdids
     
  d. sh_sam2covtab(chr) using chr.parts.rdid.bam w/ cds.readids inserted
     .. merge chr.parts.covtab,chrtab,genetab for sh_covsum(chrcovtab, genescovtab), follow on stats
     
     
=cut

  my($cdstab,$cdsreadids)= ($gnbam) ? sh_sam2covtab($OUTSH, 'cds', $gnbam, $crclassf) : (0,0);
  my($tetab,$tereadids)  = ($tebam) ? sh_sam2covtab($OUTSH, 'te', $tebam, $crclassf) : (0,0);

  # insert sh_covsum($cdstab) .. gets KU,C value independent of chrasm  
  my($sumcds)=(0);
  if($cdstab) {
    my $cdsid= basename($cdsseq,'\.(fa|cds).*'); 
    # fixme: crclassf is output of sh_annot(): ficari18tsa1cds.fa > ficari18tsa1cds.idclass
    #here? ($ucgcov)= sh_ucgcov($OUTSH, $cdstab, $cdsid, $cdsid, "", $crclassf, "norecalc");  
    ($sumcds)= sh_covsum($OUTSH, $cdstab, $cdsid, $cdsid, "", $crclassf, $genescovtab, "norecalc");  # no cds.anntab but idclass
  }
  
  my($chrcovtab);
  if(UPD21OCT) { # no readids, expect gnbam == cdsbam merge w/ tebam if avail
    # if($gnbam and $tebam) { (my $gntebam=$gnbam) =~ s/\./_te./; sh_merge_bams($gntebam, $gnbam,$tebam);  $gnbam=$gntebam; }

    #global opts here -gidprefix=$gidprefix -ridprefix=$ridprefix    
    ## gidprefix now found from crclassf ids
    if($gnbam and $crbam) {
    ($chrcovtab)= sh_sam_chr_cds2covtab($OUTSH, $crbam, $gnbam, $crclassf); 
    } else {
    ($chrcovtab)  = sh_sam2covtab($OUTSH, 'chr', $crbam, $crclassf); #?? readids still viable opt
    }
  } else {
    my $mctids= $cdstab; $mctids=~s/.covtab//; $mctids=~s/cds//; $mctids.="_cdste.readids";
    my($cdsteids) = (not $tereadids) ? $cdsreadids 
        : (not $cdsreadids) ? $tereadids
        : sh_merge_readids($OUTSH,$mctids,$cdsreadids,$tereadids);
    
    ($chrcovtab)  = sh_sam2covtab($OUTSH, 'chr', $crbam, $crclassf, $cdsteids);
  }
  
  my($crsumout)= sh_covsum($OUTSH, $chrcovtab, $asmid, $intitle, $anntable, $crclassf, $genescovtab); 
    #   genescovtab => covsum -genexcopy=$genescovtab, # UPD21JUN sam2genecov covtab == $genecdsname.genexcopy; 

  if( UPD21AUG ) { # chrplots from chrcovtab
    # gnodes_covsum.pl -plotchr -asmid $asmid  -title $outprefix -metad $sampledata
    #    -genexcopy $genescovtab -anntab $anntable -crclass $crclassf   $chrcovtab
    my($chrplots)= sh_chrcovplot($OUTSH, $chrcovtab, $asmid, $intitle, $anntable, $crclassf, $genescovtab);  
  }

  my($gsumout)=("");
  if(UPD21OCT){
    my $chrgenetab= $chrcovtab; $chrgenetab =~ s/.covtab/.genetab/;
    ($gsumout)= sh_sumgenescov8j($OUTSH, $genescovtab, $chrgenetab, $crsumout, $intitle, $cdstab); # updoct31: $cdstab
  
  } elsif( UPD21JUL ) { # revised genescov
    # sh_sumgenescov($OUTSH, $genexcopy, $chrgenetab, $lasmid, $lintitle,  $lidclassf, $lanntable) 
    #  genexcopy = gene copynum table,  primary data from gnodes_sam2genecov
    #  chrgenetab = chr x gene-cds cov table,  primary data, from updated sam2covtab:putGenetab() : GeneID ChrID Pos aCovT aCovM noCov
    my $chrgenetab= $chrcovtab; $chrgenetab =~ s/.covtab/.genetab/;
    ($gsumout)= sh_sumgenescov($OUTSH, $genescovtab, $chrgenetab, $asmid, $intitle, $crclassf); # old? opt: , $anntable
  
  } elsif( ! UPD21JUN ) { #?? maybe keep w/ update, revise ?
    ($gsumout)= sh_genescov($OUTSH, $cdstab, $chrcovtab, $asmid, $intitle, $anntable, $crclassf); 
  }
  
  close($OUTSH);  
  system("chmod +x $runsh");
  return($runsh,$chrcovtab,$crsumout,$gsumout); # ,$sumcds
}

sub gnodes_mergecov_MAIN {
  my (@covtabs)= @_;
  # my @covtabs= grep /covtab$/, @ARGV; #?? or -merge a.covtab,b.covtab,...
  
  # gnodes_pipe.pl -merge -title daphpulex_pa42r2 -chr daphpulex_pa42v2.fa  -sumdata daphplx20chrs.metad -ncpu 24  -maxmem 128gb
 
  my $nc=@covtabs;
  return "# NO covtabs to merge" unless($nc>1);
  
  my $mergetab= ($mergecov||$intitle)."_merge$nc.$VERS.covtab";  
  unless($anntable) {
    my $aname= ($chrasm) ? basename($chrasm,'\.fa.*') : $asmid;
    $anntable="$aname.anntab";
  }
  
  my $runsh="run_gnodes_merge_$intitle.sh";
  my $OUTSH; open($OUTSH,'>',$runsh);
  
  ($gnodes_setup)= sh_setup($OUTSH,$runsh);

  my($mcovtab,$mchrtab)= sh_merge_covtabs( $OUTSH, $mergetab, @covtabs);

  my($genescovtab,$ucgcovtab)=("","");
  if( UPD21JUN ) {
    $genescovtab= $covtabs[0];
    $genescovtab =~ s/.covtab//; $genescovtab.=".genexcopy"; # this can be bad name ***
  }
  
  my($crsumout)= sh_covsum($OUTSH, $mcovtab, $asmid, $intitle, $anntable, $crclassf, $genescovtab); 

  my($gsumout)=("");
  if(UPD21OCT){
    my $chrgenetab= $mcovtab; $chrgenetab =~ s/.covtab/.genetab/;
    my $mcdscov=""; #UPD21OCT31: FIXME now want $cdstab for sumgenescov
    ($gsumout)= sh_sumgenescov8j($OUTSH, $genescovtab, $chrgenetab, $crsumout, $intitle, $mcdscov); 
  } elsif( UPD21JUL ) { # revised genescov
    # sh_sumgenescov($OUTSH, $genexcopy, $chrgenetab, $lasmid, $lintitle,  $lidclassf, $lanntable) 
    #  genexcopy = gene copynum table,  primary data from gnodes_sam2genecov
    #  chrgenetab = chr x gene-cds cov table,  primary data, from updated sam2covtab:putGenetab() : GeneID ChrID Pos aCovT aCovM noCov
    my $chrgenetab= $mcovtab; $chrgenetab =~ s/.covtab/.genetab/;
    ($gsumout)= sh_sumgenescov($OUTSH, $genescovtab, $chrgenetab, $asmid, $intitle, $crclassf); # old? opt: , $anntable

  } elsif( ! UPD21JUN ) {
    my $mcdscov=""; # FIXME, merge cds.covtabs as well as chr.covtabs
    ($gsumout)= sh_genescov($OUTSH, $mcdscov, $mcovtab, $asmid, $intitle, $anntable, $crclassf); 
  }
  
  close($OUTSH);  
  system("chmod +x $runsh");
  return($runsh);
}


# run_gnodes_dnamap.sh script ============================

sub sh_setup {
  my($OUTSH,$sname)=@_;

  # need datad path: above now, add mkdir datad
  # my $datad=`pwd`; chomp($datad); # user option? make new work dir?
  # FIXME maybe here? set shell VARs for primary data here, those use to test results of each step
  # but now sh_ subs for each method are setting output names

  # FIXME2: this part should be customizable from local data : {local}gnodes_setup.sh ??
  my $localset= $ENV{gnodes_setup} || "gnodes_setup.sh";
  my $topset="";
  if( -f $localset ) {
    open(LS,$localset); $topset= join "", <LS>; close(LS);
  } else {
    $topset=<<"EOS";
## --- gnodes_setup.sh ---    
#PBS -N gnodes_pipe
#PBS -l vmem=${maxmem},nodes=1:ppn=$ncpu,walltime=${walltime}
#PBS -V

## ensure have these.. trap missing?  
module load blast samtools 
module load repeatmasker busco
# module load minimap2
# module load bwa-mem # older
bwabin=\$HOME/bio/bwa/bin/; export PATH=\$bwabin:\$PATH
## --- end gnodes_setup.sh ---    
EOS
  }
  
  my $script=<<"EOS"; 
#! /bin/bash
## $sname
## create 'gnodes_setup.sh' to add local env, paths, cluster metacomments: PBS/SBATCH/XXX ...  
$topset

export EVIGENES=$EVIGENES
export PATH=\$PATH:\$EVIGENES
export NCPU=$ncpu;
datad=$datad

if [ ! -d \$datad ]; then mkdir \$datad; fi
cd \$datad
mkdir parts

EOS

  print $OUTSH $script;
  return($script); # save $topset to find RDMAPPER ??
}

sub bad_exec { # or check_exec()
  my @exfail=(); # my %exfail;
  for my $e (@_) {
    # my $err= system("$e --help >/dev/null 2>&1"); #  --help is no good for some, use which
    my $err= system("which $e > /dev/null 2>&1");  
    push @exfail, $e if($err);
  }
  return @exfail;
}

sub readMetad {  # from gnodes_covsum.pl
  my($sampledata)= @_;
  my($nsam,$aid)=(0,0);
  my %samvals=(); my @aid=(); 
  unless($sampledata) { # want empty %samvals
    return($nsam, \%samvals, \@aid);
  }
  warn "#read sampledata=$sampledata\n" if($debug);
  open(F,$sampledata); 
  while(<F>){
    next if(/^\W/);
    my($key,$val)= (m/^(\w+)\s*=\s*(.+)$/) ? ($1,$2):(0,0);
    next unless($key);
    
    if($key eq 'pt' or $key eq 'asmid') { $aid=$val; $nsam++;
      map{ $samvals{$aid}{$_}= 0; } qw(flowcyto cytomb asmtotal atotalmb asmname );
      push @aid, $aid;
      
    } elsif($key eq 'flowcyto') { # /^flowcyto=(.+)$/)  
      my $flowcyto= $val; # flowcyto=234-391 Mb
      $samvals{$aid}{flowcyto}= $flowcyto;
      $samvals{$aid}{cytomb}  = ($flowcyto=~m/(\d+)/)?$1:0; # can be range 160-180
  
    } elsif($key eq 'glncformula') { # /^glncformula=(.+)$/ 
      $samvals{$aid}{glncformula}= $val; # glncformula=nnn Mb
      
    } elsif($key eq 'asmtotal') { # /^asmtotal=(.+)$/)   
      my $asmtotal= $val; 
      $samvals{$aid}{asmtotal} = $asmtotal;
      $samvals{$aid}{atotalmb} = ($asmtotal=~m/(\d+)/)?$1:0; # can be range 160-180
  
    # } elsif($key =~ /asmname|kulo|kuhi/) { # ($key eq 'asmname' $key eq 'kulo' or $key eq 'kuhi' )
    #   $val=~s/\s.*//; $samvals{$aid}{$key}= $val;  # /^kulo=(\S+)/)  # test hi<lo
    # } elsif($key =~ /busco|rmask|species/) {
    #   $val=~s/\s.*//; $samvals{$aid}{$key} = $val;
    } else { # save all valid key = val
      $val=~s/\s.*//; $samvals{$aid}{$key} = $val;
    }
    
  } close(F);
  push @aid, 0 unless(@aid);
  return($nsam, \%samvals, \@aid); 
}


sub sh_dnamap {
  my($OUTSH,$asm,$reads,@breads)=@_;
 
  return("") unless($asm and $reads);
  
  # RDMAPPER should be prog name? minimap2 bwa-mem  * Add bowtie2 opt
  my($DO_MIM,$DO_BWA,$BWAapp,$MIMapp)=(0,0,'bwa-mem2','minimap2');
  if(UPD21AUG) { #UPD21AUG17, check/bad_exec(mappers)
    #* use gnodes_setup.sh to find RDMAPPER ?
    my $rdxok=0;
    $RDMAPPER ||= $BWAapp;
    my($badx)= bad_exec($RDMAPPER);
    if( $badx ) { 
      my @gset= split/\n/,$gnodes_setup;  map{ s/#.*//; } @gset; 
      my($rdpath)= grep( m/$RDMAPPER/, @gset);
      if($rdpath) { $badx= 0; $rdxok=3; } # set ok ** THIS TEST FAILS
    }
    
    if( $badx ) { 
      warn "#ERR missing exec RDMAPPER=$RDMAPPER\n"; 
      if($RDMAPPER =~ m/minimap/) {
        if(! bad_exec('bwa-mem2')) { $RDMAPPER=$BWAapp='bwa-mem2'; $rdxok=2; }
        elsif(! bad_exec('bwa')) { $RDMAPPER=$BWAapp='bwa'; $rdxok=2; } 
      } elsif($RDMAPPER =~ m/bwa/ and ! bad_exec($MIMapp)){ 
        $RDMAPPER=$MIMapp; $rdxok=2; 
      }
    } else {
      $rdxok=1; 
    }
    
    $rdxok= -1 unless($rdxok);  # allow, may find later
    if($rdxok) {
      if($RDMAPPER =~ m/minimap/){ $DO_MIM=1; $MIMapp=$RDMAPPER; }
      elsif($RDMAPPER =~ m/bwa/) { $DO_BWA=1; $BWAapp=$RDMAPPER; } # bwa-mem2 or original bwa mem ?
      warn "#OK found RDMAPPER=$RDMAPPER\n" if($rdxok>1);
      warn "#ERR not-found RDMAPPER=$RDMAPPER\n" if($rdxok<0);
    }
     
  } else { # old
    if($RDMAPPER =~ m/minimap/){ $DO_MIM=1; }
    elsif($RDMAPPER =~ m/bwa/){ $DO_BWA=1; } # bwa-mem2 or original bwa mem ?
  }
  
  my $btag= ($DO_MIM) ? "_mim" : ($DO_BWA) ? "_bwa" : "_rmap";
  
  #UPD21AUG13: older samtools doesnt know --threads, requires '-' input sam param .. use that, ok for recent samt
  #old samt, but need --threads/-\@ $NCPU: my $xsam2bam="samtools view -Sb"; # add "-o $outbam -"
  my $xsam2bam="samtools view --threads \$NCPU -Sb"; # add -o $outbam
  
  my $aname= basename($asm,'\.(fa|fasta)'); # bad for teseq= chrasm.fa.out !*** '\.fa\w+ .gz
  my $nab=   basename($reads,'\.(fastq|fq)'); 
    map{ s/\.gz// } ($aname,$nab);
  my $name= $aname."_".$nab;
  my $outbam = "$name$btag.bam";
  my $bwalog = "$name$btag.log"; 

  # minimap all short reads: -x sr -N 999999 -p 0.5 << change to -p 0.80 == default ? p=2nd maps diff from 1st
  # UPD21AUG17 boost -p=0.80 or 0.90 or 0.95, sam2cov uses 0.98 for dupmin, needs test, reduces bam/speeds calc
  # orig: mimopt="-x sr -N 999999 -p 0.80 --secondary=yes"
  # UPD22APR14 : dont need *all* multimaps, adds cost in time/disk, should this be option? bwa needs change also
  my $script_mim=<<"EOM";
  mimopt="-x sr -N 9 -p 0.80 --secondary=yes"
  ( $MIMapp \$mimopt -a -t \$NCPU $asm $reads | $xsam2bam -o $outbam - ) > $bwalog 2>&1
EOM

  #above: $BWAapp= ($RDMAPPER =~ m/bwa-mem2/)? 'bwa-mem2' : 'bwa';
  my $script_bwa=<<"EOB";
  asmidx="bwax/${aname}_bwx"
  if [ ! -f \$asmidx.ann ]; then
    mkdir bwax; ( $BWAapp index -p \$asmidx $asm ) > $bwalog 2>&1
  fi
  ( $BWAapp mem -a -t \$NCPU \$asmidx $reads | $xsam2bam -o $outbam - ) > $bwalog 2>&1
EOB

  my $sh_map= ($DO_MIM) ? $script_mim : ($DO_BWA) ? $script_bwa : "echo missing read mapper";
  unless($DO_MIM or $DO_BWA) {
    warn "# ERROR: this read-mapper is not configured yet: $RDMAPPER\n# USE minimap2 or bwa|bwa-mem2\n";
  } 
  if(@breads) {
    my $nbreads= 1 + @breads;
    my $allbam= $outbam; $allbam=~s/_[12]//; $allbam=~s/$btag/_b$nbreads$btag/;

  #UPD21apr20: FIXME this way causes readids from _1/_2 to be same, spurious duplicates, one bug is that
  #   cdsasm count now ~2x of before, and 2x of CDSann when should be nearly same. Other bugs not yet found.
  #   solution? add .1/.2 suffix to read ids? maybe ok but ID number is used as integer. Change IDprefix?
  # .. cant easily fix here w/o edit read input files.. fix in sam2covtab ?
  #FIXED: s_id2 script insert:  mapper | s_id2 | samtools  -o bam
  
    my @partbams= ($outbam);
    for my $readb (@breads) {
      my $nab=   basename($readb,'\.(fastq|fq)'); map{ s/\.gz// } ($nab);
      my $name = $aname."_".$nab;
      my $bamb = "$name$btag.bam"; 
      my $logb = "$name$btag.log"; # reuse 1 log since this is serial? +> $logb
      
      my $s_id2 = ($readb =~ m/_2/) ? " perl -pe 's=\\t=/2\\t= unless(/^\\@/);' |" : "";
      push @partbams, $bamb;
      $sh_map .= 
        ($DO_MIM) ? "  ( $MIMapp \$mimopt -a -t \$NCPU $asm $readb | $s_id2 $xsam2bam -o $bamb - ) > $logb 2>&1\n"
      : ($DO_BWA) ? "  ( $BWAapp mem -a -t \$NCPU \$asmidx $readb | $s_id2 $xsam2bam -o $bamb - ) > $logb 2>&1\n"
      : "echo missing read mapper";
    }
    
    $sh_map .= "\n  samtools cat -o $allbam " . join(" ",@partbams) . "\n";
    $outbam= $allbam;
  }
  
  my $script=<<"EOS"; 
# run_gnodes1_dnamap $asm
if [ -s $outbam ]; then 
  echo reusing $outbam;
else
  echo START_dnamap  `date`
  $sh_map
  echo DONE_dnamap `date`
fi

EOS
  
  print $OUTSH $script;
  return($outbam);
}


# sub sh_dnamap_many {  #UPD21apr20 : UPD21JUN replace w/ sh_dnamap(..,@breads)
#   my($OUTSH,$asm,$reads,@breads)=@_;
# 
#   return("") unless($asm and $reads);
#   return sh_dnamap($OUTSH,$asm,$reads) unless(@breads);
#   
#   # RDMAPPER should be prog name? minimap2 bwa-mem
#   my $DO_MIM=($RDMAPPER =~ m/minimap/);
#   my $DO_BWA=($RDMAPPER =~ m/bwa/);
#   my $btag= ($DO_MIM) ? "_mim" : ($DO_BWA) ? "_bwa" : "_rmap";
#   
#   my $aname= basename($asm,'\.(fa|fasta)'); # bad for teseq= chrasm.fa.out !*** '\.fa\w+ .gz
#   my $nab=   basename($reads,'\.(fastq|fq)'); 
#     map{ s/\.gz// } ($aname,$nab);
#   my $name= $aname."_".$nab;
#   # FIXME HERE: outbam > allbam to test reusing
#   my $outbam = "$name$btag.bam";
#   my $bwalog = "$name$btag.log"; 
# 
#   my $nbreads= 1 + @breads;
#   my $allbam= $aname."_".$nab; $allbam=~s/_[12]$//;
#   $allbam .= "_b$nbreads$btag.bam";
#   my $outbam1= $outbam; $outbam=$allbam;
# 
#   # minimap all short reads: -x sr -N 999999 -p 0.5
#   my $script_mim=<<"EOM";
#   mimopt="-x sr -N 999999 -p 0.50 --secondary=yes"
#   ( minimap2 \$mimopt -a -t \$NCPU $asm $reads | samtools view --threads \$NCPU -Sb -o $outbam1 ) > $bwalog 2>&1
# EOM
# 
#   my $script_bwa=<<"EOB";
#   asmidx="bwax/${aname}_bwx"
#   if [ ! -f \$asmidx.ann ]; then
#     mkdir bwax; ( bwa-mem2 index -p \$asmidx $asm ) > $bwalog 2>&1
#   fi
#   ( bwa-mem2 mem -a -t \$NCPU \$asmidx $reads | samtools view --threads \$NCPU -Sb -o $outbam1 ) > $bwalog 2>&1
# EOB
#   
#   my $smap= ($DO_MIM) ? $script_mim : ($DO_BWA) ? $script_bwa : "echo missing read mapper";     
#   my $script=<<"EOS"; 
# # run_gnodes1_dnamap
# if [ -s $allbam ]; then 
#   echo reusing $allbam;
# else
#   echo START_dnamap  `date`
#   $smap
# EOS
#   
#   my @partbams= ($outbam1);
#   for my $readb (@breads) {
#     my $nab=   basename($readb,'\.(fastq|fq)'); map{ s/\.gz// } ($nab);
#     my $name= $aname."_".$nab;
#     my $bamb = "$name$btag.bam"; 
#     my $logb = "$name$btag.log"; 
#     push @partbams, $bamb;
#     $script .= 
#     ($DO_MIM) ? "  ( minimap2 \$mimopt -a -t \$NCPU $asm $readb | samtools view --threads \$NCPU -Sb -o $bamb ) > $logb 2>&1\n"
#   : ($DO_BWA) ? "  ( bwa-mem2 mem -a -t \$NCPU \$asmidx $readb | samtools view --threads \$NCPU -Sb -o $bamb ) > $logb 2>&1\n"
#   : "echo missing read mapper";
#   }
#   
#   $script .= "\n  samtools cat -o $allbam " . join(" ",@partbams) . "\n";;
#   $script .= "  echo DONE_dnamap `date`\n";
#   $script .= "fi\n\n";
#  
#   print $OUTSH $script;
#   return($allbam);
# }



=item sh_sam_chr_cds2covtab UPD21OCT

  $EVIGENES/genoasm/gnodes_sam2covtab8j.pl -genetab -icpu $i -ncpu $NCPU -samcpu $samcpu -debug  -crclassf=$idclass \
     -gidprefix=$gidprefix -ridprefix=$ridprefix \
     -cdsbam $cdsbam -out $covout -bam $chrbam &

  .. -samcpu opt seems not helpful, ignore.. 
  .. gidprefix can reduce mem use, gather from genes.cds seq headers? ie common prefix in all >IDPnnnn
  .. ridprefix may not be helpful
  
 UPD21OCT new sam2covtab -cdsbam cds.bam -chrbam chr.bam, no readids 
 ?? still want old sam2covtab(cdsbam) for cds.cov stats?
 
=cut

sub sh_sam_chr_cds2covtab { #UPD21OCT of sh_sam2covtab
  my($OUTSH,$chrbam,$cdsbam,$idclass)=@_;
  # opt: $dtype not used now, expect cdsbam to be merge of cds,TE if have both
  
  return("") unless($chrbam and $cdsbam); # unless($cdsbam) do old sam2covtab, no readids
  
  my $name= basename($chrbam,'\.bam'); 
  my $outtab="$name.cc$VERS.covtab";
  my $outchrtab="$name.cc$VERS.chrtab";
  
  # NOTE: opts should be same for -icpu and -merge calls
  my $opts="-nodebug "; # stop tons of debug msg from iterates
  $opts.=" -crclassf=$idclass" if($idclass);

  # UPD21OCT: if cdsbam, always add -genetab option, default after tests? for gnodes_sam2covtab8a.pl
  #   $DOGENETAB required here
  my $partset= "";
  $opts.=" -genetab ";
  (my $outgtab=$outchrtab) =~ s/chrtab/genetab/;  
  $partset .= $outgtab . '.pt* ';
  
  # if($idclass and not $gidprefix) gnodes_sam2covtab8j.pl should make gidpre from idclass.tab ids
  $opts.=" -gidprefix=$gidprefix" if($gidprefix);
  $opts.=" -ridprefix=$ridprefix" if($ridprefix);
    
  my $script= <<"EOS";
# run_gnodes1_sam_chrcds2cov $outtab
if [ -s $outtab ]; then
  echo reusing $outtab $outchrtab;
else
  echo START_sam2covtab $outtab `date`
  i=0; while [ \$i -lt \$NCPU ]; do {
    $gnodes1_sam2cov -icpu \$i -ncpu \$NCPU $opts -out $outtab -cdsbam $cdsbam -bam $chrbam  &
    sleep 5; 
    i=\$((\$i + 1));
  } done
  
  wait
 
  $gnodes1_sam2cov  -merge $opts -out $outtab  -bam $chrbam 
  mv $outtab.pt* $outchrtab.pt* $partset parts/
  echo DONE_sam2covtab  `date`
fi

EOS
  
  print $OUTSH $script;
  return($outtab); # ret($outtab,$outchrtab,$outgtab) ?
}


# run_gnodes_sam2covtab.sh script ============================

# pt=daphplx17evgt1m_cds; 
# pt=dplx20maca1r_tefams; 
# env opts="-savereadids -crclassf=dplx20cdste.idclass -minident=0.45"  bam=${pt}_SRR3090572_1_bwa.bam  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh
# env merge_readids="daphplx17evgt1m_cds_SRR3090572_1_bwa.readids dplx20maca1r_tefams_SRR3090572_1_bwa.readids"  cdstab=dplx20cdste_SRR3090572_1_bwa.readids  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh
# env opts="-crclassf=dplx20cdste.idclass" cdstab=dplx20cdste_SRR3090572_1_bwa.readids  bam=daphplx_gasm16ml_SRR3090572_1_bwa.bam  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh

sub sh_sam2covtab {
  my($OUTSH,$dtype,$ibam,$idclass,$readids)=@_;

  return("") unless($ibam);

  # FIXME readids for _1/_2 pair parts need to be uniq, add .1/.2 in sam2cov? when dupl RDID found after interval of others
  #UPD21apr20: FIXME this way causes readids from _1/_2 to be same, spurious duplicates, one bug is that
  #   cdsasm count now ~2x of before, and 2x of CDSann when should be nearly same. Other bugs not yet found.
  #   solution? add .1/.2 suffix to read ids? maybe ok but ID number is used as integer. Change IDprefix?
  # .. cant easily fix here w/o edit read input files.. fix in sam2covtab ?
  
  my $name= basename($ibam,'\.bam'); 
  my $outtab="$name.cc$VERS.covtab";
  my $outchrtab="$name.cc$VERS.chrtab";
  
  # NOTE: opts should be same for -icpu and -merge calls
  my $opts="-nodebug "; # stop tons of debug msg from iterates
  $opts.=" -crclassf=$idclass" if($idclass);

  # UPD21JUN: add -genetab option, default after tests? for gnodes_sam2covtab8a.pl
  # $opts.=" -genetab " if($DOGENETAB);
  my $partset= "";
  my $dothisgenetab= ($DOGENETAB and $dtype =~ /chr/); # not ($dtype =~ m/cds|te/)
  if($dothisgenetab){
    $opts.=" -genetab "; #?? only for dtype chr
    (my $outgtab=$outchrtab) =~ s/chrtab/genetab/;  
    if($dtype =~ /chr/){ $partset .= $outgtab . '.pt* '; }
  }

  if(UPD21OCT) { ## UPD21OCT: -savereadids/-readids obsolete for gnodes_sam2covtab8j.pl
    $readids=""; # return empty
  } else { 
    if($dtype =~ m/cds|te/){ 
      $readids="$name.readids";  $partset .= $readids . '.pt* ';
      $opts.=" -savereadids=$readids"; # was  -minident=0.45 << now too high, make option?
      }  
    elsif($dtype =~ /chr/) {
      $opts.=" -readidtab=$readids " if($readids); 
    }
  }
    
  #UPD21AUG: insert 'sleep 5;' in icpu loop to give system a bit of time in cases of huge mem overflows, like pig genome
  #Note: sam2cov -readidtab=bigdata.rids can be memory pig: 6GB for pig genome =~ 24 GB mem use per icpu task;
  # ** should revise this to add cds,te readid classes, gene ids into chr.bam as pseudo-read-group tags,
  # .. an extra pre-process but then sam2cov can avoid big mem structure for readids
  # eg. gnodes_samaddrids -readidtab=cdste.readidtab -bam chr.bam -outbam chr_cdsan.bam
  
  my $script= <<"EOS";
# run_gnodes1_sam2cov
if [ -s $outtab ]; then
  echo reusing $outtab $outchrtab;
else
  echo START_sam2covtab  `date`
  i=0; while [ \$i -lt \$NCPU ]; do {
    $gnodes1_sam2cov -icpu \$i -ncpu \$NCPU $opts -out $outtab -bam $ibam  &
    sleep 5; 
    i=\$((\$i + 1));
  } done
  
  wait
 
  $gnodes1_sam2cov  -merge  $opts -out $outtab  -bam $ibam 
  mv $outtab.pt* $outchrtab.pt* $partset parts/
  echo DONE_sam2covtab  `date`
fi

EOS
  
  print $OUTSH $script;
  return($outtab,$readids);
}

sub sh_merge_readids {
  my($OUTSH,$mergedids,@readid_parts)=@_;

  my $script=<<"EOS";
# run_gnodes1_merge_rids
if [ -s $mergedids ]; then
  echo reusing $mergedids ;
else
  echo START_gnodes1_merge_rids  `date`
  $gnodes1_sam2cov -merge readids -savereadids=$mergedids @readid_parts
  echo DONE_gnodes1_merge_rids  `date`
fi
EOS
  
  print $OUTSH $script;
  return($mergedids);
}

sub sh_merge_covtabs {
  my($OUTSH,$mergedcov,@covtab_parts)=@_;
  my $script="# run_gnodes1_merge_readids\n";
  
  my ($mergedchr,@chrparts)= map{ my $p=$_; $p=~s/.covtab//; $p.=".chrtab"; $p; } ($mergedcov,@covtab_parts);
  $script .= "$gnodes1_sam2cov -merge covtab -output=$mergedcov @covtab_parts\n";
  $script .= "$gnodes1_sam2cov -merge chrtab -output=$mergedchr @chrparts\n";
  
  print $OUTSH $script;
  return($mergedcov,$mergedchr);
}



#===== run_gnodes_annot.sh script =====================

sub sh_annot {
  my($OUTSH,$chr,$cds,$te,$idclass,$buscotsv,$genecovtab)=@_;
 
  return("") unless($chr and ($cds or $te));
  my $annout= basename($chr,'\.\w+$') . ".anntab"; 
  
  my $opts="";
  $opts.= ($debug)?" -debug":" -nodebug";
  $opts.=" -cds $cds" if($cds);
  $opts.=" -te $te" if($te);
  $opts.=" -idclass $idclass" if($idclass);
  $opts.=" -busco $buscotsv" if($buscotsv);
  $opts.=" -genecovtab $genecovtab" if($genecovtab);

  my $script=<<"EOS";
if [ -s $annout ]; then
  echo reusing $annout ;
else
  echo START_gnodes2_annotate  `date`
  $gnodes2_ann $opts -ncpu \$NCPU -chr $chr -output $annout
  echo DONE_gnodes2_annotate  `date`
fi

EOS

  print $OUTSH $script;
  # fixme: crclassf is output of sh_annot(): ficari18tsa1cds.fa > ficari18tsa1cds.idclass
  unless($idclass or not $cds){  $idclass=$cds; $idclass=~s/\.\w+$//; $idclass.=".idclass"; }
  return($annout,$idclass);  # idclass is output also, if not input
}

=item rmask.cat to te.seq

466 27.40 0.48 0.00 NC_037638.1 4243 4450 (27749750) Gypsy-24_DEl-I#LTR/Gypsy 2601 2809 (3339) m_b1s001i1

262 27.59 0.00 0.00 NC_037638.1 20221 20275 (27733925) C MuDR-1_LHu#DNA/MULE-NOF (635) 3923 3869 m_b1s001i6
   (27733925) C MuDR-1_LHu#DNA/MULE-NOF 
              ^ dingbat  extra colm for strand complement
              
 env idp=apimel perl -ne 'if(/^\d+/){ putfa() if($id and $fa);  @v=split;
 $cor=0; if($v[8] eq "C") { ($cor)= splice(@v,8,1); }
 ($cr,$cb,$ce,$te)=@v[4,5,6,8]; ($td,$tc)=split"#",$te;
 $cl=($tc=~/Unknown/)?"UNK":($tc=~/repeat|Sattelite/)?"repeat"
   :($tc=~/LTR|LINE|SINE|DNA|RC/)?"TE":"unk"; $cr=~s/\.\d//;
 $id=join"_",$cr,$cb,$ce;  $fa=""; $infa=1; } 
 elsif($infa){ ($d)=@v=split; if($d=~/^$cr/) { $s=$v[2]; $s=~s/\W//g; $fa.=$s; } } 
 END{ putfa(); } BEGIN{ $IDP=$ENV{idp}||"nada"; $IDP.="_rm"; }
 sub putfa{ $w=length($fa); $fa=~s/(.{100})/$1\n/g;
 print ">$IDP$id type=$cl; len=$w; tecl=$tc; teid=$td;\n",$fa,"\n"; } ' \
  apismel19hav31asm.fa.align > apismel19tefam.fa

=item parse_repmaskout : see gnodes_annotate.pl

=cut
 
sub sh_repeatmask {
  my($OUTSH,$chrasm,$species)=@_;
 
  return("") unless($chrasm); # guess at species for repmask? default?
  my $rmoutdir="./rmout"; my $cname=basename($chrasm);
  my $annout= $rmoutdir ."/$cname.out"; 
  my $rmlog = $rmoutdir ."/$cname.rmlog"; 

  # note rmkr appends to full chrasm name
  # want teseq.fasta output also, and te.gff, from rm -a align out?
  # rmout/
  
  # my $teseq="sss";
  # my $tegff="ggg";
  #.. rm.outs .. dont need fa.masked; use fa.cat to get align TE.seq ?
  # daphcari20chr.fa.cat.gz  daphcari20chr.fa.masked  daphcari20chr.fa.out	daphcari20chr.farm.gff	daphcari20chr.fa.tbl
  # daphcari20chr.farm.gff << bad name
  
  my $opts=" -q -xsmall";
  $opts.=" -species $species" if($species);

  my $script=<<"EOS";
if [ -s $annout ]; then
  echo reusing $annout ;
else
  echo START_repeatmasker  `date`
  mkdir $rmoutdir;
  RepeatMasker $opts -pa \$NCPU -dir $rmoutdir $chrasm > $rmlog 2>&1
  echo DONE_repeatmasker  `date`
fi

EOS

  print $OUTSH $script;
  return($annout);
}

sub sh_mrna2cds {
  my($OUTSH,$mrna,$cdsseq)=@_;
  
  unless($cdsseq) { # want .fa suffix here, not .cds
    $cdsseq=$mrna; $cdsseq =~ s/\.\w+$//; $cdsseq=~s/mrna//; $cdsseq.="cds.fa";  
  }
  (my $aaseq=$cdsseq) =~ s/\.\w+$//; $aaseq.=".aa";  # aa name= mygenes_cds.aa ?
  
  my $script=<<"EOS";
if [ -s $cdsseq ]; then
  echo reusing $cdsseq ;
else
  echo START_mrna2cds  `date`
  \$EVIGENES/cdna_bestorf.pl -nostop -cdna $mrna -outcds $cdsseq -outaa $aaseq
  echo DONE_mrna2cds  `date`
fi

EOS

  print $OUTSH $script;
  return($cdsseq,$aaseq);
}

=item buscoscan

# add sh_busco() ?? ordb=ordb_lineage
# run_BUSCO.py  -i $prot -l $ordb -o $tag$oname -m prot  --cpu $ncpu

try1:
  ordb_lineage = arthropoda from metad, not good enough for buscan
  
START_buscoscan Fri Feb 12 21:41:56 EST 2021
WARNING The dataset you provided does not contain the file dataset.cfg, likely because it is an old version. Default species (fly, eukaryota) will be used as augustus species
ERROR   Impossible to read arthropoda/   <<
DONE_buscoscan Fri Feb 12 21:42:31 EST 2021
START_gnodes2_annotate Fri Feb 12 21:42:31 EST 2021

=cut

sub sh_buscoscan {
  my($OUTSH,$cdsseq,$ordb_lineage)=@_;
 
  return("") unless($cdsseq); # guess at species for repmask? default?

  (my $aaseq=$cdsseq) =~ s/\.\w+$//; $aaseq.=".aa"; #  .cds chop bad: daphsim17evgt1cds.fa > daphsim17evgt.fa
  # unless(-f $aaseq) { }
  ## run_BU.py adds "run_" prefix
  my $buname= basename( $aaseq, ".aa");
  my $annotr= "run_" . $buname . "/full_table_" . $buname . ".tsv"; 
  my $buhmmer= "run_" . $buname . "/hmmer_output"; 
  #^^ fixme want this copied above run_xxx/ subdir, as? busco_full_table_$buname.tsv 
  # ficari18tsa1cds_buscofull_table.tsv better?
  # my $annout= "busco_full_table_".$buname.".tsv";
  my $annout=  $buname."_buscofull_table.tsv";
  
  ## dont set, busco is *supposed* to pick default eukaryota w/o -l option
  #nogo# $opts.=" -l $ordb_lineage" if($ordb_lineage);
  #this fixup may only be for IU HPC install w/ poor config file, lineage_path is correct but ;commented out
  my $optsh="";
  unless($ordb_lineage) { # missing opt bad = fail busco.sh
    $ordb_lineage="eukaryota"; # this is wrong likely
    #look at $metavals->{$asmid}{buscodb} and ENV{buscodb}
    warn "#*** WARN: Missing busco lineage option, buscodb=?, in metadata, using buscodb=$ordb_lineage *** \n";
  }
  if($ordb_lineage) {
    if($ordb_lineage =~ m,/,) { 
      $optsh="opts=-l $ordb_lineage"; 
    } else {
      $optsh=<<"EOS";
if [ X != X\$BUSCO_CONFIG_FILE ]; then
  ORDB=`grep lineage_path \$BUSCO_CONFIG_FILE | sed 's/^.*= *//;'`;
fi
if [ X != X\$ORDB ]; then opts="-l \$ORDB/$ordb_lineage"; fi
EOS
      
    }
  }

  my $script=<<"EOS";
if [ -s $annout ]; then
  echo reusing $annout ;
else
  echo START_buscoscan  `date`
  if [ ! -f $aaseq ]; then \$EVIGENES/cdna_bestorf.pl -nostop -cdna $cdsseq -outaa $aaseq ; fi
  opts=""
  $optsh
  
  run_BUSCO.py \$opts --cpu \$NCPU -m prot -i $aaseq -o $buname 
  if [ -s $annotr ]; then cp -p $annotr $annout; rm -r $buhmmer; fi
  echo DONE_buscoscan  `date`
fi

EOS

  print $OUTSH $script;
  return($annout);
}

#===== run_gnodes_ucgcov.sh script =====================


sub sh_ucgcov {
  my($OUTSH, $inbams, $covtab, $lidclassf, $norecalc)=@_;  
 
  # note: inbams can be list "a.bam b.bam c.bam .."
  if($inbams =~ m/,/){ $inbams=~s/,/ /g; }
  my($inbam1,@inbam2)= split" ",$inbams;  
  unless($covtab) {
    $covtab= $inbam1; $covtab=~s/\.bam//; $covtab=~s/.bwa//; 
    if(@inbam2){ $covtab .="n".(1+scalar(@inbam2)); }
    $covtab.="_ucg.covtab";
  }
  my $opts="";
  $opts.=" -idclassf $lidclassf" if($lidclassf);
  #x $opts.=" -asmid $lasmid" if($lasmid);
  #x $opts.=" -title $lintitle" if($lintitle);
  # $opts.=" -sumdata $sampledata" if($sampledata);
  $opts.= ($debug)?" -debug":" -nodebug";

  my $testsum= ($norecalc and $covtab)? "-a ! -s $covtab " : "";
  
  my $script=<<"EOS";
if [ -s $inbam1 $testsum]; then  
  echo START_ucgcov  `date`
  $gnodes3_ucgcov $opts -output $covtab -bam $inbams
  echo DONE_ucgcov  `date`
fi

EOS

  print $OUTSH $script;
  return($covtab);
}

sub sh_sam2genecov {
  my($OUTSH, $inbams, $covtab, $lidclassf, $norecalc)=@_;  
 
  # note: inbams can be list "a.bam b.bam c.bam .."
  if($inbams =~ m/,/){ $inbams=~s/,/ /g; }
  my($inbam1,@inbam2)= split" ",$inbams;  
  unless($covtab) {
    $covtab= $inbam1; $covtab=~s/\.bam//; $covtab=~s/.bwa//; 
    if(@inbam2){ $covtab .="n".(1+scalar(@inbam2)); }
    $covtab.=".genexcopy"; #?? what suffix
  }
  my $opts="";
  $opts.=" -idclassf $lidclassf" if($lidclassf);
  #x $opts.=" -asmid $lasmid" if($lasmid);
  #x $opts.=" -title $lintitle" if($lintitle);
  # $opts.=" -sumdata $sampledata" if($sampledata);
  $opts.= ($debug)?" -debug":" -nodebug";

  my $testsum= ($norecalc and $covtab)? "-a ! -s $covtab " : "";
  
  my $script=<<"EOS";
if [ -s $inbam1 $testsum]; then  
  echo START_genecov  `date`
  $gnodes3_sam2genecov $opts -output $covtab -bam $inbams
  echo DONE_genecov  `date`
fi

EOS

  print $OUTSH $script;
  return($covtab);
}

#===== run_gnodes_sum.sh script =====================
# env asmid=daphplx_gasm16ml covtab=daphplx_gasm16ml_SRR3090572_1_bwa.cdschr7b.covtab anntab=daphplx_gasm16ml.anntab asmdata=daphplx20chrs.metad  datad=`pwd` qsub -q debug run_gnodes_sum.sh
# $runapp $opts -asmid $asmid -anntab $anntab -sumdata $asmdata  -title $title  $covtab

sub sh_covsum {
  my($OUTSH, $covtab, $lasmid, $lintitle, $lanntable, $lidclassf,$genescovtab,$norecalc)=@_;  
  #?? want this or not? eg not for cds.covsum
  #NO: $lasmid||= $asmid; $lintitle||=$intitle; $lanntable||=$anntable; $lidclassf||=$idclassf; #global defs
  if($genescovtab eq "norecalc") {
    $norecalc= $genescovtab; $genescovtab="";
  }
  my $opts="";
  $opts.=" -asmid $lasmid" if($lasmid);
  $opts.=" -title $lintitle" if($lintitle);
  $opts.=" -anntab $lanntable" if($lanntable);
  $opts.=" -crclass $lidclassf" if($lidclassf);
  $opts.=" -genexcopy $genescovtab" if($genescovtab);
  $opts.=" -sumdata $sampledata" if($sampledata);
  $opts.= ($debug)?" -debug":" -nodebug";

  ##UPD21OCT add output name: -output xxx and return for others
  my $outsum= $covtab; $outsum=~s/\.\w+$/_sum.txt/;
  # old: chick19nc_chr_SRR3954707_b2_test8j.covtab > chick19nc_chr_test8j_SRR3954707_sum.txt
  # new: chick19nc_chr_SRR3954707_b2_test8j.covtab > chick19nc_chr_SRR3954707_b2_test8j_sum.txt
  
  #o: my $testsum= ($norecalc and $lintitle)? "-a ! -s ${lintitle}_sum.txt " : "";
  my $testsum= ($norecalc)? "-a ! -s $outsum " : "";

  my $script=<<"EOS";
if [ -s $covtab $testsum]; then  
  echo START_covsum  `date`
  $gnodes3_covsum $opts -output $outsum $covtab
  echo DONE_covsum  `date`
fi

EOS

  print $OUTSH $script;
  return($outsum);
}

=item sh_chrcovplot UPD21AUG

  gnodes_covsum.pl -plotchr -asmid $asmid  -title $outprefix -metad $sampledata
       -genexcopy $genescovtab -anntab $anntable -crclass $crclassf   $chrcovtab
       
  my($chrplots)= sh_chrcovplot($OUTSH, $chrcovtab, $asmid, $intitle, $anntable, $crclassf, $genescovtab);  

  see evigene/scripts/genoasm/gnodes_covsum2p.pl
  aweed20gnodes/
  pt=arath18tair_chr; pt=arath20max_chr; 
  env kucg=33 $evigene/scripts/genoasm/gnodes_covsum2p.pl -plotchr -debug \
   -genexcopy arath18tair1cds_SRR10178325_b2_mim.geneycopy -asmid $pt  -title ${pt}_testplota -anntab $pt.anntab  \
   -crclass arath18tair1cds.idclass -sumdata arath20asm.metad   ${pt}_SRR10178325_test8f.covtab
    
  dromel20gnodes/
  pt=drosmel6ref_chr; pt=drosmel20pi_chr;  
  env kucg=0 $evigene/scripts/genoasm/gnodes_covsum2p.pl -plotchr -debug \
    -genexcopy dromel6relt1cds_SRR11460802_b2_mim.genexcopy -asmid $pt  -title ${pt}_testplota \
   -anntab $pt.anntab  -crclass dromel6relt1cds.idclass -sumdata drosmelchr.metad ${pt}_SRR11460802_b2_mim.cc8a.covtab

=cut

sub sh_chrcovplot {
  my($OUTSH, $covtab, $lasmid, $lintitle, $lanntable, $lidclassf, $genescovtab, $norecalc)=@_;  

  my $opts="-plotchr";
  $opts.=" -asmid $lasmid" if($lasmid);
  $opts.=" -title $lintitle" if($lintitle); #? need change to output title for chrplot? or not > title.plottab default out
  $opts.=" -anntab $lanntable" if($lanntable);
  $opts.=" -crclass $lidclassf" if($lidclassf);
  $opts.=" -genexcopy $genescovtab" if($genescovtab);
  $opts.=" -sumdata $sampledata" if($sampledata);
  $opts.= ($debug)?" -debug":" -nodebug";

  my $testsum= ($norecalc and $lintitle)? "-a ! -s ${lintitle}.plottab " : "";

  # NOTE covsum -plotchr also will produce Rscript for drawing chr plots, call should be added in this script
  # my $rplot="${lintitle}.plotchr.sh"; 
  # my $addscript="if [ -x  $rplot ]; then ./$rplot; fi";
  ## dang name bug: gnodes_covsum makes this: chick19nc8f_SRR3954707.plotchr.sh NOT chick19nc8f.plotchr.sh
  ## ** gnodes_covsum should execute rplot.sh, check for Rscript?
  
  my $script=<<"EOS";
if [ -s $covtab $testsum]; then  
  echo START_chrcovplot  `date`
  $gnodes3_covsum $opts $covtab
  echo DONE_chrcovplot  `date`
fi

EOS

  print $OUTSH $script;
  return(1);
}

=item sh_sumgenescov UPD21JUL

  my $gnodes4_sumgenecov="\$EVIGENES/genoasm/gnodes_sumgenecov.pl"; #upd21jul : replaces gnodes_genescov

  pt=daphpulex_pa42v2; 
  $evigene/scripts/genoasm/gnodes_sumgenecov.pl -title ${pt}_test8f7sumgcn
  -transform=0  -genexcopy daphplx17evgt1m_cds_SRR13333791_b8_mim.geneycopy
  -chrgenetab ${pt}_SRR13333791_test8f.genetab   -asmid $pt -meta
  daphplx20chrs.metad  -idclass daphplx17evgt1m_cds.idclass

  pt=arath18tair_chr; 
  $evigene/scripts/genoasm/gnodes_sumgenecov.pl -title at18chr_test8f7sumgcn
  -transform=0  -genexcopy arath18tair1cds_SRR10178325_b2_mim.geneycopymd
  -chrgenetab ${pt}_SRR10178325_test8f.genetab    -asmid $pt -meta
  arath20asm.metad -idclass arath18tair1cds.idclass

 outputs title.sumgenetab << rename? _genesum.tab? _genecopy.sumtab?
      and title_genesum.txt << rename? _genecopy.summary.txt ?
 
  ==> aweed20gnodes/at18chr_test8f7sumgcn_genesum.txt <==
  # Chromosome Assembly Gene-Copynum Counts by Copy level
  # N Genes valid=27347, zero-cover=97, all=27444 (zero:gClass=zero or gNread<1)
   ---------------------------------------
   Genes Copynum Summary by Copy level, for nGene=27347
   aCopy=Chrasm, gCopy=Gene set, xCopy=aCopy/gCopy, xLo,xEq,xHi= xCopy<1,=1,>1
   aMiss=Missed gene reads on chrasm, % of gene reads
   CLevel   nGene   gCopy   xCopy   aCopy   aMiss           xLo,xEq,xHi
   All      27347   1.3     1.0     1.1     0.0043%,825     347,26597,403
   0        262     0.5     2.2     1.1     0.017%,4        1,0,261
   1        26200   1.0     1.0     1.0     0.00066%,88     0,26098,102
   2-9      806     3.1     0.8     2.3     0.0043%,32      319,458,29
   10-99    68      21.0    1.0     17.5    0.0099%,314     16,41,11
   99-499   8       307.2   0.1     25.9    0.022%,387      8,0,0
   500+     3       650.9   0.0     2.1     0%,0            3,0,0
   ---------------------------------------
             .. And ..
  # Moments for all genes
  Field   Median  Mean    N       StDev   Skew    Range   Correlations
  gLen    942     1351    25322   1984    0.619   120,66729
  # Moments for gCopy => [ 1.55 99.9 ] genes
  # Moments for gCopy => [ 0.66 1.55 ] genes
  # Moments for gCopy => [ 0 0.66 ] genes
 
=cut

sub sh_sumgenescov8j {
  my($OUTSH, $genexcopy, $chrgenetab, $chrsum, $lintitle, $cdscovtab)=@_;   
  
  #UPD21OCT: gnodes_sumgenecov8j.pl
  # gnodes4_sumgenecov  -title $title -genexcopy $genex -chrsum $chrsum -chrgenetab $genetab   
  #  genexcopy = gene copynum table,  primary data from gnodes_sam2genecov
  #  chrgenetab = chr x gene-cds cov table,  primary data, from updated sam2covtab:putGenetab() : GeneID ChrID Pos aCovT aCovM noCov
  #  chrsum = chr_sum.txt from covsum prior call
  #UPD21OCT31: add -cdscov xxx_cds.covtab
  
  my $opts="";
  $opts.=" -title $lintitle" if($lintitle);
  $opts.=" -genexcopy $genexcopy" if($genexcopy); # required from gnodes_sam2genecov.pl
  $opts.=" -chrgenetab $chrgenetab" if($chrgenetab); # required from gnodes_sam2covtab8e.pl
  $opts.=" -chrsum $chrsum" if($chrsum); # wanted from gnodes_covsum.pl : _sum.txt table of parts cover
  $opts.=" -cdscov $cdscovtab" if($cdscovtab); # 

  $opts.= ($debug)?" -debug":" -nodebug";
  #drop: $opts.=" -asmid $lasmid" if($lasmid);
  #drop: $opts.=" -idclass $lidclassf" if($lidclassf);
  #drop: $opts.=" -metadata $sampledata" if($sampledata);
  #drop: $opts.=" -anntab $lanntable" if($lanntable); #? drop

  my $testsum="";
  #? my $testsum= ($norecalc and $lintitle)? "-a ! -s ${lintitle}_genesum.txt " : "";

  my $script=<<"EOS";
if [ -s $chrgenetab $testsum ]; then  
  echo START_sumgenecov  `date`
  $gnodes4_sumgenecov $opts
  echo DONE_sumgenecov  `date`
fi

EOS

  print $OUTSH $script;
  return(1);
}

sub sh_sumgenescov {
  my($OUTSH, $genexcopy, $chrgenetab, $lasmid, $lintitle,  $lidclassf, $lanntable)=@_;  #, $norecalc
  
  #  genexcopy = gene copynum table,  primary data from gnodes_sam2genecov
  #  chrgenetab = chr x gene-cds cov table,  primary data, from updated sam2covtab:putGenetab() : GeneID ChrID Pos aCovT aCovM noCov
  #  idclass = gene ids w/ classes: CDS/TE/busco/UCG/other
  #  anntab = chr x annots, ?? not needed w/ gene idclass tab, all prime data is gene-id based
  
  my $opts="";
  $opts.=" -asmid $lasmid" if($lasmid);
  $opts.=" -title $lintitle" if($lintitle);
  $opts.=" -genexcopy $genexcopy" if($genexcopy); #? required from gnodes_sam2genecov.pl
  $opts.=" -chrgenetab $chrgenetab" if($chrgenetab); # required from gnodes_sam2covtab8e.pl
  $opts.=" -idclass $lidclassf" if($lidclassf);
  $opts.=" -metadata $sampledata" if($sampledata);
  $opts.=" -anntab $lanntable" if($lanntable); #? drop
  $opts.= ($debug)?" -debug":" -nodebug";

  my $testsum="";
  #? my $testsum= ($norecalc and $lintitle)? "-a ! -s ${lintitle}_genesum.txt " : "";

  my $script=<<"EOS";
if [ -s $chrgenetab $testsum ]; then  
  echo START_sumgenecov  `date`
  $gnodes4_sumgenecov $opts
  echo DONE_sumgenecov  `date`
fi

EOS

  print $OUTSH $script;
  return(1);
}

sub sh_genescov {  # obsolete for UPD21JUN
  my($OUTSH, $cdscov, $chrcov, $lasmid, $lintitle, $lanntable, $lidclassf,$norecalc)=@_;  
  
  # NOTE cdscov or chrcov may be empty, esp cdscov
  # NOTE title change for genescov? mix of cds+chr names
  # DEF out: xxx_genesum.txt, xxx.genemeans xxx.genexcopy
  
  my $opts="";
  $opts.=" -asmid $lasmid" if($lasmid);
  $opts.=" -title $lintitle" if($lintitle);
  $opts.=" -anntab $lanntable" if($lanntable);
  $opts.=" -crclass $lidclassf" if($lidclassf);
  $opts.=" -sumdata $sampledata" if($sampledata);
  $opts.= ($debug)?" -debug":" -nodebug";

  my $testsum="";
  #? my $testsum= ($norecalc and $lintitle)? "-a ! -s ${lintitle}_genesum.txt " : "";
  
  $opts.= " -cdscov $cdscov" if($cdscov);
  $opts.= " -chrcov $chrcov" if($chrcov); # assume this one is set
  
  my $script=<<"EOS";
if [ -s $chrcov $testsum]; then  
  echo START_genescov  `date`
  $gnodes4_genescov $opts
  echo DONE_genescov  `date`
fi

EOS

  print $OUTSH $script;
  return(1);
}

