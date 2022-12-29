#! /bin/bash
### env sratable=sraset.csv datad=`pwd` ncpu=16 qsub -q normal run_evgsra2genes.sh
#PBS -N evgsra2genes
#PBS -A PutAccountIdHere
#PBS -l nodes=1:ppn=16,walltime=39:55:00
#PBS -V

APPS_PATH=UDG
#. APPS_PATH=MODULE

evgapp=evgpipe_sra2genes4v.pl

if [ "X" = "X$ncpu" ]; then ncpu=8; fi
if [ "X" = "X$maxmem" ]; then maxmem=64000; fi
if [ "X" = "X$datad" ]; then echo "ERROR: missing datad=/path/to/data"; exit -1; fi
## allow dummy sratable/sraids for runsteps >= 7

if [ $APPS_PATH = "UDG" ]; then
  #ORIG# bioapps=/bio/apps
  #HOME# 
  bioapps=$HOME/bio/apps
  ## need for steps >= 7, gene set reduce, annotate, publish
  evigenes=$bioapps/evigene/scripts
  xnrbin=$bioapps/exonerate/bin
  cdhitbin=$bioapps/cdhit/bin
  ncbibin=$bioapps/ncbi/bin
  gmapbin=$bioapps/gmap/bin
  hmmerbin=$bioapps/hmmer/bin
  ## need for steps 2..5, rna data fetch, clean, assemble
  srabin=$bioapps/sratools/bin
  velobin=$bioapps/rnaseq/velvet_oases/bin  # fixme multi kmer binaries
  idbabin=$bioapps/rnaseq/idba/bin
  soapbin=$bioapps/rnaseq/soaptrans  # fixme binaries
  trinbin=$bioapps/rnaseq/trinity # fixme binaries
  kmerbin=$bioapps/khmer/scripts

  addpath=$srabin:$ncbibin:$velobin:$idbabin:$soapbin:$trinbin
  addpath=$evigenes:$xnrbin:$cdhitbin:$gmapbin:$hmmerbin:$kmerbin:$addpath
  export PATH=$addpath:$PATH
fi

if [ $APPS_PATH = "MODULE" ]; then
  module add blast sratoolkit
  module add evigene
  module add exonerate cdhit khmer
  module add gmap_gsnap hmmer
  module add velvet oases soaptrans idba trinity
fi

evopts="-NCPU $ncpu -MAXMEM $maxmem -log -debug"
## add opt -runstep start7, start at evgreduce step7
if [ "X" != "X$runsteps" ]; then evopts="$evopts -runsteps $runsteps"; fi
if [ "X" != "X$name" ]; then evopts="$evopts -runname $name"; fi
if [ "X" != "X$sratable" ]; then 
  evopts="$evopts -SRAtable $sratable"; 
elif [ "X" != "X$sraid" ]; then 
  evopts="$evopts -SRAids $sraid"; 
elif [ "X" != "X$runsteps" ]; then
  evopts="$evopts -SRAids=SRR000000";
else
  echo "ERROR: missing sratable=/path/to/data"; exit -1; 
fi

cd $datad
echo $evigenes/$evgapp  $evopts 
$evigenes/$evgapp  $evopts 

##--------------------------------------------
#i SRA2Genes Component applications currently used on PATH:
#i  path=/bio/apps/sratoolkit/bin/fastq-dump, https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
#i  path=/bio/apps/ncbi/bin/blastn,  https://blast.ncbi.nlm.nih.gov/
#i  path=/bio/apps/cdhit/bin/cd-hit-est,  https://github.com/weizhongli/cdhit/
#i  path=/bio/apps/exonerate/bin/fastanrdb, https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate
#i  path=/bio/apps/khmer/scripts/normalize-by-median.py, https://github.com/ged-lab/khmer 
#i  path=/bio/apps/ncbi/bin/vecscreen, http://ncbi.nlm.nih.gov/tools/vecscreen/
#i  path=/bio/apps/ncbi/bin/tbl2asn, http://ncbi.nlm.nih.gov/genbank/tbl2asn2/
#i  path=/bio/apps/rnaseq/velvet_oases/bin/velveth, https://www.ebi.ac.uk/~zerbino/oases/
#i  path=/bio/apps/rnaseq/idba/bin/idba_tran, https://code.google.com/archive/p/hku-idba/downloads/
#i  path=/bio/apps/rnaseq/soaptrans/SOAPdenovo-Trans-127mer, http://soap.genomics.org.cn/SOAPdenovo-Trans.html
#i  path=/bio/apps/rnaseq/trinity/Trinity, https://github.com/trinityrnaseq/trinityrnaseq
##--------------------------------------------
