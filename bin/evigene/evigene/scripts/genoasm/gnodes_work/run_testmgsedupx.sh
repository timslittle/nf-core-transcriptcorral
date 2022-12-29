#! /bin/bash
## run_testmgsedupx.sh
#SBATCH --job-name="gnodes_pipe"
#SBATCH --output="gnodes_pipe.%j.log"
#.. SBATCH --partition=debug
#SBATCH --partition=general
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=120G
#SBATCH --export=ALL

## test calc Cdupx instead of Cucg..
##  env dup=1 uniq=0 inbam=bamf/arath18tair_chr_SRR10178325_b2_bwa.bam anntab=arath18tair_chr_nocpmtdug.anntab datad=`pwd` sbatch run_testmgsedupx.sh

## samtools sort -m MemPerThread !
export NCPU=32; export TMEM=120G; MEM=3G;

if [ X = "X$inbam" ]; then echo "inbam=what?"; exit -1; fi
if [ X = "X$anntab" ]; then echo "anntab=what?"; exit -1; fi
if [ X = "X$datad" ]; then echo "datad=what?"; exit -1; fi
if [ X = "X$dup" ]; then dodup=0; else dodup=$dup; fi 

name=`basename $inbam .bam`
obam=$name.dedups.bam
gcovbed=$name.gcovbed.txt
ann=`basename $anntab .anntab`
refgenes=$ann.ucgenes.tsv

odir=mgse

module load python
sambin=$HOME/bio/apps/bin; export PATH=$sambin:$PATH;
## corrected by dgg: mgse19f/MGSE.py
mgsepy=$HOME/bio/gnomutil/mgse19f/MGSE.py
bedbin=/N/soft/rhel7/bedtools/gnu/2.26.0/bin
bedapp=genomeCoverageBed

export EVIGENES=$HOME/bio/apps/evigene/scripts
export PATH=$PATH:$EVIGENES

cd $datad
mkdir $odir

if [ ! -s $gcovbed ]; then 
echo START_gcovbed `date`

if [ ! -s $obam ]; then 
echo samtools view -F 0x100 $inbam -1 P samtools sort -m $MEM --threads $NCPU - TO $obam
samtools view -F 0x100 $inbam --threads $NCPU -1 | samtools sort -m $MEM --threads $NCPU - > $obam
if [ ! -s $obam ]; then echo samtools sort failed; exit -1; fi
fi

## mgse.py wants suffix .cov or .txt for this data, otherwise does gunzip !!!
## use  name.gcovbed.txt
echo $bedbin/genomeCoverageBed -d -split -ibam $obam  TO $gcovbed
$bedbin/genomeCoverageBed -d -split -ibam $obam > $gcovbed
if [ ! -s $gcovbed ]; then echo genomeCoverageBed failed; exit -1; fi

echo END_gcovbed `date`
fi

# make ann.ucgenes.gff  from anntab; collapse adjacent bins? if cb < lce+1 ..
# --ref = ref_genes.txt format is tsv of [ chr start end ID ]
refgenes=$ann.ucgenes.tsv
dupgenes=$ann.dupgenes.tsv
badlist=$ann.chrbad.list

if [ ! -s $refgenes ]; then

 grep busco $anntab | grep -v contam | \
 perl -ne '($cr,$cb,$ty,$ids)=split; $ce=$cb+99; $ids=~s/,.*//;  print join("\t",$cr,$cb,$ce,$ids)."\n"; ' \
  > $refgenes

 grep dupx $anntab | grep -v contam | \
 perl -ne '($cr,$cb,$ty,$ids)=split; $ce=$cb+99; $ids=~s/,.*//;  print join("\t",$cr,$cb,$ce,$ids)."\n"; ' \
  > $dupgenes

 grep contam $anntab | cut -f1 | sort -u > $badlist

fi

oname=$name.mgse
mopt=" --cov $gcovbed --ref $refgenes"
if [ $dodup = 1 ]; then
  mopt=" --cov $gcovbed --ref $dupgenes"
  oname=$name.mdupx
  if [ ! -s $dupgenes ]; then echo "ERR: no dupgenes: $dupgenes\n"; exit -1; fi
fi

if [ -s $badlist ]; then mopt="$mopt --black $badlist"; fi

echo START mgse `date`

echo python $mgsepy $mopt --name $oname --out $odir 
python $mgsepy $mopt --name $oname --out $odir 

echo END mgse `date`

