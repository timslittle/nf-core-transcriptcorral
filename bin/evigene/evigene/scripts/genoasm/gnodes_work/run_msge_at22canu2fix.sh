#! /bin/bash
## run_msge_at22canu2fix.sh : update for bad minimap -x sr, replace bedtool/Gcov w/ samtool depth
## --- gnodes_setup.sh for quartz.iu ---    
#SBATCH --job-name="gnodes_pipe"
#SBATCH --output="gnodes_pipe.%j.log"
#.. SBATCH --partition=general
#SBATCH --partition=debug
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=240G
#SBATCH --export=ALL

export NCPU=32; export TMEM=120G; MEM=3G;
# align long pacbio, test mgse cov, genosize stat

datad=/N/slate/gilbertd/chrs/daphplx/gasm20set/aweed20gnodes
# chrfa=hetzdata/at22vlr_ncbi_chrs.fa
chrfa=iatdata/at22vcan2d_chrs.fa

## pacbhifi
# reads=readvl/SRR16841688.fasta
# crname=at22vlr_ncbi_chrs_SRR16841688ph
# mimopt="-x map-pb"
## salk21 data
# reads=read22f/SRR14654461_[12].fastq.gz
# crname=at22vlr_ncbi_chrs_SRR14654461ab
# mimopt="-x sr"
## tsukuba21 data
reads=read22f/DRR214840_[12].fastq.gz
readb=read22f/DRR214842_[12].fastq.gz
crname=at22vcan2d_chrs_DRR214840cd
# mimopt="-x sr"
mimopt="-x sr --secondary=yes -p 0.95 -N 9999999"
# sr      Short single-end reads without splicing (-k21 -w11 --sr --frag=yes -A2 -B8  -O12,32
#   -E2,1  -r50 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g200 -2K50m --heap-sort=yes --secondary=no).

# at22vcan2d_dug.chrbad.list  at22vcan2d_dug.dupgenes.tsv  at22vcan2d_dug.ucgenes.tsv
refgenes=at22vcan2d_dug.ucgenes.tsv
dupgenes=at22vcan2d_dug.dupgenes.tsv
badlist=at22vcan2d_dug.chrbad.list
odir=mgse
name=$crname

gcovbed=$crname.stdcov.txt
dupcovbed=$crname.sdupcov.txt

if [ X = "X$dup" ]; then dodup=0; else dodup=$dup; fi 
if [ X = "X$datad" ]; then echo "datad=what?"; exit -1; fi

nbin=$HOME/bio/ncbi/bin; export PATH=$nbin:$PATH;
sambin=$HOME/bio/apps/bin; export PATH=$sambin:$PATH;

module load python
bwabin=$HOME/bio/bwa/bin/; export PATH=$bwabin:$PATH; # export RDMAPPER=bwa-mem2
mimbin=$HOME/bio/gnomutil/flye28/bin;  export PATH=$mimbin:$PATH; # export RDMAPPER=minimap2;

mgsepy=$HOME/bio/gnomutil/mgse19f/MGSE.py

## drop bedapp for samtools depth
# bedbin=/N/soft/rhel7/bedtools/gnu/2.26.0/bin
# bedapp=genomeCoverageBed

export EVIGENES=$HOME/bio/apps/evigene/scripts
export PATH=$PATH:$EVIGENES
cd $datad

if [ ! -s ${crname}_mim.bam ]; then
echo START_dnamap  `date`
( minimap2 $mimopt -a -t $NCPU $chrfa $reads | samtools view --threads $NCPU -Sb -o ${crname}_mim.bam - ) > ${crname}_mim.log 2>&1

if [ ! -s ${crname}_mim.bam ]; then echo "ERR: minimap $chrfa $reads to  ${crname}_mim.bam"; exit -1; fi
echo DONE_dnamap `date`
fi

if [ ! -s $gcovbed ]; then
echo START genomeCoverageBed `date`

samtools sort -m $MEM --threads $NCPU -o ${crname}_msort.bam ${crname}_mim.bam

## -H header screws mgse.py, skip
sdopt="-aa --threads $NCPU"
## this way discards 2ndayr: -G UNMAP,SECONDARY,QCFAIL,DUP
samtools depth -G 0x100 $sdopt -o $gcovbed ${crname}_msort.bam 
## this way keeps 2nday
samtools depth -g 0x100 $sdopt -o $dupcovbed ${crname}_msort.bam 

# /bin/rm ${crname}_msort.bam
if [ ! -s $gcovbed ]; then echo "ERR: fail $gcovbed"; exit -1; fi
echo DONE genomeCoverageBed `date`
fi

echo START mgse `date`
oname=$crname.mgse
mopt=" --cov $gcovbed --ref $refgenes"
if [ $dodup = 1 ]; then
  mopt=" --cov $gcovbed --ref $dupgenes"
  oname=$crname.mdupx
  if [ ! -s $dupgenes ]; then echo "ERR: no dupgenes: $dupgenes\n"; exit -1; fi
fi
if [ -s $badlist ]; then mopt="$mopt --black $badlist"; fi

echo python $mgsepy $mopt --name $oname --out $odir 
python $mgsepy $mopt --name $oname --out $odir 

echo END mgse `date`

