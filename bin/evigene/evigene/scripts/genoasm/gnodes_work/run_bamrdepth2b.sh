#! /bin/bash
## run_sambamstats.sh : samtools flagstat / stats of in.bam
## --- gnodes_setup.sh for quartz.iu ---    
#SBATCH --job-name="gnodes_pipe"
#SBATCH --output="gnodes_pipe.%j.log"
#.. SBATCH --partition=general
#SBATCH --partition=debug
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=120G
#SBATCH --export=ALL

export NCPU=32; export TMEM=120G; MEM=3G;

# version 
dv=2c

if [ X = "X$datad" ]; then echo "datad=what?"; exit -1; fi
# if [ X != "X$idclass" ]; then dostd2genecov=1; fi

if [ X = "X$bam" ]; then 
  if [ X = "X$asm" ]; then echo "asm=what?"; exit -1; fi
  if [ X = "X$reads" ]; then echo "reads=what?"; exit -1; fi
  aname=`basename $asm | sed 's/\.[a-z0-9]*$//;'`
  rname=`basename $reads | sed 's/.gz//; s/\.[a-z0-9]*$//;'`
  bname=${aname}_${rname}_mim$dv
  bam=$bname.bam

else
  bname=`basename $bam .bam`
  bname=$bname$dv
fi

if [ X = "X$mimopt" ]; then 
  mimopt="-x sr --secondary=yes -p 0.95 -N 9"
  # mimopt="-x map-pb" for pachifi
fi
# mimopt long readtypes -x map-pb  map-ont ??
# sr      Short single-end reads without splicing (-k21 -w11 --sr --frag=yes -A2 -B8  -O12,32
#   -E2,1  -r50 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g200 -2K50m --heap-sort=yes --secondary=no).
# mimopt="-x sr"
# mimopt="-x sr --secondary=yes -p 0.95 -N 9999999"
#  -p 0.95 -N 9 works, huge N slows it down too much, dont need more than a few dups to measure dupx
# try -p 0.975 ??

gcovbed=$bname.stdcov.txt
dupcovbed=$bname.sdupcov.txt

nbin=$HOME/bio/ncbi/bin; export PATH=$nbin:$PATH;
sambin=$HOME/bio/apps/bin; export PATH=$sambin:$PATH;
# bwabin=$HOME/bio/bwa/bin/; export PATH=$bwabin:$PATH; # export RDMAPPER=bwa-mem2
mimbin=$HOME/bio/gnomutil/flye28/bin;  export PATH=$mimbin:$PATH; # export RDMAPPER=minimap2;

# module load python
# mgsepy=$HOME/bio/gnomutil/mgse19f/MGSE.py

export EVIGENES=$HOME/bio/apps/evigene/scripts
export PATH=$PATH:$EVIGENES

cd $datad

if [ ! -s ${bam} ]; then
echo START_dnamap  `date`
echo minimap2 $mimopt -a -t $NCPU $asm $reads
( minimap2 $mimopt -a -t $NCPU $asm $reads | samtools view --threads $NCPU -Sb -o $bam - ) > ${bname}.log 2>&1

if [ ! -s ${bam} ]; then echo "ERR: minimap $asm $reads to  ${bam}"; exit -1; fi
echo DONE_dnamap `date`
fi

echo START samdepth `date`

samtools flagstats --threads $NCPU ${bam}  > ${bname}.fstats.txt
samtools sort -m $MEM --threads $NCPU -o ${bname}_msort.bam ${bam}

# * Diffs w/ gnodes_sam2genecov due in part to read-map quality filter
#  sam2gc skips align < $MIN_IDENT * $RDLEN, MINID=0.40 def, ie align < 60bp for 150bp read
# test these samtools depth opts:
#  -J = count D dels, as valid read-align bases, as does sam2gc
#  -Q = ? 5 .. 10 .. ? min-mapq, RDLEN dependent??
#  -q = ? base-qual, from fastq.quals, 13 in mpileup, 0 in depth
#   .. The default minimum quality value is 0 for "depth" and 13 for "mpileup".
 
## -H header screws mgse.py, so what, add back *, use w/ stmdepth.tab, sdepthtmc.tab
# sdopt1a="-aa --threads $NCPU"
# sdopt2b="-aa -J -q 10 -Q 5 --threads $NCPU"
# ^^ -q 10 -Q 5 removes most dupx, all dups class, probably -Q 5 as dup/2ndary mostly have sam.MapQ==0
sdopt="-H -aa -J -q 10 -Q 0 --threads $NCPU"
# 2c, drop -Q (0 def), ?? try samt view opt for min align, using rdlen*0.40 ?
#   -m, --min-qlen INT         ...cover >= INT query bases (as measured via CIGAR)

## 1st discards 2ndayr: -G UNMAP,SECONDARY,QCFAIL,DUP; 2nd keeps 2ndary
echo samtools depth -G 0x100 $sdopt -o $gcovbed ${bname}_msort.bam
samtools depth -G 0x100 $sdopt -o $gcovbed ${bname}_msort.bam
echo samtools depth -g 0x100 $sdopt -o $dupcovbed ${bname}_msort.bam
samtools depth -g 0x100 $sdopt -o $dupcovbed ${bname}_msort.bam

# and paste:
paste  $dupcovbed $gcovbed | cut -f1,2,3,6 > ${bname}.stmdepth.tab
#^ add in $cdscovtab, with header-name-cleaned > sdepthtmc.tab

outs="${bname}.stmdepth.tab ${bname}.fstats.txt"
# only for cds.bam ...
if [ X != "X$idclass" -a -s $idclass ]; then

rdlen=`samtools view -F 0x104 ${bam} | head -1 | cut -f10 | wc -c`;
rdlen=$(($rdlen - 1));

echo $EVIGENES/genoasm/gnodes_std2genecov.pl \
 -idclass $idclass -rdlen $rdlen  \
 -in  ${bname}.stmdepth.tab  -stats  ${bname}.fstats.txt \
 -out $bname.stdt_genexcopy
 
$EVIGENES/genoasm/gnodes_std2genecov.pl \
 -idclass $idclass -rdlen $rdlen  \
 -in  ${bname}.stmdepth.tab  -stats  ${bname}.fstats.txt \
 -out $bname.stdt_genexcopy

 outs="$outs $bname.stdt_genexcopy"
fi

echo "# OUTPUTS: "  `ls $outs`
echo "# TMP remove:  rm ${bname}_msort.bam $dupcovbed $gcovbed "
echo DONE samdepth `date`

