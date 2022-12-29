#! /bin/bash
### env rnain=bams/xxx.bam qsub -q normal bamstats.sh
#PBS -N bamintr
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=2:55:00
#PBS -o bamintr.$$.out
#PBS -e bamintr.$$.err
#PBS -V

# takes only few minutes/25Mill bam
ncpu=16

# env keepids=0 table=1 name=myrna nreads=123456000 TMPDIR=/var/tmp $sdir/samstats.pl
#no# opts="keepids=0 table=1"
opts=""

datad=$HOME/scratchn
workd=$datad/chrs/kfish/rnas/
bindir=$HOME/bio/bin
sdir=$HOME/bio/evigene/scripts/rnaseq

cd $workd/
drnabam=$rnain
outdir=$workd/bamstats
mkdir $outdir
export TMPDIR=$outdir
echo "start bamintr : `date`"  
i=0; j=1; 
for bam in $drnabam ;  do { 

  nam=`basename $bam .bam`
  #tophat# nam=`echo $bam | sed 's,/accepted_hits.bam,,;'`
  echo "$bam TO $nam.$j.intron.gff " 
  ( $bindir/samtools view $bam | env $opts name=$nam $sdir/samstats.pl > $outdir/$nam.sstats ; ) &
  
  j=$(( $j + 1 ))
  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi
} 
done
wait

echo "end bamintr: `date` " 

