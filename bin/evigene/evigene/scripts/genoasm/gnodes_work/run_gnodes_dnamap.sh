#! /bin/bash
## run_gnodes_dnamap.sh = run_bwamem
#PBS -N gnodes_dnamap
#PBS -l vmem=124gb,nodes=1:ppn=24,walltime=1:00:00
#PBS -V

# update4: map only 1 reads by default, parallel map 2, 1/2 ncpu for speed test; drop paired= opt
ncpuall=24; ncpu=24; maxmem=124000;

if [ "X" = "X$datad" ]; then echo "ERR:datad=what?"; exit -1; fi
if [ "X" = "X$asm" ]; then echo "ERR:asm=what?"; exit -1; fi
if [ "X" = "X$reads" ]; then echo "ERR:reads=what?"; exit -1; fi

lreads=$reads
if [ "X" != "X$both" ]; then 
  rreads=`echo $reads | sed 's/_1/_2/;'`; 
  if [ $lreads = $rreads ]; then rreads=""; fi
fi

module load samtools
#older: module load bwa-mem2
runbin=$HOME/bio/bwa/bin/
export PATH=$runbin:$PATH
runapp=bwa-mem2

btag="_bwa"
aname=`basename $asm | sed 's/.gz//; s/\..*//; '`
if [ "X" = "X$name" ]; then 
  nab=`echo $reads | sed 's/.gz//; s/.fastq//; s/.fq$//; s/,.*//; s/ .*//;'`
  nab=`basename $nab`
  name=$aname"_"$nab;
fi

trasmidx="bwax/${aname}_bwx"
outbam=$name$btag.bam
rname=`echo $name | sed 's/_1//; s/$/_2/;'`
routbam=$rname$btag.bam

# use bwa -a == all, best to get all multimap reads, but lots of poor-2nd-maps
bopt="-a"
if [ "X" != "X$opt" ]; then bopt=$opt; fi

cd $datad/
echo "#START " `date`

if [ ! -f "$trasmidx.ann" ]; then
  mkdir bwax
  echo $runapp index -p $trasmidx $asm
  $runapp index -p $trasmidx $asm
fi

if [ "X" != "X$rreads" ]; then
  ncpu=$(( $ncpu / 2 )); 
  echo $runapp mem $bopt -t $ncpu -o $outbam $trasmidx $lreads 
  $runapp mem $bopt -t $ncpu $trasmidx $lreads | samtools view --threads $ncpu -Sb -o $outbam &
  echo $runapp mem $bopt -t $ncpu -o $routbam $trasmidx $rreads 
  $runapp mem $bopt -t $ncpu $trasmidx $rreads | samtools view --threads $ncpu -Sb -o $routbam &
  wait
else
  echo $runapp mem $bopt -t $ncpu -o $outbam $trasmidx $lreads 
  $runapp mem $bopt -t $ncpu $trasmidx $lreads | samtools view --threads $ncpu -Sb -o $outbam
fi

echo "#DONE " `date`
