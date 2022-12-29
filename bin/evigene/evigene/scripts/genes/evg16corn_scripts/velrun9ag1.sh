#!/bin/bash
## env name=itick inpe=sraf/allpe.fa2 datad=`pwd` qsub -q normal velrun2.sh
#PBS -N velrun
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=39:55:00
#... -l vmem=96gb,nodes=1:ppn=16,walltime=19:55:00
#PBS -V

## 9ag = 9af with kmer shuffle, high only to k43? 
## 9tb: redo 9t, dapx adult9 set best asm, but modify kmers, fail below k47, add above k73
##  .. mod INSIZE, old/bad 150, new 270, split diff at 200? .. same as v9b
## 9t: redo successful 8t, remove/reduce velg covcut (=1 or none)
## .. also update velvet127s to velvet1210 

INSIZE=250
dv=9ag
if [ "X" = "X$ncpu" ]; then ncpu=8; fi
if [ "X" = "X$datad" ]; then echo "err missing datad=path/to/data"; exit -1; fi
if [ "X" = "X$name" ]; then name=ztick; fi
rund=$datad/velv$name$dv
if [ "X" = "X$inpe" ]; then echo "ERR miss inpe=sraf/allpe*.fa2"; exit -1; fi
inpe=`echo $inpe | sed "s,sraf,$datad/sraf,g;"`

noasesfork=1
##......................

export OMP_NUM_THREADS=$ncpu

#k>61# 
velbin2=$HOME/bio/velvet1210/bin2
#k<=61# 
velbin1=$HOME/bio/velvet1210/bin
velbin=$velbin2; 

kset1="93 83 73 63 57 47"
kset="87 77 67 59 53 43 35"

## for norm, ? -cov_cutoff 3 or -cov_cutoff 4, -min_pair_count 3, lower edgeFrac 0.05 ?
iopts="-ins_length $INSIZE -ins_length_sd 50"
vopts="$iopts -max_gap_count 5 " 
oopts="-scaffolding yes -min_pair_count 3 -edgeFractionCutoff 0.03 -min_trans_lgth 200 $iopts"

echo "START `date` " 
cd $datad
mkdir $rund
cd $rund/

#..... run loop for vel kmer steps
shopt -s nullglob

kseqdir=vel${dv}_seq
if [ ! -f $kseqdir/Sequences ]; then
 $velbin/velveth $kseqdir 27 -noHash -shortPaired -fmtAuto -separate $inpe
fi

for k in $kset;  do { 
  ksubdir=vel${dv}_$k
  echo "#.. start velrun $ksubdir : `date`"
  velbin=$velbin2; if [ $k -le 61 ]; then velbin=$velbin1; fi
  mkdir $ksubdir
  ln -s ../$kseqdir/Sequences $ksubdir/
  $velbin/velveth $ksubdir $k  -reuse_Sequences
  $velbin/velvetg $ksubdir $vopts -read_trkg yes 

  # $velbin/oases   $ksubdir $oopts 
  # /bin/rm $ksubdir/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}
  echo "#.. end velvetg $ksubdir : `date`"
} done
## end loop 1

## loop 2 oases, 1cpu; can we fork some at same time? enough mem?
export OMP_NUM_THREADS=2
export OMP_THREAD_LIMIT=2
i=0; 
for k in $kset;  do {
  ksubdir=vel${dv}_$k
  if [ -f $ksubdir/transcripts.fa ]; then continue; fi
  if [ -f $ksubdir/contigs.fa ]; then
    ## DONT fork here save mem ..
    velbin=$velbin2; if [ $k -le 61 ]; then velbin=$velbin1; fi
    $velbin/oases   $ksubdir $oopts 
    i=$(( $i + 1 ))
    #D# if [ $i -ge $noasesfork ]; then wait; i=0; fi
  fi
} done

wait

echo "DONE `date` " 

