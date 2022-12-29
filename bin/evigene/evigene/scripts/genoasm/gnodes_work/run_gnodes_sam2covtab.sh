#!/bin/bash
# run_gnodes_sam2covtab.sh
#PBS -N gnodes_sam2covtab
#PBS -l vmem=124gb,nodes=1:ppn=8,walltime=1:00:00
#PBS -V

ncpu=8
##upd action: -merge_readids te.readids cds.readids -out cdste.readids
# env merge_readids="dropse20t1cds_SRR11813283_1_bwa.readids dropse20chrs-families_SRR11813283_1_bwa.readids"  \
#  cdstab=dropse20cdste_SRR11813283_1_bwa.readids  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh
##.. next step:
# env opts="-crclassf=dropse20cdste.idclass" cdstab=dropse20cdste_SRR11813283_1_bwa.readids  \
#   bam=dropse20chrs_SRR11813283_1_bwa.bam  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh

## savereads opt -crclassf cds_te.idclass table for dmag7fincds_tefam20skm: classes = CDS, TE, UNK
## env opts="-savereadids -crclassf dmag7cdste.idclass -minident=0.45"  bam=dmag7cdse_SRRnnn.bam  ..
# upd to multi-cpu mode, makes ncpu part tables to be merged
# always add minalign/minaln opt for that

evigene=$HOME/bio/evigene
runapp=$evigene/scripts/genoasm/gnodes_sam2covtab.pl

if [ "X" = "X$ncpu" ]; then ncpu=8; fi
if [ "X" = "X$datad" ]; then echo "ERR:datad=what?"; exit -1; fi
if [ "X" = "X$opts" ]; then  opts=""; fi
if [ "X" != "X$merge_readids" ]; then  
   opts="$opts -merge readids";
   bam=nobam
fi
if [ "X" = "X$bam" ]; then echo "ERR:bam=what?"; exit -1; fi
if [ "X" = "X$vers" ]; then  vers=7b; fi
if [ "X" = "X$minaln" ]; then  
  echo "WARN: not set minaln=$minaln for $ncpu cpu"; 
else
  opts="$opts -minalign=$minaln"; 
fi
if [ "X" != "X$minid" ]; then  opts="$opts -minident=$minid";  fi
if [ "X" != "X$mindup" ]; then  opts="$opts -mindupid=$mindup";  fi
if [ "X" != "X$skipread" ]; then  opts="$opts -skipread=$skipread";  fi
if [ "X" != "X$chrtab" ]; then  opts="$opts -chrtab $chrtab"; fi
if [ "X" != "X$cdstab" ]; then  opts="$opts -cdstab $cdstab"; fi

module load samtools

ibam=$bam
name=`basename $ibam .bam`
outtab=$name.cdschr$vers.covtab
outchrtab=$name.cdschr$vers.chrtab

cd $datad
echo "# START `date`"
if [ ! -x $runapp ]; then echo "ERR: miss $runapp"; exit -1; fi

if [ "X" != "X$merge_readids" ]; then  
  # save to cdstab; opts="-merge readids -cdstab $cdstab"
  echo $runapp -merge readids -savereadids=$cdstab $merge_readids
  $runapp -merge readids -savereadids=$cdstab $merge_readids
  exit
fi

if [ $ncpu -gt 1 ]; then
  i=0;
  while [ $i -lt $ncpu ]; do {
    echo $runapp -icpu $i -ncpu $ncpu $opts -out $outtab -bam $bam 
    $runapp -icpu $i -ncpu $ncpu $opts -out $outtab -bam $bam  &
    i=$(($i + 1));
  } done
  
  wait

  echo $runapp -merge -out $outtab  -bam $bam 
  $runapp -merge -out $outtab  -bam $bam 
  echo /bin/rm $outtab.pt* $outchrtab.pt*

else

  echo $runapp $opts -out $outtab -bam $bam 
  $runapp $opts -out $outtab -bam $bam 

fi

wc -l $outtab; head -5 $outtab

echo "# DONE `date`"
