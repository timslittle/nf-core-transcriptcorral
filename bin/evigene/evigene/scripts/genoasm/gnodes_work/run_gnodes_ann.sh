#!/bin/bash
# run_gnodes_ann.sh
#PBS -N gnodes_ann
#PBS -l vmem=124gb,nodes=1:ppn=16,walltime=1:00:00
#PBS -V

ncpu=16
evigene=$HOME/bio/evigene
runapp=$evigene/scripts/genoasm/gnodes_annotate.pl

if [ "X" = "X$ncpu" ]; then ncpu=8; fi
if [ "X" = "X$datad" ]; then echo "ERR:datad=what?"; exit -1; fi
if [ "X" = "X$chr" ]; then echo "ERR:chr=what?"; exit -1; fi
if [ "X" = "X$opts" ]; then  opts=""; fi
if [ "X" != "X$cds" ]; then  opts="$opts -cds $cds"; fi
if [ "X" != "X$te" ]; then  opts="$opts -te $te"; fi
if [ "X" != "X$idclass" ]; then  opts="$opts -idclass $idclass"; fi
if [ "X" = "X$vers" ]; then  vers=7b; fi

# fix iucarbo perl version change:  module switch perl/5.30.1
module load blast
if [ ! -x $runapp ]; then echo "ERR: miss $runapp"; exit -1; fi

cd $datad
echo "# START `date`"

echo $runapp $opts -ncpu $ncpu -chr $chr
$runapp $opts -ncpu $ncpu -chr $chr

echo "# DONE `date`"

