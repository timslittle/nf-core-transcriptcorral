#!/bin/bash
# run_gnodes_sum.sh
#PBS -N gnodes_sum
#PBS -l vmem=64gb,nodes=1:ppn=8,walltime=1:00:00
#PBS -V

# env asmid=dropse20chrs covtab=dropse20chrs_SRR11813283_1_bwa.cdschr7b.covtab  \
#  anntab=dropse20chrs.fa.anntab  asmdata=dropse20chrs.metad datad=`pwd` qsub -q debug run_gnodes_sum.sh

ncpu=8
evigene=$HOME/bio/evigene
runapp=$evigene/scripts/genoasm/gnodes_covsum.pl

if [ "X" = "X$datad" ]; then echo "ERR:datad=what?"; exit -1; fi
if [ "X" = "X$covtab" ]; then echo "ERR:covtab=what?"; exit -1; fi
if [ "X" = "X$asmid" ]; then echo "ERR:asmid=what?"; exit -1; fi
if [ "X" = "X$opts" ]; then  opts=""; fi
if [ "X" = "X$vers" ]; then  vers=7b; fi

if [ "X" = "X$title" ]; then title=$asmid"_cov"$vers; fi
if [ "X" = "X$anntab" ]; then anntab=$asmid.anntab; fi
if [ "X" = "X$asmdata" ]; then asmdata=$asmid.data; fi

# damn iucarbo perl version change
module switch perl/5.30.1
if [ ! -x $runapp ]; then echo "ERR: miss $runapp"; exit -1; fi

cd $datad
echo "# START `date`"

echo $runapp $opts -asmid $asmid -anntab $anntab -sumdata $asmdata  -title $title  $covtab
$runapp $opts -asmid $asmid -anntab $anntab -sumdata $asmdata  -title $title  $covtab

echo "# DONE `date`"

# dropse20chrs.metad
#-------------
# pt=dropse20chrs # == asmid
# flowcyto=161-180 Mb
# glncformula=174 Mb
# asmtotal=163 Mb
# asmname=Dropse20uc
# # data dependent Ku vals, from various results..
# kulo=94  # dropse20t1cds_SRR11813283_1_bwa new busco.covtab med Covr
# kuhi=95  # dropse20chrs  CDSbusco rdacovm ave; old khi=95
