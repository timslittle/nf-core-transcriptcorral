#! /bin/bash
### env datad=path/to/data trset=myspecies_all.tr.gz qsub -q normal runtr2cds.sh
#PBS -N tr2aacds
#PBS -A PutAccountIdHere
#PBS -l nodes=1:ppn=16,walltime=39:55:00
#PBS -V

# run_s07k_tr2aacds.sh = tr2aacds2c.pl many small updates;  -aconsensus
# run_s07p_tr2cdspub.sh = more updates, exon-chain filter, includes trclass2pubset step
## runs in common test folder tr2aacds_test1908f/ for human, weed, pig, ..

if [ "X" = "X$ncpu" ]; then ncpu=8; fi
if [ "X" = "X$maxmem" ]; then maxmem=64000; fi
if [ "X" = "X$datad" ]; then echo missing datad=/path/to/data; exit -1; fi
if [ "X" = "X$trset" ]; then echo missing trset='trev19human.tr'; exit -1; fi
if [ "X" = "X$idpre" ]; then idpre=EV4m; fi
# if [ "X" = "X$name" ]; then name=`basename $trset .tr | sed 's/\.fa.*$//'`; fi
trname=`basename $trset .tr | sed 's/\.fa.*$//'`;

## $ORF_FULLvPART = $ENV{ORF_FULLvPART} || 0.85; # 0.85 gets many 5'part extensions, common gene biology
export ORF_FULLvPART=0.50
##  cds extend2utr for cdhit1 to keep 5'alts;  may be better but failed: CDSXUTR=2400,1
export CDSXUTR=900,60

##  change tr2aacds2d/asmdupfilt classify alt drop due to utrpoor/bad, major loss of ref alts
export pCDSOK=20 pCDSBAD=20 
##  -runsteps=noaadup turn off althi drops by aadup, major loss of ref alts
addopt=" -aconsensus -runsteps=noaadup"

# update evigene19sep03 to evigene19oct02
evigenes=$HOME/bio/evigene19sep03/scripts
export PATH=$HOME/bio/apps/cdhit/bin:$PATH
export fastanrdb=$HOME/bio/apps/exonerate/bin/fastanrdb
export PATH=$HOME/bio/apps/ncbi/bin:$PATH

evapp=$evigenes/prot/tr2aacds2d.pl
#use also: makeblastscore3.pl genes/blasttrset2exons2.pl rnaseq/asmrna_altreclass3c.pl

traopts="-tidy -log -debug"
if [ "X" != "X$addopt" ]; then traopts="$traopts $addopt"; fi
if [ "X" != "X$opt" ]; then traopts="$traopts $opt"; fi

cd $datad/

echo "#START `date` " 

echo '#STEP1: tr2aacds reduction of trset' 
if [ -s $trname.trclass ]; then
echo DONE $evapp -NCPU $ncpu -MAXMEM $maxmem $traopts -cdna $trset
else
echo $evapp -NCPU $ncpu -MAXMEM $maxmem $traopts -cdna $trset
$evapp -NCPU $ncpu -MAXMEM $maxmem $traopts -cdna $trset
fi

echo '#STEP2: cdsqual,  makeblastscore, blasttrset2exons, trclass2pubset, asmrna_altreclass3c' 

# NOTE: this trnrname keeps tmpfiles/ path, dont add below
trnrname=`/bin/ls tmpfiles/$trname*-self98.blastn | sed 's/-self98.blastn//;'`

echo env outcds=1 $evigenes/prot/cdsqual.sh $trnrname.cds
env outcds=1 $evigenes/prot/cdsqual.sh $trnrname.cds

echo $evigenes/makeblastscore3.pl -pIDENTMIN 99.999 -pmin 0.01 -CDSSPAN -showspan=2 -tall \
  -sizes $trnrname.cds.qual $trnrname-self98.blastn TO $trnrname-self100.btall

if [ -s $trnrname.cds.qual -a -s $trnrname-self98.blastn ]; then
  $evigenes/makeblastscore3.pl -pIDENTMIN 99.999 -pmin 0.01 -CDSSPAN -showspan=2 -tall \
  -sizes $trnrname.cds.qual $trnrname-self98.blastn > $trnrname-self100.btall
else
  echo "Fail makeblastscore3.pl input files"; exit -1;
fi

#NOTE: need trclass2pubset.pl in 2 steps here, 1st is env norealt=1 to skip altreclass, 2nd after blasttrset2exons
echo env norealt=1 $evigenes/genes/trclass2pubset.pl -onlypub -noaltdrops  -idpre $idpre -log -debug  -class $trname.trclass
env norealt=1 $evigenes/genes/trclass2pubset.pl -onlypub -noaltdrops  -idpre $idpre -log -debug  -class $trname.trclass

inbtall=$trnrname-self100.btall; 
pubids=publicset/$trname.pubids;
echo "sort -k7,7nr -k2,2 -k6,6nr -k1,1 $inbtall | env pubids=$pubids debug=1 $evigenes/genes/blasttrset2exons2.pl TO $trnrname.exontab"

if [ -s $inbtall -a -s $pubids ]; then
  sort -k7,7nr -k2,2 -k6,6nr -k1,1 $inbtall | env pubids=$pubids debug=1 $evigenes/genes/blasttrset2exons2.pl > $trnrname.exontab
else
  echo "Fail blasttrset2exons2.pl input files"; exit -1;
fi
  
# final step, integrate into trclass2pubset; see trclass2pubset:altreclass_block()
echo $evigenes/rnaseq/asmrna_altreclass3c.pl -debug -nodrops -noclasscut=60 -altrenum \
 -trclass $trname.trclass -pubids publicset/$trname.pubids  \
 -exontable $trnrname.exontab  -out publicset/$trname.pubids.realt

if [ -s $trnrname.exontab ]; then
  $evigenes/rnaseq/asmrna_altreclass3c.pl -debug -nodrops -noclasscut=60 -altrenum \
  -trclass $trname.trclass -pubids publicset/$trname.pubids  \
  -exontable $trnrname.exontab  -out publicset/$trname.pubids.realt \
  > publicset/$trname.altreclass3c.log 2>&1

  if [ -s  publicset/$trname.pubids.realt ]; then
    mv publicset/$trname.pubids publicset/$trname.pubids.old
    mv publicset/$trname.pubids.realt publicset/$trname.pubids; touch publicset/$trname.pubids.altreclass3c.done
  fi
else
  echo "Fail asmrna_altreclass3c.pl input files"; exit -1;
fi
  
echo "#DONE : `date`"
