#! /bin/bash
## run_gdebug_samsplit2cov8g.sh
## --- gnodes_setup.sh ---    
#SBATCH --job-name="gnodes_pipe"
#SBATCH --output="gnodes_pipe.%j.log"
#SBATCH --partition=debug
#SBATCH -t 1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=240G
#SBATCH --export=ALL

module load  samtools 
## --- end gnodes_setup.sh ---    

STEP=2
sam2covapp=gnodes_sam2covtab8g.pl
export DEBUG_RDFILT=1;

# ncpu = num split runs for icpu=0..ncpu-1, nparts=total out
# samcpu is used for samtools -threads=$samcpu by each icpu task; 8 * 3 = 24 

nparts=24; 
#STEP1: ncpu=8; samcpu=3;
#STEP2: ncpu=12
#STEP3: uses mem for readids, cut down ncpu
ncpu=8
#STEP4: ncpu=12 ?? uses less mem than STEP3

datad=/N/project/eugenes/chrs/daphplx/gasm20set/gnodes20f/pig20gnodes
partdir=cov9test2s

# cov9test2s/
# chrpig11c_SRR4341337_b2_bwa.pt00.bam  pig18evg4wf_t1cds_SRR4341337_b2_bwa.pt00.bam
# chrpig11c_SRR4341337_b2_bwa.pt01.bam  pig18evg4wf_t1cds_SRR4341337_b2_bwa.pt01.bam
# ..
# chrpig11c_SRR4341337_b2_bwa.pt23.bam  pig18evg4wf_t1cds_SRR4341337_b2_bwa.pt23.bam

cdsname=pig18evg4wf_t1cds_SRR4341337_b2_bwa
chrname=chrpig11c_SRR4341337_b2_bwa

# also inbam=pig18evg4wf_t1cds_SRR4341337_b2_bwa.bam
if [ X = X$inbam ]; then inbam=chrpig11c_SRR4341337_b2_bwa.bam; fi
if [ X = X$nparts ]; then nparts=24; fi
if [ X = X$ncpu ]; then ncpu=8; fi
if [ X = X$samcpu ]; then samcpu=3; fi

export EVIGENES=$HOME/bio/evigene/scripts
export PATH=$PATH:$EVIGENES
export NCPU=$ncpu

cd $datad

#-------------------------------
# STEP1: split chr + cds bam to nparts

if [ $STEP = 1 ]; then
  echo START samsplit `date`

  echo $EVIGENES/genoasm/gnodes_samaddrdid.pl -split -icpu=0 -ncpu=$ncpu -samcpu=$samcpu -npart=$nparts -bam $inbam 

  i=0; while [ $i -lt $ncpu ]; do {

    $EVIGENES/genoasm/gnodes_samaddrdid.pl -split -icpu=$i -ncpu=$ncpu -samcpu=$samcpu -npart=$nparts \
     -debug -bam $inbam  &

    sleep 5;
    i=$(($i + 1));
  } done

  wait

  echo DONE samsplit `date`
fi

#-------------------------------
# STEP2: gnodes_sam2covtab.pl cds.ptI.bam to covtab, readids

# orig:  $EVIGENES/genoasm/gnodes_sam2covtab.pl -icpu $i -ncpu $NCPU -nodebug \
#   -crclassf=pig18evg4wf_t1cds.idclass -genetab  \
#   -savereadids=pig18evg4wf_t1cds_SRR4341337_b2_bwa.readids \
#   -out pig18evg4wf_t1cds_SRR4341337_b2_bwa.cc8a.covtab \
#   -bam pig18evg4wf_t1cds_SRR4341337_b2_bwa.bam 

if [ $STEP = 2 ]; then
  echo START cds_sam2covtab $cdsname  `date`
  i=0; icpu=0; 
  while [ $i -lt $nparts ]; do {

    # FIXME: pt$i => pt0$i for i<10 !!!
    ii=$i; if [ $i -lt 10 ]; then ii="0$i"; fi
    
    ibam=$partdir/$cdsname.pt$ii.bam
    irid=$partdir/$cdsname.pt$ii.readids
    icov=$partdir/$cdsname.pt$ii.covtab
    if [ -s $icov ]; then i=$(($i + 1)); continue; fi
 
    echo gnodes_sam2covtab8g.pl -savereadids=$irid -out $icov -bam $ibam   
    
    $EVIGENES/genoasm/gnodes_sam2covtab8g.pl -ncpu 1 -debug \
     -crclassf=pig18evg4wf_t1cds.idclass -genetab  \
     -savereadids=$irid -out $icov -bam $ibam  &

    sleep 5;
    i=$(($i + 1)); icpu=$(($icpu + 1)); 
    if [ $icpu -ge $ncpu ]; then wait; icpu=0; fi

  } done
  
  wait
  echo DONE cds_sam2covtab  `date`
fi

#-------------------------------
# STEP3: gnodes_samaddrdid.pl to chr.ptI.bam with cds.ptI.readids

if [ $STEP = 3 ]; then
  echo START samaddrdid $chrname `date`
  i=0; icpu=0; 
  while [ $i -lt $nparts ]; do {

    ii=$i; if [ $i -lt 10 ]; then ii="0$i"; fi
    ibam=$partdir/$chrname.pt$ii.bam
    irid=$partdir/$cdsname.pt$ii.readids
    obam=$partdir/$chrname.pt$ii.rdid.bam
    if [ -s $obam ]; then i=$(($i + 1)); continue; fi
    
    $EVIGENES/genoasm/gnodes_samaddrdid.pl -debug -ncpu 1 -ridprefix=SRR4341337. \
       -readidtab=$irid -outbam $obam -bam $ibam  &

    i=$(($i + 1)); icpu=$(($icpu + 1)); 
    if [ $icpu -ge $ncpu ]; then wait; icpu=0; fi

  } done
  
  wait
  echo DONE samaddrdid  `date`
fi

#-------------------------------
# STEP4: gnodes_sam2covtab8g.pl -genetab -gene chr.ptI.rdid.bam  

# orig:    $EVIGENES/genoasm/gnodes_sam2covtab.pl -icpu $i -ncpu $NCPU \
#   -crclassf=pig18evg4wf_t1cds.idclass -genetab  \
#   -readidtab=pig18evg4wf_t1cds_SRR4341337_b2_bwa.readids \
#   -out chrpig11c_SRR4341337_b2_bwa.cc8a.covtab \
#   -bam chrpig11c_SRR4341337_b2_bwa.bam  
#.. and merge ...
#   $EVIGENES/genoasm/gnodes_sam2covtab.pl  -merge  -crclassf=pig18evg4wf_t1cds.idclass -genetab \
#   -readidtab=pig18evg4wf_t1cds_SRR4341337_b2_bwa.readids \
#   -out chrpig11c_SRR4341337_b2_bwa.cc8a.covtab \
#   -bam chrpig11c_SRR4341337_b2_bwa.bam 

if [ $STEP = 4 ]; then
  echo START sam2covtab8 $chrname  `date`
  i=0; icpu=0; 
  while [ $i -lt $nparts ]; do {

    ii=$i; if [ $i -lt 10 ]; then ii="0$i"; fi
    ibam=$partdir/$chrname.pt$ii.rdid.bam
    icov=$partdir/$chrname.pt$ii.covtab
    if [ -s $icov ]; then i=$(($i + 1)); continue; fi

    echo gnodes_sam2covtab8g.pl -geneidsInBam -genetab -out $icov -bam $ibam   

    $EVIGENES/genoasm/gnodes_sam2covtab8g.pl -geneidsInBam -genetab -debug -ncpu 1  \
      -crclassf=pig18evg4wf_t1cds.idclass  \
      -out $icov -bam $ibam &

    i=$(($i + 1)); icpu=$(($icpu + 1)); 
    if [ $icpu -ge $ncpu ]; then wait; icpu=0; fi

  } done
  
  wait
  echo DONE sam2covtab8  `date`
fi

if [ $STEP = 5 ]; then
  echo START sam2covtab8 merge  `date`

  # NOTE: this finds parts as:  ls $ofile.pt*
  ## NOT same as now: chrname.pt*.covtab .. 
  ## can add parts on cmdline, with -merge=covtab|chrtab|genetab

  for suf in covtab chrtab genetab; do {
    otab=$chrname.merge.$suf
    parts=$partdir/$chrname.pt*.$suf
    
   echo $EVIGENES/genoasm/gnodes_sam2covtab8g.pl -merge=$suf -genetab -debug   \
      -crclassf=pig18evg4wf_t1cds.idclass  \
      -out $otab $parts  
      
  } done
  
  echo DONE sam2covtab8 merge `date`
fi

