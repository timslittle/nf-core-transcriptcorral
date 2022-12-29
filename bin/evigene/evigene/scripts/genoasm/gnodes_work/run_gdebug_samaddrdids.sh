#! /bin/bash
## run_gdebug_samaddrdids.sh
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

STEP=1;   # samaddrdid -readids cdsname.readids -bam $chrname.bam
# STEP=2; # gnodes_sam2covtab8g  -geneidsInBam -genetab -bam $chrname.rdid.bam
# STEP=3; # gnodes_sam2covtab8g -merge -genetab -bam $chrname.rdid.bam

# pig s1,try1: 4cpu, 12/24 parts of cds.readids (2 GB) + chr.pt.bam (32 GB which is 1/2 of 40% of 167 GB total)
# slurm job stat 13 min (after readid): tcpu= 60 min, maxmem= 33 Gb
# try2: 12cpu, 12/24 parts, 2.4 hr cpu, 0.18 hr wall, 99 Gb maxmem
# try3? 12cpu, 24/24 parts?  est 3 * 2 * 33 Gb = 198 Gb?
#             JOBID PARTITION     NAME     USER    STATE       TIME TIME_LIMI  NODES NODELIST(REASON)
#            630630     debug gnodes_p gilbertd  RUNNING      10:16   1:00:00      1 c1
# sstat --format=AveCPU,AvePages,AveRSS,AveVMSize,MaxVMSize,JobID -a -j 630630
#    AveCPU   AvePages     AveRSS  AveVMSize  MaxVMSize        JobID 
#  02:23:48          0  96122696K  98538076K  98538208K 630630.batch 
#..
#.. try4/5s: 12 cpu, all orig readids (6.3 G), part4s chr.bam 33 Gb >> 213 Gb vmem used after readids
#            631150     debug gnodes_p gilbertd  RUNNING      10:38   1:00:00      1 c1
#    AveCPU   AvePages     AveRSS  AveVMSize  MaxVMSize        JobID 
#  02:24:38          0 203769732K 212736312K 212 736 312K 631150.batch 
#  08:25:14          0 203774324K 213065792K 213_065_792K 631150.batch  : 34 min runtime


nparts=24;  ncpu=8; samcpu=3;
if [ $STEP = 1 ]; then ncpu=12; fi
if [ $STEP = 2 ]; then ncpu=24; fi

# ncpu = num split runs for icpu=0..ncpu-1, nparts=total out
# samcpu is used for samtools -threads=$samcpu by each icpu task; 8 * 3 = 24 

sam2covapp=gnodes_sam2covtab8g.pl
export DEBUG_RDFILT=1;

datad=/N/project/eugenes/chrs/daphplx/gasm20set/gnodes20f/pig20gnodes

#24 parts: 
partdir2s=cov9test2s
# test4: cat 12 parts of chr.ptNN.bam and cds.ptNN.readids to one data set
partdir=cov9test4s

# cov9test2s/
# chrpig11c_SRR4341337_b2_bwa.pt00.bam  pig18evg4wf_t1cds_SRR4341337_b2_bwa.pt00.bam
# chrpig11c_SRR4341337_b2_bwa.pt23.bam  pig18evg4wf_t1cds_SRR4341337_b2_bwa.pt23.bam

cdsname=pig18evg4wf_t1cds_SRR4341337_b2_bwa
chrname=chrpig11c_SRR4341337_b2_bwa
ridprefix=SRR4341337.
gidprefix=Susscr4EV  # Susscr4EVm000001t1 cut('t1') also; leave on 'm' or only digits?
# gnodes readReadIds has quirk expects geneids to have [A-Za-z] start, other cols are class cols **
# gidprefix cut affects that; add -hasgeneids bool opt, last col in .readids == geneids

export EVIGENES=$HOME/bio/evigene/scripts
export PATH=$PATH:$EVIGENES
export NCPU=$ncpu

cd $datad

#-------------------------------

# STEP0: dataset:
# cat $partdir2s/$cdsname.pt0[0-9].readids $partdir2s/$cdsname.pt1[01].readids > $partdir/$cdsname.readids
# samtools cat -o $partdir/$chrname.bam $partdir2s/$chrname.pt0[0-9].bam  $partdir2s/$cdsname.pt1[01].bam

inreadids=$partdir/$cdsname.readids
inbam=$partdir/$chrname.bam
ridbam=$partdir/$chrname.rdid.bam
ocovtab=$partdir/$chrname.db9g.covtab

if [ $STEP = 1 ]; then
  echo START samaddrdid $chrname `date`

  i=0; icpu=0; # i == icpu here
  while [ $i -lt $ncpu ]; do {
    
    ## samaddrdid makes outbam partnames for ncpu: $outbampart =~ s/\.bam/.pt$icpu.bam/; 
    # ii=$i; if [ $i -lt 10 ]; then ii="0$i"; fi
    # obami=$partdir/$chrname.rdid.pt$ii.bam

    $EVIGENES/genoasm/gnodes_samaddrdid.pl -icpu $i -ncpu $ncpu -samcpu $samcpu -debug \
       -ridprefix=$ridprefix -gidprefix=$gidprefix \
       -hasgeneids -readidtab=$inreadids -outbam $ridbam -bam $inbam  &

    i=$(($i + 1)); 
    # icpu=$(($icpu + 1)); if [ $icpu -ge $ncpu ]; then wait; icpu=0; fi
  } done
  
  wait

  $EVIGENES/genoasm/gnodes_samaddrdid.pl -merge -bam $ridbam
  # == samtools cat -o $ridbam $partdir/$chrname.rdid.pt*.bam
  
  echo DONE samaddrdid $chrname `date`
fi


if [ $STEP = 2 ]; then

  # this test: 1 rdid.bam w/ readids, iterate over ncpu to read, write part.covtab,genetab,..
  
  echo START sam2covtab8g $chrname  `date`
  i=0; icpu=0; # i == icpu here
  while [ $i -lt $ncpu ]; do {

    echo gnodes_sam2covtab8g.pl -geneidsInBam -genetab -out $ocovtab -bam $ridbam   

    # ii=$i; if [ $i -lt 10 ]; then ii="0$i"; fi
    ## sam2cov names parts: ocovi=$partdir/$chrname.pt$ii.covtab

    $EVIGENES/genoasm/gnodes_sam2covtab8g.pl -geneidsInBam -genetab -icpu $i -ncpu $ncpu  -debug  \
      -crclassf=pig18evg4wf_t1cds.idclass  \
      -out $ocovtab -bam $ridbam &

    i=$(($i + 1)); 
    # icpu=$(($icpu + 1)); if [ $icpu -ge $ncpu ]; then wait; icpu=0; fi

  } done
  
  wait
  
  echo DONE sam2covtab8g  `date`
fi

if [ $STEP = 3 ]; then
  # merge covtab parts, add to step 2?
  echo START sam2covtab8g merge $ocovtab   `date`
  
  echo $EVIGENES/genoasm/gnodes_sam2covtab8g.pl -merge  -genetab -debug  \
      -crclassf=pig18evg4wf_t1cds.idclass  \
      -out $ocovtab  
      
  $EVIGENES/genoasm/gnodes_sam2covtab8g.pl -merge  -genetab -debug  \
      -crclassf=pig18evg4wf_t1cds.idclass  \
      -out $ocovtab  

  echo DONE sam2covtab8g merge $ocovtab   `date`

fi

