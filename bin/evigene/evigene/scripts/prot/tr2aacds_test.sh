#! /bin/bash
## RUN-LOCAL: env trset=arath_TAIR10_20101214up.cdna.gz datad=`pwd` ./tr2aacds_qsub.sh
## RUN-CLUSTER: env trset=arath_TAIR10_20101214up.cdna.gz datad=`pwd` qsub -q normal tr2aacds_qsub.sh
#PBS -N tr2cds
#PBS -l nodes=1:ppn=32,walltime=18:55:00
#PBS -o tr2cds.$$.out
#PBS -e tr2cds.$$.err
#PBS -V

function usage {
cat << 'EOT'

# TEST RUN THIS outside of PBS/qsub cluster to see if it works:
  env trset=arath_TAIR10_20101214up.cdna.gz datad=`pwd` ./tr2aacds_test.sh > & log.tr2ac1

# TEST CASE: arabidopsis TAIR10 transcripts (headers regularized, though probably not needed)
  curl -sR -o arath_TAIR10_20101214up.cdna \
    ftp://ftp.arabidopsis.org/home/tair/Sequences/blast_datasets/TAIR10_blastsets/TAIR10_cdna_20101214_updated
# get also TAIR10_pep and TAIR10_cds to compare w/ evigene output
# this takes 20 min to run on 2 cores of intel game box 

EOT
exit -1
}

ncpu=2
maxmem=1000
bapps=/camp/apps/eb/software
evigene=./evigene/scripts

#t2ac: app=cd-hit-est, path= echo MISSING_cd-hit-est
export PATH=$PATH:/camp/apps/eb/software/CD-HIT/4.8.1-GCC-10.2.0/bin
#t2ac: app=fastanrdb, path= echo MISSING_fastanrdb
export fastanrdb=/camp/apps/eb/software/Exonerate/2.4.0-foss-2016b/bin/fastanrdb
#t2ac: app=blastn, path= echo MISSING_blastn
export PATH=$PATH:/camp/apps/eb/software/BLAST+/2.11.0-gompi-2020b/bin

export makeblastdb=/camp/apps/eb/software/BLAST+/2.11.0-gompi-2020b/bin

if [ "X" = "X$trset" ]; then
  echo "Missing env trset=xxxx.tr"; usage; 
fi
if [ "X" = "X$datad" ]; then
  echo "Missing env datad=/path/to/data"; usage;
fi

cd $datad/

echo $evigene/prot/tr2aacds4.pl -tidy -NCPU $ncpu -MAXMEM $maxmem -log -cdna $trset
$evigene/prot/tr2aacds4.pl -tidy -NCPU $ncpu -MAXMEM $maxmem -log -cdna $trset

