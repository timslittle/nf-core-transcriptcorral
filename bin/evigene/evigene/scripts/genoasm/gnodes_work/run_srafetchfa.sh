#! /bin/bash
## --- gnodes_setup.sh for quartz.iu ---    
#SBATCH --job-name="gnodes_pipe"
#SBATCH --output="gnodes_pipe.%j.log"
#SBATCH --partition=general
#SBATCH -t 8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=40G
#SBATCH --export=ALL

if [ X = "X$sid" ]; then echo "ERR sid=what SRR id"; exit -1; fi
if [ X = "X$fasta" ]; then fasta=0; fi;
if [ X = "X$datad" ]; then echo "ERR datad=what path"; exit -1; fi
# subdir=readsf or ./

nbin=$HOME/bio/ncbi/bin; export PATH=$nbin:$PATH;
sambin=$HOME/bio/apps/bin; export PATH=$sambin:$PATH;
sradumpr=$HOME/bio/sratoolkit/sratoolkit300/bin/fasterq-dump

export NCPU=16
export EVIGENES=$HOME/bio/apps/evigene/scripts
export PATH=$PATH:$EVIGENES

cd $datad
mkdir readsf

if [ $fasta = 1 ]; then
  if [ ! -s readsf/${sid}.fasta.gz ]; then
  echo $sradumpr --threads $NCPU  --fasta --outdir readsf $sid
  $sradumpr --threads $NCPU  --fasta --outdir readsf/ $sid >& readsf/sradump.$sid.log 
  fi
  if [ -s readsf/${sid}.fasta ]; then
    gzip --fast readsf/${sid}.fasta
  fi
elif [ ! -s readsf/${sid}_1.fastq.gz ]; then
  echo $sradumpr --threads $NCPU  --split-3 --outdir readsf $sid
  if [ ! -s readsf/${sid}_1.fastq ]; then
    $sradumpr --threads $NCPU  --split-3 --outdir readsf/ $sid >& readsf/sradump.$sid.log 
  fi
  if [ -s readsf/${sid}_1.fastq ]; then
    gzip --fast readsf/${sid}_1.fastq &
    gzip --fast readsf/${sid}_2.fastq &
    wait
  else
    echo "ERR: missing readsf/${sid}_[12].fastq.gz"; exit -1;
  fi
fi

