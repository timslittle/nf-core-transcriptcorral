#! /bin/bash
# run_gnokmers1rdset.sh 
#SBATCH --job-name="gnodek_pipe"
#SBATCH --output="gnodek_pipe.%j.log"
#.. SBATCH --partition=debug
#SBATCH --partition=general
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=240G
#SBATCH --export=ALL

## change from many srrid to one readset (maybe many rd123_[12].fastq)..
## jellyfish make kmer histo table, then GScope.R, findGSE.R; ttag=gnok_at1u
## also maybe MGSE, but need also build, sort readmap.bam, bedtools/depth table extras...

## NOTE required param: readlen    get from reads? YES, see run_gnodesqwik_at1ksw.sh, run_rdecontam.sh
## readname.count == 1 line of readlen \t nreads \t nbases for readname_[12].fastq
# if [ "X" = "X$readlen" ]; then echo "ERR:readlen=what?"; exit -1; fi
 
export NCPU=32 MEMT=240G MEMC=3G
## last go used 240G x 24cpu w/ jellyfish, need that much?

## this kmers opt needs check: findGSE used diff sets, takes ave of successful kmers
## should skip any jellyfish kmer fail, if it does; skip or not Gscope,fiGSE that fail when jf ok
if [ X = "X$kmers" ]; then kmers="19 21 29"; fi

# if [ X = "X$intab" ]; then intab=at1k_sraruninfo1a.tsv; fi
# if [ X = "X$srrid" -a -s $intab ]; then srrid=`cut -f1 $intab | egrep '^[DES]RR'`; fi
# if [ X = "X$srrid" ]; then echo srrid=what; exit -1; fi
if [ X = "X$reads" ]; then echo "reads=what?"; exit -1; fi
if [ X = "X$readlen" ]; then echo "readlen=what?"; exit -1; fi

# if [ X = "X$maxrun" ]; then maxrun=32; fi
if [ X = "X$RM_TMP" ]; then RM_TMP=0; fi
if [ X = "X$ttag" ]; then ttag="gnokmer"; fi
if [ X = "X$datad" ]; then echo datad=what; exit -1; fi

# chrfa=arath18tair_chr.fa
# if [ X = "X$chrfa" ]; then echo chrfa=what; exit -1; fi
# asmid=`basename $chrfa .fa`
asmid=gsek
if [ X = "X$title" ]; then title=${asmid}_$ttag; fi

# kbin >> need in local dir: findGSE.R genomescope.R
module add r
kbin=$HOME/bio/gnomutil/khisto
fubin=$HOME/bio/gnomutil/bin/
jfbin=$HOME/bio/gnomutil/bin/
export PATH=$PATH:$jfbin

export EVIGENES=$HOME/bio/apps/evigene/scripts
export PATH=$PATH:$EVIGENES
## --- end gnodes_setup.sh ---    

cd $datad
echo START gnodes_kmers $title  `date`

mkdir  ${title}_tmp/
mkdir  ${title}_out/
# replace  ${title}_data/ with option
# if [ X = "X$readsf" ]; then  readsf=${title}_data; fi
# if [ ! -d $readsf ]; then  mkdir  $readsf; fi

#DROP ====== loop here for SRRid, download reads =============

irun=0;
# for sid in $srrid; do {

# if [ $irun -ge $maxrun ]; then echo "DONE gnodes_kmers $irun/$maxrun runs.. bye."; break; fi
# echo START gnodes_kmers $sid  `date`
# sinfo=`grep $sid $intab`
# echo "INFO" $sinfo

read1=`echo $reads | sed 's/ .*//;'`
# read1=$readsf/${sid}_1.fastq
# read2=$readsf/${sid}_2.fastq
readlen1=$readlen
readcnt=`echo $read1 | sed 's/_1.fastq/.count/;'`
if [ -s $readcnt ]; then readlen1=`cut -f1 $readcnt`; fi

# does jfish read gz? NO
read1noz=`echo $read1 | sed 's/.gz//g;'`
readnoz=`echo $reads | sed 's/.gz//g;'`
if [ $read1 != $read1noz ]; then 
 gunzip $reads ;
 readz=$reads; reads=$readnoz;
fi
# if [ ! -s $read1 ]; then  echo "FAIL sra data: $read1"; exit -1; fi

# reads=$readsf/${sid}_[12].fastq
readna=`basename $read1 | sed 's/\.gz//; s/\.[a-z0-9]*$//;' `
readna=`echo $readna | sed 's/_[12].*/ab/;'`
tmpf=""; outf=""

## insert here fastuniq.sh ; calc stats for both?? fu_SRR bad, _ cut in sum.txt

## kmer loop here
for kmer in $kmers; do {

  kname=$readna.jfk$kmer
  readdb=$kname.jfish

  ## Masurca uses:
  # user opt: JF_SIZE = 7_000_000_000 = 7G
  # jellyfish count -m 31 -t 32 -C -s $JF_SIZE -o k_u_hash_0 pe.cor.fa
  # export ESTIMATED_GENOME_SIZE=`jellyfish histo -t 32 -h 1 k_u_hash_0 | tail -n 1 |awk '{print $2}'`
  ##  need lower hsize,csize for higher kmer=29 .. change here, see Masurca usage
  # hsize=25G; csize=6; if [ $kmer -gt 23 ]; then hsize=15G; csize=4;  fi
  
  hsize=15G; csize=4;

if [ -s $kname.hist ]; then
  echo reusing $kname.hist;
else
  echo START jellyfish $kname  `date`
  echo jellyfish count -m $kmer  -o $readdb -t $NCPU -s $hsize -c $csize --canonical  $reads
  jellyfish count -m $kmer  -o $readdb -t $NCPU -s $hsize -c $csize --canonical  $reads

  echo DONE jellyfish $kname  `date`
  if [ ! -s $readdb ]; then
    echo "# FAIL: jellyfish count -m $kmer $reads";   continue
  fi  

  jellyfish histo $readdb > $kname.hist
  # basic genome size est
  jellyfish histo  -t $NCPU -h 1  $readdb > $kname.jf1gse.txt
  
  #? erase $readdb : YES
  rm $readdb
  outf="$outf $kname.hist $kname.jf1gse.txt"
fi
  
  ## manage outdirs for fgse, gscop, dont need all their outs, just summary table 
  # R; install.packages("pracma"); install.packages("fGarch"); q("no");
  #Rscript -e source("findGSE.R"); findGSE(histo="ERR2178784ab.jfk21.hist", sizek=21, outdir="ERR2178784ab.jfk21_fgse" );
  # GenomeScope analyzing ERR2178784ab.jfk21.hist k=21 readlen=251 outdir=ERR2178784ab.jfk21_gscop
  # Error in .External2(C_X11, paste0("png::", filename), g$width, g$height,  :    unable to start device PNG
  ## replaced png(xxx) w/ pdf(xxx)

  odir=$kname"_fgse"
  rsc='source("'$kbin'/findGSE.R"); findGSE(histo="'$kname.hist'", sizek='$kmer', outdir="'$odir'" );'
  echo Rscript -e "$rsc";
  Rscript -e "$rsc";
  # odir/v1.94.est.ERR2178784ab.jfk21.hist.genome.size.estimated.k21to21.fitted.txt
  cp -p $odir/*fitted.txt  ${title}_out/$odir.sum.txt
  tmpf="$tmpf $odir"
  
  odir=$kname"_gscop"
  echo Rscript $kbin/genomescope.R  $kname.hist $kmer $readlen1 $odir
  Rscript $kbin/genomescope.R  $kname.hist $kmer $readlen1 $odir
  cp -p $odir/summary.txt  ${title}_out/$odir.sum.txt
  tmpf="$tmpf $odir"
  
} done

echo "#TEMPS   : $tmpf ${title}_tmp/"
echo "#OUTPUTS : $outf ${title}_out/"
if [ X != "X$tmpf" ]; then mv $tmpf  ${title}_tmp/; fi
if [ X != "X$outf" ]; then mv $outf  ${title}_out/; fi
if [ $RM_TMP = 1 ]; then
  cd ${title}_tmp/; rm $tmpf; 
  cd ../
fi

#o if [ -s $read1 ]; then 
#o   gzip --fast $read1 &
#o   gzip --fast $read2 &
#o   wait
#o fi

# echo DONE gnodes_kmers $sid  `date`

#off } done
#====== NO end loop for SRRid, download reads =============

echo DONE gnodes_kmers $title  `date`

