#! /bin/bash
## run_gnokmers_atvlr.sh for vlread from run_gnokmers_at1k.sh
#SBATCH --job-name="gnodek_pipe"
#SBATCH --output="gnodek_pipe.%j.log"
#.. SBATCH --partition=debug
#SBATCH --partition=general
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=240G
#SBATCH --export=ALL

## jellyfish make kmer histo table, then GScope.R, findGSE.R; ttag=gnok_at1u
## also maybe MGSE, but need also build, sort readmap.bam, bedtools/depth table extras...

## NOTE required param: readlen    get from reads? YES, see run_gnodesqwik_at1ksw.sh, run_rdecontam.sh
## readname.count == 1 line of readlen \t nreads \t nbases for readname_[12].fastq
# if [ "X" = "X$readlen" ]; then echo "ERR:readlen=what?"; exit -1; fi
 
export NCPU=32 MEMT=120G MEMC=3G
## last go used 240G x 24cpu w/ jellyfish, need that much?

## this kmers opt needs check: findGSE used diff sets, takes ave of successful kmers
## should skip any jellyfish kmer fail, if it does; skip or not Gscope,fiGSE that fail when jf ok
if [ X = "X$kmers" ]; then kmers="19 21 29"; fi

if [ X = "X$intab" ]; then intab=at1k_sraruninfo1a.tsv; fi
if [ X = "X$srrid" -a -s $intab ]; then
  srrid=`cut -f1 $intab | egrep '^[DES]RR'`
fi
if [ X = "X$srrid" ]; then echo srrid=what; exit -1; fi
if [ X = "X$readlen" ]; then echo "readlen=what?"; exit -1; fi

if [ X = "X$maxrun" ]; then maxrun=32; fi
if [ X = "X$dofu" ]; then dofu=0; fi
if [ X = "X$RM_TMP" ]; then RM_TMP=0; fi
if [ X = "X$ttag" ]; then ttag="gnok_at1u"; fi

if [ X = "X$datad" ]; then echo datad=what; exit -1; fi

asmid="asmnone"; # `basename $chrfa .fa`
if [ X = "X$title" ]; then title=${asmid}_$ttag; fi

# need in local dir: findGSE.R genomescope.R
module add r

nbin=$HOME/bio/ncbi/bin; export PATH=$nbin:$PATH;
sambin=$HOME/bio/apps/bin; export PATH=$sambin:$PATH;
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
if [ X = "X$readsf" ]; then  readsf=${title}_data; fi
if [ ! -d $readsf ]; then  mkdir  $readsf; fi

#====== loop here for SRRid, download reads =============

irun=0;
for sid in $srrid; do {

if [ $irun -ge $maxrun ]; then echo "DONE gnodes_kmers $irun/$maxrun runs.. bye."; break; fi

echo START gnodes_kmers $sid  `date`
sinfo=`grep $sid $intab`
echo "INFO" $sinfo

# step0. check if done SRRid (sraid?); need readna,chrbname here before reads exist
readna=${sid}ab
testf=${title}_out/$readna.jfk19_fgse.sum.txt 
if [ -s $testf ]; then echo "# DONE $sid in $testf; continue.. ";  continue; fi

irun=$(($irun + 1)); 
 
read1=$readsf/${sid}.fasta
readlen1=$readlen
readcnt=`echo $read1 | sed 's/.fasta//; s/$/.count/;'`
if [ -s $readcnt ]; then readlen1=`cut -f1 $readcnt`; fi

# does jfish read gz? NO
if [ -s $read1.gz ]; then gunzip $read1.gz; fi
if [ ! -s $read1 ]; then  echo "FAIL sra data: $read1"; exit -1; fi

reads=$read1
# $readsf/${sid}_[12].fastq
readna=`basename $reads | sed 's/\.gz//; s/\.[a-z0-9]*$//;' `
readna=`echo $readna | sed 's/_[12].*/ab/;'`
tmpf=""; outf=""

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

  jellyfish histo  -t $NCPU $readdb > $kname.hist
  # basic genome size est
  jellyfish histo  -t $NCPU -h 1  $readdb > $kname.jf1gse.txt
  
  #? erase $readdb : YES
  rm $readdb
  # tmpf="$tmpf $readdb"
  outf="$outf $kname.hist $kname.jf1gse.txt"
fi
  
  ## manage outdirs for fgse, gscop, dont need all their outs, just summary table 
  # R; install.packages("pracma"); install.packages("fGarch"); q("no");
  #Rscript -e source("findGSE.R"); findGSE(histo="ERR2178784ab.jfk21.hist", sizek=21, outdir="ERR2178784ab.jfk21_fgse" );

  # GenomeScope analyzing ERR2178784ab.jfk21.hist k=21 readlen=251 outdir=ERR2178784ab.jfk21_gscop
  # Error in .External2(C_X11, paste0("png::", filename), g$width, g$height,  :    unable to start device PNG
  ## replaced png(xxx) w/ pdf(xxx)

  odir=$kname"_fgse"
  rsc='source("findGSE.R"); findGSE(histo="'$kname.hist'", sizek='$kmer', outdir="'$odir'" );'
  echo Rscript -e "$rsc";
  Rscript -e "$rsc";
  # odir/v1.94.est.ERR2178784ab.jfk21.hist.genome.size.estimated.k21to21.fitted.txt
  cp -p $odir/*fitted.txt  ${title}_out/$odir.sum.txt
  tmpf="$tmpf $odir"
  
  odir=$kname"_gscop"
  echo Rscript genomescope.R  $kname.hist $kmer $readlen1 $odir
  Rscript genomescope.R  $kname.hist $kmer $readlen1 $odir
  cp -p $odir/summary.txt  ${title}_out/$odir.sum.txt
  tmpf="$tmpf $odir"
  
  # echo DONE kmerGSE $kname  `date`

} done

echo "#TEMPS   : $tmpf ${title}_tmp/"
echo "#OUTPUTS : $outf ${title}_out/"
if [ X != "X$tmpf" ]; then mv $tmpf  ${title}_tmp/; fi
if [ X != "X$outf" ]; then mv $outf  ${title}_out/; fi
if [ $RM_TMP = 1 ]; then
  cd ${title}_tmp/; rm $tmpf; 
  cd ../
fi

#DEFER if [ -s $read1 ]; then gzip --fast $read1 ; fi

echo DONE gnodes_kmers $sid  `date`

} done
#====== end loop for SRRid, download reads =============

echo DONE gnodes_kmers $title  `date`

