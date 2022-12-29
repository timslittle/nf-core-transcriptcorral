#! /bin/bash
## run_gnodesqwik_at1k.sh
#SBATCH --job-name="gnodeq_pipe"
#SBATCH --output="gnodeq_pipe.%j.log"
#.. SBATCH --partition=debug
#SBATCH --partition=general
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=120G
#SBATCH --export=ALL

## change mimopt: -p 0.95 bad?, leave default 0.80
# run_gnodesqwik_at1knu.sh add run_fastuniq.sh

export NCPU=32 MEMT=120G MEMC=3G
if [ X = "X$TESTLOOP" ]; then TESTLOOP=0; fi

if [ X = "X$intab" ]; then intab=at1k_sraruninfo1a.tsv; fi
if [ X = "X$srrid" -a -s $intab ]; then
  srrid=`cut -f1 $intab | egrep '^[DES]RR'`
fi
if [ X = "X$srrid" ]; then echo srrid=what; exit -1; fi
if [ X = "X$maxrun" ]; then maxrun=32; fi
if [ X = "X$dofu" ]; then dofu=0; fi
if [ X = "X$RM_TMP" ]; then RM_TMP=0; fi
if [ X = "X$ttag" ]; then ttag="gnoq_at1sw"; fi

genefa=arath18tair1cds.fa
idclass=arath18tair1cds.idclass
chrfa=arath18tair_chr.fa
anntab=arath18tair_chr.anntab
metad=arath20asm.metad
gidpatt=AT;  chridpatt='Chr|NC_';

if [ X = "X$datad" ]; then echo datad=what; exit -1; fi
if [ X = "X$chrfa" ]; then echo chrfa=what; exit -1; fi
if [ X = "X$genefa" ]; then echo genefa=what; exit -1; fi
if [ X = "X$idclass" ]; then echo idclass=what; exit -1; fi
if [ X = "X$anntab" ]; then echo anntab=what; exit -1; fi
if [ X = "X$metad" ]; then echo metad=what; exit -1; fi

asmid=`basename $chrfa .fa`
genena=`basename $genefa .fa`;
if [ X = "X$title" ]; then title=${asmid}_$ttag; fi

if [ X = "X$rdlen" ]; then rdlen=0; fi
# ^?? ODD BUG got rdlen=110 instead of 76 below from samtools view col 10.. 
# this run only: env intab=at1k_gsizes_fcswens4r.ab dofu=1 RM_TMP=1 ttag=at1sw datad=`pwd` sbatch run_gnodesqwik_at1ksw.sh
# [gilbertd@h1 gnoat1k]$ grep rdlen arath18tair_chr_at1sw_out/gnodeq_pipe.623924.log | head
#   for  Nr.allcds=16324272 (38.6%), Nr.total=42288756, Lrdlen=110, MapErr=0 (0%), Nshort=0


nbin=$HOME/bio/ncbi/bin; export PATH=$nbin:$PATH;
sambin=$HOME/bio/apps/bin; export PATH=$sambin:$PATH;
sradumpr=$HOME/bio/sratoolkit/sratoolkit300/bin/fasterq-dump
fubin=$HOME/bio/gnomutil/bin/

mimbin=$HOME/bio/gnomutil/flye28/bin;  export PATH=$mimbin:$PATH; # export RDMAPPER=minimap2;
if [ X = "X$mimopt" ]; then 
  #x mimopt="-x sr --secondary=yes -N 9 -p 0.95 "
  mimopt="-x sr --secondary=yes -N 9 -p 0.80 "
  # mimopt="-x map-pb" for pachifi, map-ont ..
fi

# samtools depth 
sdopt="-H -aa -J -q 10 -Q 0 --threads $NCPU"

export EVIGENES=$HOME/bio/apps/evigene/scripts
export PATH=$PATH:$EVIGENES
## --- end gnodes_setup.sh ---    

cd $datad
echo START gnodes_qwik $title  `date`

mkdir  ${title}_tmp/
mkdir  ${title}_out/
# mkdir  ${title}_data/
# replace  ${title}_data/ with option
if [ X = "X$readsf" ]; then  readsf=${title}_data; fi
if [ ! -d $readsf ]; then  mkdir  $readsf; fi

#====== loop here for SRRid, download reads =============

irun=0;
for sid in $srrid; do {

if [ $irun -ge $maxrun ]; then echo "DONE gnodes_qwik $irun/$maxrun runs.. bye."; break; fi

echo START gnodes_qwik $sid  `date`
sinfo=`grep $sid $intab`
echo "INFO" $sinfo

# step0. check if done SRRid (sraid?); need readna,chrbname here before reads exist
readna=${sid}ab
chrbname=${asmid}_${readna}_mim 

testf=${title}_out/${chrbname}$ttag.std1a_sum.txt
if [ -s $testf ]; then echo "# DONE $sid in $testf; continue.. ";  continue; fi

testf=$readsf/sradump.$sid.log
# if [ -s $testf ]; then echo "# IN_PROGRESS? $sid in $testf; continue.. "; continue; fi

irun=$(($irun + 1)); 
if [ $TESTLOOP = 1 ]; then  echo "#TEST skip run.. "; continue; fi
 
read1=$readsf/${sid}_1.fastq
read2=$readsf/${sid}_2.fastq

if [ -s $read1.gz ]; then 
  gunzip $read1.gz &
  gunzip $read2.gz &
  wait
fi
if [ -s $read1 ]; then 
  echo reusing $read1
else
  $sradumpr --threads $NCPU --outdir $readsf --split-3 $sid >& $readsf/sradump.$sid.log 
fi
if [ ! -s $read1 ]; then  echo "FAIL sra data: $read1"; exit -1; fi

reads=$readsf/${sid}_[12].fastq
readna=`basename $reads | sed 's/\.gz//; s/\.[a-z0-9]*$//;' `
readna=`echo $readna | sed 's/_[12].*/ab/;'`

## insert here fastuniq.sh ; calc stats for both?? fu_SRR bad, _ cut in sum.txt

if [ $dofu = 1 ]; then
  outfu=`echo fu$readna | sed 's/ab$//;'`
  touch $outfu.list
  for rn in $reads; do { echo $rn >> $outfu.list; } done
  echo $fubin/fastuniq -c 0 -t q -i $outfu.list -o $readsf/${outfu}_1.fastq -p $readsf/${outfu}_2.fastq
  $fubin/fastuniq -c 0 -t q -i $outfu.list -o $readsf/${outfu}_1.fastq -p $readsf/${outfu}_2.fastq
  
  ic=`wc -l $read1 | sed 's/ .*//;'`
  oc=`wc -l $readsf/${outfu}_1.fastq | sed 's/ .*//;'`
  echo fastuniq reduction 
  echo "$ic/4; $oc/4; scale=2; 100*$oc/$ic" | bc
  rm $outfu.list
  
  # replace $reads w/ outfu for map
  oreads=$reads
  gzip --fast $reads &
  read1=$readsf/${outfu}_1.fastq
  read2=$readsf/${outfu}_2.fastq
  reads=$readsf/${outfu}_[12].fastq
  readna=`basename $reads | sed 's/\.gz//; s/\.[a-z0-9]*$//;' `
  readna=`echo $readna | sed 's/_[12].*/ab/;'`
  echo "# replaced $oreads with $reads"
fi

genebname=${genena}_${readna}_mim 
chrbname=${asmid}_${readna}_mim 
chrbam=$chrbname.bam
genebam=$genebname.bam
tmpf=""; outf=""

# always?  $rdlen = 0 FIXME rdlen: should use ave(all), and sum(all) for gnosize
if [ 1 = 1 ]; then
  frdcount=`echo $read1 | sed 's/_1.fastq/.count/;'`
  if [ ! -s $frdcount ]; then
  perl -ne 'if($i++ % 4 == 1){ $w=length; if($w>1){ $sw+= $w -1; $n++; }} END{ printf "%.0f\t%d\t%d\n",($n<1?0:$sw/$n),$n,$sw; } ' \
    $reads > $frdcount;
  cp -p $frdcount ${title}_out/
  echo -n "$frdcount "; head $frdcount
  fi
  rdlen=`cut -f1 $frdcount`;
  rdcount=`cut -f2 $frdcount`;
  rdbases=`cut -f3 $frdcount`;  
fi

if [ -s ${genebname}.stmdepth.tab ]; then 
  echo reusing $genebam;
else
  echo START cdsmap  `date`

  if [ ! -s $genebam ]; then
  ( minimap2 $mimopt -a -t $NCPU $genefa $reads | samtools view --threads $NCPU -Sb -o $genebam - ) > ${genebname}.log 2>&1
  fi
  
  samtools flagstats --threads $NCPU $genebam  > $genebname.fstats.txt
  samtools sort -m $MEMC --threads $NCPU -o ${genebname}_stemp.bam ${genebam}
  
  samtools depth -G 0x100 $sdopt -o ${genebname}.stdcov.txt ${genebname}_stemp.bam
  samtools depth -g 0x100 $sdopt -o ${genebname}.sdupcov.txt  ${genebname}_stemp.bam
  paste ${genebname}.sdupcov.txt ${genebname}.stdcov.txt | cut -f1,2,3,6 | \
    perl -pe 'if(/^#CHR/){ $_="#ChrID\tPos\tCovT\tCovM\n"; }' > ${genebname}.stmdepth.tab
    
  rm ${genebname}_stemp.bam  ${genebname}.sdupcov.txt ${genebname}.stdcov.txt
  rm $genebam
  tmpf="$tmpf ${genebname}.log  ${genebname}.stmdepth.tab"
  outf="$outf $genebname.fstats.txt"

  echo DONE cdsmap `date`
fi

if [ -s $genebname.stmdepth.tab -a ! -s $genebname.stdt_genexcopy  ]; then  
  echo START genecov  `date`
  
  $EVIGENES/genoasm/gnodes_std2genecov.pl \
   -idclass $idclass -rdlen $rdlen  \
   -in  $genebname.stmdepth.tab  -stats  $genebname.fstats.txt \
   -out $genebname.stdt_genexcopy
   
  # not at1k: ${genebname}.std1a.covtab 
  tmpf="$tmpf $genebname.stmdepth.tab"
  outf="$outf $genebname.stdt_genexcopy "
  echo DONE genecov  `date`
fi

if [ -s ${chrbname}.stmdepth.tab ]; then 
  echo reusing $chrbam;
else
  echo START chrmap  `date`

  if [ ! -s $chrbam ]; then
  ( minimap2 $mimopt -a -t $NCPU $chrfa $reads | samtools view --threads $NCPU -Sb -o $chrbam - ) > ${chrbname}.log 2>&1
  
  # done w/ reads, squeeze or erase
  gzip --fast $reads &
  fi

  samtools flagstats --threads $NCPU $chrbam  > $chrbname.fstats.txt
  samtools sort -m $MEMC --threads $NCPU -o ${chrbname}_stemp.bam ${chrbam}
  
  samtools depth -G 0x100 $sdopt -o ${chrbname}.stdcov.txt ${chrbname}_stemp.bam
  samtools depth -g 0x100 $sdopt -o ${chrbname}.sdupcov.txt ${chrbname}_stemp.bam
  
  # erase big tmp files now
  rm ${chrbname}_stemp.bam
  rm $chrbam
  tmpf="$tmpf ${chrbname}.log"

  # not at1k: $mname.scdscov.txt 
  paste ${chrbname}.sdupcov.txt ${chrbname}.stdcov.txt  | cut -f1,2,3,6 | \
    perl -pe 'if(/^#CHR/){ $_="#ChrID\tPos\tCovT\tCovM\n"; }' > ${chrbname}.stmdepth.tab

  rm ${chrbname}.sdupcov.txt ${chrbname}.stdcov.txt
  tmpf="$tmpf ${chrbname}.stmdepth.tab"
  outf="$outf $chrbname.fstats.txt"

  echo DONE chrmap  `date`
fi

if [ -s ${chrbname}.stmdepth.tab ]; then
  # stmdepth.tab to gnodes.covtab
  perl -ne \
'BEGIN{ $BN=100; } ($cr,$cb,$dat,$dam,$dgn,$dte,$dun)=split; map{ $_||=0; }($dgn,$dte,$dun); $bb=int($cb/$BN); 
if($bb > $lbb or $cr ne $lcr) { putv($lcr,$lbb) if($lcr); @sv=(0) x 4; } 
$dau=($dam == $dat) ? $dam : $dam - ($dat - $dam); $dau=0 if($dau<0); @v=($dgn,$dte,$dun,$dat,$dam,$dau); 
for $i (0..5) { $sv[$i]+=$v[$i]; } $lcr=$cr; $lbb=$bb; $lb=$cb; END{ putv($lcr,$lbb); }  
sub putv{ map{ $_=int(0.5+$_/$BN) } @sv; print join("\t",$lcr,100*(1+$lbb),@sv)."\n"; }
BEGIN{ print "#".join("\t",qw(ChrID Pos cdsCov teCov unkCov aCovT aCovM aCovU))."\n";}' \
 ${chrbname}.stmdepth.tab > ${chrbname}.std1a.covtab
  
  env rdlen=$rdlen name=$chrbname perl -ne 'my($n,@v)= split; if($n<1){ } 
elsif(m/ primary$/ and not $nrd){ $nrd=$n; } elsif(m/ primary mapped / and not $nmap){ $nmap=$n; } 
END{ $nno=$nrd-$nmap; $rdlen= ($l=$ENV{rdlen})? "readlen=$l," :""; $nam=$ENV{name}||""; 
print "#pt1.$nam $rdlen n_readid $nrd, n_nomap $nno, n_mapok $nmap\n" if($nrd); } ' \
    $chrbname.fstats.txt >> ${chrbname}.std1a.covtab
  
  tmpf="$tmpf ${chrbname}.stmdepth.tab"
  outf="$outf ${chrbname}.std1a.covtab"
fi

if [ -s ${chrbname}.std1a.covtab ]; then
  echo START covsum  `date`
  
  $EVIGENES/genoasm/gnodes_covsum.pl  -asmid $asmid -title ${title}  -anntab $anntab -crclass $idclass  -sumdata $metad  -debug \
    -genexcopy $genebname.stdt_genexcopy -output ${chrbname}$ttag.std1a_sum.txt ${chrbname}.std1a.covtab

  echo DONE covsum  `date`
fi

outf="$outf ${chrbname}$ttag.std1a*"
echo "#TEMPS   : $tmpf ${title}_tmp/"
echo "#OUTPUTS : $outf ${title}_out/"
mv $tmpf  ${title}_tmp/
mv $outf  ${title}_out/
if [ $RM_TMP = 1 ]; then
  cd ${title}_tmp/; rm $tmpf; 
  cd ../
  #and?? rm $reads
fi

echo DONE gnodes_qwik $sid  `date`

} done
#====== end loop for SRRid, download reads =============

echo DONE gnodes_qwik $title  `date`

