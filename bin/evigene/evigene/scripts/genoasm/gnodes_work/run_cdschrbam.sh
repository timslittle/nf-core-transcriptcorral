#! /bin/bash
## run_sambamstats.sh : samtools flagstat / stats of in.bam
## --- gnodes_setup.sh for quartz.iu ---    
#SBATCH --job-name="gnodes_pipe"
#SBATCH --output="gnodes_pipe.%j.log"
#.. SBATCH --partition=general
#SBATCH --partition=debug
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=120G
#SBATCH --export=ALL

export NCPU=32; export TMEM=120G; MEM=3G;

datad=/N/slate/gilbertd/chrs/daphplx/gasm20set/aweed20gnodes

if [ X = "X$gnbam" ]; then echo "gnbam=what?"; exit -1; fi
if [ X = "X$chrbam" ]; then echo "chrbam=what?"; exit -1; fi
if [ X = "X$datad" ]; then echo "datad=what?"; exit -1; fi
gnna=`basename $gnbam .bam`
chrna=`basename $chrbam .bam`
gnn=`echo $gnna | sed 's/_.*//; s/-.*//;'`
mname=${gnn}_$chrna

nbin=$HOME/bio/ncbi/bin; export PATH=$nbin:$PATH;
sambin=$HOME/bio/apps/bin; export PATH=$sambin:$PATH;
export EVIGENES=$HOME/bio/apps/evigene/scripts
export PATH=$PATH:$EVIGENES

cd $datad
echo START `date`

sdopt="--threads $NCPU"
# samtools flagstats $sdopt $gnbam > $gnna.st_fstats.txt
# samtools flagstats $sdopt $chrbam > $chrna.st_fstats.txt

## * Cat no good, neeed merge
# BUG is header, cat uses 1st only unless -h new.hdr
# samtools view -H $gnbam > $mname.hdr
# samtools view -H $chrbam >> $mname.hdr
# samtools cat $sdopt -h $mname.hdr -o $mname.bam $gnbam $chrbam

samtools merge $sdopt -n -o $mname.bam $gnbam $chrbam

# -f == fast drops 2ndary map
samtools collate -f $sdopt -o $mname.col.bam  $mname.bam

# other way, maybe want this
# samtools sort -n  $sdopt -m $MEM -o $mname.col.bam  $mname.bam

samtools flagstats $sdopt $mname.col.bam  > $mname.col_fstats.txt

## collect stats on gene x chr, AND output new bam w/ only chr ids that have gene also, for depth stat
## count reads w/ gene not chr, and vv, and neither
# samtools view $mname.col.bam | env crpat=Chr gnpat=AT perl -ne '($rid)=@v=split; if($lrid ne $rid) { } $lrid=$rid;'

echo "#TMP rm $mname.bam"
echo DONE `date`
exit;

#----------
## upd this to put chr_gene.sam for putv($inline) if(cgid =~ /$CP/ and $ngp)
# samtools view arath18tair1cds_at22vlr_ncbi_chrs-DRR214840_1_mim.col.bam | cut -f1-9 | env gpat=AT cpat=Chr perl -ne \
# 'BEGIN{ $CP=$ENV{cpat}||"[Cc]hr"; $GP=$ENV{gpat}||"g"; } ($rid,$fl,$cgid,$db,$mq,$cig)=@v=split; $rid .= ".2" if($fl & 0x80);  if($lrid ne $rid){ putv($lrid) if($lrid);  $nm=$nno=$ngp=$ncp=0; } $nno++ if($fl & 0x04);  $ngp++ if($cgid=~/$GP/); $ncp++ if($cgid =~ /$CP/); $nm++;   $lrid=$rid; END{ putv($lrid); putsum(); } sub putv{ $snm+=$nm; $srd++; $snno+=$nno; $sng+=$ngp; $snc+=$ncp; $sgnoc++ if($ngp and not $ncp); $scnog++ if($ncp and not $ngp); } sub putsum{ map{ $_ ||= 0 } ($sgnoc,$scnog);  print "srd=$srd, snm=$snm, snno=$snno, sng=$sng, snc=$snc, sgnoc=$sgnoc, scnog=$scnog\n";}' | head
# 
# srd=80316418, snm=160632836, snno=52788091, sng=35680670, snc=79643018, sgnoc=28074, scnog=43990422
#   note snno includes gene no-map where chr-hasmap, n_rdid_map total should be snc + sgnoc

#.... update
# samtools view arath18tair1cds_at22vlr_ncbi_chrs-DRR214840_1_mim.col.bam | head -220000 | \
#  env pr=1  gpat=AT cpat=Chr perl -ne \
# 'BEGIN{ $CP=$ENV{cpat}||"[Cc]hr"; $GP=$ENV{gpat}||"g"; } ($rid,$fl,$cgid,$db,$mq,$cig)=@v=split; $rid .= ".2" if($fl & 0x80);  if($lrid ne $rid){ putv($lrid) if($lrid);  $nm=$nno=$ngp=$ncp=$incr=0; } if($fl & 0x04) { $nno++; } else { $nm++; if($cgid=~/$CP/){ $ncp++; $incr=$_; } elsif($cgid =~ /$GP/) { $ngp++; } }   $lrid=$rid; END{ putv($lrid); putsum(); } sub putv{ $snm+=$nm; $srd++; $snno+=$nno; $sng+=$ngp; $snc+=$ncp; $sgnoc++ if($ngp and not $ncp); $scnog++ if($ncp and not $ngp); if($ENV{pr} and $incr and $ngp){ print $incr; }  } sub putsum{ map{ $_ ||= 0 } ($sgnoc,$scnog);  warn "# srd=$srd, snm=$snm, snno=$snno, sng=$sng, snc=$snc, sgnoc=$sgnoc, scnog=$scnog\n";}' | head
#
# srd=110000, snm=147167, snno=72833, sng=38930, snc=108237, sgnoc=168, scnog=69475
#   snm = sng + snc; 
# output chrxgene.sam | samtools sort -o chrxgene_sort.bam ; samtools depth -o chrxgene.stcov.txt chrxgene_sort.bam

