#! /bin/bash
## run_gnodesqwik.sh
## --- gnodes_setup.sh for quartz.iu ---    
#SBATCH --job-name="gnodeq_pipe"
#SBATCH --output="gnodeq_pipe.%j.log"
#.. SBATCH --partition=debug
#SBATCH --partition=general
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=120G
#SBATCH --export=ALL

export NCPU=32 MEMT=120G MEMC=3G

#OPT datad=$HOME/slatedg/chrs/daphplx/gasm20set/aweed20gnodes
#OPT genefa=arath18tair1cds.fa
#OPT idclass=arath18tair1cds.idclass
#OPT reads=read22f/DRR214840_[12].fastq.gz
#OPT chrfa=iatdata/at22vcan2d_chrs.fa
#OPT anntab=iatann/at22vcan2d_nocpmt.anntab
# ^^ instead use vv
#./arath18tair_chr_nocpmtucg.anntab ./iatann/at22vcan2d_nocpmtucg.anntab ./iatann/at22vlr_ncbi_nocpmtucg.anntab

#OPT metad=iatdata/arath22asm.metad
# FIXME: opts from idclass? or metadata?
#    gidpatt==geneIDpatterm, chridpatt==chrIDpattern
# FIXME: gather gid,chrid patts from genefa,chrfa ids ?
gidpatt='AT|g'; 
chridpatt='Chr|NC_|NW_'

if [ X = "X$datad" ]; then echo datad=what; exit -1; fi
if [ X = "X$chrfa" ]; then echo chrfa=what; exit -1; fi
if [ X = "X$genefa" ]; then echo genefa=what; exit -1; fi
if [ X = "X$idclass" ]; then echo idclass=what; exit -1; fi
if [ X = "X$reads" ]; then echo reads=what; exit -1; fi
if [ X = "X$anntab" ]; then echo anntab=what; exit -1; fi
if [ X = "X$metad" ]; then echo metad=what; exit -1; fi

asmid=`basename $chrfa .fa`
genena=`basename $genefa .fa`;
if [ X = "X$title" ]; then title=${asmid}_gqwik1a; fi
if [ X = "X$rdlen" ]; then rdlen=0; fi

readna=`basename $reads | sed 's/\.gz//; s/\.[a-z0-9]*$//;' `
readna=`echo $readna | sed 's/_[12].*/ab/;'`
genebname=${genena}_${readna}_mim 
chrbname=${asmid}_${readna}_mim 
chrbam=$chrbname.bam
genebam=$genebname.bam

nbin=$HOME/bio/ncbi/bin; export PATH=$nbin:$PATH;
sambin=$HOME/bio/apps/bin; export PATH=$sambin:$PATH;
module load r

# bwabin=$HOME/bio/bwa/bin/; export PATH=$bwabin:$PATH; # export RDMAPPER=bwa-mem2
mimbin=$HOME/bio/gnomutil/flye28/bin;  export PATH=$mimbin:$PATH; # export RDMAPPER=minimap2;
if [ X = "X$mimopt" ]; then 
  mimopt="-x sr --secondary=yes -N 9 -p 0.95 "
  # mimopt="-x map-pb" for pachifi, map-ont ..
fi

export EVIGENES=$HOME/bio/apps/evigene/scripts
export PATH=$PATH:$EVIGENES
## --- end gnodes_setup.sh ---    

cd $datad
echo START gnodes_qwik $title  `date`

# if [ -s ./rmout/$asmid.out ]; then
#   echo reusing ./rmout/$asmid.out ;
# else
#   echo START_repeatmasker  `date`
#   ## need pyth for qz,rmask
#   module load python
#   hmbin=/N/soft/rhel7/hmmer/3.3.2/bin; export PATH=$PATH:$hmbin
#   rmbin=/N/soft/rhel7/repeatmasker/4.1.1; export PATH=$PATH:$rmbin
#   rmopt="-q -xsmall -libdir $HOME/.RepeatMaskerCache/rm411libdir"
#   mkdir ./rmout;
#   RepeatMasker $rmopt  -species Arabidopsis -pa $NCPU -dir ./rmout $chrfa > ./rmout/$asmid.rmlog 2>&1
#   echo DONE_repeatmasker  `date`
# fi

# if [ -s  ${gname}_buscofull_table.tsv ]; then
#   echo reusing ${gname}_buscofull_table.tsv ;
# else
#   echo START_buscoscan  `date`
#   module load busco/3.0.2
#   geneaa=`echo $genefa | sed 's/.fa/.aa/;'`
#   if [ ! -f $geneaa ]; then $EVIGENES/cdna_bestorf.pl -nostop -cdna $genefa -outaa $geneaa ; fi
#   opts=""
#   if [ X != X$BUSCO_CONFIG_FILE ]; then
#     ORDB=`grep lineage_path $BUSCO_CONFIG_FILE | sed 's/^.*= *//;'`;
#   fi
#   if [ X != X$ORDB ]; then opts="-l $ORDB/embryophyta_odb9"; fi
#   run_BUSCO.py $opts --cpu $NCPU -m prot -i $geneaa -o $gname 
#   if [ -s run_$gname/full_table_$gname.tsv ]; then cp -p run_$gname/full_table_$gname.tsv  ${gname}_buscofull_table.tsv; rm -r run_$gname/hmmer_output; fi
#   echo DONE_buscoscan  `date`
# fi

# samtools depth 
sdopt="-H -aa -J -q 10 -Q 0 --threads $NCPU"
tmpf=""
outf=""

if [ -s $genebam -a -s ${genebname}.stmdepth.tab ]; then 
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
    
  # erase big tmp files now, genebam still needed
  rm ${genebname}_stemp.bam
  tmpf="$tmpf $genebam ${genebname}.log"
  tmpf="$tmpf ${genebname}.sdupcov.txt ${genebname}.stdcov.txt"
  outf="$outf ${genebname}.stmdepth.tab $genebname.fstats.txt"

  echo DONE cdsmap `date`
fi

# FIXME rdlen: should use ave(all), and sum(all) for gnosize
if [ $rdlen = 0 ]; then
  rdlen=`samtools view -F 0xF04 $genebam | cut -f10 | head -50000 |\
   perl -ne '($s)=split; $w=length($s);if($w>1){ $sw+=$w; $n++; } END{ printf "%.0f",($n<1?0:$sw/$n); }' `;  
fi
  
# do after cdsbam/depth.tab; 
if [ -s $genebname.stmdepth.tab -a ! -s $genebname.stdt_genexcopy  ]; then  
  echo START genecov  `date`
  
  $EVIGENES/genoasm/gnodes_std2genecov.pl \
   -idclass $idclass -rdlen $rdlen  \
   -in  $genebname.stmdepth.tab  -stats  $genebname.fstats.txt \
   -out $genebname.stdt_genexcopy
  
  # this is same for chr,cds.stmdepth.tab; should add read map stats to std.covtab
  perl -ne \
'BEGIN{ $BN=100; } ($cr,$cb,$dat,$dam,$dgn,$dte,$dun)=split; map{ $_||=0; }($dgn,$dte,$dun); $bb=int($cb/$BN); 
if($bb > $lbb or $cr ne $lcr) { putv($lcr,$lbb) if($lcr); @sv=(0) x 4; } 
$dau=($dam == $dat) ? $dam : $dam - ($dat - $dam); $dau=0 if($dau<0); @v=($dgn,$dte,$dun,$dat,$dam,$dau); 
for $i (0..5) { $sv[$i]+=$v[$i]; } $lcr=$cr; $lbb=$bb; $lb=$cb; END{ putv($lcr,$lbb); }  
sub putv{ map{ $_=int(0.5+$_/$BN) } @sv; print join("\t",$lcr,100*(1+$lbb),@sv)."\n"; }
BEGIN{ print "#".join("\t",qw(ChrID Pos cdsCov teCov unkCov aCovT aCovM aCovU))."\n";}' \
 ${genebname}.stmdepth.tab > ${genebname}.std1a.covtab

  env rdlen=$rdlen name=$genebnam perl -ne 'my($n,@v)= split; if($n<1){ } 
elsif(m/ primary$/ and not $nrd){ $nrd=$n; } elsif(m/ primary mapped / and not $nmap){ $nmap=$n; } 
END{ $nno=$nrd-$nmap; $rdlen= ($l=$ENV{rdlen})? "readlen=$l," :""; $nam=$ENV{name}||""; 
print "#pt1.$nam $rdlen n_readid $nrd, n_nomap $nno, n_mapok $nmap\n" if($nrd); } ' \
    $genebname.fstats.txt >> ${genebname}.std1a.covtab
 
  outf="$outf $genebname.stdt_genexcopy ${genebname}.std1a.covtab"
  echo DONE genecov  `date`
fi

# genedata used below: $genebname.stdt_genexcopy $genebname.std1a.covtab $genebam

if [ -s $chrbam -a -s ${chrbname}.stmdepth.tab ]; then 
  echo reusing $chrbam;
else
  echo START chrmap  `date`

  if [ ! -s $chrbam ]; then
  ( minimap2 $mimopt -a -t $NCPU $chrfa $reads | samtools view --threads $NCPU -Sb -o $chrbam - ) > ${chrbname}.log 2>&1
  fi

  samtools flagstats --threads $NCPU $chrbam  > $chrbname.fstats.txt
  samtools sort -m $MEMC --threads $NCPU -o ${chrbname}_stemp.bam ${chrbam}
  
  samtools depth -G 0x100 $sdopt -o ${chrbname}.stdcov.txt ${chrbname}_stemp.bam
  samtools depth -g 0x100 $sdopt -o ${chrbname}.sdupcov.txt ${chrbname}_stemp.bam
  
  # erase big tmp files now
  rm ${chrbname}_stemp.bam
  tmpf="$tmpf $chrbam ${chrbname}.log"
  
  #  genexchr depth calc before chr.stmdepth.tab paste
  mname=${chrbname}_$genena

  samtools merge --threads $NCPU -n -o $mname.bam $genebam $chrbam
  samtools collate -f --threads $NCPU -o $mname.col.bam  $mname.bam
  rm  $mname.bam
  tmpf="$tmpf $mname.col.bam "
  
  samtools view -H $chrbam > $mname.genexchr.sam
  
  # gnoctab=$mname.genenochr.tab: table of geneid, gposition? nreads hit gene not chr
  # FIXME: opts gpat==geneIDpatterm, cpat==chrIDpattern
  
  samtools view $mname.col.bam | \
  env gpat=$gidpatt cpat=$chridpatt gcout=$mname.genexchr.sam gnoctab=$mname.genenochr.tab name=$mname perl -ne \
'BEGIN { $CP=$ENV{cpat}||"[Cc]hr"; $GP=$ENV{gpat}||"g"; 
if($gco=$ENV{gcout}) { $GCOK=open(S,">>$gco"); }
if($gnoc=$ENV{gnoctab}) { $GNOC=open(G,">$gnoc"); } } 
$insam=$_; ($rid,$fl,$cgid,$cb,$mq,$cig)=split; $rid .= ".2" if($fl & 0x80); 
if($lrid ne $rid){ putv($lrid) if($lrid);  $nm=$nno=$ngp=$ncp=$galn=$gid=$crsam=0; } 
if($fl & 0x04) { $nno++; } else { $nm++; if($cgid=~/$CP/){ $ncp++; $crsam=$insam; } 
elsif($cgid =~ /$GP/) { $ngp++; $gid=$cgid; $galn=caln($cig); } }  $lrid=$rid; 
END{ putv($lrid); putgnoc(); putsum(); } 
sub caln{ my($c)=@_; my $al=0; while($c=~m/(\d+)[MD]/g){ $al+=$1; } return($al); }
sub putv{ $snm+=$nm; $srd++; $snno+=$nno; $sng+=$ngp; $snc+=$ncp; 
$sgnoc++ if($ngp and not $ncp); $scnog++ if($ncp and not $ngp); 
if($GCOK and $crsam and $ngp){ print S $crsam; } 
if($ngp) { if($ncp) { $sngc++; $sagc+=$galn; } else { $gnmiss{$gid}{n}++; $gnmiss{$gid}{a}+=$galn; } } } 
sub putgnoc{ my($sn,$sa,$nmis)=(0,0,0); 
for $id (sort keys %gnmiss){ my($n,$a)= map{ $gnmiss{$id}{$_}||0 } qw(n a); 
print G join("\t",$id,$n,$a)."\n"; $sn+=$n; $sa+=$a; $nmis++; }
print G join("\t","total_nochr",$sn,$sa,$nmis)."\n";
print G join("\t","total_onchr",$sngc,$sagc)."\n"; 
printf G "pcent_nochr\t%.2f\t%.2f\n",100*$sn/$sngc,100*$sa/$sagc if($sngc>0 and $sagc>0); }
sub putsum{ map{ $_ ||= 0 } ($sgnoc,$scnog);  
warn "# nread=$srd, mapped=$snm, nomap=$snno, genemap=$sng, chrmap=$snc, genenochr=$sgnoc, chrnogene=$scnog\t$ENV{name}\n"; }' 
    
  samtools sort -m $MEMC --threads $NCPU -o $mname.genexchrs.bam $mname.genexchr.sam 
  samtools depth $sdopt  -o $mname.scdscov.txt $mname.genexchrs.bam

  rm $mname.genexchr.sam     
  tmpf="$tmpf $mname.genexchrs.bam"
  outf="$outf $mname.genenochr.tab"
 
  paste ${chrbname}.sdupcov.txt ${chrbname}.stdcov.txt $mname.scdscov.txt | cut -f1,2,3,6,9 | \
    perl -pe 'if(/^#CHR/){ $_="#ChrID\tPos\tCovT\tCovM\tCovCDS\n"; }' > ${chrbname}.stmdepth.tab

  tmpf="$tmpf ${chrbname}.sdupcov.txt ${chrbname}.stdcov.txt $mname.scdscov.txt"
  outf="$outf ${chrbname}.stmdepth.tab $chrbname.fstats.txt"

  echo DONE chrmap  `date`
fi

if [ -s ${chrbname}.stmdepth.tab ]; then
  # stmdepth.tab to gnodes.covtab
  # change orig hdr: ChrID Pos CovT CovM CovU aCovT aCovM aCovU
  # to new hdr?  ChrID Pos cdsCov teCov unkCov aCovT aCovM aCovU
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
  
  outf="$outf ${chrbname}.std1a.covtab"
fi

# if [ -s $anntab ]; then
#   echo reusing $anntab ;
# else
#   echo START_gnodes2_annotate  `date`
#   $EVIGENES/genoasm/gnodes_annotate.pl -debug -ncpu $NCPU -cds $genefa -te ./rmout/$asmid.out -idclass $idclass \
#     -busco ${gname}_buscofull_table.tsv -genecovtab $genebname.stdt_genexcopy -chr $chrfa -output $anntab
#   echo DONE_gnodes2_annotate  `date`
# fi

# NOCPMITO fix
# annold=at22vcan2d_chrs.anntab
# anntab=at22vcan2d_nocpmt.anntab
# echo NOCPMITO $anntab
# 
# if [ ! -s $anntab ]; then
#  perl -p -e 'if(/^\w*ChrMT/){ $ctype="contam,MT,"; } elsif(/^\w*ChrPltd/){ $ctype="contam,Pltd,"; } 
# else { $ctype=""; } if($ctype){ ($sc,$sb,$an,$ids)=@v=split; $v[2]=$ctype.$an; 
# if($lsc eq $sc and $lsb < $sb - 100){ fillspan($sc,$lsb,$sb,$ctype); } 
# elsif($lsc ne $sc and $sb > 100){ fillspan($sc,0,$sb,$ctype); } 
# $lsc=$sc; $lsb=$sb; $_=join("\t",@v)."\n"; } $lsc=$sc;  
# sub fillspan{my($sc,$lsb,$sb,$ctype)=@_; for ($i=$lsb+100; $i<$sb; $i+=100){ 
#   print join("\t",$sc,$i,$ctype,"na")."\n"; } } ' \
#   $annold > $anntab
# 
# fi
# 
# echo "NOCPMITO contam= ",  `grep -c contam $anntab`

if [ -s ${chrbname}.std1a.covtab ]; then
  echo START covsum  `date`
  $EVIGENES/genoasm/gnodes_covsum.pl  -asmid $asmid -title ${title}  -anntab $anntab -crclass $idclass  -sumdata $metad  -debug \
    -genexcopy $genebname.stdt_genexcopy -output ${chrbname}.std1a_sum.txt ${chrbname}.std1a.covtab

  echo DONE covsum  `date`
fi

if [ -s ${chrbname}.std1a.covtab ]; then
  echo START chrcovplot  `date`
  $EVIGENES/genoasm/gnodes_covsum.pl -plotchr -asmid $asmid -title ${title} -anntab $anntab -crclass $idclass  -sumdata $metad  -debug \
    -genexcopy $genebname.stdt_genexcopy  ${chrbname}.std1a.covtab

  echo DONE chrcovplot  `date`
fi

if [ -s ${chrbname}.std1a.genetab  ]; then
  echo START sumgenecov  `date`
  $EVIGENES/genoasm/gnodes_sumgenecov.pl  -title ${title}  -debug \
   -genexcopy $genebname.stdt_genexcopy -chrgenetab ${chrbname}.std1a.genetab \
    -chrsum ${chrbname}.std1a_sum.txt -cdscov ${genebname}.std1a.covtab 

  echo DONE sumgenecov  `date`
fi

outf="$outf ${chrbname}.std1a*  ${title}*"
echo "#TEMPfiles: $tmpf"
echo "#OUTPUTS  : $outf"
# mkdir ${title}_tmp; mv $tmpf  ${title}_tmp/
# mkdir ${title}_out; mv $outf  ${title}_out/

echo DONE gnodes_qwik $title  `date`
