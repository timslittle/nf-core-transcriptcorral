#! /bin/bash
## run_rdecontam.sh : filter readset by contam.fa
#SBATCH --job-name="gnodes_pipe"
#SBATCH --output="gnodes_pipe.%j.log"
#.. SBATCH --partition=debug
#SBATCH --partition=general
#SBATCH -t 8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=120G
#SBATCH --export=ALL

## TEST effect of map to contamonly.fa vs allchrasm.fa w/ contam marked
## ie arath has MT and Cplast genes on nuclr chrs, likely some from chrs map to contamonly.fa

export NCPU=32 MEMT=120G MEMC=3G

if [ X = "X$datad" ]; then echo datad=what; exit -1; fi
if [ X = "X$contam" ]; then echo contam=what contaminants.fa ; exit -1; fi
if [ X = "X$outpre" ]; then outpre=decon; fi
if [ X = "X$nogz" ]; then nogz=0; fi

if [ X = "X$intab" ]; then intab=at1k_sraruninfo1a.tsv; fi
if [ X = "X$srrid" -a -s $intab ]; then
  srrid=`cut -f1 $intab | egrep '^[DES]RR'`
fi

useid=1
if [ X != "X$reads" ]; then  srrid=$reads; useid=0; fi
if [ X = "X$readsf" ]; then  readsf="."; fi
if [ X = "X$srrid" ]; then echo Need srrid= or reads=; exit -1; fi

sambin=$HOME/bio/apps/bin; export PATH=$sambin:$PATH;
mimbin=$HOME/bio/gnomutil/flye28/bin;  export PATH=$mimbin:$PATH;  
mimopt="-x sr --secondary=no "

cd $datad
echo START gnodes_rdecontam $title  `date`

irun=0;
for sid in $srrid; do {

if [ $useid = 1 ]; then  
  read1=$readsf/${sid}_1.fastq; 
  reads=$readsf/${sid}_[12].fastq; 
  if [ -s $read1.gz ]; then
    reads=$readsf/${sid}_[12].fastq.gz; 
  elif [ ! -s $read1 ]; then
    echo Missing $read1; continue
  fi
fi

rname=`basename $reads | sed 's/\.gz//; s/.fastq//; s/_1//;'`
oname=$outpre$rname
if [ -s ${oname}_1.fastq -o -s ${oname}_1.fastq.gz ]; then
  echo reusing ${oname}_1.fastq; continue;
fi

minimap2 $mimopt -a -t $NCPU $contam $reads | grep -v '^@' | env oname="$oname" perl -ne \
'BEGIN { $onam=$ENV{oname}||"noname"; $ridsuf=0;
$fn= $onam."_1.fastq";  if(-f $fn){ $onam=$onam."2de"; $fn= $onam."_1.fastq"; }
warn "# decontam: output to ",$onam,"_[12].fastq\n";
open(F1,">",$fn) or die "writing $fn";
$fn= $onam."_2.fastq"; open(F2,">",$fn) or die "writing $fn";  } 
($rd,$fl)=@v=split;  
unless($fl == 77 or $fl == 141){ $ncon++ unless($fl & 0x04); $nmiss1++ unless($fl & 0x08); next; }
if($ridsuf==0 and $fl == 141) { $ridsuf=($rd =~ m,/\d$,)?1:-1; }
if($ridsuf==1) { ($rdc=$rd) =~ s,/\d$,,; $newid=($rdc ne $lrdc); $lrdc= $rdc; } 
else { $newid= ($rd ne $lrd); }
if($newid){ putq($lrd,@rs) if(@rs>1); @rs=(); } 
if($fl & 0x04){ $rb=($fl & 0x80)?2:0; push @rs, @v[9,10],$rb; $nok++; } 
else { $ncon++; } $lrd=$rd; 
END { putq($lrd,@rs) if(@rs>1); 
$cnt= sprintf "%.0f\t%d\t%d\t%s",($nr<1?0:$sw/$nr),$nr,$sw,$onam;
open(FC,">$onam.count"); print FC $cnt,"\n"; close(FC); close(F1); close(F2);
warn "# decontam npairs=$npairs, ncontam=$ncon\n"; warn "# $cnt\n"; }
sub putq{ my($rid,@rv)=@_; if(@rv<5){ $nmiss1++; return; }
  if($rv[2] == 2){ @rv=@rv[3,4,5, 0,1,2]; }
  my($rs,$rq,$ra,$rsb,$rqb,$rb)=@rv; $npairs++;
  print F1 "\@$rid\n",$rs,"\n\+\n",$rq,"\n"; $nr++; $sw+=length($rs);
  print F2  "\@$rid\n",$rsb,"\n\+\n",$rqb,"\n"; $nr++; $sw+=length($rsb); 
} ' 

if [ -s ${oname}_1.fastq -a $nogz = 0 ]; then
  gzip --fast ${oname}_1.fastq &
  gzip --fast ${oname}_2.fastq &
  wait
fi

if [ $useid = 0 ]; then break; fi

} done

echo END gnodes_rdecontam $title  `date`

# SAMFLAG_nomap => 0x04, SAMFLAG_rev => 0x10, read1=>0x40? SAMFLAG_read2 => 0x80, SAMFLAG_2ndary => 0x100, SAMFLAG_suppl => 0x800
## need only two flags: 77, 141 == PAIRED,UNMAP,MUNMAP,READ1/2 : egrep ' (77|141)  ' << Unsafe tabs/spc
# 0x4d	77	PAIRED,UNMAP,MUNMAP,READ1
# 0x8d	141	PAIRED,UNMAP,MUNMAP,READ2
#?? check rid suffix /1 /2 and chomp ?

