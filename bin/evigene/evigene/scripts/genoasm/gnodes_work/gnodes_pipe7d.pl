#!/usr/bin/perl
# gnodes_pipe.pl for evigene/scripts/genoasm/ 

=item usage gnodes3_covsum

  gnodes_pipe.pl  -chr chr.fasta -cds cds.fasta -te te.fasta 
    opts: -genomedata  daphnia3genomes.data
    opts: -idclasses cdste.idclass : table of class per cds,te ID: CDS,TE,UNK ?RPT (simple repeat) 
            and modifiers: uniq,busco,duplicate, ..
  generates script run_gnodes.sh
  
=item run steps

  from dplx20gnodes.qsub.hist

  env genolist="genome/dropse20chrs.fa.gz genome/dplx20maca1refdeg.fa.gz genome/dcar20chr.fa.gz" datad=`pwd` qsub -q batch run_repmod2b.sh
  
  env asm=daphplx17evgt1m_cds.fa reads=readsf/SRR3090572_1.fastq.gz datad=`pwd` qsub -q debug run_gnodes_dnamap.sh
  env asm=dplx20maca1r_tefams.fa reads=readsf/SRR3090572_1.fastq.gz datad=`pwd` qsub -q debug run_gnodes_dnamap.sh
  env asm=daphplx_gasm16ml.fa reads=readsf/SRR3090572_1.fastq.gz datad=`pwd` qsub -q debug run_gnodes_dnamap.sh
  
  pt=daphplx17evgt1m_cds; env opts="-savereadids -crclassf=dplx20cdste.idclass -minident=0.45"  bam=${pt}_SRR3090572_1_bwa.bam  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh
  pt=dplx20maca1r_tefams; env opts="-savereadids -crclassf=dplx20cdste.idclass -minident=0.45"  bam=${pt}_SRR3090572_1_bwa.bam  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh
  env merge_readids="daphplx17evgt1m_cds_SRR3090572_1_bwa.readids dplx20maca1r_tefams_SRR3090572_1_bwa.readids"  cdstab=dplx20cdste_SRR3090572_1_bwa.readids  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh
  
  env opts="-crclassf=dplx20cdste.idclass" cdstab=dplx20cdste_SRR3090572_1_bwa.readids  bam=daphplx_gasm16ml_SRR3090572_1_bwa.bam  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh
  
  env idclass=dplx20cdste.idclass cds=daphplx17evgt1m_cds.fa te=dplx20maca1r_tefams.fa chr=daphplx_gasm16ml.fa datad=`pwd`   qsub -q debug run_gnodes_ann.sh
  
  env asmid=daphplx_gasm16ml covtab=daphplx_gasm16ml_SRR3090572_1_bwa.cdschr7b.covtab anntab=daphplx_gasm16ml.anntab asmdata=daphplx20chrs.metad  datad=`pwd` qsub -q debug run_gnodes_sum.sh
  
  #--- 2nd read set, want to merge covtab,chrtab for new sum table

=cut

use FindBin;
use lib ("$FindBin::Bin/../"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
our $EVIGENES="$FindBin::Bin/..";  
use strict;
use Getopt::Long;  

$EVIGENES =~ s,/\w+/\.\.$,,;
# export EVIGENES=/Very/long/system/path/to/bio/evigene/scripts/genoasm/..

my $VERS="7d";
my $debug=$ENV{debug}||1;
my $ncpu=$ENV{NCPU}||16; my $maxmem=$ENV{maxmem}||"64gb"; # cluster opts
   my $walltime= $ENV{walltime} || "1:00:00";
my($doREPMASK,$doBUSCO)= (1,1); # default write scripts, turn off
my $buscodb= $ENV{buscodb} || ""; # "eukaryota"; # or not?
   
my $gnodes1_dnamap=""; # local proc .. "$EVIGENES/genoasm/gnodes_dnamap.pl";
my $gnodes1_sam2cov="\$EVIGENES/genoasm/gnodes_sam2covtab.pl";
my $gnodes2_ann="\$EVIGENES/genoasm/gnodes_annotate.pl";
my $gnodes3_covsum="\$EVIGENES/genoasm/gnodes_covsum.pl";
my $gnodes3_ucgcov="\$EVIGENES/genoasm/gnodes_sam2ucgcov.pl";
my $gnodes4_genescov="\$EVIGENES/genoasm/gnodes_genescov.pl"; #upd21apr20
# add maybe: cdna_bestorf.pl ..

my $intitle= $ENV{title}||$ENV{name}||""; 
my $asmid; #no default $intitle; # intitle for output, asmid for data
my ($chrasm,$teseq,$cdsseq,$crclassf,$sampledata,$anntable,$datad,$mergecov,$species)=("") x 19;
my @reads;

my $optok= GetOptions( 
  'chrasm|assembly|genome=s',\$chrasm,
  'teseq=s',\$teseq,
  'cdsseq=s',\$cdsseq,  
  'reads=s',\@reads, # many?
  'outputdir|datad=s',\$datad,  
  'merge:s',\$mergecov,  
  'title|name=s',\$intitle,
  'asmid=s', \$asmid, # superceeds intitle for finding/using asm info
  'genomedata|sumdata|sampledata|metadata=s', \$sampledata, # optional input  << change to gnodes_pipe.cfg/config ?
  'anntable|annotations=s', \$anntable, # not an option?
  'idclass=s',\$crclassf, # table of chr/gene id => class ? want for BUSCO, TEfam, : No, rely on anntable having these
  'species=s', \$species, # repeatmasker opt only?
  'ncpu=i', \$ncpu, 'maxmem=s', \$maxmem, 'walltime=s', \$walltime,
  'REPMASK!', \$doREPMASK, # norepmask turns off
  'BUSCOdb=s', \$buscodb, 
  'debug!', \$debug, 
  );

## use metad/sampledata for opts: buscodb=  species=  repmaskdb==species=? and other params: title? ..
# check @ARGV for .fastq/.fq reads ?
if(@ARGV and my @rds= grep( /\.(fastq|fq)/, @ARGV)) { push @reads, @rds; }

unless($optok and (($chrasm and @reads) or $mergecov) ) {
  die "usage: $0 -chrasm mychr.fa -cds mygenecds.fa -teseq mytransposons.fa -reads SRRnnnn.fastq .. ";
}

sub basename { my($f,$suf)=@_; $f=`basename $f`; chomp($f); $f=~s/$suf// if($suf); return($f); }

unless($asmid){ $asmid= basename($chrasm,'\.fa.*'); } #? or not
# $asmid= $intitle if($intitle and not $asmid); # no, use chrasm name

my($nmeta,$metavals,$asmidlist)= readMetad($sampledata);
# metavals->{asmid}{keys}= vals
if($nmeta and not $species){ $species = $metavals->{$asmid}{species}||""; }
if($nmeta and not $intitle){ $intitle = $metavals->{$asmid}{asmname}||""; }
$intitle=$asmid unless($intitle);

# need datad path:
my $mydir=`pwd`; chomp($mydir); # user option? make new work dir?
if($datad) {
  unless($datad =~ m,^/, or $mydir =~ m,$datad,) { $datad="$mydir/$datad"; }
} else {
  $datad=$mydir;
}  

# add new main method: -merge_covtabs -out daphpulex_pa42v2_merge3r daphpulex_pa42v2_SRR*bwa.cc7c.covtab
# also if @reads > 1, do both: gen script from gnodes_pipe_MAIN for each @reads, then gnodes_mergecov_MAIN

if(@reads > 1) {
  #UPD21apr20: change sh_dnamap(..,@reads) to map all read sets, samtools cat -o all.bam parts.bam,
  #  then use all.bam for rest of calcs
if(1) {
  my ($reads0)= shift @reads;
  my($runsh)= gnodes_pipe_MAIN( $reads0, 0, @reads );
  warn "$runsh\n";

} else {  
  my @covtabs;
  for my $i (0..$#reads) {
    my($partsh,$chrcovtab,$sumout)= gnodes_pipe_MAIN( $reads[$i], 1+$i );
    push @covtabs, $chrcovtab;
    warn "$partsh\n";
  }
  my($runsh)= gnodes_mergecov_MAIN( @covtabs);
  warn "$runsh\n";
}
  
} elsif($mergecov) {
  my @covtabs= grep(/covtab$/, @ARGV);
  my($runsh)= gnodes_mergecov_MAIN( @covtabs );
  warn "$runsh\n";
  
} else {
  my($runsh)= gnodes_pipe_MAIN( $reads[0] );
  warn "$runsh\n";
}


#=======================================================


sub gnodes_pipe_MAIN {  
  my($reads,$it,@morereads)=@_;
  # my $reads= $reads[0]; # fixme
  
  return "# NO reads to map" unless($reads);
  if($it){ $it = "_r".$it; } else { $it=""; }
  my $runsh="run_gnodes_$intitle$it.sh"; # need read $i num
  
  my $OUTSH; open($OUTSH,'>',$runsh);
  
  sh_setup($OUTSH,$runsh);

  #? insert here call to sh_repeatmask() sh_buscoscan ?
  #  want results before sh_annot(), sh_dnamap(te) are run
  my($rmout,$buscout)=("",""); # options
  my $teann=$teseq;
  if($doREPMASK and not $teseq) {  
    my $rmaskdb= $metavals->{$asmid}{rmaskdb} || $species;
    ($rmout)= sh_repeatmask($OUTSH,$chrasm,$rmaskdb);
    $teann=$rmout; 
    # NOT $teseq= $rmout; # okay? << BUG at sh_dnamap($OUTSH,$teseq,$reads)
    ##  asmidx="bwax/ficari1asm_bwx"
    ##  bwa-mem2 index -p $asmidx ./rmout/ficari1asm.fa.out <<<
  }
  
  if($doBUSCO and $cdsseq) {  
    $buscodb ||= $metavals->{$asmid}{buscodb} || "";
    ($buscout)= sh_buscoscan($OUTSH,$cdsseq,$buscodb);
  }
  
  # this can be run independently of dnamap/sam2cov .. call before those, generates $crclassf
  my($anout,$idclassa)= sh_annot($OUTSH,$chrasm,$cdsseq,$teann,$crclassf,$buscout);
  $anntable= $anout; # is anntable option or not?
  unless($crclassf){ $crclassf= $idclassa; }
    # ? pre-make crclassf if doesnt exist and have cdsseq,teseq .. see also gnodes_annotate.pl make_idclass()
  
  #UPD21apr20: sh_dnamap_many : samtools cat -o all.bam @parts.bam, or for 1 reads does orig way
  #UPD21apr20: FIXME this causes readids _1/_2 to be same, spurious duplicates, one bug is that
  #   cdsasm count now ~2x of before, and 2x of CDSann when should be nearly same. Other bugs not yet found.
  #   solution? add .1/.2 suffix to read ids? maybe ok but ID number is used as integer. Change IDprefix?
  #  add sam flag for pair/mate?  that means more recoding w/ sam parsers
  my($gnbam,$tebam,$crbam);
  if(@morereads){
    ($gnbam)= ($cdsseq)? sh_dnamap_many($OUTSH,$cdsseq,$reads,@morereads) : 0;
    ($tebam)= ($teseq )? sh_dnamap_many($OUTSH,$teseq,$reads,@morereads) : 0;
    ($crbam)= ($chrasm)? sh_dnamap_many($OUTSH,$chrasm,$reads,@morereads) : 0;
  } else {
    ($gnbam)= ($cdsseq)? sh_dnamap($OUTSH,$cdsseq,$reads) : 0;
    ($tebam)= ($teseq )? sh_dnamap($OUTSH,$teseq,$reads) : 0;
    ($crbam)= ($chrasm)? sh_dnamap($OUTSH,$chrasm,$reads) : 0;
  }
  
  # my $ucgcovtab= $gnbam; $ucgcovtab=~ s/\.bam//; $ucgcovtab.="_ucg.covtab"; # in sh_ucgcov()
  my($ucgcovtab)= ($gnbam) ? sh_ucgcov($OUTSH, $gnbam, "",  $crclassf, "norecalc") : ("");
  #^ maybe insert before temap,chrmap
  
  my($cdstab,$cdsreadids)= ($gnbam) ? sh_sam2covtab($OUTSH, 'cds', $gnbam, $crclassf) : (0,0);
  my($tetab,$tereadids)  = ($tebam) ? sh_sam2covtab($OUTSH, 'te', $tebam, $crclassf) : (0,0);

  # insert sh_covsum($cdstab) .. gets KU,C value independent of chrasm  
  my($sumcds)=(0);
  if($cdstab) {
    my $cdsid= basename($cdsseq,'\.(fa|cds).*'); 
    # fixme: crclassf is output of sh_annot(): ficari18tsa1cds.fa > ficari18tsa1cds.idclass
    #here? ($ucgcov)= sh_ucgcov($OUTSH, $cdstab, $cdsid, $cdsid, "", $crclassf, "norecalc");  
    ($sumcds)= sh_covsum($OUTSH, $cdstab, $cdsid, $cdsid, "", $crclassf, "norecalc");  # no cds.anntab but idclass
  }
  
  my $mctids= $cdstab; $mctids=~s/.covtab//; $mctids=~s/cds//; $mctids.="_cdste.readids";
  my($cdsteids) = (not $tereadids) ? $cdsreadids 
      : (not $cdsreadids) ? $tereadids
      : sh_merge_readids($OUTSH,$mctids,$cdsreadids,$tereadids);
  
  my($chrcovtab)  = sh_sam2covtab($OUTSH, 'chr', $crbam, $crclassf, $cdsteids);

  #o: my($sumout)= sh_covsum($OUTSH, $chrcovtab); #global opts: ( $asmid, $intitle, $anntable, $sampledata )
  my($sumout)= sh_covsum($OUTSH, $chrcovtab, $asmid, $intitle, $anntable, $crclassf); 

  # add sh_genescov()
  my($gsumout)= sh_genescov($OUTSH, $cdstab, $chrcovtab, $asmid, $intitle, $anntable, $crclassf); 

  close($OUTSH);  
  system("chmod +x $runsh");
  return($runsh,$chrcovtab,$sumout,$gsumout); # ,$sumcds
}

sub gnodes_mergecov_MAIN {
  my (@covtabs)= @_;
  # my @covtabs= grep /covtab$/, @ARGV; #?? or -merge a.covtab,b.covtab,...
  
  # gnodes_pipe.pl -merge -title daphpulex_pa42r2 -chr daphpulex_pa42v2.fa  -sumdata daphplx20chrs.metad -ncpu 24  -maxmem 128gb
 
  my $nc=@covtabs;
  return "# NO covtabs to merge" unless($nc>1);
  
  my $mergetab= ($mergecov||$intitle)."_merge$nc.$VERS.covtab";  
  unless($anntable) {
    my $aname= ($chrasm) ? basename($chrasm,'\.fa.*') : $asmid;
    $anntable="$aname.anntab";
  }
  
  my $runsh="run_gnodes_merge_$intitle.sh";
  my $OUTSH; open($OUTSH,'>',$runsh);
  
  sh_setup($OUTSH,$runsh);
  my($mcovtab,$mchrtab)= sh_merge_covtabs( $OUTSH, $mergetab, @covtabs);
  #o: my($sumout)= sh_covsum($OUTSH, $mcovtab); #global opts: ( $asmid, $intitle, $anntable, $sampledata )
  my($sumout)= sh_covsum($OUTSH, $mcovtab, $asmid, $intitle, $anntable, $crclassf); 

  # add sh_genescov()
  my $mcdscov=""; # FIXME, merge cds.covtabs as well as chr.covtabs
  my($gsumout)= sh_genescov($OUTSH, $mcdscov, $mcovtab, $asmid, $intitle, $anntable, $crclassf); 

  close($OUTSH);  
  system("chmod +x $runsh");
  return($runsh);
}


# run_gnodes_dnamap.sh script ============================

sub sh_setup {
  my($OUTSH,$sname)=@_;

  # need datad path: above now, add mkdir datad
  # my $datad=`pwd`; chomp($datad); # user option? make new work dir?
  # FIXME maybe here? set shell VARs for primary data here, those use to test results of each step
  # but now sh_ subs for each method are setting output names

  # FIXME2: this part should be customizable from local data : {local}gnodes_setup.sh ??
  my $localset= $ENV{gnodes_setup} || "gnodes_setup.sh";
  my $topset="";
  if( -f $localset ) {
    open(LS,$localset); $topset= join "", <LS>; close(LS);
  } else {
    $topset=<<"EOS";
## --- gnodes_setup.sh ---    
#PBS -N gnodes_pipe
#PBS -l vmem=${maxmem},nodes=1:ppn=$ncpu,walltime=${walltime}
#PBS -V

## ensure have these.. trap missing?  
module load blast samtools 
module load repeatmasker busco
# module load bwa-mem # older
bwabin=\$HOME/bio/bwa/bin/; export PATH=\$bwabin:\$PATH
## --- end gnodes_setup.sh ---    
EOS
  }
  
  my $script=<<"EOS"; 
#! /bin/bash
## $sname
## create 'gnodes_setup.sh' to add local env, paths, cluster metacomments: PBS/SBATCH/XXX ...  
$topset

export EVIGENES=$EVIGENES
export PATH=\$PATH:$EVIGENES
export NCPU=$ncpu;

if [ ! -d $datad ]; then mkdir $datad; fi
cd $datad
mkdir parts

EOS

  print $OUTSH $script;
  return(1);
}

sub check_exec {
  my @exfail; # my %exfail;
  for my $e (@_) {
    #o: system("$e --help >/dev/null 2>&1") == 0 or fail "$e not found or failed to run";
    #o: print "$e OK\n";
    my $err= system("$e --help >/dev/null 2>&1");
    push @exfail, $e if($err);
  }
  return @exfail;
}

sub readMetad {  # from gnodes_covsum.pl
  my($sampledata)= @_;
  my($nsam,$aid)=(0,0);
  my %samvals=(); my @aid=(); 
  unless($sampledata) { # want empty %samvals
    return($nsam, \%samvals, \@aid);
  }
  warn "#read sampledata=$sampledata\n" if($debug);
  open(F,$sampledata); 
  while(<F>){
    next if(/^\W/);
    my($key,$val)= (m/^(\w+)\s*=\s*(.+)$/) ? ($1,$2):(0,0);
    next unless($key);
    
    if($key eq 'pt' or $key eq 'asmid') { $aid=$val; $nsam++;
      map{ $samvals{$aid}{$_}= 0; } qw(flowcyto cytomb asmtotal atotalmb asmname );
      push @aid, $aid;
      
    } elsif($key eq 'flowcyto') { # /^flowcyto=(.+)$/)  
      my $flowcyto= $val; # flowcyto=234-391 Mb
      $samvals{$aid}{flowcyto}= $flowcyto;
      $samvals{$aid}{cytomb}  = ($flowcyto=~m/(\d+)/)?$1:0; # can be range 160-180
  
    } elsif($key eq 'glncformula') { # /^glncformula=(.+)$/ 
      $samvals{$aid}{glncformula}= $val; # glncformula=nnn Mb
      
    } elsif($key eq 'asmtotal') { # /^asmtotal=(.+)$/)   
      my $asmtotal= $val; 
      $samvals{$aid}{asmtotal} = $asmtotal;
      $samvals{$aid}{atotalmb} = ($asmtotal=~m/(\d+)/)?$1:0; # can be range 160-180
  
    # } elsif($key =~ /asmname|kulo|kuhi/) { # ($key eq 'asmname' $key eq 'kulo' or $key eq 'kuhi' )
    #   $val=~s/\s.*//; $samvals{$aid}{$key}= $val;  # /^kulo=(\S+)/)  # test hi<lo
    # } elsif($key =~ /busco|rmask|species/) {
    #   $val=~s/\s.*//; $samvals{$aid}{$key} = $val;
    } else { # save all valid key = val
      $val=~s/\s.*//; $samvals{$aid}{$key} = $val;
    }
    
  } close(F);
  push @aid, 0 unless(@aid);
  return($nsam, \%samvals, \@aid); 
}


sub sh_dnamap {
  my($OUTSH,$asm,$reads)=@_;
 
  return("") unless($asm and $reads);
  my $btag="_bwa";
  my $aname= basename($asm,'\.(fa|fasta)'); # bad for teseq= chrasm.fa.out !*** '\.fa\w+ .gz
  my $nab=   basename($reads,'\.(fastq|fq)'); 
    map{ s/\.gz// } ($aname,$nab);
  my $name= $aname."_".$nab;
  my $outbam = "$name$btag.bam";
  my $bwalog = "$name$btag.log"; 
  
  my $script=<<"EOS"; 
# run_gnodes1_dnamap
if [ -s $outbam ]; then 
  echo reusing $outbam;
else
  echo START_dnamap  `date`
  asmidx="bwax/${aname}_bwx"
  if [ ! -f \$asmidx.ann ]; then
    mkdir bwax;
    ( bwa-mem2 index -p \$asmidx $asm ) > $bwalog 2>&1
  fi
  
  ( bwa-mem2 mem -a -t \$NCPU \$asmidx $reads | samtools view --threads \$NCPU -Sb -o $outbam ) > $bwalog 2>&1
  echo DONE_dnamap `date`
fi

EOS
  
  print $OUTSH $script;
  return($outbam);
}

sub sh_dnamap_many {  #UPD21apr20
  my($OUTSH,$asm,$reads,@breads)=@_;

  return("") unless($asm and $reads);
  return sh_dnamap($OUTSH,$asm,$reads) unless(@breads);
  
  my $btag="_bwa";
  my $aname= basename($asm,'\.(fa|fasta)'); # bad for teseq= chrasm.fa.out !*** '\.fa\w+ .gz
  my $nab=   basename($reads,'\.(fastq|fq)'); 
    map{ s/\.gz// } ($aname,$nab);
  my $name= $aname."_".$nab;
  # FIXME HERE: outbam > allbam to test reusing
  my $outbam = "$name$btag.bam";
  my $bwalog = "$name$btag.log"; 

  my $nbreads= 1 + @breads;
  my $allbam= $aname."_".$nab; $allbam=~s/_[12]$//;
  $allbam .= "_b$nbreads$btag.bam";
  my $outbam1= $outbam; $outbam=$allbam;
  
  my $script=<<"EOS"; 
# run_gnodes1_dnamap
if [ -s $outbam ]; then 
  echo reusing $outbam;
else
  echo START_dnamap  `date`
  asmidx="bwax/${aname}_bwx"
  if [ ! -f \$asmidx.ann ]; then
    mkdir bwax;
    ( bwa-mem2 index -p \$asmidx $asm ) > $bwalog 2>&1
  fi
  
  ( bwa-mem2 mem -a -t \$NCPU \$asmidx $reads | samtools view --threads \$NCPU -Sb -o $outbam1 ) > $bwalog 2>&1
EOS
  
  my @partbams= ($outbam1);
  for my $readb (@breads) {
    my $nab=   basename($readb,'\.(fastq|fq)'); map{ s/\.gz// } ($nab);
    my $name= $aname."_".$nab;
    my $bamb = "$name$btag.bam"; 
    my $logb = "$name$btag.log"; 
    push @partbams, $bamb;
    $script .= 
  "  ( bwa-mem2 mem -a -t \$NCPU \$asmidx $readb | samtools view --threads \$NCPU -Sb -o $bamb ) > $logb 2>&1\n";
  }
  $script .= "\n  samtools cat -o $allbam " . join(" ",@partbams) . "\n";;
  $script .= "  echo DONE_dnamap `date`\n";
  $script .= "fi\n\n";
 
  print $OUTSH $script;
  return($allbam);
}

# run_gnodes_sam2covtab.sh script ============================

# pt=daphplx17evgt1m_cds; 
# pt=dplx20maca1r_tefams; 
# env opts="-savereadids -crclassf=dplx20cdste.idclass -minident=0.45"  bam=${pt}_SRR3090572_1_bwa.bam  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh
# env merge_readids="daphplx17evgt1m_cds_SRR3090572_1_bwa.readids dplx20maca1r_tefams_SRR3090572_1_bwa.readids"  cdstab=dplx20cdste_SRR3090572_1_bwa.readids  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh
# env opts="-crclassf=dplx20cdste.idclass" cdstab=dplx20cdste_SRR3090572_1_bwa.readids  bam=daphplx_gasm16ml_SRR3090572_1_bwa.bam  datad=`pwd` qsub -q debug run_gnodes_sam2covtab.sh

sub sh_sam2covtab {
  my($OUTSH,$dtype,$ibam,$idclass,$readids)=@_;

  return("") unless($ibam);

  my $name= basename($ibam,'\.bam'); 
  my $outtab="$name.cc$VERS.covtab";
  my $outchrtab="$name.cc$VERS.chrtab";
  my $opts="-nodebug "; # stop tons of debug msg from iterates
  $opts.=" -crclassf=$idclass" if($idclass);
  if($dtype =~ m/cds|te/){ 
    $readids="$name.readids"; $opts.=" -minident=0.45 -savereadids=$readids"; 
    }  
  elsif($dtype =~ /chr/) {
    $opts.=" -readidtab=$readids " if($readids); 
  }
  # trap debug out to sam2cov log?   my $errlog = "$name.cc$VERS.errlog"

  my $script= <<"EOS";
# run_gnodes1_sam2cov
if [ -s $outtab ]; then
  echo reusing $outtab $outchrtab;
else
  echo START_sam2covtab  `date`
  i=0; while [ \$i -lt \$NCPU ]; do {
    $gnodes1_sam2cov -icpu \$i -ncpu \$NCPU $opts -out $outtab -bam $ibam  &
    i=\$((\$i + 1));
  } done
  
  wait
 
  $gnodes1_sam2cov  -merge -out $outtab  -bam $ibam 
  mv $outtab.pt* $outchrtab.pt* $readids.pt* parts/
  echo DONE_sam2covtab  `date`
fi
EOS
  
  print $OUTSH $script;
  return($outtab,$readids);
}

sub sh_merge_readids {
  my($OUTSH,$mergedids,@readid_parts)=@_;

  my $script=<<"EOS";
# run_gnodes1_merge_rids
if [ -s $mergedids ]; then
  echo reusing $mergedids ;
else
  echo START_gnodes1_merge_rids  `date`
  $gnodes1_sam2cov -merge readids -savereadids=$mergedids @readid_parts
  echo DONE_gnodes1_merge_rids  `date`
fi
EOS
  
  print $OUTSH $script;
  return($mergedids);
}

sub sh_merge_covtabs {
  my($OUTSH,$mergedcov,@covtab_parts)=@_;
  my $script="# run_gnodes1_merge_readids\n";
  
  my ($mergedchr,@chrparts)= map{ my $p=$_; $p=~s/.covtab//; $p.=".chrtab"; $p; } ($mergedcov,@covtab_parts);
  $script .= "$gnodes1_sam2cov -merge covtab -output=$mergedcov @covtab_parts\n";
  $script .= "$gnodes1_sam2cov -merge chrtab -output=$mergedchr @chrparts\n";
  
  print $OUTSH $script;
  return($mergedcov,$mergedchr);
}

#===== run_gnodes_annot.sh script =====================

sub sh_annot {
  my($OUTSH,$chr,$cds,$te,$idclass,$buscotsv)=@_;
 
  return("") unless($chr and ($cds or $te));
  my $annout= basename($chr,'\.\w+$') . ".anntab"; 
  
  my $opts="";
  $opts.= ($debug)?" -debug":" -nodebug";
  $opts.=" -cds $cds" if($cds);
  $opts.=" -te $te" if($te);
  $opts.=" -idclass $idclass" if($idclass);
  $opts.=" -busco $buscotsv" if($buscotsv);

  my $script=<<"EOS";
if [ -s $annout ]; then
  echo reusing $annout ;
else
  echo START_gnodes2_annotate  `date`
  $gnodes2_ann $opts -ncpu \$NCPU -chr $chr -output $annout
  echo DONE_gnodes2_annotate  `date`
fi

EOS

  print $OUTSH $script;
  # fixme: crclassf is output of sh_annot(): ficari18tsa1cds.fa > ficari18tsa1cds.idclass
  unless($idclass or not $cds){  $idclass=$cds; $idclass=~s/\.\w+$//; $idclass.=".idclass"; }
  return($annout,$idclass);  # idclass is output also, if not input
}

=item rmask.cat to te.seq

466 27.40 0.48 0.00 NC_037638.1 4243 4450 (27749750) Gypsy-24_DEl-I#LTR/Gypsy 2601 2809 (3339) m_b1s001i1

262 27.59 0.00 0.00 NC_037638.1 20221 20275 (27733925) C MuDR-1_LHu#DNA/MULE-NOF (635) 3923 3869 m_b1s001i6
   (27733925) C MuDR-1_LHu#DNA/MULE-NOF 
              ^ dingbat  extra colm for strand complement
              
 env idp=apimel perl -ne 'if(/^\d+/){ putfa() if($id and $fa);  @v=split;
 $cor=0; if($v[8] eq "C") { ($cor)= splice(@v,8,1); }
 ($cr,$cb,$ce,$te)=@v[4,5,6,8]; ($td,$tc)=split"#",$te;
 $cl=($tc=~/Unknown/)?"UNK":($tc=~/repeat|Sattelite/)?"repeat"
   :($tc=~/LTR|LINE|SINE|DNA|RC/)?"TE":"unk"; $cr=~s/\.\d//;
 $id=join"_",$cr,$cb,$ce;  $fa=""; $infa=1; } 
 elsif($infa){ ($d)=@v=split; if($d=~/^$cr/) { $s=$v[2]; $s=~s/\W//g; $fa.=$s; } } 
 END{ putfa(); } BEGIN{ $IDP=$ENV{idp}||"nada"; $IDP.="_rm"; }
 sub putfa{ $w=length($fa); $fa=~s/(.{100})/$1\n/g;
 print ">$IDP$id type=$cl; len=$w; tecl=$tc; teid=$td;\n",$fa,"\n"; } ' \
  apismel19hav31asm.fa.align > apismel19tefam.fa

=item parse_repmaskout : see gnodes_annotate.pl

=cut
 
sub sh_repeatmask {
  my($OUTSH,$chrasm,$species)=@_;
 
  return("") unless($chrasm); # guess at species for repmask? default?
  my $rmoutdir="./rmout"; my $cname=basename($chrasm);
  my $annout= $rmoutdir ."/$cname.out"; 
  my $rmlog = $rmoutdir ."/$cname.rmlog"; 

  # note rmkr appends to full chrasm name
  # want teseq.fasta output also, and te.gff, from rm -a align out?
  # rmout/
  
  # my $teseq="sss";
  # my $tegff="ggg";
  #.. rm.outs .. dont need fa.masked; use fa.cat to get align TE.seq ?
  # daphcari20chr.fa.cat.gz  daphcari20chr.fa.masked  daphcari20chr.fa.out	daphcari20chr.farm.gff	daphcari20chr.fa.tbl
  # daphcari20chr.farm.gff << bad name
  
  my $opts=" -q -xsmall";
  $opts.=" -species $species" if($species);

  my $script=<<"EOS";
if [ -s $annout ]; then
  echo reusing $annout ;
else
  echo START_repeatmasker  `date`
  mkdir $rmoutdir;
  RepeatMasker $opts -pa \$NCPU -dir $rmoutdir $chrasm > $rmlog 2>&1
  echo DONE_repeatmasker  `date`
fi

EOS

  print $OUTSH $script;
  return($annout);
}

=item buscoscan

# add sh_busco() ?? ordb=ordb_lineage
# run_BUSCO.py  -i $prot -l $ordb -o $tag$oname -m prot  --cpu $ncpu

try1:
  ordb_lineage = arthropoda from metad, not good enough for buscan
  
START_buscoscan Fri Feb 12 21:41:56 EST 2021
WARNING The dataset you provided does not contain the file dataset.cfg, likely because it is an old version. Default species (fly, eukaryota) will be used as augustus species
ERROR   Impossible to read arthropoda/   <<
DONE_buscoscan Fri Feb 12 21:42:31 EST 2021
START_gnodes2_annotate Fri Feb 12 21:42:31 EST 2021

=cut

sub sh_buscoscan {
  my($OUTSH,$cdsseq,$ordb_lineage)=@_;
 
  return("") unless($cdsseq); # guess at species for repmask? default?

  (my $aaseq=$cdsseq) =~ s/\.\w+$//; $aaseq.=".aa"; #  .cds chop bad: daphsim17evgt1cds.fa > daphsim17evgt.fa
  # unless(-f $aaseq) { }
  ## run_BU.py adds "run_" prefix
  my $buname= basename( $aaseq, ".aa");
  my $annotr= "run_" . $buname . "/full_table_" . $buname . ".tsv"; 
  my $buhmmer= "run_" . $buname . "/hmmer_output"; 
  #^^ fixme want this copied above run_xxx/ subdir, as? busco_full_table_$buname.tsv 
  # ficari18tsa1cds_buscofull_table.tsv better?
  # my $annout= "busco_full_table_".$buname.".tsv";
  my $annout=  $buname."_buscofull_table.tsv";
  
  ## dont set, busco is *supposed* to pick default eukaryota w/o -l option
  #nogo# $opts.=" -l $ordb_lineage" if($ordb_lineage);
  #this fixup may only be for IU HPC install w/ poor config file, lineage_path is correct but ;commented out
  my $optsh="";
  if($ordb_lineage) {
    if($ordb_lineage =~ m,/,) { 
      $optsh="opts=-l $ordb_lineage"; 
    } else {
      $optsh=<<"EOS";
if [ X != X\$BUSCO_CONFIG_FILE ]; then
  ORDB=`grep lineage_path \$BUSCO_CONFIG_FILE | sed 's/^.*= *//;'`;
fi
if [ X != X\$ORDB ]; then opts="-l \$ORDB/$ordb_lineage"; fi
EOS
      
    }
  }

  my $script=<<"EOS";
if [ -s $annout ]; then
  echo reusing $annout ;
else
  echo START_buscoscan  `date`
  if [ ! -f $aaseq ]; then \$EVIGENES/cdna_bestorf.pl -cdna $cdsseq -outaa $aaseq ; fi
  opts=""
  $optsh
  
  run_BUSCO.py \$opts --cpu \$NCPU -m prot -i $aaseq -o $buname 
  if [ -s $annotr ]; then cp -p $annotr $annout; rm -r $buhmmer; fi
  echo DONE_buscoscan  `date`
fi

EOS

  print $OUTSH $script;
  return($annout);
}

#===== run_gnodes_ucgcov.sh script =====================


sub sh_ucgcov {
  my($OUTSH, $inbams, $covtab, $lidclassf, $norecalc)=@_;  
 
  # note: inbams can be list "a.bam b.bam c.bam .."
  if($inbams =~ m/,/){ $inbams=~s/,/ /g; }
  my($inbam1,@inbam2)= split" ",$inbams;  
  unless($covtab) {
    $covtab= $inbam1; $covtab=~s/\.bam//; $covtab=~s/.bwa//; 
    if(@inbam2){ $covtab .="n".(1+scalar(@inbam2)); }
    $covtab.="_ucg.covtab";
  }
  my $opts="";
  $opts.=" -idclassf $lidclassf" if($lidclassf);
  #x $opts.=" -asmid $lasmid" if($lasmid);
  #x $opts.=" -title $lintitle" if($lintitle);
  # $opts.=" -sumdata $sampledata" if($sampledata);
  $opts.= ($debug)?" -debug":" -nodebug";

  my $testsum= ($norecalc and $covtab)? "-a ! -s $covtab " : "";
  
  my $script=<<"EOS";
if [ -s $inbam1 $testsum]; then  
  echo START_ucgcov  `date`
  $gnodes3_ucgcov $opts -output $covtab -bam $inbams
  echo DONE_ucgcov  `date`
fi

EOS

  print $OUTSH $script;
  return($covtab);
}
#===== run_gnodes_sum.sh script =====================
# env asmid=daphplx_gasm16ml covtab=daphplx_gasm16ml_SRR3090572_1_bwa.cdschr7b.covtab anntab=daphplx_gasm16ml.anntab asmdata=daphplx20chrs.metad  datad=`pwd` qsub -q debug run_gnodes_sum.sh
# $runapp $opts -asmid $asmid -anntab $anntab -sumdata $asmdata  -title $title  $covtab

sub sh_covsum {
  my($OUTSH, $covtab, $lasmid, $lintitle, $lanntable, $lidclassf,$norecalc)=@_;  
  #?? want this or not? eg not for cds.covsum
  #NO: $lasmid||= $asmid; $lintitle||=$intitle; $lanntable||=$anntable; $lidclassf||=$idclassf; #global defs
  
  my $opts="";
  $opts.=" -asmid $lasmid" if($lasmid);
  $opts.=" -title $lintitle" if($lintitle);
  $opts.=" -anntab $lanntable" if($lanntable);
  $opts.=" -crclass $lidclassf" if($lidclassf);
  $opts.=" -sumdata $sampledata" if($sampledata);
  $opts.= ($debug)?" -debug":" -nodebug";

  my $testsum= ($norecalc and $lintitle)? "-a ! -s ${lintitle}_sum.txt " : "";
  
  my $script=<<"EOS";
if [ -s $covtab $testsum]; then  
  echo START_covsum  `date`
  $gnodes3_covsum $opts $covtab
  echo DONE_covsum  `date`
fi

EOS

  print $OUTSH $script;
  return(1);
}

sub sh_genescov {
  my($OUTSH, $cdscov, $chrcov, $lasmid, $lintitle, $lanntable, $lidclassf,$norecalc)=@_;  
  
  # NOTE cdscov or chrcov may be empty, esp cdscov
  # NOTE title change for genescov? mix of cds+chr names
  # DEF out: xxx_genesum.txt, xxx.genemeans xxx.genexcopy
  
  my $opts="";
  $opts.=" -asmid $lasmid" if($lasmid);
  $opts.=" -title $lintitle" if($lintitle);
  $opts.=" -anntab $lanntable" if($lanntable);
  $opts.=" -crclass $lidclassf" if($lidclassf);
  $opts.=" -sumdata $sampledata" if($sampledata);
  $opts.= ($debug)?" -debug":" -nodebug";

  my $testsum="";
  #? my $testsum= ($norecalc and $lintitle)? "-a ! -s ${lintitle}_genesum.txt " : "";
  
  $opts.= " -cdscov $cdscov" if($cdscov);
  $opts.= " -chrcov $chrcov" if($chrcov); # assume this one is set
  
  my $script=<<"EOS";
if [ -s $chrcov $testsum]; then  
  echo START_genescov  `date`
  $gnodes4_genescov $opts
  echo DONE_genescov  `date`
fi

EOS

  print $OUTSH $script;
  return(1);
}
