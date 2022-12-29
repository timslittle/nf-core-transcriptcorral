#!/usr/bin/perl
# gnodes_samaddrdid.pl cut from evigene/scripts/genoasm/gnodes_sam2covtab8e.pl

=item UPD21AUG20: readReadIds() now is BIG-MEM-PIG

 UPD21AUG20: readReadIds() now is BIG-MEM-PIG with eg. pig genome data, 5/10 icpu failing in 240 GB mem
  need to rewrite to move cds-te-readid info into chr.bam rows, could also slim down readid hash w/ int ids
  .. move this call AFTER read_idclass, readChrtab  (failing samtools view -H there w/ big-pig)

=cut

use strict;
use Getopt::Long; # replace ENV w/ -opts
use constant UPD21AUG => 1; # UPD21AUG20 big-pig memory problems w/ readIds

my $PROGNAME='gnodes_samaddrdid'; # debug checkmem()
my $VERS='9a'; # 8e=revise readid meth, try to reduce mem use, misc. ; 8a=UPD21JUN prior
my $debug=1; # test
my $DEBUG_RDFILT= $ENV{DEBUG_RDFILT}||0; #UPD21AUG: DEBUG mem overflows in readfilter with big-pig data set

my($RIDPREFIX,$GIDPREFIX,$SPLITBAM,$npart, $MERGE,$NCDS, 
   $inbam,$outbam,$icpu,$ncpu,$samcpu,$cdstab, $hasgeneids)= (0) x 19;
my $SKIP_SAT='AS|XM|XO|XG|MD|YT|[a-z][a-z]'; # KEEP: NM: other? .. keep evg/gnodes eC/eG tags
my @ADD_SAT= qw(eC eG); # public info here, eC:i:1 == gene-cds Class:1|2|4|.., eG:X:id1,id2 == geneid-list, see readFilter()
my ($DROP_SEQ,$DROP_QUAL)= (0,1);
  
my $optok= GetOptions( 
  'bam=s', \$inbam, 
  'outbam=s', \$outbam, 
  'cdstab|readidtab=s',\$cdstab, # read table of cds numeric read IDs, num.aligns
  'hasgeneids!', \$hasgeneids, # for readReadIds, expect last col = geneids
    #^^ FIXME: allow multiple -readidtab => @readids option    
  'ridprefix=s',\$RIDPREFIX, 'gidprefix=s',\$GIDPREFIX,
  'icpu=i', \$icpu, 'ncpu=i', \$ncpu,
  'SPLITBAM!', \$SPLITBAM,  'npartsplit=i', \$npart, 'samcpu=i', \$samcpu,
  'merge!', \$MERGE, 
  'debug!', \$debug, 
  );

my $optinok= ($inbam and -s $inbam and ($MERGE or $SPLITBAM or $cdstab));

die "usage: gnodes_samaddrdid -bam name.chr_reads.bam -readidtab name.cdste.readids -outbam name.chr_cdsids.bam
  requires 'samtools view inbam'
  opts: -ridprefix=SRR1234 : strip readid prefix (memsave)
  -split -npartsplit=32   : split inbam to nparts, using icpu/ncpu threads, and -samcpu=4 for samtools
  -ncpu $ncpu -icpu 0..3  : subsets for parallel process, then samtools cat -o merged.bam part*.cdsids.bam
"  unless($optok and $optinok);

my %redata=();   #UPD21JUN8e: slim mem: change readid hash to one only: redata{idp}{inum}  = cdsrdid{idp}[inum] + generdid{idp}{inum}
$debug=1 if($DEBUG_RDFILT);

sub checkmem { my($ckp)=@_; 
  if(0 and $debug) { ## OFF til get something working... ? /usr/bin/ps
  # FIXME: find better checkmem, this is messy, has all icpu=[0..n] listed, ..
  # my $smc="ps -e -o pid,pcpu,pmem,stime,etime,command  | grep '$PROGNAME' | grep -v grep | head"; 
  ## NOT SEEN now: grep icpu.$icpu is likely problem.. need procid/pid of this process
  my $smc="ps -e -o pid,pcpu,pmem,stime,etime,command  | grep '$PROGNAME' | grep 'icpu.$icpu' | grep -v grep | head"; 
  my $smo=`$smc`; chomp($smo);
  warn"#checkmem[$ckp]: $smo\n"; 
  }
}

sub MAIN_stub {}

  $ncpu ||=1;  $icpu=0 if($ncpu==1);
  
  my $iname=$inbam; $iname =~ s/\.(bam|covtab|chrtab|readids)//;
  $outbam= "$iname.rdid.bam" unless($outbam); #? $iname.$VERS.covtab

  if($SPLITBAM) {
    my @res;
    warn "#$PROGNAME -split -npart $npart -icpu=$icpu/$ncpu -bam $inbam \n" if($debug);
    if(1) {
      if($npart<2){ $npart= $ncpu*2; } #??

      # my $nc= int($npart/$ncpu); my $nm= $npart % $ncpu; # if not 0, fix: extra ipart somewhere
      # for npart=24, ncpu=8 : nc=3, @ipart=($icpu, $ncpu+$icpu, 2*$ncpu+$icpu)
      # this works: $ncpu=10; $npart=24;  for($icpu=0;$icpu<$ncpu;$icpu++){  
      # my @ipart=(); for(my $i=$icpu; $i<$npart; $i+= $ncpu) { push @ipart, $i; } 
      # print "icpu=$icpu,ncpu=$ncpu,np=$npart, parts=@ipart\n"; }
    
      my @ipart=(); for(my $i=$icpu; $i<$npart; $i+= $ncpu) { push @ipart, $i; }
      @res= splitBam($inbam,$iname,$npart,@ipart); # out = iname.pt$I.bam
    
    } else {
      @res= splitBam_TRY1($inbam,$iname, $npart||$ncpu); # out = iname.pt$I.bam
    }
    exit;
  }
  
  if($MERGE) { 
    my $cmd="samtools cat -o $outbam  $iname*.pt*.bam"; # FIXME is inbam == outbam for merge?
    warn "#MERGE $cmd\n" if($debug);
    ## UPD check existence of partbams
    my $ptbam=`ls $iname*.pt*.bam`;
    my @ptbam= grep /\.bam/, split(" ",$ptbam);
    my $err=-1;
    if(@ptbam < 2 or grep(/^$outbam$/,@ptbam)) {
      warn "#MERGE $cmd\n" unless($debug);
      die "ERR: missing partbam=@ptbam or contains output:$outbam\n";    
    } else {
      $err= system($cmd);
    }
    exit $err;
  }
  
  warn "#$PROGNAME part=$icpu/$ncpu \n" if($debug);
  if($ncpu>1 and $icpu>=$ncpu) { die "# err: icpu >= ncpu : -i $icpu, -ncpu $ncpu"; }
  
  my $outbampart= $outbam; 
  if($ncpu>1) { $outbampart =~ s/\.bam/.pt$icpu.bam/; }
  
  #UPD21AUG20: readReadIds() now is BIG-MEM-PIG with eg. pig genome data, 5/10 icpu failing in 240 GB mem
  # need to rewrite to move cds-te-readid info into chr.bam rows, could also slim down readid hash w/ int ids
  
  #?? merge readid classes here? or use orig method:
  #  gnodes_sam2cov -merge readids -savereadids=$mergedids @readid_parts
  #  .. this creates class columns: readid  nrd  CDS TE UNK geneids 
  checkmem("before readReadIds");
  ($NCDS)= readReadIds($cdstab);  #?? FIXME: @readids multiple inputs possible, ie CDS.rdid, TE.rdid, UNK.rdid
  checkmem("after readReadIds");

  # ** need samtools view -h to pass header to outbam
  my $stopt = "-h"; # ($SKIPR==0)? "" : ($SKIPR==1) ? " -F 0x40" : " -F 0x80"; 
  my $swopt = "";
  if($samcpu>1) { $stopt.= " --threads $samcpu";  $swopt.= " --threads $samcpu"; }

  warn "# samtools view $stopt $inbam | samaddrdid | samtools view -Sb -o $outbampart\n" if($debug);
  if(-s $outbampart){ die "#ERR: I wont overwrite existing $outbampart\n"; }
  
  open(IN,"samtools view $stopt $inbam |") or die "ERR: samtools view $stopt $inbam";
  open(OUT,"| samtools view $swopt -Sb -o $outbampart -") or die "ERR: samtools view $swopt -Sb -o $outbampart";
  
  my ($nout,$nadd)= readFilter(*IN,*OUT);  
  
  close(IN); close(OUT);
  
  warn "gnodes_samaddrdid nadd=$nadd, nout=$nout to $outbampart\n";
  checkmem("end readFilter");
  
#-------------------------

=item splitBam($inbam,$outname,$npart, @ipart) : into n parts

  split by i-readid % nparts == ipart
  .. one way: open nparts name.part$i.bam, read inbam, write to part$i where $i = $iread % $nparts;

  .. works but runs as 1-cpu task, ie does ~ 1G/part x 20 parts in 1 hr, of 165 GB bam => needs 8hr/1cpu 
     of Nr.total=421_205_718, 1hr 1-cpu run to split np=20 gets about 10% of total, n_readid=49_000_001
     
  ?? this may not speed up any, takes same time to read thru 1 large.bam for any num of threads
  .. does samtools view --threads $ncpu affect any? test -samcpu opt?
  ** this IS much faster, TRY2 icpu/ncpu tasks to split to nparts 
      npart=24 ncpu=8 samcpu=3 (total cpu=24 mem=100G)
       => 16 GB out of 168 GB total in 15 min, vs 16 GB in 1hr of 1-cpu run TRY1
       => 65 GB out of 168 GB total in 1hr = 40%
  gnodes_pipe.626836.log : 165M read of 421M total done in 1hr = 40%
  #DERD: n_readid=165000001, nout=637097896, nout.pt1=212207400, icpu=0
  slurmstepd: error: *** JOB 626836 ON c2 CANCELLED AT 2021-08-22T15:35:05 DUE TO TIME LIMIT ***

  -- TRY3 
      npart=24 samcpu=8 ncpu=12 but run only icpu=0..3 (1/3)  (total cpu=24 mem=100G)
        to see if those fini in 1hr
     => ~19 GB out at 15 min (est from 8/24 parts),   
      #DERD: n_readid=51000001, nout=112042843, nout.pt1=55849055, icpu=0
     
  .. try ncpu, using -npart opt, icpu will write only npart/ncpu subsets?
     -nparts 24 -ncpu 8 -icpu 0..7
        each icpu writes 3 parts.bam, 
          icpu=0 => pt00.bam, pt08.bam, pt16.bam
          icpu=6 => pt06.bam, pt14.bam, pt22.bam 
          icpu=7 => pt07.bam, pt15.bam, pt23.bam [last]
  
=cut

sub splitBam {  
  my($inbam,$outname,$npart, @ipart)=@_;
  
  unless(@ipart) { my $np0=$npart-1; @ipart=(0..$np0); }
  my %ipart= map{ $_ => 1 } @ipart; 
  
  my $sopt=($samcpu>1)? "--threads $samcpu" : "";
  my($inbamh,@obamf,@obamh);
  foreach my $p (@ipart) { my $ptf=sprintf "pt%02d",$p; $obamf[$p]= "$outname.$ptf.bam"; }
  foreach my $p (@ipart) { 
    my $opipe="| samtools view $sopt -Sb -o ".$obamf[$p]." -";
    my $ok= open(my $oh, $opipe);  $obamh[$p]= $oh;
    die "#ERR: openOut($opipe)\n" unless($ok);
    warn "#ok: openOut($opipe)\n" if($debug and $ok);
  }

  my $outbamh= $obamh[0];
  open($inbamh,"samtools view -H $inbam |") or die "samtools view $inbam"; # header only here, to all obamh
  while(<$inbamh>){
    foreach my $p (@ipart) { $outbamh= $obamh[$p]; print $outbamh $_; }
  } close($inbamh);
  $outbamh= $obamh[0];
  
  my($n_readid,$nout,$nout1,$lrid,$iid,$ireadpart,$skipthis)=(0) x 9;
  open($inbamh,"samtools view $sopt $inbam |") or die "samtools view $sopt $inbam"; # data only here,  split out
  while(<$inbamh>) {

    my($rdid)= split; # my @samx= split; my($rdid)= $samx[0];
    # my($rdid,$fl,$cr,$cb,$cx,$cig,  $crp,$cxp,$px, $rseq,$rqual, @sat)= @samx;
        
    if($rdid ne $lrid) {      
      $ireadpart= ($npart > 1) ? $iid % $npart : 0;
      $iid++; $lrid=$rdid;
      $n_readid++;
      
      if($ipart{$ireadpart}) { $skipthis=0; $outbamh= $obamh[$ireadpart]; }
      else { $skipthis=1; } #? $outbamh=undef; 
      
      if($DEBUG_RDFILT and not $skipthis) {
        warn "#DERD: n_readid=$n_readid, nout=$nout, nout.pt1=$nout1, icpu=$icpu\n"
          if($n_readid % 1_000_000 == 1); # for  Nr.total=421_205_718
      }
    }
    
    unless($skipthis) {
    print $outbamh $_; 
    $nout++;  $nout1++ if($ireadpart == $ipart[0]);
    }
  }
  
  foreach my $p (@ipart) { my $h= $obamh[$p]; close($h); }  
  warn "#splitBam: npart=$npart, n_readid=$n_readid, nsam=$nout, nsam.pt1=$nout1\n" if($debug);
  return($nout, $npart, \@obamf);  
}

sub splitBam_TRY1 {  
  my($inbam,$outname,$npart)=@_;
  
  my($inbamh,@obamf,@obamh);
  for(my $p=0; $p<$npart; $p++) { my $ptf=sprintf "pt%02d",$p; $obamf[$p]= "$outname.$ptf.bam"; }
  for(my $p=0; $p<$npart; $p++) { 
    my $opipe="| samtools view -Sb -o ".$obamf[$p]." -";
    my $ok= open(my $oh, $opipe);  $obamh[$p]= $oh;
    die "#ERR: openOut(samtools -o $obamf[$p])\n" unless($ok);
    warn "#ok: openOut(samtools -o $obamf[$p])\n" if($debug and $ok);
  }

  my $outbamh= $obamh[0];
  open($inbamh,"samtools view -H $inbam |") or die "samtools view $inbam"; # header only here, to all obamh
  while(<$inbamh>){
    for(my $p=0; $p<$npart; $p++) { $outbamh= $obamh[$p]; print $outbamh $_; }
  } close($inbamh);
  $outbamh= $obamh[0];
  
  my($n_readid,$nout,$nout1,$lrid,$iid,$ireadpart)=(0) x 9;
  open($inbamh,"samtools view $inbam |") or die "samtools view $inbam"; # data only here,  split out
  while(<$inbamh>) {

    my($rdid)= split; # my @samx= split; my($rdid)= $samx[0];
    # my($rdid,$fl,$cr,$cb,$cx,$cig,  $crp,$cxp,$px, $rseq,$rqual, @sat)= @samx;
        
    if($rdid ne $lrid) {      
      $ireadpart= ($npart > 1) ? $iid % $npart : 0;
      $iid++; $lrid=$rdid;
      $n_readid++;
      $outbamh= $obamh[$ireadpart]; 
      
      if($DEBUG_RDFILT) {
        warn "#DERD: n_readid=$n_readid, nout=$nout, nout.pt1=$nout1, icpu=$icpu\n"
          if($n_readid % 1_000_000 == 1); # for  Nr.total=421_205_718
      }
    }
    print $outbamh $_; 
    $nout++;  $nout1++ if($ireadpart == 1);
  }
  
  for(my $p=0; $p<$npart; $p++) { my $h= $obamh[$p]; close($h); }  
  warn "#splitBam: npart=$npart, n_readid=$n_readid, nsam=$nout, nsam.pt1=$nout1 in $obamf[1]\n" if($debug);
  return($nout, $npart, \@obamf);  
}
  
sub readReadIds { 
  my($cdstab)= @_;
  # ?? FIXME: @readids multiple inputs possible, ie CDS.rdid, TE.rdid, UNK.rdid
  $NCDS=0;
  
  my($ok,$inh);
  if($cdstab and ($ok,$inh)= openRead($cdstab) and $ok) { 
    my $lhasgeneids=0; my @classlabel;
    # my ($isnum,$lastidp,$err)=(-1,"",0,0); 
    while(<$inh>){ 
    
      if(/^\W/) {
        if(/^#ReadID/ and not @classlabel) {           
          @classlabel = split; 
          for my $i (1..$#classlabel) { if($classlabel[$i] eq 'GeneID') { $lhasgeneids= $i; last; } }
        } # use @classlabel for output?
        next;
      }

      my @crclass=split; 
      my ($rdid,$geneids)=(0,0);

      # FIXME: readReadIds hasgeneids quirk expects geneids to have [A-Za-z] prefix, digit cols are class cols **
      # change to global hasgeneids, local lhasgeneids col = last col
      # GIDPREFIX cut affects that

      if($lhasgeneids>0) {
        ($geneids)= splice(@crclass,$lhasgeneids,1);
        if($GIDPREFIX){ $geneids=~s/$GIDPREFIX//g; } #UPD21AUG: strip GIDPREFIX like RIDPREFIX to save mem
 
      } elsif($lhasgeneids == 0) { # look for col, probly last
        if($hasgeneids) { $lhasgeneids= $#crclass; } # last col
        else {
          for my $i (1..$#crclass) { 
            if($crclass[$i] =~ /^\d/) { next; } # crclass cols are numeric
            elsif($crclass[$i] =~ /^\w/) { $lhasgeneids=$i; last; } # presume geneids like AT123, EVm123, XM_123
          }
          if($lhasgeneids == 0){ $lhasgeneids= -1; } # stop looking
        }
        if($lhasgeneids>0) {
        ($geneids)= splice(@crclass,$lhasgeneids,1);
        if($GIDPREFIX){ $geneids=~s/$GIDPREFIX//g; } #UPD21AUG: strip GIDPREFIX like RIDPREFIX to save mem
        }
      } 
      
      ($rdid)= shift @crclass; # always 1st col
      if($RIDPREFIX){ $rdid=~s/$RIDPREFIX//; }
      
      my $nd=$crclass[0]||1; 
      my $crclass=$nd;
      if( @crclass>1 ) { 
        $crclass=0; 
        $crclass |= 1 if($crclass[1]>0);  
        $crclass |= 2 if($crclass[2]>0);
        $crclass |= 4 if($crclass[3]>0);
      }

      # my($idp,$idn)= readIdNum($rdid);
      # if($lastidp and $idp ne $lastidp){ $useRIDp++; } #? dont use soft switch on hash/array
      
      ## change to int-only read ids?? save mem, needs more work to convert to ints
      ## or would this: $redata{$rdid}=xxx be smaller than that: $redata{$idp}{$idn}=xxx ?
      $redata{$rdid}= "$crclass\t$geneids"; 
      
      $NCDS++; # $lastidp=$idp;
    } close($inh);
    warn "# ncds readids=$NCDS from $cdstab \n" if $debug; 
  } 
  return($NCDS);
}

sub readFilter {  
  my($inhand,$outhand)=@_;
    
  my($n_readid,$nout,$nadd,$lrid,$iid,$ireadpart,$iskip)=(0) x 9;
  my $dopart= ($ncpu > 1);
  while(<$inhand>) {
    if(/^\@/){ print $outhand $_; next; } # sam hdr; do we want to add RG readgroup id,etc ?? value?

    # my @samx= split; 
    my($rdid,$fl,$cr,$cb,$cx,$cig,
        $crp,$cxp,$px, $rseq,$rqual, @sat)= split; # @samx;
        
    if($rdid ne $lrid) {      
      $ireadpart= ($dopart) ? $iid % $ncpu : 0;
      $iid++; $lrid=$rdid;
 
      if($ireadpart == $icpu) { 
        $n_readid++; $iskip=0;
        
        if($DEBUG_RDFILT) {
          warn "#DERD: nadd=$nadd, nout=$nout, n_readid=$n_readid, icpu=$icpu\n"
            if($n_readid % 500_000 == 1); # for  Nr.total=421_205_718
        }
      } else {
        $iskip=1;
      }
    }
    next if($iskip); # unless($ireadpart == $icpu);  # test even for ncpu=1, icpu=1  
          
    # $rseq ='*' if($DROP_SEQ); #? never? used for length(rseq)
    # below now: $rqual='*' if($DROP_QUAL); # always? qual not used
    my @osat= grep{ not m/^($SKIP_SAT):/ } @sat;

    my $ridp=$rdid; if($RIDPREFIX){ $ridp=~s/$RIDPREFIX//; }
    #? if( $fl < 0x100 ) # only need eG/eC on 1st rd map 
    if( $fl < 0x100 and my $eceg= $redata{$ridp}) {
      my($ec,$eg)= split"\t",$eceg; # always have both ec,eg from readReadIds
      $nadd++; $eg ||= 0;
      @osat= ( "eC:i:$ec", "eG:Z:$eg", @osat);
      # unshift @osat, "eG:Z:$eg"; # if($eg =~ /\w/); # always have eG:Z:0 if have eC:
      # unshift @osat, "eC:i:$ec"; # if($ec =~ /\d/);
    }
    
    print $outhand join("\t",$rdid,$fl,$cr,$cb,$cx,$cig,
                $crp,$cxp,$px, $rseq, '*', @osat)."\n"; # $rqual always = '*'
    $nout++;
  }
  
  return($nout,$nadd);  
}

sub openRead{ my($fn)=@_; my $fh;
  $fn.=".gz" unless(-s $fn or $fn=~/.gz$/);
  my $ok= (! -s $fn)? 0 : ($fn=~/.gz$/) ? open($fh,"gunzip -c $fn |") : open($fh,$fn);
  warn "#openRead($fn) = $ok\n" if $debug;
  return($ok,$fh); 
}



=item try1

START_samaddrdid Sat Aug 21 16:22:28 EDT 2021
/geode2/home/u030/gilbertd/Carbonate/bio/evigene/scripts/genoasm/gnodes_samaddrdid.pl -debug -ridprefix=SRR4341337. -readidtab=pig18evg4wf_t1cds_SRR4341337_b2_bwa.readids -bam chrpig11c_SRR4341337_b2_bwa.bam
#openRead(pig18evg4wf_t1cds_SRR4341337_b2_bwa.readids) = 1
# ncds readids=86815902 from pig18evg4wf_t1cds_SRR4341337_b2_bwa.readids 
# samtools view -h chrpig11c_SRR4341337_b2_bwa.bam | samaddrdid | samtools view -Sb -o chrpig11c_SRR4341337_b2_bwa.rdid.bam
#DERD: nadd=0, nout=0, n_readid=1, icpu=0
#DERD: nadd=9887493, nout=12479923, n_readid=500001, icpu=0
#DERD: nadd=20042205, nout=25213278, n_readid=1000001, icpu=0
#DERD: nadd=29989763, nout=37742826, n_readid=1500001, icpu=0
#DERD: nadd=40009956, nout=50359041, n_readid=2000001, icpu=0
#DERD: nadd=50191296, nout=63118890, n_readid=2500001, icpu=0

samtools view -F 260  chrpig11c_SRR4341337_b2_bwa.rdid.bam | head

SRR4341337.12	0	chr1	10047128	0	100M	*	0	0	s	*	NM:i:0	XS:i:100
SRR4341337.17	16	chr15	32807967	0	65S33M2S	*	0	0	s	*	NM:i:0	XS:i:33
SRR4341337.20	0	chr13	104534046	0	47M53S	*	0	0	s	*	eC:i:1	eG:Z:Susscr4EVm009461t1	NM:i:1 XS:i:42
SRR4341337.21	16	chr17	50962877	0	46S31M23S	*	0	0	s	*	NM:i:0	XS:i:31
SRR4341337.22	0	chr1	93124355	0	100M	*	0	0	s	*	eC:i:1	eG:Z:Susscr4EVm011357t1,Susscr4EVm016404t1,Susscr4EVm058735t1	NM:i:1	XS:i:95
SRR4341337.26	0	chr6	43865301	0	25M2I73M	*	0	0	s	*	eC:i:1	eG:Z:Susscr4EVm050552t1 NM:i:10	XS:i:50
SRR4341337.28	0	chr1	79188066	0	5S45M50S	*	0	0	s	*	NM:i:2	XS:i:35
SRR4341337.31	16	chr8	106299999	0	100M	*	0	0	s	*	eC:i:1	
  eG:Z:Susscr4EVm001606t1,Susscr4EVm018729t1,Susscr4EVm025703t1,Susscr4EVm033230t1,Susscr4EVm035459t1,Susscr4EVm043407t1,Susscr4EVm044163t1,Susscr4EVm049977t1,Susscr4EVm086073t1	
  NM:i:3	XS:i:88

=item checkmem

  memuse: perl module needs special compile: use Devel::Peak; mstat("marker1");'
  perl -e 'use Devel::Peek; mstat("one");'
     one: perl not compiled with MYMALLOC
     
  syscall: ps -e -o pid,pcpu,pmem,stime,etime,command --sort=-pmem | grep $progname | head
  PID %CPU %MEM STIME     ELAPSED COMMAND
  4701  5.2  0.7 Aug08 13-02:07:37 /usr/lpp/mmfs/bin/mmfsd
  15521  0.0  0.2 Aug19  1-19:12:39 python
  ^^^^ FIXME: find better checkmem, this is messy, has all icpu=[0..n] listed, ..

=item chr.sam

  @bwatags=qw(NM MD AS XS);
  @mimtags=qw(NM AS ms nn tp cm s1 s2 de rl);
  my $SKIP_SAT='AS|XM|XO|XG|MD|YT|[a-z]\w'; # KEEP: NM: other?

chrpig11c_SRR4341337_1_bwa.bam
SRR4341337.12	0	chr1	10047128	0	100M	*	0	0	
  AGAAACAAAGCCTTAATAAATAGCCAAATGCCTCAAACGTATCTGTTAAATCACTTCGCACAATGAACGTGCTCTTACAAACATCACCCAACTAAACTAC	
  @@@DDDDD32CFFEFFEHEGIIGIIHIGICHEHGFFGHGCD>GIIGHIIIIIIIGDHGIGE6CDC)=@(6=2;BCE>BDC(6>C5>C??B##########	
  NM:i:0	MD:Z:100	AS:i:100	XS:i:100
SRR4341337.12	256	chr3	10047128	0	100M	*	0	0	*	*	NM:i:0	MD:Z:100	AS:i:100
SRR4341337.12	256	chr9	10047128	0	100M	*	0	0	*	*	NM:i:0	MD:Z:100	AS:i:100
SRR4341337.12	256	chrY	10047128	0	100M	*	0	0	*	*	NM:i:0	MD:Z:100	AS:i:100

pig18evg4wf_t1cds_SRR4341337_1_bwa.bam
SRR4341337.10	0	Susscr4EVm095049t1	262	5	96M4S	*	0	0	
  CACTGCATTCTATCAAGCTGCATAAAATCCACTGCATTCAATCGGGATGAATTTTACCAGGCTGCATGGAATCCACTGCATTCAATCATGCTTCATGGAA	
  @@@FFFFEHHHFFHGGIGIJIIJIIJJGBFGIJJIIIIIIFGIJJIIBHE>FHIJGDC)=@GGHIEGDG7).=))7=@DEFFED>>>6@ACCC#######	
  NM:i:10	MD:Z:10A4T14C7T3A3T1C5C0T3T36	AS:i:46	XS:i:41
SRR4341337.10	256	Susscr4EVm103864t1	141	0	29H71M	*	0	0	*	*	NM:i:6 MD:Z:6G13C9T24G1G5G7	AS:i:41
SRR4341337.10	256	Susscr4EVm034071t1	249	0	60H32M8H	*	0	0	*	*	NM:i:0	MD:Z:32	AS:i:32

pig18evg4wf_t1cds_SRR4341337_1_mim.bam
SRR4341337.10	0	Susscr4EVm095049t1	262	10	96M4S	*	0	0	  s q	
  NM:i:10	ms:i:92	AS:i:92	nn:i:0	tp:A:P	cm:i:2	s1:i:2 s2:i:34	de:f:0.1042	rl:i:0
SRR4341337.10	256	Susscr4EVm118911t1	155	0	49M51S	*	0	0	*	*	
  NM:i:5ms:i:48	AS:i:48	nn:i:0	tp:A:S	cm:i:2	s1:i:34	de:f:0.102	rl:i:0
SRR4341337.22	16	Susscr4EVm011357t1	455	1	100M	*	0	0	 s q	
  NM:i:4	ms:i:164	AS:i:164	nn:i:0	tp:A:Pcm:i:3	s1:i:50	s2:i:67	de:f:0.04	rl:i:0

=cut



