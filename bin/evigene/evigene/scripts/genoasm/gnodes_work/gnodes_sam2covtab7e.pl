#!/usr/bin/perl
# sam2covtab7b.pl aka evigene/scripts/genoasm/gnodes_sam2covtab.pl
# from make_cdsxchr_covtab3h.pl, samcopyn5tab.sh

=item UPD7b 2021jan

  env opts="-savereadids -minident=0.45"  vers=7n \
   bam=dmag7fincds_tefam20skm_SRR7825549cl_1a_bwa.bam datad=`pwd` qsub -q debug run_sam2covtab7iu.sh
  add opt -crclassf cds_te.idclass table for dmag7fincds_tefam20skm: classes = CDS, TE, UNK 

  env opts="-savereadids -crclassf=dmag7cdste.idclass -minident=0.45"  bam=dmag20sk4ma_tefam_SRR7825549cl_1a_bwa.bam  datad=`pwd` \
    qsub -q debug run_sam2covtab7b.sh
  env opts="-savereadids -crclassf=dmag7cdste.idclass -minident=0.45"  bam=dmag7fincds_SRR7825549cl_1a_bwa.bam  datad=`pwd` qsub -q debug run_sam2covtab7b.sh

  env bam=dmag15nwb2asm_SRR7825549cl_1a_bwa.bam vers=7b cdstab=dmag20cdste_SRR7825549cl_1a_bwa7b.readids datad=`pwd` qsub -q debug run_sam2covtab7b.sh
  env bam=dmag19skasm_SRR7825549cl_1a_bwa.bam vers=7b cdstab=dmag20cdste_SRR7825549cl_1a_bwa7b.readids  datad=`pwd` qsub -q debug run_sam2covtab7b.sh

  # FIXMEd savereadids parts needed, ncpu overwrites same list

=item  id class table CDS/TE/UNK 

  * add more id class info, new cols?, use with annotate cds/te.seq, gff
  grep -h '>' dmag20sk4ma_tefam.fa ../dmag7finloc9b_cds.fa |  perl -ne  \
  '($id)=m/>(\S+)/; $cl="UNK"; if(/type=CDS/){ $cl=(m/Class:Transposon|Name=TE:/)?"TE":"CDS"; } 
  elsif(/teclass=(\w+)/){ $tc=$1; $cl=($tc=~/DNA|RC|LINE|SINE|LTR/)?"TE":"UNK"; } print "$id\t$cl\n"; ' \
    > dmag7cdste.idclass

=item upd7a 2021jan

  a. retool for parallel use of samtools bigdata.bam
  -- samtools -m MINALIGN for cigar min align, cut out all low qual dup aligns + noaligns

  -- parallel use by simplest way of process by read ids cut by icpu % ncpu,
        ie icpu % ncpu = 0,1,2,3,..  then process only 0th, 1st, .. read ID of modulus ncpu,
        and write full table set for ith part
        then merge all parts, part tables are additive
  from sam2covtab6c.pl of 2020nov which didnt work in parallel
  
  b. expand cdsread id table usage to table of read ids x class (cl 1,2,3)
     classes: cds-sure, te-sure, ambiguous cds-te (agt? = ambig-gene-te)
     from reads mapped to gene-cds, transposon-seq, and ambig/undef-cds-te seq
     
=item UPD7e 21apr24

 UPD7e: test opt changes, -minident 0.65+/-  -MIN_DUPIDENT 0.98+/-
  a. add 2ndary aligns (common w/ bwa) 
  b. adjust dupident counts (up? down?) MIN_DUPIDENT << new default MIN_DUPIDENT=0.999
  c. test perfect/imperfect aligns? MIN_IDENT?

=item UPD7f long low-acc reads need special align cov handling
  
  putlongmap() replaces putmap2b() for -longreads
  
  UPD7f: putmap2b/readFilter
  .. fix for shitty low accuracy long reads, ie 7000 bp long w/ 15% errors,
  .. need to cut very long cigar align into BN bin sized bytes for coverage calcs
  UPD7f:  use length(sam.rdseq)  instead of lenc ?? 
  UPD7f:  also sum all readlen for cov calcs, NSAMPLE not good enough

  trigger on read-length? ie if rdlen > 250|350|...
  
=item MINALN  use

   MINALN is critical param, 
   by default calculated from sampling input sam as $readlen * $MIN_IDENT
   for parallel, efficient use, give as command option 
     -MINALN=65 for 100 bp reads for MIN_IDENT=0.65
     -MINALN=97 for 150bp, 162 for 250 bp, etc
   
=cut

use strict;
use Getopt::Long; # replace ENV w/ -opts

my $VERS='7f';
#  '7d'; #large updates 2020.dec-2021.feb; small upd 20nov30, 6c? 7b: upd21jan > read id x chr tag class table
use constant UPD7 => 1;        
use constant UPD7e => 1;   #upd21apr24: higher dupident default (0.98 => 0.999), count 2nd suppl aligns of bwa
use constant UPD7f => 1;   #upd21apr28: for long,lowqual reads; addCigar upd should go to short rd also

my $debug=1; # test
my $MIN_IDENT= 0.40; # UPD7e was 0.65; # lo is best?
my $MIN_DUPIDENT = 0.999; #UPD7e was 0.98; # was .99/1; hi is best? or 1.0; # == ident equal to 1st/top align, lower if desired
# ^ ?? change to offby=1,2,3 .. dup read align below best align
my $ALLOW_SOFTCLIP = 0; # default? on, -nosoft turns off
my $SAVE_RDIDS= 0;  
use constant kMAXIDCLASS => 9; # idclass limit, using 3-4 now
my $CRTPAT=''; # no default, this is bad: '^[A-Za-z0-9]+'; #UPD7? TE|CDS|UNK chr class tag for SAVE_RDIDS
my $BN=$ENV{binsize}||100; 
my $SKIPR=$ENV{skipread}||0; # FIXME: skipread=0 means ignore, readmap.bam not paired

my($icpu,$ncpu,$outtab,$inbam,$inchrtab,$duptab,$cdstab,$rdidtab,$topcount,$crclassf)=(0) x 19;
$rdidtab= undef;
use constant NSAMPLE => 20000; # 500; # 500 not enough when variable like Fig.MiSeq
my( $MINALN,$nsample, $readlen, $sum_readlen)= (0) x 9; my @sreadlen; # $sminaln, $sreadlen, 
my($n_readid, $n_nomap, $n_notcr, $n_mapbad, $n_dupbad, $n_mapok, $n_partb,
   $n_intron, $n_insert, $n_delete, $n_softclip, $n_mismatch)= (0) x 19; # globals?
my($topaln,$lid,$idp,$idn,$ndi,$ncdsrd)=(0) x 9;
my $MERGE = undef; # defined($MERGE) now means do it
my ($RDIDIsNum,$RDIDIsNumerr)=(0,0);  #UPD7e need to test each read set?
my $LONGR=0; # ($inrdlen>500); # UPD7f FIXME
my $MAXSHORTRD= 500; # UPD7f, switch to LONGR method if rdlen>this

our(%crmax,%cov,%covt,%covu,%allcov,%allcovt,%allcovu,@savemap,@saveid);
#UPD7: change %cov,covt,covu to single hash of cdste-readid-hit-by-class
# FIXME: 2 uses of chrtab: input for chr, crlen, reads,  differs from output : two opt names?
# .. below outchrpart should not eq chrtab input name

my $optok= GetOptions( 
  'output|covtab=s',\$outtab, 'bam=s', \$inbam, 
  'chrtab=s',\$inchrtab, #<< change to $inchrtab to avoid fname confusion, outchrtab differs, always $outtab.chrtab
  'cdstab|readidtab=s',\$cdstab, # read table of cds numeric read IDs, num.aligns
  'savereadids:s',\$rdidtab, # write table of cds numeric read IDs, num.aligns
  'CRTPAT=s',\$CRTPAT, # use w/ savereadids
  'crclassf|idclassf=s',\$crclassf, # alt table chr => class for savereadids
  'icpu=i', \$icpu, 'ncpu=i', \$ncpu,
  'topcount=i', \$topcount, # quick sample, samview bam | top -n 10_000_000
  'minident=s',\$MIN_IDENT, 'mindupident=s', \$MIN_DUPIDENT,  
  'minalign=s',\$MINALN,  
  'softclip!',\$ALLOW_SOFTCLIP,
  'longreads!',\$LONGR,
  'binsize=i',\$BN, 'skipread=i',\$SKIPR, 
  'merge:s',\$MERGE, # UPD7, flag to merge ncpu part tables : 'merge!',
  'debug!', \$debug, 
  );

my $optinok= (($inbam and -s $inbam) or (defined $MERGE and ($outtab or @ARGV)) );

die "usage: sam2covtab -bam name.chr_reads.bam -readidtab name.cdste.readids -out name.v$VERS.covtab
  requires 'samtools view inbam'
  opts: -minident=$MIN_IDENT -mindupident=$MIN_DUPIDENT -minalign=$MINALN 
  -ncpu $ncpu -icpu 0..3  : subsets for parallel process, then -merge subsets
  -binsize=$BN  -skipread=2|1 (which of paired reads)  
  -savereadids  -idclass genes.classtab : save read_ids from CDS/TE align, for later match, with idclass types
sam2covtab -merge -out name  : merge all table parts from icpu 0..3,  
sam2covtab -merge=chrtab -out name.chrtab part1.chrtab part2.chrtab  : merge one table parts from ARGS
"  unless($optok and $optinok);



## FIXME: readIdNum() for new ids=SRRxxx.nnn and .nnn_2 or .nnn/2
## .. RDIDIsNum = 3 for idpat=prefix.nnn[_/\D]p
## Also cant expect type 3 at first,  RDIDIsNum not stable..

sub readIdNum { # UPD8a
  my($rdid)=@_;
  my($idp,$idn,$idx)=(0,0,0);
  if($RDIDIsNum == 1) { $idn= $rdid; } 
  elsif($rdid =~ /^\d/) { $idn=$rdid; $RDIDIsNum=1; }
  elsif($rdid =~ m/^(.+)(\d+)\D(\d+)$/) {
    ($idp,$idn,$idx)=($1,$2,$3); $idp.="p$idx";  # dont set RDIDIsNum
  } elsif($rdid =~ m/^(.+)(\d+)$/) {
    ($idp,$idn)=($1,$2);   # dont set RDIDIsNum
  } else {
    $RDIDIsNumerr++; ($idp,$idn)= $rdid =~ m/(.+)(\d+)/?($1,$2):($rdid,1);   
  }
    
  $idn= int($idn); #ensure, eg. srr123.456.sra style
  return($idp,$idn);   
}

sub readIdNum_OLD { # UPD7e
  my($rdid)=@_;
  my($idp,$idn)=(0,0,0);
  if($RDIDIsNum == 1) { $idn= $rdid; }
  elsif($RDIDIsNum == -1) { ($idp,$idn)=split/\./,$rdid; }
  elsif($RDIDIsNum == 0) {  
    if($rdid =~ /^\d/) { $idn=$rdid; $RDIDIsNum=1; } 
    elsif($rdid =~ /\.\d/) { ($idp,$idn)=split/\./,$rdid; $RDIDIsNum=-1; }
    else { $RDIDIsNumerr++; ($idp,$idn)= $rdid =~ m/(.+)(\d+)/?($1,$2):($rdid,1);  } #  last if($err>9); next;
  }
  $idn= int($idn); #ensure, eg. srr123.456.sra style
  return($idp,$idn);    
}

sub MAIN_stub {}

  $ncpu||=1;  $icpu=0 if($ncpu==1);
  $MIN_IDENT /= 100 if($MIN_IDENT > 1);
  $MIN_DUPIDENT /= 100 if($MIN_DUPIDENT > 1);
  
  # (my $iname=$inbam) =~ s/\.bam//;
  my $iname=$inbam || $outtab; $iname =~ s/\.(bam|covtab|chrtab|readids)//;
  $outtab= "$iname.v$VERS.covtab" unless($outtab); #? $iname.$VERS.covtab
  my $outchrtab= $outtab;  $outchrtab=~s/\.\w+$//; $outchrtab.=".chrtab";
  #o: $inchrtab= "$iname.chrtab" unless($inchrtab);  # renamed inchr, now no default input, see outchrtab
  # no default for cdstab,input readidtab
  if(defined($rdidtab)) { 
    $SAVE_RDIDS= 1;  
    $rdidtab= "$iname.readids" unless($rdidtab =~ /\w\w/); 
  }
 
  my $outpart= $outtab;
  # outchrpart should not eq chrtab input name
  my $outchrpart= $outchrtab; if($outchrpart eq $inchrtab) { $outchrpart.="out"; } # name bug: .chrtabout with -merge
  
  # gnodes_sam2covtab.pl -merge chrtab -output=ficari20r4m_merge2.7d.chrtab => chrtabout , dont want..
  # have two uses: chrtab input for read processing, chrtab output name for merges; change chrtab.in param?
  # does rdidtab need icpu parts?

  # BUG now for -merge parts, outchrpart=XXX_bwa.cc7d.covtab.ptNN, merge looks for XXX_bwa.covtab.ptNN
  # >> need outchrtab as well as outchrpart, should match but for .ptNN
  # i: $EVIGENES/genoasm/gnodes_sam2covtab.pl -icpu $i -ncpu $NCPU -nodebug  -readidtab=arath18tair1cds_ERR4586299_1_bwa.readids  -out arath18tair_chr_ERR4586299_1_bwa.cc7d.covtab -bam arath18tair_chr_ERR4586299_1_bwa.bam 
  # m: $EVIGENES/genoasm/gnodes_sam2covtab.pl  -merge -out arath18tair_chr_ERR4586299_1_bwa.cc7d.covtab  -bam arath18tair_chr_ERR4586299_1_bwa.bam

  if(defined $MERGE) { 
    # -merge : merge all table type, or -merge=covtab,chrtab,readids,.. which type only
    $MERGE="covtab,chrtab,readids" unless($MERGE=~/\w/); # all types
  }

  if($MERGE) {
    ## add @ARGV list of input parts to merge? opt to merge single table type
    ## eg. -merge=readids -out allparts.readids part1.rdid part2.rdid
    ## NOTE: this mergeparts() replaces merge(@covtabs) of generdcov4tab.pl
    ## maybe insert special call: -merge=readidsOnly -out=merged.readids set1.readids set2.readids ...
    ## if($MERGE =~ /readids/){  mergeparts( "readidlist", $outtab, @ARGV, ); exit; }
    # my @mparts= @ARGV;
    # @mparts=`ls $ofile.pt*`; chomp(@mparts); # do this in sub
    
    if(1) { #upd usage
      mergeparts("covtab", $outtab, @ARGV) if($MERGE=~/cov/); 
      mergeparts("chrtab", $outchrtab, @ARGV) if($MERGE=~/chr/); # -merge assume chrtab=output name, was outchrpart => .chrtabout bug
      mergeparts("readids", $rdidtab || "$iname.readids", @ARGV) if($MERGE=~/readid/);
      
    } else { #orig    
      mergeparts("covtab", $outtab); 
      mergeparts("chrtab", $outchrpart); 
      mergeparts("readids", $rdidtab || "$iname.readids", );
    }
    
    exit;
  }
    
  warn "#sam2covtab par: minident=$MIN_IDENT, mindupid=$MIN_DUPIDENT, softclip=$ALLOW_SOFTCLIP, BIN=$BN, part=$icpu/$ncpu \n"
    if($debug);
  if($ncpu>1 and $icpu>=$ncpu) { die "# err: icpu >= ncpu : -i $icpu, -ncpu $ncpu"; }

  my($ok,$inh);
  use constant USE_RIDp => 1; #  problem if this is big memory pig
  my @cdsrdid=(); my %cdsrdid=(); # for useRIDp == $cdsrdid{$idp}[$idnum]
  my @cdsrdclass_in=(0) x 5; my @cdsrdclass_sum=(0) x 5; my @cdsrdclass_miss=(0) x 5; 
  my ($NCDS,$useRIDp)=(0,USE_RIDp);  
  
  #UPD7: revise here to read multi-col read id table: rdid, class1, class2, class3; classes are 0/1 flag
  if($cdstab and ($ok,$inh)= openRead($cdstab) and $ok) { 
    my ($isnum,$lastidp,$err)=(-1,"",0); 
    while(<$inh>){ 
      if(/^\W/){
        if(/^#ReadID/){ my($lrd, @classlabel)=split; } # use @classlabel for output?
        next;
      }
      ##o my($rdid,$nd)=split; 
      ##UPD7: change input to multicol for classes
      ## ($nd,@crclass) may be counts for 1,2,3 read classes, zero for no class, 1+ for yes
      ## for @crclass>0 use 0/1+ as flag to set class bit
      ## $crclass = bit0 | bit1 | bit2 ? ie 1/2/4 or 1+2,1+4,1+2+4 = 7
      ##-- FIXME for diff @crclass cols: nmap | nmap CDS | nmap CDS TE | nmap CDS TE UNK ..
      ## set z=1 if @crclass>1, crclass |= 1,2,4 if($crclass[1,2,3]>0); 
      #ReadID	nmap	CDS	TE
      # SRR7825548.10000011	1	1	0
      # SRR7825548.10000013	1	0	1
      # SRR7825548.10000014	1	1	0
      # SRR7825548.10000024	1	1	0
      # SRR7825548.10000038	2	2	0

      my($rdid,@crclass)=split; 
      my $nd=$crclass[0]||1; 
      my $crclass=$nd;
      if(UPD7 and @crclass>1) { 
        $crclass=0; 
        $crclass |= 1 if($crclass[1]>0);  
        $crclass |= 2 if($crclass[2]>0);
        $crclass |= 4 if($crclass[3]>0);
        #? $crclass |= 8 if($crclass[4]>0); # how many cols?
        ## old        
        # my $z = (@crclass > 3) ? 1 : 0; # FIXME: crclass[0] may be nhit total, may want cols 1,2,3 if @crclass > 3
        # $crclass |= 1 if($crclass[$z+0]>0);  
        # $crclass |= 2 if($crclass[$z+1]>0);
        # $crclass |= 4 if($crclass[$z+2]>0);
        
        $cdsrdclass_in[1]++ if($crclass & 1);
        $cdsrdclass_in[2]++ if($crclass & 2);
        $cdsrdclass_in[3]++ if($crclass & 4);
      }
      # FIXME: may be full read id, want only SRR.num, 
      # FIXME: have mixed SRRsets, use cdsrdid{idp}[idn] 
      my($idp,$idn)=(0,0);
if(UPD7e) {
      ($idp,$idn)= readIdNum($rdid);
} else {
      if($isnum == 1) { $idn=$rdid; }
      elsif($isnum == 0) { ($idp,$idn)=split/\./,$rdid; }
      elsif($isnum == -1) {
        if($rdid =~ /^\d/) { $idn=$rdid; $isnum=1; } 
        elsif($rdid =~ /\.\d/) { ($idp,$idn)=split/\./,$rdid; $isnum=0; }
        else { $err++; last if($err>9); next; }
      }
}
      if($lastidp and $idp ne $lastidp){ $useRIDp++; } #? dont use soft switch on hash/array
      #o: if(USE_RIDp){ $cdsrdid{$idp}[$idn]=$nd; } else { $cdsrdid[$idn]=$nd; }
      #UPD7: replace val nd w/ cdsrdid class val
      if(USE_RIDp) { 
        if(UPD7e){ 
          if($cdsrdid{$idp}[$idn]){ $idp.="b"; }# FIX for sam cat of _1/_2 reads same id
        }
        $cdsrdid{$idp}[$idn]=$crclass; 
      } else { $cdsrdid[$idn]=$crclass; }
      
      $NCDS++; $lastidp=$idp;
    } close($inh);
    warn "# ncds readids=$NCDS from $cdstab \n" if $debug; 
  } else { warn "# warn: no cdstab: $cdstab\n" if $debug; }

  #---- new crlen from bam hdr ----
  my($NOK,$hasread,$haslen,%crok,%cread,%crlen,%crdtype)=(0,0,0);  
  #o my(%crclass,@crclasses); # for SAVE_RDIDS
  my($nidclass,$idclassh,$idclasslist)=(0); # renamed crclass
  
  if($inchrtab and ($ok,$inh)= openRead($inchrtab) and $ok) { 
    while(<$inh>){ next if(/^\W/); 
      my($cr,$crlen,$nread)=split;  #? UPD7 add crclass here? but chrtab now is output read counts table, no class info
      $cread{$cr}= $nread||0; $hasread++ if($nread>0);
      $crlen{$cr}= $crlen; $haslen++ if($haslen>0);  # may be only list of chr, crlen w/o nread
    } close($inh); 
  } else {
    open(IN,"samtools view -H $inbam |") or die "ERR: samtools view -H $inbam";
    while(<IN>) {
      #@SQ  SN:chrid LN:chrlen @SQ	SN:dmag20ug4d_sc000001	LN:1885227
      if(/^\@SQ/) {
        my($cr)= m/SN:(\S+)/; my($crlen)= (m/LN:(\d+)/)?$1:0;
        if($cr and $crlen > 0){ $crlen{$cr}= $crlen; $haslen++; }
      } elsif(/^\@/) {
        ;
      } else {
        last; # -H means only header ?
      }
    } close(IN);
  }

  #UPD7: my $CRTPAT='^[A-Za-z0-9]+'; #? TE|CDS|UNK chr class tag
  # ^ change this CRTPAT match to chrtag{cr} => tag hash lookup, make chrtag{} at top
  # if($SAVE_RDIDS and $cr ne $lcr){ my($crt)= ($cr =~ m/($CRTPAT)/)?$1:0; $crtlist{$crt}++; }
  # UPD7b : $CRTPAT usage is bad, unless caller sets it just right... change to no default,
  #       and create empty (UNK) class ?
  # * Need to check idclassh, cant allow > kMAXIDCLASS => 3..9 ? or process mem blows up
  
  #UPD7b: always call read_idclass, check kMAXIDCLASS and sam hdr IDs ~ $CRTPAT 
  ($nidclass,$idclassh,$idclasslist)= read_idclass($crclassf,0,\%crlen); # cds,te class by id, $idclass->{id} = class

  #UPD7a: no go
  # ($nidclass,$idclassh,$idclasslist)= ($crclassf) ? read_idclass($crclassf) : (0); # cds,te class by id, $idclass->{id} = class
  # if($nidclass==0 and $haslen){
  #   my(%crclass,%crcl);
  #   for my $cr (sort keys %crlen) {
  #     my($crt)= ($CRTPAT and $cr =~ m/($CRTPAT)/)? $1 : 'UNK'; 
  #     $crclass{$cr}= $crt;
  #   }
  #   unless(%crclass){ $crcl{'UNK'} = 0; }
  #   else { map{ $crcl{$crclass{$_}}++; } keys %crclass; } # or values %crclass > not uniq
  #   my @crclasses= sort keys %crcl;
  #   $nidclass= scalar keys %crclass; 
  #   $idclassh= \%crclass; # idclassh->{id}= class
  #   $idclasslist= \@crclasses;
  # }
  
  ## ORIG  
  # if($crclassf and ($ok,$inh)= openRead($crclassf) and $ok) { 
  #   while(<$inh>){ next if(/^\W/); 
  #     my($cr,$crclass)=split;   
  #     $crclass{$cr}= $crclass || 0;  
  #   } close($inh); 
  #   
  # } elsif($haslen) {
  #   for my $cr (sort keys %crlen) {
  #     my($crt)= ($cr =~ m/($CRTPAT)/)?$1:0; $crclass{$cr}= $crt;
  #   }
  # }
  # 
  # #bad not keys = id: @crclasses= sort keys %crclass; # UPD7: should be short, standard set of classes: CDS, TE, UNK ?
  # my %crcl=();
  # unless(%crclass){ $crcl{'UNK'} = 0; }
  # else { map{ $crcl{$crclass{$_}}++; } keys %crclass; } # or values %crclass > not uniq
  # @crclasses= sort keys %crcl;

   
  if(($hasread or $haslen) and ($ncpu>1 and $icpu<$ncpu)) {
    my @cr; 
    
    ## UPD7: replace this icpu part filter: crok and NOK for simpler readid modulus ncpu
    if($hasread>$ncpu) { @cr=sort{$cread{$b}<=>$cread{$a} or $a cmp $b} keys %cread; }
    else { @cr=sort{$crlen{$b}<=>$crlen{$a} or $a cmp $b} keys %crlen; }
    for my $i (0..$#cr){ if($icpu == ($i % $ncpu)){ $crok{$cr[$i]}=1; $NOK++; } }
    
    $outpart="$outtab.pt$icpu";
    $outchrpart= "$outchrpart.pt$icpu";
    $rdidtab .= ".pt$icpu" if($SAVE_RDIDS and $rdidtab); # rdidtabpart ?
    warn "# output part $icpu/$ncpu nchr=$NOK to $outpart\n" if($debug);
  }

  ## inbam must be read-ordered, ie all multimaps together, only 1 of paired reads counted now

  # my $sflag= 0x04 + (($SKIPR==1) ? 0x40 : 0x80); # 0x04 == no map DONT skip nomap, count
  #? other samt filter opts: -q minqual
  my $stopt = ($SKIPR==0)? "" : ($SKIPR==1) ? " -F 0x40" : " -F 0x80"; 
  my $STMINFILT=0;
  if($MINALN > 0) { $stopt .= " -m $MINALN"; $STMINFILT=1; } #upd7
  
  warn "# samtools view $stopt $inbam\n" if($debug);
  if($topcount != 0) {
    my $pheadtail= ($topcount < 0) ? "tail -n$topcount" : "head -n$topcount";
    open(IN,"samtools view $stopt $inbam | $pheadtail |") or die "ERR: samtools view $stopt $inbam";
  } else {
    open(IN,"samtools view $stopt $inbam |") or die "ERR: samtools view $stopt $inbam";
  }
  
  my ($ntmap)= readFilter(*IN,$STMINFILT);  
  close(IN);
  
  putcovtab( $outpart);
  putchrtab( $outchrpart); 
  putrdids($rdidtab) if($SAVE_RDIDS); 

  # UPD7: add merge parts code, separate main procedure or program

#-------------------------

=item add merge readid class tables

  -- input 2+ cdste_dnamap.readids tables w/ mapto class cols, 
     should be same columns, from opts="-savereadids -crclassf cds_te.idclass" bam=xxx.bam run_sam2covtab7iu.sh 
  .. output merged cdste_dnamap.readids, adds counts in columns w/ same readid
  .. reuse mergeparts(), run after merge parts, ie separate invoke
  
    #!/bin/bash
    ## mkreadids7b.sh
    #PBS -N mergerdids
    #PBS -l vmem=64gb,nodes=1:ppn=4,walltime=1:00:00
    #PBS -V
  
    rdt=dmag20sk4ma_tefam_SRR7825549cl_1a_bwa.readids
    rdc=dmag7fincds_SRR7825549cl_1a_bwa.readids
    rdo=dmag20cdste_SRR7825549cl_1a_bwa7b.readids
    
    if [ "X" = "X$datad" ]; then echo "ERR:datad=what?"; exit -1; fi
    cd $datad
    
    sort -k1,1 $rdt $rdc | perl -ne \
    'if(/^\W/) { print if(/#/ and not $hd++); } else { ($rd,@v)=split; if($lrd eq $rd){ 
      for $i (0..$#v){ $lv[$i]+=$v[$i]; } $nv++; } else{ putv() if($lrd); @lv=@v; $nv=1; } $lrd=$rd; } 
      END{ putv(); } sub putv{ print join("\t",$lrd,@lv)."\n"; @lv=(); }  ' \
      > $rdo

=cut

sub mergeparts {
  my($mergeflag,$ofile,@pt)=@_;
  
  unless(@pt) { # need @pt > 1 ?
    @pt=`ls $ofile.pt*`; chomp(@pt); 
  }
  my $npt=@pt;
  if(@pt<1){ warn "ERR: no parts to merge for $ofile\n"; return 0; } # need @pt > 1 ?
  rename($ofile, "$ofile.old") if(-s $ofile);
  
  #  also merge $rdidtab= "$iname.readids" parts from SAVE_RDIDS
  my $mtype= ($mergeflag and $mergeflag =~ /readids/)? 1 : 2;
    
  # table merge is adding nums in columns, need to read each part, keep all rows, 1st col = id key
  # there may be diff id keys in parts, tab is hash on ids x @val cols
  my %tab=(); my($nrow,$ncol,$ntin,$it)=(0,0,0);
  my(@thead,@tcomm);
  for my $pt (@pt){
    open(IN,$pt) or next; 
    ++$it; my($pti)= ($pt =~ m/\.(pt\w+)/)?$1:"i$it";
    while(<IN>){ 
      if(/^\W/) { 
        # chrtab has #totals, covtab has #info, preserve? at least col header
        if(/^#(ChrID|ReadID)/) { push @thead,$_ unless(@thead); }
        elsif(/^#/ and /n_mapok|total_locs|total_reads|cds_reads|minident=/){ s/#/#$pti./; push @tcomm,$_; }
        next;
      }
      # FIXME: sum tab[i] not for chr len, position .. are these same for all pt? pick max ?
      # covtab: chrid Pos sum-vals
      # chrtab: chrid chrlen sum-vals

      # FIXME2: hash key needs to be id.loc not id alone! for covtab, ok for chrtab?
      my($id, $cb, @vals);
      if($mtype == 1) { ($id,@vals)= split; $cb=1; } else { ($id, $cb, @vals)=split; }
      for my $i (0..$#vals) { $tab{$id}{$cb}[$i] += $vals[$i]; $ncol=$i if($i>$ncol);  }
      $ntin++;
      
      #... bad here..
      # my($id, @vals)=split; $ntin++;
      # my $tl= $tab{$id}[0]||0;
      # if($vals[0] > $tl){ $tab{$id}[0]= $vals[0]; } # non-sum: Pos/chrlen 
      # for my $i (1..$#vals) { $tab{$id}[$i] += $vals[$i]; $ncol=$i if($i>$ncol);  }
      
    } close(IN);
  }
       
  open(OUT,">$ofile"); 
  print OUT @thead if(@thead);
  #bad: for my $id (sort keys %tab) { my @vals= @{$tab{$id}}; print OUT join("\t",$id,@vals)."\n"; $nrow++; }
  for my $id (sort keys %tab) { 
    for my $cb (sort{ $a <=> $b } keys %{$tab{$id}} ) {
      my @vals= @{$tab{$id}{$cb}}; 
      if($mtype == 1) { print OUT join("\t",$id,@vals)."\n"; } # no cb for mtype == 1
      else { print OUT join("\t",$id,$cb,@vals)."\n"; }
      $nrow++; 
      }
    }
  print OUT @tcomm if(@tcomm);
  close(OUT);
  
  warn "#sam2covtab merge: nparts=$npt, nrows=$nrow, ncols=$ncol, ntin=$ntin to $ofile\n";
  return($nrow);
}

=item try1a upd

** NOTE: n_readid is total of view -m minaln bam, all other counts are for icpu part
   .. maybe want n_readid for icpu part only, so all are additive?
   .. Also, nomap == 0 with -m minaln use
   
head dmag19skasm_SRR7825549cl_1a_bwa.cdschr7a.covtab.pt1
#cdsxchr_covtab par: minident=0.65, mindupid=0.98, softclip=0, BIN=100, part=1/8 
#dmag19skasm_SRR7825549cl_1a_bwa.cdschr7a.covtab.pt1 n_readid 40419628, n_nomap 0, n_notcr 0, 
#  n_mapok 7997661, n_mapbad 4321056, n_dupbad 16887553, 
#  n_mismatch 26422783, n_insert 4655486, n_delete 3637119, n_softclip 0, n_intron 0
#ChrID	Pos	CovT	CovM	CovU	aCovT	aCovM	aCovU
dapmag19sk:LG2	100	0	0	0	1	1	1
dapmag19sk:LG2	200	0	0	0	4	4	4
dapmag19sk:LG2	300	0	0	0	9	9	9
dapmag19sk:LG2	400	0	0	0	12	12	12
dapmag19sk:LG2	500	0	0	0	10	10	10
dapmag19sk:LG2	600	0	0	0	6	5	5
dapmag19sk:LG2	700	0	0	0	8	7	7
...
dapmag19sk:scaffold8526_cov134	100	0	0	0	2	2	2
dapmag19sk:scaffold8526_cov134	200	0	0	0	2	2	2
dapmag19sk:scaffold4640_cov102	100	0	0	0	1	1	1
dapmag19sk:scaffold4640_cov102	200	0	0	0	2	2	2


head dmag19skasm_SRR7825549cl_1a_bwa.cdschr7a.chrtab.pt1
# nmap,nmult,nuniq = locations, nmread,nuniq,nomap = read counts
# n_notcr = 0
#ChrID        	chrlen  	nmap   	nuniq  	nmult  	nmread	nomap
dapmag19sk:LG2	16359456	3126127	479851	2646276	650830	0
dapmag19sk:LG1	14067088	1322635	382735	939900	455601	0
dapmag19sk:LG3	11088946	1167344	298312	869032	353424	0
dapmag19sk:LG7	10157464	5137049	453155	4683894	854496	0
dapmag19sk:LG5	10124675	2210566	252371	1958195	302225	0
...
# total_locs	122901750	26387247	3595539	22791708	4905372	0
# total_reads	40419628	mapt:4905372,12.14%	uniq:3595539,8.90%	mult:1309833,3.24%	nomap:0,0.00%	notcr:0,0.00%
---
#sam2covtab merge: nparts=8, nrows=4195, ncols=21, ntin=9841926 to dmag19skasm_SRR7825549cl_1a_bwa.cdschr7a.covtab 
#sam2covtab merge: nparts=8, nrows=4187, ncols=6, ntin=33319 to dmag19skasm_SRR7825549cl_1a_bwa.cdschr7a.chrtab

=cut


sub readFilter {  
  my($inhand,$STMINFILT)=@_;
  # my $sflag= 0x04 + (($SKIPR==1) ? 0x40 : 0x80);  
  # inhand == open(IN,"samtools view -F $sflag $inbam |") 
  my($lid, $iid, $ireadpart, $inrdlen, $nokid,$ntmap,$nmap,$crtlist)=(0) x 9;
  my @saver;
 
  if($SAVE_RDIDS){
    push @saveid, join("\t", '#ReadID','nmap',@$idclasslist); # add col header for @crclass
  }
  use constant rfQUICKEN => 1; 
  while(<$inhand>) {
    my @samx= split; 
    my($id,$fl,$cr,$cb,$cx,$cig)=@samx;
    
    if($id ne $lid) {
      if($lid) {
        if(UPD7f and $LONGR){  
        ($nmap,$crtlist)= ($nokid) ? putlongmap( $lid, $inrdlen, $nokid, \@saver) : (0,0); 
        } elsif(UPD7f) {
        ($nmap,$crtlist)= ($nokid) ? putshortmap( $lid, $inrdlen, $nokid, \@saver) : (0,0); 
        } else { # old
        ($nmap,$crtlist)= ($nokid) ? putmap2b( $lid, $nokid, \@saver) : (0,0); 
        }
        $ntmap += $nmap;
        if($SAVE_RDIDS and $nmap>0) { push @saveid, "$lid\t$nmap\t$crtlist"; }  
      }
      
      $lid=$id; $nokid=0; @saver=();      
      $ireadpart= ($ncpu > 1) ? $iid % $ncpu : 0;
      $iid++; 
      if($ireadpart == $icpu) { # always count new ids, but for nomap, or $n_readid++ unless($fl & 0xF00); # all 2nd maps
        $n_readid++;
        $inrdlen= length($samx[9]); # read seq, may be '*' or other placeholder
        $LONGR=1 if($inrdlen > $MAXSHORTRD); # FIXME; my $MAXSHORTRD= 500;
        $sum_readlen += $inrdlen; # new global counter, ave(rdlen)= sum_readlen/n_readid
      }
      #^^ UPD21: dont count unless($ireadpart == $icpu)
    }
    next unless($ireadpart == $icpu);  # test even for ncpu=1, icpu=1  
    
    # limit filters here .. need to do in putmap() w/ all aligns/read *

    my $badmap= 0;
    if($fl & 0x4){ $n_nomap++; $badmap=1; } ## next; #? count n_readid here
    if(rfQUICKEN and not $STMINFILT and not $LONGR) { # speed up some, skip low qual quick, MINALN set in putmap2b() 
      unless(UPD7e){ if($fl & 0x800){ $n_partb++;  next; } } # 2nd part map from bwa, dont count as dup map
      if($MINALN>0) {
        my $alen=0; while($cig =~ m/(\d+)M/g) { $alen += $1; }
        if($alen < $MINALN) { $n_mapbad++; $badmap=2; } ## next; 
      } 
    }    
    # UPD7maybe: for nomap,badmap, check/count cds_read match if flagged, for accurate total & num missing
    # @cdsrdclass_miss[1,2,3] counterpart of cdsrdclass_sum == hit reads
    if($badmap>0) { 
      if(UPD7 and $NCDS and $fl < 0x100) { 
        # FIXME over-counting same read 2x w/ badmap 2nd maps here, not so hits: only ($fl < 0x100): nomap and 1st map
        # my($idp,$idn)=split/\./,$id; 
        my($idp,$idn)= readIdNum($id); # UPD7e
        my $cdsread=0;
        $idn=int($idn); # ensure for srr123.456.sra style
        if(USE_RIDp){ $cdsread= $cdsrdid{$idp}[$idn]||0; } else { $cdsread=$cdsrdid[$idn]||0; } 
        if($cdsread>0) { my @clbit= (0,1,2,4); 
          $cdsrdclass_miss[4]++; # == any rdclass
          for my $i (1,2,3) { $cdsrdclass_miss[$i]++ if ($cdsread & $clbit[$i]); }
        } else { 
          $cdsrdclass_miss[0]++;  # this is for this-rd not in cdsrd set, dont need
        }
      } 
      next;
    }
    
    my($nmi)= (m/NM:i:(\d+)/)?$1:0;    # UPD7f saver: add inreadlen if > 9
    my $saver= [$cr,$cb,$cig,$nmi,$fl]; #? add $fl, check  0x800 supple from bwa == 2nd half map, maybe 
    push @saver, $saver;
    $nokid++ ;
    #old6: $nokid++ unless($NOK and not $crok{$cr}); # FIXME bug need this in putmap2b() also for @saver
        
  } close($inhand);

  if(UPD7f and $LONGR){  
  ($nmap,$crtlist)= ($nokid) ? putlongmap( $lid, $inrdlen, $nokid, \@saver) : (0,0); 
  } elsif(UPD7f) {
  ($nmap,$crtlist)= ($nokid) ? putshortmap( $lid, $inrdlen, $nokid, \@saver) : (0,0); 
  } else { # old
  ($nmap,$crtlist)= ($nokid) ? putmap2b( $lid, $nokid, \@saver) : (0,0); 
  }
    
  $ntmap+= $nmap;
  if($SAVE_RDIDS and $nmap>0) { push @saveid, "$lid\t$nmap\t$crtlist"; }  #NotE: crtlist is counts in @crclasses columns
  $n_readid++ if($ireadpart == $icpu); # always count new ids, but for nomap
  return($ntmap);
}


# gnodes_sam2covtab.lrinc.pl

sub putlongmap {
  my($rdid, $rdlen, $inmap, $saveset)=@_;
  our(%crmax,%cov,%covt,%covu,%allcov,%allcovt,%allcovu); # ,@savemap,@saveid

  my $cdsread=  0; # ($NCDS) ? $cdsrdid[$idn] : 0;
  if($NCDS){ 
    my($idp,$idn)= readIdNum($rdid); # UPD7e
    if(USE_RIDp){ $cdsread= $cdsrdid{$idp}[$idn]||0; } else { $cdsread=$cdsrdid[$idn]||0; } 
  } 
  
  my(@thismap,%crl,%crtlist);
  my($nmap,$topaln,$lcr)=(0) x 9; # ndi   
  
  for my $saver (@$saveset) {
    my($cr,$cb,$cig,$nmi,$sfl)= @$saver;

    my($alen,$lenc,$cend,$softclip)= addCigar($cr,$cb,$cig,$sfl,$nmi,$rdlen,$inmap,$cdsread);  
    #  my($alen,$lenc,$cend,$softclip)= smokeCigar($cig,$cb);
    
    if($nmi>0){ $alen -= $nmi; $n_mismatch += $nmi; } 
    if($ALLOW_SOFTCLIP and $softclip>0) { $lenc -= $softclip; $n_softclip++; }
    
    my $mapbad=0;    
    if($mapbad) {  
    
    } elsif($topaln == 0) { 

      $nsample++; 
      if($nsample < NSAMPLE) {
        push @sreadlen, $lenc; 
      } elsif($nsample == NSAMPLE or $readlen<1) {
        @sreadlen = sort{$b <=> $a} @sreadlen; $readlen= $sreadlen[ int($nsample/2) ];
        my $areadlen= ($n_readid<2) ? $readlen : int(0.5 + $sum_readlen/$n_readid); # UPD7f
        # if($areadlen != $readlen){ } #some warning, use which?
        #UPD7f.no: $MINALN= int( $readlen * $MIN_IDENT); # was $MINALN= $sminaln / $nsample; 
        $MINALN= 1 if($MINALN<1); #?? or not
        warn "# readlen=$readlen, avelen=$areadlen, minaln=$MINALN from $nsample\n" if($debug);
      }

      $topaln=$alen;  $crdtype{$cr}{'mread'}++; 
      if($SAVE_RDIDS and $cr ne $lcr){ my($crt)= $idclassh->{$cr} || 0; $crtlist{$crt}++; }
    } else {
      # $mapbad=2 if( $mapbad or $alen < $topaln * $MIN_DUPIDENT);  # count? n_dupbad   
    }
    
    if($mapbad){ 
      if($mapbad==2){ $n_dupbad++; } else { $n_mapbad++; }  
    } else { 
      $nmap++; 
      $n_mapok++; # num maps/read id, count after mapbad, before crok
    }

    $lcr= $cr;      
  }  
  
  if($SAVE_RDIDS) {
    my $crtlist = join"\t", map{ $crtlist{$_}||0 } @$idclasslist; #<<< FIXME, need counts/chrclass for output cols
    return($nmap,$crtlist); 
  }

  if(UPD7 and $NCDS) { 
    my @clbit= (0,1,2,4);
    if($cdsread>0) { 
      $cdsrdclass_sum[4]++; # == hit any rdclass
      for my $i (1,2,3) { $cdsrdclass_sum[$i]++ if ($cdsread & $clbit[$i]); }
    } else { 
      $cdsrdclass_sum[0]++; # this rd not in any cds-rdclass
    }  
  }

  return($nmap,0);
}

  
=item addCigar update
  
  my($alen,$lenc,$cend,$softclip)= 
    addCigar($cr,$cb,$cigar,$fl,$nmi,$rdlen,$nmap,$cdsread);  
    
  fill global tables: %crmax,%cov,%covt,%covu,%allcov,%allcovt,%allcovu
    
  -- need to count M align per base for very long cig from long/poor reads: nMnDnMnI... 1000s long
  -- fill in cov bin tables from M/base count
  -- replace smokeCigar w/ addCigar()  for short/hiqual and long reads

=cut

sub addCigar { 
  my($cr,$cb,$cigar,$fl,$nmi,$rdlen,$nmap,$cdsread)=@_;
   # $cr,$cb,$cig,$sfl,$nmi,$rdlen,$inmap,$cdsread
  our(%crmax,%cov,%covt,%covu,%allcov,%allcovt,%allcovu); # ,@savemap,@saveid
  my($alen,$lenc,$cend,$nlen,$indel,$softclip)= (0) x 9;
  # cend == ci, alen == cm, bi == w # $cm+=$w; 

  # $cb--; $cend=$cb; # FIXME: tb/cend is 1-base, need 0-base for [cend+i] below
  $cend=0;
  my @alnb= (0) x $rdlen;
  while($cigar =~ m/(\d+)([A-Z])/g) { 
    my($bi,$bt)=($1,$2); # my($w,$t)=($1,$2); 
 
    # unless($bt eq 'H' or $bt eq 'S' or $bt eq 'I') {  $nlen+=$bi; }#? include both I/D indels for N? not I
     
    if($bt eq 'H') { 
      $lenc += $bi; $bi=0; 
    } elsif($bt eq 'S') {
      $softclip += $bi; $lenc += $bi; $bi=0; 
    } elsif($bt eq 'M') {
      for(my $i=0;$i<$bi;$i++){ $alnb[$cend+$i]++; }
      $alen+= $bi; $lenc += $bi;  $cend += $bi;  
    } elsif($bt eq 'N') { 
      $n_intron++;  $cend += $bi; # cend not changed for HISP; 
    } elsif($bt eq 'I') {
      #x $n_insert++; 
      $n_insert += $bi; #UPD7f change for base count eq n_mismatch += nmi
      $lenc += $bi; $bi=0; $indel++;
    } elsif($bt eq 'D') {  # P also?
      #x $n_delete++; 
      $n_delete += $bi; #UPD7f change for base count
      $cend += $bi; $indel++;
    } elsif($bt eq 'P') {  # what is P now?
      $cend += $bi;
    } else {
      # unknown here, what? record?
    }
  }

  my $rv=1; my $uv=0; my $rty='zero';
  if($nmap==1) { $uv=1; $rty='uniq' } 
  elsif($nmap>1) { $rv= 1.0/$nmap; $rty='mult'; }   
  
  # my($cb,$ce)=($cb,$cb+$cend); #?  my($cr,$cb,$ce)= @$cbe; 
  my $be= int(($cb+$cend)/$BN);   # my $bb=int($cb/$BN);
  $crdtype{$cr}{$rty} += $nmap;
  $crmax{$cr}=$be if($be>$crmax{$cr});
  my $bdepth= 1.0/$BN; # per base per bin depth 
  my $rvdepth= $rv * $bdepth;
  my $uvdepth= $uv * $bdepth;
  
  # if($SAVE_RDIDS and $cr ne $lcr){ my($crt)= $idclassh->{$cr} || 0; $crtlist{$crt}++; }
  for(my $i=0; $i<=$cend; $i++) {
    next unless($alnb[$i]>0);
    # my $bdepth=$alnb[$i] or next; # only 1 or 0 now
    # $bdepth /= $BN; # should be 1/100 fixed val
    my $bi= int( ($cb+$i)/$BN );
    $allcov{$cr}[$bi]  += $rvdepth; 
    $allcovu{$cr}[$bi] += $uvdepth; 
    $allcovt{$cr}[$bi] += $bdepth;  
    if($cdsread > 0) {
      $covt{$cr}[$bi] += $bdepth if($cdsread & 1); # 1st col
      $cov{$cr}[$bi]  += $bdepth if($cdsread & 2); # 2nd col
      $covu{$cr}[$bi] += $bdepth if($cdsread & 4); # 3rd col .. more cols? later
    }
  }
      
  return($alen,$lenc,$cend,$softclip);#?
}

# gnodes_sam2covtab.cigdepth7f.pl

=item UPD7f addCigarDepth 
  
  -- for multipmap reads, count depth per read base pos, then collapse to cov bins with per-base depth
  
  my @adepth= (0) x $rdlen;
  my $nmap= @$saveset; # minus badmap, esp dups < topalign*mindupident
  
  for my $rdmap (@$saveset) {
    my($cr,$cb,$cigar,$nmi,$sfl)= @$rdmap;
    my($alen,$lenc,$cend,$softclip)= 
      addCigarDepth( \@adepth, $cr,$cb,$cigar,$sfl,$nmi,$rdlen,$nmap,$cdsread);  
  }
    
  fill global tables: %crmax,%cov,%covt,%covu,%allcov,%allcovt,%allcovu

=cut

sub addCigarDepth { 
  my( $addToBins, $aligndepth, $cr,$cb,$cigar,$rdlen,$nmap,$cdsread)=@_;
  our(%crmax,%cov,%covt,%covu,%allcov,%allcovt,%allcovu); # ,@savemap,@saveid
  my($alen,$lenc,$cend,$nlen,$indel,$softclip)= (0) x 9;
  # cend == ci, alen == cm, bi == w # $cm+=$w; 
 
  my @thisdepth= (0) x $rdlen;
  while($cigar =~ m/(\d+)([A-Z])/g) { 
    my($bi,$bt)=($1,$2);
     
    if($bt eq 'H') { 
      $lenc += $bi; $bi=0; 
    } elsif($bt eq 'S') {
      $softclip += $bi; $lenc += $bi; $bi=0; 
    } elsif($bt eq 'M') {
      for(my $i=0;$i<$bi;$i++){ $thisdepth[$cend+$i]=1; } #  $aligndepth->[$cend+$i]++; 
      $alen+= $bi; $lenc += $bi;  $cend += $bi;  
    } elsif($bt eq 'N') { 
      $n_intron++;  $cend += $bi; # cend not changed for HISP; 
    } elsif($bt eq 'I') {
      $n_insert += $bi; #UPD7f change for base count eq n_mismatch += nmi
      $lenc += $bi; $bi=0; $indel++;
    } elsif($bt eq 'D') {  # P also?
      $n_delete += $bi; #UPD7f change for base count
      $cend += $bi; $indel++;
    } elsif($bt eq 'P') {  # what is P now?
      $cend += $bi;
    } else {
      # unknown here, what? record?
    }
  }

  if($addToBins) {
    # @$aligndepth should be sum of thisdepth over mappings, use for rvcov depth

    my $rv=1; my $uv=0; my $rty='zero';
    if($nmap==1) { $uv=1; $rty='uniq' } 
    elsif($nmap>1) { $rv= 1.0/$nmap; $rty='mult'; }   
    
    my $be= int(($cb+$cend)/$BN);  
    $crdtype{$cr}{$rty} += $nmap;
    $crmax{$cr}=$be if($be>$crmax{$cr});

    for(my $i=0; $i<=$cend; $i++) {
      next unless($thisdepth[$i]); # 0 or 1
      my $adepth= $aligndepth->[$i]||1; # 1..nmap
      my $bi= int( ($cb+$i)/$BN );
      
      my ($rvd,$uvd)=(0,0);
      if($adepth == 1){ 
        $uvd= $rvd= 1.0/$BN; ## $bdepth;  # == 1/BN
        $allcovu{$cr}[$bi] += $uvd; 
      } else { 
        $rvd= (1.0/$adepth) / $BN; #?? 1/$BN; # or $bdepth; ??? rvd here is problem, for multimap always 1/BIN?
      }
      
      $allcov{$cr}[$bi]  += $rvd; 
      
      #bad? my $bdepth= $adepth/$BN; #? BUG Here? map2b() uses bspan/BN, using adepth is wrong?
      my $bdepth= 1.0/$BN; # map2b() equiv $bspan/$BN, allcovt counts reads at this (cr,cb+i) base / BINSIZE
      $allcovt{$cr}[$bi] += $bdepth;  
        # addToBins call for each map of rd, but adepth includes each, ie N x N
        # .. no, this is ok: each map of rd has diff cr,cb .. but this new calc has some much larger total-dupl vals,
        # .. ** map2b() sets bdepth=bspan/BN, ignoring nmap/adepth .. want that again? probably
        # .. where rvd is near same
        
      if($cdsread > 0) {
        $covt{$cr}[$bi] += $bdepth if($cdsread & 1); # 1st col
        $cov{$cr}[$bi]  += $bdepth if($cdsread & 2); # 2nd col
        $covu{$cr}[$bi] += $bdepth if($cdsread & 4); # 3rd col .. more cols? later
      }
    }

  }
  
  # FIXME: addCigarDepth returns cend == readlen, 0-base not cb-base, change from prior
  return($alen,$lenc,$cend,$softclip, \@thisdepth);#?
}

=item UPD7f putshort examle changed covtab

  OLD daphpulex_pa42v2_SRR3090572_b2_bwa.cc7d.covtab
#ChrID	Pos	CovT	CovM	CovU	aCovT	aCovM	aCovU
dpulex19pa42b:scaffold_1	100	0	0	0	28	7	0
dpulex19pa42b:scaffold_1	200	0	0	0	76	27	0
dpulex19pa42b:scaffold_1	300	0	0	0	101	34	0
dpulex19pa42b:scaffold_1	400	0	0	0	83	27	0
dpulex19pa42b:scaffold_1	500	0	0	0	67	22	0
   NEW cc7fdup9hi.covtab
#ChrID	Pos	CovT	CovM	CovU	aCovT	aCovM	aCovU
dpulex19pa42b:scaffold_1	100	0	0	0	1	0	0
dpulex19pa42b:scaffold_1	200	0	0	0	66	7	1
dpulex19pa42b:scaffold_1	300	0	0	0	113	13	1
dpulex19pa42b:scaffold_1	400	0	0	0	92	12	1
dpulex19pa42b:scaffold_1	500	0	0	0	32	0	0   << Odd result aCovT>0, aCovM == 0 (ie < 0.5)
..

OLD: grep 'scaffold_1    21.00   '  daphpulex_pa42v2_SRR3090572_b2_bwa.cc7d.covtab
#ChrID	Pos	CovT	CovM	CovU	aCovT	aCovM	aCovU
dpulex19pa42b:scaffold_1	21000	1980	0	0	3081	167	0
dpulex19pa42b:scaffold_1	21100	4302	0	0	6621	344	2
dpulex19pa42b:scaffold_1	21200	6270	0	0	7861	408	2
dpulex19pa42b:scaffold_1	21300	7397	0	0	7600	386	0
dpulex19pa42b:scaffold_1	21400	5995	0	0	7501	382	0
dpulex19pa42b:scaffold_1	21500	3024	0	0	5381	280	0
dpulex19pa42b:scaffold_1	21600	919 	0	0	2355	132	0
dpulex19pa42b:scaffold_1	21700	2856	0	0	3228	180	1
dpulex19pa42b:scaffold_1	21800	7900	0	0	8611	481	2
dpulex19pa42b:scaffold_1	21900	8606	0	0	11796	690	2

NEW: grep 'scaffold_1    21.00   '  daphpulex_pa42v2_SRR3090572_b2_bwa.cc7fdup9hi.covtab
#ChrID	Pos	CovT	CovM	CovU	aCovT	aCovM	aCovU
dpulex19pa42b:scaffold_1	21000	1648	0	0	5896	78	54
dpulex19pa42b:scaffold_1	21100	7870	0	0	26129	139	47
dpulex19pa42b:scaffold_1	21200	7635	0	0	16885	110	41
dpulex19pa42b:scaffold_1	21300	29235	0	0	30304	140	41
dpulex19pa42b:scaffold_1	21400	50463	0	0	71071	247	23
dpulex19pa42b:scaffold_1	21500	36218	0	0	64631	210	3
dpulex19pa42b:scaffold_1	21600	6075	0	0	16417	73	13
dpulex19pa42b:scaffold_1	21700	3894	0	0	4582	154	115
dpulex19pa42b:scaffold_1	21800	40107	0	0	45000	450	178
dpulex19pa42b:scaffold_1	21900	48713	0	0	64963	690	286

OLD: grep 'scaffold_1    9.00    '  daphpulex_pa42v2_SRR3090572_b2_bwa.cc7d.covtab
dpulex19pa42b:scaffold_1	9000	49	0	0	58	21	0
dpulex19pa42b:scaffold_1	9100	21	0	0	31	9	0
dpulex19pa42b:scaffold_1	9200	10	0	0	19	5	0
dpulex19pa42b:scaffold_1	9300	23	0	0	28	9	0
dpulex19pa42b:scaffold_1	9400	38	0	0	49	17	0
dpulex19pa42b:scaffold_1	9500	51	0	0	66	21	0
dpulex19pa42b:scaffold_1	9600	54	0	0	67	23	0
dpulex19pa42b:scaffold_1	9700	60	0	0	84	28	0
dpulex19pa42b:scaffold_1	9800	31	0	0	82	27	0
dpulex19pa42b:scaffold_1	9900	14	0	0	81	23	0

NEW: grep 'scaffold_1    9.00    '  daphpulex_pa42v2_SRR3090572_b2_bwa.cc7fdup9hi.covtab
dpulex19pa42b:scaffold_1	9000	121	0	0	132	20	1
dpulex19pa42b:scaffold_1	9100	35	0	0	47	12	3
dpulex19pa42b:scaffold_1	9200	11	0	0	21	11	9
dpulex19pa42b:scaffold_1	9300	39	0	0	43	14	11
dpulex19pa42b:scaffold_1	9400	61	0	0	69	23	16
dpulex19pa42b:scaffold_1	9500	57	0	0	63	32	30
dpulex19pa42b:scaffold_1	9600	32	0	0	38	38	38
dpulex19pa42b:scaffold_1	9700	31	0	0	43	41	41
dpulex19pa42b:scaffold_1	9800	23	0	0	77	42	26
dpulex19pa42b:scaffold_1	9900	32	0	0	159	35	5
=cut

=item UPD7f putshortmap 
  
  .. uses addCigarDepth() twopass on read saveset to calc 
      1st. valid aligns among multimaps, 2nd. cov depth/valids
  .. replaces putmap2b
  
=cut

sub putshortmap { 
  my($rdid, $rdlen, $nokid, $saveset)=@_;
  our(%crmax,%cov,%covt,%covu,%allcov,%allcovt,%allcovu); # ,@savemap,@saveid

  my($idp,$idn)= readIdNum($rdid); # UPD7e
  my $cdsread=  0; # ($NCDS) ? $cdsrdid[$idn] : 0;
  if($NCDS){ if(USE_RIDp){ $cdsread= $cdsrdid{$idp}[$idn]||0; } else { $cdsread=$cdsrdid[$idn]||0; } } 
  
  my(@thismap,@aligndepth,%crl); @aligndepth=();
  my($nmap,$topaln,$mindupaln,$lcr)=(0) x 9; # ndi   
  my $nmapin= @$saveset;
  use constant { ADD1ST => 0, ADD2BINS => 2 };
  
  for my $saver (@$saveset) {
    my($cr,$cb,$cigar,$nmi,$sfl)= @$saver;
    
    # FIXME: addCigarDepth returns cend == readlen, 0-base not cb-base, change from prior
    my($alen,$lenc,$cend,$softclip, $thisdepth)= 
      addCigarDepth(ADD1ST, [], $cr,$cb,$cigar,$rdlen,$nmapin,$cdsread); # nmap changes in loop
    #o: my($alen,$lenc,$cend,$softclip)= smokeCigar($cig,$cb);
    
    if($nmi>0){ $alen -= $nmi; $n_mismatch += $nmi; } 
    if($ALLOW_SOFTCLIP and $softclip>0) { $lenc -= $softclip; $n_softclip++; }
    
    my $mapbad=0;    
    if($MINALN>0) { 
      $mapbad=1 if($alen < $MINALN); 
    } else { # sample lenc
      my $mina= $lenc * $MIN_IDENT; 
      $mapbad=1 if($alen < $mina);
      if($topaln == 0) { # sample only 1st map
      $nsample++; push @sreadlen, $lenc; 
      if($nsample >= NSAMPLE) {
        @sreadlen = sort{$b <=> $a} @sreadlen; $readlen= $sreadlen[ int($nsample/2) ];
        my $areadlen= ($n_readid<2) ? $readlen : int(0.5 + $sum_readlen/$n_readid); # UPD7f
        # if($areadlen != $readlen){ } #some warning, use which?
        $MINALN= int( $readlen * $MIN_IDENT); # was $MINALN= $sminaln / $nsample; 
        $MINALN= 1 if($MINALN<1);
        warn "# readlen=$readlen, avelen=$areadlen, minaln=$MINALN from $nsample\n" if($debug);
        }
      }
    }

    if($mapbad) { } # added, maybe 2nd map better?
    elsif($topaln == 0) { 
      $topaln=$alen; $mindupaln= $topaln * $MIN_DUPIDENT;
      $crdtype{$cr}{'mread'}++; 
    } else {
      $mapbad=2 if( $mapbad or $alen < $mindupaln); # == $topaln * $MIN_DUPIDENT # count? n_dupbad   
    }
    
    if($mapbad){ 
      if($mapbad==2){ $n_dupbad++; } else { $n_mapbad++; }  
    } else { 
      $nmap++; $n_mapok++; # num maps/read id, count after mapbad, before crok
      # here, need to skip not crok, but count in nmap
      for(my $i=0; $i<$cend; $i++) { $aligndepth[$i] += $thisdepth->[$i]; }
      #OLD: push @thismap, [$cr,$cb,$cend]; #UPD7
      push @thismap, $saver; # NEW for ADD2BINS
    }

    $lcr= $cr;      
  } # @saveid
  
  # $nmap == @thismap;
  $lcr=0; my(%crtlist);
  for my $saver (@thismap) {
    my($cr,$cb,$cigar,$nmi,$sfl)= @$saver;
    
    my($alen,$lenc,$cend,$softclip, $thisdepth)= 
      addCigarDepth(ADD2BINS, \@aligndepth, $cr,$cb,$cigar,$rdlen,$nmap,$cdsread); 

    if($SAVE_RDIDS and $cr ne $lcr){ my($crt)= $idclassh->{$cr} || 0; $crtlist{$crt}++; }
    $lcr=$cr;
  }
    
  if($SAVE_RDIDS) {
    my $crtlist = join"\t", map{ $crtlist{$_}||0 } @$idclasslist; #<<< FIXME, need counts/chrclass for output cols
    return($nmap,$crtlist);
  } else {
    if(UPD7 and $NCDS) { 
      my @clbit= (0,1,2,4);
      if($cdsread>0) { 
        $cdsrdclass_sum[4]++; # == hit any rdclass
        for my $i (1,2,3) { $cdsrdclass_sum[$i]++ if ($cdsread & $clbit[$i]); }
      } else { 
        $cdsrdclass_sum[0]++; # this rd not in any cds-rdclass
      }  
    }
    return($nmap,0);
  }
}

sub smokeCigar {  
  my($cigar,$cbeg)= @_;
  my($alen,$lenc,$cend,$softclip)= (0) x 9;
  $cend= $cbeg;

  while($cigar =~ m/(\d+)([A-Z])/g) { 
    my($bi,$bt)=($1,$2);
    if($bt eq 'M') {
      $alen+= $bi; $lenc += $bi;  $cend += $bi;  
    } elsif($bt eq 'N') { 
      $n_intron++;  $cend += $bi; # cend not changed for HISP; 
    } elsif($bt eq 'H') { 
      $lenc += $bi; $bi=0; 
    } elsif($bt eq 'S') {
      $softclip += $bi; $lenc += $bi; $bi=0; 
    } elsif($bt eq 'I') {
      #x $n_insert++; 
      $n_insert += $bi;#UPD7f change for base count
      $lenc += $bi; $bi=0; 
    } elsif($bt eq 'D') {  # P also?
      #x $n_delete++; 
      $n_delete += $bi; #UPD7f change for base count, eq n_mismatch += nmi
      $cend += $bi;
    } elsif($bt eq 'P') {  
      $cend += $bi;
    } else {
      # unknown here, what? record?
    }
  }
  
  return($alen,$lenc,$cend,$softclip);
}

sub putmap2b {
  my($rdid, $nokid, $saveset)=@_;
  our(%crmax,%cov,%covt,%covu,%allcov,%allcovt,%allcovu); # ,@savemap,@saveid

  # $n_readid++; # count before this
  #x my($idp,$idn)=split/\./,$rdid; 
  my($idp,$idn)= readIdNum($rdid); # UPD7e
  my $cdsread=  0; # ($NCDS) ? $cdsrdid[$idn] : 0;
  if($NCDS){ if(USE_RIDp){ $cdsread= $cdsrdid{$idp}[$idn]||0; } else { $cdsread=$cdsrdid[$idn]||0; } } 
  #UPD7: change meaning of cdsread => crclass: 0,1,2,4,... bit flags
  
  my(@thismap,%crl);
  my($nmap,$topaln,$mindupaln,$lcr)=(0) x 9; # ndi   
  # comput nmap, cdsread from @saveid all
  for my $saver (@$saveset) {
    my($cr,$cb,$cig,$nmi,$sfl)= @$saver;
    # my $saver= [$cr,$cb,$cig,$nmi,$fl]; #? add $fl, check  0x800 supple from bwa == 2nd half map, maybe 
    
    unless(UPD7e) { # upd7e count 2nd part maps, important variable?
    if($sfl & 0x800) { # 2nd part map from bwa, dont count as dup map
      $n_partb++;  next;  # skip for now ?? or add align to 1st part if same chr eq $lcr
    }
    }
    
    my($alen,$lenc,$cend,$softclip)= smokeCigar($cig,$cb);
   
    #UPD to this, replacing @thismap, %allcov table fill below
    # my($alen,$lenc,$cend,$softclip)= 
    #   addCigar($cr,$cb,$cigar,$fl,$nmi,$rdlen,$nmap,$cdsread); 
    
    if($nmi>0){ $alen -= $nmi; $n_mismatch += $nmi; } 

    if($ALLOW_SOFTCLIP and $softclip>0) { $lenc -= $softclip; $n_softclip++; }
    
    #UPD add/change to mapbad= alen < $minaln; set minaln from read-size * MIN_IDENT sampled here?
    # my $mapbad= ( $lenc < 1 || $alen < $lenc * $MIN_IDENT ) ? 1 : 0;
    my $mapbad=0;    
    if($MINALN>0) { 
      $mapbad=1 if($alen < $MINALN); 
    } else { # sample lenc
      my $mina= $lenc * $MIN_IDENT; 
      $mapbad=1 if($alen < $mina);
      if($topaln == 0) { # sample only 1st map
      $nsample++; push @sreadlen, $lenc; 
      if($nsample >= NSAMPLE) {
        # $readlen= int(0.5 + $sreadlen/ $nsample );
        @sreadlen = sort{$b <=> $a} @sreadlen; $readlen= $sreadlen[ int($nsample/2) ];
        my $areadlen= ($n_readid<2) ? $readlen : int(0.5 + $sum_readlen/$n_readid); # UPD7f
        # if($areadlen != $readlen){ } #some warning, use which?
        $MINALN= int( $readlen * $MIN_IDENT); # was $MINALN= $sminaln / $nsample; 
        $MINALN= 1 if($MINALN<1);
        warn "# readlen=$readlen, avelen=$areadlen, minaln=$MINALN from $nsample\n" if($debug);
        }
      }
    }

    if($mapbad) { } # added, maybe 2nd map better?
    elsif($topaln == 0) { 
      $topaln=$alen; $mindupaln= $topaln * $MIN_DUPIDENT;
      $crdtype{$cr}{'mread'}++; 
    } else {
      $mapbad=2 if( $mapbad or $alen < $mindupaln); # == $topaln * $MIN_DUPIDENT # count? n_dupbad   
      #? change $topaln * $MIN_DUPIDENT to $topaln - $MIN_DUPOFFBY == 1,2 ; offby1,2 ? 
    }
    
    if($mapbad){ 
      if($mapbad==2){ $n_dupbad++; } else { $n_mapbad++; }  
    } else { 
      $nmap++; $n_mapok++; # num maps/read id, count after mapbad, before crok
      # here, need to skip not crok, but count in nmap
      push @thismap, [$cr,$cb,$cend]; #UPD7
      #o6: unless($NOK and not $crok{$cr}) { push @thismap, [$cr,$cb,$cend]; }
    }

    $lcr= $cr;      
  } # @saveid
  
  #UPD3d: dont use precalc ndup, count up from inputs, same id eq lid, collect then store at next id
  my $rv=1; my $uv=0; my $rty='zero'; # nmap == @thismap
  if($nmap==1) { $uv=1; $rty='uniq' } elsif($nmap>1) { $rv= 1.0/$nmap; $rty='mult'; }   #? skip if read not counted here?
  # for my $c (keys %crl){ $crdtype{$c}{$rty} += $nmap; } 
  
=item vers3h
  -- correct bin readcov value: adjust for base-level span in bin, ie
    binb=200, bine=300, and cb=220, ce=279, then cov values * cw=1+279-220=60,
      covt[binb] += cw; covu[bbin] += uv * cw; covm[bbin] += rv * cw;
    AND adjust output by /BN bin size, or else cw*BN here to get per-base read depth in bins
=cut

  $lcr=0; my(%crtlist);
  for my $cbe (@thismap) {
    my($cr,$cb,$ce)= @$cbe; 
    # my($cr,$cb,$cig,$nmi)= @$saver; # NOT this
    $crdtype{$cr}{$rty} += $nmap;
    my $bb=int($cb/$BN); my $be=int($ce/$BN);  
    $crmax{$cr}=$be if($be>$crmax{$cr});
    #UPD7: 
    # my $CRTPAT='^[A-Za-z0-9]+'; #? TE|CDS|UNK chr class tag
    # ^ change this CRTPAT match to chrtag{cr} => tag hash lookup, make chrtag{} at top
    #o: if($SAVE_RDIDS and $cr ne $lcr){ my($crt)= ($cr =~ m/($CRTPAT)/)?$1:0; $crtlist{$crt}++; }
    if($SAVE_RDIDS and $cr ne $lcr){ my($crt)= $idclassh->{$cr} || 0; $crtlist{$crt}++; }

    $lcr=$cr;
    for(my $i=$bb; $i<=$be; $i++) { 
      my $ib= 1 + $BN*$i; my $ie=$BN + $ib - 1; #  1..100, 201..300
      my $bspan= 1 + (($ce < $ie)? $ce : $ie) - (($cb > $ib)? $cb : $ib); #? use this int val, or fract bdepth?
      my $bdepth= $bspan/$BN; # 1 for cb,ce outside bin, <1 for inside .. this is per-base depth, for all of bin
      my $rvd= $bdepth * $rv; my $uvd= $bdepth * $uv;
      $allcov{$cr}[$i]  += $rvd; 
      $allcovu{$cr}[$i] += $uvd; 
      $allcovt{$cr}[$i] += $bdepth;  

    #UPD7: change meaning of cdsread => crclass: 0,1,2,4,... bit flags    
      if($cdsread > 0) {
        #upd7: replace these w/ $crcov{$cr}{$crclass}[$i] += $bdepth;
        if(UPD7) {
          $covt{$cr}[$i] += $bdepth if($cdsread & 1); # 1st col
          $cov{$cr}[$i]  += $bdepth if($cdsread & 2); # 2nd col
          $covu{$cr}[$i] += $bdepth if($cdsread & 4); # 3rd col .. more cols? later
        } else {
          $cov{$cr}[$i]  += $rvd; 
          $covu{$cr}[$i] += $uvd; 
          $covt{$cr}[$i] += $bdepth;  
        }
      }
    }
  }

  #o: return $nmap;
  if($SAVE_RDIDS) {
    #? my $crtlist = join",", sort{ $crtlist{$b} <=> $crtlist{$a} } keys %crtlist;
    my $crtlist = join"\t", map{ $crtlist{$_}||0 } @$idclasslist; #<<< FIXME, need counts/chrclass for output cols
    #? my $crtlist = join",", sort keys %crtlist;
    return($nmap,$crtlist);
  } else {
  
    if(UPD7 and $NCDS) { 
      my @clbit= (0,1,2,4);
      if($cdsread>0) { 
        $cdsrdclass_sum[4]++; # == hit any rdclass
        for my $i (1,2,3) { $cdsrdclass_sum[$i]++ if ($cdsread & $clbit[$i]); }
      } else { 
        $cdsrdclass_sum[0]++; # this rd not in any cds-rdclass
      }  
    }
      
    return($nmap,0);
  }
}


# FIXME for duptab.gz at least
sub openRead{ my($fn)=@_; my $fh;
  $fn.=".gz" unless(-s $fn or $fn=~/.gz$/);
  my $ok= (! -s $fn)? 0 : ($fn=~/.gz$/) ? open($fh,"gunzip -c $fn |") : open($fh,$fn);
  warn "#openRead($fn) = $ok\n" if $debug;
  return($ok,$fh); 
}

sub read_idclass {
  my($crclassf,$allclasses,$crlenh)=@_; 
  my($nid,$ncl,$ok,$inh)=(0,0,0);
  my %crclass=();   my @crclass=();
  $allclasses||=0;
  
  if($crclassf and ($ok,$inh)= openRead($crclassf) and $ok) { 
    while(<$inh>){ next if(/^\W/); 
      my($cr,$crclass,@clx)=split;  # may have more columns .. keep all ie CDS,BUSCO 
      if($allclasses and @clx){ $crclass=join(" ",$crclass,@clx); }
      $crclass{$cr}= $crclass || 0;  
    } close($inh); 
    $nid= scalar(keys %crclass);
  }
  #.. insert other way from sam hdr ids and CRTPAT?
  if($nid==0 and $CRTPAT and ref($crlenh)){
    for my $cr (sort keys %$crlenh) {
      my($crt)= ($cr =~ m/($CRTPAT)/)? $1 : 'UNK'; 
      $crclass{$cr}= $crt;
    }
    $nid= scalar(keys %crclass);
  }
  
  unless($nid>0) { @crclass=('UNK'); }
  else {
    my %crcl=(); # make list of values = classes
    map{ $crcl{$crclass{$_}}++; } keys %crclass;  # or values %crclass > not uniq
    @crclass= sort keys %crcl; # sort by count?
  }
  
  $ncl= @crclass; 
  if($ncl > kMAXIDCLASS) { # problem, cancel..
    warn "#ERR too many idclasses n=$ncl from nid=$nid of $crclassf or sam ids x CRTPAT='$CRTPAT' \n";
    $nid=0; %crclass=(); @crclass=();
  }
  warn "# read nid=$nid, nclass=$ncl from $crclassf\n" if($debug);
  return($nid,\%crclass,\@crclass);
}  

sub putcovtab {
  my($outpart)= @_;
  warn "# output $outpart\n" if($debug);
  rename($outpart,"$outpart.old") if(-s $outpart);
  open(OUT,">$outpart");
  my $areadlen= ($n_readid<2) ? $readlen : int(0.5 + $sum_readlen/$n_readid); # UPD7f
  #x my $readleno= ($areadlen != $readlen)?"$readlen,avelen:$areadlen" : $readlen;
  #OR make ave 1st .. better for size calcs when readsize varies much
  my $readleno= ($areadlen != $readlen)?"$areadlen,median:$readlen" : $readlen;
  # fixme; n_nomap not in n_readid now; UPD7 now nomap are back in n_readid counts
  my $sst= "#$outpart n_readid $n_readid, n_nomap $n_nomap, n_notcr $n_notcr, n_mapok $n_mapok, n_mapbad $n_mapbad, n_dupbad $n_dupbad, n_mismatch $n_mismatch, n_insert $n_insert, n_delete $n_delete, n_softclip $n_softclip, n_intron $n_intron\n";
  print OUT "#cdsxchr_covtab par: minident=$MIN_IDENT, mindupid=$MIN_DUPIDENT, softclip=$ALLOW_SOFTCLIP, readlen=$readleno, minaln=$MINALN, BIN=$BN, part=$icpu/$ncpu \n";
  print OUT $sst; warn $sst if($debug);
  
  # FIXME: drop cds CovT CovM CovU ; only aCovT aCovM aCovU here???
  my @ocols=qw(ChrID Pos CovT CovM CovU aCovT aCovM aCovU);
  print OUT "#".join("\t",@ocols)."\n";
  # FIXME not %cov == cdscov may be empty

  # FIXME: %crlen may be empty; use crdtype or cmax
  my @cr = keys %crmax;  
  @cr= sort{ $crmax{$b}<=>$crmax{$a} } @cr; 
  
  for my $cr (@cr) { 
    my $crmax= $crmax{$cr}||0; # ?? int( $crlen{$cr} / $BN)
    for(my $i=0; $i<=$crmax; $i++) { 
          
      # vers3h: all these now are fractional? want int of all
      my @cdsv= ($NCDS) ? ( $covt{$cr}[$i], $cov{$cr}[$i], $covu{$cr}[$i]) : (0,0,0);
      my @allv= ($allcovt{$cr}[$i], $allcov{$cr}[$i], $allcovu{$cr}[$i]);
      @ocols= map{ (defined $_) ? int( 0.5 + $_) : 0 } ( @cdsv, @allv );
     
      my $ib=$BN*(1+$i);  # < end of bin, or $ib= ($i==0) ? 1 : $BN * $i;  == start of bin 
      print OUT join("\t",$cr,$ib,@ocols)."\n";   
      } 
    }
  close(OUT); 
}

sub putrdids { 
  my($rdidtab)= @_; # or $outh
  warn "# rdidtab $rdidtab\n" if($debug);
  rename($rdidtab,"$rdidtab.old") if(-s $rdidtab);
  open(OT,">$rdidtab"); 
  #see above:  push @saveid, join("\t", '#ReadID','nmap',@$idclasslist); # add col header for @crclass
  for my $sid (@saveid) { 
    # FIXME for mixed readsets, remove or not SRRprefix. from rdnum here
    # CRTPAT adds 3rd column: rdid, nmap, chrtaglist << chrtaglist must be columns of counts/@$idclasslist
    print OT $sid,"\n";   # push @saveid, "$lid\t$nmap" 
    }
  close(OT);
}

sub putchrtab { 
  my($chrtab)= @_; # or $outh
  # my $chrtab=$ENV{chrtab}||"chrsizeread.tab"; 
  warn "# chrtab $chrtab\n" if($debug);
  rename($chrtab,"$chrtab.old") if(-s $chrtab);
  open(OT,">$chrtab"); 
  my $areadlen= ($n_readid<2) ? $readlen : int(0.5 + $sum_readlen/$n_readid); # UPD7f
  # my $readleno= ($areadlen != $readlen)?"$readlen,avelen:$areadlen" : $readlen;
  #OR make ave 1st .. better for size calcs when readsize varies much
  my $readleno= ($areadlen != $readlen)?"$areadlen,median:$readlen" : $readlen;

  print OT "#cdsxchr_chrtab: nmap,nmult,nuniq = locations, nmread,nuniq,nomap = read counts, readlen=$readleno, part=$icpu/$ncpu\n";
  #x print OT "# n_notcr = $n_notcr\n" if($NOK);
  print OT "#".join("\t",qw(ChrID chrlen nmap nuniq nmult nmread nomap))."\n";

  my($tcw,$trt,$tru,$trd,$trm)= (0) x 9;
  # FIXME: %crlen may be empty; use crdtype or crmax
  # also for icpu parts dont want all of crlen, use crmax == found cr
  my @cr= keys %crmax;  
  @cr= sort{ $crmax{$b}<=>$crmax{$a} } @cr;  
  for my $c (@cr) {
    my $cw= $crlen{$c} || 0; unless($cw){ $cw= $BN* $crmax{$c}; }
    my ($ru,$rd,$rm)=map{ $crdtype{$c}{$_}||0 } qw(uniq mult mread); 
    my $rt=$ru+$rd;
    $tcw+=$cw; $trt+= $rt; $tru+=$ru; $trd+=$rd; $trm+=$rm;
    print OT join("\t",$c,$cw,$rt,$ru,$rd,$rm,0)."\n"; 
  } 

  # add $n_notcr where? similar to n_nomap
  my $otloc= join("\t","# total_locs",$tcw,$trt,$tru,$trd,$trm,$n_nomap) . "\n";  
  print OT $otloc; warn "#$chrtab:$otloc" if($debug);  

  #o: my $tnread=$trm + $n_nomap + $n_notcr; ## need n_notcr if use n_nomap
  # $n_readid + $n_nomap != $tnread .. NOT, $trm is <= n_readid due to mismaps
  #ob: my $tnread= $n_readid + $n_nomap;  
  my $tnread= $n_readid; # UPD7: nomap now back in n_readid
  my $tmm= $trm - $tru;
  my ($pru,$pmm,$prm,$pnomap,$pnotcr)= map{ sprintf "%.2f%%", 100*$_/$tnread; } ($tru,$tmm,$trm,$n_nomap,$n_notcr);
  my $onotcr=""; #x  ($NOK)? "notcr:$n_notcr,$pnotcr" : "";
  my $otread= join("\t","# total_reads",$tnread,"mapt:$trm,$prm","uniq:$tru,$pru","mult:$tmm,$pmm", "nomap:$n_nomap,$pnomap", "readlen=$readlen")."\n";  #,$onotcr
  print OT $otread;  warn "#$chrtab:$otread" if($debug); 
 
   
  #UPD7: add sum stats for cdsreads .. crclasses
  if(UPD7 and $NCDS) { 
    # NOTE: cdsrdclass_in is all IDs mapped for another run, maybe not same input read ids for this map run (chrasm), 
    #  instead count unmapped input cds-class read ids = cdsrdclass_miss
    
    # fixme cdsany: 1,2,3 sums are not exclusive, will have some double/triple counts of same read
    # .. add bit 4 == any rd class: cdshitany = $cdsrdclass_sum[4]; cdsmissall= $cdsrdclass_miss[4] ?
    # my ($cdsany,$cdsanyin,$cdsmiss)=(0,0,0); 
    # for my $i (1,2,3) { 
    #  $cdsany +=  $cdsrdclass_sum[$i]; $cdsmiss += $cdsrdclass_miss[$i];  
    #  $cdsanyin += $cdsrdclass_sum[$i] + $cdsrdclass_miss[$i]; # was $cdsrdclass_in[$i]; 
    #} 
    
    my $cdsanyhit= $cdsrdclass_sum[4]||0; 
    my $cdsallmiss= $cdsrdclass_miss[4]||0; 
    my $cdsno= $cdsrdclass_sum[0]||0; #< no cdsrd class hit ## my $cdsnoin= $cdsrdclass_in[0]; #?
    my $cds123  = join(",",@cdsrdclass_sum[1,2,3]);
    my $cds123miss= join(",",@cdsrdclass_miss[1,2,3]);
    #x my $cds123in= join(",",@cdsrdclass_in[1,2,3]); #<< should this be $cdsrdclass_sum[$i] + $cdsrdclass_miss[$i]
    
    my $otread= join("\t","# cds_reads",$NCDS,"hit/miss:$cdsanyhit/$cdsallmiss","cl123hit/miss:$cds123/$cds123miss", "no_cdsrd:$cdsno")."\n";   
    print OT $otread;  warn "#$chrtab:$otread" if($debug); 
  }
  
  close(OT); 
} 

=item run_sam2covtab.sh cluster script

  #!/bin/bash
  # run_sam2covtab7b.sh
  #PBS -N evg_sam2covtab
  #PBS -l vmem=124gb,nodes=1:ppn=8,walltime=1:00:00
  #PBS -V
  
  ncpu=8
  runapp=./sam2covtab7b.pl
  # relocate to runapp=$evigene/scripts/genoasm/gnodes_sam2covtab.pl
  module load samtools
  
  if [ "X" = "X$ncpu" ]; then ncpu=8; fi
  if [ "X" = "X$datad" ]; then echo "ERR:datad=what?"; exit -1; fi
  if [ "X" = "X$bam" ]; then echo "ERR:bam=what?"; exit -1; fi
  if [ "X" = "X$vers" ]; then  vers=7b; fi
  if [ "X" = "X$opts" ]; then  opts=""; fi
  if [ "X" = "X$minaln" ]; then  
    echo "WARN: not set minaln=$minaln for $ncpu cpu"; 
  else
    opts="$opts -minalign=$minaln"; 
  fi
  if [ "X" != "X$minid" ]; then  opts="$opts -minident=$minid";  fi
  if [ "X" != "X$mindup" ]; then  opts="$opts -mindupid=$mindup";  fi
  if [ "X" != "X$chrtab" ]; then  opts="$opts -chrtab $chrtab"; fi
  if [ "X" != "X$readids" ]; then  opts="$opts -readidtab $readids"; fi
  # NOTE -cdstab old name for -readidtab inputs
  
  ibam=$bam
  name=`basename $ibam .bam`
  outtab=$name.cdschr$vers.covtab
  outchrtab=$name.cdschr$vers.chrtab
  
  cd $datad
  echo "# START `date`"
  if [ ! -x $runapp ]; then echo "ERR: miss $runapp"; exit -1; fi
  
  if [ $ncpu -gt 1 ]; then
    i=0;
    while [ $i -lt $ncpu ]; do {
      echo $runapp -icpu $i -ncpu $ncpu $opts -out $outtab -bam $bam 
      $runapp -icpu $i -ncpu $ncpu $opts -out $outtab -bam $bam  &
      i=$(($i + 1));
    } done
    
    wait
  
    echo $runapp -merge -out $outtab  -bam $bam 
    $runapp -merge -out $outtab  -bam $bam 
    echo /bin/rm $outtab.pt* $outchrtab.pt*
  
  else
  
    echo $runapp $opts -out $outtab -bam $bam 
    $runapp $opts -out $outtab -bam $bam 
  
  fi
  
  wc -l $outtab; head -5 $outtab
  
  echo "# DONE `date`"
  
=cut
