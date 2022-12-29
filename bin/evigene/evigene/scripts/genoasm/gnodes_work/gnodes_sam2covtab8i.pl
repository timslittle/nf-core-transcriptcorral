#!/usr/bin/env perl
# gnodes_sam2covtab8a.pl aka evigene/scripts/genoasm/gnodes_sam2covtab.pl
# from make_cdsxchr_covtab3h.pl, samcopyn5tab.sh

=item UPD21AUG27, vers=8i major gene data change

  -- big-pig genome data exposes memory piggery of cds-gene linkage data,
    mem use blows up w/ very high copy (TE)genes in Susscrofa EVm gene set
      on top of very large, repeat/TE rich chr asm (3 Gb pig =~ human)

  -- two major changes:
    
    a. replace cds.readids table ( 10x bigger for pig than chick)
       with cds.bam reading.  Assumption here is that read id order is same in chr.bam and cds.bam,
       so all cds reads can be gotten by iterative sam view cds.bam along side sam view chr.bam
    
    b. restructure genetab hash table, removing covgenebin{geneid}{chrid}{locbin} hash that can grow hugh w/ 3 Gb genomes
       replace w/ cds
        
        $covgenebin{$gid}{$cr}{$bi}{'covm'} => $covgenebin{$gid}{$gbin}[icovm,icovt,icovu,icovz]
          old: {$cr}{$bi} can become huge for pig genome,
          new: {gid}{gbin} is limited to gene sizes, ie 10 .. 100 bins but for few rare, uniq large genes (TITIN =~ 100k, 1k bins)
          .. this upd lacks cr/cloc info, should keep some of that cds/geneid x chrloc 
        $covgene{gid}{sums} remains same?
        genetab output differs somewhat (col 1,2 = geneid, geneloc) but should match cds.covtab in
          useful way
          
   putGenetab8i: table col changes; besides sumgene/readgene non-locus rows, want some kind of GeneID x ChrID
   possible: GeneID:[abc]bin ChrID CPos covs... where GeneID:[abc]bin are pre-selected representative bins/gene,
    select from cds.covtab, median/25%/75% coverage bins? avoid extremes
   for gene 2-9 dupx, fully dupx will have 2-9 diff chr:loc bins for each [abc]
   for genes skew/partial dups, have variable counts at gene:[abc], 
   for xhigh-copy (TE) genes, need some limit on chr:loc bins, to keep mem use in bounds
          
  -- possible changes:
    .. could test larger binsize (ie 1000 from 100) for large genomes, to reduce mem
    
                  
=item UPD21JUN

  -- add gene ids to readid methods, 
    .. save/collect gene ids from cds x dna,
    .. use cds-geneids-rdids into chr x dna to count depth/copynum of genes via that mapping
  -- add gene covtab output
  -- update/cleanup
    
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

my $VERS='8i'; # 8i: remove cds.readids, add iterative cds.bam x chr.bam read map pairings
# 8e=revise readid meth, try to reduce mem use, misc. ; 8a=UPD21JUN prior
#  '7d'; #large updates 2020.dec-2021.feb; small upd 20nov30, 6c? 7b: upd21jan > read id x chr tag class table

use constant UPD7 => 1;        
use constant UPD7e => 1;   #upd21apr24: higher dupident default (0.98 => 0.999), count 2nd suppl aligns of bwa
use constant UPD7f => 1;   #upd21apr28: for long,lowqual reads; addCigar upd should go to short rd also
use constant UPD21JUN => 1; # 8a,upd21jun14+
use constant UPD21JUN8e => 1; # 8e,upd21jun20+

use constant UPD21AUG => 1; # UPD21AUG20 big-pig memory problems w/ readIds
  # UPD21AUG27 : major big-pig changes: cdsbam read instead of readIds, genecovbins struct change
use constant GENECOVBINS => 5; # 3 == 8c 8b upd for genescov      
use constant GENEBINFIX => 0; # getCdsRead drop 1, 2 after tests
  # GENECOVBINS 4 == reinstate median cov tabulation, xCopy from sums is subject to skew by partial highly duplic reads in cds
  # .. extreme case in drosmel where many CDS genes have this; see in sam2genecov genexcopy outputs, C.nz >> C.M
  # ^^ not sure this is useful yet, solution may be to use C.nz from genexcopy table instead of C.M as gene cn val
  #  for comparison to chr-cds-read copynum vals from sum of cds-reads mapped (covt/covm) =~ C.nz
use constant { SAMFLAG_nomap => 0x04, SAMFLAG_rev => 0x10, SAMFLAG_2ndary => 0x100, SAMFLAG_suppl => 0x800 }; 
  
my $debug=1; # test
my $DEBUG_RDFILT= $ENV{DEBUG_RDFILT}||0; #UPD21AUG: DEBUG mem overflows in readfilter with big-pig data set
my $DEBUG_GENELOC= $ENV{DEBUG_GENELOC}||0; #UPD21AUG: DEBUG mem overflows in readfilter with big-pig data set
my $DEBUG_GBINLIM= $ENV{DEBUG_GBINLIM}||0; #UPD21AUG: DEBUG genetab over-counts
#d warn "#  $0 V=$VERS\n"; for my $db (qw(DEBUG_RDFILT DEBUG_GENELOC DEBUG_GBINLIM)) { warn "# $db=$ENV{$db}\n"; }
use constant TRY_GBINLIM => 4; $DEBUG_GBINLIM= TRY_GBINLIM;
# see TRY_GBINLIM == 5 tests of geneloc x chrloc, now TRY_GENELOC == 2? add %covgenecrloc
use constant TRY_GENELOC => 2; $DEBUG_GENELOC= TRY_GENELOC;

my $MIN_IDENT= 0.40; # UPD7e was 0.65; # lo is best?
my $MIN_DUPIDENT = 0.98; # 7f= 0.999; #UPD7e was 0.98; # was .99/1; hi is best? or 1.0; # == ident equal to 1st/top align, lower if desired

# ^ ?? change to offby=1,2,3 .. dup read align below best align
my $ALLOW_SOFTCLIP = 0; #  
my $SAVE_RDIDS= 0;  
use constant kMAXIDCLASS => 9; # idclass limit, using 3-4 now
my $CRTPAT=''; # no default, this is bad: '^[A-Za-z0-9]+'; #UPD7? TE|CDS|UNK chr class tag for SAVE_RDIDS
my $BN=$ENV{binsize}||100; 
my $SKIPR=$ENV{skipread}||0; # FIXME: skipread=0 means ignore, readmap.bam not paired

use constant NSAMPLE => 20000; ## 500 not enough when variable like Fig.MiSeq
use constant MIN_MINALN => 30; ## require this min align
use constant MAX_MINALNCDS => 120; ## max for shorter cds aligns, w/ intron breaks

my($icpu,$ncpu,$outtab,$inbam,$cdsbam,$inchrtab,$duptab,$topcount,$crclassf,$samcpu,$cdstab,$rdidtab)=(0) x 19;
$rdidtab= undef;
my( $MINALN,$nsample, $readlen, $sum_readlen)= (0) x 9; 
my @sreadlen; # $sminaln, $sreadlen, 
my($n_readid, $n_nomap, $n_notcr, $n_mapbad, $n_dupbad, $n_mapok, $n_partb,
   $n_intron, $n_insert, $n_delete, $n_softclip, $n_mismatch)= (0) x 19; # globals?
my($topaln,$lid,$idp,$idn,$ndi,$ncdsrd)=(0) x 9;
my $MERGE = undef; # defined($MERGE) now means do it
my ($RIDPREFIX,$GIDPREFIX,$ncutGIDPRE,$RDIDIsNum,$RDIDIsNumerr)=(0) x 9;  #UPD7e need to test each read set?
  # $RIDPREFIX obsolete v8i; $RDIDIsNum,$RDIDIsNumerr obsolete v8f
my $LONGR=0; # ($inrdlen>500); # UPD7f FIXME
my $MAXSHORTRD= 500; # UPD7f, switch to LONGR method if rdlen>this
my $outgenetab = undef; my $only_genetab=0; # UPD21JUN option
my $SAVE_GENECOV= 0;  
my $GeneIdsInBam= $ENV{'GeneIdsInBam'} || 0; # UPD21AUG20 ** Drop for cdsbam reading

our( %crmax, # chr max loc found in chr.bam; see also %crlen should be ~same
  %cov,%covt,%covu, # UPD21AUG20 how are these used now?
  %allcov,%allcovt,%allcovu, # primary chr cov data
  @savemap, @saveid, # UPD21AUG20 drop usage
  %covgene, %covgenebin, # UPD21AUG20 change struct
  %covgenecrloc,  # TRY_GENELOC == 2 UPD21SEP04
);

#UPD7: change %cov,covt,covu to single hash of cdste-readid-hit-by-class
# FIXME: 2 uses of chrtab: input for chr, crlen, reads,  differs from output : two opt names?
# .. below outchrpart should not eq chrtab input name

my $optok= GetOptions( 
  'output|covtab=s',\$outtab, 
  'bam=s', \$inbam, 
  'cdsbam=s', \$cdsbam, 
  'chrtab=s',\$inchrtab, #<< change to $inchrtab to avoid fname confusion, outchrtab differs, always $outtab.chrtab
#o:  'cdstab|readidtab=s',\$cdstab, # read table of cds numeric read IDs, num.aligns
#o:  'savereadids:s',\$rdidtab, # write table of cds numeric read IDs, num.aligns
#o:  'GeneIdsInBam!',\$GeneIdsInBam, #UPD21AUG
  'ridprefix=s',\$RIDPREFIX, 'gidprefix=s',\$GIDPREFIX,
  'genetab:s',\$outgenetab, # write table of genes coverage, using cds-readids x chr read map cov
  'onlygenetab!',\$only_genetab,
  'CRTPAT=s',\$CRTPAT, # use w/ savereadids
  'crclassf|idclassf=s',\$crclassf, # alt table chr => class for savereadids
  'icpu=i', \$icpu, 'ncpu=i', \$ncpu,  'samcpu=i',$samcpu,
  'topcount=i', \$topcount, # quick sample, samview bam | top -n 10_000_000
  'minident=s',\$MIN_IDENT, 'mindupident=s', \$MIN_DUPIDENT,  
  'minalign|MINALN=s',\$MINALN,  
  'softclip!',\$ALLOW_SOFTCLIP,
  'longreads!',\$LONGR,
  'binsize=i',\$BN, 'skipread=i',\$SKIPR, 
  'merge:s',\$MERGE, # UPD7, flag to merge ncpu part tables : 'merge!',
  'debug!', \$debug, 
  );

my $optinok= (($inbam and -s $inbam) or (defined $MERGE and ($outtab or @ARGV)) );

die "usage: sam2covtab -bam name.chr_reads.bam -cdsbam name.cdste.bam -out name.v$VERS.covtab
  requires 'samtools view inbam'
  opts: -minident=$MIN_IDENT -mindupident=$MIN_DUPIDENT -minalign=$MINALN 
  -ncpu $ncpu -icpu 0..3  : subsets for parallel process, then -merge subsets
  -binsize=$BN  -skipread=2|1 (which of paired reads)  
  -idclass genes.classtab :  
sam2covtab -merge -out name  : merge all table parts from icpu 0..3,  
sam2covtab -merge=chrtab -out name.chrtab part1.chrtab part2.chrtab  : merge one table parts from ARGS
"  unless($optok and $optinok);

# readChrtab() globals
my($NOK,$haslen,%crok,%crlen,%crdtype)=(0,0,0); #8e.drop: $hasread,,%cread 
  
# global readid data
use constant USE_RIDp => 1; #  DROP; problem if this is big memory pig
my ($NCDS,$hasCdsBam,$useRIDp)=(0,0,USE_RIDp);  
my %cdslen=(); # from cds.bam header, or table of id/len ??
#obs: my %cdsrdid=(); # for useRIDp == $cdsrdid{$idp}[$idnum]
#obs: my %generdid=(); #UPD21JUN: mem overflow problem with huge drospse SRR read set (CDepth=160 )
my %redata=();   #UPD21JUN8e: slim mem: change readid hash to one only: redata{idp}{inum}  = cdsrdid{idp}[inum] + generdid{idp}{inum}
my @cdsrdclass_in=(0) x 5; my @cdsrdclass_sum=(0) x 5; my @cdsrdclass_miss=(0) x 5; 
  
my($nidclass,$idclassh,$idclasslist)=(0); # read_idclass globalse; renamed crclass
my($MINALNCDS, $lcdsreadid, $cdsiid, $lastcdsrow) = (0) x 4; # getCdsRead globals, was our()
$MINALNCDS= $MINALN;

sub cleangid{ my($g)=@_; if($GIDPREFIX and $g =~ s/$GIDPREFIX//){ $ncutGIDPRE++; } $g =~ s/:/_/g; return $g; }

sub MAIN_stub {}

  $ncpu||=1;  $icpu=0 if($ncpu==1);
  $MIN_IDENT /= 100 if($MIN_IDENT > 1);
  $MIN_DUPIDENT /= 100 if($MIN_DUPIDENT > 1);
  
  my $iname=$inbam || $outtab; $iname =~ s/\.(bam|covtab|chrtab|genetab|readids)//;
  $outtab= "$iname.v$VERS.covtab" unless($outtab); #? $iname.$VERS.covtab
  my $outchrtab= $outtab;  $outchrtab=~s/\.\w+$//; $outchrtab.=".chrtab";

  # no default for cdstab,input readidtab
#   if(defined($rdidtab)) { 
#     $SAVE_RDIDS= 1;  $rdidtab= "$iname.readids" unless($rdidtab =~ /\w\w/); 
#   }
 
   # UPD21JUN: outgenetab Option with input cdstab|readidtab
  if(defined $outgenetab) {
    $SAVE_GENECOV= ($SAVE_RDIDS)? 0 : 1; #not both
    unless($outgenetab =~ m/\w/) { 
      $outgenetab= $outtab; $outgenetab=~s/\.\w+$//; $outgenetab.=".genetab"; 
    } 
  }
  
  my $outpart= $outtab;
  my $outchrpart= $outchrtab; if($outchrpart eq $inchrtab) { $outchrpart.="out"; } # name bug: .chrtabout with -merge
  my $outgenepart= $outgenetab; # may be undef
  
  if(defined $MERGE) { 
    # -merge : merge all table type, or -merge=covtab,chrtab,readids,.. which type only
    unless($MERGE=~/\w/){ $MERGE=($only_genetab)?"genetab":"covtab,chrtab,genetab,readids"; }# all types
    
    mergeparts("covtab", $outtab, @ARGV) if($MERGE=~/cov/); 
    mergeparts("chrtab", $outchrtab, @ARGV) if($MERGE=~/chr/); # -merge assume chrtab=output name, was outchrpart => .chrtabout bug
    mergeparts("genetab", $outgenetab, @ARGV) if($MERGE=~/gene/);   
    #o: mergeparts("readids", $rdidtab || "$iname.readids", @ARGV) if($MERGE=~/readid/);    
    exit;
  }
    
  warn "#sam2covtab par: minident=$MIN_IDENT, mindupid=$MIN_DUPIDENT, softclip=$ALLOW_SOFTCLIP, BIN=$BN, part=$icpu/$ncpu \n"
    if($debug);
  if($ncpu>1 and $icpu>=$ncpu) { die "# err: icpu >= ncpu : -i $icpu, -ncpu $ncpu"; }

=item UPD21AUG20: readReadIds() now is BIG-MEM-PIG

  #UPD21AUG20: readReadIds() now is BIG-MEM-PIG with eg. pig genome data, 5/10 icpu failing in 240 GB mem
  # need to rewrite to move cds-te-readid info into chr.bam rows, could also slim down readid hash w/ int ids
  # .. move this call AFTER read_idclass, readChrtab  (failing samtools view -H there w/ big-pig)
  # memuse: perl module needs special compile: use Devel::Peak; mstat("marker1");'
  # perl -e 'use Devel::Peek; mstat("one");'
  #    one: perl not compiled with MYMALLOC
  # syscall: ps -e -o pid,pcpu,pmem,stime,etime,command --sort=-pmem | grep $progname | head
  # PID %CPU %MEM STIME     ELAPSED COMMAND
  # 4701  5.2  0.7 Aug08 13-02:07:37 /usr/lpp/mmfs/bin/mmfsd
  # 15521  0.0  0.2 Aug19  1-19:12:39 python

=item GeneIdsInBam

  GeneIdsInBam replaces readReadIDs mem pig .. 
  
  gnodes_samaddrdid.pl -readids cds.readids -bam chrasm.bam -outbam chrasm.rdid.bam
     rdid .... seq qual  eC:i:1  eG:X:geneid1,geneid2 
     
  samtools view -F 260  chrpig11c_SRR4341337_b2_bwa.rdid.bam | head
  
  SRR4341337.12	0	chr1	10047128	0	100M	*	0	0	s	*	NM:i:0	XS:i:100
  SRR4341337.17	16	chr15	32807967	0	65S33M2S	*	0	0	s	*	NM:i:0	XS:i:33
  SRR4341337.20	0	chr13	104534046	0	47M53S	*	0	0	s	*	eC:i:1	eG:Z:Susscr4EVm009461t1	NM:i:1 XS:i:42
  SRR4341337.21	16	chr17	50962877	0	46S31M23S	*	0	0	s	*	NM:i:0	XS:i:31
  SRR4341337.22	0	chr1	93124355	0	100M	*	0	0	s	*	eC:i:1	eG:Z:Susscr4EVm011357t1,Susscr4EVm016404t1,Susscr4EVm058735t1	NM:i:1	XS:i:95
  SRR4341337.26	0	chr6	43865301	0	25M2I73M	*	0	0	s	*	eC:i:1	eG:Z:Susscr4EVm050552t1 NM:i:10	XS:i:50

=cut

  #obsolete# ($NCDS)= readReadIds($cdstab); # ret=$NCDS
  
  ($haslen)= readChrtab($inchrtab, $inbam); # NOTE: this reads chr-len from inbam headers @SQ ID: LEN: unless inchrtab has it

  #UPD7b: always call read_idclass, check kMAXIDCLASS and sam hdr IDs ~ $CRTPAT 
  ($nidclass,$idclassh,$idclasslist)= read_idclass($crclassf,0,\%crlen); # cds,te class by id, $idclass->{id} = class
  
  # BUG if haslen==0 but ncpu,icpu>1 .. reset one, need part tags
  if(($ncpu>1 and $icpu<$ncpu)) { # ($hasread or $haslen) and 
  
    ## ?? crok{} no longer used, icpu makes parts from readid index, icpu == iread % ncpu
    if($haslen) {
      # if($hasread>$ncpu) { @cr=sort{$cread{$b}<=>$cread{$a} or $a cmp $b} keys %cread; }  else { }
      my @cr=sort{$crlen{$b}<=>$crlen{$a} or $a cmp $b} keys %crlen; 
      for my $i (0..$#cr){ if($icpu == ($i % $ncpu)){ $crok{$cr[$i]}=1; $NOK++; } }
    }
    
    $outpart="$outtab.pt$icpu";
    $outchrpart= "$outchrpart.pt$icpu";
    if($outgenepart) { $outgenepart= "$outgenepart.pt$icpu"; } # outgenepart may be undef
    $rdidtab .= ".pt$icpu" if($SAVE_RDIDS and $rdidtab); # rdidtabpart ?
    warn "# output part $icpu/$ncpu nchr=$NOK to $outpart\n" if($debug);
  } 

  #--------- moved from above ------
  #UPD21AUG27: readReadIds obsolete
  #UPD21AUG20: readReadIds() now is BIG-MEM-PIG with eg. pig genome data, 5/10 icpu failing in 240 GB mem
  # need to rewrite to move cds-te-readid info into chr.bam rows, could also slim down readid hash w/ int ids
  # .. move this call AFTER read_idclass, readChrtab  (failing samtools view -H there w/ big-pig)

  # if($cdstab and not $GeneIdsInBam) {
  #   my $gidinbam= "$cdstab.inbam";  #FIXUP: look for file flag GeneIdsInBam
  #   if( -f $gidinbam ) { $GeneIdsInBam=1; }
  #   else { ($NCDS)= readReadIds($cdstab); }
  # }
  
  #----------------
  
  ## inbam must be read-ordered, ie all multimaps together, only 1 of paired reads counted now

  #? other samt filter opts: -q minqual
  my $STMINFILT=0;
  my $stopt = ($SKIPR==0)? "" : ($SKIPR==1) ? " -F 0x40" : " -F 0x80"; # drop SKIPR option?
  if($MINALN > 0) { $stopt .= " -m $MINALN"; $STMINFILT=1; } #upd7
  if($samcpu>1){ $stopt.=" --threads $samcpu"; }
  
  warn "# samtools view $stopt $inbam\n" if($debug);
  if($topcount != 0) {
    my $pheadtail= ($topcount < 0) ? "tail -n$topcount" : "head -n$topcount";
    open(IN,"samtools view $stopt $inbam | $pheadtail |") or die "ERR: samtools view $stopt $inbam";
  } else {
    open(IN,"samtools view $stopt $inbam |") or die "ERR: samtools view $stopt $inbam";
  }
  
  $hasCdsBam=0;
  if($cdsbam and -s $cdsbam) { # UPD21AUG27: minalign/MINALN different for cds
    $hasCdsBam= 1;
    open(CDSBAM,"samtools view -H $cdsbam |");     # read hdr geneids, len
    while(<CDSBAM>) {
      if(/^\@SQ/){
        my($cr)= m/SN:(\S+)/; my($crlen)= (m/LN:(\d+)/)?$1:0;
        if($cr and $crlen > 0){ $cr= cleangid($cr); $cdslen{$cr}= $crlen; $NCDS++; }
      }
    } close(CDSBAM);
    
    $stopt =~ s/ \-m $MINALN// if($STMINFILT);
    warn "# samtools view $stopt $cdsbam\n" if($debug);
    open(CDSBAM,"samtools view $stopt $cdsbam |") or die "ERR: samtools view $stopt $cdsbam";
  }
  
  my ($ntmap)= readFilter(*IN,$STMINFILT);  #? add *CDSBAM
  close(IN); if($hasCdsBam != 0) { close(CDSBAM); }
  
  putcovtab( $outpart);
  putchrtab( $outchrpart); 

  putGenetab8i($outgenepart) if($SAVE_GENECOV); # GENECOVBINS==5, simplified, merge-able UPD21AUG
#  if( GENECOVBINS >= 5) {
#   putGenetab8i($outgenepart) if($SAVE_GENECOV); # simplified, merge-able 
#   } elsif( GENECOVBINS == 4) {
#   putGenetab8e($outgenepart) if($SAVE_GENECOV); # simplified, merge-able UPD21JUN8e
#   } elsif( GENECOVBINS == 3) {
#   putGenetab8c($outgenepart) if($SAVE_GENECOV); # simplified, merge-able
#   } elsif( GENECOVBINS == 2 ) {
#   putGenetab8b($outgenepart) if($SAVE_GENECOV); # tab8b mergable, tab8a not merge-able 
#   } elsif(  GENECOVBINS == 1 ) {
#   putGenetab8a($outgenepart) if($SAVE_GENECOV); #  UPD21JU, may be undef (GENECOVBINS == 1)
#   }
  
  #obs: putrdids($rdidtab) if($SAVE_RDIDS); 

# END
#-------------------------

## FIXME: readIdNum() for new ids=SRRxxx.nnn and .nnn_2 or .nnn/2
## .. RDIDIsNum = 3 for idpat=prefix.nnn[_/\D]p
## Also cant expect type 3 at first,  RDIDIsNum not stable..

# sub readIdNum { # UPD8a UPD21AUG: Obsolete
#   my($rdid)=@_;
#   my($idp,$idn,$idx)=(0,0,0);
#   if($RDIDIsNum == 1) { $idn= $rdid; } 
#   elsif($rdid =~ /^\d/) { $idn=$rdid; $RDIDIsNum=1; }
#   elsif($rdid =~ m/^(.+)(\d+)\D(\d+)$/) {
#     ($idp,$idn,$idx)=($1,$2,$3); $idp.="p$idx";  # dont set RDIDIsNum
#   } elsif($rdid =~ m/^(.+)(\d+)$/) {
#     ($idp,$idn)=($1,$2);   # dont set RDIDIsNum
#   } else {
#     $RDIDIsNumerr++; ($idp,$idn)= $rdid =~ m/(.+)(\d+)/?($1,$2):($rdid,1);   
#   }
#   $idn= int($idn); #ensure, eg. srr123.456.sra style
#   return($idp,$idn);   
# }

# sub readReadIds { 
#   my($cdstab)= @_;
#   my($ok,$inh);
#   $NCDS=0;
#   #UPD7: revise here to read multi-col read id table: rdid, class1, class2, class3; classes are 0/1 flag
#   if($cdstab and ($ok,$inh)= openRead($cdstab) and $ok) { 
#     my ($isnum,$lastidp,$err,$hasgeneids)=(-1,"",0,0); 
#     while(<$inh>){ 
#       if(/^\W/){
#         if(/^#ReadID/){ my($lrd, @classlabel)=split; } # use @classlabel for output?
#         next;
#       }
# 
#       #o: my($rdid,@crclass)=split; 
#       # UPD21JUN: include gene-IDs column for readids
#       my @crclass=split; 
#       my ($rdid,$geneids)=(0,0);
#       if($hasgeneids == 0) { # look for col, probly last
#         for my $i (1..$#crclass) { 
#           if($crclass[$i] =~ /^\d/) { next; } # crclass cols are numeric
#           elsif($crclass[$i] =~ /^\w/) { $hasgeneids=$i; last; } # presume geneids like AT123, EVm123, XM_123
#         }
#         if($hasgeneids == 0){ $hasgeneids= -1; } # stop looking
#       } 
#       
#       if($hasgeneids>0) {
#         ($geneids)= splice(@crclass,$hasgeneids,1);
#         if($GIDPREFIX){ $geneids=~s/$GIDPREFIX//g; } #UPD21AUG: strip GIDPREFIX like RIDPREFIX to save mem
#       }
#       ($rdid)= shift @crclass; # always 1st col
#       if($RIDPREFIX){ $rdid=~s/$RIDPREFIX//; }
#       
#       my $nd=$crclass[0]||1; 
#       my $crclass=$nd;
#       if( @crclass>1 ) { 
#         $crclass=0; 
#         $crclass |= 1 if($crclass[1]>0);  
#         $crclass |= 2 if($crclass[2]>0);
#         $crclass |= 4 if($crclass[3]>0);
#         #? $crclass |= 8 if($crclass[4]>0); # how many cols?        
#         if(0) { #UNUSED  cdsrdclass_in
#         $cdsrdclass_in[1]++ if($crclass & 1);
#         $cdsrdclass_in[2]++ if($crclass & 2);
#         $cdsrdclass_in[3]++ if($crclass & 4);
#         }
#       }
# 
#       # my($idp,$idn)= readIdNum($rdid);
#       # if($lastidp and $idp ne $lastidp){ $useRIDp++; } #? dont use soft switch on hash/array
#       
#       ## change to int-only read ids?? save mem, needs more work to convert to ints
#       ## or would this: $redata{$rdid}=xxx be smaller than that: $redata{$idp}{$idn}=xxx ?
#       $redata{$rdid}= "$crclass\t$geneids"; # UPD21AUG
#       
#       # if(UPD21JUN8e) {
#       #   $redata{$idp}{$idn}= "$crclass\t$geneids";
#       # } elsif(USE_RIDp) { 
#       #   if(UPD7e){ if($cdsrdid{$idp}[$idn]){ $idp.="b"; } } # FIX for sam cat of _1/_2 reads same id
#       #   $cdsrdid{$idp}[$idn]=$crclass; 
#       #   $generdid{$idp}{$idn}=$geneids if($geneids);  #UPD21Jun21: changed [idn] to {idn} for sparse geneids, mem oflow
#       # } 
#       # #old: else { $cdsrdid[$idn]=$crclass; }
#       
#       $NCDS++; # $lastidp=$idp;
#     } close($inh);
#     warn "# ncds readids=$NCDS from $cdstab \n" if $debug; 
#   } 
# 
#   return($NCDS);
# }

sub readChrtab {
  my($inchrtab, $inbam)=@_;
  my($ok,$inh);
  if($inchrtab and ($ok,$inh)= openRead($inchrtab) and $ok) { #skip inchrtab and use only sam -H sizes ??
    while(<$inh>){ next if(/^\W/); 
      my($cr,$crlen,$nread)=split;  #? UPD7 add crclass here? but chrtab now is output read counts table, no class info
      #drop 8e: $cread{$cr}= $nread||0; $hasread++ if($nread>0);
      $crlen{$cr}= $crlen; $haslen++ if($crlen>0);  # may be only list of chr, crlen w/o nread
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
  return($haslen);
}


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
  if(@pt<1){ warn "Warn: no parts to merge for $ofile\n"; return 0; } # need @pt > 1 ?
  rename($ofile, "$ofile.old") if(-s $ofile);
  
  # is genetab like readids? ie no cb bin location column?
  my $mtype= ($mergeflag and $mergeflag =~ /readids/)? 1 : 2;
  my $dolocgene=0; # only if DEBUG_GENELOC ?
  if( GENECOVBINS >= 5) {
     $mtype=2 if($mergeflag =~ /gene/); # new genetab 2 fixed cols 0,1,2 = GeneID, GeneBin
     $dolocgene=1 if($mergeflag =~ /gene/); # 5: genetab locgene needs special merge
  } elsif( GENECOVBINS == 4) {
     $mtype=3 if($mergeflag =~ /gene/); # genetab 3 fixed cols 0,1,2 = GeneID, ChrID, ChrBin
  } elsif( GENECOVBINS == 3) { # GENECOVBINS == 4 may add back fixed col3 ChrBin
     $mtype=2 if($mergeflag =~ /gene/); # genetab 2 fixed cols 0,1 = GeneID, class    
  } elsif( GENECOVBINS == 2) {
     $mtype=3 if($mergeflag =~ /gene/); # genetab 3 fixed cols 0,1,2 = GeneID, ChrID, ChrBin
  }
  
  #x my $mgeneids= (defined $outgenetab and $SAVE_RDIDS)?1:0; # only for cds.readids, not cds.covtab** GENECOVBINS >= 1
  my $mgeneids= 0; #UPD21AUG.drop: (defined $outgenetab and $SAVE_RDIDS and $mergeflag =~ /readids/)?1:0; # GENECOVBINS >= 1

  our(%cloc);
  sub get_locgene { # dolocgene
    my($id)=@_; our(%cloc); my @oloc; 
    use constant MAX_locgene => 9;
    return(0) unless(exists $cloc{$id});
    my @loc=sort{ $cloc{$id}{$b} <=> $cloc{$id}{$a} } keys %{$cloc{$id}}; 
    my $nloc=@loc; 
    for(my $i=0; $i<$nloc; $i++) { 
      my $lc=$loc[$i]; my $nr=$cloc{$id}{$lc} || 0; 
      last if($i >= MAX_locgene or $nr <= 9); # nr <= MIN_READS
      push @oloc, "$lc\t$nr\t0"; # print "$id\tlocgene\t$lc\t$nr\t0\n"; 
    } 
    return($nloc, @oloc);
  }

     
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
        if(/^#(ChrID|ReadID|GeneID)/) { push @thead,$_ unless(@thead); }
        elsif(/^#/ and /n_mapok|total_locs|total_reads|cds_reads|minident=/){ s/#/#$pti./; push @tcomm,$_; }
        next;
      }

      my($id, $cb, @vals);
      if($mtype == 3) { my $cid; ($id,$cid,$cb,@vals)= split; $id.="\t$cid"; }
      elsif($mtype == 1) { ($id,@vals)= split; $cb=1; } 
      else { ($id, $cb, @vals)=split; }
      
#      ## FIXME merge=readids now has GeneID col to merge, not += additive
#       if($mgeneids) {
#         my $ie=$#vals;
#         my $gids= pop(@vals); # always last col?
#         $tab{$id}{$cb}[$ie] .= "$gids," if($gids =~ /\w\w/);
#       }

      if($dolocgene and $cb =~ m/locgene/){ my($cr,$cl,$nr)=@vals; $cloc{$id}{"$cr\t$cl"}+=$nr; } 
      elsif($dolocgene and $vals[0] =~ m/^cr:/){ 
        my $crcb= shift @vals; $cb .= "\t$crcb"; # TRY_GENELOC == 2
        for my $i (0..$#vals) { $tab{$id}{$cb}[$i] += $vals[$i]; }
        } 
      else { for my $i (0..$#vals) { $tab{$id}{$cb}[$i] += $vals[$i]; $ncol=$i if($i>$ncol);  } }
      $ntin++;
    } close(IN);
  }
  
  # GENECOVBINS == 4, -merge=genetab, want median of chr,pos vals .. calc here?

  open(OUT,">$ofile"); 
  print OUT @thead if(@thead);
  for my $id (sort keys %tab) { 
    # sort genetab col2: readgene,sumgene at top before 100,200,.. ok?

    for my $cb (sort{ $a <=> $b } keys %{$tab{$id}} ) {
      my @vals= @{$tab{$id}{$cb}}; 
#       if($mgeneids) {
#         my $gids= $vals[-1]; # always last col? 
#         if($gids =~ /\w\w/){ my %gids= map{ $_ => 1 } split ",",$gids;  
#           $gids=join",",sort keys %gids; $vals[-1]=$gids; #x push @vals,$gids; 
#           }
#       }
      if($mtype == 1) { print OUT join("\t",$id,@vals)."\n"; } # no cb for mtype == 1
      else { print OUT join("\t",$id,$cb,@vals)."\n"; } # ok for mtype=3 where id == gid\tcrid
      $nrow++; 
      }

    if($dolocgene) {
      my($nl,@oloc)= get_locgene($id);
      if($nl){ for my $oloc (@oloc) { print OUT join("\t",$id,"locgene",$oloc)."\n"; } }
    }

  }
  print OUT @tcomm if(@tcomm);
  close(OUT);
  
  warn "#sam2covtab merge: nparts=$npt, nrows=$nrow, ncols=$ncol, ntin=$ntin to $ofile\n";
  return($nrow);
}


# ( $cdsread, $geneids, $genebins)= getCdsRead_UPD( $rdid, $iid);
sub getCdsRead {  # == getCdsRead_UPD
  my( $rdid, $iid)= @_; # , $crfl,$cr,$cb from chr.bam row cols[1..4]
  my( $cdsread, $geneids, $genebins)=(0) x 9;
  my( $moredata,$cdsrow)=($hasCdsBam > 0,''); 
  
  return ( $cdsread, $geneids, $genebins) unless($moredata);
  if($lastcdsrow) { $cdsrow= $lastcdsrow; $lastcdsrow=''; $lcdsreadid=0; $cdsiid--; }
  else { $cdsrow= <CDSBAM>; unless($cdsrow){ $moredata=0; $hasCdsBam= -1;  } } 
  
  my @thisrows=();
  while($moredata) {
    # my($drdid,$dfl,$dgid,$dgb,$dgx,$dcig,@dsamx)= split" ",$cdsrow;    
    my($drdid,$dfl)= split" ",$cdsrow;    
    
    #?speedup? if($dfl >= 256) { $cdsrow= <CDSBAM>; unless($cdsrow){ $moredata=0; $hasCdsBam= -1; last; } next; }
    my $firtsmap= ($dfl < 256);
    my $thisid= ($drdid eq $rdid)?1:0;
    $cdsiid++ if($firtsmap); # always 1st map here?
    
    if($thisid) {
      $moredata=0; 
      my $badmap=($dfl & 0x4)?1:0;
      push @thisrows, $cdsrow unless($badmap);
      $cdsrow='';
      while(<CDSBAM>){ # collect all dupmaps here, read to next rid even if badmap/nodups
        my($erdid)= split" ",$_;                 
        if($erdid ne $drdid) { $lastcdsrow=$_;  $cdsiid++; $moredata=0; last; } 
        else { push @thisrows, $_; } # never for dup:($efl & 0x4)
      }
    } elsif( $cdsiid > $iid) { 
      $lastcdsrow= $cdsrow; $cdsrow=''; $moredata=0;   
    } else {
      $cdsrow= <CDSBAM>; unless($cdsrow){ $moredata=0; $hasCdsBam= -1;  }
    }
  }

  if(@thisrows) {
    $cdsread=1; 
    my $firstm= 1; my(@gbins,%gids);
    my $minalncds= $MINALNCDS || MIN_MINALN; # $MINALN || int( $readlen * $MIN_IDENT);  
    for my $ro (@thisrows) {
      my($grdid,$gfl,$gid,$gb,$gx,$gcig)= split" ",$ro;
      my $alen=0; while($gcig =~ m/(\d+)M/g) { $alen += $1; }
      my $cend=$alen; while($gcig =~ m/(\d+)[DN]/g) { $cend += $1; } # add to gend bin
      if($ro =~ m/NM:i:(\d+)/){ $alen -= $1; }
      
      last if($alen < $minalncds); # use GENEMIN_LAST way
      # use constant GENEMIN_LAST => 1; 
      # test1:aweed_test8ipbwa vs inbwa_genesum, no strong effect for next vs last; next slightly slower, may have *fewer* dups
      # FIXME: last/next needs test, a. last drops dup cds cov? b. next slower?
      # if($alen < $minalncds){ if(GENEMIN_LAST) { last; } else { next; } } # sam mult aligns *should* be ordered by +align 
      
      $gid= cleangid($gid); 
      $gids{$gid}++; 
      
      #FIXME data bug from too much genebins, see addCigarDepth
      # 1a .. try genespan: gid:gb:ge base-level, let addCigar figure cover span? not tested, skip
      # 2b .. try 1 bin only? or n-full bins, n= int(0.9+aln/BN)
      #
      # use constant GENEBINFIX => 0; # 1, 2, 0
      # if(GENEBINFIX == 1) {  # drop this way
      #   my $ge= $gb+$cend;
      #   push @gbins, "$gid:$gb:$ge";
      # } elsif(GENEBINFIX == 2) {   #? use this way, or orig/0
      #   my $nbin = int(0.49 + $alen/$BN); 
      #   my $gbin = int($gb/$BN); 
      #   for(my $i=0; $i<$nbin; $i++) { my $ib=$gbin + $i; push @gbins, "$gid:$ib"; } 
      #   #bad: for(my $i=0; $i<$nbin; $i++) { my $ib=$gbin + $i*$BN; push @gbins, "$gid:$ib"; } # *BN is wrong !!
      # } else  # use orig
      
      # UPD21SEP13 ??? $gb-$cend when SAM REV flag  $gfl & 0x10 .. no, think not, it is rseq that is revcompl to align;
      #  gb..gb+cend is span on chr/gene seq that revcomp(rdseq) aligns to, so target seq is always fwd
      my $gbin = int($gb/$BN); 
      my $gebin= ($cend < $BN) ? $gbin : int(($gb+$cend)/$BN);   # alen should be cend from cigar parse
      for(my $ib=$gbin; $ib<=$gebin; $ib++) { push @gbins, "$gid:$ib"; }
      
      if($firstm) { $minalncds= $MIN_DUPIDENT * $alen; $firstm=0; } 
    }
    
    $geneids= join",",sort keys %gids;
    $genebins= join",",@gbins; #?? SEE addCigarDepth: data bug does this have dup entries? should be only 1 of each gid:ib per read maps
    #? my %gbins= map{ $_,1 } @gbins; $genebins= join",", sort keys %gbins;
  }
  
  return($cdsread,$geneids, $genebins);
}


#FIXME: read <CDSBAM> is hanging or slowing greatly... bad algo, debug
# sub  getCdsRead_BAD { #UPD21AUG27 .. match rdid from cdsbam, other ways?
#   my( $crdid, $criid, $crfl,$cr,$cb)= @_; # from chr.bam row cols[1..4]
#   
#   # this needs some readFilter() tests, ie min-align, min-dup-align, for cds, same as chr? 
#   my($cdsread,$geneids, $genebins)=(0) x 9;
#   return ($cdsread, $geneids, $genebins) unless($hasCdsBam > 0);
#   
#   #global my >> our($MINALNCDS, $lcdsreadid, $cdsiid, $lastcdsrow);
#   my ($thisid,$firstmap,$minaln,$minduplen,$ireadpart,$iskip)= (0) x 9;
#   my $minalncds =  $MINALN || int( $readlen * $MIN_IDENT);  # readFilter() sets MINALN,$readlen
#   $minalncds = 30 if($minalncds<30);
# 
#   our (@gid,@gloc);
#   sub retCdsRead {
#     our (@gid,@gloc);
#     if(@gid) {  
#       $cdsread||=1;
#       my %gi= map{ $_,1 } @gid;  $geneids= join",", sort keys %gi;
#       my %gb= map{ $_,1 } @gloc; $genebins= join",", sort keys %gb;
#       # $geneids= join",", @gid; #? filter out dupls %gid=map{ $_,1 } @gid; $geneids=join",",sort keys %gid; ? order
#       # $genebins=join",", @gloc; 
#     }
#     return($cdsread,$geneids, $genebins);
#   }
#     
#   my $moredata=($hasCdsBam > 0)?1:0;
#   while($moredata) { # <CDSBAM>
#     my $cdsrow;
#     if($lastcdsrow) { $cdsrow= $lastcdsrow; $lastcdsrow=''; $lcdsreadid=0; $cdsiid--; }
#     else { $cdsrow= <CDSBAM>; unless($cdsrow){ $moredata=0; $hasCdsBam= -1; return retCdsRead(); } }
#     
#     my($rdid,$fl,$gid,$gb,$gx,$cig,@samx)= split" ",$cdsrow; 
#     my $badmap= ($fl & 0x4)?1:0;
#        
#     if($rdid ne $lcdsreadid) {
# 
#       $thisid= ($rdid eq $crdid)?1:0;
# 
#       $lcdsreadid=$rdid;      
#       $firstmap=1; # == topaln below
#       $ireadpart= ($ncpu > 1) ? $cdsiid % $ncpu : 0; # is this right for both cds.bam, chr.bam?
#       # $iskip= ($ireadpart == $icpu)?0:1; #? maybe not, want thisid and chr-iskip?
#       # pass chr.bam iskip? or ireadpart= criid % ncpu
#       
#       $cdsiid++; # this should be kept == chr.iid, expect/require same num/order of readids in both
#       if($DEBUG_RDFILT) {
#         warn "#DERD: cds.bam cdsiid=$cdsiid, criid=$criid, cdsrd=$rdid, crrd=$crdid, icpu=$icpu\n"
#             if($criid % 100_000 == 1);  
#       }
#         
#       if($thisid) {
#         # ok collect all mappings; @gid = @gloc = ();
#       } elsif($cdsiid < $criid) {  
#         # read more
#       } elsif($cdsiid > $criid) { 
#         # stop, return any result, but preserve last @samx read for next call
#         # this is only stop/return point until end<CDSBAM> ?
#         $lastcdsrow= $cdsrow;  # BUT could infloop if criid doesn't inc
#         return retCdsRead();
#       } else {
#         # read all mappings of this rdid, then return gids, genebins to caller
#       }
# 
#       # if($ireadpart == $icpu) { # always count new ids, but for nomap, or $n_readid++ unless($fl & 0xF00); # all 2nd maps
#       # }
#       
#     } else {
#       $firstmap=0;
#     }
#     
#     if($thisid and not $badmap) # fixme: not_thisread=(rdid ne crdid) and badmap=($fl & 0x4 or alen < minaln) differ
#     {
#       my $okmap= 1; 
#       my $alen=0; while($cig =~ m/(\d+)M/g) { $alen += $1; }
#       if($firstmap) {  $okmap= 0 if($alen < $minalncds); $minduplen= $MIN_DUPIDENT * $alen; }
#       else { $okmap=0 if($alen < $minduplen); }
#       if($okmap) {
#         $cdsread=1; # FIXME this is bit flag for CDS|TE|UNK data type, now only CDS ;
#         my $gbin = int($gb/$BN);  # end of bin,
#         my $gebin= int(($gb+$alen)/$BN);   # alen should be cend from cigar parse
#         if($GIDPREFIX){ $gid=~s/$GIDPREFIX//; }
#         push @gid, $gid; 
#         push @gloc, "$gid:$gbin"; # for i ($gbin .. $gebin) ???
#         # for(my $ib=$gbin; $ib<=$bebin; $ib++) { push @gloc, "$gid:$ib"; }
#       }
#     }
#     
#   } # moredata
#   
#   return retCdsRead();
# }


sub readFilter {  
  my($inhand,$STMINFILT)=@_;
  # inhand == open(IN,"samtools view $inbam |") 
  my($lid, $iid, $ireadpart, $inrdlen, $nokid,$ntmap,$nmap,$no_cdsread,$crtlist)=(0) x 19;
  my @saver;
  my($cdsread,$geneids, $genebins)=(0,0,0); # UPD21AUG, move cdsread/geneids to bam rdid rows as tags: eC:i:1, eG:Z:gid1,gid2
  if($MINALNCDS < MIN_MINALN){ $MINALNCDS= MIN_MINALN; } elsif($MINALNCDS > MAX_MINALNCDS) { $MINALNCDS= MAX_MINALNCDS; } 
  #not this, needs calc below: $MINALN= MIN_MINALN if($MINALN < MIN_MINALN); 

#   if($SAVE_RDIDS){
#     my @ch=@$idclasslist; push @ch,'GeneID' if($outgenetab);
#     push @saveid, join("\t", '#ReadID','nmap',@ch); # add col header for @crclass
#   }

  my $iskip= 0;
  while(<$inhand>) {
    my @samx= split; 
    my($rdid,$fl,$cr,$cb,$cx,$cig)=@samx;
    
    if($rdid ne $lid) {
    
      if($lid and $nokid > 0) {
        if($LONGR){  
        ($nmap,$crtlist)=  putLongmap( $lid, $inrdlen, $nokid, \@saver, $cdsread,$geneids, $genebins); 
        } else {
        ($nmap,$crtlist)=  putShortmap( $lid, $inrdlen, $nokid, \@saver, $cdsread,$geneids, $genebins); 
        }
        $ntmap += $nmap;
        #o: if($SAVE_RDIDS) { push @saveid, "$lid\t$nmap\t$crtlist"; }  
      }
      
      $lid=$rdid; $nokid=0; @saver=();      
      $ireadpart= ($ncpu > 1) ? $iid % $ncpu : 0;
      $iskip= ($ireadpart != $icpu);
      $iid++; 
      
      unless($iskip) {
       ($cdsread, $geneids, $genebins)=(0,0,0);  # update only when rdid ne lid
       if($hasCdsBam > 0) { ($cdsread, $geneids, $genebins)= getCdsRead( $rdid, $iid); }
      
        # always count new ids, but for nomap, or $n_readid++ unless($fl & 0xF00); # all 2nd maps
        $n_readid++;
        $inrdlen= length($samx[9]); # read seq, may be '*' or other placeholder
        $LONGR=1 if($inrdlen > $MAXSHORTRD); # FIXME; my $MAXSHORTRD= 500;
        $sum_readlen += $inrdlen; # new global counter, ave(rdlen)= sum_readlen/n_readid
        
        if($DEBUG_RDFILT) {
          warn "#DERD: n_readid=$n_readid, ntmap=$ntmap, iid=$iid, ciid=$cdsiid, cdsids=$cdsread,$geneids, icpu=$icpu\n"
            if($n_readid % 100_000 == 1); # for  Nr.total=421_205_718
        }
      }
    }

    next if($iskip);  # test even for ncpu=1, icpu=1  
    if( GENECOVBINS>=4 and $only_genetab ) { #? after ireadpart == icpu, was before
      unless($geneids) { $no_cdsread++; next;  } 
    }
    
    # limit filters here .. need to do in putmap() w/ all aligns/read *
    my $badmap= 0;
    my($nmi)= (m/NM:i:(\d+)/)?$1:0;    # UPD7f saver: add inreadlen if > 9
    if($fl & 0x4){ $n_nomap++; $badmap=1; } ## next; #? count n_readid here
    elsif(not $STMINFILT and not $LONGR) { #rfQUICKEN:  speed up some, skip low qual quick, MINALN set in putmap2b() 
      if($MINALN>0) {
        my $alen=0; while($cig =~ m/(\d+)M/g) { $alen += $1; } $alen -= $nmi if($nmi>0);
        if($alen < $MINALN) { $n_mapbad++; $badmap=2; } ## next; 
      } 
    } 
       
    if($badmap>0) { 
      putBadmap($badmap, $rdid, $inrdlen, [$cr,$cb,$cig,0,$fl], $cdsread,$geneids, $genebins );

    } else {
      my $saver= [$cr,$cb,$cig,$nmi,$fl]; #? add $fl, check  0x800 supple from bwa == 2nd half map, maybe 
      push @saver, $saver;
      $nokid++ ;
    }
  } close($inhand);

  if($nokid < 0) {
    ($nmap,$crtlist)=(0,0);
  } elsif($LONGR){  
    ($nmap,$crtlist)= putLongmap( $lid, $inrdlen, $nokid, \@saver, $cdsread,$geneids, $genebins); 
  } else {
    ($nmap,$crtlist)= putShortmap( $lid, $inrdlen, $nokid, \@saver, $cdsread,$geneids, $genebins);
  }  
    
  $ntmap+= $nmap;
  # if($SAVE_RDIDS and $nmap>0) { push @saveid, "$lid\t$nmap\t$crtlist"; }  #NotE: crtlist is counts in @crclasses columns
  
  $n_readid++ unless($iskip); # if($ireadpart == $icpu); # always count new ids, but for nomap
  return($ntmap);
}


sub putBadmap {
  my($badmaptype, $rdid, $rdlen, $saveset, $cdsread,$geneids, $genebins)=@_;
   
  my($cr,$cb,$cig,$nmi,$fl)= @$saveset;
  # only process bad map if is CDS read and 1st align to chr: $fl < 0x100
  
  if( $fl < 0x100) { # (UPD21AUG) 
  
    # 'badmap' == 'covz' is cds-read missing chr align
    if( $SAVE_GENECOV and $genebins) {
      my @glocs= split",",$genebins; 
      my $lgid=0;    
      for my $gloc (@glocs){
        my($gid,$gi)=split":",$gloc;
        $covgenebin{$gid}{$gi}{'covz'} += 1;
        unless($gid eq $lgid) {
        $covgene{$gid}{'badmap'} += 1; # maybe add rdlen? not $bdepth; 
        }
        $lgid=$gid;
      }
    }
  
    if($cdsread>0) { my @clbit= (0,1,2,4); 
      $cdsrdclass_miss[4]++; # == any rdclass
      for my $i (1,2,3) { $cdsrdclass_miss[$i]++ if ($cdsread & $clbit[$i]); }
    } else { 
      $cdsrdclass_miss[0]++;  # this is for this-rd not in cdsrd set, dont need
    }
  }
    
}

sub putLongmap {
  my($rdid, $rdlen, $inmap, $saveset, $cdsread,$geneids, $genebins)=@_;
  our(%crmax,%cov,%covt,%covu,%allcov,%allcovt,%allcovu); # ,@savemap,@saveid

  if($SAVE_GENECOV and $geneids){ my @gids= split",",$geneids; $geneids= \@gids; } #UPD21AUG
  if($SAVE_GENECOV and $genebins){ my @gids= split",",$genebins; $genebins= \@gids; } #UPD21AUG
  
  # FIXME: *** putLongmap()  replace addCigar() w/ addCigarDepth()  ***
  
  my(@thismap,%crids,%crtlist);
  my($nmap,$topaln,$lcr)=(0) x 9; # ndi   
  
  for my $saver (@$saveset) {
    my($cr,$cb,$cig,$nmi,$sfl)= @$saver;

    my($alen,$lenc,$cend,$softclip)= addCigar($cr,$cb,$cig,$sfl,$nmi,$rdlen,$inmap,$cdsread,$geneids, $genebins);  
    # ** this should be addCigarDepth( 0|1, .. ) as per shortread ??
    
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
        $MINALN= MIN_MINALN if($MINALN < MIN_MINALN); #?? or not
        warn "# readlen=$readlen, avelen=$areadlen, minaln=$MINALN from $nsample\n" if($debug);
      }

      $topaln=$alen;  $crdtype{$cr}{'mread'}++; 
      #o if($SAVE_RDIDS and $cr ne $lcr){ my($crt)= $idclassh->{$cr} || 0; $crtlist{$crt}++; }
    } else {
      # $mapbad=2 if( $mapbad or $alen < $topaln * $MIN_DUPIDENT);  # count? n_dupbad   
    }
    
    if($mapbad){ 
      if($mapbad==2){ $n_dupbad++; } else { $n_mapbad++; }  
    } else { 
      $nmap++; 
      $n_mapok++; # num maps/read id, count after mapbad, before crok
      $crids{$cr}++; # outgenetab
    }

    $lcr= $cr;      
  }  

  if(GENECOVBINS and ref($geneids) and $nmap>0) {
    for my $gid (@$geneids){ 
      $covgene{$gid}{'nreads'}++;  # == mread above in crdtype
      $covgene{$gid}{'nmaps'} += $nmap;  
    }
  } 
        
#   if($SAVE_RDIDS) {
#     my $crtlist = join"\t", map{ $crtlist{$_}||0 } @$idclasslist; #<<< FIXME, need counts/chrclass for output cols
#     $crtlist.="\t".join(",",sort keys %crids) if($outgenetab); #SAVE_RDIDS & -genetab means save gene id == cr here
#     return($nmap,$crtlist); 
#   }

  if(UPD21AUG  or $NCDS) { 
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
  my($cr,$cb,$cigar,$fl,$nmi,$rdlen,$nmap,$cdsread,$geneids, $genebins)=@_;
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

  if(GENECOVBINS and $geneids) {
    #?? need align/readlen factor for these?
    my $pal= ($cend<1) ? 0: $alen/$cend; my $prv= $pal * $rv;
    for my $gid (@$geneids){  
      $covgene{$gid}{'allcovm'} += $prv;
      $covgene{$gid}{'allcovt'} += $pal; 
      }
  }  
  
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

  UPD21AUG30: addCigarDepth_NU() revise genebins    

=cut

sub addCigarDepth_OLD { 
  my( $addToBins, $aligndepth, $cr,$cb,$cigar,$rdlen,$nmap,$cdsread,$geneids, $genebins)=@_;
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

    # upd GENECOVBINS, drop covgenebin as not helpful
    # TEST 999_2 this way allcovm,t == nreads,nmaps .. otherway below these are higher .. why? readlen/BN = 1.5 maybe dif
    # if(GENECOVBINS > 999_2 and $geneids) {
    #   #?? need align/readlen factor for these?
    #   my $pal= ($cend<1) ? 0: $alen/$cend; my $prv= $pal * $rv;
    #   for my $gid (@$geneids){ # for now same as 8a
    #     $covgene{$gid}{'allcovm'} += $prv;
    #     $covgene{$gid}{'allcovt'} += $pal; 
    #     }
    # }
    
    my $be= int(($cb+$cend)/$BN);  
    $crdtype{$cr}{$rty} += $nmap;
    $crmax{$cr}=$be if($be>$crmax{$cr});
    my $bdepth= 1.0/$BN; # map2b() equiv $bspan/$BN, allcovt counts reads at this (cr,cb+i) base / BINSIZE

    for(my $i=0, my $ifirst=1; $i<=$cend; $i++, $ifirst=0) {
      next unless($thisdepth[$i]); # 0 or 1
      my $bi= int( ($cb+$i)/$BN );
      
      my $adepth= $aligndepth->[$i]||1; # 1..nmap
      my ($rvd,$uvd) = (0,0); # rdv = total rmap / dup maps
      if($adepth > 1) { 
        $rvd= $bdepth/$adepth; # was (1/adepth)/BN
      } elsif($adepth == 1){ 
        $uvd= $rvd= $bdepth;
        $allcovu{$cr}[$bi] += $uvd;   # uniq read depth
      } 
      $allcov{$cr}[$bi]  += $rvd; 
      $allcovt{$cr}[$bi] += $bdepth;  

        # addToBins call for each map of rd, but adepth includes each, ie N x N
        # .. no, this is ok: each map of rd has diff cr,cb .. but this new calc has some much larger total-dupl vals,
        # .. ** map2b() sets bdepth=bspan/BN, ignoring nmap/adepth .. want that again? probably
        # .. where rvd is near same


=item FIXME gene-loc cov counts too high      
    
  FIXME: gi bin loc may not overlap cb+i .. ***
  problem is v8i change to gene-locs vs chr-locs has nearly doubled  cov[tm] vals
    -- may be due to read partly over bins are adding uncovered bin extra, ie ~50% on ave?
    -- maybe pick 1 gene bin only, bin at ave(bb,be)? in getCdsRead
    
  sam2covtab8f
   grep g10178782t1  drosmel6ref_chr_SRR11460802ab_test8f.genetab
  dromel6:g10178782t1	NT_033779.5	16705500	3	3	0
  dromel6:g10178782t1	NT_033779.5	16705600	23	23	0
  dromel6:g10178782t1	NT_033779.5	16705700	17	17	0
  dromel6:g10178782t1	NT_033779.5	16705800	3	3	0
  dromel6:g10178782t1	sumgene	1	43	43	0
  
  New sam2covtab8i
  grep g10178782t1  test8if/drosmel6ref_chr_SRR11460802ab_test8id.genetab
  dromel6:g10178782t1	readgene	43	43	43	0 << same as old sumgene
  dromel6:g10178782t1	sumgene	47	47	47	0
  dromel6:g10178782t1	100	36	36	36	0   << ~50% higher than old max
  dromel6:g10178782t1	200	29	29	29	0   <<  higher

  # ?? change genebins val to genespan base level?  gid:gb:ge
  # ?? need to match read-index here: 0..cend
  # if($cbi >= $gb and $cbi <= $ge) { add cov }
  # NO: next unless($gi == $bi); gi == gene-loc, bi == chr-loc

getCdsRead:
      FIX==0: for(ib=gbin; ib<=gebin; ib++) push @gbins, gid:ib;
      if(GENEBINFIX == 1) {
        my $ge= $gb+$cend;
        push @gbins, "$gid:$gb:$ge";
      } elsif(GENEBINFIX == 2) {
        my $nbin = int(0.3 + $alen/$BN); #? 0.49 maybe
        my $gbin = int($gb/$BN); 
        for(my $i=0; $i<$nbin; $i++) { my $ib=$gbin + $i*$BN; push @gbins, "$gid:$ib"; }
      } 
      
=cut

      if(GENECOVBINS >=5 and $genebins) {
        my($lgid)=(0);

        for my $gloc (@$genebins) { # for now same as 8a
          my($gid,$gi)=split":",$gloc;
if(GENEBINFIX == 1) {
          my($gidx,$gb,$ge)=split":",$gloc; $ge=$gb if($ge<$gb);
          if($gb + $i <= $ge) {
          my $di= int( ($gb+$i)/$BN ); # like bi?
          $covgenebin{$gid}{$di}{'covm'} += $rvd; # or? $bdepth; #? maybe want rdv, udv as per allcov ?
          $covgenebin{$gid}{$di}{'covt'} += $bdepth;
          $covgenebin{$gid}{$di}{'covu'} += $uvd;
          }         
} else { # orig == 0 and GENEBINFIX == 2         
          $covgenebin{$gid}{$gi}{'covm'} += $rvd; # or? $bdepth; #? maybe want rdv, udv as per allcov ?
          $covgenebin{$gid}{$gi}{'covt'} += $bdepth;
          $covgenebin{$gid}{$gi}{'covu'} += $uvd;
}          
          if($gid ne $lgid) { 
          $covgene{$gid}{'allcovm'} += $rvd;
          $covgene{$gid}{'allcovu'} += $uvd;
          $covgene{$gid}{'allcovt'} += $bdepth; # was 'alldepth'
          }
          $lgid=$gid;
        }
      }  
      
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


sub addCigarDepth { # == addCigarDepth_NU 21AUG30 .. seems okay now
  my( $addToBins, $aligndepth, $cr,$cb,$cigar,$rdlen,$nmap,$cdsread,$geneids, $genebins, $imap)=@_;
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

  # addCigarDepth_NU2 : okay now
  if($addToBins) {
    my $rv=1; my $uv=0; my $rty='zero';
    if($nmap==1) { $uv=1; $rty='uniq' } 
    elsif($nmap>1) { $rv= 1.0/$nmap; $rty='mult'; }   

    my $be= int(($cb+$cend)/$BN);  
    $crdtype{$cr}{$rty} += $nmap;
    $crmax{$cr}=$be if($be>$crmax{$cr});
    
    my($cds1,$cds2,$cds3)=(0) x 3;
    if($cdsread > 0) {
     $cds1= ($cdsread & 1);
     $cds2= ($cdsread & 2);
     $cds3= ($cdsread & 4);
    }
        
    my($bi0,@biv)=(-1);
    my $bdepth= 1.0/$BN; 
    for(my $i=0; $i<=$cend; $i++) {
      next unless($thisdepth[$i]); # 0 or 1      
      my $adepth= $aligndepth->[$i]||1; # 1..nmap
      my $bi= int( ($cb+$i)/$BN );
      push @biv, $bi if($bi > $bi0); $bi0= $bi; # try4
      
      my ($rvd,$uvd) = (0,0); # rdv = total rmap / dup maps
      if($adepth > 1) { 
        $rvd= $bdepth/$adepth; # was (1/adepth)/BN
      } elsif($adepth == 1){ 
        $uvd= $rvd= $bdepth;
        $allcovu{$cr}[$bi] += $uvd;   # uniq read depth
      } 
      $allcovt{$cr}[$bi] += $bdepth;  
      $allcov{$cr}[$bi]  += $rvd; 
      
      if($cdsread > 0) {
        # these are all total depth counts, where CDS/TE/UNK exist
        $covt{$cr}[$bi] += $bdepth if($cds1); # 1st col
        $cov{$cr}[$bi]  += $bdepth if($cds2); # 2nd col
        $covu{$cr}[$bi] += $bdepth if($cds3); # 3rd col .. more cols? later
      
#        # TRY_GBINLIM < 4 is obsolete ?
#         if(TRY_GBINLIM < 4 and GENECOVBINS >=5 and $genebins) {
#           my($lgid,$nbin)=(0,0);
# 
#           my ($acovt,$acovm,$acovu)=(0,0,0); # DEBUG_GBINLIM
#           if($DEBUG_GBINLIM){ # try2, set eq to acovt,m,u counts
#             ($acovt,$acovm,$acovu)= ($allcovt{$cr}[$bi], $allcov{$cr}[$bi], $allcovu{$cr}[$bi]);
#           }
# 
#           for my $gloc (@$genebins) { # for now same as 8a
#             my($gid,$gi)=split":",$gloc;
#             
#               # orig == 0 and GENEBINFIX == 2      
#               # COUNT BUG still, this is counting more reads/bin than above allcov{cr}[bi]
#               # .. presumably 2+ genebins/gi bins overlap one bi bin
#               # .. possible fix, limit gloc gi to one count, BUT only for adepth == 1 ? limit to adepth counts?
#               # .. may need to interpolate bi to closest gi, or get 1st gi bin overcount, 2nd undercount
#               # ** last if nbin >= adepth pushes donw overcount, but not enough
#               # ** try2 covgenebin.gid.gi.covt,m,u <= allcov.cr.bi.covt,m,u ; ie keep identical to covtab
#               # .. try2 likely bad for dupx genes as 2nd maps will have diff/lower? allcovt,m,u counts at diff locus
#               # .. try3: covgenebin= max(covgenebin.gid.gi, acov) ? boosts counts, overcount probably
#               
#             #try1: 
#             last if(TRY_GBINLIM == 1 and $DEBUG_GBINLIM and $nbin >= $adepth); #? only for same gid?
#             $nbin++;
#             
#             if(TRY_GBINLIM > 1 and $DEBUG_GBINLIM) { #try2,3
#               
#               if(0) { # try3
#               # try3: looks partly right, partly wrong, overcounts, but some are better than t2
#               # note: also changed to GENEBINFIX => 0, extending genebin spans some
#               sub cmax{ my($h,$d,$i,$t,$v)=@_;  $h->{$d}{$i}{$t}=$v if($v > $h->{$d}{$i}{$t}); }
#               cmax( \%covgenebin, $gid,$gi,'covt',$acovt); # try3:
#               cmax( \%covgenebin, $gid,$gi,'covm',$acovm);
#               cmax( \%covgenebin, $gid,$gi,'covu',$acovu);
#               } else { # try2 ok, undercount?
#               $covgenebin{$gid}{$gi}{'covt'} = $acovt;
#               $covgenebin{$gid}{$gi}{'covm'} = $acovm; 
#               $covgenebin{$gid}{$gi}{'covu'} = $acovu;
#               }
#               
#               if($gid ne $lgid) { 
#               $covgene{$gid}{'allcovt'} = $acovt;  
#               $covgene{$gid}{'allcovm'} = $acovm;
#               $covgene{$gid}{'allcovu'} = $acovu;
#               }
#             
#             } else { # TRY_GBINLIM == 0, DEBUG_GBINLIM == 0
#               $covgenebin{$gid}{$gi}{'covt'} += $bdepth;
#               $covgenebin{$gid}{$gi}{'covm'} += $rvd; 
#               $covgenebin{$gid}{$gi}{'covu'} += $uvd;
#               
#               if($gid ne $lgid) { 
#               $covgene{$gid}{'allcovt'} += $bdepth;  
#               $covgene{$gid}{'allcovm'} += $rvd;
#               $covgene{$gid}{'allcovu'} += $uvd;
#               }
#             }
#             
#             $lgid=$gid;
#           }
#         }  # TRY_GBINLIM < 4 obsolete?
        
      } # has cdsread
      
    } # read align span 0..cend

    if(TRY_GBINLIM >= 4 and $DEBUG_GBINLIM and $cdsread > 0 and $genebins ) # GENECOVBINS >=5 and 
    { 
      my($lgid,$nbin)=(0,0);
      my $nbiv=@biv; my $bi= $biv[ int($nbiv/2) ];
      my($acovt,$acovm,$acovu)= ($allcovt{$cr}[$bi], $allcov{$cr}[$bi], $allcovu{$cr}[$bi]);

      my $crbi = 0;
      if( TRY_GENELOC == 2 ) { 
        (my $ccr=$cr) =~ s/:/_/g; # cr id cant have : here
        my $mi= 10 * int(0.5 + $bi/10); #? is change from bi a problem? count over larger span ?
        $crbi= "$ccr:$mi";  #  bigger bin? int($bi/10) == 1000
      }
      
      ##? how to match chr bi to gene gi? simple way, pick 1 of each 
      # ** ^^ pick 1 maybe not good enough for dupx genes, 
      ## addCigar here has 1-read maps, likely 1 or 2 bins/read tho dupx has diff cr/loci maps per read
      # for my $bi (@biv) {
      #  my($acovt,$acovm,$acovu)= ($allcovt{$cr}[$bi], $allcov{$cr}[$bi], $allcovu{$cr}[$bi]);
      # }
      
      for my $gloc (@$genebins) { # for now same as 8a
        my($gid,$gi)=split":",$gloc;
        if($gid ne $lgid) { # pick 1
        $covgenebin{$gid}{$gi}{'covt'} = $acovt;
        $covgenebin{$gid}{$gi}{'covm'} = $acovm; 
        $covgenebin{$gid}{$gi}{'covu'} = $acovu;
        
        if( TRY_GENELOC == 2 ) { #TEST 21SEP04
        if($imap < 10) { # ? use $imap to limit excess crbi ?
        $covgenecrloc{$gid}{$gi}{$crbi} = $acovm; # +=1 for nreads?   $acovm == covdepth ? NOT += acovm
        }
        # ? check for overloaded {crbi}, limit to 9? vals per gid.gi
        
        }
        
        $covgene{$gid}{'allcovt'} = $acovt;  
        $covgene{$gid}{'allcovm'} = $acovm;
        $covgene{$gid}{'allcovu'} = $acovu;
        }
        $lgid=$gid;
      }
      
    }
    
  } # addToBins

  # FIXME: addCigarDepth returns cend == readlen, 0-base not cb-base, change from prior
  return($alen,$lenc,$cend,$softclip, \@thisdepth);#?
}


=item addCigarDepth_NU1_addToBins

  # drop this way .. slower ??
  
  if($addToBins) {
    # @$aligndepth should be sum of thisdepth over mappings, use for rvcov depth

    my $rv=1; my $uv=0; my $rty='zero';
    if($nmap==1) { $uv=1; $rty='uniq' } 
    elsif($nmap>1) { $rv= 1.0/$nmap; $rty='mult'; }   

    my $be= int(($cb+$cend)/$BN);  
    $crdtype{$cr}{$rty} += $nmap;
    $crmax{$cr}=$be if($be>$crmax{$cr});

    my $bdepth= 1.0/$BN; 
    my $cds1= ($cdsread & 1);
    my $cds2= ($cdsread & 2);
    my $cds3= ($cdsread & 4);
    
    my (%biv,@tdv,@rdv,@udv);
    for(my $i=0; $i<=$cend; $i++) {
      next unless($thisdepth[$i]); # 0 or 1
      my $bi= int( ($cb+$i)/$BN );
      $biv{$bi}++; # or push @biv, $bi if($bi > $bi0); $bi0= $bi;
      
      my $adepth= $aligndepth->[$i] || 1; # 1..nmap
      my ($rvd,$uvd) = (0,0); # rdv = total rmap / dup maps
      if($adepth > 1) { 
        $rvd= $bdepth/$adepth; # was (1/adepth)/BN
      } elsif($adepth == 1){ 
        $uvd= $rvd= $bdepth;
        $udv[$bi] += $uvd;   # uniq read depth
      } 
      $rdv[$bi] += $rvd; 
      $tdv[$bi] += $bdepth;  
    }
    
    for my $bi (sort{ $a<=>$b } keys %biv) {
      my $tdvi= $tdv[$bi]||0;
      my $rdvi= $rdv[$bi]||0;
      my $udvi= $udv[$bi]||0;
     
      $allcovt{$cr}[$bi] += $tdvi;  
      $allcov{$cr}[$bi]  += $rdvi; 
      $allcovu{$cr}[$bi] += $udvi;   # uniq read depth

      # these are all total depth counts, where CDS/TE/UNK exist
      $covt{$cr}[$bi] += $tdvi if($cds1); # 1st col
      $cov{$cr}[$bi]  += $tdvi if($cds2); # 2nd col
      $covu{$cr}[$bi] += $tdvi if($cds3); # 3rd col .. more cols? later
      
      if(GENECOVBINS >=5 and $genebins) {
        my($lgid)=(0);
        for my $gloc (@$genebins) { # for now same as 8a
          my($gid,$gi)=split":",$gloc;
          
          # for GENEBINFIX == 0,2 : not GENEBINFIX == 1
          $covgenebin{$gid}{$gi}{'covt'} += $tdvi;
          $covgenebin{$gid}{$gi}{'covm'} += $rdvi;
          $covgenebin{$gid}{$gi}{'covu'} += $udvi;

          if($gid ne $lgid) { 
          $covgene{$gid}{'allcovt'} += $tdvi;
          $covgene{$gid}{'allcovm'} += $rdvi;
          $covgene{$gid}{'allcovu'} += $udvi;
          }
          $lgid=$gid;
        }
      }
    
    }
  }
  
=cut

    
=item UPD7f putShortmap 
  
  .. uses addCigarDepth() twopass on read saveset to calc 
      1st. valid aligns among multimaps, 2nd. cov depth/valids
  .. replaces putmap2b
  
=cut

sub putShortmap { 
  my($rdid, $rdlen, $nokid, $saveset, $cdsread,$geneids, $genebins)=@_;
  our(%crmax,%cov,%covt,%covu,%allcov,%allcovt,%allcovu); # ,@savemap,@saveid

  if($SAVE_GENECOV and $geneids){ my @gids= split",",$geneids; $geneids= \@gids; } #UPD21AUG
  if($SAVE_GENECOV and $genebins){ my @gids= split",",$genebins; $genebins= \@gids; } #UPD21AUG

  my(@thismap,@aligndepth,%crids); @aligndepth=();
  my($nmap,$topaln,$mindupaln,$lcr)=(0) x 9; # ndi   
  my $nmapin= @$saveset;
  use constant { ADD1ST => 0, ADD2BINS => 2 };
  use constant GENEBIN => 2000; # DEBUG_GENELOC: should be option, boost default to 5k ?
  my $GENEBIN2BN=  GENEBIN / $BN;
  
  for my $saver (@$saveset) {
    my($cr,$cb,$cigar,$nmi,$sfl)= @$saver;
    
    # FIXME: addCigarDepth returns cend == readlen, 0-base not cb-base, change from prior
    my($alen,$lenc,$cend,$softclip, $thisdepth)= 
      addCigarDepth(ADD1ST, [], $cr,$cb,$cigar,$rdlen,$nmapin,$cdsread,$geneids, $genebins); # nmap changes in loop
    
    if($nmi>0){ $alen -= $nmi; $n_mismatch += $nmi; } 
    if($ALLOW_SOFTCLIP and $softclip>0) { $lenc -= $softclip; $n_softclip++; }
    
    my $mapbad=0;    
    if($MINALN>0) { 
      $mapbad=1 if($alen < $MINALN); 
    } else { # sample lenc
      my $mina= $lenc * $MIN_IDENT; $mina= MIN_MINALN if($mina < MIN_MINALN);
      $mapbad=1 if($alen < $mina);
      if($topaln == 0) { # sample only 1st map
      $nsample++; push @sreadlen, $lenc; 
      if($nsample >= NSAMPLE) {
        @sreadlen = sort{$b <=> $a} @sreadlen; $readlen= $sreadlen[ int($nsample/2) ];
        my $areadlen= ($n_readid<2) ? $readlen : int(0.5 + $sum_readlen/$n_readid); # UPD7f
        # if($areadlen != $readlen){ } #some warning, use which?
        $MINALN= int( $readlen * $MIN_IDENT); # was $MINALN= $sminaln / $nsample; 
        $MINALN= MIN_MINALN if($MINALN < MIN_MINALN); # set higher? 10? 30?
        $MINALNCDS= $MINALN; # for now, for getCdsRead(), maybe set below chr minaln due to intron/ends misses, maxmincds?
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
      push @thismap, $saver; # NEW for ADD2BINS
      $crids{$cr}++; #outgenetab
    }

    $lcr= $cr;      
  } # @saveid
  
  $lcr=""; my(%gbins,%crtlist);  my $imap=0;
  for my $saver (@thismap) {
    my($cr,$cb,$cigar,$nmi,$sfl)= @$saver;
    # add imap param? imap/nmap ..
    my($alen,$lenc,$cend,$softclip, $thisdepth)= 
      addCigarDepth(ADD2BINS, \@aligndepth, $cr,$cb,$cigar,$rdlen,$nmap,$cdsread,$geneids, $genebins, ++$imap); 

if($DEBUG_GENELOC) {
    my $cgbin= int( 0.5 +  int(0.5 + $cb/GENEBIN) * $GENEBIN2BN); # GENEBIN == 2000 converted to BN loc
    $gbins{$cr}{$cgbin}++; # this way sums to readgene aCovT,  split into cgbins
    # adding += 1/nmap should sum to readgene aCovM; readgene aCovU when nmap == 1
    #? and/or add $alen;  
}    
    #o if($SAVE_RDIDS and $cr ne $lcr){ my($crt)= $idclassh->{$cr} || 0; $crtlist{$crt}++; }
    $lcr=$cr;
  }

  if(GENECOVBINS and ref($geneids) and $nmap>0) {
    #?? add/replace w/ @$genebins => covgene{gid}{gi}{nreads}++
    my $ngid= @$geneids;
    for my $gid (@$geneids){ 
      $covgene{$gid}{'nreads'}++;  # == mread above in crdtype
      $covgene{$gid}{'nmaps'} += $nmap;  
      $covgene{$gid}{'nuniq'} ++ if($nmap == 1); # and $ngid == 1  #??
      
if($DEBUG_GENELOC) {
      for my $cr (sort keys %gbins) {
        for my $cb (sort keys %{$gbins{$cr}}) {
          my $vc= $gbins{$cr}{$cb}; # == readgene.aCovT
          #? add equiv for readgene.aCovM ?
          $covgene{$gid}{'crloc'}{$cr}{$cb} += $vc; # vc? or +1, may become huge w/ big-pig data
        }
      }  
}        
    }
  } 
      
#   if($SAVE_RDIDS) {
#     my $crtlist = join"\t", map{ $crtlist{$_}||0 } @$idclasslist; #<<< FIXME, need counts/chrclass for output cols
#     $crtlist.="\t".join(",",sort keys %crids) if($outgenetab); #SAVE_RDIDS & -genetab means save gene id == cr here
#     return($nmap,$crtlist);
#   } 
  
  if(UPD21AUG  or $NCDS) { 
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

# sub putrdids { 
#   my($rdidtab)= @_; # or $outh
#   warn "# rdidtab $rdidtab\n" if($debug);
#   rename($rdidtab,"$rdidtab.old") if(-s $rdidtab);
#   open(OT,">$rdidtab"); 
#   #see above:  push @saveid, join("\t", '#ReadID','nmap',@$idclasslist); # add col header for @crclass
#   for my $sid (@saveid) { 
#     # FIXME for mixed readsets, remove or not SRRprefix. from rdnum here
#     # CRTPAT adds 3rd column: rdid, nmap, chrtaglist << chrtaglist must be columns of counts/@$idclasslist
#     print OT $sid,"\n";   # push @saveid, "$lid\t$nmap" 
#     }
#   close(OT);
# }

sub putGenetab8i {
  my($chrtab)= @_; # == genetab
  warn "# putgenetab $chrtab\n" if($debug);
  rename($chrtab,"$chrtab.old") if(-s $chrtab);
  open(OUT,">$chrtab"); 

  print OUT "# genetab version=$VERS\n";
  print OUT "# cols for GeneID, GPos: aCovT aCovM aCovU aCovZero = read map counts/bin as per chr.covtab\n";
  print OUT "# cols for sumgene,  1: allCovT allCovM allCovU noCov\n"; # gave=allCovT/nbins
  print OUT "# cols for readgene, 1: nrdmaps, nreads, nuniq, nomap\n";
  
  # putGenetab8i: table col changes; besides sumgene/readgene non-locus rows, want some kind of GeneID x ChrID
  # possible: GeneID:[abc]bin ChrID CPos covs... where GeneID:[abc]bin are pre-selected representative bins/gene,
  #  select from cds.covtab, median/25%/75% coverage bins? avoid extremes
  # for gene 2-9 dupx, fully dupx will have 2-9 diff chr:loc bins for each [abc]
  # for genes skew/partial dups, have variable counts at gene:[abc], 
  # for xhigh-copy (TE) genes, need some limit on chr:loc bins, to keep mem use in bounds
  
  my @ocols=qw(GeneID GPos aCovT aCovM aCovU aCovZ); # aCovZ == noCov as per 8e, aCovU uniq rd cov as per other covtab 
    # adds GeneID 1st? or last? for merge better 1st
  print OUT "#".join("\t",@ocols)."\n";
  
  my $areadlen= ($n_readid<2) ? $readlen : int(0.5 + $sum_readlen/$n_readid); # UPD7f

  #UPD21AUG: Fixme add back GIDPREFIX if chopped by cleangid(), for agreement w/ other tables, 
  my $addGIDPRE= ($GIDPREFIX and $ncutGIDPRE>0);
  
  my @gid= sort keys %covgene; # do we know gene lengths?
  my $ngene= @gid; 
  my($tgmaplen,$tallmap,$tallmapm,$tallmapu,$tnmaps,$tnuniq,$tnreads,$tnomap,$tnodepth,$ttnb,$tmedcovt,$tmedcovm)=(0) x 19; 
  my(@gmed,@gave);
  for my $gid (@gid) {
    my $gidout= ($addGIDPRE)? $GIDPREFIX.$gid : $gid;
    my $alldepth= $covgene{$gid}{'allcovt'}||0; # merge-able count
    my $alldepthm= $covgene{$gid}{'allcovm'}||0; # merge-able count
    my $alldepthu= $covgene{$gid}{'allcovu'}||0; # merge-able count
    my $nreads= $covgene{$gid}{'nreads'}||0; # merge-able count
    my $nmaps = $covgene{$gid}{'nmaps'}||0; # merge-able count
    my $nuniq = $covgene{$gid}{'nuniq'}||0; # merge-able count
    my $nomap = $covgene{$gid}{'badmap'}||0;   # merge-able count  #? badmap == allcovz
 
    my $ncrloc=0; my @vcrloc=();
 if($DEBUG_GENELOC) {
    # NEED Mergable table: 
    #  geneid 'locene' vcol1 vcol2 vcol3 vcol4
    #  geneid locgene1  chr1  ib1  nread1  0
    #  geneid locgene9  chr9  ib9  nread9  0
    if(exists $covgene{$gid}{'crloc'}) {
    my @cr= sort keys %{$covgene{$gid}{'crloc'}};
    for my $cr (@cr) { 
      my @cb= sort{ $a<=>$b } keys %{$covgene{$gid}{'crloc'}{$cr}}; 
      #? sort output by read count cnrd, pick top 9
      for my $cb (@cb) { 
        my $cnrd= $covgene{$gid}{'crloc'}{$cr}{$cb}||0;
        my $ib=$BN*(1+$cb);
        push @vcrloc, join("\t",$cr,$ib,$cnrd,0); # output vcols
        $ncrloc++; last if($ncrloc>=9);
        }
      last if($ncrloc>=9);
    }
    }
}
    
    my $tnb=0; my(@covt,@covm);
    # my @cr= sort{ $crmax{$b}<=>$crmax{$a} } keys %{$covgenebin{$gid}};  
    my @ibins= sort{$a <=> $b} keys %{$covgenebin{$gid}}; # {$cr}
    for my $i (@ibins) {
      # my @allv= ($allcovt{$cr}[$i], $allcov{$cr}[$i], $allcovu{$cr}[$i]); # reuse allcov bin counts
      # my $covm=$covgenebin{$gid}{$i}{'covm'}; # GENECOVBINS == 5
      # my $covt=$covgenebin{$gid}{$i}{'covt'}; # 
      # my $covu=$covgenebin{$gid}{$i}{'covu'}; # 
      # my $covz=$covgenebin{$gid}{$i}{'covz'}; # == badmap/nomap

      my @ocols= map{ $covgenebin{$gid}{$i}{$_}||0 } qw( covt covm covu covz );
      @ocols= map{ int( 0.5 + $_) } @ocols; # ( $covt,$covm,$covu,$covz );
      my $ib=$BN*(1+$i);  # < end of bin, or $ib= ($i==0) ? 1 : $BN * $i;  == start of bin 
      print OUT join("\t",$gidout, $ib, @ocols)."\n";  $tnb++;
      
      # try for gene bin x chr loc rows? or add to ocols?
      # need merge-able 3 value tuple: crid/crbin/count, following geneid/genebin, or crid:crbin as 1 col
      #  geneid/gbin/cr:crbin/count ? or even gid:gbin/cr:cbin/count for 3 cols, merge = sum(count)
      # ? belongs w/ locgene/$DEBUG_GENELOC
      # ** Only useful when locgene, ncrloc > 1 
      
      if(TRY_GENELOC == 2 and exists $covgenecrloc{$gid}{$i}) { 
        # key == cr:crbin, val == covm/acovm = covdepth, and/or nreads count?
        # ? how to make this mergeable, like above crloc?  cr/ib/nrd 
        # * add 'cr:' prefix to col3 for merge 
        # FIXME: cr:cb , cb is MBN bin val, output $ib=$MBN*(1+$cb); 
        #  if( TRY_GENELOC == 2 ) { my $mi= 10 * int(0.5 + $bi/10); $crbi= "$cr:$mi"; } #  bigger bin? int($bi/10) == 1000
        my @crmb= sort keys %{ $covgenecrloc{$gid}{$i} }; #? sort by +cv count?
        for my $crmb (@crmb) {
          my $cv= $covgenecrloc{$gid}{$i}{$crmb}||0; $cv=int( 0.5 + $cv );
          my($cr,$cb)= split/:/,$crmb; $cb= $BN*(1+$cb); 
          print OUT join("\t",$gidout, $ib, "cr:$cr:$cb", $cv)."\n";   # needs merge fix  
        }
      }
      
    }
   
    #? replace gave w/ nodepth = nomap * readlen ? is merge-able; NOTE /BN as alldepth is for BN sizes
    my $nodepth= $areadlen/$BN * $nomap;
    ($alldepth,$alldepthm,$alldepthu,$nodepth)= map{ int(0.5+ $_ ) } ($alldepth,$alldepthm,$alldepthu,$nodepth); 
    $tallmap+=$alldepth; $tallmapm+=$alldepthm; $tnodepth += $nodepth; $tallmapu+=$alldepthu;
    $tnmaps+=$nmaps;  $tnuniq+= $nuniq; $tnreads+=$nreads; $tnomap+=$nomap; 

    ## dont need both these sum,read count .. sum is better stat
    print OUT join("\t",$gidout, 'sumgene', $alldepth, $alldepthm, $alldepthu, $nodepth)."\n";  # note: alldepth > nreads by >=10% 
    print OUT join("\t",$gidout, 'readgene', $nmaps, $nreads, $nuniq, $nomap)."\n";  #upd21jul04: was off, want for match to other tabs
 if($DEBUG_GENELOC) {
    for(my $i=1; $i<=$ncrloc; $i++)  {
      print OUT join("\t",$gidout, "locgene$i", $vcrloc[$i-1])."\n"; #? locgene1..locgene9 ?
    } 
}
  }

  my $Cave = ($tallmapm < 1)? 0: sprintf"%.1f",$tallmap/$tallmapm;
  my $pMiss= ($tnreads  < 1)? 0: sprintf"%.3f",100*$tnomap/$tnreads;
  
  print OUT join("\t",'total','sumgene',$tallmap, $tallmapm, $tallmapu, $tnodepth)."\n"; #  
  print OUT join("\t",'total','readgene', $tnmaps, $tnreads, $tnuniq, $tnomap)."\n";  #   
  close(OUT); 
  warn "#genecov: ngene=$ngene, C.ave=$Cave, missed=$pMiss% in  $chrtab\n" if($debug);   
  return($ngene);
  
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

__END__

sub putGenetab8e { # if(GENECOVBINS == 4)
  my($chrtab)= @_; # == genetab
  warn "# putgenetab $chrtab\n" if($debug);
  rename($chrtab,"$chrtab.old") if(-s $chrtab);
  open(OUT,">$chrtab"); 
  
  print OUT "# cols for ChrID, Pos: CovT CovM zero = read map counts/bin as per chr.covtab\n";
  print OUT "# cols for sumgene,  1: allCovT allCovM noCov\n"; # gave=allCovT/nbins
  print OUT "# cols for readgene, 1: nrdmaps, nreads, nomap\n";
  my @ocols=qw(GeneID ChrID Pos aCovT aCovM noCov); # adds GeneID 1st? or last? for merge better 1st
  print OUT "#".join("\t",@ocols)."\n";

  my $areadlen= ($n_readid<2) ? $readlen : int(0.5 + $sum_readlen/$n_readid); # UPD7f
  
  my @gid= sort keys %covgene; # do we know gene lengths?
  my $ngene= @gid; 
  my($tgmaplen,$tallmap,$tallmapm,$tnmaps,$tnreads,$tnomap,$tnodepth,$ttnb,$tmedcovt,$tmedcovm)=(0) x 9; 
  my(@gmed,@gave);
  for my $gid (@gid) {

    my $alldepth= $covgene{$gid}{'allcovt'}||0; # merge-able count
    my $alldepthm= $covgene{$gid}{'allcovm'}||0; # merge-able count
    my $nreads= $covgene{$gid}{'nreads'}||0; # merge-able count
    my $nmaps = $covgene{$gid}{'nmaps'}||0; # merge-able count
    my $nomap = $covgene{$gid}{'badmap'}||0;   # merge-able count
 
    my $tnb=0; my(@covt,@covm);
    my @cr= sort{ $crmax{$b}<=>$crmax{$a} } keys %{$covgenebin{$gid}};  
    for my $cr (@cr) {
      my @ibins= sort{$a <=> $b} keys %{$covgenebin{$gid}{$cr}};
      for my $i (@ibins) {
        # my @allv= ($allcovt{$cr}[$i], $allcov{$cr}[$i], $allcovu{$cr}[$i]); # reuse allcov bin counts
        my $covm=$covgenebin{$gid}{$cr}{$i}{'covm'}; # GENECOVBINS == 4
        my $covt=$covgenebin{$gid}{$cr}{$i}{'covt'};
        my @ocols= map{ (defined $_) ? int( 0.499 + $_) : 0 } ( $covt,$covm,0 );
        my $ib=$BN*(1+$i);  # < end of bin, or $ib= ($i==0) ? 1 : $BN * $i;  == start of bin 
        print OUT join("\t",$gid, $cr, $ib, @ocols)."\n";  $tnb++;
        # push @covt,$ocols[0]; push @covm, $ocols[1]; 
      }
    }
   
    #? replace gave w/ nodepth = nomap * readlen ? is merge-able; NOTE /BN as alldepth is for BN sizes
    my $nodepth= $areadlen/$BN * $nomap;
    ($alldepth,$alldepthm,$nodepth)= map{ int(0.5+ $_ ) } ($alldepth,$alldepthm,$nodepth); 
    $tallmap+=$alldepth; $tallmapm+=$alldepthm; $tnodepth += $nodepth;
    $tnmaps+=$nmaps; $tnreads+=$nreads; $tnomap+=$nomap; 
    # $ttnb+=$tnb; $tgmaplen+= $gmaplen;
    # $tmedcovt+= $medcovt; $tmedcovm+=$medcovm;

    ## dont need both these sum,read count .. sum is better stat
    print OUT join("\t",$gid, 'sumgene', 1, $alldepth, $alldepthm, $nodepth)."\n";  # note: alldepth > nreads by >=10% 
    print OUT join("\t",$gid, 'readgene', 1, $nmaps, $nreads, $nomap)."\n";  #upd21jul04: was off, want for match to other tabs
  }

  my $Cave = ($tallmapm < 1)? 0: sprintf"%.1f",$tallmap/$tallmapm;
  my $pMiss= ($tnreads  < 1)? 0: sprintf"%.3f",100*$tnomap/$tnreads;
  
  print OUT join("\t","total",'sumgene',$tallmap,$tallmapm,$tnodepth)."\n"; # was gavmed?
  print OUT join("\t",'total', 'readgene', $tnmaps, $tnreads, $tnomap)."\n";    
  close(OUT); 
  warn "#genecov: ngene=$ngene, C.ave=$Cave, missed=$pMiss% in  $chrtab\n" if($debug);   
  return($ngene);
}

sub putGenetab8c { # if(GENECOVBINS == 3)
  my($chrtab)= @_; # == genetab
  warn "# putgenetab $chrtab\n" if($debug);
  rename($chrtab,"$chrtab.old") if(-s $chrtab);
  open(OUT,">$chrtab"); 

  # FIXME: drop cds CovT CovM CovU ; only aCovT aCovM aCovU here???
  
  # print OUT "# cols for ChrID: aCovT aCovM aCovU read map counts/bin as per chr.covtab\n";
  print OUT "# cols for sumgene : allCovT, allCovM, noCov\n"; # gave=allCovT/nbins
  print OUT "# cols for readgene: nrdmaps, nreads, nomap\n";
  my @ocols=qw(GeneID class aCovT aCovM noCov); # adds GeneID 1st? or last? for merge better 1st
  print OUT "#".join("\t",@ocols)."\n";

  my $areadlen= ($n_readid<2) ? $readlen : int(0.5 + $sum_readlen/$n_readid); # UPD7f
  
  my @gid= sort keys %covgene; # do we know gene lengths?
  my $ngene= @gid; 
  my($tgmaplen,$tallmap,$tallmapm,$tnmaps,$tnreads,$tnomap,$tnodepth,$ttnb,$tmedcovt,$tmedcovm)=(0) x 9; 
  my(@gmed,@gave);
  for my $gid (@gid) {

    my $alldepth= $covgene{$gid}{'allcovt'}||0; # merge-able count
    my $alldepthm= $covgene{$gid}{'allcovm'}||0; # merge-able count
    my $nreads= $covgene{$gid}{'nreads'}||0; # merge-able count
    my $nmaps = $covgene{$gid}{'nmaps'}||0; # merge-able count
    my $nomap = $covgene{$gid}{'badmap'}||0;   # merge-able count
    # my $gmaplen= $BN * $tnb; #? skip this val, use tnb
    # my $gave= ($tnb<1)?0: sprintf"%.2f",$alldepth/$tnb; # tnb, not gmaplen= BN*tnb
    # push @gave, $gave;
    
    #? replace gave w/ nodepth = nomap * readlen ? is merge-able; NOTE /BN as alldepth is for BN sizes
    my $nodepth= $areadlen/$BN * $nomap;
    ($alldepth,$alldepthm,$nodepth)= map{ int(0.5+ $_ ) } ($alldepth,$alldepthm,$nodepth); 
    # $ttnb+=$tnb; $tgmaplen+= $gmaplen;
    $tallmap+=$alldepth; $tallmapm+=$alldepthm; $tnodepth += $nodepth;
    $tnmaps+=$nmaps; $tnreads+=$nreads; $tnomap+=$nomap; 
    # $tmedcovt+= $medcovt; $tmedcovm+=$medcovm;

    ## dont need both these sum,read count .. sum is better stat
    print OUT join("\t",$gid, 'sumgene', $alldepth, $alldepthm, $nodepth)."\n";   # gave is not merge-able
    print OUT join("\t",$gid, 'readgene', $nmaps, $nreads, $nomap)."\n";    
  }

  # @gave=sort{$b <=> $a} @gave; my $gavmed= $gave[ int($ngene/2)];
  # @gmed=sort{$b <=> $a} @gmed; my $gmedmed= $gmed[ int($ngene/2)];

  my $Cave = ($tallmapm < 1)? 0: sprintf"%.1f",$tallmap/$tallmapm;
  my $pMiss= ($tnreads < 1)? 0: sprintf"%.3f",100*$tnomap/$tnreads;
  
  print OUT join("\t","total",'sumgene',$tallmap,$tallmapm,$tnodepth)."\n"; # was gavmed?
  print OUT join("\t",'total', 'readgene', $tnmaps, $tnreads, $tnomap)."\n";    
  close(OUT); 
  warn "#genecov: ngene=$ngene, C.ave=$Cave, missed=$pMiss% in  $chrtab\n" if($debug);   
  return($ngene);
}

sub putGenetab8b { # if(GENECOVBINS == 2)
  my($chrtab)= @_; # == genetab
  warn "# putgenetab $chrtab\n" if($debug);
  rename($chrtab,"$chrtab.old") if(-s $chrtab);
  open(OUT,">$chrtab"); 

=item gene covtab bins
# PROBLEM: this may be wanted only for merged chr.covtab bins,
#   merge of genetab gets wacky intron/end dips in cov due to few of part map having reads in those bins

#GeneID 	ChrID   	Pos	aCovT	aCovM	aCovU
AT1G01010t1	NC_003070.9	3700	1	1	1
AT1G01010t1	NC_003070.9	3800	12	12	12
AT1G01010t1	NC_003070.9	3900	26	26	26
AT1G01010t1	NC_003070.9	4000	23	23	23
AT1G01010t1	NC_003070.9	4100	9	9	9
AT1G01010t1	NC_003070.9	4200	12	12	12
AT1G01010t1	NC_003070.9	4300	12	12	12
AT1G01010t1	NC_003070.9	4400	2	2	2

=cut

  # FIXME: drop cds CovT CovM CovU ; only aCovT aCovM aCovU here???
  
  print OUT "# cols for ChrID: aCovT aCovM aCovU read map counts/bin as per chr.covtab\n";
  print OUT "# cols for sumgene : nbins, allCovT, allCovM, noCovT\n"; # gave=allCovT/nbins
  print OUT "# cols for readgene: nbins, nrdmaps, nreads, nomap\n";
  my @ocols=qw(GeneID ChrID Pos aCovT aCovM aCovU); # adds GeneID 1st? or last? for merge better 1st
  print OUT "#".join("\t",@ocols)."\n";

  my $areadlen= ($n_readid<2) ? $readlen : int(0.5 + $sum_readlen/$n_readid); # UPD7f
  
  my @gid= sort keys %covgene; # do we know gene lengths?
  my $ngene= @gid; 
  my($tgmaplen,$tallmap,$tallmapm,$tnmaps,$tnreads,$tnomap,$tnodepth,$ttnb,$tmedcovt,$tmedcovm)=(0) x 9; 
  my(@gmed,@gave);
  for my $gid (@gid) {
    my $tnb=0; my(@covt,@covm);
    my @cr=sort{ $crmax{$b}<=>$crmax{$a} } keys %{$covgenebin{$gid}};  
    ## print alldepth, nomap as GeneID alldepth 1 ... aCovT x x, same col format as cr,ib
    for my $cr (@cr) {
      my @ibins= sort{$a <=> $b} keys %{$covgenebin{$gid}{$cr}};
      for my $i (@ibins) {
        # my @cdsv= ($NCDS) ? ( $covt{$cr}[$i], $cov{$cr}[$i], $covu{$cr}[$i]) : (0,0,0);
        my @allv= ($allcovt{$cr}[$i], $allcov{$cr}[$i], $allcovu{$cr}[$i]); # reuse allcov bin counts
        @ocols= map{ (defined $_) ? int( 0.5 + $_) : 0 } ( @allv );
        my $ib=$BN*(1+$i);  # < end of bin, or $ib= ($i==0) ? 1 : $BN * $i;  == start of bin 
        print OUT join("\t",$gid, $cr, $ib, @ocols)."\n";  $tnb++;
        push @covt,$ocols[0]; push @covm, $ocols[1]; 
      }
    }

    ## add 'nreads' 'nmaps'
    my $alldepth= $covgene{$gid}{'allcovt'}||0; # merge-able count
    my $alldepthm= $covgene{$gid}{'allcovm'}||0; # merge-able count
    my $nreads= $covgene{$gid}{'nreads'}||0; # merge-able count
    my $nmaps = $covgene{$gid}{'nmaps'}||0; # merge-able count
    my $nomap = $covgene{$gid}{'badmap'}||0;   # merge-able count
    my $gmaplen= $BN * $tnb; #? skip this val, use tnb
    my $gave= ($tnb<1)?0: sprintf"%.2f",$alldepth/$tnb; # tnb, not gmaplen= BN*tnb
    push @gave, $gave;

    #?? print median val from all @cr.@ibins ? for allcovt, allcov ?
    my $hnb=int($tnb/2);
    @covt= sort{$b<=>$a}@covt; my $medcovt= $covt[$hnb];
    @covm= sort{$b<=>$a}@covm; my $medcovm= $covm[$hnb];
    push @gmed, $medcovm;
    
    ($alldepth,$alldepthm,$medcovt,$medcovm)= map{ sprintf"%.0f",$_ } ($alldepth,$alldepthm,$medcovt,$medcovm); 
    $ttnb+=$tnb; $tgmaplen+= $gmaplen;
    $tallmap+=$alldepth; $tallmapm+=$alldepthm; 
    $tnmaps+=$nmaps; $tnreads+=$nreads; $tnomap+=$nomap; 
    $tmedcovt+= $medcovt; $tmedcovm+=$medcovm;
    
    #?? fixme? gmaplen, tnb are variable over icpu parts so dont add on merge, unless fixed there
    #? replace gave w/ nodepth = nomap * readlen ? is merge-able
    my $nodepth= $areadlen * $nomap; $tnodepth += $nodepth;
    print OUT join("\t",$gid, 'sumgene', $tnb, $alldepth, $alldepthm, $nodepth)."\n";   # gave is not merge-able
    print OUT join("\t",$gid, 'readgene', $tnb, $nmaps, $nreads, $nomap)."\n";    
    #?? print OUT join("\t",$gid, 'medgene', $tnb, $medcovt, $medcovm, 0)."\n";   # gave is not merge-able
  }

  @gave=sort{$b <=> $a} @gave; my $gavmed= $gave[ int($ngene/2)];
  @gmed=sort{$b <=> $a} @gmed; my $gmedmed= $gmed[ int($ngene/2)];

  print OUT join("\t","total",'sumgene',$tgmaplen,$tallmap,$tallmapm,$tnodepth)."\n"; # was gavmed?
  print OUT join("\t",'total', 'readgene', $ttnb, $tnmaps, $tnreads, $tnomap)."\n";    
  # print OUT join("\t","total",'medgene',$ngene,$gmedmed,$tmedcovt,$tmedcovm)."\n";  
  close(OUT); 
  warn "#genecov: ngene=$ngene, C.ave=$gavmed in  $chrtab\n" if($debug);   
  return($ngene);
}


sub putGenetab8a { # if(GENECOVBINS == 1)
  my($chrtab)= @_; # == genetab
  warn "# putgenetab $chrtab\n" if($debug);
  rename($chrtab,"$chrtab.old") if(-s $chrtab);
  open(OT,">$chrtab"); 

# FIXME swap chr tables for gene tables
#          $covgenebin{$gid}{$cr}[$bi] += $bdepth; 
#          $covgene{$gid}{'alldepth'} += $bdepth; # also badmap, counts reads? or rdlen?

# FIXME: merge for this table: NOT sum all columns.. need max(Galign), sum(Nmap,Nomap), 
# .. merge: recalc(Cmed,Cave); Cave= sum(Nmap) / max(Galign) ??  Cmed = ?? cant do merge w/o more data
# ** change to output all of counts per gid,chr,cbin for merge, ie like chr.covtab but w/ added GeneID level

  print OT "#".join("\t",qw(GeneID Galign Cmed Cave Nmap Nomap Nitem))."\n";
  my @gid= sort keys %covgene; # do we know gene lengths?
  my ($smaplen,$tallmap,$tnomap)=(0) x 9; my(@gmed,@gave);
  for my $gid (@gid) {
    my $allmap=$covgene{$gid}{'alldepth'}||0;
    my $nomap= $covgene{$gid}{'badmap'}||0;
    my @cr=sort keys %{$covgenebin{$gid}};  
    
    # START_sam2covtab Mon Jun 14 20:30:38 EDT 2021
    # putgenetab arath18tair_chr_SRR10178325_test8a.ts8a.genetab
    # Can't use string ("270.159999999899") as a HASH ref while "strict refs" in use at gnodes_sam2covtab8a.pl line 1016.
    # DONE_sam2covtab Tue Jun 15 02:12:52 EDT 2021
    #x @cr = grep{ $crmax{$_} and not m/alldepth|badmap/ } @cr;
    
    my ($tnb,@bcov)= (0);
    for my $cr (@cr){ 
      my @ibins= sort{$a <=> $b} keys %{$covgenebin{$gid}{$cr}};
      #? can we parse tandem gene spans from ibins gaps ? need to know more of cds-seq size, 
      my $nb=@ibins; $tnb+= $nb;
      for my $i (@ibins){ push @bcov, $covgenebin{$gid}{$cr}{$i}; }
      }
      
    @bcov=sort{$b <=> $a} @bcov;
    my $gmaplen= $BN * $tnb;
    my $gmed= $bcov[int($tnb/2)];
    my $gave= ($gmaplen<1)?0: sprintf"%.2f",$allmap/$tnb; # * $allmap/$gmaplen is /$BN too small
    ($gmed,$allmap)= map{ sprintf "%.1f",$_ } ($gmed,$allmap);
    
    # PROBLEM: Most of these vals are not merge-able over read-subsets for ncpu, only allmap?
    # .. would need to save chr x cbin list of counts for all genes to merge.
    # this: $covgenebin{$gid}{$cr}{$i} could be added to asm_rd.covtab as gene ids per bin, then use in merge for
    #   similar genetab output?
    
    # is nomap == rdcount or rdlen sum ?
    print OT join("\t",$gid,$gmaplen,$gmed,$gave,$allmap,$nomap, $tnb)."\n";
    
    push @gmed,$gmed; push @gave,$gave;
    $smaplen+=$gmaplen; $tallmap+=$allmap; $tnomap+=$nomap;
  }
  my $ngene= @gid; my $nh=int($ngene/2);
  @gmed=sort{$b <=> $a} @gmed; @gave=sort{$b <=> $a} @gave; 
  my $gmedmed= $gmed[ $nh ]; my $gavmed= $gave[ $nh];
  print OT join("\t","#total",$smaplen,$gmedmed,$gavmed,sprintf("%.0f",$tallmap),$tnomap,$ngene)."\n";
  close(OT); 
  warn "#$chrtab: ngene=$ngene, C.ave=$gavmed\n" if($debug);   
  return($ngene);
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
