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

my $VERS='8j'; # 8i: remove cds.readids, add iterative cds.bam x chr.bam read map pairings
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

use constant UPD21SEP => 1; # test revised cds/genetab binning, fixed cds-points to meaure chr loci at: 100,200,..,len/100 base points
  #? replace genetab hashes w/ arrays for gid.gbin, use input cds.bam hdr id/length > fixed bins
use constant UPD21DEC => 1; # CDSBAM classed w/ idclass table CDS,TE,UNK
use constant UPD22JAN => 1; # if(TEST_FRACTIONS)  #UPD.TEST 22JAN02 
my $TEST_FRACTIONS= UPD22JAN; #default now, debug opt to turn off? $ENV{TEST_FRACTIONS}||0;
   $TEST_FRACTIONS= 0 if($ENV{NO_FRACTIONS});
use constant UPD22MAR => 1; #  pairmap, same read ids fwd/rev, fix revid
   
use constant { SAMFLAG_nomap => 0x04, SAMFLAG_rev => 0x10, SAMFLAG_read2 => 0x80, SAMFLAG_2ndary => 0x100, SAMFLAG_suppl => 0x800 }; 
my $isPAIRMAP=0; #UPD22MAR, dont need this opt, auto-correct _rev ids
  
my $debug=$ENV{debug}||0; # test
my $DEBUG_RDFILT= $ENV{DEBUG_RDFILT}||0; #UPD21AUG: DEBUG mem overflows in readfilter with big-pig data set
# my $DEBUG_GENELOC= $ENV{DEBUG_GENELOC}||0; #UPD21AUG: DEBUG mem overflows in readfilter with big-pig data set
# my $DEBUG_GBINLIM= $ENV{DEBUG_GBINLIM}||0; #UPD21AUG: DEBUG genetab over-counts

use constant TRY_GBINLIM => 4; 
my $DEBUG_GBINLIM= TRY_GBINLIM;
use constant TRY_GENELOC => 2; 
my $DEBUG_GENELOC= TRY_GENELOC;
my $DEBUG_CRCB= 0; # $ENV{DEBUG_CRCB}||0; # UPD21SEP test for genetab gid.gb.cr.cb locus points ** CRCB=1 result is worse 

my $MIN_IDENT= 0.40; # UPD7e was 0.65; # lo is best?
my $MIN_DUPIDENT = 0.98; # 7f= 0.999; #UPD7e was 0.98; # was .99/1; hi is best? or 1.0; # == ident equal to 1st/top align, lower if desired

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
my $SAVE_GENECOV= 0;  # now == hasCdsBam
#obs: my $GeneIdsInBam= $ENV{'GeneIdsInBam'} || 0; # UPD21AUG20 ** Drop for cdsbam reading
my $USE_CDSBAMCLASS=1; # UPD21DEC, maybe make default, depends on idclass table values (CDS,TE,UNK)

my $ALN_SUB_NMI = (defined $ENV{nonmi}) ? 0 : 1; # UPD21NOV: TEST for hetero adjust, dont subtract NM:i:nn  from align val


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
  'usecdsclass!',\$USE_CDSBAMCLASS, # requires -idclass/crclassf and hasCdsBam
  'nmi!',\$ALN_SUB_NMI, # -nonmi turns off
  'binsize=i',\$BN, 'skipread=i',\$SKIPR, 
  'merge:s',\$MERGE, # UPD7, flag to merge ncpu part tables : 'merge!',
#  'PAIRMAP!', \$isPAIRMAP, #UPD22MAR
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
my($NOK,$haslen,%crlen,%crdtype)=(0,0,0); #8e.drop: %crok, $hasread,,%cread 
  
# global readid data
use constant USE_RIDp => 1; #  DROP; problem if this is big memory pig
my ($NCDS,$hasCdsBam,$useRIDp)=(0,0,USE_RIDp);  
my @cdsrdclass_in=(0) x 5; my @cdsrdclass_sum=(0) x 5; my @cdsrdclass_miss=(0) x 5; 
# my %cdslen=(); # now @agenesize; from cds.bam header, or table of id/len ??

our( %crmax, # chr max loc found in chr.bam; see also %crlen should be ~same
  %cov, %covt, %covu, # UPD21AUG20 how are these used now?
  %allcov, %allcovt, %allcovu, # primary chr cov data
  # @savemap, @saveid, # UPD21AUG20 drop usage
  %covgene, %covgenebin, # UPD21AUG20 change struct
  %covgenecrloc,  # TRY_GENELOC == 2 UPD21SEP04
  #UPD21SEP ^v replace w/ @covgene[igene], .. arrays indexing geneids, to save mem, one hash geneindex{gid} = igene++
  %geneindex, # indices from gene id/size table
     @agenesize, # replace  %cdslen
     @agenebins, @ageneclass, #UPD21DEC class  
  @agenecov, @agenebincov, @agenecrloc, # cover counts from bams
);

#UPD7: change %cov,covt,covu to single hash of cdste-readid-hit-by-class
# FIXME: 2 uses of chrtab: input for chr, crlen, reads,  differs from output : two opt names?
# .. below outchrpart should not eq chrtab input name
  
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
 
  $hasCdsBam= ($cdsbam and -s $cdsbam)?1:0;
  # FIXME: -genetab opt still needed for -merge w/o -cdsbam ; should look for .genetab.pt* also
  if($hasCdsBam or defined $outgenetab) { # was (defined $outgenetab)
    $SAVE_GENECOV=  1;
    unless($outgenetab =~ m/\w\w/) {  $outgenetab= $outtab; $outgenetab=~s/\.\w+$//; $outgenetab.=".genetab"; }    
  }
  
  my $outpart= $outtab;
  my $outchrpart= $outchrtab; if($outchrpart eq $inchrtab) { $outchrpart.="out"; } # name bug: .chrtabout with -merge
  my $outgenepart= $outgenetab; # may be undef
  
  if(defined $MERGE) { 
    # -merge : merge all table type, or -merge=covtab,chrtab,readids,.. which type only
    unless($MERGE=~/\w/){ $MERGE=($only_genetab)?"genetab":"covtab,chrtab,genetab,readids"; }# all types
    if($MERGE=~/gene/ and not $outgenetab =~ m/\w/) { # fix if genetab.pt* but no -genetab or -cdsbam flag
       $outgenetab= $outtab; $outgenetab=~s/\.\w+$//; $outgenetab.=".genetab"; 
       #? or fail: unless(-s "$outgenetab.pt1") { $MERGE=~s/genetab//; }
    } 
    
    mergeparts("covtab", $outtab, @ARGV) if($MERGE=~/cov/); 
    mergeparts("chrtab", $outchrtab, @ARGV) if($MERGE=~/chr/); # -merge assume chrtab=output name, was outchrpart => .chrtabout bug
    mergeparts("genetab", $outgenetab, @ARGV) if($MERGE=~/gene/);   
    #o: mergeparts("readids", $rdidtab || "$iname.readids", @ARGV) if($MERGE=~/readid/);    
    exit;
  }
    
  my $topline= "#sam2covtab par: minident=$MIN_IDENT, mindupid=$MIN_DUPIDENT, softclip=$ALLOW_SOFTCLIP, BIN=$BN, part=$icpu/$ncpu \n";
  $topline =~ s/BIN/hetz.ALN_SUB_NMI=off, BIN/ unless($ALN_SUB_NMI);
  warn $topline  if($debug);
  if($ncpu>1 and $icpu>=$ncpu) { die "# err: icpu >= ncpu : -i $icpu, -ncpu $ncpu"; }

=item UPD21AUG20: readReadIds() now is BIG-MEM-PIG
 gone...
 
=item GeneIdsInBam is obsolete
  GeneIdsInBam replaces readReadIDs mem pig .. 
  
=cut

  
  ($haslen)= readChrtab($inchrtab, $inbam); # NOTE: this reads chr-len from inbam headers @SQ ID: LEN: unless inchrtab has it

  #UPD7b: always call read_idclass, check kMAXIDCLASS and sam hdr IDs ~ $CRTPAT 
  ($nidclass,$idclassh,$idclasslist)= read_idclass($crclassf,0,\%crlen); # cds,te class by id, $idclass->{id} = class
  $USE_CDSBAMCLASS=0 unless($nidclass and $hasCdsBam); #UPD21DEC
  
  # BUG if haslen==0 but ncpu,icpu>1 .. reset one, need part tags
  if(($ncpu>1 and $icpu<$ncpu)) { # ($hasread or $haslen) and 
  
    $outpart="$outtab.pt$icpu";
    $outchrpart= "$outchrpart.pt$icpu";
    if($outgenepart) { $outgenepart= "$outgenepart.pt$icpu"; } # outgenepart may be undef
    $rdidtab .= ".pt$icpu" if($SAVE_RDIDS and $rdidtab); # rdidtabpart ?
    warn "# output part $icpu/$ncpu nchr=$NOK to $outpart\n" if($debug);
  } 

  ## inbam must be read-ordered, ie all multimaps together, only 1 of paired reads counted now

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
  
  if($nidclass and not $GIDPREFIX) { getGidPrefix(); } # GIDPREFIX from %$idclassh ids
  if($hasCdsBam) { # == ($cdsbam and -s $cdsbam)  # UPD21AUG27: minalign/MINALN different for cds
    my $IGENE=0; my $skipcdsid=0; # note NCDS == IGENE
    open(CDSBAM,"samtools view -H $cdsbam |");     # read hdr geneids, len
    while(<CDSBAM>) {
      if(/^\@SQ/){
        my($gid)= m/SN:(\S+)/; my($crlen)= (m/LN:(\d+)/)?$1:0;
        if($gid and $crlen > 0){ 
          my $ogid=$gid; $gid= cleangid($gid); 
          # UPD21DEC.18: reinstate idclass option to classify CDSBAM geneids: TE,UNK,0 instead of CDS
          my $cdsread=0;
          if($USE_CDSBAMCLASS){ # UPD21DEC
            my($crt)= $idclassh->{$ogid} || $idclassh->{$gid} || 0; 
            $cdsread |= 1 if($crt =~ /CDS/i);
            $cdsread |= 2 if($crt =~ /TE/i);
            $cdsread |= 4 if($crt =~ /UNK/i);
            unless($cdsread){ $skipcdsid++; next; } #??? skip silently?
          } else {
            $cdsread=1; # default unless classified
          }
          my $igene= ++$IGENE; # or use ++$NCDS here
          $geneindex{$gid}= $igene;
          $agenesize[$igene]= $crlen; # $cdslen{$cr}= $crlen; 
          $ageneclass[$igene]= $cdsread;# UPD21DEC
          my @gbins;
          for( my $i=100; $i<$crlen; $i+=100) { push @gbins, $i; } #? reuse $BN here or other? readsize?
          if($crlen < 140) { push @gbins, 70; } # short ones, want this?
          elsif($crlen < 180) { push @gbins, 140; }
          elsif($crlen < 220) { push @gbins, 170; }
          $agenebins[$igene]= join(' ',@gbins);
          $NCDS++; 
          }
      }
    } close(CDSBAM);
    warn "# CDSBAM nid=$NCDS, nskip=$skipcdsid in $cdsbam\n" if($debug);
   
    $stopt =~ s/ \-m $MINALN// if($STMINFILT);
    warn "# samtools view $stopt $cdsbam\n" if($debug);
    open(CDSBAM,"samtools view $stopt $cdsbam |") or die "ERR: samtools view $stopt $cdsbam";
  }
  
  my ($ntmap)= readFilter(*IN,$STMINFILT);  #? add *CDSBAM
  close(IN); if($hasCdsBam != 0) { close(CDSBAM); }
  
  putcovtab( $outpart);
  putchrtab( $outchrpart); 

  putGenetab8j($outgenepart) if($SAVE_GENECOV); # GENECOVBINS==5, simplified, merge-able UPD21AUG
#  putGenetab8i($outgenepart) if($SAVE_GENECOV); # GENECOVBINS==5, simplified, merge-able UPD21AUG


# END
#-------------------------

sub _max{ return($_[0] < $_[1])?$_[1]:$_[0]; }
sub _min{ return($_[0] > $_[1])?$_[1]:$_[0]; }
  
sub getGidPrefix {
  if($nidclass and not $GIDPREFIX) { # ,$idclassh
    my @gids=sort keys %{$idclassh};
    my $nh=int(@gids/2); my($da,$db,$dc)= @gids[0,$nh,-1];
    my($dl)= _min(length($da), _min(length($db),length($dc)));
    while($dl > 2) {
    my $ip=0; 
    for(my $i= $dl-2; $i--; $i>1) {
      if(substr($da,0,$i) eq substr($db,0,$i) and substr($da,0,$i) eq substr($dc,0,$i)) {
        $ip=$i; last;
      }
    }
    if($ip>0){ 
      my $gidp= substr($da,0,$ip); my @nop= grep{ not m/^$gidp/ } @gids; 
      unless(@nop){ $GIDPREFIX= $gidp; last; } 
    }
    $dl= $ip;
    }
  }
  return $GIDPREFIX;
}
  

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


  # FIXME: dang ':' in chrid bug, removed from crcovat, but orig in .covtab
  # dapmag19sk_LG7 in genetab, dapmag19sk:LG7 in covtab
  our($ncrcovat,%crcovat)=(0);
  sub read_crcovat { 
    my($crcovtab)= @_; my $nr=0;
    open(CCT,$crcovtab) or return 0;
    while(<CCT>){
      next if(/^\W/); my @v=split; my($cr,$cb,@acov)=@v[0,1,5,6,7];
      $cr =~ s/:/_/g; # cr id cant have : here, FIXME covtab contains orig cr:id, genetab not
      if($crcovat{$cr}{$cb}){ $nr++; $crcovat{$cr}{$cb}= join("\t",@acov); }
    } close(CCT);
    return($nr);
  }
  
  # table merge is adding nums in columns, need to read each part, keep all rows, 1st col = id key
  # there may be diff id keys in parts, tab is hash on ids x @val cols
  my %tab=(); my($nrow,$ncol,$ntin,$nid,$nidmiss,$it)=(0) x 9;
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
        if(TRY_GENELOC >= 2) { # up to TRY_GENELOC == 3?
          my($cc,$cr,$cb)=split":",$crcb; $ncrcovat++; $crcovat{$cr}{$cb}++;
        }
        # UPD21SEP21: genetab < insert covtab data for crcb: 
        # a.collect crcb hash, crcovat{cr}{cb}=1 or push @genecrat{cr}{cb}, gid? gid,gb ?
        # b. read covtab: cov{cr}{cb}=vals[5,6,7] == crcov{$cr}{$cb}[3,4,5] 
        # c. write genetab inserting cov{cr}{cb} cols
        } 
      else { for my $i (0..$#vals) { $tab{$id}{$cb}[$i] += $vals[$i]; $ncol=$i if($i>$ncol);  } }
      $ntin++;
    } close(IN);
  }
  
  # GENECOVBINS == 4, -merge=genetab, want median of chr,pos vals .. calc here?
  my($ncrcovread,$ncrcovset)=(0,0);
  if(TRY_GENELOC >= 2 and $ncrcovat>0) { #? and $dolocgene
    ($ncrcovread)= read_crcovat($outtab); # FIXME outtab == global chr.covtab file
  }
  
  open(OUT,">$ofile"); 
  print OUT @thead if(@thead);
  for my $id (sort keys %tab) { 
    # sort genetab col2: readgene,sumgene at top before 100,200,.. ok?
    my $nic=0;
    for my $cb (sort{ $a <=> $b } keys %{$tab{$id}} ) {
      my @vals= @{$tab{$id}{$cb}}; 
      
      if($ncrcovread and $cb =~ m/cr:/) { # cb == genepos\tcr:crid:crpos
        my($gpos,$cc,$cr,$crb)= split/[\t:]/, $cb;
        if(my $covs= $crcovat{$cr}{$crb}) {
          unshift @vals, $covs; $ncrcovset++; #?? ok
        } else {
          # ?? skip cr: output if no crcovat val
        }
      }
      
#       if($mgeneids) {
#         my $gids= $vals[-1]; # always last col? 
#         if($gids =~ /\w\w/){ my %gids= map{ $_ => 1 } split ",",$gids;  
#           $gids=join",",sort keys %gids; $vals[-1]=$gids; #x push @vals,$gids; 
#           }
#       }

      if($mtype == 1) { print OUT join("\t",$id,@vals)."\n"; } # no cb for mtype == 1
      else { print OUT join("\t",$id,$cb,@vals)."\n"; } # ok for mtype=3 where id == gid\tcrid
      $nrow++; $nic++; 
      }
    $nid++; $nidmiss++ if($nic == 0);

    if($dolocgene) {
      my($nl,@oloc)= get_locgene($id);
      if($nl){ for my $oloc (@oloc) { print OUT join("\t",$id,"locgene",$oloc)."\n"; } }
    }

  }
  print OUT @tcomm if(@tcomm);
  close(OUT);

  $nid="$nid,idmiss=$nidmiss" if($nidmiss>0);
  warn "#sam2covtab merge: nparts=$npt, nrows=$nrow, ncols=$ncol, ntin=$ntin, nid=$nid to $ofile\n";
  return($nrow);
}

sub midbin{ my($cb,$ce)=@_; return int( ($cb+$ce)/(2*$BN) ); }

sub getCdsRead {  # == getCdsRead_UPD
  # ( $cdsread, $geneids, $genebins)= getCdsRead_UPD( $rdid, $iid);
  my( $rdid, $iid)= @_; # , $crfl,$cr,$cb from chr.bam row cols[1..4]
  my( $cdsread, $geneids, $genebins)=(0) x 9;
  my( $moredata,$cdsrow)=($hasCdsBam > 0,''); 
  
  return ( $cdsread, $geneids, $genebins) unless($moredata);
  if($lastcdsrow) { $cdsrow= $lastcdsrow; $lastcdsrow=''; $lcdsreadid=0; $cdsiid--; }
  else { $cdsrow= <CDSBAM>; unless($cdsrow){ $moredata=0; $hasCdsBam= -1;  } } 
  # my $fixidPAIRMAP=0; # use from readFilter
  
  my @thisrows=();
  while($moredata) {
    # my($drdid,$dfl,$dgid,$dgb,$dgx,$dcig,@dsamx)= split" ",$cdsrow;    
    my($drdid,$dfl)= split" ",$cdsrow;    

    #UPD22MAR: use SAMFLAG_read2 NOT SAMFLAG_rev to fix srid from pairmap, same id for fwd/rev reads needs changing *** other gnodes 2
    if($dfl & SAMFLAG_read2) { $drdid .= "/2" unless($drdid =~ m,/2,); } # always add to _rev id dont need:  $isPAIRMAP
    
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
        my($erdid,$erfl)= split" ",$_;                 
        if($erfl & SAMFLAG_read2) { $erdid .= "/2" unless($erdid =~ m,/2,); }
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
    # $cdsread=1; #UPD21dec dont set till have one
    my $firstm= 1;
    my(@gbins,@gids,%gids);
    my $minalncds= $MINALNCDS || MIN_MINALN; # $MINALN || int( $readlen * $MIN_IDENT);  
    for my $ro (@thisrows) {
      my($grdid,$gfl,$gid,$gb,$gx,$gcig)= split" ",$ro;
      my $alen=0; while($gcig =~ m/(\d+)M/g) { $alen += $1; }
      my $cend=$alen; while($gcig =~ m/(\d+)[DN]/g) { $cend += $1; } # add to gend bin
      if($ALN_SUB_NMI and $ro =~ m/NM:i:(\d+)/){ $alen -= $1; }
      
      last if($alen < $minalncds); # use GENEMIN_LAST way
       
      $gid= cleangid($gid); 
      if( $gids{$gid}++) { next; } # ignore dup maps to same gene.cds
      my $igene= $geneindex{$gid} or next;
      push @gids, $igene; #  $gid; # replace w/ igene?
      
      # NOTE: this way can return @gids but no/less @gbins
      # agenebins[ig] == 100,200,300 .. but for shorties w/ 70/140 start/end
      my @allbins= split " ", $agenebins[$igene]; # this is All gene bins, pick overlaps to this read
      my @overbins;
      my $ge= $gb + $cend;
      # FIXME here? allbins are base-position, want gbins as array index [ 0..nbin]
      # for my $ab (split" ",$allbins)
      
# if($DEBUG_CRCB == 1) { # NOT helpful
#       # pick gbin 1 only this read overlaps most? what of longish reads that fully span 2+ 100b/BN bins?
#       my $gm= $BN * midbin($gb,$ge); # midbin() returns BN bin val, $BN* for base val here
#       my $lab=0;
#       for my $i (0..$#allbins) {
#         my $ab= $allbins[$i];
#         if($gm >= $lab and $gm <= $ab) {  push @overbins, $i; last; } 
#         $lab= $ab;
#       }
#       
# } else {
      for my $i (0..$#allbins) {
        my $ab= $allbins[$i];
        if($ab >= $gb and $ab <= $ge) { push @overbins, $i; } # i-array not ab-base
      }
# }      
      if(@overbins) {
        my $gbins= join(" ",@overbins);
        push @gbins, "$igene:$gbins";  # eg == 1:1 2 3, 222:2 3, 321:3 
      }
      
      if($firstm) { $minalncds= $MIN_DUPIDENT * $alen; $firstm=0; } 
           
      ##  gb..gb+cend is span on chr/gene seq that revcomp(rdseq) aligns to, so target seq is always fwd
      # my $gbin = int($gb/$BN); 
      # my $gebin= ($cend < $BN) ? $gbin : int(($gb+$cend)/$BN);   # alen should be cend from cigar parse
      # for(my $ib=$gbin; $ib<=$gebin; $ib++) { push @gbins, "$gid:$ib"; }
    }

  #UPD21DEC, use idclassh->{gid} =~ /CDS|TE|UNK/ to set type flags of cdsread 1|2|4
  ## cdsread use in putShortmap/readFilt:    $cds1= ($cdsread & 1); $cds2= ($cdsread & 2); $cds3= ($cdsread & 4);
  if($USE_CDSBAMCLASS){ # UPD21DEC
    $cdsread=0; #UPD21DEC, set  1|2|4 class bits, may have all/none
    for my $igene (@gids) { $cdsread |= $ageneclass[$igene]; } # note @gids is geneindex igene
  } else {
    $cdsread=1 if(@gids);
  } 
  
use constant CDSRETREF => 1;    
if(CDSRETREF) {
    # FIXME: change ret to \@gids, \@gbins ; dont need strings
    $geneids= \@gids;
    $genebins= \@gbins;
} else {
    $geneids= join",",@gids; # sort keys %gids;
    $genebins= join",",@gbins; #?? SEE addCigarDepth: data bug does this have dup entries? should be only 1 of each gid:ib per read maps
}
  }
  
  return($cdsread,$geneids, $genebins);
}




sub readFilter {  
  my($inhand,$STMINFILT)=@_;
  # inhand == open(IN,"samtools view $inbam |") 
  my($lid, $iid, $ireadpart, $inrdlen, $nokid,$ntmap,$nmap,$no_cdsread,$crtlist)=(0) x 19;
  my @saver;
  my($cdsread,$geneids, $genebins)=(0,0,0); # UPD21AUG, move cdsread/geneids to bam rdid rows as tags: eC:i:1, eG:Z:gid1,gid2
  if($MINALNCDS < MIN_MINALN){ $MINALNCDS= MIN_MINALN; } elsif($MINALNCDS > MAX_MINALNCDS) { $MINALNCDS= MAX_MINALNCDS; } 
  #not this, needs calc below: $MINALN= MIN_MINALN if($MINALN < MIN_MINALN); 

  my $testAlign= (not $STMINFILT and not $LONGR)?1:0;
  my $fixidPAIRMAP=0;  # can auto-detect PAIRMAP from srid and flags: srid fwd == srid rev
  
  my $iskip= 0;
  while(<$inhand>) {
    my @samx= split; 
    my($rdid,$fl,$cr,$cb,$cx,$cig)=@samx;
    
    #UPD22MAR: use SAMFLAG_read2 to fix srid from pairmap, same id for fwd/rev reads needs changing *** other gnodes 2
    #auto-detect: if($isPAIRMAP == -1 and ($rdid eq $lid) and ($fl & SAMFLAG_read2) and not($lastfl  & SAMFLAG_rev)) { $isPAIRMAP=1; }
    # but: if($isPAIRMAP == -1 and not($rdid eq $lid)) .. turn off ?
    # only for PAIRMAP? ($fl & SAMFLAG_read2), ie always test ids when have _rev flag
    if($fl & SAMFLAG_read2){
      if($fixidPAIRMAP > 0) { $rdid .= "/2"; }
      elsif($fixidPAIRMAP==0) { $fixidPAIRMAP=($rdid =~ m,/2,)?-1:1; $rdid .= "/2" if($fixidPAIRMAP > 0); }
    }
    # if($isPAIRMAP and $fl & SAMFLAG_read2) { #?? dont need isPAIRMAP opt, always add /2 to _rev id
    #   if($fixidPAIRMAP==0) { $fixidPAIRMAP=($rdid =~ m,/2,)?-1:1; }
    #   if($fixidPAIRMAP >0) { $rdid .= "/2"; }
    # }
   
    
    if($rdid ne $lid) {
    
      if($lid and $nokid > 0) {
        if($LONGR){  
        ($nmap,$crtlist)=  putLongmap( $lid, $inrdlen, $nokid, \@saver, $cdsread,$geneids, $genebins); 
        } else {
        ($nmap,$crtlist)=  putShortmap( $lid, $inrdlen, $nokid, \@saver, $cdsread,$geneids, $genebins); 
        }
        $ntmap += $nmap;
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
        if($inrdlen > $MAXSHORTRD) { # FIXME; my $MAXSHORTRD= 500;
          $LONGR=1; $testAlign= 0;
        }
        
        $sum_readlen += $inrdlen; # new global counter, ave(rdlen)= sum_readlen/n_readid
        
        if($DEBUG_RDFILT) { # CDSRETREF $geneids is ref(@) now 
          warn "#DERD: n_readid=$n_readid, ntmap=$ntmap, iid=$iid, ciid=$cdsiid, cdsids=$cdsread, icpu=$icpu\n"
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
    my($nmi)= ($ALN_SUB_NMI and m/NM:i:(\d+)/)?$1:0;    # UPD7f saver: add inreadlen if > 9
    if($fl & 0x04){ $n_nomap++; $badmap=1; } ## next; #? count n_readid here
    elsif($testAlign) { #rfQUICKEN:  speed up some, skip low qual quick, MINALN set in putmap2b() 
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
      my $lgid=0;    
      my @glocs= (CDSRETREF) ? @$genebins : split(",",$genebins); #UPD21SEP now  igene:100,200,300..
      for my $gloc (@glocs) {
        # UPD21SEP
        my($igene,$gbins)= split":",$gloc; # no igene dups in glocs
        my @gbins= split" ",$gbins;
        for my $gi (@gbins) {
          $agenebincov[$igene][$gi]{'covz'} += 1; #??
        }
        $agenecov[$igene]{'badmap'} += 1; # badmap? or covz

        ## old
        # my($gid,$gi)=split":",$gloc;
        # $covgenebin{$gid}{$gi}{'covz'} += 1;
        # unless($gid eq $lgid) {
        # $covgene{$gid}{'badmap'} += 1; # maybe add rdlen? not $bdepth; 
        # }
        
        $lgid=$igene;
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

  my @geneids; if($geneids) { if(CDSRETREF){ @geneids= @$geneids; } else { @geneids=split",",$geneids } }
  my @genebins; if($genebins) { if(CDSRETREF){ @genebins= @$genebins; } else { @genebins= split",",$genebins; } }
  # if($SAVE_GENECOV and $geneids){ my @gids= split",",$geneids; $geneids= \@gids; } #UPD21AUG
  # if($SAVE_GENECOV and $genebins){ my @gids= split",",$genebins; $genebins= \@gids; } #UPD21AUG
  
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

# FIXME UPD21SEP
  if(GENECOVBINS and ref($geneids) and $nmap>0) {
    for my $gid (@$geneids){ 
#       $covgene{$gid}{'nreads'}++;  # == mread above in crdtype
#       $covgene{$gid}{'nmaps'} += $nmap;  
    }
  } 
        

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

#UPD21SEP : FIXME > acovgene[igene] ...
  if(GENECOVBINS and $geneids) {
    #?? need align/readlen factor for these?
    my $pal= ($cend<1) ? 0: $alen/$cend; my $prv= $pal * $rv;
    for my $gid (@$geneids){  
#       $covgene{$gid}{'allcovm'} += $prv;
#       $covgene{$gid}{'allcovt'} += $pal; 
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

# sub addCigarDepth_OLD { 
#   my( $addToBins, $aligndepth, $cr,$cb,$cigar,$rdlen,$nmap,$cdsread,$geneids, $genebins)=@_;
#   our(%crmax,%cov,%covt,%covu,%allcov,%allcovt,%allcovu); # ,@savemap,@saveid
#   my($alen,$lenc,$cend,$nlen,$indel,$softclip)= (0) x 9;
#   # cend == ci, alen == cm, bi == w # $cm+=$w; 
#  
#   my @thisdepth= (0) x $rdlen;
#   while($cigar =~ m/(\d+)([A-Z])/g) { 
#     my($bi,$bt)=($1,$2);
#      
#     if($bt eq 'H') { 
#       $lenc += $bi; $bi=0; 
#     } elsif($bt eq 'S') {
#       $softclip += $bi; $lenc += $bi; $bi=0; 
#     } elsif($bt eq 'M') {
#       for(my $i=0;$i<$bi;$i++){ $thisdepth[$cend+$i]=1; } #  $aligndepth->[$cend+$i]++; 
#       $alen+= $bi; $lenc += $bi;  $cend += $bi;  
#     } elsif($bt eq 'N') { 
#       $n_intron++;  $cend += $bi; # cend not changed for HISP; 
#     } elsif($bt eq 'I') {
#       $n_insert += $bi; #UPD7f change for base count eq n_mismatch += nmi
#       $lenc += $bi; $bi=0; $indel++;
#     } elsif($bt eq 'D') {  # P also?
#       $n_delete += $bi; #UPD7f change for base count
#       $cend += $bi; $indel++;
#     } elsif($bt eq 'P') {  # what is P now?
#       $cend += $bi;
#     } else {
#       # unknown here, what? record?
#     }
#   }
# 
#   if($addToBins) {
#     # @$aligndepth should be sum of thisdepth over mappings, use for rvcov depth
# 
#     my $rv=1; my $uv=0; my $rty='zero';
#     if($nmap==1) { $uv=1; $rty='uniq' } 
#     elsif($nmap>1) { $rv= 1.0/$nmap; $rty='mult'; }   
# 
#     
#     my $be= int(($cb+$cend)/$BN);  
#     $crdtype{$cr}{$rty} += $nmap;
#     $crmax{$cr}=$be if($be>$crmax{$cr});
#     my $bdepth= 1.0/$BN; # map2b() equiv $bspan/$BN, allcovt counts reads at this (cr,cb+i) base / BINSIZE
# 
#     for(my $i=0, my $ifirst=1; $i<=$cend; $i++, $ifirst=0) {
#       next unless($thisdepth[$i]); # 0 or 1
#       my $bi= int( ($cb+$i)/$BN );
#       
#       my $adepth= $aligndepth->[$i]||1; # 1..nmap
#       my ($rvd,$uvd) = (0,0); # rdv = total rmap / dup maps
#       if($adepth > 1) { 
#         $rvd= $bdepth/$adepth; # was (1/adepth)/BN
#       } elsif($adepth == 1){ 
#         $uvd= $rvd= $bdepth;
#         $allcovu{$cr}[$bi] += $uvd;   # uniq read depth
#       } 
#       $allcov{$cr}[$bi]  += $rvd; 
#       $allcovt{$cr}[$bi] += $bdepth;  
# 
#         # addToBins call for each map of rd, but adepth includes each, ie N x N
#         # .. no, this is ok: each map of rd has diff cr,cb .. but this new calc has some much larger total-dupl vals,
#         # .. ** map2b() sets bdepth=bspan/BN, ignoring nmap/adepth .. want that again? probably
#         # .. where rvd is near same
# 
# 
# =item FIXME gene-loc cov counts too high      
#     
#   FIXME: gi bin loc may not overlap cb+i .. ***
#   problem is v8i change to gene-locs vs chr-locs has nearly doubled  cov[tm] vals
#     -- may be due to read partly over bins are adding uncovered bin extra, ie ~50% on ave?
#     -- maybe pick 1 gene bin only, bin at ave(bb,be)? in getCdsRead
#     
#   sam2covtab8f
#    grep g10178782t1  drosmel6ref_chr_SRR11460802ab_test8f.genetab
#   dromel6:g10178782t1	NT_033779.5	16705500	3	3	0
#   dromel6:g10178782t1	NT_033779.5	16705600	23	23	0
#   dromel6:g10178782t1	NT_033779.5	16705700	17	17	0
#   dromel6:g10178782t1	NT_033779.5	16705800	3	3	0
#   dromel6:g10178782t1	sumgene	1	43	43	0
#   
#   New sam2covtab8i
#   grep g10178782t1  test8if/drosmel6ref_chr_SRR11460802ab_test8id.genetab
#   dromel6:g10178782t1	readgene	43	43	43	0 << same as old sumgene
#   dromel6:g10178782t1	sumgene	47	47	47	0
#   dromel6:g10178782t1	100	36	36	36	0   << ~50% higher than old max
#   dromel6:g10178782t1	200	29	29	29	0   <<  higher
# 
#   # ?? change genebins val to genespan base level?  gid:gb:ge
#   # ?? need to match read-index here: 0..cend
#   # if($cbi >= $gb and $cbi <= $ge) { add cov }
#   # NO: next unless($gi == $bi); gi == gene-loc, bi == chr-loc
# 
# getCdsRead:
#       FIX==0: for(ib=gbin; ib<=gebin; ib++) push @gbins, gid:ib;
#       if(GENEBINFIX == 1) {
#         my $ge= $gb+$cend;
#         push @gbins, "$gid:$gb:$ge";
#       } elsif(GENEBINFIX == 2) {
#         my $nbin = int(0.3 + $alen/$BN); #? 0.49 maybe
#         my $gbin = int($gb/$BN); 
#         for(my $i=0; $i<$nbin; $i++) { my $ib=$gbin + $i*$BN; push @gbins, "$gid:$ib"; }
#       } 
#       
# =cut
# 
#       if(GENECOVBINS >=5 and $genebins) {
#         my($lgid)=(0);
# 
#         for my $gloc (@$genebins) { # for now same as 8a
#           my($gid,$gi)=split":",$gloc;
# if(GENEBINFIX == 1) {
#           my($gidx,$gb,$ge)=split":",$gloc; $ge=$gb if($ge<$gb);
#           if($gb + $i <= $ge) {
#           my $di= int( ($gb+$i)/$BN ); # like bi?
#           $covgenebin{$gid}{$di}{'covm'} += $rvd; # or? $bdepth; #? maybe want rdv, udv as per allcov ?
#           $covgenebin{$gid}{$di}{'covt'} += $bdepth;
#           $covgenebin{$gid}{$di}{'covu'} += $uvd;
#           }         
# } else { # orig == 0 and GENEBINFIX == 2         
#           $covgenebin{$gid}{$gi}{'covm'} += $rvd; # or? $bdepth; #? maybe want rdv, udv as per allcov ?
#           $covgenebin{$gid}{$gi}{'covt'} += $bdepth;
#           $covgenebin{$gid}{$gi}{'covu'} += $uvd;
# }          
#           if($gid ne $lgid) { 
#           $covgene{$gid}{'allcovm'} += $rvd;
#           $covgene{$gid}{'allcovu'} += $uvd;
#           $covgene{$gid}{'allcovt'} += $bdepth; # was 'alldepth'
#           }
#           $lgid=$gid;
#         }
#       }  
#       
#       if($cdsread > 0) {
#         $covt{$cr}[$bi] += $bdepth if($cdsread & 1); # 1st col
#         $cov{$cr}[$bi]  += $bdepth if($cdsread & 2); # 2nd col
#         $covu{$cr}[$bi] += $bdepth if($cdsread & 4); # 3rd col .. more cols? later
#       }
#       
#     }
# 
#   }
#   
#   # FIXME: addCigarDepth returns cend == readlen, 0-base not cb-base, change from prior
#   return($alen,$lenc,$cend,$softclip, \@thisdepth);#?
# }


# sub addCigarDepth { # == addCigarDepth_NU 21AUG30 .. seems okay now
#   my( $addToBins, $aligndepth, $cr,$cb,$cigar,$rdlen,$nmap,$cdsread,$geneids, $genebins, $imap)=@_;
#   our(%crmax,%cov,%covt,%covu,%allcov,%allcovt,%allcovu); # ,@savemap,@saveid
#   my($alen,$lenc,$cend,$nlen,$indel,$softclip)= (0) x 9;
#   # cend == ci, alen == cm, bi == w # $cm+=$w; 
#  
#   my @thisdepth= (0) x $rdlen;
#   while($cigar =~ m/(\d+)([A-Z])/g) { 
#     my($bi,$bt)=($1,$2);
#      
#     if($bt eq 'H') { 
#       $lenc += $bi; $bi=0; 
#     } elsif($bt eq 'S') {
#       $softclip += $bi; $lenc += $bi; $bi=0; 
#     } elsif($bt eq 'M') {
#       for(my $i=0;$i<$bi;$i++){ $thisdepth[$cend+$i]=1; } #  $aligndepth->[$cend+$i]++; 
#       $alen+= $bi; $lenc += $bi;  $cend += $bi;  
#     } elsif($bt eq 'N') { 
#       $n_intron++;  $cend += $bi; # cend not changed for HISP; 
#     } elsif($bt eq 'I') {
#       $n_insert += $bi; #UPD7f change for base count eq n_mismatch += nmi
#       $lenc += $bi; $bi=0; $indel++;
#     } elsif($bt eq 'D') {  # P also?
#       $n_delete += $bi; #UPD7f change for base count
#       $cend += $bi; $indel++;
#     } elsif($bt eq 'P') {  # what is P now?
#       $cend += $bi;
#     } else {
#       # unknown here, what? record?
#     }
#   }
# 
#   # addCigarDepth_NU2 : okay now
#   if($addToBins) {
#     my $rv=1; my $uv=0; my $rty='zero';
#     if($nmap==1) { $uv=1; $rty='uniq' } 
#     elsif($nmap>1) { $rv= 1.0/$nmap; $rty='mult'; }   
# 
#     my $be= int(($cb+$cend)/$BN);  
#     $crdtype{$cr}{$rty} += $nmap;
#     $crmax{$cr}=$be if($be>$crmax{$cr});
#     
#     my($cds1,$cds2,$cds3)=(0) x 3;
#     if($cdsread > 0) {
#      $cds1= ($cdsread & 1);
#      $cds2= ($cdsread & 2);
#      $cds3= ($cdsread & 4);
#     }
#         
#     my($bi0,@biv)=(-1);
#     my $bdepth= 1.0/$BN; 
#     for(my $i=0; $i<=$cend; $i++) {
#       next unless($thisdepth[$i]); # 0 or 1      
#       my $adepth= $aligndepth->[$i]||1; # 1..nmap
#       my $bi= int( ($cb+$i)/$BN );
#       push @biv, $bi if($bi > $bi0); $bi0= $bi; # try4
#       
# #UPD21SEP : change this?
#       my ($rvd,$uvd) = (0,0); # rdv = total rmap / dup maps
#       if($adepth > 1) { 
#         $rvd= $bdepth/$adepth; # was (1/adepth)/BN
#       } elsif($adepth == 1){ 
#         $uvd= $rvd= $bdepth;
#         $allcovu{$cr}[$bi] += $uvd;   # uniq read depth
#       } 
#       $allcovt{$cr}[$bi] += $bdepth;  
#       $allcov{$cr}[$bi]  += $rvd; 
#       
#       if($cdsread > 0) {
#         # these are all total depth counts, where CDS/TE/UNK exist
#         $covt{$cr}[$bi] += $bdepth if($cds1); # 1st col
#         $cov{$cr}[$bi]  += $bdepth if($cds2); # 2nd col
#         $covu{$cr}[$bi] += $bdepth if($cds3); # 3rd col .. more cols? later
#               
#       } # has cdsread
#       
#     } # read align span 0..cend
# 
# 
# =item UPD21SEP Obsolete TRY_GBINLIM=4,  %covgenebin   
#   
# # UPD21SEP add chr-loci counts agenecovloc?
# 
#     next unless($imap < 10);
#     my $nbiv=@biv; my $bi= $biv[ int($nbiv/2) ];
#     my($acovt,$acovm,$acovu)= ($allcovt{$cr}[$bi], $allcov{$cr}[$bi], $allcovu{$cr}[$bi]);
#     (my $ccr=$cr) =~ s/:/_/g; # cr id cant have : here
#     my $crbi= "$ccr:$bi";  #??
#     
#     for my $gloc (@$genebins) {
#       my($igene,$gbins)= split":",$gloc; # no igene dups in glocs
#       my @gbins= split" ",$gbins;
#       for my $gi (@gbins) {
#         $agenecrloc[$igene][$gi]{$crbi} = $acovm; #?? = nmaps
#       }
#     }
#       
# #UPD21SEP : drop this way
# 
#     if(TRY_GBINLIM >= 4 and $DEBUG_GBINLIM and $cdsread > 0 and $genebins ) # GENECOVBINS >=5 and 
#     { 
#       my($lgid,$nbin)=(0,0);
#       my $nbiv=@biv; my $bi= $biv[ int($nbiv/2) ];
#       my($acovt,$acovm,$acovu)= ($allcovt{$cr}[$bi], $allcov{$cr}[$bi], $allcovu{$cr}[$bi]);
# 
#       my $crbi = 0;
#       if( TRY_GENELOC == 2 ) { 
#         (my $ccr=$cr) =~ s/:/_/g; # cr id cant have : here
#         my $mi= 10 * int(0.5 + $bi/10); #? is change from bi a problem? count over larger span ?
#         $crbi= "$ccr:$mi";  #  bigger bin? int($bi/10) == 1000
#       }
#       
#       for my $gloc (@$genebins) { # for now same as 8a
#         my($gid,$gi)=split":",$gloc;
#         if($gid ne $lgid) { # pick 1
#         $covgenebin{$gid}{$gi}{'covt'} = $acovt;
#         $covgenebin{$gid}{$gi}{'covm'} = $acovm; 
#         $covgenebin{$gid}{$gi}{'covu'} = $acovu;
#         
#         if( TRY_GENELOC == 2 ) { #TEST 21SEP04
#         if($imap < 10) { # ? use $imap to limit excess crbi ?
#         $covgenecrloc{$gid}{$gi}{$crbi} = $acovm; # +=1 for nreads?   $acovm == covdepth ? NOT += acovm
#         }
#         # ? check for overloaded {crbi}, limit to 9? vals per gid.gi
#         
#         }
#         
#         $covgene{$gid}{'allcovt'} = $acovt;  
#         $covgene{$gid}{'allcovm'} = $acovm;
#         $covgene{$gid}{'allcovu'} = $acovu;
#         }
#         $lgid=$gid;
#       }
#       
#     }
# =cut
#     
#   } # addToBins
# 
#   # FIXME: addCigarDepth returns cend == readlen, 0-base not cb-base, change from prior
#   return($alen,$lenc,$cend,$softclip, \@thisdepth);#?
# }

sub addCigarDepth_NUNU { # == addCigarDepth_NU 21AUG30 .. seems okay now
  my( $addToBins, $aligndepth, $cr,$cb,$cigar,$rdlen,$nmap)=@_;
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

  return($alen,$lenc,$cend,$softclip, \@thisdepth)
    unless($addToBins);

  if($addToBins) {            
    my($bi0,$ib,@biv,@covs)=(-1,-1);
    my $bdepth= 1.0/$BN; 
    for(my $i=0; $i<=$cend; $i++) {
      next unless($thisdepth[$i]); # 0 or 1      
      my $adepth= $aligndepth->[$i]||1; # 1..nmap
      my $bi= int( ($cb+$i)/$BN );
      if($bi > $bi0) { push @biv, $bi ; $bi0= $bi; $ib++; }
      
      my ($rvd,$uvd) = (0,0); # rdv = total rmap / dup maps
      if($adepth > 1) { 
        $rvd= $bdepth/$adepth; # was (1/adepth)/BN
      } elsif($adepth == 1){ 
        $uvd= $rvd= $bdepth;
      }       
      $covs[$ib]{'covt'} += $bdepth;  
      $covs[$ib]{'covm'} += $rvd; 
      $covs[$ib]{'covu'} += $uvd;   # uniq read depth
    }
    return($alen,$lenc,$cend, \@biv, \@covs); # diff ret()
  } # addToBins

  # FIXME: addCigarDepth returns cend == readlen, 0-base not cb-base, change from prior
  return($alen,$lenc,$cend,$softclip, \@thisdepth);#?
}


    
=item UPD7f putShortmap 
  
  .. uses addCigarDepth() twopass on read saveset to calc 
      1st. valid aligns among multimaps, 2nd. cov depth/valids
  .. replaces putmap2b
  
=cut

sub putShortmap { 
  my($rdid, $rdlen, $nokid, $saveset, $cdsread, $geneids, $genebins)=@_;
  our(%crmax,%cov,%covt,%covu,%allcov,%allcovt,%allcovu); # ,@savemap,@saveid

  # FIXME: getCdsRead() should return \@geneids so dont reparse array each read
  my @geneids=(); if($geneids) { if(CDSRETREF) { @geneids= @$geneids; } else {  @geneids= split",",$geneids; } }
  my @genebins=(); if($genebins) { if(CDSRETREF) { @genebins= @$genebins; } else {  @genebins= split",",$genebins;  } }
  # if($SAVE_GENECOV and $geneids){ my @gids= split",",$geneids; $geneids= \@gids; } #UPD21AUG
  # if($SAVE_GENECOV and $genebins){ my @gids= split",",$genebins; $genebins= \@gids; } #UPD21AUG

  my(@thismap,@aligndepth,%crids); @aligndepth=();
  my($nmap,$topaln,$mindupaln,$lcr)=(0) x 9; # ndi   
  my $nmapin= @$saveset;
  use constant { ADD1ST => 0, ADD2BINS => 2 };
  use constant GENEBIN => 2000; # DEBUG_GENELOC: should be option, boost default to 5k ?
  my $GENEBIN2BN=  GENEBIN / $BN;
  
  for my $saver (@$saveset) {
    my($cr,$cb,$cigar,$nmi,$sfl)= @$saver;
    
    # FIXME: addCigarDepth returns cend == readlen, 0-base not cb-base, change from prior
    # my($alen,$lenc,$cend,$softclip, $thisdepth)= 
    #   addCigarDepth(ADD1ST, [], $cr,$cb,$cigar,$rdlen,$nmapin,$cdsread,$geneids, $genebins); # nmap changes in loop
    my($alen,$lenc,$cend,$softclip, $thisdepth)= 
      addCigarDepth_NUNU(ADD1ST, [], $cr,$cb,$cigar,$rdlen,$nmapin);
      
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

  my($cds1,$cds2,$cds3)=(0) x 3;
  if($cdsread > 0) {
     $cds1= ($cdsread & 1);
     $cds2= ($cdsread & 2);
     $cds3= ($cdsread & 4);
  }

  for my $saver (@thismap) {
    my($cr,$cb,$cigar,$nmi,$sfl)= @$saver;

    # my($alen,$lenc,$cend,$softclip, $thisdepth)= 
    #  addCigarDepth(ADD2BINS, \@aligndepth, $cr,$cb,$cigar,$rdlen,$nmap,$cdsread,$geneids, $genebins, ++$imap); 

    my($alen,$lenc,$cend,$biva,$covs)= 
      addCigarDepth_NUNU(ADD2BINS, \@aligndepth, $cr,$cb,$cigar,$rdlen,$nmap);
    ++$imap;

    my $rty= ($nmap==1) ? 'uniq' : ($nmap>1) ? 'mult' : 'zero';
    $crdtype{$cr}{$rty} += $nmap;
    my $be= int(($cb+$cend)/$BN);  
    $crmax{$cr}=$be if($be>$crmax{$cr});
    
    my $nb=@$biva; my $acovm=0;
    for(my $i=0; $i<$nb; $i++) {
      my $bi= $biva->[$i]; # N: bi == int( ($cb+$i)/$BN ) for i=0..cend
      my($covt,$covm,$covu)= map{ $covs->[$i]{$_}||0 } qw( covt covm covu );  

      $allcovt{$cr}[$bi] += $covt;  
      $allcov{$cr}[$bi]  += $covm; 
      $allcovu{$cr}[$bi] += $covu;   # uniq read depth
      $acovm=  $allcov{$cr}[$bi] if($i==0); # max(acovm, new)?
      
      if($cdsread > 0) {
        # these are all total depth counts, where CDS/TE/UNK exist
        $covt{$cr}[$bi] += $covt if($cds1); # 1st col
        $cov{$cr}[$bi]  += $covt if($cds2); # 2nd col
        $covu{$cr}[$bi] += $covt if($cds3); # 3rd col .. more cols? later
      } # has cdsread
    }


    if(@genebins and $imap < 10) {
    (my $ccr=$cr) =~ s/:/_/g; # cr id cant have : here, FIXME covtab contains orig cr:id, genetab not
    my $bi= int( $cb/$BN); # use this, $DEBUG_CRCB == 0 
    
    # if($DEBUG_CRCB == 1) { # NOT helpful
    #   $bi= midbin($cb,$cb+$cend);
    # } elsif($DEBUG_CRCB == 0) {
    #   $bi= int( $cb/$BN); # was $biva->[0]; .. same val 
    # } elsif($DEBUG_CRCB == 2) {
    #   $bi= $cb; # test unbinned base point align, more or less agreement from read aligns?
    # } else {
    #   $bi= $biva->[0]; #orig
    # }
    
      # could use $ib=$BN*(1+$cb) .. tabout does that calc
      # acovm should be same as  sum($covs->[0]{'covm'}) .. try += covm here instead?
    my $crbi= "$ccr:$bi";  #??
    for my $gloc (@genebins) {
      my($igene,$gbins)= split":",$gloc; # no igene dups in glocs
      my @gbins= split" ",$gbins;
      for my $gi (@gbins) {
        $agenecrloc[$igene][$gi]{$crbi} = $acovm; #?? = nmaps
      }
    }
   }
    
# if($DEBUG_GENELOC) {
#     my $cgbin= int( 0.5 +  int(0.5 + $cb/GENEBIN) * $GENEBIN2BN); # GENEBIN == 2000 converted to BN loc
#     $gbins{$cr}{$cgbin}++; # this way sums to readgene aCovT,  split into cgbins
#     # adding += 1/nmap should sum to readgene aCovM; readgene aCovU when nmap == 1
#     #? and/or add $alen;  
# }    

    #o if($SAVE_RDIDS and $cr ne $lcr){ my($crt)= $idclassh->{$cr} || 0; $crtlist{$crt}++; }
    $lcr=$cr;
  }

  if(@geneids and $nmap>0) { # GENECOVBINS
    for my $igene (@geneids) {
      # UPD21SEP
      $agenecov[$igene]{'nreads'} ++;
      $agenecov[$igene]{'nmaps'}  += $nmap;
      $agenecov[$igene]{'nuniq'}  ++ if($nmap == 1);
    }

    for my $gloc (@genebins) { # @genebins may be empty w/ @geneids
      # UPD21SEP
      my($igene,$gbins)= split":",$gloc; # no igene dups in glocs
      my @gbins= split" ",$gbins;
      for my $gi (@gbins) {
        #FIXME: $gi == base val, want $gi/BN ? for [] array index 0..n
        $agenebincov[$igene][$gi]{'covm'} ++; #?? = nreads
        $agenebincov[$igene][$gi]{'covt'} += $nmap; #?? = nmaps
        $agenebincov[$igene][$gi]{'covu'} ++ if($nmap == 1); #??
      }
    }
    
  } 
    
## old GENECOVBINS
#     for my $gid (@$geneids){ 
#       $covgene{$gid}{'nreads'}++;  # == mread above in crdtype
#       $covgene{$gid}{'nmaps'} += $nmap;  
#       $covgene{$gid}{'nuniq'} ++ if($nmap == 1); # and $ngid == 1  #??
#       
# if($DEBUG_GENELOC) {
#       for my $cr (sort keys %gbins) {
#         for my $cb (sort keys %{$gbins{$cr}}) {
#           my $vc= $gbins{$cr}{$cb}; # == readgene.aCovT
#           #? add equiv for readgene.aCovM ?
#           $covgene{$gid}{'crloc'}{$cr}{$cb} += $vc; # vc? or +1, may become huge w/ big-pig data
#         }
#       }  
# }        
#     }
    
      
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
if(UPD22JAN and $TEST_FRACTIONS) { #UPD.TEST 22JAN02 
      # rounding aCovM to int may bias high dupx spans, with largish NCPU=32, dropping cov-fractions/32 to zero
      # where merged val would be up to ~16, ie bad for cases of UCG <= 50, bad for daphmag w/ large dupx spans.
      # Test compromise 0.00 digits, should work ok in merge (sum row vals), elsewhere
      @ocols= map{ (defined $_) ? int( 0.5 + 100 * $_)/100 : 0 } ( @cdsv, @allv );
} else { # orig      
      @ocols= map{ (defined $_) ? int( 0.5 + $_) : 0 } ( @cdsv, @allv );
}     
      my $ib=$BN*(1+$i);  # < end of bin, or $ib= ($i==0) ? 1 : $BN * $i;  == start of bin 
      print OUT join("\t",$cr,$ib,@ocols)."\n";   
      } 
    }
  close(OUT); 
}

sub putGenetab8j {
  my($chrtab)= @_; # == genetab
  warn "# putgenetab $chrtab\n" if($debug);
  rename($chrtab,"$chrtab.old") if(-s $chrtab);
  open(OUT,">$chrtab"); 

  print OUT "# genetab version=$VERS\n";
  print OUT "# cols for GeneID, GPos: aCovT aCovM aCovU aCovZero = read map counts/bin as per chr.covtab\n";
  print OUT "# cols for sumgene,  1: allCovT allCovM allCovU noCov\n"; # gave=allCovT/nbins
  print OUT "# cols for readgene, 1: nrdmaps, nreads, nuniq, nomap\n";
  
  my @ocols=qw(GeneID GPos aCovT aCovM aCovU aCovZ); # aCovZ == noCov as per 8e, aCovU uniq rd cov as per other covtab 
    # adds GeneID 1st? or last? for merge better 1st
  print OUT "#".join("\t",@ocols)."\n";
  
  my $areadlen= ($n_readid<2) ? $readlen : int(0.5 + $sum_readlen/$n_readid); # UPD7f
  my $addGIDPRE= ($GIDPREFIX and $ncutGIDPRE>0);

  my $ngene= scalar(@agenesize); # scalar keys %geneindex; # 
  my($tgmaplen,$tallmap,$tallmapm,$tallmapu,$tnmaps,$tnuniq,$tnreads,$tnomap,$tnodepth,$ttnb,$tmedcovt,$tmedcovm)=(0) x 19; 
  my(@gmed,@gave,@igeneid);
  my $igmax=0; 
  for my $gid (keys %geneindex){ my $ig= $geneindex{$gid}; $igeneid[$ig]= $gid; $igmax=$ig if($ig>$igmax); }
  for(my $igene=1; $igene <= $igmax; $igene++) { # igene=1 is 1st; ugh (@agenesize)
    my $glen= $agenesize[$igene];
    # my $allbins= $agenebins[$igene];
    my $gid= $igeneid[$igene];
    next unless($gid); #? any missing ig indices?
    my $gidout= ($addGIDPRE)? $GIDPREFIX.$gid : $gid;

    my $nreads= $agenecov[$igene]{'nreads'}||0; # merge-able count
    my $nmaps = $agenecov[$igene]{'nmaps'}||0; # merge-able count
    my $nuniq = $agenecov[$igene]{'nuniq'}||0; # merge-able count
    my $nomap = $agenecov[$igene]{'badmap'}||0;   # merge-able count  #? badmap == allcovz
    $tnmaps+=$nmaps;  $tnuniq+= $nuniq; $tnreads+=$nreads; $tnomap+=$nomap; 

    ## dont need both these sum,read count .. sum is better stat .. NO, readcount best accuracy
    print OUT join("\t",$gidout, 'readgene', $nmaps, $nreads, $nuniq, $nomap)."\n";  #upd21jul04: was off, want for match to other tabs
    #o print OUT join("\t",$gidout, 'sumgene', $alldepth, $alldepthm, $alldepthu, $nodepth)."\n";  # note: alldepth > nreads by >=10% 

    my $tnb=0; my(@covt,@covm);
    my @allbins= split" ", $agenebins[$igene]; # this is All gene bins, pick overlaps to this read
    for(my $i=0; $i<@allbins; $i++) {
      my $gi= $allbins[$i];
      my @ocols= map{ $agenebincov[$igene][$i]{$_}||0 } qw( covt covm covu covz );

if(UPD22JAN and $TEST_FRACTIONS) { #UPD.TEST 22JAN02 
      @ocols= map{ (defined $_) ? int( 0.5 + 100 * $_)/100 : 0 } @ocols; # ( $covt,$covm,$covu,$covz );
} else { # orig      
      @ocols= map{ int( 0.5 + $_) } @ocols; # ( $covt,$covm,$covu,$covz );
}
      print OUT join("\t",$gidout, $gi, @ocols)."\n";  $tnb++;
    
      if(exists $agenecrloc[$igene][$i]) {
        #  $agenecrloc[$igene][$gi]{$crbi} = $acovm; #?? = nmaps
        my @crcb= sort keys %{ $agenecrloc[$igene][$i] }; #? sort by +cv count?
        
        # insert here? join adjacent crcb locs .. this helps
        my ($lcr,$lcb,$lcrcb,$cbi)=(0) x 4; my %crdbcov=();
        for my $crcb (@crcb) {
          my $cv= $agenecrloc[$igene][$i]{$crcb}||0; 
          my($cr,$cb)=split":",$crcb;
          if( $lcr eq $cr and abs($cb - $lcb) <= 2) { 
            $crdbcov{$lcrcb} += $cv;
          } else {
            $lcb= $cb; $crdbcov{$crcb} += $cv;
          }
          $lcr=$cr; $lcrcb=$crcb;
        }
        @crcb= sort keys %crdbcov;
             
        for my $crcb (@crcb) {
          # my $cv= $agenecrloc[$igene][$i]{$crcb}||0; 
          my $cv= $crdbcov{$crcb};
if(UPD22JAN and $TEST_FRACTIONS) { #UPD.TEST 22JAN02 
          $cv=int( 0.5 + 100 * $cv)/100;
} else { # orig      
          $cv=int( 0.5 + $cv );
}
          my($cr,$cb)= split/:/,$crcb; $cb= $BN*(1+$cb); 
          print OUT join("\t",$gidout, $gi, "cr:$cr:$cb", $cv)."\n";   # needs merge fix  
        }
      
      }
    }
  }
  
  # my $Cave = ($tallmapm < 1)? 0: sprintf"%.1f",$tallmap/$tallmapm;
  my $pMiss= ($tnreads  < 1)? 0: sprintf"%.3f",100*$tnomap/$tnreads;
  
  #o print OUT join("\t",'total','sumgene',$tallmap, $tallmapm, $tallmapu, $tnodepth)."\n"; #  
  print OUT join("\t",'total','readgene', $tnmaps, $tnreads, $tnuniq, $tnomap)."\n";  #   
  close(OUT); 
  warn "#genecov: ngene=$ngene, missed=$pMiss% of nrd=$tnreads in  $chrtab\n" if($debug);   #  C.ave=$Cave,
  return($ngene);
}


=item old putGenetab8i

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
        my @crcb= sort keys %{ $covgenecrloc{$gid}{$i} }; #? sort by +cv count?
        for my $crcb (@crcb) {
          my $cv= $covgenecrloc{$gid}{$i}{$crcb}||0; $cv=int( 0.5 + $cv );
          my($cr,$cb)= split/:/,$crcb; $cb= $BN*(1+$cb); 
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

=cut


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



