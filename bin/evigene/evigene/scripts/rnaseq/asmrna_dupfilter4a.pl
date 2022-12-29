#!/usr/bin/env perl
# asmrna_dupfilter.pl

=item about
  
  Filter out hi identity duplicate transcripts from transcript assembly,
  using megablast AFTER identifying best proteins.
  
  Principle is that over-assembling transcript reads with many assembly options
  produces a subset of accurate assemblies in a superset of crappy ones.
  
  Best-protein detection is done on superset (longest orf + cds/utr quality),
  then clustered (cd-hit) to reduce to "best" subset by protein quality.
  
  However cd-hit protein filtering retains nearly identical transcripts assembled
  from same reads, but differ in protein enough to avoid cd-hit cluster.
  
  This step identifies high identity transcripts from bestcd subset, from megablast align,
  then marks/removes the hi-id fraction, keeping longest-protein/longest-tr of each
  hi-id cluster.

=item usage

 pt=myproject
 $evigene/scripts/rnaseq/asmrna_dupfilter4.pl -aasize $pt.aa.qual  -aconsensus $pt.consensus \
   -blast $pt.selfcds.blastn
   -outeqtab $pt.alntab -outclass $pt.trclass

=item inputs

  inputs: 
    protein sizes, transcript sizes, as faCount table (id, size); see aacount for gappy proteins

   tr.megablast output from:
   makeblastdb -in $pt.tr -dbtype nucl -out $pt
   blastn -db $pt -query $pt.tr -task megablast -evalue 1e-19 -perc_identity 99 -dust no -outfmt 7 -out $pt.mblast

  outputs: 
    table of same/subset transcripts, rough equiv of cd-hit-est
    bestset.tr, filtered: 1.bestsetaa.ids only, 2.remove trsame subset (like cd-hit-est, but diff methods)

=item FIXMEs

  2013.aug : IS_CDSALIGN ORIENT or is IMPORTANT : need to know when alt-cds are reversed
  	patch in sumblastpart: $or add to bspans
  	
	## FIXME asmrna_dupfilter.pl: check for ID mismatches == too many zeros in aa.count /tr.count / blastn
	## eg. tables have diff prefixes like 'litova:vel' vs 'litovavel'

	** -eqgene genemap.eqgene use is problem .. maybe fixed, need
	   converts all main class to altmap, loses main link, and
		 causes further samemap-locus loci to be created via NOMAIN patch.
	
	2013.sep:
#m2t: ERR: trclass inline:Funhe2Exy3m004019t1_G2        drop    altmap  Fungr1EG3m004332t1      99/48           0,0,pflag:3
#m2t: ERR: trclass inline:Funhe2Exy3m149508t1_G2        drop    altmap  Funhe2Eq7m074376t1      99/35           0,0,pflag:3
## these are from new asmrna_dupfilter using xxx.eqgene table, should not be in trclass ..
## kfish2evg367mixx.trclass6:296 kfish2evg367mixx.trclass8:821 kfish2evg367mixx.trclass9:821
	grep  Funhe2Exy3m149508t1_G2 kfish2evg367mixx.eqgene
	Funhe2Eq7m074376t1      Funhe2Eq7m074376t1      Funhe2Exy3m149508t1_G2/35.55    1095sc:48614-48822:.
	Funhe2Exy3m149508t1_G2  noid    Funhe2Eq7m074376t1/35.55,Funhe2E6bm091492t1/35.39       1095sc:48708-48907:.
	Funhe2E6bm091492t1      Funhe2E6bm091492t1      Funhe2Exy3m149508t1_G2/35.39    1095sc:48614-48785:.

=item FIXME 201402 drop.mains problem

  201402: 
  AAPART cutoff is problem, with UTRBAD/POOR, causing drop.mains of uniq orthologs 
  AAPART & AAMIN, AAMINBAD, AAMINPOO  need fix to prevent drop uniq/valid main/noclass genes, 
  including recover perfectfrag drops of complete/utrgood trasm for AAPART-bad drop.mains
  part of problem that shorter-AA are sometimes best trasm
  
  update for drop.mains, perfectfrag replacements, needs input list from tr2aacds 
    of perfectfrag/dups from fastanr/cdhit prior steps
      
=item revise classifier/output

  FIXME: make OUTSPANTAB=1 default, see below alnclass.sh as better/faster classifier.
  
  update: fold in tested methods for tr class assignment, to replace outclusters
  -- newer class assignment: using %ident (>98) and %align (>50%)
          from  aabugs4qual/tsaevg/cdsidentclass.sh

  shorter alignments of identity are tested as valid criteria for alternate-transcripts
  using genome-mapped transcripts, hi-ident + part-align gives *mostly* alt-tr, some paralogs,
  depending on species & freq of hi-ident paralogs.
              
  update: add refblast  input table for added classification, per
            aabugs4/aabugs4qual/tsaevg/classparts.pl
        tall4 format: Refid Trid  Bits Iden Algn Rlen Tlen, where Ref/Tr may be swapped columns

  update: add 2ndary alt-tr table, in aa.qual format?,
      from alt-tr called by trasm soft, but excluded by cdhit/other filtering as not distinct
      these all should have same geneid + tNNN suffix as primary input tr (defined in -aasize -blast trself.blast)
      ID matching only used, plus aa.qual, to decide whether to keep/drop
      - input to classifier may be 2nd OUTSPANTAB, generated from 1st pass of this and 2nd-alts.aa.qual table.

=item add UTR BAD/POOR filters

  problem now is that main class accumulates utrbad/poor genes, many w/ other utrbad alternates,
  so as more trasm sets are added to this filtering, main class (and some of others) increase w/ junk.
  
  use input aasizes/aaqual scores
  when reading cdhit/blast aligns (esp cdhit clusters), where feasible swap/replace UTRBAD top/first gene 
    if equivalent utrok gene is in cluster/align (ie. utrok/utrbad have same prot/cds size, or ok is nearly same),
    esp. if utrok is also aacomplete.    
  per evigene/scripts/prot/aabest3.sh

=item add ncRNA classing 2016.02
  
  -- asmrna_dupfilter can do 1st pass, likely want tr2genome mapping for intron measures
  -- putative ncRNA can be pulled from tr2aacds dropset:
        drop.(main|noclass), maybe drop.altmid and frag/part classes
  -- use only subset with no CDS overlap (ie not althi, maybe not altmid)
    .. but will need added full tr/cdna overlap test      
  -- maybe select from only utrbad/utrpoor subset
  -- firmest ncrna class, w/o special tests, would include introns, ie tr2genome exon parts alignments

=item add classifyFullTr 2016.05 

=item UPD1912
  -- adds PHETERO, reduces all identity classing params by 1..3 % to account for heterozygous trasm
  -- cds-stranded align now, removes -sense problem (any reason for -sense CDS aligns for locus class?)
      calling tr2aacds now +strand align only,  blastn -strand plus  self.cds

=item AACONS consensus classing 2018.05

  evigene/scripts/prot/tr2aacds_aaconsensus_item.txt
  evigene/scripts/prot/tr2aacds2c.pl 
  evigene/scripts/rnaseq/asmrna_dupfilter3c.pl

  Consensus criterion is an important addition, to retain transcripts where
  the coding sequence (protein) is supported by identical assemblies across
  assembler methods.
  
  Adam Voshall brought this to my attention, with this thesis,
  https://digitalcommons.unl.edu/computerscidiss/145
  "Consensus Ensemble Approaches Improve De Novo Transcriptome Assemblies"
  
  showing that Evigene tr2aacds (and asmrna_dupfilter3 classifier) reduces,
  drastically, the number of perfect protein matches to a reference gene
  set, versus the number of perfects from individual assemblers.  The
  reduction is due to longest-ORF often not being the "true" coding
  sequence.

  There is a puzzle here, as extensive measure by myself, others, show this
  longest-ORF criteria does return greatest average alignment, and most
  genes with high alignment, using protein or CDS aligned to reference
  genes.  The difference, by digging thru discrepancies, is mostly at the
  99+% level, where perfect means 100% identity, a harder criteria to
  reach. Longest-ORF selection by Evigene is pulling rare? assemblies of
  genes that are slightly longer (e.g. 1..10 aminos, or <=1% longer) than
  "true" reference gene.  More work to determine if these are artifacts to
  be dropped or valid, if slightly different transcripts.

  tr2aacds should retain all valid, or likely valid, coding transcripts, removing
  only truely redundant or fragment models.  To that end, AACON consensus scoring
  and classification is added.  
  
  tr2aacds will attempt to score input transcripts for CDS-consensus, using
  ID patterns that distinguish assembler sources (as added by evigene trformat.pl, 
  but customizable), or/and by SRA2Genes with its knowledge of separate assemblies
  in  trsets/.  Consensus is measured as identical coding sequences, produced by 2+
  different assembly methods (programs, or kmer setting, or data slices).
  
  "Fastanrdb all.cds > nr.cds", as now used for assembly reduction, is also used now for
  consensus detection.
  
  Input to asmrna_dupfilter3c.pl for this is a table of
    tr_ID <tab> consensus_score
  where any score > 0 means the tr_ID will be retained in okayset, suitably classed
  (generally as an alternate to longer ORF).  This is similar to how homology table
  scores are used, 
    tr_ID <tab> ref_ID <tab> homology_score

  This command usage is then 
    tr2aacds2c.pl -aconsensus[=acons_options] to enable, 
      default mode measures consensus in trset.nrcds, writes trset.consensus table
      
    asmrna_dupfilter3c.pl -aconsensus trset.consensus  -anames trset.names (blastp scores + refnames) 
     (asmrna_dupfilter is called by tr2aacds)
    
=item add cd-hit-est input alignments

  maybe replace cd-hit-aa with cd-hit-cds(est) for both aa-based and nt-based equivalence/reduction
  of trasm sets.  cd-hit-aa has problem of lumping paralogs w/ silent subsititutions, want to keep
  those in early filtering.  cd-hit-cds/est gives approx same classes of same-locus vs diff-locus as
  blastn or blat (depending on parameters).
    cd-hit-est -c 0.99 -G 0 -aS $pALIGN -l 150 -d 0 -i cacao11pub3ig.cds
        pIDENT=0.99 is good as w/ others; pALIGN=0.25 .. 0.50
        
  cacao11pub3ig.class3
  # alternate classifier, replace both cd-hit aa and mblast-tr ? 
  cacao11g: 44403 cds, 29451 loci, 7630 loci have alts,
  blastn (99% id, >=50%? align): alt class=16179 alts, 484 diff locus; diffclass: noclass=noalts; altmid=236 alts, 1156 diff
      false-alts: 484; false-loci: 236 of 44403
  cdhit25: 29893 clusters, 730 loci have alts in diff clusters; 295 diff loci in same cluster;
      false-alts: 295; false-loci: 730 of 44403
        -- about same as above self.mblast (same loci?)
  cdhit50: 30057 clusters, 841 loci alts diff cluster; 252 diff loci in same cluster
    -- cd-hit pALIGN=0.33 may be good choice;
     
=item updated for blat psl input

  blat -out=psl -t=dna -minIdentity=98 -maxIntron=0 kfish1cds.tr kfish1cds.tr kfish1cds.selfblat
  blat -out=psl -t=dna  -minScore=99  -minIdentity=99 -maxIntron=0 kfish1cds.tr kfish1cds.tr kfish1cds.selfblat99
  #  -minScore=99 or such to reduce fragment aligns, 30 default, nmatch - nmiss
  : tests show blat and megablast give nearly same results; blat can be very slow vs mblast for large tr set
  
  
=item see also
  
  evigene/scripts/prot/aaqual.sh
  evigene/scripts/makeblastscore.pl
  evigene/scripts/rnaseq/asmrna_equalgenes.pl
  
=item author
  
  don gilbert, gilbertd near indiana edu, 2012
  part of EvidentialGene, evigene/scripts/

=cut

use FindBin;
use lib ("$FindBin::Bin/../","$FindBin::Bin"); # 201405 add; assume evigene/scripts/ layout: main, prot/, rnaseq/, ...

use strict;
use warnings;
use Getopt::Long;
#add# use cdna_evigenesub; # replacing local subs 201405

# subs were in asmrna_dupfilter3_overlocsub.pm; now in asmrna_overlocusclass4a.pl
# sub classifyFullTr;
# sub overlocusclass;

use constant VERSION => '2019.12.14'; # now version 4, 

# UPD1912v4:   drop old/unused code, fixed identityclass() NOMAIN bug again; still testing
# UPD1912a small fixes + pHETERO adjust/lower all identity params by 1-3% TEST
#  '2019.08.25'; #UPD1908 class improvements from bestfiltcompare tests, lose valid althi esp, aacons-frag/part NOT maybeok
use constant UPD1912 => 1;
use constant UPD1908 => 1;# thru UPD1911
use constant AACONS => 1;  # aaconsensus: keep tr with consensus of identical aa across asm types
use constant CODEPOT1607 => 1; # add from above codepot work

# '2018.05.12'; # aaconsensus update
# '2016.07.09'; # BAD_GAPS param added, 25% default is too high
# '2016.05.27'; # classifyFullTr() / overLocusclass()
# '2016.02.11'; # MINBASEOVER/TINYALN changes
# '2014.05.17'; #'2014.02.21'; # test/fix drop.main.utrbad problems missed orthogenes
# '2013.09.07'; # 01'; # 08.09'; #08.07; 03.25; 24; # '2013.02.21';

my $debug= $ENV{debug} || 0;
my $OUTSPANTAB= 1;  # make default until replace outclusters()
my $OVSLOP=6; # FIXME: need overlapslop ~ < 0.01 of max part span ; or < 0.02..0.05 of min part span?
my $pctOVSLOP=$ENV{pctover} || 0.02; # need opt?
   $pctOVSLOP=$pctOVSLOP/100 if($pctOVSLOP > 0.9);

# classifier options
our($AAMIN,$AAPART,$AAMINBAD, $AAMINPOO, $AADUP_IDENT, $OK_CDSUTR, $BAD_CDSUTR, $noRESET_CDSUTR, $BAD_GAPS, 
    $TINYALN, $IS_CDSALIGN,$ALTFRAG,$PHIALN,$NHIALN,$PHI,$PMID,$PLOW, $PHETERO);

# 201402: these 4 AAsize cutoffs become options, need to reset to cure drop.main.utrbad bug
# also reset BAD_CDSUTR,OK_CDSUTR options to ignore input aaqual if user sets limits: noRESET_CDSUTR
# FIXME2: need adjust pCDS/utrpoor calc for aaqual "utrorf" cases that now have full mrna size but should be split
$AAMIN =$ENV{aamin}||30; #was 40;   # for aacomplete, utrok
$AAPART=$ENV{aapart}||100; # was 100; # for aapartial # 201402: THIS aapart cut is problem, drop.mains uniq orthologs here
      # 201402: update for drop.mains, perfectfrag replacements, need input list of perfectfrag/dups from fastanr/cdhit prior steps
      # 201402: AAPART & AAMIN, AAMINBAD, AAMINPOO  need fix to prevent drop uniq/valid main/noclass genes, 
      # including recover perfectfrag drops of complete/utrgood trasm for AAPART-bad drop.mains
      # part of problem that shorter-AA are sometimes best trasm
$AAMINBAD=$ENV{aaminbad}||60; #was 200;  # for utrbad class
$AAMINPOO=$ENV{aaminpoo}||60; #was 100;  # for utrpoor class
## 2014 new opts for asmrna_dupfilter2:
## tcas4/tribol beetle: 2014 ncbi has 10 prot under 40aa, 5 are NP curated/refseq; 
## export aamin=30;  export aapart=120; export aaminbad=90; export aaminpoo=60

my $OK_CDSUTR_default= $ENV{pcdspoor} ||60; ##  # too high? changes w/ cds-size
my $BAD_CDSUTR_default= $ENV{pcdsbad} ||30; ## % CDS/trlen
$OK_CDSUTR= 0;
$BAD_CDSUTR= 0; # % CDS/trlen
$noRESET_CDSUTR= 1; # flag to ignore input tqual/aaqual utrbad

my $MINUTR=$ENV{minutr}||300; # ~ 300b "fixed" average utr sizes, maybe too low, 
# see sample of org gene set pcds for small aa .. per species min-utr, ranging from 400b insect .. 700b fish .. 1000b mouse
# should adjust %cds bad for cds size, for large cds> 10k, pcds >=90%, for small cds < 300, pcds may be 33%
# should replace utrbad/poor flags with coding potential flags, adding some cds-seq cp calcs.

# classifier should also downweight qual by more than aaspan-aagaps as aagaps damage homol scores
# eg:  score = aanogapsize - aagaps ? reduce qual by ngaps
#old# $BAD_GAPS= 25;  # % gaps in AA; 160709: make this opt, 25% too high, 5% maybe best default?
## BUT review of homol sez gappy prots can still be best/only locus rep.. mostly want to reduce qual score of gappy prots.
$BAD_GAPS= $ENV{aagapmax} || $ENV{BAD_GAPS} || 10; # was 25; # need opt? 'aagapmax' ?

## add UTR BAD/POOR filters, per evigene/scripts/prot/aabest3.sh
#  if(UPD1908)  bump up AADUP_IDENT to 99.0? 
$AADUP_IDENT= (UPD1908)? 99 : 98; # %ident, option for aacluster ident drops
$PHETERO = (UPD1912) ? ($ENV{PHETERO}||0) : 0;

$TINYALN = 20; # %align, was 25; was $ENV{mina}||50; ignore less than this == MINAL
  ## UPD 2016.02: -tinyaln 35  bad (sometimes), need opts and diff test: tiny by bases not pct overlap
  ## ie. cds-exon overlap should be valid for >= tinybases, eg. 90? need to sample ref species alt exon sizes
  ## UPD 2016.02: reinstate MINALIGN, for tiny by bases not pct overlap, or reuse TINYALN for bases overlap 
  ## ?? UseTINYALNBASES only for IS_CDSALIGN ?
use constant UseTINYALNBASES => 1;   
my $TINYALNBASES= $ENV{MINALIGN}||90; # was $MINALIGN= 90; # REUSE//NOT USED NOW; use just TINYALN; change to several levels
my $MINCDS = $ENV{MINCDS} || 90; # or 3*MINAA  
my $TINYALNCDSLEN= 3*$TINYALNBASES;
my $MIN_EQGENE= $ENV{MINEQGENE}||15; # %equalgene cds-align, was 33; # global/opt..

$ALTFRAG= 0.5;
$PHIALN=$ENV{ahi} || 98; # was 99; was 65 !!;  DROP THIS; use mina?
$NHIALN=$ENV{nhi} || 9; # what? for short cds cancel few base diff < PHIALN
#  my $tinyalndiff= ($aln - $alnmax < $NHIALN) ? 1 : 0 # TEST3 add

## tr-self %identity levels to classify alt-tr
## <90% ident prob too low to class alt-tr; use pmid=95 plow=90 ??
$PHI = $ENV{phi} ||99; 
$PMID= $ENV{pmid}||90; 
$PLOW= $ENV{plow}||80;  

$IS_CDSALIGN=1; #  $ENV{cdsw}||0; # UPD1912 now default, use -noCDSALIGN or -RNAALIGN
my $IS_RNAALIGN= ! $IS_CDSALIGN; # alias noCDSALIGN
   # WHICH default?  tralign tests better than cds, but cds for other, eg. cds-cdhits
my $SORTEDALIGNTAB=0; # debug input: sort -k2,2nr -k7,7nr -k6,6n evg2anofunz4c.alntab
my $OIDFIX=undef;
my $DOOVERLOCUS=0;

my ($aasizes,$trsizes,$blatpsl,$blastab,$lastz,$bcdhit,$aablast,$aanames,$aacdhit,$outeqtab,$outclass,
    $dupids,$logfile,$head,$eqgene)= (0) x 20;
our($aconsensus, %aconsensus)= (undef); # AACONS

my $optok= GetOptions(
  "aasizes=s", \$aasizes, 
  "trsizes=s", \$trsizes, 
  "ablastab=s", \$aablast,   # this is traa-refaa.blastp
  "anames=s", \$aanames,   # variant aablast used also for naming
  "aconsensus:s", \$aconsensus,   # for AACONS
 		## FIXME: allow -ablastab to be .names table, for tr2aacds, ie decide by file name
  "acdhit=s", \$aacdhit,    # this is traa-self.cdhit.clstr

  "blastab=s", \$blastab,    # this is tr-self.blastn, -CDSALIGN for cds-self.blastn
  #UPD1912.off# "lastz=s", \$lastz, # lastz general format; was -blastz option -blast[ab] conflict **
  #UPD1912.off# "blat=s", \$blatpsl, # tr-self.blat;  other input format opts here..?  -informat=xxx
  #UPD1912.off# "bcdhit=s", \$bcdhit, #  tr-self.cdhits.clstr; 
  "dupids=s", \$dupids, #  

  #UPD1912.off# "eqmap|eqgene=s", \$eqgene,    # 130901: mapping equivalence table, of $evigene/equalgene.pl 
  
  "aligntab|outeqtab=s", \$outeqtab,  ## this should have aliases: -outalntab and -inalntab, or just -aligntab ?  
  "outclass=s", \$outclass,   # 2nd out: outclasstab
  "sortedaligntab!", \$SORTEDALIGNTAB, 

  # 201402 option update: $AAMIN,$AAPART,$AAMINBAD, $AAMINPOO
  "AAMIN=i", \$AAMIN,   
  "AAPARTMIN=i", \$AAPART,   
  "AABADMIN=i", \$AAMINBAD,   
  "AAPOOMIN=i", \$AAMINPOO,   
  "aagapmax=i", \$BAD_GAPS,   

  "TINYALN|MINALIGN=i", \$TINYALN,  # pTINYALN ? FIXME .. both pctMINALN, basesMINALN
  "pCDSOK=i", \$OK_CDSUTR,   # CDSOKUTR
  "pCDSBAD=i", \$BAD_CDSUTR, # CDSBADUTR 
  "pHeterozygosity=i", \$PHETERO, # UPD1912 
 
  "ALTFRAG|fragment=s", \$ALTFRAG,  # pALTFRAG ?
  "OUTSPANTAB!", \$OUTSPANTAB, # now fixed=1 ? drop opt
  "CDSALIGN!", \$IS_CDSALIGN,  "RNAALIGN!", \$IS_RNAALIGN, #UPD1912 alias -noCDSALIGN 
  #UPD1912.off# "overlocus!", \$DOOVERLOCUS, 
  "debug!", \$debug, 
  );


(my $evapp=$0) =~ s,^.*/,,;
warn "# EvidentialGene $evapp, VERSION ",VERSION,"\n"; # if($debug); # change to loggit() ?
#o my $hasbalign= ($blastab or $lastz or $blatpsl or $bcdhit) ? 1 : 0;
my $hasbalign= ($blastab) ? 1 : 0;

die "usage:  asmrna_dupfilter.pl -aasize=name.aa.qual  -blast=name.cdsalign.blastn  
 input cds|rna alignment: -blast=name.blastn | -aligntab=name.aligntab
 opts: 
   -AAMIN=$AAMIN
   -ablastab=name-refaa-blast.table 
   -aconsensus=name.consensus
   -pHeterozygosity=2 : use heteroz. identity reduction
  more options, see source.
" unless($optok and $aasizes and ($hasbalign or $outeqtab)); ## and $trsizes < dont need if aasizes=aa.qual

# UPD1912.off: -blat=name.blatpsl | -bcdhit=name.cdhit.clstr
# -trsize=name.tr.count : dont need w/ aa.qual table
#UPD1912 maybe move param checks to adjustClassifyParams()

$noRESET_CDSUTR=($OK_CDSUTR>0 or $BAD_CDSUTR>0)?0:1; # user option set, ignore input aaqual 201402 update
$OK_CDSUTR ||=$OK_CDSUTR_default;
$BAD_CDSUTR||=$BAD_CDSUTR_default;
$BAD_CDSUTR= $OK_CDSUTR if($BAD_CDSUTR > $OK_CDSUTR);
if($IS_RNAALIGN) { $IS_CDSALIGN= 0; } # alias

if($TINYALN>0) { $TINYALN= 100*$TINYALN if($TINYALN<1); }
$ALTFRAG= $ALTFRAG/100 if($ALTFRAG>1); # prop not percent
 
# should be our() shared vars ..
our(%aasize, %sizeval, %trsize, %aaqual, %oids, %oidsof, %cdsoff); # all from readSizes(); one hash? == sizeval
our(%aablast, %aablastref,$naabl); $naabl=0;
our(%aacluster,%aaclustermain); # globals for readAAcdhit
our(%better, %outrows, %validids, %bspans); ## %validids was %outids; 

our(%dupids, %dupfirst, $ndupdrop); $ndupdrop=0;
use constant DUPFILTER1 => 1; # tests ok
use constant DUPFILTER2 => 0; # bad?


sub adjustClassifyParams { # UPD1912
  map{ $_= $AAMIN if($AAMIN > $_); } ($AAPART, $AAMINBAD, $AAMINPOO);
  $MINCDS= 3*$AAMIN if(3*$AAMIN > $MINCDS);
  if($PHETERO > 9 or $PHETERO < 0.001) { 
    warn "# asmrna_dupfilter: ignoring pHETERO=$PHETERO out of range 0.001 to 9 \n" unless($PHETERO == 0); # loggit
    $PHETERO= 0; }
  if($PHETERO) {
    warn "# asmrna_dupfilter: pHETERO=$PHETERO reduced classify identities \n" if($debug); # loggit
    map{ $_ -= $PHETERO; } ($PHI, $PMID, $PLOW, $PHIALN, $AADUP_IDENT);
  }
}


sub MAINstub {}

my $OUTH= *STDOUT;

adjustClassifyParams(); # UPD1912

readSizes();
readDupIds($dupids)   if($dupids);
 		##  allow -ablastab to be .names table, for tr2aacds, ie decide by file name; was '.names'
if($aablast and $aablast =~ /\.name/ and not $aanames) { $aanames= $aablast; $aablast=""; }
($naabl)= readAAnametab($aanames) if($aanames); # prefer this now?
($naabl)= readAAblast($aablast)   if($aablast and not $naabl);  # one or other of aablast,aanames

readAAcdhit($aacdhit) if($aacdhit); # also correctAAcluster()

#UPD 201805 : require $aconsensus to be defined? 
# alt consfile input: (my $nrcds = $blastab) =~ s/.self.*blastn/.cds/;
use constant READ_AACONS_FROM_NRCDS => 1; # use or not alt cons input from trset.nrcds ?
my($nconsensus,$consfile)= (AACONS) ? readConsensus($aconsensus,$blastab) : (0); 

use constant EQGENE_OVERRIDES_ALN => 0;
use constant EQGENE_CHANGES_NOALN => 0; # not quite right yet, 160217, tho new eqgene data useful
## $eqgenes->{$tid}{$sid} >= $MIN_EQGENE alignment
## FIXME: need to adjust eqgenes to account for poor gmapping, 
##   paralogs map to same locus, poorly, 
##   but blastn/any align says they are different loci
## BUT also have perfect gmap eqgenes, not scored same by blastn, should be called same locus dups

my(%eqflag); # set in identityclass, output?
my($eqgenes,$nineqgene,$neqgene,$gmapqual,$neqalts,$neqexon,$eqexons)
    = ($eqgene) ? readEqualGene($eqgene) : (undef,0,0,0,0,0); # UPD1912 drop? eqgenes removed for now


# FIXME 201405, allow reuse $outeqtab even if hasbalign.. for tr2aacds restarts
if( ($ENV{useouteqtab} or $ENV{usealigntab} ) and -s $outeqtab ) { $hasbalign=0; }

my $nbalign=0;
if($hasbalign) {
  if($outeqtab) {
    rename($outeqtab,"$outeqtab.old") if( -f $outeqtab );
    open(OUTH,">$outeqtab") or die "write $outeqtab";
    $OUTH= *OUTH;
  }
  
  # allow input of OUTSPANTAB instead of regenerate?
  if($blastab) { ($nbalign)= readblasttab($blastab); }
  #o elsif($lastz) { ($nbalign)= readlastz($lastz); } # UPD1912 off
  #o elsif($bcdhit) { ($nbalign)= readcdhit($bcdhit); }
  #o elsif($blatpsl) { ($nbalign)= readblatpsl($blatpsl); } # detect from input table ??
  if($outeqtab) { close($OUTH); $OUTH=undef; } # *STDOUT ?
  warn "# readalign: nids=$nbalign to $outeqtab\n" if $debug;
}

if($outeqtab) {
  my $infile=$outeqtab;   # infile == outeqtab ? STDIN?
  my $insorted=$SORTEDALIGNTAB; # add input opt for sorted table ?? test classing bugs
  unless($infile and -f $infile) { 
    die "ERR: missing input align table -aligntab $outeqtab";
  } elsif($infile =~ /\.gz/) { 
    die "ERR: cant use gzipped align table -aligntab $outeqtab";
  }

  # >> set this BEFORE correctAAcluster from read{blast|cd|blat} : $validids{$id}
  if($aacdhit) {
    my $havevalid= scalar(%validids)?1:0;
    $havevalid= readIdsFromAlnTab($infile) unless( $havevalid); 
    correctAAcluster($havevalid); # update %aacluster,%aaclustermain
  }

  $OUTH= *STDOUT; 
  if($outclass) {
    rename($outclass,"$outclass.old") if( -f $outclass );
    open(OUTC,">$outclass") or die "write $outclass";
    $OUTH= *OUTC;
  }
  
#   if($DOOVERLOCUS) { # UPD1912 off
#   $IS_CDSALIGN=0; # always?
#   overlocusclass($OUTH,$infile,$insorted, $IS_CDSALIGN);
#   } else {  

  identityclass($OUTH,$infile,$insorted); 
#  }

  if($outclass) { close($OUTH); $OUTH=undef; } # *STDOUT ?
}


#.................................

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }
sub bint { my $b=shift; return ($b<0)?0 : ($b=~/e\+/)? int($b) : $b; }

sub openRead { # in cdna_evigenesub.pm
  my($fna, $nostdin)= @_; my($ok,$hin)= (0,undef);
  $ok= ($fna =~ /\.gz$/) ? open($hin,"gunzip -c $fna|") 
  	 : ($fna =~ /stdin|^-/ and not $nostdin) ? *STDIN 
  	 : open($hin,$fna);  
  # loggit(1,"ERR: openRead $fna") unless($ok);
	die "ERROR: openRead $fna" unless($ok);
  return ($ok,$hin);
}



sub readSizes {  # see also cdna_evigenesub.pm getAaQual()

  my($naa,$ntr,$nerr,$ok,$inh)=(0) x 10;
  if($aasizes) {  
    ## fix for aacount gaps: id,size,gaps : NOT NOW, aa.qual: id,size-gaps,gaps,..
    ## drop faCount? require aa.qual here?

    # ($ok,$inh)= openRead($aasizes,1);
    open(F,$aasizes) or die "FAIL: read $aasizes ..."; 
    # if($aasizes =~ /count|qual/) { 
    # } else { open(F,"faCount $aasizes |") or die "FAIL: faCount $aasizes ..."; } # drop this..

## FIXME 201405: Maybe replace w/ cdna_evigenesub:getAaQual() -- somewhat different %AAQUALH data
## getAaQual id => $alen,$pctcds,$acv,$aqual1
## FIXMEd for cds.qual not aa.qual, col2 == cds-size, col4=aasize,qual, col5=trsize
## FIXME 2016.05, cds.qual now *may* have "Code/Noncode,.." column after tlen, before offs, 
## use constant CODEPOT1607 

    my $iscds=0;
    while(<F>) { 
      next if(/^#/ or /^total/); 
      my @v=split; # aa.qual cols; gap is removed from alen
      my($id,$alen,$gap,$aqual,$tlen,$offs,$oids);
      my $codepot=0; #  Code/Nonc/Unknown col maybe there..
      if($v[5] =~ /^(Code|Noncode|Unknown)/) {
        ($id,$alen,$gap,$aqual,$tlen,$codepot,$offs,$oids)=@v; 
      } else {
        for my $i (0..$#v) { if($v[$i] =~ /^(Code|Noncode|Unknown)/) { $codepot=$v[$i]; splice(@v,$i,1); last; } }
        ($id,$alen,$gap,$aqual,$tlen,$offs,$oids)= @v;  # ,offs,oids may be missing. and Code/Nonc
      }
      $offs||=0; $oids||=0;
      
      my $cdlen=0;
      unless($alen =~ /^\d/) { $nerr++; next; }
      if($aqual) { 
        my($aqlen)= $aqual=~m/(\d+)/; 
        if($aqlen < $alen and $aqlen*3 >= $alen-3) { $cdlen=$alen; $alen= int($alen/3); $iscds++; } #cdslen not alen
        elsif($aqlen == 0 and $tlen == 0 and $iscds) { $cdlen=$alen; $alen= int($alen/3);  } # datatable bug
        else { $cdlen=$alen*3; }
        $aqual .= "-gapbad" if($gap>0 and (100*$gap/($alen+$gap) > $BAD_GAPS)); # add qual flag for gap/(alen+gap) > MAXGAP:  
        $aaqual{$id}= $aqual; 
      }
      if($tlen =~ /^\d/) { $trsize{$id}= $tlen; $ntr++; } 
      if($offs =~ /^\d/) { $offs=~s/:.$//; $cdsoff{$id}= $offs; }  
      if($oids) { # capture, for fixups?
        my @oids=grep{ /^\w/ and not($_ eq "na" or $_ eq $id) } split",",$oids;
        for my $d (@oids) { $oids{$d}=$id; $aasize{$d}=$alen unless($aasize{$d});} 
        $oidsof{$id}=$oids if(@oids);  
      }       
      $aasize{$id}=$alen; $naa++; 
      ## all in one?
      @{$sizeval{$id}}{ qw(aasize cdsize gap aaqual trsize cdsoff oids codepot) }
        = ($alen,$cdlen,$gap,$aqual,$tlen,$offs,$oids,$codepot);

      
      } close(F); 
  }
  
  if($trsizes) {
  if($trsizes =~ /^aaqual/ or $ntr>0) { 
    # got above;now default
  } elsif($trsizes =~ /^aasize|^cdssize/) { # for tr == cds
    foreach my $id (keys %aasize) { $trsize{$id}= 3*$aasize{$id}; }
  } else {  
    ## drop faCount? expected use is aa.qual w/ trsize column
    # ($ok,$inh)= openRead($trsizes,1);
    open(F,$trsizes) or die "FAIL: read $trsizes ..."; 
    # if($trsizes =~ /count|qual/) { 
    # } else { open(F,"faCount $trsizes |") or die "FAIL: faCount  $trsizes ..."; }
    while(<F>) { next if(/^#/ or /^total/); my($id,$al)=split; $trsize{$id}=$al; $ntr++; } close(F); 
  }
  }
 
 warn "# readSizes: naa=$naa; ntr=$ntr\n" if $debug;
 return($naa,$ntr); 
}


# readConsensus($nrcds); # FIXME: read hdr of evg3weed1rcnrcd1.cds given blastab=evg3weed1rcnrcd1-self98.blastn
# FIXME2: aconsens table 1st ID mismatch w/ input_nrcd1.cds ids; 2nd way, reading input.cds agree= works right

sub readConsensus { 
  my($aconsensus,$blastab)=@_; # globals
  my($ncons,$consfile)=(0,"none");

  if(defined $aconsensus) {
    return($ncons,$consfile) if($aconsensus =~ m/^(no|0)/);  
  } else { #? require -acons or not
    return($ncons,$consfile);
  }
  
  if($aconsensus and -f $aconsensus) {
    open(F,$aconsensus); $consfile=$aconsensus;
    while(<F>){ 
      next if(/^\W/); 
      #NO# my($id,$ag)=split; $aconsensus{$id}=$ag if($ag); 
      # fixme: need @agids to match input nrcd1.cds
      my($id,$agscore,$agsrc,$agids)= split;
      if($agscore) {
        $aconsensus{$id}=$agscore; $ncons++;
        my @agb= split",",$agids;
        map{ $aconsensus{$_}=$agscore; } @agb;
      }
    } close(F);

  } elsif(READ_AACONS_FROM_NRCDS) { # turn off this test hack?
    my $nrcds = $blastab; # first hack..works
    $nrcds=~s/.self.*blastn//; $nrcds.=".cds";
    
    if(-f $nrcds) { 
      open(F,$nrcds); $consfile=$nrcds;
      while(<F>){ 
        if(/^>(\S+)/){ my $id=$1; #? read all equal ids on >hdr?
          if(my($ag)=m/agree=(\w+)/){ $aconsensus{$id}=$ag; $ncons++; } 
        }
      }  close(F); 
    }
  }
  
  warn "# readConsensus n=$ncons from $consfile\n" if($debug);
  return($ncons,$consfile);
}



# fix from fastanrdb all.cds > all_nr.cds >> hdr has cds-identical ids; should have prefiltered this
sub readDupIds {
  my($infile)= @_;
  my($nids,$ndups,$inh)=(0,undef);
  %dupids= %dupfirst=();
	# ($ok,$inh)= openRead($infile,1);
  open(IN,$infile) or die "reading $infile"; $inh=*IN; 
  while(<$inh>) {  
    s/^>//; # if from fastanrdb
    my  @dupids= split; # grep /\w/ or any other?
    next unless(@dupids>1); #? or record all ids?
    my $firstid= $dupids[0];
    if($dupids{$firstid}) {
      my $nextdup= $firstid;
      $firstid= $dupids{$nextdup};
      shift @dupids;
      push @{$dupfirst{$firstid}}, @dupids;  
    } else {
      $dupfirst{$firstid}= \@dupids;  $nids++; 
    }
    map{ $dupids{$_}= $firstid } @dupids;
    $ndups += @dupids;
  } close($inh);  
  warn "# readDupIds: nfirst=$nids; ndups=$ndups\n" if $debug;
  return ($nids,$ndups);
}


=item readEqualGene

Table of map equivalences, add as adjunct to align tab of blast/lastz, which have omissionn mistakes ~10%?
	$evigene/scripts/equalgene.pl -in kf2mixx.main.gff -over kf2mixx.main.gff > kf2mixx.main.eqgene
* need to change main/over ids to input oids somewhere ..
Mainid                  Oid                     Overlapids
Funhe2Exy3m110884t1     Fungr1EG3m041003t1      Funhe2Exy3m129932t1/C97.72,Funhe2Exy3m110083t1/82.77    10098sc:638046-638339:.
Funhe2Exy3m125954t1     Fungr1EG3m043564t1      na      479sc:197743-197931:.
Funhe2Exy3m091394t1     Fungr1EG3m037315t1      Funhe2Exy3m076422t1/I100,Funhe2Exy3m083869t1/89.90,Funhe2Exy3m044035t1_C1/87.88,Funhe2Exy3m063866t1/87.87,Funhe2Exy3m054191t1_C1/87.83,Funhe2Exy3m059312t1/87.81,Funhe2Exy3m025662t1/87.79,Funhe2Exy3m060328t1_C1/86.87,Funhe2Exy3m050942t1/86.86,Funhe2Exy3m052474t1/86.86,Funhe2Exy3m085750t1/86.86,Funhe2Exy3m050513t1/86.84,Funhe2Exy3m040128t1/86.82,Funhe2Exy3m052598t1/86.78,Funhe2Exy3m052599t1/86.78,Funhe2Exy3m049962t1/85.85,Funhe2Exy3m054654t1/85.85,Funhe2Exy3m123553t1/57.83,Funhe2Exy3m026282t1_C1/53.87,Funhe2Exy3m090395t1_C2/0.97    1967sc:10172-10679:-

FIXME: need to adjust eqgenes to account for poor gmapping, 
  paralogs map to same locus, poorly, 
  but blastn/any align says they are different loci
BUT also have perfect gmap eqgenes, not scored same by blastn, should be called same locus dups
	-- add mapqual to eqgene.table ? need it per geneID
	-- or separate table?

Fixme: dont output these to trclass
#m2t: ERR: trclass inline:Funhe2Exy3m149508t1_G2        drop    altmap  Funhe2Eq7m074376t1      99/35           0,0,pflag:3

grep  Funhe2Exy3m149508t1_G2 kfish2evg367mixx.eqgene
Funhe2Eq7m074376t1      Funhe2Eq7m074376t1      Funhe2Exy3m149508t1_G2/35.55    1095sc:48614-48822:.
Funhe2Exy3m149508t1_G2  noid    Funhe2Eq7m074376t1/35.55,Funhe2E6bm091492t1/35.39       1095sc:48708-48907:.
Funhe2E6bm091492t1      Funhe2E6bm091492t1      Funhe2Exy3m149508t1_G2/35.39    1095sc:48614-48785:.
	
=cut

=item eqgene classifier

  - do this in readEqualGene, mark which overlap alts are bad, which ok
  - also need to regard mapqual align, Split values to decide
  
  need separate classifier to handle various eqgene attributes, decide which tr/alts are bad/good
  eg eqgene classifier for this case: good g131, bad g453t5
  .. for this case need to know that g453t5 <mismap> g453t1,2,3,4 .. count overlap/alt/locus? drop outliers? use Split info?
ok Anofunz4gEVm000131t1	noid	Anofunz4gEVm000452t5/19	KB668936:299087-306372:-	99a,99i,7287l,3x	0
ok Anofunz4gEVm000131t2	noid	Anofunz4gEVm000452t5/33	KB668936:299737-303981:-	98a,99i,4269l,2x	0
ok Anofunz4gEVm000131t3	noid	Anofunz4gEVm000452t5/58	KB668936:299632-301584:-	70a,99i,1809l,2x	0
ok Anofunz4gEVm000452t1	noid	na	KB668900:3696-9724:-	99a,100i,4410l,8x	0
ok Anofunz4gEVm000452t2	noid	na	KB668900:3696-9724:-	99a,100i,4410l,8x	0
ok Anofunz4gEVm000452t3	noid	na	KB668900:3696-9724:-	96a,100i,4329l,8x	0
ok Anofunz4gEVm000452t4	noid	na	KB668900:3696-8875:-	100a,100i,3372l,9x	0
bad Anofunz4gEVm000452t5	noid	Anofunz4gEVm000131t1/56,Anofunz4gEVm000131t2/56,Anofunz4gEVm000131t3/43	KB668936:300530-301921:-	85a,100i,2469l,4x,Spl:29%,KB668900	0
ok Anofunz4gEVm000452t6	noid	na	KB668900:6592-8636:-	79a,100i,948l,3x	0

  g131,3/3 alts are over 1/6 g453 alts
  -- should this indicate g453t5 is not-much-over other g453t alts?
  
=cut

sub readEqualGene {
  my($infile)= @_;
  my($nids,$nov,$nxeq, $nne)=(0) x 9; 
  my (%eqgenes,%eqexons,%neqalts,%gmapqual); # extended eqgene.tab cols
		  #* change this eqgene/altmap classing,
		  # NO cds-align b/n td,qd means something, poor mapping paralogs end up here
		  # .. need (a) gmap qual score (align,ident,split), and
		  # ..      (b) paralog flag reverse of eqgene, for classed-alts that *dont* gmap same locus
	
		  # TEST1602: ($eqgenes,$nineqgene,$neqgene,$gmapqual,$neqalts)= readEqualGene($eqgene) extended table

  # open(IN,$infile) or die "reading $infile"; $inh=*IN; 
	my($ok,$inh)= openRead($infile,1);
  while(<$inh>) {  
  	next unless(/^\w/);
		my @v=split"\t";
  	my($mid,$oid,$overids,$loc,$mapqual,$altnotover)=@v;
		# extended overeqgene from altclassed-cds-genome.blastn calc: add alt_notover column for paralog classing
		
		next if($mid =~ /_G\d+$/); # dup map id syntax
    $mid =~ s/_C\d+$//; #what of split maps? change mapqual?
    my $pmid= $oids{$mid}||$mid;
    
    # new bleqgene mapqual: 85a,100i,2469l,4x,Spl:29%,KB668900 == %align,%ident,bases,exons,Spl:split
		$gmapqual{$pmid}=$gmapqual{$mid}="$mapqual\t$loc"; # or $loc\t$mapqual
		
    if($altnotover and $altnotover ne "na") {
      my @nealt= split",",$altnotover; $nne++;
      map{ my $pa= $oids{$_}||$_; 
        $neqalts{$pmid}{$pa}= $neqalts{$mid}{$_}=1; 
        } @nealt;
    }
    
		next if($overids eq "na" or not $overids);
    $nids++; my $jov=0; 
		my @ovd= grep{$_ ne "na"} split",",$overids;
		foreach my $ov (@ovd) {
			my($ovd,$cx)=split"/", $ov; 
			next if($ovd =~ /_G\d+$/); # dup map id syntax
		  $ovd =~ s/_C\d+$//;  ## what of split maps: id _C[12] ? skip or chop _C?
			$cx=~s/^[IC]//; #dont.care# $cx="$cx.$cx" unless($cx=~/\.\d/);

			# overeqcdsloc.pl test syntax,may change, adds exon-equal:  
			#  Anofunz4jEVm000002t47		Anofunz4jEVm000002t1/100xe100,..
			# problem cases, mis-align due to missing exon align of true alt, should aln<100 reduce xe value?
			# Anofunz4hEVm000070t2	Anofunz4hEVm000070t1/87xe100  	KB668689:297566-312085:-	87a<<,99i,10536l,8x
			# Anofunz4hEVm000070t1	Anofunz4hEVm000070t2/89xe88   	KB668689:297566-312085:-	100a,99i,10314l,8x
			
			my($xeq)= ($cx=~s/xe(\d+)//)?$1:0; 
			
			#old# $cx.=".100" if($cx =~ /^I100/); $cx=~s/^[IC]//; 
			my($ca,$xa)=split /[\.]/,$cx;  $ca||=0;  $xa||=0;


			if($ca >= $MIN_EQGENE) { 
  		  my $povd= $oids{$ovd}||$ovd;
			  ## ca is not reciprocal equal now: mid-larger over ovd-smaller, ca is large portion of ovd
			  ## mid-smaller over ovd-larger, ca is small portion of ovd
			  $jov++;  $eqgenes{$pmid}{$povd}= $eqgenes{$mid}{$ovd}= $ca; # pct of mid aligned to ovd
			  if($xeq) { $nxeq++; $eqexons{$pmid}{$povd}= $eqexons{$mid}{$ovd}= $xeq; } # new, test
			  }
		}
		$nov++ if($jov);
	}
	warn "# readEqualGene: nequal=$nov/$nids\n" if $debug;
	return (\%eqgenes,$nids,$nov,\%gmapqual,\%neqalts,$nxeq,\%eqexons);
	#oold#return (\%eqgenes,$nids,$nov);
}
			

=item readAAnametab / readAAblast

  maybe revise this to also use trasm.names table, computed for publicset naming,
  has essentially same info w/ best ref per tr.
  names format:
  TrID   Name   Align_score  RefID  RepID (uniprot)
  Funhe2Exx4m000455t12    CDD: Na_trans_assoc, Sodium ion transport-associated    100%,245/230,1997       CDD:219069      pfam06512
  Funhe2Exx4m000455t12    Sodium channel protein type 5 subunit alpha     100%,2057/2016,1997     RefID:UniRef50_Q14524   RepID:SCN5A_HUMAN
  Align_score = val%,nalign/nref,ntr
  
  make aablast.tab from query-ref.blastp
   a. (best)?
   cat evg2anofunz4g-bugsref.aa.blastp | \
     env blsum=1 aa="$query.aa.qual,$ref.aa.qual" ISCORE=4 nst=3 $evigene/scripts/blast2bestgenes.pl \
     > evg2anofunz4g-bugs6ref.aa.bltab
     # bltab: Query   Source     Bits    Iden    Algn    Qlen    Slen
   b. older, mostly same for cols 1-7
    $evigene/scripts/makeblastscore3.pl -tall -aasize "$query.aa.qual,$ref.aa.qual" evg2anofunz4g-bugsref.aa.blastp \
     > evg2anofunz4g-bugsref.aa.btall
   
     
=cut 

sub readAAnametab
{
  my($aanametab)= @_;
  my $swapids=0; my $naabl=0; my $nblerr=0;
  my %bscore;
  open(F, $aanametab) or die "FAIL: read $aanametab";
  while(<F>) { next unless(/^\w/); 
    my($td,$name,$alnscore,$rd,@more)=split"\t"; 
    # alnscore format expected: 72%,3270/4555,3282 ;  may be '72' or '72%' only
    #x my($apct,$aln,$refsize,$tsize)= $alnscore =~ m=^(\d+)%,(\d+)/(\d+),(\d+)=;
    my($apct,$aln)= $alnscore =~ m=^(\d+)%,(\d+)\b=; # partial match ok
    unless($aln){ ($aln)= $alnscore =~ m/^(\d+)/; }  # only care about relative score to other td
    next unless($aln and $aln>0);

    my $bscore= $aln;
    $naabl++; 
    $rd =~ s/^RefID://; # or CDD: or other
    unless($OIDFIX or $naabl>9) { 
      $OIDFIX=1 if(%oids and ($oids{$td} or $oids{$rd}));
    }
    if($OIDFIX){ $td= $oids{$td}||$td; $rd= $oids{$rd}||$rd; }
    # local $^W = 0;   # no warnings; no help
    unless($aablast{$td} and $bscore{$td} > $bscore) {
      $aablast{$td}="$bscore,$rd";  $bscore{$td}= $bscore;
      unless( $aablastref{$rd} and $bscore{$rd} > $bscore ) {
        $aablastref{$rd}= "$bscore,$td";  $bscore{$rd}= $bscore;
        $aablastref{$td}= "$bscore,$rd";   # revhash?
      }  
    } else {
      unless( $aablastref{$rd} and $bscore{$rd} > $bscore ) {
        $aablastref{$rd}= "$bscore,$td";  $bscore{$rd}= $bscore;
        $aablastref{$td}= "$bscore,$rd";   # revhash?
      }  
    }
    
  } close(F); 
 warn "# readAAnametab: naabl=$naabl\n" if $debug;
 return($naabl); 
}


sub readAAblast 
{  
  my($aablast)= @_;
  my $swapids=0; my $naabl=0; my $nblerr=0;
  my %bscore;
  if($aablast =~ /\.names$/) { return readAAnametab($aablast); }
  
  ## precheck format before sort of possibly very large file ..
  use constant nCHECK => 29;
	my($ok,$inh)= openRead($aablast,1);
  # open(F,$aablast) or  die "FAIL: read $aablast ..."; 
  while(<$inh>) { next unless(/^\w/); 
    ## FIXME: allow orig blastp table format, not postprocess .tall4 version ?
    ## at least check for blastp table
    my($td,$rd,@v)=split"\t"; $naabl++; # expect: trid refid bits ident align .. may be refid,trid instead

    if(%oids) { 
      $OIDFIX=1 if($oids{$td} or $oids{$rd});   # is aasize{oid} ok here
      if($OIDFIX){ $td= $oids{$td}||$td; $rd= $oids{$rd}||$rd; }
    }
    if($swapids==1) {
      ($td,$rd)=($rd,$td);
    } elsif($swapids==0) {
      if(@v >= 12 and $v[8] =~ /^\d/) { $nblerr++; } # blast.tab
      else {
      if($aasize{$td}) { $swapids= -1; }
      elsif($aasize{$rd}) { $swapids= 1; }
      else { $nblerr++; }
      }
    }
    
    last if($naabl > nCHECK);
  } close($inh);
  if($nblerr>2 or $swapids==0) { 
    die "ERR: expect table of aablast scores: trid refid bitscore identity align ..\n"
      ." $nblerr trids from aasize not found in -aablast=$aablast\n";
  }
  
  $naabl=$nblerr= 0;  
  #?? add -aablastsorted option? sort wastes time if not needed.
  #old# open(F,"sort -k3,3nr $aablast |") or die "FAIL: read $aablast ..."; 
	# ($ok,$inh)= openRead($aablast,1);
	
	## FIXME sort for swapids, td first, refd 2nd when tied scores, prob not much effect
	## FIXME TEST1602: -k3 bitscore sort may be flaky, opt try -k4 ident or -k5 align
	my $sortord=($swapids==1)?'-k3,3nr -k5,5nr -k2,2 -k1,1' : '-k3,3nr -k5,5nr -k1,1 -k2,2';
  # open(F,"sort $sortord $aablast |") or die "FAIL: read $aablast ..."; 
  $ok= ($aablast =~ /\.gz$/) ? open($inh,"gunzip -c $aablast | sort $sortord |") 
        : open($inh,"sort $sortord $aablast |");  
  
  while(<$inh>) { next unless(/^\w/); 
    ## FIXME: allow orig blastp table format, not postprocess .tall4 version ?
    ## at least check for blastp table
    my($td,$rd,@v)=split"\t"; $naabl++; # expect: trid refid bits ident align .. may be refid,trid instead
    
    if($OIDFIX) { # and %oids %oid2pubid; aablast, others?
      $td= $oids{$td}||$td; $rd= $oids{$rd}||$rd;
    }
    
    if($swapids==1) {
      ($td,$rd)=($rd,$td);
    }
    
    my $bscore= $v[0]; # bits, ident, algn; want choice of? ident or algn maybe better
    # local $^W = 0;   # no warnings; no help
    ## maybe fixme: missing some uniq blastref{rd} here? pull aablastref out of aablast{} ?
    unless($aablast{$td} and $bscore{$td} >= $bscore) {
      $aablast{$td}="$bscore,$rd";  $bscore{$td}= $bscore;
      unless( $aablastref{$rd} and $bscore{$rd} >= $bscore ) {
        $aablastref{$rd}= "$bscore,$td";  $bscore{$rd}= $bscore;
        $aablastref{$td}= "$bscore,$rd";   # revhash?
      }  
    } else {
      unless( $aablastref{$rd} and $bscore{$rd} >= $bscore ) {
        $aablastref{$rd}= "$bscore,$td";  $bscore{$rd}= $bscore;
        $aablastref{$td}= "$bscore,$rd";   # revhash?
      }  
    
    }
    
  } close($inh); 
 warn "# readAAblast: naabl=$naabl\n" if $debug;
 return($naabl); 
}


use constant SPANSUM => 1;
use constant ADDSPANS1605 => 1;

sub putspans {
  my($lq)= @_;
  # fixme: save no-match lq ids:  lq, aq, wq, na.... or self-score ?
  my $nmatch=0;
  # fixme: sort output bspans by tidn/taln? $bspans{$b}->[0]->[4] = xbit
  # fixme? throw away dupids here? should be in same bspans align cluster
  
  # our(%dupids, %dupfirst);
  my $dupfirst=""; 
if(DUPFILTER2) {    
  if($dupids) { $dupfirst= $dupids{$lq} || ""; }
}
  
  foreach my $lt (sort keys %bspans) {
    next if($lt eq $lq); # is lt eq lq allowed here?
    
    ## move this dupid filter before into readblastab, readcdhit ? see active DUPFILTER1 
if(DUPFILTER2) {    
    if( $dupids and $dupids{$lt} ) {
      if($dupfirst eq $dupids{$lt}) { $ndupdrop++; next; }
      $dupfirst= $dupids{$lt};
    }
}
    
    my @bspans= @{$bspans{$lt}};
    my ($tbit,$taln,$tidn,$aq,$at,$wq,$wt,$torient,$xbm,$xem,$tbm,$tem)= (0) x 19;
    foreach my $sp (@bspans) {
      my($xb,$xe,$tb,$te,$xbit,$aln,$aident,$or)= @$sp; # 2013.aug: IS_CDSALIGN add $or
      $tbit += $xbit; $taln+= $aln; $tidn+= $aident; 
      ##$torient+=$or; ## weight by aln so tiny -or dont throw it off?
      $torient += $aln * $or; ## weight by aln so tiny -or dont throw it off?
if(ADDSPANS1605) {
      $xbm=$xb if($xbm==0 or $xb<$xbm); $xem=$xe if($xe>$xem);
      $tbm=$tb if($tbm==0 or $tb<$tbm); $tem=$te if($te>$tem);
}
      }
    # my $mis= $taln - $tidn; # dont need this; replace w/ tidn
    if($OUTSPANTAB) { 
      # add Qsize,Tsize  cds/trlen ?
      $aq= $aasize{$lq}||0; $aq *=3;  # 3*convert to cds-size
      $at= $aasize{$lt}||0; $at *=3;
      $wq= $trsize{$lq}||0; 
      $wt= $trsize{$lt}||0;  
      ## NO: add here?? $aaqual{$lq} ; $aaqual{$lt}; 
      ## FIXME taln may be rel CDS-len (aq,at) or TR-len (wq,wt)
      my $alnmax= _max(1, ($IS_CDSALIGN) ? _min($aq,$at) : _min($wq,$wt) );
      
      ## UPD 2016.02: -tinyaln 35  bad (sometimes), need opts and diff test: tiny by bases not pct overlap
      ## eg. TINYALNBASES= 90; MIN_PCTOVERLAP=10;
  if(UseTINYALNBASES) {
      if($alnmax >= $MINCDS) { # never skip tiny prots, MINCDS =~ 90 bases ?
        my $minbaseover= ($alnmax >= $TINYALNCDSLEN)? $TINYALNBASES : $alnmax * 0.10; ##$MIN_PCTOVERLAP;
        next if($taln < $minbaseover);
      } 
  } else { # old, pct-TINYALN problem    
      next if( (100 * $taln / $alnmax ) < $TINYALN); # skip trival matches
  }      
      $nmatch++;
      
      # 2013.aug: IS_CDSALIGN add $torient here? only care if $tor < 0; add sign to one of these cols?
      ## cant use taln, will screw up sort -n by maxalign; use tbit, unused now in scoring
      if($IS_CDSALIGN and $torient<0) { $tbit= -$tbit; } #? what followon problems does this cause?
if(ADDSPANS1605) {
      print $OUTH join("\t",qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits Qspan Tspan))."\n" unless($head++);
      print $OUTH join("\t",$lq,$aq,$wq,$lt,$at,$wt,$taln,$tidn,$tbit, "$xbm-$xem","$tbm-$tem")."\n"; 
} else {      
      print $OUTH join("\t",qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits))."\n" unless($head++);
      print $OUTH join("\t",$lq,$aq,$wq,$lt,$at,$wt,$taln,$tidn,$tbit)."\n"; 
}
      $validids{$lq}++; $validids{$lt}++;
      } 
    else { 
      puts($lq,$lt,$taln,$tidn);   $nmatch++; # DROP this old version
    }
  } 

if(DUPFILTER2) {  
  if($nmatch==0 and $dupfirst and $dupids{$lq} and $dupids{$lq} ne $dupfirst) { $nmatch=-1; $ndupdrop++; }
}
  if($nmatch==0) {
    if($OUTSPANTAB) { 
      my ($tbit,$taln,$tidn,$aq,$at,$wq,$wt)= (0) x 19; 
      my $lt="self"; # or lq, or use blast-self scores?
      $aq= $aasize{$lq}||0; $aq *=3;  # 3*convert to cds-size
      $wq= $trsize{$lq}||0; 
      $at=$aq; $wt=$taln=$tidn=$tbit= $wq; # or aq
if(ADDSPANS1605) {
      print $OUTH join("\t",qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits Qspan Tspan))."\n" unless($head++);
      print $OUTH join("\t",$lq,$aq,$wq,$lt,$at,$wt,$taln,$tidn,$tbit, "1-$wq","1-$wt")."\n"; 
} else {      
      print $OUTH join("\t",qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits))."\n" unless($head++);
      print $OUTH join("\t",$lq,$aq,$wq,$lt,$at,$wt,$taln,$tidn,$tbit)."\n"; 
}
      $validids{$lq}++; $validids{$lt}++;
    } else {
    
    }
  }
  
  %bspans=();
  return($nmatch); # ,$ndupdrop
}



=item identity class constants

=cut

use constant { kAATINY => 1, kAAUTRBAD => 2, kAAUTRPOOR => 4, kAADUP => 8, kAAGAPS => 16, 
              kNONCODE => 32, }; # kNONCODE == cds.qual Noncode flag
              # kDUPEXONS => 64 ?? for cullExonEq
## eor flags to clear others
use constant NOTTINY => kAADUP + kAAGAPS + kNONCODE; # - kAATINY - kAAUTRPOOR - kAAUTRBAD - kAAGAPS
use constant NOTPOORBAD => NOTTINY + kAATINY; # kAATINY + kAADUP + kAAGAPS + kNONCODE; # - kAAUTRPOOR - kAAUTRBAD - kAAGAPS
use constant NOTPOOR => NOTPOORBAD + kAAUTRBAD; # kAATINY + kAADUP + kAAUTRBAD + kAAGAPS + kNONCODE; # - kAAUTRPOOR

=item identity align defaults

$AAMIN =$ENV{aamin}||30; #was 40;   # for aacomplete, utrok

$TINYALN = 20; # %align, was 25; was $ENV{mina}||50; ignore less than this == MINAL
$TINYALNBASES= $ENV{MINALIGN}||90; # was $MINALIGN= 90; # REUSE//NOT USED NOW; use just TINYALN; change to several levels
$TINYALNCDSLEN= 3*$TINYALNBASES;
$MINCDS = $ENV{MINCDS} || 90; # or 3*MINAA  

$ALTFRAG= 0.5;

$PHIALN=$ENV{ahi} || 98; # was 99; was 65 !!;  DROP THIS; use mina?
$NHIALN=$ENV{nhi} || 9; # what? for short cds cancel few base diff < PHIALN

$PHI = $ENV{phi} ||99; 
$PMID= $ENV{pmid}||90; 
$PLOW= $ENV{plow}||80;  

$OK_CDSUTR_default= 60;
$BAD_CDSUTR_default= 30; # % CDS/trlen

=cut

=item  overlocusclass, overlocus_eqgene

  temp in asmrna_dupfilter3_overlocsub.pm
  now in asmrna_overlocusclass4a.pl

=cut

#this works WITH our(%vars) # 
sub overlocusclass4a_INCLUDE { }
#o BEGIN{ require "asmrna_overlocusclass4a.pl"; } # test replacing _PRE4a subs below.


#========================= identityclass UPD4a ============================

#UPD1912 this works WITH our(%vars) # TEST code in idclass4.pl, MOVE into this code for pub uses
sub identityclass4a_INCLUDE { }
#======== 20.01.14 moved asmrna_identityclass4a.pl back to  asmrna_dupfilter4a.pl ===========================
# BEGIN{ require "asmrna_identityclass4a.pl"; } # test replacing _PRE4a subs below.

=item identityclass

  from  aabugs4qual/tsaevg/cdsidentclass.sh
  input from above putspans, sorted by ^Qclen, ^Align, vTtlen:
     qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits) 

  update classes per aabugs4qual/tsaevg/classparts.pl
  .. need blastp/blast table input for this..
  .. need 3-4 final categories:  keep-main, keep-alts, trash-fragments (alt,noalt) 
     .. keep includes partial prot where good enough.
     .. trash-alts maybe (not fragment but not distinct enough)
     .. add poor cut to main,alt,noclass per below a,b qualities.
     
  cl_main.alttr and cl_alts, keep if: 
  a. homology unique or best (near best for alt),  or
  b. aasize >= 45% of main and aacomplete? and not aasame-as-main?
  c. main == fragment if aasize < MINAAFULL or (aasize < MINAAPART and partial/utrpoor)

  identityclass parameters ... see above now
    $TINYALN=$ENV{mina}||50; 
    $IS_CDSALIGN=$ENV{cdsw}||0; # default using tralign, tests better than cds
    $PHIALN=$ENV{ahi}||98; # was 65; # use mina?
    $PHI=$ENV{phi}||99; $PMID=$ENV{pmid}||90; $PLOW=$ENV{plow}||80; 
    $ALTFRAG: add isfrag pct; now 50 (0.5)
  
=cut

sub identityclass {
  my($outh, $infile, $insorted)= @_;

  use constant TEST3 => 1; # 13aug09 test fixes to alt/main classing # UPD1912: all accepted now
  use constant TEST1602 => 1; ## UPD 2016.02 problem w/ dup equal mains, equivalence after 1st see both as mains, 1 to be alt,
  use constant TEST1603 => 0; #UPD1912 off, only for now-removed eqgene measures, was 1; ## test use here eqgenes/eqexons??
  
  # our(%validids,$IS_CDSALIGN);
	my $havevalidids= scalar(%validids)?1:0; ## need this?

  ### infile is putspans() alntab: 
  ###  Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits
  ###    1     2     3   4     5     6     7     8    9
  ### sort ^Qclen, ^Align, vTtlen  : orig
  ### sort ^Qclen, ^Align, vQtlen, vTtlen, =Qid, =Tid : now
  ### maybe ^Ident before =Tid ?  after ^Align ?
  ### should it be ^Qclen, =Qid, ^Align, vTtlen :  no want top Aligns first
  ### added vQtlen before vTtlen, else get utrbad before good ***
  # my $ALNSORTORD='-k2,2nr -k7,7nr -k6,6n '; #orig
  # my $ALNSORTORD='-k2,2nr -k7,7nr -k6,6n -k1,1 -k4,4'; #add IDs to order ties

  my $ALNSORTORD='-k2,2nr -k7,7nr -k3,3n -k6,6n -k1,1 -k4,4';  #add vQtlen/3 before vTt/6 IDs to order ties
  # unless($IS_CDSALIGN) { #..not sure what is right, want cds to play role in best choice
  #   $ALNSORTORD='-k2,2nr -k7,7nr -k3,3n -k6,6n -k1,1 -k4,4';  # test same as IS_CDSALIGN  
  # }
  
  my($inh, %class,%bestmatch,%ismain, %mainof, %havepair); 
  if(ref($infile)) { $inh= $infile; } # STDIN,.. sort ??
  else {
    my $tmpdir=$ENV{TMPDIR};
    $tmpdir= './' unless($tmpdir and $tmpdir ne "/tmp"); ## { $tmpdir=`pwd`; chomp($tmpdir); }
    if($insorted) { open(IN,$infile) or die "reading $infile"; $inh=*IN; }
    else { open(IN,"sort -T $tmpdir $ALNSORTORD $infile |") or die "reading $infile"; $inh=*IN; }
  }
  
  use constant UPD1912doswapc => 0; # 0, turn this off, likely bug
  use constant UPD1912fixNOMAIN => 1; # use DEBUGnomain == $ENV{fixnomain} for testing
  use constant ANTI1912c => 1; ## antisense checks
  my $PAL_ANTIMIN= 60; #ANTI1912c opt? needs test w/ refs
  # ref_arath: antisense pAlign for same gene, diff gene == 60 is good
  #  sameg.pal	n:146 av:75 med:65 81 90;  quartiles: .25, .50, .75
  #  diffg.pal	n:333 av:66 med:44 70 91
  
  my $DEBUGnomain= $ENV{fixnomain}||0; # or UPD1912fixNOMAIN
  warn "# UPD1912doswapc => ",UPD1912doswapc," UPD1912fixNOMAIN => ",UPD1912fixNOMAIN," ANTI1912c => ",ANTI1912c," fixnomain=",$DEBUGnomain,"\n";
    
  my($lastd)=("");
  while(<$inh>) {  # input sorted align table
    next if(/^Qid|^\W/); chomp; 
    my @v= split"\t"; 
    my($qd,$qc,$qw,$td,$tc,$tw,$aln,$iden,$bits)= @v; 
    my($isfrag,$aclass,$alnmax,$pid,$pal,$samesize)= (0) x 10;
    
    my $isanti= ($IS_CDSALIGN and $bits < 0)?1:0;
    my $antiflag= ($isanti) ? "/-sense" : "/."; ## append to bestmatch ??
    $isfrag= $aclass="";
    $samesize=($tc == $qc)?1:0; # OR ? abs($tc - $qc) < $NHIALN
    my $swapqt=($tc > $qc)?1:0;
    
    if(AACONS and UPD1912doswapc) {
      if($nconsensus and $samesize) {
        #UPD1912 : Is this swapqt a problem, for dup loci, almost identicals? changes their sort order some what randomly
        # .. Could replace this aacons boost by adding .acons as decimal to alntab Qclen, Tclen so sort will favor tc.tcon > qc.nocon
        my $tdcon= $aconsensus{$td}||0; my $qdcon=$aconsensus{$qd}||0;
        $swapqt=1 if($tdcon > $qdcon);
      }
    }
    if($swapqt){ ($qd,$td)=($td,$qd); ($qc,$qw,$tc,$tw)=($tc,$tw,$qc,$qw); } # swap to greater cds?

    # eqgene classifier : drop for now
		# if(EQGENE_CHANGES_NOALN && $lastd && $lastd ne $qd) {}
				
		$lastd=$qd; # here?
    if($td eq "self" or $qd eq $td) { # putspan fix for no-align cases
      $td=$qd; # $aclass="noclass"; << should be this, but ??
      $bestmatch{$qd}="$td,100/100/self1" unless($bestmatch{$qd});
      next;  # can lose eqgene if no further td/qd ..
    }
    
    my($qsize,$tsize)= ($IS_CDSALIGN) ? ($qc,$tc) : ($qw,$tw);
    # note: alnmax is min(qsize,tsize) not max
    $alnmax= ($qsize>$tsize and $tsize>0)?$tsize:($qsize>0)?$qsize:$tsize;  
    $isfrag= ($tsize < $ALTFRAG*$qsize)?"frag":"";
    
    $pid= ($aln<1)?0: int(0.5+ 100*$iden/$aln); $pid=100 if($pid>100);
    $pal= ($alnmax<1)?0 : int(0.5+ 100*$aln/$alnmax); $pal=100 if($pal>100);
    
    my $tinyalndiff= ((($alnmax - $aln) < $NHIALN)) ? 1 : 0; # TEST3 add, use w/ $pal

    ##UPD1912 drop eqgenes parts for now
		#u $havepair{$qd}{$td}= $pal; #UPD1912:not used ?? part of EQGENE; change this to pal align val for altpar TEST1603?
		#u if($neqgene>0) {}
		
    if(ANTI1912c and $isanti) { # UPD1912 maybe .. TEST .. NOT HERE, need bestmatch for output trclass
      $antiflag= "/." if($pal < $PAL_ANTIMIN); # turn off flag
      #   next if($pal < 60); # force separate loci for rev-aa, partial cds align??
      #    # treat as NO ALIGN ? skip before bestmatch NO cant skip bestmatch
    } 

    my $pidalnval="$pid/$pal$antiflag"; # old: $pid/$pal
    my $tisaltpar=0; # treat like $skiptinyaln for now
    if(TEST1603) {
      $tisaltpar=1 if($pidalnval=~/altpar\d/); #?  and $qclass // not altparx\d
    }
    if(TEST1603) { # UPD1912:off, this may be sort of working right now..
    $bestmatch{$qd}="$td,$pidalnval" unless($bestmatch{$qd} and not($bestmatch{$qd}=~/altpar\d/));
    $bestmatch{$td}="$qd,$pidalnval" unless($bestmatch{$td} and not($bestmatch{$td}=~/altpar\d/)); 
    } else {    
    $bestmatch{$qd}="$td,$pidalnval" unless($bestmatch{$qd});
    $bestmatch{$td}="$qd,$pidalnval" unless($bestmatch{$td});  
    }
        
    my $qclass= $class{$qd}||""; my $tclassD= $class{$td}||"";
    #D warn "#dnomain1 $qd,$qc $td,$tc pid=$pid, pal=$pal, al=$aln, qcla=$qclass, tcla=$tclassD\n" if($DEBUGnomain);

    if($samesize and $qclass =~ /alt/ and $ismain{$td}) { # check/stop alt1 <> alt2 assigns
      next if($bestmatch{$qd} =~ /$td,/); # problem here? for many equal bestmatch
    }

    #..........    
    ## reclassify? althi > althi100 for 2ndary alts; althi100 = unneeded duplicates : TEST
    
    my $skiptinyaln=0;
    if(UseTINYALNBASES) {
      if($alnmax >= $MINCDS) {  
        my $minbaseover= ($alnmax >= $TINYALNCDSLEN)? $TINYALNBASES : $alnmax * 0.10; ##$MIN_PCTOVERLAP;
        $skiptinyaln=($aln < $minbaseover)?1:0;
      } else {
        $isfrag= "frag" unless($isfrag); #or "tiny" ? # dont skip assign alt/frag to tiny cds, but maybe always set isfrag ?
      }
    } else {  
      $skiptinyaln= ($pal < $TINYALN)?1:0;
    }
  
    # set alt aclass here, or not.
    if( $skiptinyaln or $tisaltpar) { } ## defer: $aclass="noalign$isfrag";  
    elsif( $tc == 0 and $qc > 0 and $pid >= $PLOW ) { $aclass="frag0aa"; } # tc==0 special case "frag0"
    elsif( $pid >= $PHI ) { 
      $aclass= ($isfrag)?"parthi":"althi"; # Primary classifier; change parthi to althipart ? althi$isfrag ?
      $aclass .= "1" if($pal >= $PHIALN or $tinyalndiff); # althi1 == hi-align + hi-ident, nearly same (identicals gone here)
      }
    elsif( $pid >= $PMID ) { $aclass="altmid$isfrag"; }
    elsif( $pid >= $PLOW ) { $aclass="altlo$isfrag"; }
    # else { $aclass=""; } # otherwise this pair are not aligned well enough to call same locus

    if(ANTI1912c and $isanti) { # UPD1912 maybe .. TEST .. NOT WORKING ??
      $aclass="" if($pal < $PAL_ANTIMIN); # force separate loci for rev-aa, partial cds align??
      $pidalnval =~ s/.sense/./ unless($aclass);
      ## $antiflag= ($isanti) ? "/-sense" : "/."; 
    } 
    
=item samesize/samebestmatch 
  problem of 2+mains w/ several essentially identicals, 
  after 2 mains created, then row linking them needs to break tie
  
  * problems losing main links with large num 910 $samesize pairs .. sort order isn't enough
  
=cut

    if(my $ocl= $class{$td}) {  
      if($ocl =~ m/^alt|^part/) { 
        $class{$td}= $aclass if($aclass eq "althi1"); 
        next;   # TEST1602  ok 
        # UPD1912, check this, NOMAIN major group now, bad
        # maybe want to drop into FINDMAIN below, add qd main if missing, but dont change old td class to aclass?
        #no if($DEBUGnomain == 2) {  ## UPD1912fixNOMAIN .. TEST 
        #no  $aclass=$ocl unless($aclass eq "althi1"); # and continue into FINDMAIN
        #no } else { next;  }  # 
      } elsif($ocl =~ m/^main/) { 
        if($qclass =~ m/^main/ or $tisaltpar) { 
          # skip to below aclass set, test5, yes, do this also?
          } 
        else { next; } # * YES, need this : test4
      } 
    }
      
    if($aclass) { 
      my $qmain=$qd; my $attr= $pidalnval; ## "$pid/$pal$antiflag";
      ## use constant FINDMAINID => 1;  
      if($class{$qd}) {   
        my $qnext= $qd; 
        my $more= ($bestmatch{$qnext})?1:0;
        
        my %qdid=( $qd => 1, $td => 1);
        while($more) { 
          $more=0;
          my($qbest,$qpal)=split",",$bestmatch{$qnext};
          
          if($qbest) {
          $qbest=0 if($ismain{$qnext} and not $ismain{$qbest}); # TEST1603 want always? maybe bad?    
                
          #n if($DEBUGnomain) { 
          #n# this way adds more /qmain tags, doesnt resolve alt1<=>alt2 lost mains, for samesize near identicals
          #n$qbest=0 if($class{$qnext} eq "main" and  $class{$qbest} ne "main"); # TEST 1912
          #n} else {
          #n $qbest=0 if($ismain{$qnext} and not $ismain{$qbest}); # TEST1603 want always? maybe bad?          
          #n}
          }
          
          if($qbest and not $qdid{$qbest}) {
            $qnext= $qbest; $qdid{$qbest}=1; 
            $more=1 if($class{$qnext} and $bestmatch{$qnext});
          }
        }
        if($qnext ne $qd) { $qmain= $qnext; 
          unless(UPD1912fixNOMAIN) { $attr="$pidalnval/$qmain"; } # $DEBUGnomain
          } 
      }
      
      $ismain{$qmain}++;   
      $class{$td}=$aclass; #? problem here when class{td} == main set before, should not reset to alt
      
      # TEST1602 this way works.  if($bestmatch{$td}) must be true from above
      if(not $bestmatch{$td} or $bestmatch{$td} =~ /^$qd/){ $bestmatch{$td}="$qd,$attr"; }
      #old if($bestmatch{$td} =~ /^$qd/) { $bestmatch{$td}="$qd,$attr"; }
      #old elsif(not $bestmatch{$td}) { $bestmatch{$td}="$qd,$attr"; }

      ## add to prevent 2+ circular alta <> altb <> altc with no main ..
      $class{$qmain}="main" unless($class{$qmain}); # should always be:  $qc >= $tc//$samesize..
      if(UPD1912fixNOMAIN) { $mainof{$td}= $qmain unless($mainof{$td}); }
      #D warn "#dnomain2 qmain:$qmain $td,$tc  qmcla=$class{$qmain}, tcla=$class{$td}\n" if($DEBUGnomain);
    }
      
  } close($inh);

  # UPD1912fixNOMAIN/DEBUGnomain .. findmain() may work, does on simple test case w/ many samesize dups
  # should it be used for other than %class keys? probably not, mainof only set for %class'ed
  sub findmain { 
    my($td,$mainof)=@_;
    my $nextd=$td; my(%mdid);
    while (my $md= $mainof->{$nextd}) { last if($mdid{$md}++); $nextd=$md; }
    return ($nextd eq $td) ? "" : $nextd;
  }
  
  # END:  print trclass table; add more fields to output: aaqual, aablast, tr,aa sizes?
  { my($q,$pal,$c,$d);
  foreach $d (sort keys %class) {
    ($q,$pal)=split",",$bestmatch{$d}; $c= $class{$d}; # alts mostly here in class list
    if(UPD1912fixNOMAIN) { my($md)= findmain($d, \%mainof); $pal .= "/$md" if($md and $md ne $q); }
    my @cla= classifytr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  foreach $d (sort grep{ not($class{$_}) } keys %ismain) { # mains, not in class list
    ($q,$pal)=split",",$bestmatch{$d}; $c= "main";    
    #d if(UPD1912fixNOMAIN) { my($md)= findmain($d, \%mainof); $pal .= "/$md" if($md and $md ne $q); }
    my @cla= classifytr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  foreach $d (sort grep{ not($class{$_} or $ismain{$_}) } keys %bestmatch) { # unclassified by shared aligns
    ($q,$pal)=split",",$bestmatch{$d}; $c= "noclass";  
    #d if(UPD1912fixNOMAIN) { my($md)= findmain($d, \%mainof); $pal .= "/$md" if($md and $md ne $q); }
    my @cla= classifytr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  }  
  
}


=item classifytr

  classifier of locus primary, alternate, fragment, redundant 
  using identityclass() collection of overlapping transcripts (CDS or full tr)  
  focused on CDS qualities, maybe should be classifyCDS() ..
  
=cut


sub classifytr {
  my($tid,$cla,$qid,$pidal)= @_;
  $qid||="0"; $pidal||="0";

  # use constant { kAATINY => 1, kAAUTRBAD => 2, kAAUTRPOOR => 4, kAADUP => 8, kAAGAPS => 16, kNONCODE => 32, };
  # use constant NOTPOORBAD => kAATINY + kAADUP + kAAGAPS; # - kAAUTRPOOR - kAAUTRBAD - kAAGAPS
  # use constant NOTPOOR => kAATINY + kAADUP + kAAUTRBAD + kAAGAPS; # - kAAUTRPOOR
  
  my($pidi,$pali)= $pidal =~ m/^(\d+).(\d+)/;
  my($isanti)= $pidal =~ m/\-sense/ ? 1 : 0;
  my $aw= $aasize{$tid} || 0;
  my $tqual= $aaqual{$tid} || ""; #? parse for "aasize,pcds,aaqual"
  my $tw= $trsize{$tid} || 1; 
  my $pcds= int(300*$aw/$tw);  
  my $butr= ($tw - 3*$aw); # okutr if butr<=300p, modify BAD_CDSUTR for small aw
  my $ispoor= ($aw < $AAMIN or ($tqual =~ m/partial/ and $aw < $AAPART))?kAATINY:0;
  my $tbits= $aablast{$tid} || "0,0";
  my($tbscore,$tbref)=split",",$tbits;
  my $aacons= (AACONS) ? ($aconsensus{$tid}||0) : 0;
  my $mapqloc= ($nineqgene) ? $gmapqual->{$tid} : 0;  
  my($mapqual,$maploc)= ($mapqloc)?(split"\t",$mapqloc):(0,0);
  # new mapqual: 85a,100i,2469l,4x,Spl:29%,KB668900 == %align,%ident,bases,exons,Spl:split

  my($tbrscore,$tbrefbest)= ($tbref) ? (split",", $aablastref{$tbref}) : (0,0); 
  my $rbits2= $aablastref{$tid} || "0,0";
  my($tbrscore2,$tbref2)= split",", $rbits2; # this is revhash tid => score,ref refbest

  my $aaclus= $aacluster{$tid} || "0,0";  
  my($aamainid,$aamainpi)= split",", $aaclus; # ="$mainid,$pi";
  if($aamainid eq $tid) {  $aamainid=""; $aamainpi=0; }
  ## FIXME: aamain can be bad/missing; fix before this, need cds/tr align info to know if aamainid is good.
  
  if($tqual =~ m/gapbad/) { $ispoor |= kAAGAPS; }  # gaps discrepancy w/ aaqual utrbad and pcds utrbad 
  if($butr >= $MINUTR) {  
    if($tqual =~ m/utrbad/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRBAD; } 
    elsif($tqual =~ m/utrpoor/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRPOOR; }
    else { $ispoor |= ($pcds < $BAD_CDSUTR)? kAAUTRBAD: ($pcds < $OK_CDSUTR)? kAAUTRPOOR : 0; }
  }

  my $tcode=0;
  if(CODEPOT1607) {
    $tcode= $sizeval{$tid}{'codepot'}||"";  
    if($tcode and $tcode =~ /^Noncode/) {  
    $tcode=0 if($tbscore>0);
    $tcode=0 if($aw > 139 or ($aw > 99 and not $ispoor)); #?? which? many utrorf/utrbad for aa>99
    $ispoor |= kNONCODE if($tcode and $tcode =~ /^Noncode/);  
    }
  }

  if($aamainpi > $AADUP_IDENT) { # TEST aacluster added info on dup prots; AADUP_IDENT=98 default
    $cla.="a2" unless($cla=~/hi[012]/);  
    $ispoor |= kAADUP if($cla =~ /althi/);  
  }
  unless( $ispoor & kAADUP ) {
    if($pidal =~ m/altmap\d+xeq/ and $cla =~ /althi1|part/) {
      $ispoor |= kAADUP;   
    }
  }
  my $aadupflag= ($ispoor & kAADUP and $aamainid) ? ",aadup:$aamainid" : "";  

  # CHECKME: adding aablast kept 40k more okay, all althi/ahia2 + 2k parthi
  # .. are these true uniq aablast or just althi aadups ?
  my $MINBLASTSCORE= 60; # aablast only? bitscore always?
  my $keepdrop="";
  unless( $tbscore == 0 or $tbits=~/^0,0/) {
    my $risbest= (($tbrefbest eq $tid) or ($tbscore >= $tbrscore))?1:0;# was ($tbrefbest eq $tid), allow 2+ same score
    my $risgood= ($tbrscore > $MINBLASTSCORE and $tbscore >= 0.90 * $tbrscore)?1:0;
    my $risok= ($tbref2 and $tbrscore2 >= $MINBLASTSCORE)?1:0;
    if($risok and $risgood) { $risgood=0 if($tbrscore2 > $tbrscore); } #?
    $ispoor = $ispoor & NOTTINY if($risgood); ## FIXME1712: clear kAATINY for $risok/risgood
    my $isaadup=($ispoor & kAADUP);
    if($risbest) { 
      $tbits.=",refbest"; $keepdrop.="okay:refbest,"; $ispoor=0 unless($isaadup); }
    elsif($risgood) { # was 3rd
      $tbits.=",refgood"; unless($isaadup){ $keepdrop.="okay:refgood,"; $ispoor=0;} } # maybe2 ok ?
    elsif($risok) {  # was 2nd ? why superceed refgood? 2nd ref match?
      $tbits.=",refok"; unless($isaadup){ $keepdrop.="okay:refok,"; $ispoor=0;} }
  }
  
  if($ispoor > kAATINY and $cla !~ /althi|parthi/) {  
    if(UPD1908) { } # maybe rescue althi also here
    if($aw >= $AAMINBAD) { $ispoor = $ispoor & NOTPOORBAD; } ##  ^ (kAAGAPS+kAAUTRPOOR+kAAUTRBAD);  ##?? ^ XOR bad
    elsif($aw >= $AAMINPOO) { $ispoor = $ispoor & NOTPOOR; } ##  ^ (0+kAAUTRPOOR);
  }

  # cdsantisense, try to guess true cases from low align, not NONCODE, aacons, other ?
  if(UPD1912 and $isanti and ($ispoor or $cla =~ /parthi/)) {
    my $antiok= 0;
    $antiok += 1 if($pali < 75);  $antiok += 2 if($pali < 55);
    $antiok += $aacons if($aacons > 1);
    if($antiok) { 
      $ispoor = $ispoor & (kAATINY + kAADUP + kAAGAPS + kNONCODE);
      if($cla =~ /parthi/ and $antiok > 2) { $cla="altmidfrag"; } # change to not part? "altmid$isfrag";  
    }
  }
  
  if(AACONS) {
    if($aacons > 0) {
      if(UPD1908 and $cla =~ /parthi|frag0aa/) { } # maybe also skip aacons rescue of poor-noclass
      else { $keepdrop .="okay:aacons$aacons,"; }
      $aadupflag .= ",aacons:$aacons";  #? remove aadup:mainid or not
    }
  }

  # FIXME1805: rescue short/part cds with good/best match to refgenes
  if($cla =~ /parthi|frag0aa/) {
    $keepdrop.= "drop";
    
  } elsif($cla =~ /althi/) {
    if(UPD1908) {  
        if($ispoor == kNONCODE) { $cla =~ s/a2//; $cla =~ s/hi1/hi/; $cla= $cla."nc"; $keepdrop.="okay"; } ## dont rescue NC alt if other poor qual
        $keepdrop.= (($ispoor & NOTPOORBAD) or ($ispoor and $pali >= 95))?"drop":"okay";
    } else {
        $keepdrop.= ($ispoor)?"drop":"okay";  # ispoor vs main size?
    }

  } elsif($cla =~ /main/) {
    if(UPD1908) { if($ispoor & kNONCODE) { $cla =~ s/a2//; $cla= $cla."nc"; $keepdrop.="okay"; } } # rescue NC main if other ispoor
    $keepdrop.= ($ispoor)?"drop":"okay";   
    
  } elsif($cla =~ /noclass/) {
    if(UPD1908) { if($ispoor == kNONCODE) { $cla= $cla."nc"; $keepdrop.="okay"; } } # dont rescue NC noclass if other ispoor
    $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop? ditto to altmid
    
  } else { # other altmid/low 
    $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop this class, ispoor maybe uniq prot: maybeok?
  }
  
  my $okay;
  if($keepdrop =~ /drop/) {
    if($keepdrop =~ /okay/) { $okay= "maybeok"; } else { $okay= "drop"; }
  } else {
    $okay= "okay";
  }
  
  $tbits.= $aadupflag if($aadupflag); # defer, AFTER refbest/good/..
  $tbits.= ",chrmap:$mapqual,$maploc" if($maploc and $mapqual);
  $tbits.= ",pflag:".int($ispoor); # DEBUG info, but keep
  $tbits="aaref:$tbits" unless($tbits=~/^0,0/); #FIXME: add tag to tbits output
  
  # if(defined $eqflag{$tid}) { # TEST1603
  #   my $eqfl= $eqflag{$tid}{$qid}||""; 
  #   if($eqfl) { $eqfl="$qid/$eqfl,"; }
  #   my @q= grep{ $_ ne $qid } sort keys %{$eqflag{$tid}}; 
  #   $eqfl .= join",",map{ "$_/".$eqflag{$tid}{$_} }@q;  
  #   $tbits.= ",feq:$eqfl";
  # }
  
  return (wantarray) ? ($tid,$okay,$cla,$qid,$pidal,$tqual,$tbits) : $okay;
}



#p4o #========================= PRE4a ============================
#p4o 
#p4o =item classifytr
#p4o 
#p4o   classifier of locus primary, alternate, fragment, redundant 
#p4o   using identityclass() collection of overlapping transcripts (CDS or full tr)  
#p4o   focused on CDS qualities, maybe should be classifyCDS() ..
#p4o   
#p4o =cut
#p4o 
#p4o sub classifytr_PRE4a {
#p4o   my($tid,$cla,$qid,$pidal)= @_;
#p4o 
#p4o   # use constant { kAATINY => 1, kAAUTRBAD => 2, kAAUTRPOOR => 4, kAADUP => 8, kAAGAPS => 16, kNONCODE => 32, };
#p4o   # use constant NOTPOORBAD => kAATINY + kAADUP + kAAGAPS; # - kAAUTRPOOR - kAAUTRBAD - kAAGAPS
#p4o   # use constant NOTPOOR => kAATINY + kAADUP + kAAUTRBAD + kAAGAPS; # - kAAUTRPOOR
#p4o   $qid||="0"; $pidal||="0";
#p4o 
#p4o =item UPD1908 classifytr drop valid alts problem
#p4o   
#p4o     drops of 'true' reference genes, human, plant, are mainly 'althi' with 2 pflag qualities:
#p4o       aadup (>= 99% ident from cdhit.aa cluster)
#p4o       utrpoor/bad
#p4o     large portion have pia = 100/100/ but signif num have 100/99..80..small %align-global, should treat those as valids?  
#p4o     maybe drop aadup criteria, or at least balance it w/ other rescue criteria: pi/pa < 100/100 ? 
#p4o     add criteria of total non-ident cds-bases? some threshold (10? 100?) keeps alt?
#p4o 
#p4o     pflag, see above constants 2=utrbad 4=utrpoor 8=aadup
#p4o     see also FIXME 201402 drop.mains problem
#p4o       AAPART cutoff is problem, with UTRBAD/POOR, causing drop.mains of uniq orthologs 
#p4o       AAPART & AAMIN, AAMINBAD, AAMINPOO  need fix to prevent drop uniq/valid main/noclass genes, 
#p4o   
#p4o     drop.parthi has similar loss of valid ref alts, adjust parthi criteria?
#p4o       frag/part == $tsize < $ALTFRAG*$qsize; ALTFRAG = 0.50 default, should retest
#p4o       
#p4o     drop valid alts via perffrag: this happens when true alt is perfect substring of longer alt,
#p4o       not sure this is likely alternating-exon gene structure, but may be,
#p4o       more likely is experimentally found alternate start site in same exon structure.
#p4o       
#p4o     for ref_arath16ap.trclass
#p4o       5167 drop.althi, 3695 have aadup, 1445/1472 have utrpoor/bad of no-aadup
#p4o       5601 okay.althi, 0 have aadup, 1062 have utrpoor/bad
#p4o 
#p4o     for bevgc_evg3weed1rd.trclass
#p4o       85208 drop.althi, 41846 have aadup, 30681/43362 have utrpoor/bad, 14032/43362 aapartial of no-aadup
#p4o       106526 okay.althi, 0 have aadup, 32274 have utrpoorbad, 22851 aapartial
#p4o 
#p4o     for ref_human.trclass
#p4o       26715 drop.althi, 14346 aadup, 12305/12369 have utrpoor/bad of no-aadup
#p4o       22563 okay.althi, 0 aadup, 5975 have utrpoor/bad
#p4o     
#p4o     for trev19human.trclass
#p4o       126070 drop.althi, 45301 aadup, 37126/80769 utrpoorbad of noaadup, 35811/80769 aapartial of noaadup
#p4o       138448 okay.althi, 0 aadup, 39146 utrpoorbad, 30744 aapartial
#p4o       
#p4o =cut
#p4o   
#p4o use constant ALTPOOR1607 => 1;
#p4o =item ALTPOOR1607
#p4o    9jul16: alt class bug: drop.althi when <30% align to main, this should be valid alt due to small align
#p4o    problem ispoor from utrpoor, but is valid unique aahoref, possible paralog not alt
#p4o    cornhi8mtrinLocDN35783c0g3t2    drop    althi   cornhi8m9agvelvk35Loc6006t1     100/32/.        262,54%,complete-utrpoor        0,0,pflag:4
#p4o    test solution: ignore ispoor when $pal < 75% ?
#p4o =cut
#p4o     
#p4o   ## pidal == pident/palign,otherattr ==  "$pid/$pal$antiflag";
#p4o   my($pidi,$pali)= $pidal =~ m/^(\d+).(\d+)/;
#p4o   
#p4o   my $aw= $aasize{$tid} || 0;
#p4o   my $tqual= $aaqual{$tid} || ""; #? parse for "aasize,pcds,aaqual"
#p4o   ## add these utr class quals : this should be in tqual string: aasize, pcds, aaqual
#p4o   my $tw= $trsize{$tid} || 1; 
#p4o   my $pcds= int(300*$aw/$tw);  
#p4o   my $butr= ($tw - 3*$aw); # okutr if butr<=300p, modify BAD_CDSUTR for small aw
#p4o   
#p4o   my $ispoor= ($aw < $AAMIN or ($tqual =~ m/partial/ and $aw < $AAPART))?kAATINY:0;
#p4o 
#p4o   my $tbits= $aablast{$tid} || "0,0";
#p4o   my($tbscore,$tbref)=split",",$tbits;
#p4o   #FIXME: add tag to tbits output:  $tbits="aaref:$tbits" unless($tbits=~/^0,0/);
#p4o   
#p4o   my $aacons=0;
#p4o   if(AACONS) {
#p4o     $aacons= $aconsensus{$tid}||0;
#p4o     ## $ispoor= 0 if($aacons>0); # FIXME;
#p4o   }
#p4o   
#p4o   #TEST1602 add on tbits, like aablast ref info? equalgene mapqual info 
#p4o   my $mapqloc= ($nineqgene) ? $gmapqual->{$tid} : 0;  
#p4o   my($mapqual,$maploc)= ($mapqloc)?(split"\t",$mapqloc):(0,0);
#p4o   # $gmapqual{$mid} == "$mapqual\t$loc"; 
#p4o   # new bleqgene mapqual: 85a,100i,2469l,4x,Spl:29%,KB668900 == %align,%ident,bases,exons,Spl:split
#p4o 
#p4o   # change to aablastref{tid} ?
#p4o   my($tbrscore,$tbrefbest)= ($tbref) ? (split",", $aablastref{$tbref}) : (0,0); 
#p4o   my $rbits2= $aablastref{$tid} || "0,0";
#p4o   my($tbrscore2,$tbref2)= split",", $rbits2; # this is revhash tid => score,ref refbest
#p4o 
#p4o   my $aaclus= $aacluster{$tid} || "0,0";  
#p4o   my($aamainid,$aamainpi)= split",", $aaclus; # ="$mainid,$pi";
#p4o   if($aamainid eq $tid) {  $aamainid=""; $aamainpi=0; }
#p4o   ## FIXME: aamain can be bad/missing; fix before this, need cds/tr align info to know if aamainid is good.
#p4o   
#p4o   ## main class should use aasize, aapartial and utrbad/poor with AAMINBAD, AAMINPOO
#p4o   ## FIXME: pcds bad here? comes from aw - gaps effects; separate gapbad from utrbad
#p4o   ## 201402: maybe ignore tqual utrbad/poor and use only pcds with these cutoffs? dont know tqual cutoffs
#p4o   if($tqual =~ m/gapbad/) { $ispoor |= kAAGAPS; }  # gaps discrepancy w/ aaqual utrbad and pcds utrbad 
#p4o   if($butr >= $MINUTR) { # ~ 300b "fixed" average utr sizes
#p4o     # should adjust %cds bad for cds size, for large cds> 10k, pcds >=90%, for small cds < 300, pcds may be 33%
#p4o     if($tqual =~ m/utrbad/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRBAD; } 
#p4o     elsif($tqual =~ m/utrpoor/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRPOOR; }
#p4o     else { $ispoor |= ($pcds < $BAD_CDSUTR)? kAAUTRBAD: ($pcds < $OK_CDSUTR)? kAAUTRPOOR : 0; }
#p4o     # if($pcds < $BAD_CDSUTR) { $ispoor |= kAAUTRBAD; } elsif($pcds < $OK_CDSUTR) { $ispoor |= kAAUTRPOOR; } 
#p4o   }
#p4o 
#p4o   my $tcode=0;
#p4o if(CODEPOT1607) {
#p4o   $tcode= $sizeval{$tid}{'codepot'}||""; # $tcode =~ /^Noncode/ set what? kAATINY+kAAUTRBAD
#p4o   if($tcode and $tcode =~ /^Noncode/) { # /^N/ enough
#p4o   $tcode=0 if($tbscore>0);
#p4o   # $tcode=0 if($aw > $AAPART);# UPD1908 TESTING,not very accurate but use to separate too many tiny aa
#p4o   $tcode=0 if($aw > 139 or ($aw > 99 and not $ispoor)); #?? which? many utrorf/utrbad for aa>99
#p4o   $ispoor |= kNONCODE if($tcode and $tcode =~ /^Noncode/);  
#p4o   }
#p4o }
#p4o 
#p4o   ## TEST1603: $pidal =~ m/altmap\d+xeq/ ; xeq for althi1 means redundant == AADUP
#p4o   # problem cases, xeq but not same due to poor align, diff aa size, etc. : want to keep both these, t2 is valid alt
#p4o   # Anofunz4hEVm000070t1	okay	althi	Anofunz4hEVm000070t2	100/52/./altmap87xeq	3437,61%,complete	aaref:4110,AGAP000222-PB,refgood,chrmap:100a,99i,10314l,8x,KB668689:297566-312085:-,pflag:0
#p4o   # Anofunz4hEVm000070t2	okay	main	Anofunz4hEVm000070t1	100/52/./altmap87xeq	3511,52%,complete-utrpoor	aaref:4064,AGAP000222-PB,refgood,chrmap:87a,99i,10536l,8x,KB668689:297566-312085:-,pflag:0
#p4o   
#p4o   if(UPD1908) { 
#p4o     # ?? need finer filter for althi than current ($aamainpi > $AADUP_IDENT)
#p4o     # bump AADUP_IDENT from 98 to 99?  needs testing
#p4o     # losing too many valid alts, but still have too many invalids
#p4o   }
#p4o   if($aamainpi > $AADUP_IDENT) { # TEST aacluster added info on dup prots; AADUP_IDENT=98 default
#p4o     ## check that aamainid is main ??
#p4o     $cla.="a2" unless($cla=~/hi[012]/); # hi1/a2 redundant info
#p4o     $ispoor |= kAADUP if($cla =~ /althi/); # only for /althi|parthi/
#p4o   }
#p4o   unless( $ispoor & kAADUP ) {
#p4o     if($pidal =~ m/altmap\d+xeq/ and $cla =~ /althi1|part/) {
#p4o       $ispoor |= kAADUP;  # add new bitflag kDUPEXONS ?
#p4o     }
#p4o     # V4: eqgene altmap100 fixup, these have identical CDS set..
#p4o     # if($pidal =~ m/altmap(\d+)/) { my $am=$1; 
#p4o     #   if($am >= $AADUP_IDENT and $cla =~ /althi1|part/) {
#p4o     #     $cla.="a2" unless($cla=~/hi[012]/); # hi1/a2 redundant info
#p4o     #     #? $tbits.=",aadup:$aamainid";
#p4o     #     $ispoor |= kAADUP;  # add new bitflag kDUPEXONS ?
#p4o     #   }
#p4o     # }
#p4o     # if(($maploc =~ /\.icdup/) and ($cla =~ /alt|part/)) { $ispoor |= kAADUP; }  #? kDUPEXONS
#p4o   }
#p4o   my $aadupflag= ($ispoor & kAADUP and $aamainid) ? ",aadup:$aamainid" : ""; # defer, AFTER refbest/good/..
#p4o   # below add to  $tbits.= $aadupflag
#p4o 
#p4o   # TEST3 : add? ispoor for cla == althi1 and pid/pal == 100/100  ie identicals?? NO, pal is align/shortlen
#p4o   # NO, not this, althi1 can be good alt, need align/mainlen stat also. or use only AADUP score as now.
#p4o   
#p4o   # CHECKME: adding aablast kept 40k more okay, all althi/ahia2 + 2k parthi
#p4o   # .. are these true uniq aablast or just althi aadups ?
#p4o   my $MINBLASTSCORE= 60; # aablast only? bitscore always?
#p4o   my $keepdrop="";
#p4o   unless( $tbscore == 0 or $tbits=~/^0,0/) {
#p4o     my $risbest= (($tbrefbest eq $tid) or ($tbscore >= $tbrscore))?1:0;# was ($tbrefbest eq $tid), allow 2+ same score
#p4o     my $risgood= ($tbrscore > $MINBLASTSCORE and $tbscore >= 0.90 * $tbrscore)?1:0;
#p4o     my $risok= ($tbref2 and $tbrscore2 >= $MINBLASTSCORE)?1:0;
#p4o     if($risok and $risgood) { $risgood=0 if($tbrscore2 > $tbrscore); } #?
#p4o     $ispoor = $ispoor & NOTTINY if($risgood); ## FIXME1712: clear kAATINY for $risok/risgood
#p4o     ## ispoor & kAADUP should not be cleared here???
#p4o     my $isaadup=($ispoor & kAADUP);
#p4o     if($risbest) { 
#p4o       $tbits.=",refbest"; $keepdrop.="okay:refbest,"; $ispoor=0 unless($isaadup); }
#p4o     elsif($risgood) { # was 3rd
#p4o       $tbits.=",refgood"; unless($isaadup){ $keepdrop.="okay:refgood,"; $ispoor=0;} } # maybe2 ok ?
#p4o     elsif($risok) {  # was 2nd ? why superceed refgood? 2nd ref match?
#p4o       $tbits.=",refok"; unless($isaadup){ $keepdrop.="okay:refok,"; $ispoor=0;} }
#p4o   }
#p4o   
#p4o       # 201402: update for drop.mains, perfectfrag replacements, need input list of perfectfrag/dups from fastanr/cdhit prior steps
#p4o       # 201402: AAPART & AAMIN, AAMINBAD, AAMINPOO  need fix to prevent drop uniq/valid main/noclass genes, 
#p4o   ## maybe ignore utrbad for main .. otherwise can miss true longutr/ncutr genes; but refbest/good will keep those w/ homology
#p4o   ## keep largeaa,utrpoor for main, noclass, altlow; but drop smallaa,utrpoor
#p4o   if($ispoor > kAATINY and $cla !~ /althi|parthi/) {  
#p4o     if(UPD1908) { } # maybe rescue althi also here
#p4o     if($aw >= $AAMINBAD) { $ispoor = $ispoor & NOTPOORBAD; } ##  ^ (kAAGAPS+kAAUTRPOOR+kAAUTRBAD);  ##?? ^ XOR bad
#p4o     elsif($aw >= $AAMINPOO) { $ispoor = $ispoor & NOTPOOR; } ##  ^ (0+kAAUTRPOOR);
#p4o   }
#p4o 
#p4o   if(AACONS) {
#p4o     # $aacons= $aconsensus{$tid}||0;
#p4o     # $ispoor= 0 if($aacons>0); # NO;
#p4o     if($aacons > 0) {
#p4o       if(UPD1908 and $cla =~ /parthi|frag0aa/) { } # maybe also skip aacons rescue of poor-noclass
#p4o       else { $keepdrop .="okay:aacons$aacons,"; }
#p4o       $aadupflag .= ",aacons:$aacons";  #? remove aadup:mainid or not
#p4o       ## $ispoor & kAADUP : remove or not? NOT, has other uses/meaning
#p4o     }
#p4o   }
#p4o 
#p4o   # FIXME1805: rescue short/part cds with good/best match to refgenes
#p4o   # .. fragment/part class may include distinct loci, matching ref genes, but also w/ align to longer locus
#p4o   # .. prior  keepdrop= okay:refbest should keep; check for losses
#p4o   if($cla =~ /parthi|frag0aa/) {
#p4o     $keepdrop.= "drop";
#p4o     # not frag? if($ispoor & kNONCODE) { $cla= $cla."nc"; $keepdrop.="okay"; } 
#p4o     
#p4o   } elsif($cla =~ /althi/) {
#p4o if(UPD1908) {  
#p4o     # valid althi loss: maybe bump ($pali > 75) to ($pali > 95) ? needs test
#p4o     if($ispoor == kNONCODE) { $cla =~ s/a2//; $cla =~ s/hi1/hi/; $cla= $cla."nc"; $keepdrop.="okay"; } ## dont rescue NC alt if other poor qual
#p4o     $keepdrop.= (($ispoor & NOTPOORBAD) or ($ispoor and $pali >= 95))?"drop":"okay";
#p4o } elsif(ALTPOOR1607) {
#p4o     $keepdrop.= (($ispoor & NOTPOORBAD) or ($ispoor and $pali > 75))?"drop":"okay";
#p4o } else {
#p4o     $keepdrop.= ($ispoor)?"drop":"okay";  # ispoor vs main size?
#p4o }
#p4o 
#p4o   } elsif($cla =~ /main/) {
#p4o     # $keepdrop.= "okay"; # always ok if has alts? NO, 
#p4o     # Fixme: keeping all "main" class gives ever expanding trsets with added trasm
#p4o     # should apply aaqual, aaref criteria to main + its althi; 
#p4o 
#p4o     # $ispoor = $ispoor & kAATINY; # clear bits 2,4,.. # now above
#p4o     if(UPD1908) { if($ispoor & kNONCODE) { $cla =~ s/a2//; $cla= $cla."nc"; $keepdrop.="okay"; } } # rescue NC main if other ispoor
#p4o     $keepdrop.= ($ispoor)?"drop":"okay";   
#p4o 
#p4o     ## 21feb update: dmag5icd35all_cde35.class4
#p4o     ## class7 maindrops: utrpoor/bad:48144, then partial:19493, gapbad:2866,  397 are aalong+complete
#p4o     ## class6 maindrops: utrpoor/bad:57636, then partial:19815, gapbad:3178,  1075 are aalong+complete
#p4o     # dmag4vel4ibnk31Loc11665t1       drop    main    dmag4vel4ifik31Loc13364t1       99/99   102,71%,complete        0,0,pflag:4
#p4o     # dmag4vel4ibnk31Loc1189t4        drop    main    dmag4vel4ipak31Loc21101t1       100/35  364,62%,complete        0,0,pflag:18
#p4o     ## pflag:18 = kAAGAPS + kAAUTRBAD
#p4o     #..
#p4o     ## class4 maindrops: utrpoor/bad:60742, then partial:19737, but 3000 are aalong+complete, WHY drop?
#p4o     ## class5 maindrops: utrpoor/bad:56971, then partial:18807, but 1788 are aalong+complete (have gaps, so utrbad/poor now)
#p4o     # dmag4vel4ibnk31Loc10846t1       drop    main    dmag4vel4ipak31Loc13986t1       99/94   116,73%,complete        0,0,pflag:4
#p4o     # dmag4vel4ibnk31Loc11089t1       drop    main    dmag4vel4ifik31Loc18962t1       99/68   129,70%,complete        0,0,pflag:4
#p4o     # >dmag4vel4ibnk31Loc10846t1 is aacluster:main, not aatiny, not utrpoor, not kAADUP, not althi;
#p4o     # .. its alt is poor: partial5-utrbad, 37aa; did AAcluster fix mangle main?  pflag:4 == kAAUTRPOOR, miscalc from $pcds ???
#p4o     # ** GAPS in aa above; 116aa,73%pcds is wrong if -gaps removed. BUT gaps not removed from nt size, so pcds off by that.
#p4o     #  dmag4vel4ibnk31Loc10846t1 68aa,48gap,475nt,pcds=3*68/475=43%;  pcds= 3*68/(475-3*48) = 62%
#p4o     # .. so drop.main/aacomplete is ~correct calc given aagaps : probably dont want utrbad flag here, but a gapbad may be useful
#p4o     # .. 48gaps/116aa is not good prot; add flag gapbad == kAAGAPS above in readAAqual ?
#p4o 
#p4o     
#p4o   } elsif($cla =~ /noclass/) {
#p4o     # UPD1908 maybe modify aacons rescue for poor-noclass, as with parthi/frag
#p4o     if(UPD1908) { if($ispoor == kNONCODE) { $cla= $cla."nc"; $keepdrop.="okay"; } } # dont rescue NC noclass if other ispoor
#p4o     $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop? ditto to altmid
#p4o     
#p4o   } else { # other altmid/low 
#p4o     ## code err here, in $cla?: Can't call method "other" without a package or object reference
#p4o     ##? if(UPD1908) { if($ispoor == kNONCODE) { $cla= "$cla"."nc"; $keepdrop.="okay"; } } dont rescue NC noclass if other
#p4o     $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop this class, ispoor maybe uniq prot: maybeok?
#p4o   }
#p4o   
#p4o   my $okay;
#p4o   if($keepdrop =~ /drop/) {
#p4o     if($keepdrop =~ /okay/) { $okay= "maybeok"; } else { $okay= "drop"; }
#p4o   } else {
#p4o     $okay= "okay";
#p4o   }
#p4o   
#p4o   ## AACONS: add flag here, like aadup:
#p4o   $tbits.= $aadupflag if($aadupflag); # defer, AFTER refbest/good/..
#p4o   $tbits.= ",chrmap:$mapqual,$maploc" if($maploc and $mapqual);
#p4o   $tbits.= ",pflag:$ispoor"; # DEBUG info, but keep
#p4o   $tbits="aaref:$tbits" unless($tbits=~/^0,0/); #FIXME: add tag to tbits output
#p4o   
#p4o   if(defined $eqflag{$tid}) { # TEST1603
#p4o     my $eqfl= $eqflag{$tid}{$qid}||""; 
#p4o     # debug output all eqflag{tid}{otherid}
#p4o     if($eqfl) { $eqfl="$qid/$eqfl,"; }
#p4o     my @q= grep{ $_ ne $qid } sort keys %{$eqflag{$tid}}; 
#p4o     $eqfl .= join",",map{ "$_/".$eqflag{$tid}{$_} }@q;  
#p4o     $tbits.= ",feq:$eqfl";
#p4o   }
#p4o   
#p4o   return (wantarray) ? ($tid,$okay,$cla,$qid,$pidal,$tqual,$tbits) : $okay;
#p4o }
#p4o 
#p4o 
#p4o =item identityclass
#p4o 
#p4o   ## UPD1912: evigene/scripts/rnaseq/asmrna_identityclass4a.pl
#p4o 
#p4o   from  aabugs4qual/tsaevg/cdsidentclass.sh
#p4o   input from above putspans, sorted by ^Qclen, ^Align, vTtlen:
#p4o      qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits) 
#p4o 
#p4o   update classes per aabugs4qual/tsaevg/classparts.pl
#p4o   .. need blastp/blast table input for this..
#p4o   .. need 3-4 final categories:  keep-main, keep-alts, trash-fragments (alt,noalt) 
#p4o      .. keep includes partial prot where good enough.
#p4o      .. trash-alts maybe (not fragment but not distinct enough)
#p4o      .. add poor cut to main,alt,noclass per below a,b qualities.
#p4o      
#p4o cl_main.alttr and cl_alts, keep if: 
#p4o   a. homology unique or best (near best for alt),  or
#p4o   b. aasize >= 45% of main and aacomplete? and not aasame-as-main?
#p4o   c. main == fragment if aasize < MINAAFULL or (aasize < MINAAPART and partial/utrpoor)
#p4o 
#p4o FIXME/Check : aablast score has dropped/poor a set of 653 uniq ref (Nref1) matches vs prior clidtab/class1
#p4o   see env swapids=1 refblast=aaeval/refdebl/refde-.tall3 ./classparts.pl trsets/$pt.clid2tab
#p4o Partition               Nclass  Nmatch  Nrefa   Nref1   Bits    Iden    Algn    Rlen    Tlen
#p4o locust1best5.cl_poor    68122   1331    3485    653     113     61      170     361     389
#p4o 
#p4o # # FIXME2: maybe ignore utrbad for main .. otherwise can miss true longutr/ncutr genes; 
#p4o # # but refbest/good will keep those w/ homology
#p4o # cat $pt.class3 | grep 'drop.main' | sort -k6,6nr | head
#p4o # fungrvel4ik53Loc89t53   drop    main    fungrvel3bk43Loc142t24  99/66   1408,27%,complete-utrbad        187.2,DRERI:ENSDARG00000086955
#p4o #.. ^only aa equiv is also drop:
#p4o #  fungrvel3bk43Loc142t24  drop    althi   fungrvel4ik53Loc89t53   99/66   1373,32%,complete-utrbad        181.4,DRERI:ENSDARG00000086955
#p4o #
#p4o # fungrvel4k25Loc3571t1   drop    main    fungrvel4k29Loc22240t1  99/100  668,28%,complete-utrbad 0,0
#p4o # fungrvel4k25Loc15515t1  drop    main    fungrvel2k35Loc34556t1  100/100 580,24%,complete-utrbad 261,DRERI:ENSDARG00000068192
#p4o # .. other DR92 genes:
#p4o #  fungrvel3bk35Loc16459t1 okay    main    fungrvel3bk29Loc13518t1 99/92   876,82%,complete        543,DRERI:ENSDARG00000068192,refbest
#p4o #  fungrvel3bk29Loc13518t1 okay    althi   fungrvel3bk35Loc16459t1 99/92   836,61%,complete        572,DRERI:ENSDARG00000068192,refbest
#p4o #  fungrvel4k25Loc15515t1  drop    main    fungrvel2k35Loc34556t1  100/100 580,24%,complete-utrbad 261,DRERI:ENSDARG00000068192
#p4o #  fungrvel4k35Loc25356t1  okay    althi   fungrvel4k25Loc15515t1  100/65  363,99%,partial3        256.3,DRERI:ENSDARG00000068192
#p4o #
#p4o # fungrvel3bk29Loc5288t1  drop    main    fungrvel4k25Loc5314t2   99/87   565,23%,partial5-utrbad 0,0
#p4o # fungrvel4k25Loc3268t8   drop    main    fungrvel4ik53Loc14576t1 100/58  526,24%,complete-utrbad 396,DRERI:ENSDARG00000038737
#p4o #...    
#p4o     
#p4o =cut
#p4o 
#p4o sub identityclass_PRE4a {
#p4o   my($outh, $infile, $insorted)= @_;
#p4o 
#p4o   use constant TEST1602 => 1; ## UPD 2016.02 problem w/ dup equal mains, equivalence after 1st see both as mains, 1 to be alt,
#p4o   use constant TEST3 => 1; # 13aug09 test fixes to alt/main classing
#p4o 
#p4o 	my $havevalidids= scalar(%validids)?1:0; ## need this?
#p4o 
#p4o ## 2016.feb : problems now classing  alts => main/noclass, tho blastab infile has proper hi identity scores
#p4o ##   found w/ genom-map eqgene, some essentially identical tr classed as sep loci ..
#p4o ##   repeated runs of this on new oksets reduces this alt>loc misclass, but not entirely gone.
#p4o 
#p4o ## 2013.aug : IS_CDSALIGN ORIENT in alntab, add -sign/antisense to one field
#p4o ### -aln == reverse align ?? for CDS antisense problem
#p4o ### -bits == reverse align ?? bits not used here for scoring.. best choice other than adding column
#p4o ### -aln affects here sort: -k7,7nr; use other field Ident?
#p4o 
#p4o #... see above now
#p4o #   $TINYALN=$ENV{mina}||50; 
#p4o #   $IS_CDSALIGN=$ENV{cdsw}||0; # default using tralign, tests better than cds
#p4o #   $PHIALN=$ENV{ahi}||98; # was 65; # use mina?
#p4o #   $PHI=$ENV{phi}||99; $PMID=$ENV{pmid}||90; $PLOW=$ENV{plow}||80; 
#p4o #   # $ALTFRAG: add isfrag pct; now 50 (0.5)
#p4o   
#p4o ### infile is putspans() alntab: Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits
#p4o ###  sort ^Qclen, ^Align, vTtlen: cat $atab | grep -v '^Qid' | sort -k2,2nr -k7,7nr -k6,6n 
#p4o ### should it be ^Qclen, ^QID, ^Align, vTtlen : to keep qids together? no want top Aligns first
#p4o ### * add vQtlen before vTtlen, else get utrbad before good ***
#p4o   # my $ALNSORTORD='-k2,2nr -k7,7nr -k6,6n '; #orig
#p4o   # my $ALNSORTORD='-k2,2nr -k7,7nr -k6,6n -k1,1 -k4,4'; #add IDs to order ties
#p4o 
#p4o   my $ALNSORTORD='-k2,2nr -k7,7nr -k3,3n -k6,6n -k1,1 -k4,4';  #add vQtlen/3 before vTt/6 IDs to order ties
#p4o   unless($IS_CDSALIGN) { # sort tlen, not clen 
#p4o     # ..not sure what is right, want cds to play role in best choice
#p4o     # .. input blast align k7 is for tlen, not clen, but want to choose best by long clen > long tlen,talign
#p4o     #t1: $ALNSORTORD='-k3,3nr -k7,7nr -k2,2nr -k6,6nr -k1,1 -k4,4';  # ^Qtlen,^Align,^Qclen,vTtlen,Qid,Tid
#p4o     $ALNSORTORD='-k2,2nr -k7,7nr -k3,3n -k6,6n -k1,1 -k4,4';  # test same as IS_CDSALIGN  
#p4o   }
#p4o   
#p4o   my($inh, %class,%bestmatch,%ismain, %havepair); 
#p4o   if(ref($infile)) { $inh= $infile; } # STDIN,.. sort ??
#p4o   else {
#p4o     ## FIXME fail /tmp sort use, what bug? using TMPDIR.. use -T ./
#p4o     my $tmpdir=$ENV{TMPDIR};
#p4o     $tmpdir= './' unless($tmpdir and $tmpdir ne "/tmp"); ## { $tmpdir=`pwd`; chomp($tmpdir); }
#p4o     if($insorted) { open(IN,$infile) or die "reading $infile"; $inh=*IN; }
#p4o     else { open(IN,"sort -T $tmpdir $ALNSORTORD $infile |") or die "reading $infile"; $inh=*IN; }
#p4o     #orig#else { open(IN,"sort $ALNSORTORD $infile |") or die "reading $infile"; $inh=*IN; }
#p4o     
#p4o     ## UPD 1602: maybe revise to drop sort, expect per-qd, top-hit sort order of blastn output?
#p4o     ## then change class assignments as new qd x td warrants?
#p4o   }
#p4o   
#p4o   my($lastd)=("");
#p4o   while(<$inh>) { # maybe local table, sorted, or from file
#p4o     next if(/^Qid|^\W/); chomp; 
#p4o     my @v= split"\t"; 
#p4o     my($qd,$qc,$qw,$td,$tc,$tw,$aln,$iden,$bits)= @v; 
#p4o     my($isfrag,$aclass,$alnmax,$pid,$pal,$samesize)= (0) x 10;
#p4o     
#p4o 		## 2013.aug : IS_CDSALIGN ORIENT in alntab, add -sign/antisense to one field: not aln, iden?, bits?
#p4o     ## FIXME: *** 98/100/-sense/PitaEaR000975t24 << buggers, mixes up users of this table.
#p4o     ## all piad entries should have pd and sense/asense fields, placeholder where needed '.' ?
#p4o     
#p4o     #old# my $antiflag= ($IS_CDSALIGN and $bits < 0) ? "/-sense" : ""; ## append to bestmatch ??
#p4o     my $antiflag= ($IS_CDSALIGN and $bits < 0) ? "/-sense" : "/."; ## append to bestmatch ??
#p4o     		## bestmatch="id,99/89/-sense" for antisense?
#p4o     		## ^^ is this new flag causing problems?  Besthit is appended /after, parsed where?
#p4o     
#p4o     #AACONS.2019 : special case of samesize  consensus(td) = rare (qd), swap  
#p4o     # NEEDS TEST .. if($aconsensus{$td} > $aconsensus{$qc}) .. swap qd,td
#p4o     
#p4o     $isfrag= $aclass="";
#p4o     $samesize=($tc == $qc)?1:0; # OR tinyalndiff? abs($tc - $qc) < $NHIALN
#p4o     my $swapqt=($tc > $qc)?1:0;
#p4o     use constant UPD1912doswapc => 0; # 0, turn this off, likely bug
#p4o     if(AACONS and UPD1912doswapc) {
#p4o       if($nconsensus and $samesize) {
#p4o         #UPD1912 : Is this swapqt a problem, for dup loci, almost identicals? changes their sort order some what randomly
#p4o         my $tdcon= $aconsensus{$td}||0; my $qdcon=$aconsensus{$qd}||0;
#p4o         $swapqt=1 if($tdcon > $qdcon);
#p4o       }
#p4o     }
#p4o     if($swapqt){ ($qd,$td)=($td,$qd); ($qc,$qw,$tc,$tw)=($tc,$tw,$qc,$qw); } # swap to greater cds?
#p4o 
#p4o     if(EQGENE_CHANGES_NOALN && TEST1602) {
#p4o       ## neqalts data not yet clean enough for this reclassing
#p4o       ## ? just ignore cds-align row for those marked as non-alts? 
#p4o       ## set any vals? class{},ismain{},havepair{} ?>
#p4o       #?? next if(ref $neqalts and $neqalts->{$qd}{$td}	);  
#p4o 
#p4o       # if(TEST1602 && $lastd && $lastd ne $qd) {	# reverse, neqalts	OR below/above: skiptonext if($neqalts and $neqalts->{$qd}{$td}	);	
#p4o       # if(ref $neqalts) { 
#p4o       #   my @ted= sort{ $aasize{$b} <=> $aasize{$a} or $a cmp $b } keys %{$neqalts->{$lastd}}; 
#p4o       #   foreach my $te (@ted) { 
#p4o       #     if($havepair{$lastd}{$te}) { # break pairing alt class
#p4o       #       my $laclass="mainmap"; #? mainparalog?
#p4o       #       if( $class{$te} =~ /^alt/) {  $class{$te}= $laclass; }
#p4o       #       #?? elsif( $class{$lastd} =~ /^alt/) {  $class{$lastd}= $laclass; }
#p4o       #       $havepair{$lastd}{$te}=0; 
#p4o       #       $havepair{$te}{$lastd}=0;
#p4o       #     } 
#p4o       #   }
#p4o       # }
#p4o       # }				  
#p4o     }
#p4o     
#p4o 
#p4o =item eqgene classifier
#p4o 
#p4o   - do this in readEqualGene, mark which overlap alts are bad, which ok
#p4o   - also need to regard mapqual align, Split values to decide
#p4o   
#p4o   need separate classifier to handle various eqgene attributes, decide which tr/alts are bad/good
#p4o   eg eqgene classifier for this case: good g131, bad g453t5
#p4o   .. for this case need to know that g453t5 <mismap> g453t1,2,3,4 .. count overlap/alt/locus? drop outliers? use Split info?
#p4o ok Anofunz4gEVm000131t1	noid	Anofunz4gEVm000452t5/19	KB668936:299087-306372:-	99a,99i,7287l,3x	0
#p4o ok Anofunz4gEVm000131t2	noid	Anofunz4gEVm000452t5/33	KB668936:299737-303981:-	98a,99i,4269l,2x	0
#p4o ok Anofunz4gEVm000131t3	noid	Anofunz4gEVm000452t5/58	KB668936:299632-301584:-	70a,99i,1809l,2x	0
#p4o ok Anofunz4gEVm000452t1	noid	na	KB668900:3696-9724:-	99a,100i,4410l,8x	0
#p4o ok Anofunz4gEVm000452t2	noid	na	KB668900:3696-9724:-	99a,100i,4410l,8x	0
#p4o ok Anofunz4gEVm000452t3	noid	na	KB668900:3696-9724:-	96a,100i,4329l,8x	0
#p4o ok Anofunz4gEVm000452t4	noid	na	KB668900:3696-8875:-	100a,100i,3372l,9x	0
#p4o bad Anofunz4gEVm000452t5	noid	Anofunz4gEVm000131t1/56,Anofunz4gEVm000131t2/56,Anofunz4gEVm000131t3/43	KB668936:300530-301921:-	85a,100i,2469l,4x,Spl:29%,KB668900	0
#p4o ok Anofunz4gEVm000452t6	noid	na	KB668900:6592-8636:-	79a,100i,948l,3x	0
#p4o 
#p4o   g131,3/3 alts are over 1/6 g453 alts
#p4o   -- should this indicate g453t5 is not-much-over other g453t alts?
#p4o   
#p4o =cut
#p4o 		  
#p4o 		if(EQGENE_CHANGES_NOALN && $lastd && $lastd ne $qd) { 
#p4o 		  ## revise lastd class if warranted, when have new qd
#p4o 		  #* change this altmap class: NO cds-align b/n td,qd means something, poor mapping paralogs end up here
#p4o 		  # .. need (a) gmap qual score (align,ident,split), and
#p4o 		  # ..      (b) paralog flag reverse of eqgene, for classed-alts that *dont* gmap same locus
#p4o 		  # TEST1602: ($eqgenes,$nineqgene,$neqgene,$gmapqual,$neqalts)= readEqualGene($eqgene) extended table
#p4o 
#p4o 			if($neqgene>0 && $eqgenes->{$lastd}) {
#p4o 				my @ted= sort{ $aasize{$b} <=> $aasize{$a} or $a cmp $b } keys %{$eqgenes->{$lastd}}; 
#p4o       
#p4o         ## FIXMEd remove dupl map ids _G2,3.., from eqgene table: Funhe2Exy3m149508t1_G2 ..
#p4o         ## FIXME2: this turns all main class into altmap, need to keep 1 as main or otherwise flag main/altmap
#p4o         ## check each lastd/@ted for ismain{id}, class{id}, and aasize{id} : pick one to keep as main
#p4o 			  
#p4o 			  ## FIXME3: Still bad, 160217, replaces combined locus groups' main w/ althi1/altmap w/o reciprocal alt to main
#p4o 			  
#p4o 			  ##?? decide here w/ mapquals to reclass? or should equalgene table be revised for lowqual cases?
#p4o 			  ## mapqual now == gmapalign%,nnn,location
#p4o 			  # my $ldmapqual= (ref $gmapqual)? $gmapqual->{$lastd} : ""; # TEST1602
#p4o 			  # skip_reclass if($ldmapqual < $MINGMAPQUAL);
#p4o 			  
#p4o 				my $miss=0;
#p4o 				foreach my $te (@ted) { 
#p4o           unless($havepair{$lastd}{$te}) { 
#p4o             next if($te =~ /_G\d+$/); # dup map id syntax, dropped from eqgene table now?
#p4o             next if($havevalidids and not $validids{$te}); #?? fix for _Gnnn, or not? 
#p4o             ## should do in readEqualGene but dont have %validids then
#p4o 
#p4o  			      # my $temapqual= (ref $gmapqual) ? $gmapqual->{$te} : ""; # TEST1602
#p4o  			      # next if($temapqual < $MINGMAPQUAL);
#p4o  			      
#p4o             my $palmap= $eqgenes->{$lastd}{$te}||0;  # val is te align to lastd, not reverse
#p4o             my $palmapt= $eqgenes->{$te}{$lastd}||0;   
#p4o             # swap ids for largest pal ??
#p4o             if($palmap > 0) {
#p4o               $miss++; 
#p4o               #bad?# my $pidalnval="99/$palmap";  # .$fakeantiflag ??
#p4o               my $pidalnval="$palmap";  
#p4o              
#p4o               ## this is wrong way, te here is (often) main gene, lastd is fragment overlapping 
#p4o               ## maybe NOT bestmatch{lastd}, val is te align to lastd
#p4o               ## maybe replace bestmatch{te} if palmap > bestmatch{te} ?
#p4o 
#p4o               $bestmatch{$lastd}="$te,$pidalnval" unless($bestmatch{$lastd});
#p4o 
#p4o               #BAD??# $bestmatch{$te}="$lastd,$pidalnval" unless($bestmatch{$te});
#p4o               
#p4o               my $laclass="altmap"; # this may be replacing main wrongly, need mainmap? test class{te} first?
#p4o               #bad# $class{$te}= $laclass; #?? is this going to work
#p4o               $class{$lastd}= $laclass; 
#p4o               $havepair{$lastd}{$te}++; $havepair{$te}{$lastd}++; #??
#p4o               }  
#p4o             } 
#p4o 					}
#p4o 				if($miss) { } #.. reclass $lastd ??
#p4o 			} # eqgene
#p4o 		}
#p4o 		
#p4o 				
#p4o 		$lastd=$qd; # here?
#p4o     if($td eq "self" or $qd eq $td) { # putspan fix for no-align cases
#p4o       $td=$qd; # $aclass="noclass"; << should be this, but ??
#p4o       $bestmatch{$qd}="$td,100/100/self1" unless($bestmatch{$qd});
#p4o       next;  # can lose eqgene if no further td/qd ..
#p4o     }
#p4o 
#p4o     #FIXME: have 2nd perfect matches, e.g qd1 >> td but qd2 == td; need that class to drop dups.
#p4o     # ** BETTER: Remove these before align/cluster; fastanrdb on .cds, .aa?
#p4o     #  asmrna5/cdsx/alldmag5x.cds n=4093166; alldmag5x.nrcds n=1782663; 501494 have dups (some many-many)
#p4o     # .. info is in blastn aligns, not in cdhit clusters (no 2ndary aligns).
#p4o     # .. use alntab ident column, if cdssize1 = cdssize1 = ident-3 (stop), identicals
#p4o     # eg dmag4vel4xbxk55Loc9866t1  dmag4vel4xfik55Loc11404t1 dmag4vel4xpak25Loc1083t9 : cds-ident alts of main dmag4vel4ibxk55Loc9119t1
#p4o     
#p4o     ## add here?? $aaqual{$qd} ; $aaqual{$td}; .. maybe do this below, END
#p4o     #  my($qqual,$tqual)= ($aaqual{$qd},$aaqual{$td});
#p4o     #  my($qbits,$tbits)= ($aablast{$qd},$aablast{$td});
#p4o 
#p4o     ## FIXME: check qw-main pCDS for UTRBAD/POOR; dont call tw frag if qc/qw is UTRBAD/POOR
#p4o     ## FIXME1405: sometimes? use input $aaqual{$qd}, ignore OK_CDSUTR/BAD_CDSUTR unless option..
#p4o     ## my $maw = $aasize{$qd} || 0;
#p4o     ## my $mqv = aaqualscore($aaqual{$qd});  #? parse for "aasize,pcds,aaqual"
#p4o 
#p4o     if(0) {
#p4o     ## NOT YET USED: add these utr class quals?
#p4o     my $qcds= ($qw>0) ? 100*$qc/$qw : $OK_CDSUTR;
#p4o     my $tcds= ($tw>0) ? 100*$tc/$tw : $OK_CDSUTR;
#p4o     my $qutrbad= ($qcds >= $OK_CDSUTR)?0:1; # ($qcds > $BAD_CDSUTR)?1: 2;
#p4o     my $tutrbad= ($tcds >= $OK_CDSUTR)?0:1; # ($tcds > $BAD_CDSUTR)?1: 2;
#p4o     }
#p4o     
#p4o     # FIXME: when qc == tc, can assign althi1 = althi2 instead of main1 = althi2
#p4o     my($qsize,$tsize)= ($IS_CDSALIGN) ? ($qc,$tc) : ($qw,$tw);
#p4o     # note: alnmax is min(qsize,tsize) not max
#p4o     $alnmax= ($qsize>$tsize and $tsize>0)?$tsize:($qsize>0)?$qsize:$tsize;  
#p4o     $isfrag= ($tsize < $ALTFRAG*$qsize)?"frag":"";
#p4o     
#p4o     #o if($IS_CDSALIGN) { $alnmax=($qc>$tc and $tc>0)?$tc:($qc>0)?$qc:$tc;  $isfrag= ($tc < $ALTFRAG*$qc)?"frag":""; } 
#p4o     #o else { $alnmax=($qw>$tw and $tw>0)?$tw:($qw>0)?$qw:$tw;  $isfrag= ($qutrbad==0 and $tw < $ALTFRAG*$qw)?"frag":"";  }
#p4o     ## if($qc < $MINCDS) { $class{$qd}.="tiny"; } elsif($qc < $MINCDSPART and $ispartial<needqual){ ..}
#p4o     
#p4o     $pid= ($aln<1)?0: int(0.5+ 100*$iden/$aln); $pid=100 if($pid>100);
#p4o     $pal= ($alnmax<1)?0 : int(0.5+ 100*$aln/$alnmax); $pal=100 if($pal>100);
#p4o     #o my $palq= ($qsize<1)?0 : int(0.5+ 100*$aln/$qsize); $palq=100 if($palq>100);
#p4o     #o my $palt= ($tsize<1)?0 : int(0.5+ 100*$aln/$tsize); $palt=100 if($palt>100);
#p4o     #^ add $palq = $aln/$qc; $palt= $aln/$tc; or palmax = aln/alnmin 
#p4o     
#p4o     ##WrongWayPeachy#my $tinyalndiff= (TEST3 && (($aln - $alnmax) < $NHIALN)) ? 1 : 0; # TEST3 add, use w/ $pal
#p4o     my $tinyalndiff= (TEST3 && (($alnmax - $aln) < $NHIALN)) ? 1 : 0; # TEST3 add, use w/ $pal
#p4o     #NOT this,see samesize# my $samealn= ($tinyalndiff && abs($qc - $tc) < $NHIALN) ? 1 : 0; # TEST3 add2; move up to IS_CDS ??BAD calc?
#p4o 
#p4o ## PROBLEM using eqgenes .. wont get here unless td x qd is in blastn.alntab ; some are NOT.
#p4o ## check all $eqgenes->{$td} ?? should have self match.
#p4o ##* Need to change this altmap class: NO cds-align b/n td,qd does mean something, poor mapping paralogs end up here
#p4o ## 13Sep01 : add eqgene map info for missed alignments : pal adjustment
#p4o ## 13sep06 : maybe not here .. trust align score if it exists for pair? at least needs flag if changes pal.
#p4o 
#p4o use constant TEST1603 => 1; ## test use here eqgenes/eqexons??
#p4o 
#p4o 		$havepair{$qd}{$td}= $pal; #?? change this to pal align val for altpar TEST1603?
#p4o 
#p4o =item altpar problem case
#p4o 
#p4o -- may need align distance tree/cluster for altpars, and decision w/ other data on where to cut tree to separate loci
#p4o -- part of problem is at alt findmain, iteration should stop when palign drops a lot
#p4o    it t1/t4 is wrong main for many of these (t2,t3 and children)
#p4o    palign drops from >60% to <20% for t2,t3,..
#p4o    
#p4o -- there are 2 loci at least, diff gmap, diff aaref, alts of 2nd remain linked to 1st tho.
#p4o egrep '^Anofunz4hEVm004829t.    ok' evg24m2banofun.tgclass3     
#p4o >>locus1
#p4o >> t1 should be main not t4, aafull, refbest, 100align, 
#p4o >> is t4 artifact? may be mashup of loc1,loc2 as it aligns well to both sets
#p4o Anofunz4hEVm004829t1	okay	althi	Anofunz4hEVm004829t4	99/68/./altmap65xeq/Anofunz4hEVm004829t5	509,95%,complete	aaref:972,AGAP002866-PA,refbest,chrmap:100a,100i,1530l,2x,KB669169:1413398-1415008:-,pflag:0
#p4o Anofunz4hEVm004829t4	okay	main	Anofunz4hEVm004829t5	100/98/./altmap58	554,96%,partial3	aaref:897,AGAP002866-PA,refgood,chrmap:65a,98i,1662l,2x,KB669169:1413434-1414585:-,pflag:0
#p4o Anofunz4hEVm004829t5	okay	althi1	Anofunz4hEVm004829t4	100/98/./altmap58	490,84%,complete	aaref:885,AGAP002866-PA,refgood,chrmap:100a,100i,1430l,3x,KB669169:1413398-1415008:-,pflag:0
#p4o 
#p4o >>locus2
#p4o Anofunz4hEVm004829t2	okay	main	Anofunz4hEVm004829t4	99/3/./paralt15	509,89%,complete	aaref:928,AGAP002865-PA,refbest,chrmap:100a,100i,1530l,2x,KB669169:1375626-1377211:-,pflag:0
#p4o Anofunz4hEVm004829t8	okay	althi	Anofunz4hEVm004829t4	99/3/./paralt17	433,90%,complete	aaref:779,AGAP002865-PA,chrmap:100a,98i,1302l,2x,KB669169:1375626-1376983:-,pflag:0
#p4o >> ? locus2b w/ crossmap ?
#p4o Anofunz4hEVm004829t7	okay	althi1	Anofunz4hEVm004829t4	99/3/./paralt18	433,92%,complete	aaref:778,AGAP002865-PA,chrmap:100a,98i,1302l,2x,KB669169:1375626-1412196:-,pflag:0
#p4o Anofunz4hEVm004829t9	okay	althi1	Anofunz4hEVm004829t4	99/3/./paralt18	433,90%,complete	aaref:772,AGAP002865-PA,chrmap:100a,99i,1302l,2x,KB669169:1375626-1412196:-,pflag:0
#p4o Anofunz4hEVm004829t6	okay	altmid	Anofunz4hEVm004829t4	99/3/./paralt18	433,79%,complete	aaref:775,AGAP002865-PA,chrmap:100a,98i,1302l,2x,KB669169:1410838-1412196:-,pflag:0
#p4o 
#p4o >> ? locus3 : belongs w/ Anofunz4hEVm004760t23 Anofunz4hEVm005570t1 Anofunz4hEVm005571t1 and others n=29 same aaref, all over KB668962:33384-35764:-
#p4o >> but some have low align, may be 4th? paralog loc?
#p4o Anofunz4hEVm004829t3	okay	althi	Anofunz4hEVm004829t4	100/3/./paralt16	508,83%,complete	aaref:820,AGAP012957-PA,chrmap:84a,99i,1527l,6x,KB668962:33634-35764:-,pflag:0
#p4o 
#p4o =cut
#p4o 
#p4o 
#p4o 		if($neqgene>0) { 
#p4o if(TEST1603) {
#p4o 			my $palmap= $eqgenes->{$qd}{$td}||0; # note this is pct of td aligned to qd
#p4o       if($palmap>0) { 
#p4o 			  my $xeq= $eqexons->{$qd}{$td}||0; # only for ($neqexon>0)
#p4o 			  my $XPHI = 95; my $XPLO= 3;
#p4o 			  # xeq care about (a) xeq>= identity == not alt but redundant, 
#p4o 			  #   (b) xeq <= noxover, not alt but paralog maybe, (c) middle = usual alt
#p4o         
#p4o         #* set aln=0/TINYALNBASES when pal=0/TINYALN
#p4o         #*? change 'paralt' to 'altpar' to avoid other evg parse problems?
#p4o         if($neqexon>0 and $xeq >= $XPHI) { 
#p4o           if($palmap < $PHI) { } #?? no change or yes
#p4o           $antiflag .="/altmap$palmap"."xeq"; 
#p4o           $eqflag{$td}{$qd}="altmapxe$palmap.$xeq"; # for report, which way?
#p4o           $pal=$palmap if($palmap >= $PHI and $pal < $PMID);  # FIXME bad annot
#p4o           }
#p4o         elsif($neqexon>0 and $xeq < $XPLO) { # problems here.. dont reset pal yet
#p4o           $antiflag .="/altparx$pal";  #later? $pal=3; $aln=13;
#p4o           $eqflag{$td}{$qd}="altparx$palmap.$xeq"; # for report
#p4o           }
#p4o         else { 
#p4o           ## skip this case for tid =~ qid already alts of same locus, dont really need now?
#p4o           (my $qg=$qd)=~s/t\d+$//;  unless($td=~/^$qg/) {
#p4o           $antiflag .="/altmap$palmap";
#p4o           $eqflag{$td}{$qd}="altmap$palmap.$xeq"; # for report
#p4o           $pal=$palmap if($palmap >= $PHI and $pal < $PMID);  # FIXME bad annot
#p4o           }
#p4o           }
#p4o         }
#p4o       elsif($pal>0 and exists( $eqgenes->{$td}) and exists( $eqgenes->{$qd})) { # defined or exists??
#p4o         ## bad here for nomap alts ** why isnt defined/exists working? for ids not in eqgene
#p4o         if($gmapqual->{$qd} and $gmapqual->{$td}) {
#p4o           my $revmap= $eqgenes->{$td}{$qd}||0; 
#p4o           ## unmapped gene cases are problem here.. need what? defined ($eqgenes->{td}) 
#p4o           $antiflag .="/altpar$pal" if($revmap<3);  # $pal=3; $aln=13; #? change later?
#p4o           $eqflag{$td}{$qd}="altpar$pal.$palmap.$revmap"; # for report
#p4o         } else {
#p4o           #lots# warn "#DBG: bad eqgene $qd,$td\n" if($debug);
#p4o         }
#p4o         }
#p4o } else {			
#p4o 			my $palmap= $eqgenes->{$qd}{$td}||0; 
#p4o 			if($palmap>0) { $antiflag .="/altmap$palmap"; 
#p4o   			if(EQGENE_OVERRIDES_ALN && $palmap > $pal) { $pal=$palmap; } ##? not this change
#p4o 			}
#p4o }			
#p4o 		}
#p4o 
#p4o     my $pidalnval="$pid/$pal$antiflag"; # old: $pid/$pal
#p4o     my $tisaltpar=0; # treat like $skiptinyaln for now
#p4o     if(TEST1603) {
#p4o       $tisaltpar=1 if($pidalnval=~/altpar\d/); #?  and $qclass // not altparx\d
#p4o     }
#p4o     if(TEST1603) { # this may be sort of working right now..
#p4o     $bestmatch{$qd}="$td,$pidalnval" unless($bestmatch{$qd} and not($bestmatch{$qd}=~/altpar\d/));
#p4o     $bestmatch{$td}="$qd,$pidalnval" unless($bestmatch{$td} and not($bestmatch{$td}=~/altpar\d/)); 
#p4o      # prob here for altpar ?? need to replace td>qd
#p4o     } else {    
#p4o     ## .. change these to use palq, palt : palign per q,t size ?
#p4o     $bestmatch{$qd}="$td,$pidalnval" unless($bestmatch{$qd});
#p4o     $bestmatch{$td}="$qd,$pidalnval" unless($bestmatch{$td});  
#p4o         #^^ maybe should replace best{td}= new qd if new qd.pid >> oldqd.pid ?
#p4o     }
#p4o     #new# unless($bestmatch{$qd}) { $bestmatch{$qd}="$qd,100/100/self1"; } # always need bestmatch{qd} entry
#p4o         
#p4o # Maybe problem here making too many fragment alts: a,b near identical alts of c, b <= a subset,
#p4o # but both classed/bestmatch to c main.  b should be dropped as fragment/99equiv of a, but dont see that due
#p4o # to b<c precedence.
#p4o 
#p4o     unless(TEST3) {
#p4o     if($class{$td}) { next; }  # TEST3: yes, defer this to afer aclass, change for althi1 
#p4o     }    
#p4o     
#p4o     my $qclass= $class{$qd}||"";
#p4o     if($samesize and $qclass =~ /alt/ and $ismain{$td}) { # check/stop alt1 <> alt2 assigns
#p4o       next if($bestmatch{$qd} =~ /$td,/); # problem here? for many equal bestmatch
#p4o     }
#p4o 
#p4o     #..........    
#p4o       ## reclassify? althi > althi100 for 2ndary alts; althi100 = unneeded duplicates : TEST
#p4o       ##if($class{$td} eq "althi" and $pid>99 and $pal>=99) { $class{$td}= "althi100b"; }
#p4o       # ^ this probably wont happen, 2ndary althi, 1st will also be 100 pi/pa
#p4o       # ** problem using cdhits clusters here, it doesnt give 2ndary aligns as blastn does
#p4o       # .. so partial align to big tr may also be perfect align to smaller, valid alt tr
#p4o       ## next;   # problem above for qd1,qd2 > td, need to keep qd2 with bestmatch{qd}
#p4o     
#p4o     # UPD 1602: change TINYALN use for MINALIGN bases
#p4o     my $skiptinyaln=0;
#p4o     if(UseTINYALNBASES) {
#p4o       if($alnmax >= $MINCDS) { # recall alnmax is min(qsize,tsize); never skip tiny prots, MINCDS =~ 90 bases ?
#p4o         my $minbaseover= ($alnmax >= $TINYALNCDSLEN)? $TINYALNBASES : $alnmax * 0.10; ##$MIN_PCTOVERLAP;
#p4o         $skiptinyaln=($aln < $minbaseover)?1:0;
#p4o       } else {
#p4o         $isfrag= "frag" unless($isfrag); #or "tiny" ? # dont skip assign alt/frag to tiny cds, but maybe always set isfrag ?
#p4o       }
#p4o     } else { # old, pct-TINYALN is a problem for large cds  
#p4o       $skiptinyaln= ($pal < $TINYALN)?1:0;
#p4o     }
#p4o   
#p4o     # set alt aclass here, or not.
#p4o     if( $skiptinyaln or $tisaltpar) { } ## defer: $aclass="noalign$isfrag";  
#p4o     elsif( $tc == 0 and $qc > 0 and $pid >= $PLOW ) { $aclass="frag0aa"; } # tc==0 special case "frag0"
#p4o     elsif( $pid >= $PHI ) { 
#p4o       $aclass= ($isfrag)?"parthi":"althi"; # Primary classifier; change parthi to althipart ? althi$isfrag ?
#p4o       $aclass .= "1" if($pal >= $PHIALN or $tinyalndiff); # althi1 == hi-align + hi-ident, nearly same (identicals gone here)
#p4o       }
#p4o     elsif( $pid >= $PMID ) { $aclass="altmid$isfrag"; }
#p4o     elsif( $pid >= $PLOW ) { $aclass="altlo$isfrag"; }
#p4o       ## old $pid >= $PHI       
#p4o       # old0:  $aclass= ($pal<$PHIALN or $isfrag)?"parthi":"althi";
#p4o     	# TEST3: add min-basediff here to PHIALN tests, so shorty tr dont take over w/ minor base diff
#p4o     	# NO, see samesize: add $samealn flag, althi0 or althi2 ?? for ~identical size and ~perfect align
#p4o       # if($pid >= 99) { ## PHI == 99 default
#p4o       #   #old2# if($samealn) {$aclass .= "2"; } elsif # bad calls here, pal lowish 80%
#p4o       #   #old3.use this# if($pal >= $PHIALN or $tinyalndiff) { $aclass .= "1"; }  # "100"; use or not? PHIALN reused: was 99 here
#p4o       #   #old1 # $aclass .= "1" if($pid > 99 and $pal >= 99);  # "100"; use or not? PHIALN reused: was 99 here
#p4o       #   }
#p4o 
#p4o 
#p4o =item samesize/samebestmatch 
#p4o   problem of 2+mains w/ several essentially identicals, 
#p4o   after 2 mains created, then row linking them needs to break tie
#p4o 
#p4o egrep '^Anofunz4hEVm00398[2345]t1       ' $pt.tgclass | egrep -v 't[23456789]|t..       ' 
#p4o Anofunz4hEVm003982t1	okay	main	Anofunz4hEVm003985t1	100/100/./altmap98	629,94%,complete	aaref:1108,AGAP001965-PA,refbest,chrmap:98a,100i,1890l,5x,KB669169:464434-466658:+,pflag:0
#p4o Anofunz4hEVm003983t1	okay	main	Anofunz4hEVm003984t1	100/100/./altmap98	629,94%,complete	aaref:1101,AGAP001965-PA,refgood,chrmap:98a,99i,1890l,5x,KB669169:464434-466658:+,pflag:0
#p4o Anofunz4hEVm003984t1	okay	althi1	Anofunz4hEVm003983t1	100/100/./altmap98	629,94%,complete	aaref:1100,AGAP001965-PA,refgood,chrmap:98a,98i,1890l,5x,KB669169:464434-466658:+,pflag:0
#p4o Anofunz4hEVm003985t1	okay	althi1	Anofunz4hEVm003982t1	100/100/./altmap98	629,94%,complete	aaref:1108,AGAP001965-PA,refbest,chrmap:98a,100i,1890l,5x,KB669169:464434-466658:+,pflag:0
#p4o 
#p4o egrep '^Anofunz4hEVm00398[2345]t1       ' evg24mergeanofun.tgalntab | egrep -v 't[23456789]|t.. ' | sort -k2,2nr -k7,7nr -k6,6n
#p4o Anofunz4hEVm003982t1	1890	1996	Anofunz4hEVm003985t1	1890	1996	1888	1883	3602  # main1
#p4o Anofunz4hEVm003985t1	1890	1996	Anofunz4hEVm003982t1	1890	1996	1888	1883	3602  # altof main1
#p4o Anofunz4hEVm003983t1	1890	1996	Anofunz4hEVm003984t1	1890	1996	1886	1885	3621  # main2
#p4o Anofunz4hEVm003984t1	1890	1996	Anofunz4hEVm003983t1	1890	1996	1886	1885	3621  # altof main2
#p4o Anofunz4hEVm003983t1	1890	1996	Anofunz4hEVm003985t1	1890	1996	1841	1840	3535  # main2 = altof main1
#p4o Anofunz4hEVm003985t1	1890	1996	Anofunz4hEVm003983t1	1890	1996	1841	1840	3535  # altof main1 = main2
#p4o Anofunz4hEVm003982t1	1890	1996	Anofunz4hEVm003983t1	1890	1996	1811	1811	3483  # main1 = main2, break tie
#p4o Anofunz4hEVm003983t1	1890	1996	Anofunz4hEVm003982t1	1890	1996	1811	1811	3483  .. more, dont redo break
#p4o Anofunz4hEVm003982t1	1890	1996	Anofunz4hEVm003984t1	1890	1996	1800	1800	3462
#p4o Anofunz4hEVm003984t1	1890	1996	Anofunz4hEVm003982t1	1890	1996	1800	1800	3462
#p4o Anofunz4hEVm003984t1	1890	1996	Anofunz4hEVm003985t1	1890	1996	1800	1799	3456
#p4o Anofunz4hEVm003985t1	1890	1996	Anofunz4hEVm003984t1	1890	1996	1800	1799	3456
#p4o 
#p4o =cut
#p4o 
#p4o   if(TEST3) {
#p4o     if(my $ocl= $class{$td}) {  
#p4o     
#p4o   if(TEST1602) {  
#p4o       # test2: next is bad here, for class{td} = class{qd} == main both
#p4o       #* this helps, reduces main: 10632/17263, +noclass: 3964/4141, but adds those to okay.alt: 90048/60193
#p4o       #* need to fiddle w/ alt drops now
#p4o       ## maybe should not NEXT here, but set ismain{qd},class{qdmain} ??
#p4o       
#p4o       #** this may be LOSING valid mains, set in prior row, now its alt is resetting main to althi1 ..
#p4o       if($ocl =~ m/^alt|^part/) { $class{$td}= $aclass if($aclass eq "althi1"); next; } # test3
#p4o       elsif($ocl =~ m/^main/) { 
#p4o         # **  problem of 2+mains w/ several essentially identicals,  break tie here **
#p4o         if($qclass =~ m/^main/ or $tisaltpar) { 
#p4o           # skip to below aclass set, test5, yes, do this also?
#p4o           } 
#p4o         else { next; } # * YES, need this : test4
#p4o       } 
#p4o   } else {
#p4o       $class{$td}= $aclass if($ocl =~ m/^alt/ and ($aclass eq "althi1")); # is this right?  or  $aclass eq "althi2"??
#p4o       next;   
#p4o   } 
#p4o     }
#p4o   }    
#p4o 
#p4o     ## 1602 update: $ismain{$qd} and $ismain{$td} already called, later row sez they are equal also, reclass 1.
#p4o     ## .. but below changes 2nd main to alt: $class{$td}=$aclass
#p4o 
#p4o     if($aclass) { 
#p4o       ## FIXME ismain{qd} >> qd may be alt of other main, .. follow chain here?  keep qd as bestmatch? both?
#p4o       my $qmain=$qd; my $attr= $pidalnval; ## "$pid/$pal$antiflag";
#p4o 
#p4o       ## FIXME: this still leaves some alts w/o final main id link; including circular alt1 <=> alt2 links
#p4o       ## fix for (TEST1603) altpars, need to stop findmain iter when palign drops
#p4o       # my $isaltpar=$tisaltpar; if(TEST1603) { $isaltpar= $tisaltpar;}  ## == ($pidalnval=~/altpar/)?1:0; # 
#p4o       
#p4o       use constant FINDMAINID => 1;  
#p4o       # this is right now: if(FINDMAINID) ..
#p4o       if($class{$qd}) {  #  =~ /alt|part/
#p4o         my $qnext= $qd; 
#p4o         my $more= ($bestmatch{$qnext})?1:0;
#p4o         my %qdid=( $qd => 1, $td => 1);
#p4o         while($more) { 
#p4o           $more=0;
#p4o           my($qbest,$qpal)=split",",$bestmatch{$qnext};
#p4o           
#p4o           if(TEST1603) {
#p4o             $qbest=0 if($ismain{$qnext} and not $ismain{$qbest}); # want always? maybe bad?
#p4o  		        #? my $qtpal= $havepair{$qnext}{$td} || $havepair{$td}{$qnext}; #?? change this to pal align val for altpar TEST1603?
#p4o  		        #? $qbest=0 if($tisaltpar and defined $qtpal and $qtpal < $PLOW); # not good?
#p4o           }
#p4o           
#p4o           if($qbest and not $qdid{$qbest}) {
#p4o             $qnext= $qbest; $qdid{$qbest}=1; 
#p4o             $more=1 if($class{$qnext} and $bestmatch{$qnext});
#p4o           }
#p4o         }
#p4o         if($qnext ne $qd) { $qmain= $qnext; $attr="$pidalnval/$qmain"; } #was $attr.="/$qmain";
#p4o         ## WARNING: now attr has /-sense or /. before /qmain; 
#p4o         ##   evgmrna2tsa2.pl needs to know this field's structure, /qmain at end esp.
#p4o       }
#p4o       
#p4o       $ismain{$qmain}++;   
#p4o       $class{$td}=$aclass; #? problem here when class{td} == main set before, should not reset to alt
#p4o if(TEST1602) {   
#p4o       if(1) { # this way works.  if($bestmatch{$td}) must be true from above
#p4o       if($bestmatch{$td} =~ /^$qd/) { $bestmatch{$td}="$qd,$attr"; }
#p4o       elsif(not $bestmatch{$td}) { $bestmatch{$td}="$qd,$attr"; }
#p4o       #? elsif($qmain ne $qd and $bestmatch{$td} !~ m,/,) { $bestmatch{$td} .= "/$qmain"; }
#p4o       }    
#p4o       # if(0) { # test16: this is BAD, bestmatch{td} always set here, usually to qd w/o updated attr
#p4o       #   $bestmatch{$td}="$qd,$attr" unless($bestmatch{$td}); #?keep 1st bestmatch? NOT this test, always have bmatch
#p4o       # }
#p4o       #^ this change if already bestmatch has effect, is it adding /qmain above, or is it not qd ?  test this:
#p4o       # if($ismain{$td}) { 
#p4o       #  #? delete $ismain{$td}; #? yes or no? prevent both aligned tr as main, but class{td}=alt does that
#p4o       # }
#p4o } else {      
#p4o       $bestmatch{$td}="$qd,$attr"; # orig way 
#p4o }
#p4o if(TEST3) {
#p4o       ## add to prevent 2+ circular alta <> altb <> altc with no main ..
#p4o       $class{$qmain}="main" unless($class{$qmain}); # should always be:  $qc >= $tc//$samesize..
#p4o }
#p4o     }
#p4o       
#p4o   } close($inh);
#p4o 
#p4o   
#p4o   # END:  print trclass table; add more fields to output: aaqual, aablast, tr,aa sizes?
#p4o   # FIXME: here, elsewhere create ID main,alt table, with new numeric IDs, old/cur ids, main/alt num
#p4o   { my($q,$pal,$c,$d);
#p4o   foreach $d (sort keys %class) {
#p4o     ($q,$pal)=split",",$bestmatch{$d}; $c= $class{$d}; 
#p4o     my @cla= classifytr( $d, $c, $q, $pal);
#p4o     print $outh join("\t",@cla)."\n"; 
#p4o     }
#p4o   foreach $d (sort grep{ not($class{$_}) } keys %ismain) {
#p4o     ($q,$pal)=split",",$bestmatch{$d}; $c= "main";    
#p4o     my @cla= classifytr( $d, $c, $q, $pal);
#p4o     print $outh join("\t",@cla)."\n"; 
#p4o     }
#p4o   foreach $d (sort grep{ not($class{$_} or $ismain{$_}) } keys %bestmatch) {
#p4o     ($q,$pal)=split",",$bestmatch{$d}; $c= "noclass";  
#p4o     my @cla= classifytr( $d, $c, $q, $pal);
#p4o     print $outh join("\t",@cla)."\n"; 
#p4o     }
#p4o   }  
#p4o   
#p4o }


# pre-read all of infile ids into validids, before idenityclass
sub readIdsFromAlnTab {
  my($infile, $insorted)= @_;
  my($nids,$inh)=(0,undef);
  # ($ok,$inh)= openRead($infile,1);
  open(IN,$infile) or die "reading $infile"; $inh=*IN; 
  while(<$inh>) {  
    next if(/^Qid|^\W/); chomp; 
    my($qd,$qc,$qw,$td,$tc,$tw,$aln,$iden,$bits)= split"\t"; 
    
    #?? use dupids here to mark as not validids ?
    if( $dupids ) {
      $validids{$qd}++ unless( $dupids{$qd} and $qd ne $dupids{$qd});    
      $validids{$td}++ unless( $dupids{$td} and $td ne $dupids{$td});    
    } else {
      $validids{$qd}++; $validids{$td}++;
    }
    $nids++;
  } close($inh); # rewind if ref() ???
  return $nids;
}


sub readblasttab {
  my ($bother)= @_;
  
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq,$nids)= (0) x 10;
  my($ok,$fh)= openRead($bother);
  # my $fh;
  # if($bother =~ /\.gz/) { open(F,"gunzip -c $bother |") or die $bother; $fh=*F; }
  # elsif($bother =~ /stdin|^-/) { $fh= *STDIN; }
  # else { open(F,$bother) or die $bother; $fh=*F; }
  
  %bspans=();
  my %dupskipids=(); my $dupskipspan=0;
  
  ## FIXMEd: self-only matches need  recording via putspans; should be but dupdrop is droping firstid also

  while(<$fh>) { 
    unless(/^\w/) { next; } 
    #  if(/^# Query/) { ## Unused info
    #  #Notnow# puts($lq,$lt,$sa,$sm) if($lt); 
    #  ($qd)=m/Query:\s+(\S+)/; $wq=(m/len=(\d+)/)?$1:0; }
      
    my @v= split; 
    my($q,$t,$bits,$pctident,$aln,$mis,$gap,@bspan)= @v[0,1,-1,2,3,4,5, 6,7,8,9]; # 6-9 =  q. start, q. end, s. start, s. end, 

    # FIXME: need to parse align parts before summing;
    # only(SPANSUM)  now
    if($lq and $q ne $lq) {  
      putspans($lq) unless($dupskipspan); $nids++; %bspans=(); $dupskipspan=0; 
    }

  ## FIXME: self-only matches need  recording via putspans
    my $dupskip=0;
if(DUPFILTER1) {    ## DF1 active; see also inactive DUPFILTER2
    # checkme: is dupids mainid always in blast set? if not dupskip drops useful items
    # ** DUPID list has ids not in blast/cluster set...
    #   add dupskipids{id}=mainid and check later in validids
    if( $dupids ) { 
      if( $dupids{$q} and $q ne $dupids{$q} ) { 
        $dupskipspan=$q; ## if($q eq $t); # flag to skip putspan
        $dupskipids{$q}++; $dupskip++;   
      } elsif( $q ne $t and $dupids{$t} and $t ne $dupids{$t}) { 
        $dupskipids{$t}++; $dupskip++; 
      }  
      if($dupskip) { $ndupdrop++; }  # next; below, not here  WRONG now, need to skip self putspan unless this dup is main 
    }
}   
    ## UNUSED now: qd, wq : FIXME save ids of no-match
    if($t eq $q) { } ## { $qd=$q unless($qd); $wq=$aln unless($wq); }
    elsif($dupskip==0 and $dupskipspan == 0) { 
 
     $bits= bint($bits);
      # only(SPANSUM)  
      my $aident= _max(0,$aln-$mis-$gap); # other way to calc: $aident = $pctident * $aln;
      #new: my $aident= int(0.5 + $aln*$pctident/100); # better?
      sumblastpart( $q, $t, $bits,$aln,$aident, @bspan); # if($bits > $pMINLOW * $bmax or $bspans{$t}); 
    } 
    $lq=$q; $lt=$t;  
  } close($fh);
  
  #only(SPANSUM)
  putspans($lq) unless($dupskipspan); $nids++; $dupskipspan=0;

  if($ndupdrop>0) { $ndupdrop=scalar(keys %dupskipids); }
  warn "# readblasttab: nids=$nids; ndupdrop=$ndupdrop\n" if $debug;
  return($nids);
}

## 2011.aug BUG here, need to test sb-se outside tb-te spans also
## 2013.aug : IS_CDSALIGN ORIENT or is IMPORTANT : need to know when alt-cds are reversed
sub sumblastpart {
  my( $q, $t, $bits,$aln,$aident, $qb,$qe,$sb,$se) = @_;
  my $or=0; # 2013.aug : ORIENT problem.  if both here are reversed, or=0; if only 1, or=-1
  if($qb > $qe) { ($qb,$qe)= ($qe,$qb); $or=-1; }
  if($sb > $se) { ($sb,$se)= ($se,$sb); $or=($or<0)?0:-1; } #was $or--
  unless($bspans{$t}) { 
    $bspans{$t}=[]; push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or]); 
    return; }
  my $ov=0;
  ## 2011oct overslop pct fix
  my $qlen=1+$qe-$qb; my $slen=1+$se-$sb;
  my $qslop= _max($OVSLOP, int($pctOVSLOP*$qlen));
  my $sslop= _max($OVSLOP, int($pctOVSLOP*$slen));
  
  foreach my $sp (@{$bspans{$t}}) {
    my($xb,$xe,$tb,$te,$xbit)= @$sp;
    if($qe < $xb or $qb > $xe) { }
    elsif($qe > $xe and $qb >= $xe - $qslop) { }
    elsif($qb < $xb and $qe <= $xb + $qslop) { }
    else { $ov=1; last; }
    ## add 2011.aug
    if($se < $tb or $sb > $te) { }
    elsif($se > $te and $sb >= $te - $sslop) { }
    elsif($sb < $tb and $se <= $tb + $sslop) { }
    else { $ov=1; last; }
  }  
  ## IS_CDSALIGN add $or to bspans, problems?
  unless($ov) { push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or]); }
}



# change mainaa/subaa in %aacluster, using other input info: aaqual/size-gaps, cds/tr align ids
sub aaqualscore  # in cdna_evigenesub.pm
{
  my($mqual)= @_;
  my $mqv=0; $mqual ||="missing";
  if($mqual =~ /complete/) { $mqv += 2; } elsif($mqual =~ /partial[35]/) { $mqv += 1; }
  if($mqual =~ /utrbad/) { $mqv -= 2; } elsif($mqual =~ /utrpoor/) { $mqv -= 1; }
  if($mqual =~ /gapbad/) { $mqv -= 1; } # or -2?
  return $mqv; # range is now -3..+2, was -2 .. +2
}


sub correctAAcluster
{
  my($havevalidids)= @_;
  # % aaclustermain == hash{mainid}[subids]
  # % aacluster == hash{eachid} = "mainid,pctident"
  my $nreset=0; my %didmain;
  foreach my $mid (sort keys %aaclustermain) {
    next if($didmain{$mid}); # probably dont need. 
    my $mainid= $mid;
    my @cids= @ { $aaclustermain{$mid} }; # has all ids incl mainid?
    my $maw = $aasize{$mid} || 0;
    my $mqv = aaqualscore($aaqual{$mid});  #? parse for "aasize,pcds,aaqual"
    ## sort cids by aasize? or need to check thru all?
    # @cids = sort{$aasize{$b} <=> $aasize{$a}} @cids;
    my $reset=0;  my @goodids=(); my $ninval=0; 
    if($havevalidids and not $validids{$mainid}) { $maw=0; $mqv=-9; $mainid=0;  $ninval++; }
    foreach my $id (@cids) {
      ## also check each id is valid for tr/cds align, drop invalids including current mainid
      if($havevalidids and not $validids{$id}) { 
        ## if($id eq $mainid) { $maw=0; $mqv=-9; $mainid=0; } # above now: main could be last id.. fixme
        $ninval++; next; # drop from aacluster
      }
      push @goodids, $id; 
      next if($id eq $mainid);
      my $aw= $aasize{$id} || 0; 
      my $qv= aaqualscore($aaqual{$id});
      if($aw > $maw and $qv >= $mqv - 1) {  # reset main; ignore qv if aw >>> maw ?
        ($mainid,$maw,$mqv)=($id,$aw,$qv); $reset++; 
      } elsif($qv > $mqv and ($aw > 0.98*$maw)) {
        ($mainid,$maw,$mqv)=($id,$aw,$qv); $reset++;       
      }
    }
    
    if( @goodids == 0 and $ninval >= @cids) {
      foreach my $id (@cids) { delete $aacluster{$id}; }
      $aaclustermain{$mid}= [];
      $nreset++; $reset=0; 
    }
    if($reset) {
      foreach my $id (@goodids) {
        my($oldmain,$pi)=split",",$aacluster{$id};
        $aacluster{$id}= "$mainid,$pi";
        }
      $aaclustermain{$mainid}= \@goodids;
      $nreset++;
    }
    $didmain{$mainid}++;
  }
  warn "# correctAAcluster nreset=$nreset\n" if($debug); # got nreset=212165 for nclust=232463 TOO HIGH?
  return($nreset);
}


sub readAAcdhit {
  my ($bother)= @_;
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq)= (0) x 10;
  my($ok,$fh)= openRead($bother);

  %aacluster=(); # global
  our @cluster=();
  my $iscdhit=0;
  my($cli,$ncl,$nalt,$mainid,$mainlen,$nerr)=(0)x10;

  # FIXME: use aaqual utrbad/poor for main/alt with pi=100; reset main if utrbad and alt utrok
  sub putclus { 
    my($mainid,$cluster1)=@_; ## our @cluster;
    if($mainid) {  
      $aaclustermain{$mainid}=[] unless($aaclustermain{$mainid});
      foreach my $cl (@$cluster1) { 
      my($id,$pi)=split",",$cl; $aacluster{$id}="$mainid,$pi"; 
      push @{$aaclustermain{$mainid}}, $id;
      } 
    }
  }
  
  while(<$fh>) { 
    if(/^>/) { 
      putclus($mainid,\@cluster) if(@cluster); @cluster=(); $mainid=0;
      ($cli)= m/(\d+)/;  $ncl++; #   m/Cluster\s*(\d+)/;  fixed or not?
      
    # elsif(/^(\d\w*)\s+(\d+)(..), >(.+)\.\.\. (.+)$/)  ## BADD patt ??
    } elsif(/^(\d+)/) {
      my $i=$1;
      m/(\d+)(nt|aa), >(.+)\.\.\. (.+)$/;  
      my($tlen,$typ,$tid,$pinfo)=($1,$2,$3,$4); 
      ## FIXME: merge.clstr:  1.f1  9999aa, >id... at 100  << .f1,.f2 added // DROP this?
      ## FIXME2: new merge fnum at end of line now;
      my @pmore; ($pinfo,@pmore)= split /\t/, $pinfo;
      #cdhit perls: /(aa|nt), >(.+)\.\.\./
      my $pi;
      my $ismain=($pinfo =~ /\*/)?1:0;  # FIXME: drop /$/; new merge fnum at end of line now;
      # FIXME: merge.aa.clster bad; lacks some main * ; skip those for now
      unless($tid) {
        $nerr++;
      } elsif($ismain) {
        $mainid= $tid; $mainlen=$tlen; $pi=100;
        push @cluster, "$tid,$pi";
      } else {
        # pi for cdhit-est: at 1:492:892:1383/+/100.00%
        # pi for cdhit-aa : at 99.76%   <<<<<
        $pinfo =~ s/at //;  $pinfo =~ s/^\D+//; #? always \digit start
        ($pi)= $pinfo =~ m/([\d\.]+)%/; $pi=~s/\.00//; $pi||=0;
        push @cluster, "$tid,$pi";  $nalt++;
      }
    } elsif(/^\w/) {
      $nerr++; # warn/die not cdhit format ..
    }
  }
  
  putclus($mainid,\@cluster) if(@cluster); ### putspans($mainid) if($mainid); %bspans=();
  warn "# readAAcdhit: nclust=$ncl, nalt=$nalt, nerr=$nerr \n" if $debug;
  return($ncl); 
}



__END__


=item readlastz: lastz align general format

test case:
  aabugs4/tsaevgc/daphmag5xbest5/dmag5xau13c2011f
  pt=dmag5xau13c2011

$evigene/scripts/rnaseq/asmrna_dupfilter2.pl -debug -CDSALIGN -tinyaln 35 \
  -aasize inputset/$pt.aa.qual -acdhit tmpfiles/${pt}_cd90.aa.clstr \
  -lastz tmpfiles/$pt-self97.lastz.gz  \
  -outeq tmpfiles/$pt.alnlztab2 -outclass $pt.trclasslz2

$bg/mb/galn/bin/lastz \
 --identity=100 --coverage=25 --step=10 --seed=match15 --notransition --exact=20 --match=1,5 \
 --ambiguous=n --nochain --nogapped  --format=general 'altset1.okboth.cds[multiple]' altset1.okboth.cds

#score  name1              strand1 size1 zstart1 end1  name2               strand2 size2 zstart2 end2 identity idPct  coverage covPct

779     dmvel4xpak25Loc11378t4   +  1083    0    779   dmvel4xpak25Loc11378t10  +   1083  0    779  779/779 100.0%  779/1083  71.9%
1083    dmvel4xpak25Loc11378t10  +  1083    0    1083  dmvel4xpak25Loc11378t10  +   1083  0    1083 1083/1083 100.0%  1083/1083 100.0%
770     dmvel4xpak25Loc11378t5   +  978     208  978   dmvel4xpak25Loc11378t10  +   1083  313  1083 770/770 100.0%  770/978 78.7%
251     dmvel4xpak25Loc11378t1   +  684     0    251   dmvel4xpak25Loc11378t10  +   1083  399  650  251/251 100.0%  251/684 36.7%
330     dmvel4xpak25Loc11378t1   +  684     354  684   dmvel4xpak25Loc11378t10  +   1083  753  1083 330/330 100.0%  330/684 48.2%
324     dmvel4xpak25Loc11378t7   +  753     429  753   dmvel4xpak25Loc11378t10  +   1083  753  1077 324/324 100.0%  324/753 43.0%
    ...
455     dmvel4xpak25Loc11378t3   +  753     0    455   dmvel4xpak25Loc11378t7   +   753   0    455  455/455 100.0%  455/753 60.4%
753     dmvel4xpak25Loc11378t7   +  753     0    753   dmvel4xpak25Loc11378t7   +   753   0    753  753/753 100.0%  753/753 100.0%
456     dmvel4xpak25Loc11378t1   +  684     222  678   dmvel4xpak25Loc11378t7   +   753   297  753  456/456 100.0%  456/684 66.7%
324     dmvel4xpak25Loc11378t10  +  1083    753  1077  dmvel4xpak25Loc11378t7   +   753   429  753  324/324 100.0%  324/753 43.0%
324     dmvel4xpak25Loc11378t5   +  978     648  972   dmvel4xpak25Loc11378t7   +   753   429  753  324/324 100.0%  324/753 43.0%

=cut

sub readlastz {
  my ($bother)= @_;
  
  my $islzg=0;
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq, $nids)= (0) x 10;
  my($ok,$fh)= openRead($bother);
  # my $fh;
  # if($bother =~ /\.gz/) { open( $fh,"gunzip -c $bother |") or die $bother; }
  # elsif($bother =~ /stdin|^-/) { $fh= *STDIN; }
  # else { open($fh,$bother) or die $bother;  }
  %bspans=();

  ## TEST: is lastz output file sorted properly? LOOKS OK //name2 should all be grouped, does split-run undo that?
  
  my @hd; # $islzg=1 ;
  while(<$fh>) { 
    if(!$islzg and m/^#score\tname1/ and /\tcoverage/) { $islzg=1; chomp; s/^#//; @hd=split"\t"; } # should do but want some slack here?
    next unless(/^\d/);     
    chomp; my @v= split"\t";
    unless(@v==15 and $islzg) { die "# ERR: doesnt look like lastz general format I know: hd=",@hd," val=",@v,"\n" ; }
    
    ## allow subset columns?
    #score  name1   strand1 size1 zstart1 end1 name2   strand2 size2   zstart2 end2  identity idPct   coverage   covPct
    my( $lzscore, ## is lz score = count of ident bases? no,
      $tid, $tor, $tsize, $tstart, $tend,
      $qid, $qor, $qsize, $qstart, $qend,
      $nida, $pid, $ncovb, $pcov)= @v;
    $qstart++; $tstart++; # move to 1-origin

    ## NOTE: my lastz out has name2 as first-order == qid, name1 == tid
    if($lq and $qid ne $lq) {  
      putspans($lq);  $nids++; %bspans=();
    }

    if($tid eq $qid) { 
      #?? add:  $trsize{$qid} = $qsize unless($trsize{$qid});
      # $qd=$qid unless($qd); $wq=$aln unless($wq); 
    } else { 
      my($aident,$na)= split "/",$nida;
      my($aln,$nb)= split "/",$ncovb;
      my $bits=  $lzscore; # what? maybe 2*aident? 
      
# if(1) { ## SPANSUM
      ## *?* Need this; lastz hsp as for blastn can overlap lots ..
      sumblastpart( $qid, $tid, $bits,$aln,$aident, $qstart,$qend,$tstart,$tend); 
# } else {      
#       $bspans{$tid}=[]  unless(defined $bspans{$tid}); 
#       push( @{$bspans{$tid}}, [$qstart,$qend,$tstart,$tend,$bits,$aln,$aident]);  
# }
    } 
    ($lq,$lt)= ($qid,$tid);
  } close($fh);
  
  putspans($lq);  $nids++;   
  warn "# readlastz: nids=$nids\n" if $debug;  # ; ndupdrop=$ndupdrop
  return($nids);
}



=item cd-hit(est) align cluster format

  has align-span, percent ident, enough for use.
  
  cacao11pub3ig_cde25.cds.clstr
  >Cluster 0
  0       16221nt, >Thecc1EG029098t1... *
  >Cluster 1
  0       15408nt, >Thecc1EG019010t1... at 1:15408:1:15411/+/99.98%
  1       15411nt, >Thecc1EG019010t2... *
  >Cluster 2
  0       14343nt, >Thecc1EG007738t1... *
  1       10578nt, >Thecc1EG007738t2... at 1:10574:1831:12404/+/100.00%
  >Cluster 4
  0       12732nt, >Thecc1EG015810t1... at 1:12732:304:13035/+/100.00%
  1       13035nt, >Thecc1EG015810t2... *
  2       12504nt, >Thecc1EG015810t3... at 74:12503:148:12584/+/99.94%
  3       12717nt, >Thecc1EG015810t4... at 1:12717:304:13035/+/99.88%
  >Cluster 15
  0       9804nt, >Thecc1EG034527t1... *
  1       9804nt, >Thecc1EG034527t2... at 1:9804:1:9804/+/100.00%
  2       9804nt, >Thecc1EG034527t3... at 1:9804:1:9804/+/100.00%
  3       7512nt, >Thecc1EG034527t4... at 1:7512:2287:9804/+/99.92%
  >Cluster 21
  0       8751nt, >Thecc1EG006991t1... *
  1       8304nt, >Thecc1EG006991t2... at 2894:8304:3336:8751/+/99.83%
  2       8178nt, >Thecc1EG006991t3... at 2894:8178:3336:8630/+/99.74%
  3       8391nt, >Thecc1EG006991t4... at 3239:8377:3356:8494/+/100.00%

  dmag4vel4xfi_cde60.cds.clstr
  >Cluster 8633
  0       1383nt, >dmag4vel4xfik65Loc51t4... *
  1       1383nt, >dmag4vel4xfik75Loc49t1... at 1:1383:1:1383/+/100.00%
  >Cluster 8634
  0       1383nt, >dmag4vel4xfik65Loc96t9... *
  1       1383nt, >dmag4vel4xfik65Loc96t15... at 1:1383:1:1383/+/99.42%
  2       1383nt, >dmag4vel4xfik65Loc96t18... at 1:1383:1:1383/+/99.86%
  3       684nt, >dmag4vel4xfik81Loc3813t5... at 1:684:463:1146/+/99.12%
  4       864nt, >dmag4vel4xfik85Loc1522t3... at 1:864:379:1242/+/100.00%
  5       864nt, >dmag4vel4xfik85Loc1522t5... at 1:864:379:1242/+/99.31%
  6       774nt, >dmag4vel4xfik91Loc1049t2... at 1:774:373:1146/+/99.22%
  >Cluster 8635
  0       873nt, >dmag4vel4xfik45Loc1t9452... at 1:873:70:942/+/99.89%
  1       873nt, >dmag4vel4xfik45Loc1t9453... at 1:873:70:942/+/99.89%
  2       483nt, >dmag4vel4xfik45Loc1t9455... at 1:483:901:1383/+/100.00%
  3       483nt, >dmag4vel4xfik45Loc1t9479... at 1:483:901:1383/+/100.00%
  4       492nt, >dmag4vel4xfik55Loc412t36... at 1:492:892:1383/+/100.00%
  5       1383nt, >dmag4vel4xfik65Loc127t19... *
  6       240nt, >dmag4vel4xfik65Loc127t21... at 1:240:1144:1383/+/100.00%
  7       1383nt, >dmag4vel4xfik75Loc138t9... at 1:1383:1:1383/+/99.86%

=cut

sub readcdhit {
  my ($bother)= @_;
  
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq)= (0) x 10;
  my($ok,$fh)= openRead($bother);
  # my $fh;
  # if($bother =~ /\.gz/) { open(F,"gunzip -c $bother |") or die $bother; $fh=*F; }
  # elsif($bother =~ /stdin|^-/) { $fh= *STDIN; }
  # else { open(F,$bother) or die $bother;  $fh=*F; }
  %bspans=();
  
  my $iscdhit=0;
  my($cli,$ncl,$nalt,$mainid,$mainlen,$nerr,$hasdupid)=(0)x10;
  while(<$fh>) { 
    if(/^>/) { 

    ## move this dupid filter before into readblastab, readcdhit ?
if(DUPFILTER1) {    
    if( $hasdupid > 0 ) {
      foreach my $lt (sort keys %bspans) {
        if( my $dpmain= $dupids{$lt} ) {
          if($bspans{$dpmain} and $lt ne $dpmain) {  # ok to drop ..
            if($lt eq $mainid) { $mainid= $dpmain; }
            delete $bspans{$lt}; $ndupdrop++; 
          }
        }
      }
    }  
}
      putspans($mainid) if($mainid); %bspans=(); $hasdupid=0;
      ($cli)= m/(\d+)/;  $ncl++; #   m/Cluster\s*(\d+)/;  fixed or not?
      
    # elsif(/^(\d\w*)\s+(\d+)(..), >(.+)\.\.\. (.+)$/) # problem: >(ID\.1)\.\.\.
    } elsif(/^(\d+)/) { 
      my $i=$1;
      m/(\d+)(nt|aa), >(.+)\.\.\. (.+)$/; # problem: >(ID\.1)\.\.\.
      my($tlen,$typ,$tid,$pinfo)=($1,$2,$3,$4); 
      ## FIXME: merge.clstr:  1f1  9999aa, >id... at 100  << .f1,.f2 added // DROP this?
      #cdhit perls: /(aa|nt), >(.+)\.\.\./
      my $pi;
      my $ismain=($pinfo =~ /\*$/)?1:0;
      $hasdupid++ if($dupids and $dupids{$tid}); # need to read all cluster then drop dups
      unless($tid) {
        $nerr++;
      } elsif($ismain) {
        $mainid= $tid; $mainlen=$tlen; $pi=100;
        push( @{$bspans{$tid}}, [1,$tlen,1,$tlen,$tlen,$tlen,$tlen]); 
      } else {
        # pi for cdhit-est: at 1:492:892:1383/+/100.00%
        # pi for cdhit-aa : at 99.76%
        $pinfo =~ s/at //;  $pinfo =~ s/^\D+//; #? always \digit start
        ($pi)= $pinfo =~ m/([\d\.]+)%/; $pi=~s/\.00//; $pi||=0;
        my($qstart,$qend,$tstart,$tend)= (0) x 4;
        my @pinfo= split( m=/=, $pinfo);
        my @aln= split/:/,$pinfo[0]; #? multiple align segs or 1 only?
        my $aor=$pinfo[1]; # dont need?
        # 249nt, >dmag5vel5xco1k75Loc13984t1... at 249:1:421:669/-/100.00% << revalign
        if(@aln>3) { ($qstart,$qend,$tstart,$tend)= @aln; # what if @aln>4 ?
          ($qstart,$qend)= ($qend,$qstart) if($qstart>$qend); } 
        my $aln= 1+$qend-$qstart;
        my $aident= int($aln * $pi/100);
        my $bits= $aln; # aident? 
        push( @{$bspans{$tid}}, [$qstart,$qend,$tstart,$tend,$bits,$aln,$aident]); # add tlen?
        $nalt++;
      }
      # $iscdhit=2; 
    } elsif(/^\w/) {
      $nerr++; # warn/die not cdhit format ..
    }
  }
  
if(DUPFILTER1) {    
    if( $hasdupid > 0 ) {
      foreach my $lt (sort keys %bspans) {
        if( my $dpmain= $dupids{$lt} ) {
          if($bspans{$dpmain} and $lt ne $dpmain) {  # ok to drop ..
            if($lt eq $mainid) { $mainid= $dpmain; }
            delete $bspans{$lt}; $ndupdrop++; 
          }
        }
      }
    }  
}
  putspans($mainid) if($mainid); %bspans=();
  warn "# readcdhit: nclust=$ncl, nalt=$nalt, ndupdrop=$ndupdrop, nerr=$nerr \n" if $debug;
  return($ncl+$nalt); 
}

sub readblatpsl {
  my ($bother)= @_;
  
  my $ispsl=0;
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq, $nids)= (0) x 10;
  my($ok,$fh)= openRead($bother);
  # my $fh;
  # if($bother =~ /\.gz/) { open(F,"gunzip -c $bother |") or die $bother; $fh=*F; }
  # elsif($bother =~ /stdin|^-/) { $fh= *STDIN; }
  # else { open(F,$bother) or die $bother;  $fh=*F; }
  %bspans=();
  
  $ispsl=1 ;
  while(<$fh>) { 
    ## $ispsl=1 if m/^psLayout/; # should do but want some slack here?
    next unless(/^\d/);     
    chomp; my @v= split"\t";
    if(@v==22) { shift(@v); } # ucsc psl starts with a 'bin' field?
    unless(@v==21 and $ispsl){ die "# error: doesnt look like psl format I know" ; }
    my( $mat, $mis, $rep_matches, 
      $qgapw, $tgapw, $orient,
      $qid, $qsize, $qstart, $qend,
      $tid, $tsize, $tstart, $tend,
      $blocksizes, $qstarts, $tstarts,
      )= @v[0..2, 5,7, 8..16, 18..20];
    $qstart++; $tstart++; # move to 1-origin

    ## FIXME: allow for blat split rows per query, happens sometimes.., use bspans as for blast
    if($lq and $qid ne $lq) {  
      putspans($lq); $nids++; %bspans=();
    }

    if($tid eq $qid) { 
      # $qd=$qid unless($qd); $wq=$aln unless($wq); 
    } else { 
      # my $pmat= 100 * $mat / $qsize;
      # next if($pmat < $TINYALN); # NOT NOW with bspan parts
      my $aident= _max(0,$mat-$mis);  # not sure -gap right there, align should be +gap, ident = match
      my $aln= $mat + $qgapw; # is this right? align span includes query gaps, and mismatches
      my $bits=  $mat; # what? maybe 2*mat ? 2*aident? 
      #** should split blocks in qstarts, tstarts to sum up q/t start,end
      ## but dont need see putspans:       $tbit += $xbit; $taln+= $aln; $tidn+= $aident;
      
      $bspans{$tid}=[] unless(defined $bspans{$tid}); ## [$qstart,$qend,$tstart,$tend,$bits,$aln,$aident]; 
      # putspans($qid);  # only 1 pairmatch/span per psl row.
      push( @{$bspans{$tid}}, [$qstart,$qend,$tstart,$tend,$bits,$aln,$aident]); 
    } 
    ($lq,$lt)= ($qid,$tid);
  } close($fh);
  
  putspans($lq); $nids++;  
  return($nids);
}

