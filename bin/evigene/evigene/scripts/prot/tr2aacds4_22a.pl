#!/usr/bin/perl
# tr2aacds.pl

=item ABOUT

  EvidentialGene tr2aacds.pl
   
  This program converts mRNA assembly sets of too many transcripts, to best open
  reading frame coding sequences, with filtering of identical coding sequences, 
  and alternate CDS transcripts found by high identity subset alignments. 

=item EXAMPLE

  evigene/scripts/prot/tr2aacds.pl -log -cdna transcripts.fa -NCPU 8 -MAXMEM 16000 

  Input is one file of all transcripts.fasta[.gz]
  
  Output is folder of non-redundant transcripts, proteins and CDS,
  classified by gene locus, with tables of gene/transcript qualities,
  and a separate drop folder of redundant sequences.

  See RUN SCRIPT below for a cluster script.

=item REQUIREMENTS

  Current software required:
  EvidentialGene source code from http://eugenes.org/EvidentialGene/evigene/
            or http://sourceforge.net/projects/evidentialgene/
            packaged by release as evigeneDATE.tar          
  blastn, makeblastdb of NCBI C++ BLAST (tested ncbi2227 and later)
          from  https://blast.ncbi.nlm.nih.gov/ (Altschul et al. 1990)
  cd-hit, cd-hit-est, of https://github.com/weizhongli/cdhit/ or http://cd-hit.org/ (Li & Godzik 2006)
  fastanrdb, of exonerate package from https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate (Slater & Birney 2005)
    rapid and reliable non-redundant fasta sequence database
  These need to be in PATH or ENV. 

=item More ABOUT 

  This pipeline script comes out of many tests and refinements for
  selecting best transcript assemblies from a superset of many, many
  assemblies of same RNAseq data.  The focus on CDS of mRNA transcripts
  allows for several of those refinements in balancing too few/too long mistakes,
  and too many subset/nearly same assemblies.  The goal is a biologically
  valid transcript assembly set.

  The input mRNA assembly set is presumed a large collection from many
  assemblers, data slices and options from the same transcript source.
  Typically 1+ million input transcripts are processed, and reduced to a 
  biologically realistic set of 20k - 50k main plus similar number of
  alternate transcripts.
  
  CDS sequence is primary quality of mRNA transcripts, after trying/discarding 
  AA sequence (for those missed silent codon paralog changes).  CDS and AA
  are guessed as ~longest ORF of transcript (with sizes options for complete/partial).
  Long UTRs are also scanned for ORFs, of joined genes (common w/ hi-express neighbors,
  some assemblers).   
  
  CDS are classified by identity/alignment to each other, with fastanrdb (identicals), 
  cd-hit-est (ident. fragments) and blastn (local hi-identity alternates).
  High identity redundant assemblies are discarded, alternates classified
  by hi-identity local alignments are classified and CDS-duplicative ones
  discarded.   A final okay set contains main transcripts with informative
  alternates, and unique CDS (no alternate or low identity).  Minimum CDS/AA sizes
  are used. Protein completeness/partial, and UTR-poor qualities also
  score and filter excess assemblies.
  
  Transcript as whole is ignored due to many UTR misassemblies and difficulty
  in assessing quality, but for CDS/UTR ratio as quality measure.  Test show
  selection of longest CDS-ORF has strong +correlation with highest 
  protein orthology score.

  Use of high-identity CDS filtering (nrdb and cdhit-est at 100% identity, then
  blastn exon-level identity >=98%), balances te removal of redundant CDS 
  with retention of valid alternate, paralogs and clonal variant transcripts.
  
  This reduction can be run on subset transcript assemblies, then
  re-run on okay outputs to combine.  This script currently uses
  -NCPU for parallelization on clusters.  The slowest computes are
  for blastn-self of a large set.   Both use NCPU efficiently on one compute 
  node w/ many cores.  

=item SEE ALSO SRA2Genes

  Please also consider using SRA2Genes of EvidentialGene, which includes this tr2aacds
  as an important component. Read more in  evigene/docs/EvidentialGene_howto.txt 
  and  evigene/docs/evgpipe_sra2genes.help.txt
  
  tr2aacds is one component in a full gene reconstruction pipeline program.
  It has a limited function of efficiently reducing large, over-assemblies of genes, using 
  only intrinsic evidence (self-referential) of these sequences.
  
  A full gene reconstruction pipeline program like SRA2Genes has many more operations, and is more
  complex with choices of input gene evidence, measurements of protein homology, 
  contamination, and added outputs.  As named, this SRA2Genes fetches RNA-seq input data 
  from public repositories like NCBI SRA, to over-assemble, reduce, classify and analyze.  
  It also operates on your own RNA-seq or/and transcript assemblies.
  
  
=item PIPELINE ALGORITHM of tr2aacds 

  Stage 1 reduction, steps 0 to 5, are original tr2aacds algorithm.
  Stage 2 reduction, adds classifications steps, for version 4, tr2aacds4.
  
  Stage1 reduction
  0. collect input transcripts in one FASTA file, e.g. asm_name.tr, 
     and compute best ORFs (aa and CDS sequences), working mostly on CDS sequences

  1. perfect redundant removal:  fastanrdb  input.cds > input_nr.cds

  2. perfect fragment removal: cd-hit-est -c 1.0 -l $MINCDS ..

  3. blastn, basic local align hi-ident subsequences for alternate tr.

  4. classify main/alternate cds, okay & drop subsets, using evigene/rnaseq/asmrna_dupfilter2.pl
     .. merges alignment table, protein-quality and identity, to score okay-main, ok-alt, and drop sets.

  5. make final output files from outclass: okay-main, okay-alts, drops 
      okayset is for public consumption, drops for data-overload enthusiasts (may contain valids).

  Stage2 reduction (more details below)
  
  revised steps 0..5: 
    a. reset default options of Stage1 to retain more alternates/paralogs, both true and false.
    b. measure additional gene qualities, coding potential, assembler-replication (consensus).
    c. revise S2, perfect fragment removal, to include UTR in comparison to keep valid 5'alternate starts.
    
  6. preliminary, unreduced, gene-classified data set, via trclass2pubset.pl
  
  7. alternate exon pattern analysis, measure constituative, variable and fragment alternate patterns
  
  8. reclassify with added measures, for alternate/paralog/fragment and shortaa/noncoding, via trclass2pubset.pl
     Outputs of gene/locus classified sequences and tables, in okayset/ folder to preserve usage.  
     Stage1 reduction is preserved in okayset1st/

=item  VERSION 2016.07

  -- update classifier version
  -- utrbad,gapbad and related mods: problems dropping valid short aa, flip side is too many tiny frags
  
=item  VERSION 2017.12

  -- minor updates for evgpipe_sra2genes
  
  -- add opt to switch classifier opts
     between many-small-noblastp vs few-small-unless-blastp-hit 
     e.g. for 2-stage reduction, so don't miss the few, short ortholog genes
      -smallclass='keeptiny|droptinynoho' OR -tinynohomology=keep|drop
     ** test first, direct w/ asmrna_dupfilter3.pl ** 
     asmrna_dupfilter3.pl min-size opts (all ENV settings)
  $AAMIN =$ENV{aamin}||30; #was 40;   # for aacomplete, utrok
  $AAPART=$ENV{aapart}||100; # was 100; # for aapartial # 201402: THIS aapart cut is problem, drop.mains uniq orthologs here
  $AAMINBAD=$ENV{aaminbad}||60; #was 200;  # for utrbad class
  $AAMINPOO=$ENV{aaminpoo}||60; #was 100;  # for utrpoor class
  tr2aacds: should use ONLY aamin or mincds not both
  my $MINCDS = $ENV{MINCDS} || 90; # what? #maketraa3.sh: $AAMIN=40; $AMINPOO=100; $AMINBAD=200

=item VERSION 2018.06.18

  -- tr2aacds3.pl (orig tr2aacds2c.pl), majorish version change
  -- added, tested AACONS replication/consensus over assemblers (AACONS below)
   
=item VERSION 2019.08 .. 2019.12 

  -- tr2aacds4.pl (orig tr2aacds2d.pl), major revision of reduction filter algorithm
  -- in conjunction with reference comparison tests of 3 reduce filter types: 
      longest CDS, longest RNA, most RNA-seq.  See public doc doi: 10.1101/829184 
  -- old algo measurement + classify details changed
  -- new algo measure + class details added
  -- updates marked as UPD1908 below
  
  -- PHETERO/pHeterozygosity added UPD1912, to handle heterozyg. samples (eg. outbreed deer ticks),
      by reducing align identity cut-offs, by 1..9 %, which lumps more nearly-same CDS,
      and more high-identity transcripts are dropped.  Should be used only where heterozygosity
      is measured, values of 1 or 2 should work reasonably well (as for Ixodes)
     
=item UPD1908 explained

Recovery of accurate alternates and paralogs was lower than expected, 
and had a drop in number of valid transcripts per alt, par classes.
dgg examined implementation details using these comparisons,
and identified specific problems:
  
  p1. Perfect fragment removal was removing valid alternate 5'starts,
  by ignoring 5'UTR sequence where alternate exons are spliced to
  produce shorter, but fragment-identical coding sequences.  This is
  qualitatively different form of alternate splicing that the more
  common exon-alternation in CDS.  The updated implementation adds
  5'UTR for measure and reduction of perfect-fragment CDS alternates;
  UTR-only alternate splicing is still ignored by evigene methods.
  
  p2. Transcript aberrations includes a measure of percent CDS, where
  transcripts with long UTR, e.g. less than 50% CDS, were reduced by
  that measure.  However, coding genes with longer UTR are now a
  common occurrence in reference gene sets, increasing since 2010, as
  evidence of long-UTR transcript assemblies has been incorporated in
  reference gene sets.  Now this reduction weight is only applied to
  transcripts with less than 20% CDS.
  
  p3. Alternates with nearly identical proteins, as opposed to CDS,
  were reduced.  This is a mistake based on reference gene sets, with
  valid CDS alternates that produce near or identical proteins.
  
  p4. The ORF calculator option to distinguish complete versus 5'
  partial proteins was adjusted for better agreement with reference
  gene sets, increasing portion of complete proteins, with a related
  result that fewer CDS alternates, at 5' partial end, are retained.

These corrections have the general result of retaining more alternate
models.  New transcript set qualities have been examined and found
useful additions to this "longest CDS and valid alternates"
self-referential classifier:

  q1. Coding potential calculation, to help distinguish random
  sequence that has codon chains but is less likely to be valid
  protein.  This is common for the shorter putative CDS, where random
  results with large samples contain such codon chains.  Various
  coding potential calculations are known (eg Kang et al 2017), but
  are not perfect classifiers as many known reference proteins score
  as having no coding potential.  Evigene's implementation uses coding
  potential calculation as a modifier quality, to reduce low scoring
  transcripts in absence of other positive qualities (ie homology).  A
  subclass 'nc' for non-coding potential is added to the
  classifications and transcript measures.
  
  q2. Replication by assemblers of the same coding sequence
  transcripts is a form of technical validation, where replicated
  forms have higher likelihood as accurate than un-replicated
  alternate models [cite Voshall thesis].  This also is not a perfect
  classifier, as with many complex genes, a single assembly
  parameter/data set often produces the most accurate model.  Another
  type of replication score, across populations or individuals of a
  species, will be distinguished, analogous to biological and
  technical replication approaches used in gene expression analyses.
  
  q3. Alternate exon splicing patterns are assessed from
  self-alignment of gene alternates. This is an approximate measure
  that compares to alternate exon splicing found with chromosome
  mapped transcripts.  Valid alternates typically share common or
  constitutive exons, invalid models do not share these, or contain a
  fragment of the shared exons.  This exon splice pattern measure
  helps discriminate valid and invalid models of alternates.
    See also public doc doi: 10.1101/829184

=item UPD1908, 2019.08..11, tr2aacds2d.pl and asmrna_dupfilter3c.pl

  * improved best perfdup and perffrag choices, tr2aacds2c nrcheckaaqual,
    ie no longer picks shortest utr of cds-dups, but median utr size
  -- fewer perfect-fragments; modified perfect-fragment test to include 5'UTR in test of identical CDS fragments, 
      as this test was dropping most valid alternate start exons with splices in 5'UTR (proper name is what?)
  -- modified keep/drop settings, esp. keep more althi, high identity alternates as these are measured 
    to contain many biologically valid alternates.
  
  * added assembly replication quality score, keeps alternates and noclass when assemblers agree.
    replicate/agree/consensus evidence update, doesnt rescue frag/partial for agree>1  

  * bestorf ORF_FULLvPART default 0.85 should be reduced to 0.50 or lower for accurate uses,
    ref gene sets have few aapartials due largely to rarely calling long partial vs short full orfs.
  
  * aaqual table now includes Code/Noncode column from bestorf codepot qual.
  
  * added 'nc' class tag to tr2aacds classes, for pot. noncoding, 
    ignores 'noncode' tag if aasize >= 140, or aasize > 99 and not poor 
    keeps all 'nc' subclass of main/alt/noclass, but drops nc part/frag/utrorf,
    
  * modified default parameters: AADUP  99% from 98%, 
  * later removed this AADUP test entirely, removed 'a2' sublcass tag,
   as causing loss of valid alt/par w/o useful redundancy reduction.

   * tr classifier: reduced cost of long UTR (ie. aberrant or joined-gene score) as such are now a common
   occurence in reference gene sets (due to more evidence of long-UTR in transcript assemblies); 
   params pCDSOK=pCDSBAD=20 (min % of mRNA as CDS is okay quality)

=item UPD1908, 2nd stage reduction filter

  * overall, many more transcripts are retained in okayset as potentially valid, but wanting further
   evidence tests to reduce to less redundant set.  This 2nd stage reduction filter is done now in
   trclass2pubset.pl, with two variants: self-referential only (using only input trasm sequences),
   and adding external evidence of more complete SRA2Genes/Tr2Genes, with cross-species homology,
   contaminant checks, chromosome map checks.
   
  * in keeping value of  self-referential only reduction filter, this 2nd stage reduction should always
  be used w/ tr2aacds, now is included as step-7 of MAIN pipeline.  Test proc in "prot/tr2aacds4_stg2a.pl"
  as tr2aacds4_stage2b($outclass,$cdsblast);  That calls several evg single-method scripts, plus
  trclass2pubset that makes publicset/ sequences, with gene-tr-evidence table.

  * 2nd stage reduction filter adds new self-ref quality algorithm (genes/blasttrset2exons): 
  alternate/paralog exon splice site measurement, w/o chromosomes,
  by alt transcript alignments per gene locus, assessment of common/constituative exons shared across
  alternates, and quality scoring of alternates as to common/unique exon splice patterns.  This is
  an analog of chromosome-mapped exon splice pattern analysis.  It agrees with chr-map results enough to be 
  useful to help decide if alternate models are redundant/fragment models, 
  as tested with reference and assembled genes.
     
=item RENAME tr2aacds to what?

  'tr2aacds' is funny name for biologists, others, and 2019 update changes
  enough to warrant new name, but also keep tr2aacds as it functions like old vers.
  name1: tr2genes ? but this is only part of sra2genes omnibus pipe..
     tr2codinggenes ?   
     tr2genes_classifier ? no does more than class, and see asmrna_dupfilter, asmrna_altreclass
  name2: tr2genes_self_referential_filter is too wordy
         tr2genes_selfref ? tr2geneset ?
  evgpipe_tr2genes.pl should be symlink to evgpipe_sra2genes.pl,
     ie same program w/ diff start points.
  
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

=item RUN SCRIPT run_tr2aacds4.sh 

  This script can be used, or modified, to execute tr2aacds4 on a cluster batch system,
  or from Unix command line w/o batch system.  
  NOTE you need to configure required software in path evgapps=$HOME/bio/apps/ ,
  or equivalent w/ edits.  Required are 
    evigene/scripts/ cdhit/bin/  exonerate/bin/ ncbi/bin
  
  #! /bin/bash
  #----------------------------------------
  ### env trset=tevg9arathap.tr datad=`pwd` ./run_tr2aacds4.sh
  ### env trset=myspecies_all.tr.gz datad=path/to/data qsub -q normal run_tr2aacds4.sh
  #PBS -N tr2aacds
  #PBS -A PutAccountIdHere
  #PBS -l nodes=1:ppn=8,walltime=15:55:00
  #PBS -V
  
  if [ "X" = "X$ncpu" ]; then ncpu=8; fi
  if [ "X" = "X$maxmem" ]; then maxmem=32000; fi
  if [ "X" = "X$datad" ]; then echo "datad=?"; exit -1; fi
  if [ "X" = "X$trset" ]; then echo "trset=what?"; exit -1; fi
  
  evgapps=$HOME/bio/apps
  evigenes=$evgapps/evigene/scripts
  export PATH=$evgapps/cdhit/bin:$PATH
  export fastanrdb=$evgapps/exonerate/bin/fastanrdb
  export PATH=$evgapps/ncbi/bin:$PATH
  
  evapp=$evigenes/prot/tr2aacds4.pl
  traopts="-log"
  #x traopts="-debug -log"
  if [ "X" != "X$opt" ]; then traopts="$traopts $opt"; fi
  
  cd $datad/
  echo "#START `date` " 
  echo $evapp -NCPU $ncpu -MAXMEM $maxmem $traopts -cdna $trset
  $evapp -NCPU $ncpu -MAXMEM $maxmem $traopts -cdna $trset
  echo "#DONE : `date`"
  #----------------------------------------

=cut


=item AUTHOR
  
  don gilbert, gilbertd near indiana edu, 2013 .. 2019 for v4
  part of EvidentialGene, http://arthropods.eugenes.org/EvidentialGene/
  
  =======================================================================================================
  
=cut

use constant VERSION =>  '2022.04.05'; # '2022.01.20'; #  adds resolve_revaa
use constant UPD22a => 1; # 2022.04.05 upd test for -stranded fix
use constant UPD20j_resolve_revaa => 1;
use constant UPD1912 => 1; # PHETERO test
use constant UPD1908 => 1; use constant UPD1911 => 1; 
use constant AACONS => 1;  # aaconsensus: keep tr with consensus of identical aa across asm types

## '2019.12.07'; # '2019.11.26'; ## rev 4a = 11.20,  4b = 11.26, 4c=12.04
## '2019.08.30'; 4prelim = 2d 
## '2018.06.18'; # aaconsensus update: uncertain value; utrorf bugfix
## '2017.12.21'; 
## '2017.12.01' # minor upd for sra2genes, + manysmall/fewsmall+blastp opt
## 2016.07 utrbad,gapbad and related mods: problems dropping valid short aa, flip side is too many tiny frags
## '2016.05.03'; ## blast ncpu fail retry
## '2014.05.25'; #.15 ; rescue nrdup/frag discards using aaqual (complete+utrgood)
## Vs '2013.07.27'; # '04.15'; #04.07; 03.14'; # '.03.11'; '.03.06';
## 11mar.FIXME: need blastn -ungapped otherwise miss perfect match of parts ; e.g. alt-tr half=perfect, other=imperf
## 14mar: blast_ncpu using fasplit(query.fa,ncpu) instead of blastn -num_tasks ..
## 07apr: add -ablastab homolog scores for classifier
## 15apr: lastz not ready to replace blastn; testing option..
## 27jul: lastz added, maybe ready/better than blastn, maybe not .. compare more w/ blastn
##-----------------------------------

use FindBin;
use lib ("$FindBin::Bin/../"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
# this script is now in evigene/scripts/prot/..

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname);

use cdna_evigenesub; # replacing local subs 201405
# use cdna_proteins; #?

# cdna_evigenesub globals:
our $EVIGENES= $ENV{EVIGENES} || "$FindBin::Bin/..";  
our $EGAPP='tr2aacds';  
our $EGLOG='t2ac';
our $dryrun=0;  
our $DEBUG= $ENV{debug}|| 0;
our $EVGLOGH;
our %EVG_SEQSUFFIX; # UPD1911 cdna_evigenesub.pm required seq file suffices
# our (%genenames, %genedbxref, %genenamepct, %namedgenes, %cddnames, %pubids); # cdna_evigenesub.pm globals

use constant default_MINCDS => 90;
my $MINCDS = $ENV{MINCDS} || default_MINCDS;  
my $AAMIN  = int($MINCDS/3); # UPD1912, use only as option, keep sync w/ MINCDS

my @default_cdsblastopt= (98, 1e-19); # UPD1912, PHETERO test,  change if default
my $CDSBLAST_IDENT= 98; # 98 + -ungapped ; was 95; # or 98; # or 95? ** drop this back to 95 default; at least for high okay count species
my $CDSBLAST_EVALUE= 1e-19; # is this ok? will skip short exons of 100% ident .. how short?
## .. may be getting high okay main/noclass gene count from mixed strain/heterozygotes? need some checking.
## refinealtpar uses CDSBLAST_EVALUE=1e-5 
my $PHETERO= $ENV{PHETERO}||0; # UPD1912: should this be option? as for asmdupfilter -pHETERO=$phet

my $NCPU= 1;
my $MAXMEM= 2000; # was 1000, in Mb always 
my $ORGANISM=$ENV{ORGANISM}||$ENV{species}||"No_name"; # want for selfref pubids ? short Spenam Genspe Nonam

## UPD 2016.02: -tinyaln 35  bad (sometimes), need opts and diff test: tiny by bases not pct overlap
our $okayclass  ='main|noclass';
our $okaltclass ='alt|part';
our $dropclass  ='drop|cull';
# problem class: maybeok.althi1 (trclass); some have best refaa quals, but this is large subset
# maybeok counted as okay in sum tables.


# my $prefix = $ENV{prefix} || ""; # what?
my ($runsteps,$debug,$tidyup,$USE_LASTZ,$noutrorf,$reorient)= (0) x 9; # pm: $dryrun, 
my ($logfile,$aaseq,$aasize,$cdnaseq,$cdsseq,$aacdseq,$aaclstr,$aablast)= (undef) x 20; 
my (@okayset,@dropset,@inputset,@tmpfiles,@erasefiles); # tidyup file sets
my ($USE_AACONS, $aconsensus, )= (0,undef); # AACONS
my $smallclass=undef; # upd1712
my $smallclassmax= defined($ENV{'SMALLMAXAA'}) ? $ENV{'SMALLMAXAA'} : 100; # ?? or 120;
if(UPD1908) { $tidyup=1; } # new defaults; below $USE_AACONS=1;
$reorient= $ENV{DO_RESOLVESENSE}||0; # uses DID_RESOLVESENSE
my $RNAIsStranded= undef; # UPD22a: -strand[edrna]:yes|no|auto auto==maybe|dunno|test_it
## added opt, as -runsteps nofrag : skip_step_nofrag .. or as -skipstep='nofrag|step2,...' ? or as -stepopt='skipnofrag|skip2' ?
## like sra2genes  -runsteps=(no|yes)(frag|aadup|stage2) ie noaadup,nofrag,nostage2
## 201405: replace w/ subs: openloggit
## 201609: need -noutrorf to pass on to cdna_bestorf.pl; could do this? no, opt not allowed in test -x


my @saveopt= @ARGV;  
my $optok= GetOptions( 
  "cdnaseq|mrnaseq|trinput=s", \$cdnaseq, ## \@input,  # one only for this?
  "aaseq|aainput:s", \$aaseq,
  "cdsseq|cdsinput:s", \$cdsseq,
  "logfile:s", \$logfile,
  "runsteps=s", \$runsteps, # upd1805  
  "aconsensus:s", \$aconsensus,   # for AACONS, this has optional params, "-acons" alone means default consensus method
  "ablastab|anames=s", \$aablast,   # option traa-refaa.blastp.tall4 table for asmrna_dupfilter2.pl classifier
                # now also allow evigene .names table(trid, name, alignval, refid, ..)
  "tinynohomology|smallclass!", \$smallclass, # is this bool or val ? -notiny/-nosmall ?
   # upd1712 ^^ -smallclass='keeptiny|droptinynoho' OR -tinynohomology=keep|drop
  "MINCDS=i", \$MINCDS,  #? switch to MINAA? .. allow that as option at least
  "MINAA=i", \$AAMIN, # -AAMIN opt clash w/ -aa[seq]  
  "pHeterozygosity=i", \$PHETERO,  #UPD1912? 
  "SMALLAAMAX=i", \$smallclassmax,  # UPD1908 : should be CDS if using MINCDS
  "CDSBLAST_IDENT=s", \$CDSBLAST_IDENT, "CDSBLAST_EVALUE=s", \$CDSBLAST_EVALUE,  
  "NCPU=i", \$NCPU, "MAXMEM=i", \$MAXMEM,  
  "ORGANISM|species=s", \$ORGANISM, # UPD1911  
  #UPD1908:obsolete# "uselastz|lastz!", \$USE_LASTZ, # instead of blastn for cds-align
  "tidyup!", \$tidyup, 
  "noutrorf!", \$noutrorf,  # 201609
  "reorient!", \$reorient,  # UPD20j_resolve_revaa.. will be option-less default
  "strandedrna:s", \$RNAIsStranded,   # UPD22a
  "dryrun|n!", \$dryrun,
  "debug!", \$DEBUG, ## was $debug, ## $DEBUG now, dont need both
);


## temp fix for bad input cdnaseq: -s require file exists (stdin??) and loggit
die "EvidentialGene tr2aacds.pl VERSION ",VERSION,"
  convert large, redundant mRNA assembly set to best protein coding sequences, 
  filtering by quality of duplicates, fragments, and alternate transcripts.
  See http://eugenes.org/EvidentialGene/about/EvidentialGene_trassembly_pipe.html
Usage: tr2aacds.pl -cdnaseq transcripts.fasta[.gz] 
options: 
  -NCPU=$NCPU  -MAXMEM=$MAXMEM.Mb : number of CPU, max memory in megabytes (recommended)
  -MINAA=$AAMIN  : minimum valid protein size (or -MINCDS=$MINCDS),
      set -MINAA as high in 20..120 as reasonable for species, to avoid excess tiny coding-genes 
  -pHeterozygosity=[0..9] : reduce percent identities for heterozygous organism sample (default 0),
      this lowers alternate/paralog identity cutoffs and classes, and drops more high-identity transcripts
  -logfile       : write log of operations  (recommended); 
  -species=Homo_sapiens : define species, now used for public IDs only
  -ablastab=blastp_table : use blast scores to keep some transcripts (blastp traa x refdb -outfmt 7; see docs)
  -aconsensus    : use replication over assemblers as quality (default) or -aconsensus=no (dont use), 
  -runsteps=list : skip or do certain pipeline steps (nofrag|yesaadup|nostage2, see docs/code)
  -strandedrna   : auto|yes|no; default auto == maybe|dunno, tests cds-orient, reliable but must recalc cds if stranded
  -reorient      : resolve mixed strand alternates from evidence of map strand and code potentials
  -noutrorf      : dont look for ORF in long UTRs, from fused genes (ie 2nd longest ORF)
  -debug         : more output to logfile
" unless($optok and $cdnaseq and -s $cdnaseq); 
#was   -MINCDS=$MINCDS     : minimum coding sequence size

## need retesting before common use
#  -[no]smallclass: option to keep small ORFs (>=20 aa), or skip small ORFs (<90 aa), default is neither
#  -[no]tidyup    : relocate outputs nicely (default $tidyup); -dryrun : list operations w/o action


openloggit( $logfile, $cdnaseq); 
loggit(1, "EvidentialGene tr2aacds.pl VERSION",VERSION);
loggit(1, "CMD: tr2aacds.pl ",@saveopt);
loggit(1, "ERR: unused arguments:",@ARGV) if(@ARGV>0); # do something w/ remaining @ARGV .. warn

$debug= $DEBUG; # dont need both
$tidyup= 1 unless($dryrun); #was ||$debug;  default on unless debug|dryrun ?

# see above: $USE_AACONS=0; or 1 default?  
if(defined $aconsensus) {
  if($aconsensus=~/^(no|0)/) { $USE_AACONS=0; }
  else { $USE_AACONS=1; }
} else {
  if(UPD1908) { $USE_AACONS=1; $aconsensus=""; } # UPD1908: new default 
}

if($MINCDS != 3*$AAMIN) { # UPD1912, AAMIN only as option, keep sync w/ MINCDS
  if($MINCDS == default_MINCDS) { $MINCDS= 3*$AAMIN; }
  else { $AAMIN= int($MINCDS/3); } 
}

if(UPD1911) {
  # require("prot/tr2aacds4_stg2a.pl"); # now inline below 
  # #  tr2aacds4_stage2b($outclass,$cdsblast);

  ## need these new defaults in stage1 (also stage2), used in sub asmdupfilter_cds()
  unless($ENV{'asmrna_altreclass'}) {
  $ENV{'asmrna_altreclass'}="$EVIGENES/rnaseq/asmrna_altreclass4.pl"; # for trclass2pubset, fix this how?
  ## see below  findevigeneapp("rnaseq/asmrna_altreclass4.pl") ..
  }
  $ENV{ORGANISM}= $ORGANISM if($ORGANISM); # pass this way if from option
  
  ## ALTFRAG cutoff is question, old def == 0.50
  my %updef= ( pCDSOK => 20, pCDSBAD => 20, ALTFRAG => 0.50, CDSXUTR =>  "900,60");
  if($MINCDS > default_MINCDS) {  # > 99; nosmall retest
    # my $aam= int($MINCDS/3); # now == $AAMIN
    my %aaenv=(aamin => $AAMIN, aapart => $AAMIN, aaminbad => $AAMIN, aaminpoo => $AAMIN); 
    for my $uk (keys %aaenv) { $updef{$uk}= $aaenv{$uk}; }
  }
  for my $uk (keys %updef) { $ENV{$uk}= $updef{$uk} unless($ENV{$uk}); } # no zeros

  ## TEST classify results of reduced identities for heterozygous trasm sets
  if(UPD1912) { 
    # my @default_cdsblastopt= (98, 1e-19); # UPD1911, PHETERO test, change only if default
    # my $PHETERO = $ENV{PHETERO}||0; # now option
    # also add -runsteps=yesaadup, if default
    if($PHETERO>0 and $PHETERO<9) { 
      $CDSBLAST_IDENT -= $PHETERO if($CDSBLAST_IDENT == $default_cdsblastopt[0]); 
      $CDSBLAST_EVALUE = 1e-14    if($CDSBLAST_EVALUE == $default_cdsblastopt[1]); 
    } else {
      $PHETERO=0; # cancel out-of-range choice
    }
  }

  # UPD1911 for [no]smallclass .. what? user opt?  MINCDS enough?
  #?? change smallclass to filter out SMALLMAX >= smallaa > MINAA unless( has_aacon or has_homol or has_evidence)
  #? add as new reduction step b/n reduce nr and reduce frags
  
  ## NCPU,MAXMEM  .. add calc for sensible NCPU/asm script, dont want 100 cpu x 50 GB mem calls
  use constant MEMperCPU => 2000; # default MAXMEM = 2000, change? 8cpu x 2g = 16g total ..
  use constant baseCPU   => 8;  
  $NCPU=1 if($NCPU<1); $MAXMEM= MEMperCPU if($MAXMEM<MEMperCPU);
  if($NCPU > baseCPU) {
    my $cpm= int($MAXMEM/$NCPU); my $ncpum= $NCPU;
    if($cpm < MEMperCPU) { $ncpum= int(0.5 + $MAXMEM/MEMperCPU); $ncpum= baseCPU if($ncpum < baseCPU); } 
    if($ncpum < $NCPU) { loggit(1,"Reset NCPU=$ncpum due to lower MAXMEM=$MAXMEM"); $NCPU=$ncpum; }
    }
}


my $APPblastn=    findapp("blastn"); 
my $APPmakeblastdb= findapp("makeblastdb");
# my $APPlastz="echo MISSING_lastz"; # Obsolete UPD1908 
#  	 $APPlastz= findapp("lastz") if ($USE_LASTZ); # nodebug $debug; 
# 	 ## lastz 04.15 testing still; 2013.03.24 : replace blastn default for basic local align: better at finding perfect local aligns
my $APPfastanrdb= findapp("fastanrdb");
my $APPcdhitest=  findapp("cd-hit-est"); 
my $APPcdhit=     findapp("cd-hit");  
## .. these should call findevigeneapp() to warn/die if missing
my $APPcdnabest= findevigeneapp("cdna_bestorf.pl"); # allow ENV/path substitutions? ? UPD 2016.07? or cdna_proteins.pm
my $APPtraa2cds= findevigeneapp("prot/traa2cds.pl");
my $APPaaqual=    findevigeneapp("prot/aaqual.sh"); # replace w/ cdna_evigenesub:makeAaQual()

# ACONS: asmrna_dupfilter3c.pl is proper version, for now, fix here or not?
# should update to asmrna_dupfilter4 ..
# my $APPtrdupfilter= findevigeneapp("rnaseq/asmrna_dupfilter3c.pl"); # OLD
my $APPtrdupfilter= findevigeneapp("rnaseq/asmrna_dupfilter4.pl"); # UPD to asmrna_dupfilter4.pl

my @UPD4APPS;
if(UPD1908) {
  # tr2aacds4_stage2b uses makeblastscore.pl rnaseq/asmrna_altreclass3c.pl prot/cdsqual.sh 
  #  genes/blasttrset2exons2.pl genes/trclass2pubset.pl
  my($addapp);
  ## FIXME: asmrna_altreclass3c should be aliased .. asmrna_altreclass4 ?
  ## FIXME: ditto blasttrset2exons2 > blasttrset2exons
  if( $addapp= findevigeneapp("rnaseq/asmrna_altreclass4.pl")) { push @UPD4APPS, $addapp; }
  if( $addapp= findevigeneapp("makeblastscore.pl")) { push @UPD4APPS, $addapp; }
  if( $addapp= findevigeneapp("prot/cdsqual.sh")) { push @UPD4APPS, $addapp; }
  if( $addapp= findevigeneapp("genes/blasttrset2exons2.pl")) { push @UPD4APPS, $addapp; }
  if( $addapp= findevigeneapp("genes/trclass2pubset.pl")) { push @UPD4APPS, $addapp; }
}

## FAIL at this point if any apps missing?
#-------------------------------------

=item FIXME 201402 drop.mains problem

  FIXME 201402: step5. need perfect_frag/dup info input to step4 asmdupfilter_cds, to replace poor qual drop.mains
  asmdupfilter_cds: AAPART cutoff is problem, with UTRBAD/POOR, causing drop.mains of uniq orthologs 
  update for drop.mains, perfectfrag replacements, needs input list from tr2aacds 
    of perfectfrag/dups from fastanr/cdhit prior steps
    
  asmdupfilter_cds update to relax AAMINBAD,.. dropping, tests suggest options for this
    -AABADMIN=60 -AAPOOMIN=60 -AAPARTMIN=80 -AAMIN=40 -pCDSOK=50 
  will about double main/noclass okay count, but keep most of shortish/utrpoor but thru orthogenes,
  not much change to alt classing/count (still dropping poorish alts).
  
=item FIXME 201402 tblastn-refaa screening needed for orthogene completeness
 
  /bio/bio-grid/aabugs4/evg6c/reft/catfishmiss, reft/outz/
  blp1catfishevg8aa-refvert_kogkfish.blastp.gz  :  blastp refvert.aa x evg8.aa  : adds 300 orthomiss
  blp1catfishevg8tr-refvert_kogkfish.tblastn.gz :  tblastn refvert.aa x evg8.tr : adds 113 orthomiss

  orthogene blastp x evg8.aa,tblastn x evg8.tr
  evg8.aa recovers 300 orthogene misses of evg1, of 415  catfish1evg.nohitmissed.ids
  evg8.tr adds 62 of evg1.misses (mostly short-revaa?), for 362 of 415 in all input catfishasm.aa
    *and* adds 51 of 544 orthogene miss of all catfishasm.aa, mostly/all? shorter revaa
    for 413 orthogene found, of 959 orthogene misses (catfish1evg.nohitmissed.ids catfish1asm3.nohit.ids)
    * some of tr hits span apparent inframe stop codons, need way to check on those.
  -- together these add 148 missing common ortho groups of Catfish=0, min taxa=9, moving catfish up a few levels
    in completeness, below/near median, in kfish2/prot/fish11c/fish11gor3/fish11gor3-orthomcl-gcommon.tab
  ** need tblastn refaa x alltrasm.tr as early step in tr2aacds (option..)
    -- in 3 substeps: 
      ra0: prefilter alltrasm.tr w/ cdhit-est -p 0.99? 
            dont need all tr, but collect all IDs to reapply tblastn results, this may reduce reduncancy
      ra1: if -inrefaa, tblastn -q refaa -db alltrasm1.tr.nsq ..
      ra2: if -inreftblastn or ra1, retabulate reftblastn for alltrasm IDs, 
            reftbln include align-spanbe for revised bestorf to catch revaa, utrorf and inframestops, 
                    refIDs and scores for asmdupfilter classing, replacing blastp input
      ra3: if -inreftabn or ra2, calc bestorf .aa,.cds adding revaa and alignspan/inframestop info for bestorfs,
            calc new mrna using utrorf, revaa, inframestop info, splitting utrorf mrna parts
            
=cut

sub failIfEmpty { 
  my($nrecout,$label,$seqfile)=@_;
  if($nrecout<1) { loggit(LOG_DIE, "ERR: failed ", $label."=",$seqfile,"nrec=",$nrecout); } 
  else { loggit(0, $label."=",$seqfile,"nrec=",$nrecout); }
}


sub MAIN_start {}

MAIN: {
  # final output is classify table from asmrna_dupfilter2
  # and seq filesets classified: perfect dups, perfect fragments, alt partial dups : keep,drop class

  $ENV{TMPDIR}= './' unless($ENV{TMPDIR} and $ENV{TMPDIR} ne '/tmp');  # mostly for big sorts 
  
  ## set TMPDIR to localdir always? dont use problematic perl mods
  # unless($ENV{TMPDIR}) {  ## FIXME sort TMPDIR, cant use /tmp
  #  require Cwd; if(my $cdir=Cwd::getcwd()) { $ENV{TMPDIR}=$cdir; } 
  #  # 201506: perl fail cwd2 .. require Cwd::getcwd()
  #  #o# if(my $cdir=FindBin::cwd2()) { $ENV{TMPDIR}=$cdir; } 
  #  }

 
  my($cdsseqnr,$cdsseqnrcd1,$cdsblast,$outclass,$outaln); # intermediate/output files
  my(%perfect_dups,%perfect_frags,$smallaa_dropids);
  
  loggit(0, "BEGIN with cdnaseq=",$cdnaseq,"date=",`date`);
  # fail unless -s $cdnaseq; .. is STDIN ok for this?
  
  my %steprun=(); # see s2g:  per {stepnum} options
  # UPD1908: ? make default:  -runsteps=noaadup, reverse w/ -runsteps=okaadup ?
  # if(UPD1908){ unless($runsteps) { $runsteps="noaadup"; } }# not good enough, other -runstep= cancels default
  # $steprun{7} == tr2aacds4_stage2b
  if($runsteps) {
    my @rs= split /[,;\s|]/, $runsteps; 
    for my $rs (@rs) {
      if($rs =~ m/^(skip|no|yes|do)(\w+)/i) { my ($noyes,$s)=($1,$2);
      if($s=~/frag|^2\b/){ $s=2; } # nofrag | defrag | fragremoval ..
      elsif($s=~/(aadup|aacluster|41)/) { $s=41; }
      elsif($s=~/^(stage2|7)/) { $s=7; }
      $steprun{$s}= ($noyes =~ m/yes|do/)?'yes':'no';
      }
    }
  }
  if(UPD1908) { 
    if(UPD1912 and $PHETERO) { 
    $steprun{41}='yes' unless(defined $steprun{41});  # PHETERO default
    }
    $steprun{41}='no' unless(defined $steprun{41});  # this way ensures default
    }
  
  if(defined $smallclass) { #upd1712 .. UPD1908 turn this off?? replace w/ other variant
    #   asmrna_dupfilter3.pl min-size opts (all ENV settings)
    #   tr2aacds MINCDS: should use ONLY aamin or MINCDS not both
    # MINCDS or aamin is primary? MINCDS limits bestorf calls
    my $MINAA= int($MINCDS/3);
    my %aaenv=(aamin => $MINAA, aapart=>100, aaminbad=>60, aaminpoo=>60); # current defaults in asmrna_dupfilter3
    if($smallclass) {
      %aaenv=(aamin =>20, aapart=>40, aaminbad=>20, aaminpoo=>20); 
      $MINCDS= 3 * $aaenv{aamin};
    } else { # nosmallclass ; require aablast or not?
      %aaenv=(aamin =>90, aapart=>120, aaminbad=>120, aaminpoo=>90); 
      $MINCDS= 3 * 20; # FIXME: NOT for -nosmallclass, need smalls from get_bestorf + aablastp
    }
    #X $MINCDS= 3 * $aaenv{aamin}; # FIXME: NOT for -nosmallclass, need smalls from get_bestorf + aablastp
    
    my $log="MINCDS=$MINCDS,"; # used for get_bestorf(), elsewhere?
    for my $k (keys %aaenv) { my $v=$aaenv{$k}; $ENV{$k}=$v; $log.="$k=$v,"; }
    loggit(1,"set smallclass=",$smallclass," to keep/drop tiny-aa:",$log); 
  }
    
  
  #  0ra. FIXME 201402 insert tblastn-ref.aa-asm.tr here, need for get_bestorf asm.aa,cds
  # UPD1908: ORF_FULLvPART of cdna_bestorf is important param, "standard" orf callers call ~all complete prots
  #     maybe change default ORF_FULLvPART from 0.85 to 0.50 or 0.33 ..
  # UPD1908: cdna_bestorf -codepotential new option, codepot=Noncode in classifier to filter tiny cds
  
  # 0. make/collect input asm_name.{tr,aa,cds}, working mostly on .cds here
  # FIXMEd: parallelize this for NCPU by splitting input cdnaseq to NCPU parts; cat parts> whole.
  ($cdsseq,$aaseq) = get_bestorf($cdnaseq,$aaseq,$cdsseq, $RNAIsStranded);
  
  # FIXME1606: fail if(facount(cdsseq)==0), add below also, fail_if(facount(outfile)==0)
  failIfEmpty( facount($cdsseq),"bestorf_cds=",$cdsseq);
 
  # FIXME: Check DUP IDs in aaseq, FAIL if found ..
  ($aasize)= make_aaqual($aaseq); # save global hash in $AAQUALH= getAaQual($aasize) ??
  if(my $ndup= fadupids($aaseq)) { loggit(LOG_DIE,"ERR: $ndup duplicate ids in $aaseq\n"); }
  
  #UPD20apr:
  # if($RNAIsStranded) ... redo get_bestorf(option=forward) ? use with -reorient
  #   # our $aaqualh= getAaQual($aasize); # global hash $AAQUALH, read once  
  #   my($ismrna,$nfwd,$nrev,$nt)= isStrandedRNA(undef,$aasize);
  #  -strandedrna : auto-detect from aaqual, dont need user opt?
  # see below UPD20j_resolve_revaa -reorient,  trclass_resolve_strandmix
  if(UPD22a) {
    if(defined $RNAIsStranded and $RNAIsStranded =~ /1|yes/i) {
      loggit( 0, "isStrandedRNA=$RNAIsStranded, option -stranded=yes \n");
    } elsif(defined $RNAIsStranded and $RNAIsStranded =~ /0|no/i) {
      loggit( 0, "isStrandedRNA=$RNAIsStranded, option -stranded=no \n");
    } else {
      # my($ismrna,$nfwdstr,$nrevstr,$ntstr)= isStrandedRNA($AAQUALH); #? always test?
      my($aaqualhash)= getAaQual($aasize); #??? ,$aaqualhash); # cdna_evigenesub;
      my($ismrna,$nfwdstr,$nrevstr,$ntstr)= isStrandedRNA($aaqualhash,$aasize);

      loggit( 0, "isStrandedRNA=$ismrna, f:$nfwdstr/r:$nrevstr\n");
      if($ismrna) {
        $RNAIsStranded=1;
        loggit( 0, "recalculating CDS for stranded mRNA inputs\n");
        rename($aaseq,"$aaseq.old"); rename($cdsseq,"$cdsseq.old"); rename($aasize,"$aasize.old");
        push @tmpfiles, "$aaseq.old","$cdsseq.old","$aasize.old"; 
        ($cdsseq,$aaseq) = get_bestorf($cdnaseq,$aaseq,$cdsseq, $RNAIsStranded);
        ($aasize)= make_aaqual($aaseq); 
        my $strflag= makename($cdnaseq,".is_stranded"); # if(defined $cdsseq and not $cdsseq);
        runcmd("touch $strflag");
        push @inputset, $strflag; # keep flag w/ input.cdna/cds
      }
    }
    
  }
  
  # these subtasks will check existence of output file before run task
  # 1. nonredundant removal:  fastanrdb  input.cds 
  ($cdsseqnr)= nonredundant_cds($cdsseq);
  failIfEmpty( facount($cdsseqnr),"nonredundant_cds=",$cdsseqnr);
  my $cdsseqnr_orig= $cdsseqnr; # keep name for counts
  
  ## 1.1. reassign redundant cds from aaqual = pCDS to reduce utrbad/utrpoor set.
  ## .. need only change header ID in $cdsseqnr ; fastanrdb puts all redundant ids on header.
  ## Note this rewrites $cdsseqnr file, to same name
  ## Note2 nrbest also tabulates agreeset/consensus/replication .. add opt to merge subset agree tables
  my($nagree,$agreeset) = (0,{});
  my($nhashomol,$homolset) = (0,{});

  if(UPD1912) {
    my($maind,@agrees)= getFileset('./','agree$|consadd$|aacons$');   # consensus
    for my $agf (@agrees) { my($nadd,$agreesetr)= add_consensus_idset($agf, $agreeset); $nagree+= $nadd; }
  }
  
  my($nbest,$ndups)= nonredundant_reassignbest($cdsseqnr, $aasize, $agreeset);
  loggit(0,"nonredundant_reassignbest=",$nbest,"of",$ndups); 

  # 1.2.  smallaa-class reduction step, require evidence to keep, only if have evidence at this step
  if(UPD1908) { 
    my($maind,$agf)= getFileset('./','consensus$'); 
    if($agf) {
      ($nagree,$agreeset) = add_consensus_idset($agf, $agreeset);  
    } else {
      ($nagree,$agreeset) = consensus_idset($cdsseqnr,'agreeonly'); # read headers # BAD for add_cons
    }
    
    ($nhashomol,$homolset) = homology_idset($aablast) if($aablast); # read headers  

    # UPD1908:?? new smallclass reduction step b/n  nonredundant_cds() and nofragments_cds() 
    #  remove SMALLMAXAA >= smallaa >= MINAA unless( has_evidence( aahomol | aaconsensus | other) )
    
    #upd20ja: trclass2pubset.pl uses $ENV{aaminnoclass}||60; for asmrna_altreclass -noclasscut=60
    if($smallclassmax) { # match aaminnoclass to this, but if no evidence, use lower cut?
      my $nca= $ENV{aaminnoclass}||0;
      if($nagree or $nhashomol) { $ENV{aaminnoclass}= $smallclassmax if($nca < $smallclassmax); } # 100 default here, 60 in trclass2pubset .. which ?
      # else { $ENV{aaminnoclass}= $smallclassmax if($nca > 0 and $nca < $smallclassmax); }
    }
    if($smallclassmax and ($nagree or $nhashomol)) {
      my($nok, $ncut, $cdsseqnrUPD, $smalldropidh)= 
        evd_smallclassfilter( $cdsseqnr, $aasize, $smallclassmax, $nagree, $agreeset, $nhashomol,$homolset);
      loggit(0,"evd_smallclassfilter[$smallclassmax aa] cut:",$ncut,"kept:",$nok, "on evidence of nagree:",$nagree,"nhomol:",$nhashomol); 
      if($nok and $ncut) {
        push @tmpfiles, $cdsseqnr; #? or rename in smallclassfilter
        $cdsseqnr= $cdsseqnrUPD;
        $smallaa_dropids= $smalldropidh;
        }  
      } 
  }

    
  ## 1.9. extend cds into utr to detect cds-end alts
  ## .. xutr doesnt seem useful in blastn align tab and dupfilter classing ..
  ## .. change to use -extend2utr, cdsseqnrxu only in nofragments_cds/cdhit ?
  ## OR postprocess cdhit.clstr and test perfdups for valid 5'alt utr? another align step
  
  my $cdsseqnrxu= $cdsseqnr;
  if(UPD1908) { 
    # TEST align CDS extended into utr, for steps 2/fragcds, 3/blastn, to detect mainly 5'start-alts
    # need only cdsnr ids here .. separate aaseq
    my $CDSXUTR= $ENV{CDSXUTR}||'900,60'; # tested default .. was 120;
    my $cdsnrids= faidlist($cdsseqnr,'',"ashash"); 
    my $aaseqnr = faextract($aaseq, makename($cdsseqnr,".aa"), $cdsnrids,1);
    push @tmpfiles, $aaseqnr;
    ($cdsseqnrxu) = make_cdsFromTrAA($cdnaseq, $aaseqnr, makename($cdnaseq,"nrxu.cds"),"-extend2utr $CDSXUTR");

    ## BUG: lost consensus val, agree=nnn, from $cdsseqnr headers ; USED to rescue shorter cds in nofragments_cds 
    #above# ($nagree,$agreeset) = consensus_idset($cdsseqnr,'agreeonly'); # read headers  
    if($nagree) {
      my $tmpxu="$cdsseqnrxu.tmp"; open(FX,$cdsseqnrxu); open(OX,">$tmpxu");
      while(<FX>){ if(/^>(\S+)/){ if(my $ag=$agreeset->{$1}){ s/$/ agree=$ag;/; } } print OX $_; }
      close(FX); close(OX); rename($tmpxu,$cdsseqnrxu); 
      }    
    push @tmpfiles, $cdsseqnrxu;
    loggit(0,"cds -extend2utr ",$CDSXUTR,"bp in",$cdsseqnrxu); 
  }

  ## FIXME1702: fastanrdb makes long hdr >id1 id2 id3 id3 ... id20000 ... and blast whines about that
  ## see below %perfect_dups = redundant_idset(),   %perfect_frags= fragment_idset()
  
  
  
  # 2. perfect fragment removal  : clusterize : $NCPU, $MAXMEM
  # 2.FIXME: some perf-frag maybe/are aabest (complete) vs perf-longest (aapoor), reassignbest cdhit.clstr + aaqual, 
  # 2.1805: add -nonofrag, ie skip nofragments_cds() step on request (may be dropping valid homol aa)
  my $cdsnrids= {}; # hash ref now, or id file?
  my ($fragidsetdata);
  
  if($steprun{2} eq "no") { # if($skip_step_nofrag) ;; -runsteps='skipnofrag|skip2'
    #O# $cdsseqnrcd1= (UPD1908)?$cdsseqnrxu:$cdsseqnr;
    $cdsseqnrcd1= $cdsseqnr; # UPD1908 skip xu here ?
    $fragidsetdata= "$cdsseqnrcd1.clstr"; # no such data
    $cdsnrids= faidlist($cdsseqnrcd1,'',"ashash"); 
    loggit(0, "skip step2. nofragments_cds=",$cdsseqnrcd1);
  } else {  
    ($cdsseqnrcd1,$cdsnrids)= nofragments_cds( (UPD1908)?$cdsseqnrxu:$cdsseqnr );
    failIfEmpty( facount($cdsseqnrcd1),"nofragments_cds=",$cdsseqnrcd1);

    # UPD1908 replace xu here with cds-only subset for blastn/further work ?? try xu only in nofragments_cds
    $fragidsetdata= "$cdsseqnrcd1.clstr"; # preserve name for fragment_idset()
    # below: %perfect_frags= fragment_idset("$cdsseqnrcd1.clstr"); # read headers after rebest, for dupclass_fileset
    if(UPD1908) {
    my $cdsnr2c = faextract($cdsseq, makename($cdsseqnr,"cd1x.cds"), $cdsnrids,1);
    push @tmpfiles, $cdsnr2c;
    $cdsseqnrcd1= $cdsnr2c; # BUT cant change file name for xxxcd1.cds.clstr
    }
  }

  # UPD1908: add step 2.05; filter out excessive short things that appear noncoding by cdsqual.pl tests 
  # combine w/ ok_noncode set methods, use similar reduction classing (perfrag,altblast) on ok_nc_rna set
  # add cdna_bestorf measure: longorf diff from other 5 frame orfs, large diff for good orf, small for poor?
  # also include in followup prot homol tests to ensure no mistakes in ncrna class
  #   ($cdsseqcodeok,$cdsqual) = noncode_shortcut($cdsseqnr); OR noncode_class($cdsseqnr, $cdnaseq)
  # and add class/count info on those
  #  %noncode_shortcuts= noncodecut_idset($cdsqual);
    
  # 2.1. FIXME need to check the cdsseqnrcd1.clstr clusters for useful alternate tr (?) 
  # slightly shorter CDS can be valid alts OR best main, use aaqual to rescue
  #  nofrag_reassignbest("$cdsnrseq.clstr", ...); see nofragments_cds()
  # FIXME.201405 reassignbest replace old 2: nr_reass and nfrag_reass .. or not, correct other 2 instead 
  #   my($nbest,$ndups)= reassignbest_aaqual($cdsseqnrcd1,$aasize);
  # * this affects other parts, redundant_idset & frag_idset below ..
  
  # UPD1908 drop lastz_cds
  # if($USE_LASTZ) {
  #   # 3z. lastz replace blastn for alignments of hi-ident subsequences ; tested but unsure result more valid than blastn
  #   ($cdsblast)= lastz_cds($cdsseqnrcd1);  # asmrna_dupfilter2.pl can now parse lastz.general.format
  #   loggit(0, "lastz_cds=",$cdsblast);
  # } else 
  
	{
    # 3. blastn alignments of hi-ident subsequences : clusterized  : $NCPU, $MAXMEM
    ($cdsblast)= blastn_cds($cdsseqnrcd1);
    loggit(0, "blastn_cds=",$cdsblast);
  }


  #  dropping valid homol prots, ie shorter than relative, but cluster with at 90% ..
  # 4.1 cdhit -c 0.9 -i $aaseq -o $aaseqcd .. for aaseqcd.clstr input to dupfilter ..
  #  .. filter uses aaclstr only to drop aa-hiident alts as uninformative; eg 1/2 of 13k althi.
  # FIXME 4.1 : reduce aaseq to $cdsseqnrcd1 id subset, otherwise large inputs become very slow here..
  # .. extract id lists from each reduce step?, at least nofragments_cds
  # 4.1,1805: FIXME add -runsteps="xxx,skipaacluster|skip4.1?" : cdhit -c 0.90 seq.aa  ends up 

  if($steprun{41} eq "no") { #  -runsteps='skip(aadup|aacluster|4.1)'
    ($aacdseq,$aaclstr)=(undef,undef);
    loggit(0, "skip step4.1 aadup clustering");
  } else {
    my $aaseqnr= makename($cdsseqnrcd1,".aa");
    $aaseqnr= faextract($aaseq,$aaseqnr,$cdsnrids,1);
    push @tmpfiles, $aaseqnr;
    ($aacdseq,$aaclstr)= aacluster($aaseqnr); # but use only aaseq IDs from cdsseqnrcd1  
  }
   
  # 4. classify main/alternate cds, okay & drop subsets, using evigene/rnaseq/asmrna_dupfilter2.pl
  ($outclass,$outaln)= asmdupfilter_cds($cdsblast,$aasize,$aaclstr);
  loggit(0, "asmdupfilter_cds=",$outclass);
  # FIXME? here, or in asmrna_dupfilter2.pl, create ID main,alt table from outclass/trclass, 
  #    with new numeric IDs, old/cur ids, main/alt num
  ## Better leave to other script, evgrna2genbanktsa.pl, merging trclass, naming table, other annots ?
    
  # 5. make final output files from outclass: okay-main, okay-alts, drops 
  # 5. fixme: add hdr classinfo from drops: %perfect_dups, %perfect_frags, 
  # FIXME 201402: step5. need perfect_frag/dup info input to step4 asmdupfilter_cds, to replace poor qual drop.mains
  # FIXME: fragment_idset nrcheckaaqual better changes clstr main ids..
  
  ## redundant_idset (UPD1908) ? $cdsseqnrxu : $cdsseqnr << here want orig cdsnr not nrxu 
  %perfect_dups = redundant_idset($cdsseqnr_orig); # read headers after rebest, for dupclass_fileset

  # $smallaa_dropids = ; # see above
  
  %perfect_frags= fragment_idset($fragidsetdata); #was ("$cdsseqnrcd1.clstr"); # read headers after rebest, for dupclass_fileset

  #new1712# 
  my $OUTH= (defined $EVGLOGH) ? $EVGLOGH : *STDOUT;
  my $OUTB= undef; my $oksumt= open(OUTBT,">$outclass.sum.txt");  $OUTB=*OUTBT if($oksumt);

  asmdupclass_sum($OUTH,$OUTB,$outclass,$aasize,\%perfect_dups,\%perfect_frags,$smallaa_dropids);  # add info from %perfect_dups , %perfect_frags
  close($OUTB) if($oksumt);
  
  my @outfiles= asmdupclass_fileset($outclass,$cdnaseq,$aaseq,$cdsseq,\%perfect_dups,\%perfect_frags,$smallaa_dropids); 
  loggit(0, "asmdupfilter_fileset=", @outfiles);

  ## turn off tidyup unless( -s $outclass);
  if($tidyup and -s $outclass) {
    loggit(0, "tidyup output folders: okayset dropset inputset tmpfiles");
    ## tidyup needs tod/basename($fn) ?? use perl:rename($fn,"$tod/$tf") ? log nfiles moved not each filename?
    sub tidyup{ my($tod,@td)= @_; mkdir($tod); foreach my $fn (@td) { my $tf=basename($fn); runcmd("mv $fn $tod/$tf") if(-f $fn); } }
    tidyup("okayset",@okayset);  
    tidyup("dropset",@dropset);  
    tidyup("inputset",@inputset);  
    tidyup("tmpfiles",@tmpfiles);  
    my @rmlist;
    foreach my $fn (@erasefiles) { 
      if(-f $fn) { unlink($fn); push @rmlist,$fn; } 
      elsif(-d $fn and $fn =~ m/split/) { rmdir($fn); push @rmlist,$fn; } # only tidy splits now.. leave other subdirs
      ## elsif(-d $fn and $fn ne "." and $fn !~ m,^/,) { rmdir($fn); push @rmlist,$fn; }
      } 
    if(@rmlist) { my $nrm=@rmlist; my $rml=join" ",@rmlist[0..4]; loggit(0,"tidyup erase: n=$nrm, $rml .."); } 
  }


=item step 7, tr2aacds4_stage2b, 2nd stage filter
  
  from sra2genes SCRIPT_evgpubset_UPD1911()
  call AFTER tidyup moves file sets, using okayset/ tmpfiles/,
  creates publicset/  want that renamed? as okayset_stage2?  orig okayset_stage1?
  
  require("tr2aacds4_stg2a.pl");
  @res= tr2aacds4_stage2b( $trclass, $cdsblast, TIDY2STAGE);
  
=cut

  if(UPD1908) {  # add user option for tr2aacds4_stage2b ?
    unless($steprun{7} eq "no") {    
    my $cdb=($tidyup) ? "tmpfiles/$cdsblast" : $cdsblast; # tidyup moves cdsblast
    use constant TIDY2STAGE => 1; # debug
    my @stg2res= tr2aacds4_stage2b($outclass, $cdb, TIDY2STAGE);
    }
  }
  
=item Note step 6 is obsolete
  step 6, genome-map quals, is now obsolete (below)
=cut

  loggit(0, "DONE at date=",`date`);
  loggit(0, "======================================"); # log may append runs
}

#-------------------------------------------------

=item sub tr2aacds4_stage2b, 2nd stage filter

  code from evigene/scripts/prot/tr2aacds4_stg2a.pl
  
=cut

sub tr2aacds4_stage2b {
  my($trclass,$cdsblast,$tidy2stg)= @_;
  $tidy2stg ||= 0;
  
  #UPD1911: cdsqual,  makeblastscore, blasttrset2exons, trclass2pubset, asmrna_altreclass3c
  $ENV{debug}=$debug; #set for callees
  
  unless ( -s $trclass ) { loggit(1,"ERR: missing -class $trclass"); return(-1); }
  unless ( -s $cdsblast ) { loggit(1,"ERR: missing cdslbast $cdsblast"); return(-1); }

  ## Above at start
  # $ENV{'asmrna_altreclass'}="$EVIGENES/rnaseq/asmrna_altreclass3c.pl"; # for trclass2pubset, fix this how?
  # my %updef= ( pCDSOK => 20, pCDSBAD => 20, ALTFRAG => 0.30, CDSXUTR =>  "900,60");
  # for my $uk (keys %updef) { $ENV{$uk}= $updef{$uk} unless($ENV{$uk}); } # no zeros
  
  # NOTE: assume this may be rerun, with some intermediate data, reuse ??
  # NOTE: this trnrname keeps tmpfiles/ path, dont add below
  
  my $idprefix=$ENV{idprefix} || ""; # fixme, from ORGANISM / species
  # my $ORGANISM=$ENV{ORGANISM}||$ENV{species}||"No_name"; # want for selfref pubids ? short Spenam Genspe Nonam
  $ENV{ORGANISM}= $ORGANISM if($ORGANISM); # pass this way if from option
  $ENV{EVGFLAGS}||=""; $ENV{EVGFLAGS} .="tr2aacds4_stage2b"; # for callees
  
  my $trname= $trclass; $trname=~s/.trclass$//; # .gz?
  #x my $trnrname = makename($cdsblast,"",".self\d+.blastn"); # ugh: tevg9arathapnrcd1x-self98.blastn_out-self100.btall
  (my $trnrname = $cdsblast) =~ s/.self\d+.blastn.*//;
  #     my($tmpd,$trnrcds)= getFileset('tmpfiles',"trnameXXXX.cds"); .. too many .cds variants, use cdsblast name
  # tmpfiles/tevg9arathapnrcd1x.cds    tmpfiles/tevg9arathapnrxu.cds
  # tmpfiles/tevg9arathapnrxucd1.cds   tmpfiles/tevg9arathapnr.cds

  ## FIXME: use found evgapps .. as found above 
  ## $eap="makeblastscore.pl"; if( $addapp= findevigeneapp($eap) ){ $UPD4APPS{$eap}= $addapp; }
  ## $cmd= $UPD4APPS{"makeblastscore.pl"};
  ## $cmd.= " -pIDENTMIN 99.999 -pmin 0.01 -CDSSPAN -showspan=2 -tall -sizes $trnrname.cds.qual $cdsblast > $inbtall"; 
  # if( $addapp= findevigeneapp("rnaseq/asmrna_altreclass4.pl")) { push @UPD4APPS, $addapp; }
  # if( $addapp= findevigeneapp("makeblastscore.pl")) { push @UPD4APPS, $addapp; }
  # if( $addapp= findevigeneapp("prot/cdsqual.sh")) { push @UPD4APPS, $addapp; }
  # if( $addapp= findevigeneapp("genes/blasttrset2exons2.pl")) { push @UPD4APPS, $addapp; }
  # if( $addapp= findevigeneapp("genes/trclass2pubset.pl")) { push @UPD4APPS, $addapp; }

  
  my($cmd,$opts,$err);
  $cmd="env outcds=1 $EVIGENES/prot/cdsqual.sh $trnrname.cds"; #?? err if missing cds
  $err= runcmd($cmd);

  my $inbtall="$trnrname-self100.btall"; 
  if ( -s $inbtall ) {
    $cmd=""; $err= 0; loggit(0, "Reusing $inbtall" );
  } else {
    $cmd="$EVIGENES/makeblastscore.pl -pIDENTMIN 99.999 -pmin 0.01 -CDSSPAN -showspan=2 -tall -sizes $trnrname.cds.qual $cdsblast > $inbtall";
    $err= runcmd($cmd);
    #? return ($err) if($err); # loggit(LOG_DIE, "Fail $cmd, err:",$err) if ($err);
  }

  #NOTE: need trclass2pubset.pl in 2 steps here, 1st is env norealt=1 to skip altreclass, 2nd after blasttrset2exons
  $opts="-onlypub -norealt -noaltdrops -log -debug"; # make -norealt an option in trclass2pubset
  $opts .=" -idpre $idprefix" if($idprefix);
  $cmd="$EVIGENES/genes/trclass2pubset.pl $opts -class $trclass"; # was env norealt=1 ..
  my $pubset1cmd= $cmd; #? rerun after trclass_resolve_strandmix
  $err= runcmd($cmd);
  return ($err) if($err); # loggit(LOG_DIE, "Fail $cmd, err:",$err) if ($err);
  
  # assume here pubids is always made new, so always recalc exontab
  my $pubids="publicset/$trname.pubids"; # made by trclass2pubset

  ## UPD20jan: trclass_resolve_strandmix.pl .. after 1st trclass2pubset ? before blasttrset2exons ?
  ## *** redo cdsqual, intball for cds and strand changes **
  ## .. resolution method: use pubids gene groups, for each group 2+ trs, check trclass for -sense flag,
  ## .. if -sense flag, check all aaqual codepot, pick one sense for largest codepot, with caveats? 
  ## .. quals aaref* if have it, aasize/complete/aaqual? aacons? no
  ## .. then swap revaa for -sense trs, changing okayset seqs
  #above use constant UPD20j_resolve_revaa => 0; # also $reorient
  
  if(UPD20j_resolve_revaa and $reorient) {
    $opts="-calccodepot -fragcheck -log -out"; #! fragcheck needs ncbibin/blastn, not essential
    #UPD20apr: -nomaybe user opt, no fwd+rev ambiguous, pick 1.
    #   : auto-detect from aaqual, dont need user opt?
    $opts .= " -strandedrna=$RNAIsStranded" if($RNAIsStranded); # UPD22a, only if RNAIs=1|0|yes|no
    $opts .= " -nomaybe" if($ENV{reor_nomaybe}); # nomaybeor ?
    # should become defaults: -calccodepot X -changesonly 
    $cmd="$EVIGENES/genes/trclass_resolve_strandmix.pl $opts -pubids $pubids -class $trclass"; 
    $err= runcmd($cmd); 
    $ENV{DID_RESOLVESENSE}=1; # *** flag for altreclass ***
      ## collect info from trname.reorient.log
        #reor: restrand_seqs n=615
        #reor: rewriteTrClass human18ncx.trclass to human18ncx.trclass.upreor,  nup= 847
    my $nreort = -1; my $reorlog="$trname.reorient.log"; 
    if(-f $reorlog) { 
      my $reort=`tail $reorlog | egrep 'restrand_seqs|nup='`; 
      ($nreort)= ($reort =~ m/restrand_seqs n=(\d+)/) ? $1 : -1;
      loggit(0, $reort);
      }      
    if($nreort and not $err) {
      # trclass_resolve_strandmix makes new trclass, adds okayset/name.okreor.aa,cds,idtab
      if( -d 'publicset'){ rename('publicset','publicset_orold'); } # push @erasefiles, 'publicset_orold'; #?
      $err= runcmd($pubset1cmd); # Need to rerun trclass2pubset -onlypub 
    }
  }

  if ( -s $inbtall and -s $pubids ) {
    $cmd= "sort -k7,7nr -k2,2 -k6,6nr -k1,1 $inbtall | env pubids=$pubids debug=1 $EVIGENES/genes/blasttrset2exons2.pl > $trnrname.exontab";
    $err= runcmd($cmd); # return ($err) if($err); 
  }
  
  ## final step, integrate into trclass2pubset; see trclass2pubset:altreclass_block()
  ## add -NOALTDROPS, same as -nodrops for  asmrna_altreclass3c, this turns drops to culls, separate files

  $opts="-noaltdrops -exontab $trnrname.exontab  -log -debug ";
  $opts .=" -idpre $idprefix" if($idprefix);
  $cmd="$EVIGENES/genes/trclass2pubset.pl $opts -class $trclass";
  $err= runcmd($cmd);
  return ($err) if($err); # loggit(LOG_DIE, "Fail $cmd, err:",$err) if ($err);

  # maybe here, rename okayset/ to okayset_stage1/ ?
  #  rename publicset/  to okayset/ or okayset_stage2/ ?   
  #  for compat, symlink publicset/name.{aa,cds,mrna}_pub.fa to name.okay.{aa,cds,mrna}, skip okalt.
  #  for clarity, symlink name.pubids to name.genetable.tsv or .txt
  #  recreate name.trclass.sum.txt table of classes, for stage2? or merge into genesummary.txt ?

  if($tidy2stg) {
    loggit(0,"tidy okayset => okayset1st, stage2 reduction => okayset");
    rename("okayset","okayset1st");
    my($pubd,@pubseq);
    
    ## FIXME: prepubset $trname.{mrna,cds,aa} need to move .. these are eq okayset1st/*.tr,cds,aa ; erase?
    ($pubd,@pubseq)= getFileset('publicset','mrna$|cds$|aa$',$pubd,$trname); 
    for my $ps (grep { not /_pub|_cull/ }  @pubseq) {
      (my $tps=$ps) =~ s/publicset/tmpfiles/;
      rename($ps,$tps);
    }
    
    ($pubd,@pubseq)= getFileset('publicset','fa',$pubd,$trname);  
    for my $ps (grep/_pub.fa|_cull.fa/, @pubseq) {
      my ($ok,$oks)= (0,$ps);
      if($oks =~ s/(\w+)_pub.fa/okay.$1/ ) { $ok=1; }
      elsif($oks =~ s/(\w+)_cull.fa/cull.$1/ ) { $ok=1; }
      rename($ps, $oks) if($ok); # symlink($ps,$oks); #?? or rename them, and culls also      
    }

    rename("publicset","okayset"); #  reserve pubset for later named annots
  }
  
  return(0); #??
}


=item OBSOLETE step6 genome-map reclassing : see SRA2Genes

use constant GMAP_OPTION6 => 0;
if(GMAP_OPTION6) {
  ## UPD 2016.02 : add option step 6, ?? after tidyup moves to okayset
  # 6. genome-map reclassing, *after* 1st outclass using okayset
  # make separate pipeline * use after evgmrna2tsa2.pl publicset, for proper alt matching
  my $chrasm=0;
  sub blastn_cds_genome{};
  sub cdsgmap_maketables{};
  if($chrasm) {
    # 6.1. blastn okayset/my.cds to chrasm
    my($cdsgmapblast)= blastn_cds_genome($chrasm, @okayset);
    # 6.2. tabulate blast > cdsgmap.tall > cdsgmap.equalgene
    my @xx= cdsgmap_maketables($cdsgmapblast);
    # 6.4. classify main/alternate cds, okay & drop subsets, using evigene/rnaseq/asmrna_dupfilter2b.pl
    # combine both cdsblastab + cdsgmapblast.eqgene sources of overlap classing
    ($outclass,$outaln)= asmdupfilter_cds($cdsgmapblast,$aasize,$aaclstr);
    # 6.5: new asmdupclass_sum(), asmdupclass_fileset(), ..
  }
}

=cut  
  

#-------------------------------------

sub nonredundant_cds
{
  my($cdsseq)=@_;

  my $cdsnrseq = makename($cdsseq,"nr.cds");
  my $cmd="$APPfastanrdb $cdsseq > $cdsnrseq";
  # should be: fastanrdb -i -f $cdsseq ; -i == --ignore-case
  unless(-s $cdsnrseq) {
  my $runerr= runcmd($cmd);
  }
  push @tmpfiles, $cdsnrseq;
  return($cdsnrseq);
}


  ## 1.1. reassign redundant cds from aaqual = pCDS to reduce utrbad/utrpoor set.
  ## .. need only change header ID in $cdsseqnr ; fastanrdb puts all redundant ids on header.

# %perfect_dups= redundant_idset($cdsseqnr); # read headers after rebest, for dupclass_fileset
sub redundant_idset
{
  my($cdsnrseq)=@_;
  my %pdup=();
  open(F,$cdsnrseq) or return %pdup; 
  while(<F>) { if(/^>/) { s/>//; 
    my($d1,@dp)=split; 
    if(AACONS) { @dp= grep{ not m/agree=/ } @dp; }
    map{ $pdup{$_}=$d1; }@dp; } } close(F);
  return %pdup;
}

#  my($nagree,$agreeset) = consensus_idset($cdsseqnr); # read headers after rebest, for dupclass_fileset
sub consensus_idset # (AACONS)
{
  my($cdsnrseq,$flag)=@_;
  my %agree=(); my $nagree=0;
  my $retall=($flag=~/only/)?0:1; # agreeonly
  open(F,$cdsnrseq) or return 0;  # my ($ok,$hin)= openRead($cdsnrseq);
  while(<F>) { if(/^>/) { s/>//; 
    my($d1,@dp)=split; 
    # my($ag)= grep /agree=/, @dp;
    my($ag)= (m/agree=(\w+)/)?$1:0;    
    $agree{$d1}= $ag if($ag or $retall); 
    $nagree++ if($ag);
    # @dp= grep{ not m/agree=/ } @dp; }
    } 
  } close(F);
  return ($nagree, \%agree);
}

sub add_consensus_idset
{
  my($agreefile, $agreelist)=@_; 
  my($ncons,$consfile)=(0,"none");
  $agreelist={} unless(ref $agreelist);

  if($agreefile and -f $agreefile) {
    open(F,$agreefile);  
    while(<F>) { 
      next unless(/^\w\S+\s+\d/);  
      my($id,$agscore,$agsrc,$agids)= split;
      if($agscore > 1) {
        $agreelist->{$id} += $agscore; $ncons++;
        if($agids) { map{ $agreelist->{$_} += $agscore; } split(",",$agids); }
      }
    } close(F);
  }  
  loggit(0,"add_consensus_idset n=$ncons from $agreefile"); 
  return($ncons,$agreelist);
}

sub homology_idset # UPD1908: collect trids from aablast table: trid <tab> name or refid <tab> ...
{
  my($ablastab,$flag)=@_;
  my %hashomol=(); my $nhasho=0;
  open(F,$ablastab) or return 0;  # my ($ok,$hin)= openRead($cdsnrseq);
  while(<F>) { if(/^\w/) { 
    my($d1,@dp)=split; # ignore other cols, score
    $hashomol{$d1}= 1; 
    $nhasho++;
    } 
  } close(F);
  return ($nhasho, \%hashomol);
}

sub consensusof { 
  my($nrtest, $dmain, @d2)= @_;
  my($cscore,$cmain,$ctop,$cset)= (0) x 9;
  $cmain= $dmain;
if($USE_AACONS) {

  # assume for now trasm IDs are from evg trformat.pl with common patterns:
  #  subsetprefix | asmmethod | kKMER | [Ll]oc NNNN tIII where locusNNNN alternate tIII are chopped off
  # FIXME: need to use $aconsensus param as options for how to parse ids
  # eg: acons='(trin|velv|idba|soap|pacb)' is source pattern
  # or: acons='(trin|idba|pacb|(?:velv|soap)\w*k\d+)' for kmer setting
  # or: acons='^(\w+)' .. etc
  
  my $userparse= ($aconsensus =~ m/\w/)?1:0;
  # my $nokmer= 1; 
  my $nokmer= (($aconsensus =~ m/kmer/) or $ENV{'conkmer'})?0:1; #UPD19
  my @dsrc=();
  for my $d ($dmain,@d2) {
    my($s,$sdef)=("",$d);
    if(1) { # default parse
      unless($sdef=~s/[Ll]oc.*$//) { # need checks, if no LocNNN chomp numbers
        $sdef=~s/utrorf//; $sdef=~s/_//g;
        while($sdef=~s/(\d)\D\d+$/$1/) { ; } # AT2G1234.1 > AT2G1234; XM_1234.2 > XM_1234; TrinDN1234_c0_g2_i1 TrinDN1234c1g2t1  > TrinDN; 
        $sdef=~s/\d+$//; 
        }
      $sdef=~s/k\d+$// if($nokmer or $sdef=~/idba/); # fixme for asmbler tags
      $s= $sdef;
    } 
    if($userparse) { ($s)= ($d=~m/($aconsensus)/)?$1:$sdef; }
    push @dsrc, $s if($s=~/\w/);
  }
    
  my %dsrc=();  map{ $dsrc{$_}++ } @dsrc;
  @dsrc= sort{ $dsrc{$b}<=>$dsrc{$a} or $a cmp $b } keys %dsrc; # dont care for sort by count?
  if(@dsrc>1) {
    $cscore= @dsrc; # number of sources, for 1st test
    $ctop= $dsrc[0]; 
    $cset=join",", sort @dsrc; # alpha sort for display
    $cmain= $dmain; # for now
    #x ($cmain)= grep /$ctop/, ($dmain,@d2); # dont want this, or needs testing
  }
}  
  return($cscore,$cmain,$ctop,$cset);
}


## $nbetter += nrcheckaaqual('nrdup|nrfrag',$aaqualh,\%better,$mainid,@altids) : sets hash %better and $nbetter
## revise nrcheckd to use same sub w/ nofrag_reassignbest 
use constant fPCDS_DIFF => 1; # what? should change any small diff?
use constant dPCDS_DIFF => 5; # what? should change any small diff?

## $nbetter += nrcheckaaqual('nrdup', $aaqualh, \%better, @dp);} 
## $nbetter += nrcheckaaqual('nrfrag', $aaqualh, \%better, $mainid, @fragids); 
sub nrcheckaaqual 
{
  my($nrtest, $aaqualh, $better, $dmain, @d2)=@_;  
  my $dofrag= ($nrtest =~ /frag/)?1:0;
  my $aqmain= $aaqualh->{$dmain} or return 0; # warn ERR ?
  my($awm,$apm,$acm,$aqm)=split",",$aqmain;  
  our($pCDSpoor,$pCDSbad); # in cdna_proteins.pm

  my ($isbetter,$okmain,$awmMIN, $amMIN, $amMINd2 )= (0) x 9;
if(UPD1908) {
  $okmain=($acm < kAAQUAL_MAX-1  or $dmain =~ /utrorf/)?0:1; ## or $aqmain=~/utrbad|utrpoor|partial/)?0:1;
  $amMIN  = 0.97; $awmMIN= 0.97 * $awm;  # drop awmMIN, awm changes in loop
  $amMINd2= 0.93; # for larger aaqual diff
  #UPD1908: $acm < kAAQUAL_MAX  too stringent to kick out longest orf, allow utrpoor at least
  #UPD1908: also bump 0.90 to 0.95 for shortest frag replacemtn
  #UPD1908: ** need to respect cdna_bestorf:ORF_FULLvPART option with awmMIN and acmain/acfrag =~ partial
} else {
  $okmain=($acm < kAAQUAL_MAX  or $dmain =~ /utrorf/)?0:1; ## or $aqmain=~/utrbad|utrpoor|partial/)?0:1;
  $awmMIN= 0.90 * $awm; # orig, too low/short
}
  
  return 0 if($dofrag and $okmain);

  my $idbest= $dmain; 
  ## sort @d2 by best aaqual ? using aaqualh acm scores ?
  my($medianutr)=(0,0);
  unless($dofrag) { 
    my @pcds= sort{ $b <=> $a } map{ my($aw,$ap)= split",", $aaqualh->{$_}; $ap||0; } @d2;
    my $n2= int( @pcds / 2); 
    $medianutr= $pcds[$n2]; # $medianid = $d2[$n2]; #  id may have poor qual.. dont use
  }
  foreach my $d2 (@d2) { 
    my $aqfrag= $aaqualh->{$d2} or next;
    my($awf,$apf,$acf,$aqf)=split",", $aqfrag; # BUG: $aqmain; 
    #o my $pdif= $apf - $apm; # ap[fm] is percent coding, ie. utrok/poor/bad score .. shouldn't weight too much
     
    if($dofrag) {
      my $okfrag= ($acf < kAAQUAL_MAX or $d2 =~ /utrorf/) ? 0 : 1;
if(UPD1908) {
      ### $acf has CDSpoor value : $okfrag= 0 if($apf < $pCDSpoor); #UPD1908
      if($okfrag) { # drop: $pdif >= fPCDS_DIFF
        my $cdif = $acf - $acm; # aaqualscore() for complete/utrbad/gapbad, range is -3..+3,
        my $fok= ($cdif >= 2) ? ($awf >= $awm * $amMINd2) : ($awf >= $awm * $amMIN);
        if($fok) { $idbest= $d2;  $isbetter=1;
          ($awm,$apm,$acm,$aqm)=($awf,$apf,$acf,$aqf);
        }
      }
} else {     
      if($okfrag and $awf >= $awmMIN ) {  # drop: and $pdif >= fPCDS_DIFF
        $idbest= $d2; $apm= $apf;  $isbetter=1; 
      }
}      
    } else { 
      ## cds perfect-duplicate, so aaquals are same, only utrs differ
      #o my $okdup= ($d2 =~ /utrorf/)?0:1; #upd1905
      
      ##UPD1908 use of dPCDS_DIFF filter problem: selects *shortest* UTR for same CDS,
      ## instead want median? or longest < aberrant? ie pCDS in 60% to 90% range, but depends on other qual
      #o my $pdif= $apf - $apm; # ap[fm] is percent coding, ie. utrok/poor/bad score .. shouldn't weight too much

      my $okdup= ($acf < $acm or $d2 =~ /utrorf/)?0:1; #upd1905
      $okdup= 0 unless($okdup and abs($apf - $medianutr) < abs($apm - $medianutr)); #?? pick median utr size

      if( $okdup ) { # and $pdif >= dPCDS_DIFF
        $idbest= $d2; $isbetter=1; 
        ($awm,$apm,$acm,$aqm)=($awf,$apf,$acf,$aqf);
        }  
    }
  }
  
  my $dobett= ($idbest ne $dmain)?1:0;
  my $addflag= "";
  if($USE_AACONS) {  # add header flag in better->{dmain}= newhdr
    if($nrtest=~m/,(agree=\w+)/) { $addflag=$1; $dobett=2; }
     ##   $nrflag= "$nrdup,agree=$cscore"; #?
  }
  
  if($dobett) { ## ($idbest ne $dmain) ; # if($isbetter)
    $isbetter=1; 
    # if(AACONS and $dofrag) .. DONT delete dmain if it has agree= flag ??
    if($dofrag) { 
      if($idbest ne $dmain) {
        delete $better->{$dmain} unless($isbetter); # UPD1908: dont delete longer, if shorter rescued only by agree= val
        }
      $better->{$idbest}= $dmain;
      }
    else { my @da= grep { $_ ne $idbest } ($dmain,@d2);  $better->{$dmain}= join(" ",$idbest,$addflag,@da);  }
    #o else { @d2= grep { $_ ne $idbest } @d2;  $better->{$dmain}= join(" ",$idbest,$dmain,@d2);  }
  }

  return $isbetter;
}

sub nonredundant_reassignbest
{
  my($cdsnrseq, $aaqual, $agreeset)=@_;
  my ( %better, %aconsensus ); # %aq,
  my ($nbetter,$nrec,$nrtot,$naddcons)=(0) x 9;
  my $flagfile= "$cdsnrseq.isbest";
  unless( -s $cdsnrseq and -s $aaqual ) {
    loggit(1,"ERR: nonredundant_reassignbest missing cdsnr:$cdsnrseq or aaqual:$aaqual"); 
    return($nbetter,$nrec);
  }
  return($nbetter,$nrec) if( -f $flagfile or $dryrun);
  
#   my $cdsnrhdr = makename($cdsnrseq,".hdr");
#   my $cmd="grep '^>' $cdsnrseq > $cdsnrhdr"; # dont really need .hdr file; combine below
  
  runcmd("touch $flagfile"); push @tmpfiles, $flagfile; #??
  
  our $aaqualh= getAaQual($aaqual); # $qualh # may be global hash $AAQUALH, read once  
  # open(F,$aaqual); while(<F>) { my($id,$aw,$gp,$aq,$tw)=split; $aq{$id}="$aw,$aq" if($aq); } close(F);
  
  #   # sub nrcheckaaqual == revised nrcheckd to use same sub w/ nofrag_reassignbest 
  
  my($ok,$hin)= openRead($cdsnrseq); # open(F,$cdsnrseq); 
  #o# while(<$hin>) { if(/^>/){ s/>//; my @dp=split; nrcheckd(@dp) if(@dp>1); } } close($hin);
  while(<$hin>) { if(/^>/) { 
    s/>//; my @dp=split; $nrtot++;
    
    if(@dp>1) {
      my $nrflag= 'nrdup';
      
      if($USE_AACONS) {
      my($cscore,$cmain,$ctop,$cset)= consensusof('nrdupcds',@dp); # if(AACONS)
      if(ref $agreeset){ map{ $cscore += $agreeset->{$_}||0; } @dp;  } #UPD1912
      #? alter aaqualh->{$dp[0]} with cscore ? 
      # * add consensus flag in cdsnrseq header, dont change aaqual better hdr order tho
      if($cscore) {
        $nrflag= "nrdup,agree=$cscore"; #?
        my $cids= join",",@dp; # YES, need these as nofrag rewrites nr.cds headers
        $aconsensus{ $dp[0] }= join"\t",$cscore,$cset,$cids;  
        }
      }
      
      $nrec++; 
      $nbetter += nrcheckaaqual($nrflag, $aaqualh, \%better, @dp);
    } else {
      if($USE_AACONS and ref $agreeset) { #UPD1912
        my $cids= $dp[0];
        my $cscore = $agreeset->{$cids}||0;   # FIXME: cset == prefix of cid .. from consensusof()??
        if($cscore) { $aconsensus{ $dp[0] }= join"\t",$cscore,$cids,$cids; $naddcons++; }
        ## ... this is a problem to fix, below: 
        #   my $maxcon= int($nrtot * 0.40);  # NOT $nrec; dunno how many to allow
        #t2ac: add_consensus_idset n=1051546 from ./evg4567hetfix.agree
        #t2ac: too many cds have consensus to be usable: 813627 of 1555374 (52%); skipping consensus..
      }
    }
    } 
  } close($hin);

  # AACONS fixme here? write table trset.consensus for asmrna_dupfilter3c.pl
  if(AACONS and $USE_AACONS) {
    my $maxcon= int($nrtot * 0.40);  # NOT $nrec; dunno how many to allow
    my @cids= sort keys %aconsensus;
    my $ncons= @cids; my $nconstest= $ncons;
    if($naddcons) { #UPD1912 maxcon for addcons/agreeset
      # $maxcon += $naddcons; #? or pct of naddcons? # or $ncons - $naddcons?
      $nconstest= $ncons - $naddcons;
    }
    
    my $cname= makename($cdnaseq,".consensus"); #? global? need below for asmrna_dupfilter3c.pl
    #t2ac: consensus transcripts n=87908 tabled in evg3weed1rc.consensus
    #^ FIXME: avoid rewriting prior consensus .. need opt?
    
=item too many cons bug, nrec wrong base

#t2ac: CMD= /home/ux455375/bio/exonerate/bin/fastanrdb evg3weed1rc.cds > evg3weed1rcnr.cds
#t2ac: nonredundant_cds== evg3weed1rcnr.cds nrec= 767122
#t2ac: CMD= touch evg3weed1rcnr.cds.isbest
#t2ac: getAaQual: naa=2016751 in evg3weed1rc.aa.qual, val1 arab6roo3dn4msoapk27loc170215t1= 509,85,2,partial3
#t2ac: too many tr marked consensus to be usable: 87908 of 246815 ; skipping consensus..
        nrec=246815 is wrong total, 35.6% of those w/ dups, vs 11.4% of all nr.cds
#t2ac: nonredundant_reassignbest= 130619 of 246815
    # grep '^>' evg3weed1rcnr.cds | wc -l = 767122
    # grep '^>' evg3weed1rcnr.cds | grep agree= | wc = 87908  
    
=cut

    if($nconstest > $maxcon) { 
      my $pcons= ($nrtot>0)?int(100*$ncons/$nrtot):0;
      loggit(1,"too many cds have consensus to be usable: $ncons of $nrtot ($pcons%); skipping consensus..");  

    } else {
      #* put list of all nr ids into trset.consensus? col3, Yes, need these later
      if(-s $cname) { rename($cname, "$cname.old"); }
      my $ok= open(B,">$cname");
      for my $cid (@cids) { print B "$cid\t",$aconsensus{$cid},"\n"; }
      close(B);
      # push @tmpfiles, $cname; #? or leave in main folder
      loggit(0,"consensus transcripts n=$ncons tabled in $cname"); 
    }
  } 
  
  if($nbetter>0) { # rewrite $cdsnrseq headers 
    my $ok= open(B,">$cdsnrseq.best");
    if($ok) { 
      ($ok,$hin)= openRead($cdsnrseq); 
      while(<$hin>) { if(/^>(\S+)/) { if(my $hdr= $better{$1}) { s/>.*$/>$hdr/; } } print B $_; } 
      close(B); close($hin); 
      my $cmd="mv $cdsnrseq $cdsnrseq.old; mv $cdsnrseq.best $cdsnrseq";
      system($cmd); #? runcmd($cmd);
      push @tmpfiles, "$cdsnrseq.old";
    }
  }
   
  return($nbetter,$nrec); # ,$ncons what?
}

=item 1.1 reassign redundant cds using aaqual  

  grep '^>' shrimt1trin1nr.cds > shrimt1trin1nr.cds.hdr
  grep '  ' shrimt1trin1nr.cds.hdr | cat shrimt1trin1.aa.qual - | perl -ne \
'if(/^>/) { s/>//; my @dp=split; checkd(@dp) if(@dp>1); } elsif(/^(\w+)\t/) { ($id,$aw,$gp,$aq,$tw)=split; 
$aq{$id}="$aw,$aq"; } END{ warn "done ncheck=$ncheck nbad=$nbad nbetter=$nbet\n"; } 
sub checkd{ my($d1,@d2)=@_; ($aw1,$aww1,$pc1,$aq1)=split",",$aq{$d1}; $ncheck++; $nbad++ if($pc1<1); 
foreach $d (@d2) { ($aw,$aww,$pc,$aq)=split",",$aq{$d}; $nbad++ if($pc<1); $pdif=$pc - $pc1; 
if($pdif>9) {  print "$d:$aw,$pc,$aq\tbetter\t$d1:$aw1,$pc1,$aq1\n"; $nbet++; } }}' 

done ncheck=15630 nbad= nbetter=8389  # all diff
done ncheck=15630 nbad= nbetter=2292  # pdif > 9

  must changes:  good vs utrbad
shrimt1trin1loc21564c0t5:1385,81%,complete      better  shrimt1trin1loc21564c0t2:1385,47%,complete-utrbad
shrimt1trin1loc9285c0t1:389,65%,complete        better  shrimt1trin1loc8886c0t1:389,33%,complete-utrbad
shrimt1trin1loc468087c0t1:42,50%,complete-utrpoor       better  shrimt1trin1loc229219c0t1:42,31%,complete-utrbad
shrimt1trin1loc43780c0t2:95,54%,complete-utrpoor        better  shrimt1trin1loc43780c0t1:95,37%,complete-utrbad

  maybe changes: same utr-qual but minor pCDS improvement
shrimt1trin1loc21426c0t4:555,83%,complete       better  shrimt1trin1loc21426c0t3:555,66%,complete
  = 15% pcds improve
shrimt1trin1loc65367c0t1:55,30%,complete-utrbad better  shrimt1trin1loc36046c0t1:55,7%,complete-utrbad
  = 23% pcds improve -- still bad but 30% >> 7% cds  
  
  ignore small changes?
shrimt1trin1loc58089c0t1:523,73%,complete       better  shrimt1trin1loc45529c0t1:523,70%,complete
  
=cut
 
=item smallclassfilter

  UPD1908:  new smallclass reduction step b/n  nonredundant_cds() and nofragments_cds() 
   remove SMALLMAXAA >= smallaa >= MINAA unless( has_evidence( aahomol | aaconsensus | other) )

  my($nok, $ncut, $cdsseqnrUPD)= 
     evd_smallclassfilter( $cdsseqnr, $aaqual, $smallclassmax, $agreeset, $homolset);

  SEE ALSO asmrna_altreclass noclasscut (0 def in altreclass, 60 def as aaminnoclass in trclass2pubset.pl caller)
  # caller trclass2pubset.pl uses $ENV{aaminnoclass}||60; for asmrna_altreclass -noclasscut=60

=cut
    
sub evd_smallclassfilter {
  my( $cdsseqnr, $aaqual, $smallclassmax, $nagree, $agreeset, $nhashomol, $homolset)= @_;
  my( $nin, $nok, $ncut, $ninevid)= (0) x 9;

  unless( -s $cdsseqnr and -s $aaqual) {
    loggit(1,"ERR: evd_smallclassfilter missing cdsnr:$cdsseqnr or aaqual:$aaqual"); 
    return(0,0);
  }
  # return($nok, $ncut) if( $dryrun);
  
  my $nevid= $nagree + $nhashomol;
  my($aaqualh,$naa)= getAaQual($aaqual); ## may be global AAQUALH, read once?
  
  if($naa > 20000 and ($nevid < 0.10 * $naa and $nevid < 20000)) {
    loggit(1,"NOTE: evd_smallclassfilter not enough evidence nev=$nevid for naa=$naa"); 
    return(0,0);
  }
  
  my $cdsnrids= faidlist($cdsseqnr,'',"ashash"); 
  my (%cdsfiltok,%smalldropid);
  for my $id (sort keys %$cdsnrids) {
    my $ok=1; $nin++; 
    my $hasevd= ($agreeset->{$id} or $homolset->{$id}); # count to ensure enough to decide?
    $ninevid++ if($hasevd);
    my $aqr= $aaqualh->{$id}; # my($awf,$apf,$acf,$aqf)=split",",$aqr; 
    my($awf)= ($aqr)? (split",",$aqr) : 0;    
    if($awf > 0 and $awf < $smallclassmax) {
      $ok=0 unless($hasevd); # ($agreeset->{$id} or $homolset->{$id});
    }
    if($ok){ $nok++; $cdsfiltok{$id}= $ok; } else { $ncut++; $smalldropid{$id}=$awf||1; }
  }

  return(0,0,$cdsseqnr) unless($nok and $ncut);
  if( $ninevid < 0.10 * $nin) {
    loggit(1,"NOTE: evd_smallclassfilter not enough evidence nev=$ninevid for ncds=$nin"); 
    return(0,0,$cdsseqnr);
  }
  
  my $cdssmalfilt = faextract($cdsseqnr, makename($cdsseqnr,"sm.cds"), \%cdsfiltok, 1); # 1 == no update
  push @tmpfiles, $cdssmalfilt; 
  return($nok,$ncut,$cdssmalfilt,\%smalldropid);
}

  # 2. perfect fragment removal
  # 2.FIXME: some perf-frag maybe/are aabest (complete) vs perf-longest (aapoor), reassignbest cdhit.clstr + aaqual, 
sub nofragments_cds
{
  my($cdsseq)=@_;

  my $cdsnrseq = makename($cdsseq,"cd1.cds");
  my $cdlog = makename($cdsseq,"cd1.log");
  ## *** clusterize this .. -T NCPU works ok
  ## cd-hit-est -M 1500  -c 1.00 -d 0;  -l  length of throw_away_sequences, default 10
  ## .. redo opts so can use ENV replacements ..
  my $opts=" -c 1.00"; 
  $opts.=" -T $NCPU" if($NCPU>1);
  $opts.=" -M $MAXMEM" if($MAXMEM>1); # M in Megabytes
  $opts.=" -l ".($MINCDS-1) if($MINCDS>10); 
  if(my $cdopt=$ENV{CDHITESTOPT}) { $opts .= " ".$cdopt; }
  
## pogonus1all3nr.cds : cd-hit-est failing w/ bizzare mem request; maybe cds > 1.2G file too big?
## opt here to fasplit(2+), then cd-hit-est-2d -i pt1.cds -i2 pt2.cds, combine parts..

  # my $cdsnrids = makename($cdsnrseq,".ids"); # idfile-only or hash?
  # my $idops="noupdate";  
  my $cdsnrids = {}; # idfile-only or hash?
  my $idops="ashash"; # ashash <=> noupdate for idfile-only
  my $cmd="$APPcdhitest $opts -d 0 -i $cdsseq -o $cdsnrseq 1> $cdlog 2>&1"; # log to ...
  unless(-s $cdsnrseq) {
    my $runerr= runcmd($cmd); # $idops="update";
    unlink("$cdsnrseq.bak.clstr") if (-f "$cdsnrseq.bak.clstr");
    # push @erasefiles,"$cdsnrseq.bak.clstr";

    my($nbetter,$nrec)= nofrag_reassignbest("$cdsnrseq.clstr", $cdsnrseq, $cdsseq,  $aasize);
    loggit(0,"nofrag_reassignbest=",$nbetter,"of",$nrec,$cdsnrseq); 
    # $cdsnrseq= $cdsnrbest if($cdsnrbest); ## sub renames new best to old cdsnrseq .. dont need this
  }
  
  #file# $cdsnrids= faidlist($cdsnrseq,$cdsnrids,$idops); # returns old if exists..
  $cdsnrids= faidlist($cdsnrseq,'',"ashash"); 
  # FIXME: grep '^>' $cdsnrseq > $cdsnrseq.hdr for later ref
  push @tmpfiles, $cdsnrseq, "$cdsnrseq.clstr", $cdlog; #, $cdsnrids ;
  return($cdsnrseq,$cdsnrids);
}

sub nofrag_reassignbest
{
  my( $cdhitclstr, $cdsnrseq, $cdsinseq, $aaqual)=@_;
  my ( %better ); # %aq, 
  my ($nbetter,$nrec, $ntest, $nerr)=(0) x 9;
  my $flagfile= "$cdsnrseq.isbest";
  unless( -s $cdsnrseq and -s $aaqual and -s $cdhitclstr) {
    loggit(1,"ERR: nofrag_reassignbest missing cdsnr:$cdsnrseq, $cdhitclstr or aaqual:$aaqual"); 
    return($nbetter,$nrec);
  }
  return($nbetter,$nrec) if( -f $flagfile or $dryrun);

  ## 20140525: correct from evigene/scripts/prot/testnrcd1best.pl  

  runcmd("touch $flagfile"); push @tmpfiles, $flagfile; #??
  
  my $aaqualh= getAaQual($aaqual); # ret: ($qualh,$sizeh,$naa) : $qualh # may be global hash, read once?
  my %aqv; ## for sort only ?
  for my $id (keys %$aaqualh) { 
    my $aqr= $aaqualh->{$id}; my($awf,$apf,$acf,$aqf)=split",",$aqr; 
    $aqv{$id}=$awf * $apf * $acf; 
    }
  
  my(@fragids,$mainid,$mainlen);  # %fragids, > better
  open(F, $cdhitclstr); ## or return($nbetter,$nrec);
  while(<F>) { # cd-hit clstr format
  
    if(/^>/) { # new cluster, decide if last mainid is good, or replace by 2nd best ?
      $nrec++; 
      ## can this be replaced w/ shared sub nrcheckd() from above ??
      ## nrcheckd($mainid,@altids) : sets hash %better and $nbetter

      if($mainid and @fragids) { $nbetter += nrcheckaaqual('nrfrag', $aaqualh, \%better, $mainid, @fragids); } 
           
      @fragids=(); $mainid=0; 
       
    } elsif(/^(\d+)/) { my $i=$1;
      m/(\d+)(nt|aa), >(.+)\.\.\. (.+)$/;  
      my($tlen,$typ,$tid,$pinfo)=($1,$2,$3,$4); 
      my $ismain=($pinfo =~ /\*/)?1:0;  # FIXME: drop /$/; new merge fnum at end of line now;
      unless($tid) {
      } elsif($ismain) {
        $mainid= $tid; $mainlen=$tlen; # $pi=100;  
        $better{$mainid}= $mainid; # always set
      } else {
        push @fragids, $tid; #save tlen in hash for sort?; need to collect all, may not have mainlen yet.
      }
    } 
  } close(F);
  if($mainid and @fragids) { $nbetter += nrcheckaaqual('nrfrag', $aaqualh, \%better, $mainid, @fragids); $mainid=0; } 

  if(0 and AACONS) {
    #?? dont need this 2nd read thru $cdsinseq . do below
    my($nagree,$agreeset) = consensus_idset($cdsinseq,'agreeonly'); # read headers after rebest, for dupclass_fileset
    if($nagree>0) {
      my @ida= sort keys %$agreeset; # returns all ids? or just with consensus?
      for my $id (@ida) { $better{$id} = $agreeset->{$id} unless($better{$id}); $nbetter++; }
    }
  }
     
   ## FIXME: maybe dont replace main w/ shorter/aabetter frag, but add aabetter to cdsnrseq ?
   ## NOTE: %better must be prefilled w/ orig cdsnrseq IDs to replace.
  if($nbetter>0) { # replace cds longest seq w/ better frag seq
    my($ok,$nok,$hin)= (0,0); ## check nok == nrec ??
    ($ok,$hin)= openRead($cdsinseq); # open(INF,$cdsinseq)
    $ok= open(OUTF,">$cdsnrseq.best") if($ok);
    if($ok) { $ok=0;
      while(<$hin>) { 
        if(/^>(\S+)/) { 
          my $id=$1; $ok= $better{$id}||""; 
          $nok++; 
          if($USE_AACONS) {
            my($ag)= (m/agree=(\w+)/)?$1:0;    
            # $ok.=" agree=$ag" if($ag); #old# $ok++ if($ag);
            $ok.=" agree=$ag" if($ag and $ok ne $id); #old# $ok++ if($ag);
          }
          s/$/ replacecdtop=$ok/ if($ok and $ok ne $id); # annotate change
        } 
        print OUTF $_ if($ok);
      }  close(OUTF); close($hin); 
      
      if($nok < $nrec) { # ERROR?? skip rename
        loggit(1,"ERR: nofrag_reassignbest fewer than input: reassign=$nok, orig=$nrec for $cdsnrseq.best"); 
        push @tmpfiles, "$cdsnrseq.best";
      } else {
        my $cmd="mv $cdsnrseq $cdsnrseq.old; mv $cdsnrseq.best $cdsnrseq";
        system($cmd);  
        push @tmpfiles, "$cdsnrseq.old";
      }
    }
  }
   
  return($nbetter,$nrec); # what?
}

=item FIXME.201405 new reassignbest_aaqual ? or not

  FIXME.201405 new reassignbest replace old 2: nr_reass and nfrag_reass
  nofrag_reassignbest() fails, nr_reassignbest() not always best..
  #   my($nbest,$ndups)= reassignbest_aaqual($cdsseqnrcd1,$aasize);

  mainaltbetter.info
  -- this does better job at replacing nr and nofrag results using aaqual (aa-complete, no utrbad/poor)
  #nbetter: nmain=98921; mainok=1551; altbet=9728; noaltbet=87642; 

=cut

# sub reassignbest_aaqual  # DROP FOR NOW, revise other 2 ..


# %perfect_frags= fragment_idset("$cdsseqnrcd1.clstr"); # read headers after rebest, for dupclass_fileset
# FIXME: fragment_idset nrcheckaaqual better changes clstr main ids..
sub fragment_idset
{
  my($cdhitclstr)=@_;
  my(%fragids,@fragids,$mainid);
  open(F, $cdhitclstr) or return %fragids;
  while(<F>) { # cd-hit clstr format
    if(/^>/) { if($mainid and @fragids) { map{ $fragids{$_}= $mainid } @fragids; } @fragids=(); $mainid=0; }
    elsif(/^(\d+)/) { my $i=$1;
      m/(\d+)(nt|aa), >(.+)\.\.\. (.+)$/;  
      my($tlen,$typ,$tid,$pinfo)=($1,$2,$3,$4); 
      my $ismain=($pinfo =~ /\*/)?1:0;  # FIXME: drop /$/; new merge fnum at end of line now;
      unless($tid) {
        # $nerr++;
      } elsif($ismain) {
        $mainid= $tid; ## $mainlen=$tlen; $pi=100; # push @cluster, "$tid,$pi";
      } else {
        push @fragids, $tid;
      }
    }
  } close(F);
  return %fragids;
}

sub aacluster
{
  my($aaseq)=@_;
  return () unless(-s $aaseq and $aaseq !~ /\.gz$/);
  ## 1805: make -c 0.90 param? ie -c 0.95 may be safer .. or 0.98 or 1.0 
  
  my $aacdseq = makename($aaseq,"_cd90.aa");
  my $cdlog   = makename($aaseq,"_cd90.log");
  ## cd-hit  -M 1500  -c 0.9 -d 0;  -l  length of throw_away_sequences, default 10
  my $opts=""; 
  $opts.=" -T $NCPU" if($NCPU>1);
  $opts.=" -M $MAXMEM" if($MAXMEM>1); # M in Megabytes
  $opts.=" -l ".int($MINCDS/3) if($MINCDS>30); 

  my $cmd="$APPcdhit $opts -c 0.90 -d 0 -i $aaseq -o $aacdseq 1> $cdlog 2>&1"; # log to ...
  unless(-s $aacdseq) {
  my $runerr= runcmd($cmd);
  unlink("$aacdseq.bak.clstr") if (-f "$aacdseq.bak.clstr");
  # push @erasefiles,"$aacdseq.bak.clstr";
  }
  my $aaclstr= "$aacdseq.clstr";
  push @tmpfiles, $aacdseq, $aaclstr, $cdlog ;
  return($aacdseq,$aaclstr); ## want only "$aacdseq.clstr"
}


sub asmdupclass_sum
{
  #o# my($OUTH,$outclass,$aasize)=@_; 
  my($OUTH,$OUTB,$outclass,$aasize,$perfect_dups,$perfect_frags,$smallaa_dropids)=@_; 
  my($nt,%tab,%idclass);
  
  #UPD1712: write to separate file: trname.tr2aacds.info or .sum.txt	
  #UPD1712: add stats of drops: %perfect_dups, %perfect_frags, 
  use constant kPERFSTAT => 1;
  
  # open(OC,"cat $outclass | cut -f2,3 | sort | uniq -c |") or warn "ERR: asmdupclass_sum $outclass";
  unless( open(OC,$outclass) ) { loggit(1,"ERR: asmdupclass_sum missing: $outclass"); return -1; }

  print $OUTH  "# Class Table for $outclass \n";
  if($OUTB){ print $OUTB  "# Class Table for $outclass \n" ; }
  while(<OC>) {
    s/maybeok/okay/; s/altmidfrag/altmfrag/; # was amfrag
    my($id,$ac,$cl)=split; 
    for my $ct ($cl,"total") { $tab{$ct}{$ac}+=1; } $nt+=1;
    # my $fclass= ($ac =~ /drop/) ? "drop" : ($cl =~ /^alt/) ? "okalt" : "okay"; # is this bad default? problem was "amfrag" not altmidfrag!
    ## dang, add "part*" to alts,
    my $fclass= "none"; 
    if($ac =~ /^drop/) { $fclass="drop"; } elsif($ac =~ /^okay/) { $fclass= ($cl =~ /^alt|^part/)?"okalt":"okay"; } #|^amfrag
    $idclass{$id}= $fclass; 
  } close(OC);

if(kPERFSTAT) {
  my($droppdup,$mainpdup,$droppfrag,$dropsmall,$mainpfrag)=(0) x 9; 
  if(ref $perfect_dups) { $droppdup= scalar(keys %$perfect_dups); }
  if(ref $smallaa_dropids) { $dropsmall= scalar(keys %$smallaa_dropids); }
  if(ref $perfect_frags) { $droppfrag= scalar(keys %$perfect_frags); }
    # my %md=(); foreach my $pid (keys %$perfect_frags) { my $mid= $perfect_frags->{$pid}; $droppfrag++;  $md{$mid}++; } 
    # $mainpfrag= scalar(keys %md);  # ignore?
  if($droppdup or $droppfrag or $dropsmall) {
    $tab{'perfdupl'}{'drop'} += $droppdup;  $tab{'total'}{'drop'} += $droppdup;
    $tab{'perffrag'}{'drop'} += $droppfrag; $tab{'total'}{'drop'} += $droppfrag;
    $tab{'smallorf'}{'drop'} += $dropsmall; $tab{'total'}{'drop'} += $dropsmall;
    $nt += $droppdup + $droppfrag + $dropsmall;
  }
}  

  $nt||=1;
  my @ac=qw(okay drop); my @pac=map{ '%'.$_ } @ac;
  printf $OUTH "%-9s\t","class"; print $OUTH join("\t",@pac,@ac)."\n";
  if($OUTB){ printf $OUTB "%-9s\t","class"; print $OUTB join("\t",@pac,@ac)."\n"; }
  foreach my $cl (sort keys %tab) { 
    my @pv=(); my @nv=();
    foreach my $ac (@ac) { 
      my $n= $tab{$cl}{$ac}||0; push @nv,$n; 
      my $pf=$n/$nt; 
      my($fn,$fd)= ($pf<0.005)?(10000,100):(1000,10);
      my $p= int($fn*$pf)/$fd; push @pv,$p; 
    }
    if($cl eq "total"){ print $OUTH (("-") x 45); print $OUTH "\n"; } ;
    if($OUTB and $cl eq "total"){ print $OUTB (("-") x 45); print $OUTB "\n"; } ;
    printf $OUTH "%-9s\t",$cl; print $OUTH join("\t",@pv,@nv)."\n"; 
    if($OUTB){ printf $OUTB "%-9s\t",$cl; print $OUTB join("\t",@pv,@nv)."\n"; }
  } 
  print $OUTH (("=") x 45); print $OUTH "\n";
  if($OUTB){ print $OUTB (("=") x 45); print $OUTB "\n"; }
  
  ##  aastat for top1k of okay idclass;  evigene/scripts/prot/aastat.sh
  if(-s $aasize) {
    my $top=1000; 
    our(@aw,@nn,@aq); my($nok,$n,$sw,$sn,$nt,$nc)=(0) x 9; 
    #?? add same stats for okalt, drop sets? : @$aw{$class} @$nn{class}, $nok{class}
    ##? count aaqual: complete+utrok+gapok vs partial/utrbad/gapbad 
    open(AA,"sort -k2,2nr $aasize |"); # FIXME: sort TMPDIR=localdir
    while(<AA>) { next unless(/^\w/); my($id,$aw,$nn,$aqual)= split; $nt++;
      if($idclass{$id} =~ /okay/) { $nok++; push @aw,$aw; push @nn,$nn; push @aq,$aqual; } 
    } close(AA);
    
    # stat top1k and allokay
    sub aastat { 
      my ($OUTH,$OUTB,$name,$n)= @_; 
      $n=@aw if($n > @aw); my($sw,$sn,$nc)=(0,0,0); 
      for my $i (0..$n-1) { $sw+=$aw[$i]; $sn+=$nn[$i]; $nc++ if($aq[$i]=~/complete/); }
      my $aw=int($sw/$n); my $an=int(10*$sn/$n)/10;  my($mx,$md,$mi)= @aw[0,int($n/2),$n-1];
      #o#print $OUTH  "$name\t n=$n; average=$aw; median=$md; min,max=$mi,$mx; sum=$sw; gaps=$sn,$an\n";     
      print $OUTH  "$name\t n=$n; average=$aw; median=$md; min,max=$mi,$mx; nfull=$nc; sum=$sw; gaps=$sn,$an\n";     
      if($OUTB){print $OUTB  "$name\t n=$n; average=$aw; median=$md; min,max=$mi,$mx; nfull=$nc; sum=$sw; gaps=$sn,$an\n";}
      # as per aastat.sh, add nfull=$nc
    }
    print $OUTH  "# AA-quality for okay set of $aasize (no okalt): all and longest $top summary \n";
    if($OUTB){ print $OUTB  "# AA-quality for okay set of $aasize (no okalt): all and longest $top summary \n";}
    aastat($OUTH,$OUTB,"okay.top",$top);
    aastat($OUTH,$OUTB,"okay.all",$nok) unless($nok<$top);
    # bug for cacao3all7.aa.qual: okay.all n=51551 when only n=31341 are -drop -alt in cacao3all7.trclass
    # problem was "^amfrag" vs ^altmidfrag above .. no dupids in aa.qual; cacao3all7.aa n=3397790;
  }
  
}


# @outfiles= asmdupclass_fileset($outclass,$cdnaseq,$aaseq,$cdsseq); 
# fixme: add hdr classinfo from drops: %perfect_dups, %perfect_frags, 

# our %EVG_SEQSUFFIX = ( 'tr' => 'tr|cdna', 'aa' => 'aa|pep', 'cds' => 'cds', 
#                        'mrna' => 'mrna', 'ncrna' => 'ncrna'); # UPD1911
                      
sub asmdupclass_fileset
{
  my($outclass,$cdnaseq,$aaseq,$cdsseq,$perfect_dups,$perfect_frags,$smallaa_dropids)=@_;
  my($nt,%idclass,%classinfo,@outfiles,%outfiles);

  # UPD1908: maybe add 'oknoc', putative noncode, to okay okalt class set
  # UPD1911: cdnaseq.fa|.other suffix mixup >> resolve by changing all subdir suffices to .tr
  
  my %seqset= ( 'tr' => $cdnaseq, 'aa' => $aaseq, 'cds' => $cdsseq);
if(UPD1911) {
  for my $sf (qw( tr aa cds)) {
    my $inf = $seqset{$sf};
    my($suf)= $inf =~ m,(\.[\.\w]+)$,; $suf =~ s/\.gz//;
    my $sok= $EVG_SEQSUFFIX{$sf} || $sf;
    unless($suf and $suf =~ m/\.($sok)/) { $suf=".$sf"; }
    
    foreach my $class (qw(okay okalt drop)) {  # 
      my $cname= makename($inf,".$class$suf");
      push @outfiles, $cname if($dryrun or -s $cname);  
    }
  }
 
} else {
  # check/skip for existing fileset
  foreach my $inf ($cdnaseq,$aaseq,$cdsseq) 
  {
    my($suf)= $inf =~ m,(\.[\.\w]+)$,;   $suf =~ s/\.gz//;
    foreach my $class (qw(okay okalt drop)) {  # 
      my $cname= makename($inf,".$class$suf");
      push @outfiles, $cname if($dryrun or -s $cname);  
    }
  }
}
  return @outfiles if(@outfiles > 3);

  if(ref $perfect_dups) { foreach my $pid (keys %$perfect_dups) { 
    my $mainid= $perfect_dups->{$pid}; $classinfo{$pid}= "perfectdup,drop,match:$mainid"; } }
  if(ref $perfect_frags) { foreach my $pid (keys %$perfect_frags) { 
    my $mainid= $perfect_frags->{$pid}; $classinfo{$pid}= "perfectfrag,drop,match:$mainid"; } }
  if(ref $smallaa_dropids) { foreach my $pid (keys %$smallaa_dropids) { 
    my $aalen= $smallaa_dropids->{$pid}; $classinfo{$pid}= "smallaa,drop,size:$aalen"; } }

  unless( open(OC,$outclass) ) { loggit(1,"ERR: asmdupclass_fileset missing: $outclass"); return -1; }
  while(<OC>) {
    next unless(/^\w/); 
    s/maybeok/okay/; # s/altmidfrag/amfrag/; 
    my($id,$ac,$cl,$bestid,$pia,$aaqual,$flags)=split; 
    my $fclass= ($ac =~ /drop/) ? "drop" : ($cl =~ /^alt|^part/) ? "okalt" : "okay";
    # UPD1908: $fclass="oknoc" if($ac =~ /okay/ and $cl =~ /nc$/);
    $idclass{$id}= $fclass; $nt++;
    my $match= ($bestid eq $id)?"":",match:$bestid,pct:$pia";
    # fixme: match bestid == thisid ; ie no match for noclass...
    if(AACONS) { # add consensus from trclass.flags to hdr: aacons:3
      if( my($aacons)= $flags=~m/(aacons:\w+)/ ) { $match .= ",$aacons"; }
    }
    
    $classinfo{$id}= "$cl,$ac$match; aalen=$aaqual;";
    # .. above classinfo for perfectdups from 1.fastanrdb, 2. perfectfragments cdhitest
  } close(OC);
  if($nt < 2) { loggit(1,"ERR: asmdupclass_fileset classes=$nt"); return -1; }
  
  @outfiles= (); %outfiles= ();

if(UPD1911) {
  for my $sf (qw( tr aa cds)) {
    my $inf = $seqset{$sf};
    my($suf)= $inf =~ m,(\.[\.\w]+)$,; $suf =~ s/\.gz//;
    my $sok= $EVG_SEQSUFFIX{$sf} || $sf;
    unless($suf and $suf =~ m/\.($sok)/) { $suf=".$sf"; }

    foreach my $class (qw(okay okalt drop)) {
      my $cname= makename($inf,".$class$suf");
      if(-s $cname) { rename($cname,"$cname.old"); } # or what?
      my $chandle=undef;
      my $okw= open( $chandle, ">$cname");
      if($okw) {
        push @outfiles, $cname;
        $outfiles{"$inf.$class"}= $chandle;
        push @okayset, $cname unless($class eq "drop");
        push @dropset, $cname if($class eq "drop");
      } else {
        loggit(1,"ERR: writing $cname");
      }
    }    
  }

} else {
  foreach my $inf ($cdnaseq,$aaseq,$cdsseq) {
    
    my($suf)= $inf =~ m,(\.[\.\w]+)$,;  $suf =~ s/\.gz//;
    
    foreach my $class (qw(okay okalt drop)) {
      my $cname= makename($inf,".$class$suf");
      if(-s $cname) { rename($cname,"$cname.old"); } # or what?
      my $chandle=undef;
      my $okw= open( $chandle, ">$cname");
      if($okw) {
        push @outfiles, $cname;
        $outfiles{"$inf.$class"}= $chandle;
        push @okayset, $cname unless($class eq "drop");
        push @dropset, $cname if($class eq "drop");
      } else {
        loggit(1,"ERR: writing $cname");
      }
    }
  }
}
  
  foreach my $inf ($cdnaseq,$aaseq,$cdsseq) {
    my $okin= 0; 
    if($inf =~ /\.gz$/) { $okin= open(IN,"gunzip -c $inf|"); }
    else { $okin= open(IN,$inf); }
    loggit(1,"ERR: asmdupclass_fileset reading: $inf") unless($okin); 
    my $chandle=undef;
    my $CHECKUTRORF=($inf eq $cdnaseq)?1:0; # 1806 bugfix: missing cdna(utrorf) for okay aa.utrorf only
    ## BUT insure no re-rev : UPDATE.maybe: revcomp(cdnaseq) if aaseq strand=- ??? or instead traa2cds -trout
    while(<IN>) {
      if(/^>(\S+)/) { my $id=$1; 
        my $class= $idclass{$id} || "drop"; # NOTE: 1,2 perfect dups are not in idclass; dropped already
        my $isuo=0; if($CHECKUTRORF and $class eq "drop") {
          if(my $clutr= $idclass{$id.'utrorf'}) { $class= $clutr; $isuo=1; }
        }
        $chandle= $outfiles{"$inf.$class"}; # fail if missing?
        ## FIXME: want classinfo also for 1. perfdups, 2.perffrags, need ids from those files
        if(my $classinfo= $classinfo{$id}) { # if($isuo) id.utrorf; $CHECKUTRORF?
          $classinfo =~ s/ aalen=\S+// if(/aalen=/);
          s/evgclass=\S+//; s/$/ evgclass=$classinfo/; 
          }
        loggit(1,"ERR: missing output id=$id; $inf.$class") unless(ref $chandle);        
      }
      print $chandle $_ if(ref $chandle); ## bad chandle
    } close(IN);
  }

## add aa.qual per okayset, dropset .. pull from aaseq, not orig.aa.qual, add evgclass= info
## grep -v drop $pt.trclass | cut -f1,2 | sed 's/okay//' | ggrep -F -f - inputset/$pt.aa.qual > okayset/$pt.aa.qual

  foreach my $inf ($cdnaseq,$aaseq,$cdsseq) {
    foreach my $class (qw(okay okalt drop)) {
      my $ohand= $outfiles{"$inf.$class"};
      close($ohand) if(ref $ohand);
    }
  }
  
  return(@outfiles);
}



sub asmdupfilter_cds
{
  my($cdsblast,$aasize,$aacdhit)=@_;

  my $insuf='aa.qual|aa.size|aa|fasta';
  my $outclass= makename($aasize,".trclass",$insuf);
  my $outaln  = makename($aasize,".alntab",$insuf); # change to .cds.alntab ?
  my $aflog   = makename($aasize,".adupfilt.log",$insuf);
  ## my $aacdhit = makename($aasize,".aa.clstr");
  my $aaclstropt= ( $aacdhit and -s $aacdhit ) ? "-acdhit $aacdhit" : "";
  my $dbg=($debug)? "-debug" : "";

	## opt for cdsblast ==  lastz : $USE_LASTZ global, file=name-self97.lastz
	my $cdsblastop= ( $cdsblast =~ m/blast/ ) ? "-blastab $cdsblast"
			: ( $cdsblast =~ m/lastz/ or $USE_LASTZ ) ? "-lastz $cdsblast" : "-blastab $cdsblast";

  
	##  opt for -aablast=$aablast, table of trid refid bitscore identity align
	##  for homology scoring to preserve from drops..
	## FIXME: allow -ablastab to be .names table in asmrna_dupfilter2.pl, ie decide by file name
	## add opt trdupfilter  -aanames trset.names == table($td,$name,$alnscore,$rd,), used also for  namegenes
	## -ablastab "trset.names" works same way (names suffix aware)
  my $aablastopt= ( $aablast and -s $aablast ) ? "-ablastab $aablast" : "";

  if(AACONS) {
    #? for noAACONS, add -aconsensus no 
    my $cname= makename($cdnaseq,".consensus"); #? global? need for asmrna_dupfilter3c.pl
    if($USE_AACONS) {
      if(-f $cname) { $aablastopt .= " -aconsensus $cname"; }
      else { $aablastopt .= " -aconsensus "; } # this way uses trset.nrcds headers .. skip or no?
    } else {
      $aablastopt .= " -aconsensus no"; # or leave off for back-compatability
    }
  }
  
  if(UPD1908) { # add user/env opts to asmdupfilter this way?
    if(my $pv=$ENV{pCDSOK}){ $aablastopt .= " -pCDSOK=$pv"; } # d:60 .. too high for alt drops
    if(my $pv=$ENV{pCDSBAD}){ $aablastopt .= " -pCDSBAD=$pv"; } # d:30
    if(my $pv=$ENV{ALTFRAG}){ $aablastopt .= " -ALTFRAG=$pv"; } # d:0.5, alt/part class
  }
  if(UPD1912) { # now opt $PHETERO # my $PHETERO= $ENV{PHETERO}||0; 
    if($PHETERO>0 and $PHETERO<9) { $aablastopt .= " -pHETERO=$PHETERO";} 
    }
  

  ##... simplest..
  # $evigene/scripts/rnaseq/asmrna_dupfilter2.pl -debug -CDSALIGN -tinyaln 35 
  #  -aasize $pt.aa.qual -blast slf95-$pt*blastn -outeq $pt.baln5 -outclass $pt.bclass5 > & log.dupfltb5.$pt
  ##... more inputs.. skip -dupids
  # $evigene/scripts/rnaseq/asmrna_dupfilter2.pl -debug -CDSALIGN -tinyaln 35 -aasize $pt.aa.qual 
  #  -acdhit $pt.aa.clstr -dupids $pt.dupids -blast outz/sd-slf95-$pt*blastn.gz 
  #  -outeq $pt.baln5c -outclass $pt.bclass5c
  
  ## UPD 2016.02: -tinyaln 35  bad (sometimes), need opts and diff test: tiny by bases not pct overlap
  ## ie. cds-exon overlap should be valid for >= tinybases, eg. 90? need to sample ref species alt exon sizes
  ## my $MINALIGN= 90; # NOT USED NOW; use just TINYALN; change to several levels

  ## UPD 1602: removed -tinyaln 35
  ## UPD1912: dupfilter4, -CDSALIGN is default, -PHETERO is option here and dupfilt, dont need ENV{PHETERO}
  my $cmd="$APPtrdupfilter $dbg  -aasize $aasize -CDSALIGN $cdsblastop $aablastopt $aaclstropt"
    ." -outeqtab $outaln -outclass $outclass >$aflog 2>&1";
    
  unless(-s $outclass) {
  my $runerr= runcmd($cmd);
  }
  push @tmpfiles, $outaln, $aflog; ## outclass is main output file?
  return($outclass,$outaln);
}
  
sub make_aaqual
{
  my($aaseq)=@_;
  
  # UPD1908: read new attr: codepot=Noncode/Code/.. 
  # fixed minor update bug: this always makes new, should instead reuse old aasize..
  my $aasize= makeAaQual($aaseq,1); # use constant makeAaQual_SETHASH=>1;
  
  # my $aasize= makename($aaseq,".aa.qual"); 
  # my $cmd="$APPaaqual $aaseq"; ## copy sub here/cdna_evigenesub.pm
  # unless(-s $aasize) { my $runerr= runcmd($cmd); }
  push @inputset, $aasize;
  return($aasize);
}


# sub fadupids #  in cdna_evigenesub.pm
# sub facount #  in cdna_evigenesub.pm
# sub faidlist #  in cdna_evigenesub.pm
# sub faextract #  in cdna_evigenesub.pm
# sub fasize #  in cdna_evigenesub.pm
# sub fasplit #  in cdna_evigenesub.pm


#  make_bestorf_ncpu: parallelized for NCPU by splitting input cdnaseq to NCPU parts; cat parts> whole.
sub make_bestorf_ncpu
{
  my($cmdin, $npart,$cdnaseq,$cdsseq,$aaseq)=@_;
  #o# my($npart,$cdnaseq,$cdsseq,$aaseq)=@_;
  return unless(-s $cdnaseq);
  
  my $MINAA=int($MINCDS/3);
  my $csize= fasize($cdnaseq);
  my $spldir= makename($cdnaseq,"_split/");
  my $splsize= 1 + int($csize/$npart);
  
  mkdir($spldir); # dryrun?
  my @splset= fasplit( $cdnaseq, $spldir, $npart, $splsize); 
  my (@cdsset,@aaset);
  my $icpu= 0; 
  foreach my $cdna1 (@splset) {
    # my $aa1 = makename($cdna1,".aa");  # BUG fixed? makename chopped icpu index; always add here? not a conflict?
    # my $cds1= makename($cdna1,".cds");
    my $aa1  = makename($cdna1,".aa$icpu"); 
    my $cds1 = makename($cdna1,".cds$icpu"); 
    my $cmd1= $cmdin . " -cdna $cdna1 -aaseq $aa1 -cdsseq $cds1";
    #o my $cmd1="$APPcdnabest -nostop -minaa=$MINAA -cdna $cdna1 -aaseq $aa1 -cdsseq $cds1";
    #o $cmd1 =~ s/nostop/nostop -noutrorf / if($noutrorf); # .=" -noturorf" if($noutrorf); ##201609 NO NOT AT END

    push @aaset, $aa1; push @cdsset, $cds1;
    my $pid= forkcmd($cmd1);    
    if(++$icpu > $npart) { while (wait() != -1) { }; $icpu= 0; }
  }
  ## THIS FAILED TO WAIT for pids ...
  ## not enough: wait(); # waitpid(@pids); ??
  while (wait() != -1) { };
  
  my $cmd; #? cat ok or do perl read/print ?
  cat_splitset($aaseq,\@aaset); # ret: ($runerr, $nok, \@ofail);
  cat_splitset($cdsseq,\@cdsset);
  # $cmd= "cat ".join(' ',@aaset)." > $aaseq";   runcmd($cmd);
  # $cmd= "cat ".join(' ',@cdsset)." > $cdsseq"; runcmd($cmd);
  
  #?? are these useful parts to keep for later?
  #was ($debug||$dryrun) # upd1712: default nodebug
  push @erasefiles, @splset, @aaset, @cdsset, $spldir; # which?
  
  return($cdsseq,$aaseq);
}


sub make_bestorf
{
  my($cdnaseq,$rnaIsStranded)=@_;
  my $MINAA=int($MINCDS/3);
  my $aaseq = makename($cdnaseq,".aa");
  my $cdsseq= makename($cdnaseq,".cds");
  # FIXME: parallelize this for NCPU by splitting input cdnaseq to NCPU parts; cat parts> whole.

  #o# my $cmd="$APPcdnabest -nostop -minaa=$MINAA -cdna $cdnaseq -aaseq $aaseq -cdsseq $cdsseq";
  my $cmdin  = "$APPcdnabest -nostop -minaa=$MINAA";
     $cmdin .= " -noutrorf" if($noutrorf); # .=" -noturorf" if($noutrorf); ##201609 NO NOT AT END
     $cmdin .= " -codepot" if(UPD1908); # -codepotential new option, always? 
     $cmdin .= " -ostrand=fwd" if($rnaIsStranded and $rnaIsStranded=~/yes|1/i); #UPD22a
  
  unless(-s $cdsseq) {
    if($NCPU>1 and not $dryrun) { ($cdsseq,$aaseq)= make_bestorf_ncpu($cmdin, $NCPU,$cdnaseq,$cdsseq,$aaseq); }
    else { runcmd($cmdin . " -cdna $cdnaseq -aaseq $aaseq -cdsseq $cdsseq"); } # FIX20sep: ',' => '.'
  }
  return($cdsseq,$aaseq);
}


sub make_cdsFromTrAA
{
  my($cdnaseq,$aaseq,$cdsout,$makeopt)=@_;
  $makeopt||="";
  $cdsout= makename($cdnaseq,".cds") unless($cdsout);
  my $cmd="$APPtraa2cds $makeopt -cdna $cdnaseq -aaseq $aaseq -cdsseq $cdsout";
  unless(-s $cdsout) {
    my $runerr= runcmd($cmd);
  }
  return($cdsout);
}

sub get_bestorf
{
  my($cdnaseq,$aaseq,$cdsseq,$rnaIsStranded)=@_; #UPD22a: $rnaIsStranded
  $aaseq = makename($cdnaseq,".aa") unless($aaseq); # if(defined $aaseq and not $aaseq);
  $cdsseq= makename($cdnaseq,".cds") unless($cdsseq); # if(defined $cdsseq and not $cdsseq);
  if( $cdsseq and -s $cdsseq ) {
    # have already ..  
  } elsif( $aaseq and -s $aaseq ) {
    ($cdsseq) = make_cdsFromTrAA($cdnaseq,$aaseq);
  } else {
    ($cdsseq,$aaseq) = make_bestorf($cdnaseq,$rnaIsStranded);
  }
  push @inputset, $cdsseq,$aaseq ;
  return($cdsseq,$aaseq);
}

=item UPD1912 blastn_cds -strand plus or both ?

   UPD1912a: blastn -strand plus  .. gets rid of -sense problematic CDS aligns
   * tested stranded align years back, but kept both w/ -sense classify : WHY?  
   maybe used as poor cds assembly qual? but get many valid cds genes w/ partial antisense aligns
   
   * WHY == antisense transcripts, mostly artifact, some real, are "properly" classed into same locus
      as fwd transcripts.  But this has problems, real antisense tr are considered separate loci,
      and some are located outside of fwd locus as other loci, some are best conserved genes
      But antisense *proteins* are not known biologically, vs noncode antisense, and are
      artifacts of various trasm errors (eg. gaps spanning rev stop codon), and biological nonsense rev prots
      
   UPD1912b: 
   * need to add new classifier of antisense tr, combine into self-mrna align measure,
   add class meaures of anti-cds, utr-orf (either strand), mixed/ambiguous self-mrna-part-aligns
   following primary reduction using cds -strand align.  NO,
   
   UPD1912c:
   * strand problem is self-evid lacks enough to decide which of fwd/rev protein is biological,
     revert to cds-align both strands, with -sense flagged for classifier, later resolution
   * need external evid: prot homol or/and chr-align strand (from introns)
  
=cut 

sub blastn_cds_opts
{
  my $opts=""; 

  ## -dust no/yes ? does it make a diff?  -ungapped may need to be option..
  ## 11mar.FIXME: need blastn -ungapped otherwise miss perfect match of parts ; e.g. alt-tr half=perfect, other=imperf
  ## 16feb.FINALLY: "-ungapped -xdrop_ungap 4" cures problem of non-local basic-local-alignments, -xdrop 20 default allows extension into low-ident aligns

  ## 16feb : TEST blastn -ungapped -xdrop_ungap 4 << reduce xdrop to prevent extending identical exon aligns to nonident alt exons
  ## .. this may be cure for miscalling alts as new loci.  needs tests of alts vs paralogs vs hetero transcripts
  ## default xdrop_ungap 30 means extending ident align to non-ident w/ drop in bitscore of 30 ? low val reduces this

  ## NO UPD1912: blastn -strand plus  .. gets rid of -sense problematic CDS aligns ** NO GOOD, turns -sense into many new loci,

  if(my $blopt=$ENV{BLASTNOPT}) { 
    $opts .= " $blopt"; 
    $opts .= " -perc_identity $CDSBLAST_IDENT" unless($opts=~/perc_identity/);
    $opts .= " -evalue $CDSBLAST_EVALUE" unless($opts=~/evalue/);
    }
  else { 
    $opts .= " -ungapped -xdrop_ungap 4 -dust no"; 
    #no $opts .= " -strand plus" if(UPD1912); #NO, but maybe as option? 
    $opts .= " -perc_identity $CDSBLAST_IDENT -evalue $CDSBLAST_EVALUE";
    }  
  
  ## caller add: $opts.=" -num_threads $NCPU" if($NCPU>1); # -num_threads isnt effective
  ##nothere# my $cmd="$APPblastn -task megablast $opts -outfmt 7 -db $cdsdb -query $cdsseq -out $cdsbltab";

  return $opts; 
}

sub blastn_cds
{
  my($cdsseq)=@_;
  if($NCPU>1 and not $dryrun) { return blastn_cds_ncpu($NCPU,$cdsseq); }
 
  ## *** clusterize this .. NCPU works ok .. NOT!
  my $cdsdb= makename($cdsseq,"_db");
  my $cdsbltab= makename($cdsseq,"-self$CDSBLAST_IDENT.blastn");
  my $blog= makename($cdsdb,".log","xxxsuf");
  return($cdsbltab) unless(-s $cdsseq);

## 11mar.FIXME: need blastn -ungapped otherwise miss perfect match of parts ; e.g. alt-tr half=perfect, other=imperf
## 16feb.FINALLY: "-ungapped -xdrop_ungap 4" cures problem of non-local basic-local-alignments, -xdrop 20 default allows extension into low-ident aligns
## 14mar : Very Sloooowwwww for daphmag, cluster not using NCPU effectively w/ -num_threads (e.g get 2x speedup for 32 cpu)
##   .. redo this w/ fasplit(query) and ncpu x 1-cpu parts.

## 16feb : TEST blastn -ungapped -xdrop_ungap 4 << reduce xdrop to prevent extending identical exon aligns to nonident alt exons
## .. this may be cure for miscalling alts as new loci.  needs tests of alts vs paralogs vs hetero transcripts
## default xdrop_ungap 30 means extending ident align to non-ident w/ drop in bitscore of 30 ? low val reduces this

  my $fmtcmd="$APPmakeblastdb -in $cdsseq -dbtype nucl -out $cdsdb -logfile $blog";
  
  # my $opts=""; 
  ## -dust no/yes ? does it make a diff?  -ungapped may need to be option..
  # if(my $blopt=$ENV{BLASTNOPT}) { $opts .= " ".$blopt; }
  # else { $opts.= " -ungapped -dust no "; } #  -evalue $CDSBLAST_EVALUE ??
  # $opts .= " -perc_identity $CDSBLAST_IDENT -evalue $CDSBLAST_EVALUE";

  my $opts= blastn_cds_opts();
  
  $opts.=" -num_threads $NCPU" if($NCPU>1); # -num_threads isnt effective
  my $cmd="$APPblastn -task megablast $opts -outfmt 7 -db $cdsdb -query $cdsseq -out $cdsbltab";
    
  unless(-s $cdsbltab) {
  runcmd($fmtcmd);
  my $runerr= runcmd($cmd);
  }
  push @erasefiles, "$cdsdb.nsq","$cdsdb.nin","$cdsdb.nhr";
  push @tmpfiles, $cdsbltab,$blog;
  return($cdsbltab);    
}

  # parallelize blastn by splitting query.fa to NCPU parts; blastn -num_threads not effective.
sub blastn_cds_ncpu
{
  my($npart,$cdsseq)=@_;
  
  my $cdsdb= makename($cdsseq,"_db");
  my $cdsbltab= makename($cdsseq,"-self$CDSBLAST_IDENT.blastn");
  my $blog= makename($cdsdb,".log","xxxsuf");
  return($cdsbltab) unless(-s $cdsseq);
  
  # my $opts=""; 
  # if(my $blopt=$ENV{BLASTNOPT}) { $opts .= " ".$blopt; }
  # else { $opts.= " -ungapped -dust no "; } #  -evalue $CDSBLAST_EVALUE ??
  # $opts .= " -perc_identity $CDSBLAST_IDENT -evalue $CDSBLAST_EVALUE";
 
  my $opts= blastn_cds_opts();

  ## FIXME 1604: check output parts for failed runs, rerun if found and few .. eg. memoverload sometimes
  
  my $blcmd0="$APPblastn -task megablast $opts -outfmt 7 -db $cdsdb "; # add parts: -query $cdsseq -out $cdsbltab
  my $fmtcmd="$APPmakeblastdb -in $cdsseq -dbtype nucl -out $cdsdb -logfile $blog";

  unless(-s $cdsbltab) {
    runcmd($fmtcmd);
    push @erasefiles, "$cdsdb.nsq","$cdsdb.nin","$cdsdb.nhr";

    my $csize= fasize($cdsseq);
    my $spldir= makename($cdsseq,"_blsplit/");
    my $splsize= 1 + int($csize/$npart);

    mkdir($spldir); # dryrun?
    my @splset= fasplit( $cdsseq, $spldir, $npart, $splsize); 
    my (@bloset);
    my $icpu= 0; my $ipart=0;
    foreach my $cds1 (@splset) {
      $ipart++;
      my $cdsbltab1= makename($cds1,"-self$CDSBLAST_IDENT.blastn$ipart");
      my $cmd1= $blcmd0 . " -query $cds1 -out $cdsbltab1";  # add parts: -query $cdsseq -out $cdsbltab
      push @bloset, $cdsbltab1;
      my $pid= forkcmd($cmd1);    
      if(++$icpu > $npart) { while (wait() != -1) { }; $icpu= 0; }
      else { sleep(2); } # fork/wait bug on old solaris unix, needs delay here maybe?
    }
    while (wait() != -1) { };
    
  ## FIXME 1604 here: check blast outputs for failed runs, report, rerun if found and few
  use constant FIXME_blastn_cds_ncpu_fails => 1;
  my $nfail=0;
  if(FIXME_blastn_cds_ncpu_fails) {  
    my $FAILDIE= ($debug) ? LOG_WARN : LOG_DIE;
    # FIXME 1606: need @$ofailset, @rerunset update to @bloset
    my($runerr,$nok,$failset,$ofailset)= blastn_collate_result($npart, $cdsbltab, \@bloset,\@splset);
    if($runerr and ref $failset) {
      $nfail= @$failset; ## nfail == runerr now
      if($nfail >= $npart/2) { 
        loggit($FAILDIE,"ERR=blastn_cds_ncpu nfail=$nfail, nok=$nok, fatal cmd=",$blcmd0); 
      } elsif($nfail>0) {
        loggit(LOG_WARN,"ERR=blastn_cds_ncpu, nfail=$nfail, nok=$nok, rerun cmd=",$blcmd0); 
        $icpu= 0; $ipart=0;
        my @rerunset=();
        foreach my $cds1 (@$failset) {
          $ipart++;
          my $cdsbltab1= makename($cds1,"-self$CDSBLAST_IDENT.blastn$ipart");
          my $cmd1= $blcmd0 . " -query $cds1 -out $cdsbltab1";  # add parts: -query $cdsseq -out $cdsbltab
          push @rerunset, $cdsbltab1; # same output will rewrite fails
          my $pid= forkcmd($cmd1);    
          if(++$icpu > $npart) { while (wait() != -1) { }; $icpu= 0; }
        }
        while (wait() != -1) { };

        # reuse @bloset or merge @rerunset ? BUG HERE using @$ofailset again from @bloset, diff names
        my %ofailn= map{ $_=>1 } @$ofailset;
        my @bloset2= grep{ not $ofailn{$_} } @bloset; 
        push @bloset2, @rerunset;
        
        my($runerr2,$nok2,$failset2)= blastn_collate_result($npart, $cdsbltab, \@bloset2,\@splset);
        if($runerr2 and ref $failset2) {
          $nfail= @$failset2;
          loggit($FAILDIE,"ERR=blastn_cds_ncpu rerun nfail=$nfail, nok=$nok2, fatal cmd=",$blcmd0); 
        }
      }
    }
    
  } else { # orig cat output parts
    cat_splitset($cdsbltab,\@bloset);
    # my $cmd= "cat ".join(' ',@bloset)." > $cdsbltab"; 
    # my $runerr= runcmd($cmd); #??
  }
  
    if($nfail>0) { #  and $debug
      # dont erase..
    } else { #upd1712, was $debug||$dryrun
      push @erasefiles, @splset, @bloset, $spldir;  #? do always?
    } 
    
  }

  push @tmpfiles, $cdsbltab,$blog;
  return($cdsbltab);    
}

sub blastn_collate_result
{
  my($npart, $cdsbltab, $bloset, $splset)=@_;
  my $nin= @$splset; 
  my $nok= 0;
  # $OKLINE = "# BLAST processed 3833 queries";
  my (@qfail,@ofail);
  if(-s $cdsbltab) { } # save? leave partial result.
  my $ok= open(O,'>', $cdsbltab);
  for my $i (0..$nin-1) {
    my $outf= $bloset->[$i];
    my $qseq= $splset->[$i];
    my $isend= 0;
    if($outf and -s $outf) {
      my $lastl="";      
      $ok= open(I,$outf); while(<I>) { $lastl=$_; print O $_; } close(I);
      $isend++ if( $lastl =~ m/^# BLAST processed/); 
    }
    if($isend>0) { $nok++; } 
    else { push @qfail, $qseq; push @ofail, $outf; }
  }
  close(O);
  # if($runerr) { loggit(1,"ERR=$err ",$cmd[0]); } # ..
  
  my $runerr= $nin - $nok;
  my $failset= \@qfail; # add ofail to remove from bloset..
  return ($runerr, $nok, $failset, \@ofail);
}

__END__

=item UPD1908 drop lastz_cds

  out of date code ; ubiq blastn is working well enough, 
  
=cut

# sub lastz_cds
# {
#   my($cdsseq)=@_;
#   if($NCPU>1 and not $dryrun) { return lastz_cds_ncpu($NCPU,$cdsseq); }
#  
#   my $cdsbltab= makename($cdsseq,"-self$CDSBLAST_IDENT.lastz");
#   my $blog= makename($cdsbltab,".log","lastz");
#   return($cdsbltab) unless(-s $cdsseq);
# 
#   my $opts=""; 
#   if(my $blopt=$ENV{LASTZOPT}) { $opts = $blopt; }
# 	else { 
# 	$opts="--identity=$CDSBLAST_IDENT --coverage=20 --strand=plus --step=10 --seed=match12 --filter=nmatch:90 --gfextend --notransition --exact=20 --match=1,5 --nochain --nogapped --ambiguous=n";
# 	}
#   ##my $cmd="$APPlastz \"$cdsseq\[multiple\]\" $cdsseq  $opts --format=general --output=$cdsbltab"; 
#   my $cmd= $APPlastz .' "' . $cdsseq . '[multiple]" ' . $cdsseq . ' ' . $opts .' --format=general --output='.$cdsbltab; 
# 
#   unless(-s $cdsbltab) {
#   my $runerr= runcmd($cmd);
#   }
#   push @tmpfiles, $cdsbltab, $blog;
#   return($cdsbltab);    
# }
# 
# sub lastz_cds_ncpu
# {
#   my($npart,$cdsseq)=@_;
#   
#   my $cdsbltab= makename($cdsseq,"-self$CDSBLAST_IDENT.lastz");
#   my $blog= makename($cdsbltab,".log","lastz");
#   return($cdsbltab) unless(-s $cdsseq);
#   
#  	## lastz opts tested; not sure how many are needed yet. PI,PA main ones.
#  	## FIXME: ** before use, this needs more tests w/ other options; got wacky results e.g. self partial aligns, not full. 
# 	## PI=97 == CDSBLAST_IDENT; PA=20;
# 	## lseedopt="--step=10 --seed=match12"
# 	## lopt="--identity=$PI --coverage=$PA --filter=nmatch:90 --gfextend --notransition --exact=20 --match=1,5 --nochain --nogapped --ambiguous=n"
#   ## $bbin/lastz "$trf[multiple]" $qfile $lseedopt $lopt --format=general --output=$onam.lzout  &
#   ##
#   ## ?? add --strand since we have stranded CDS-seq
#   ## --strand=plus   search + strand only (matching strand of query spec)
#   ## --coverage=$CDSLASTZ_COVER ??
#   ## maybe not? --exact=20 --match=1,5 .. problems for snp,indel? 
#   my $opts=""; 
#   if(my $blopt=$ENV{LASTZOPT}) { $opts = $blopt; }
# 	else { 
# 	$opts="--identity=$CDSBLAST_IDENT --coverage=20 --strand=plus --step=10 --seed=match12 --filter=nmatch:90 --gfextend --notransition --exact=20 --match=1,5 --nochain --nogapped --ambiguous=n";
# 	}
# 	
# 	## lastz target[xxx] query  options : is recommended form, does option order mater?
#   ##my $blcmd0="$APPlastz \"$cdsseq\[multiple\]\" QUERYHERE $opts --format=general "; # per part :  --output=$onam.lzout $qfile
#   my $blcmd0= $APPlastz .' "' . $cdsseq . '[multiple]" QUERYHERE ' . $opts .' --format=general '; # per part :  --output=$onam.lzout $qfile
# 
#   unless(-s $cdsbltab) {
#     my $csize= fasize($cdsseq);
#     my $spldir= makename($cdsseq,"_blsplit/");
#     my $splsize= 1 + int($csize/$npart);
# 
#     mkdir($spldir); # dryrun?
#     my @splset= fasplit( $cdsseq, $spldir, $npart, $splsize); 
#     my (@bloset);
#     my $icpu= 0; my $ipart=0;
#     foreach my $cds1 (@splset) {
#       $ipart++;
#       my $cdsbltab1= makename($cds1,"-self$CDSBLAST_IDENT.lastz$ipart");
#       my $cmd1= $blcmd0 . " --output=$cdsbltab1";  $cmd1 =~ s/QUERYHERE/$cds1/;  
#       push @bloset, $cdsbltab1;
#       my $pid= forkcmd($cmd1);    
#       if(++$icpu > $npart) { while (wait() != -1) { }; $icpu= 0; }
#     }
#     while (wait() != -1) { };
#     
#     cat_splitset($cdsbltab,\@bloset);
#     # my $cmd= "cat ".join(' ',@bloset)." > $cdsbltab"; 
#     # my $runerr= runcmd($cmd); #??
# 
#     push @erasefiles, @splset, @bloset, $spldir;  
#   }
# 
#   push @tmpfiles, $cdsbltab,$blog;
#   return($cdsbltab);    
# }


# sub findapp  #  in cdna_evigenesub.pm
# sub findevigeneapp  #  in cdna_evigenesub.pm
# sub runcmd #  in cdna_evigenesub.pm
# sub forkcmd #  in  cdna_evigenesub.pm
# sub makename #  in cdna_evigenesub.pm




=item test1

  /bio/bio-grid/aabugs4/crusts/shrimpt/best1test
  log.tr2aacds1                     shrimt1trin1.drop.cds             shrimt1trin1nr.cds
  log.tr2aacds2                     shrimt1trin1.drop.tr              shrimt1trin1nrcd1.blastn
  log.tr2aacds3                     shrimt1trin1.okalt.aa             shrimt1trin1nrcd1.cds
  log.tr2aacds4                     shrimt1trin1.okalt.cds            shrimt1trin1nrcd1.cds.bak.clstr
  shrimt1trin1.aa                   shrimt1trin1.okalt.tr             shrimt1trin1nrcd1.cds.clstr
  shrimt1trin1.aa.qual              shrimt1trin1.okay.aa              shrimt1trin1nrcd1.log
  shrimt1trin1.adupfilt.log         shrimt1trin1.okay.cds             shrimt1trin1nrcd1_db.nhr
  shrimt1trin1.alntab               shrimt1trin1.okay.tr              shrimt1trin1nrcd1_db.nin
  shrimt1trin1.cds                  shrimt1trin1.tr.gz                shrimt1trin1nrcd1_db.nsq
  shrimt1trin1.drop.aa              shrimt1trin1.trclass

  >> cleanup
  dropset/               log.tr2aacds2          okayset/               tmpfiles/
  inputset/              log.tr2aacds3          shrimt1trin1.tr.gz
  log.tr2aacds1          log.tr2aacds4          shrimt1trin1.trclass
  
  OUTPUT seq:
    shrimt1trin1.{okay,okalt}.{tr,cds,aa} >> okayset/
    shrimt1trin1.drop.{tr,cds,aa} >> dropset/

  INPUT/made seq: >> inputset/
    shrimt1trin1.aa
    shrimt1trin1.aa.qual
    shrimt1trin1.cds
    shrimt1trin1nr.cds
    shrimt1trin1nrcd1.cds
      
  TEMP files to erase >> tmpfiles/ # (move to tmpfile subdir?)
    shrimt1trin1nrcd1_db.{nsq,nhr,nin} : drop
    shrimt1trin1nrcd1.blastn  : probably drop  
    shrimt1trin1nrcd1.cds.bak.clstr : always erase junk
    shrimt1trin1nrcd1.log : drop
    shrimt1trin1nrcd1.cds.clstr : probably drop; leave for inspection?
    shrimt1trin1.adupfilt.log : leave for inspection?
    shrimt1trin1.alntab : for inspect, reruns?
    
=cut


# =item OLD script tr2aacds_qsub.sh
# 
#   #! /bin/bash
#   ### env trset=allstrimpt1.tr datad=`pwd` qsub -q normal tr2aacds_qsub.sh
#   #PBS -N tr2cds
#   #PBS -l nodes=1:ppn=32,walltime=18:55:00
#   #PBS -o tr2cds.$$.out
#   #PBS -e tr2cds.$$.err
#   #PBS -V
# 
#   ncpu=30
#   maxmem=50000
#   evigene=$HOME/bio/evigene/scripts
#   #t2ac: app=cd-hit-est, path= echo MISSING_cd-hit-est
#   #v45bad32kseq#export PATH=$HOME/bio/cdhit/bin:$PATH
#   export PATH=$HOME/bio/cdhit461/bin:$PATH
#   #t2ac: app=fastanrdb, path= echo MISSING_fastanrdb
#   export fastanrdb=$HOME/bio/exonerate/bin/fastanrdb
#   #t2ac: app=blastn, path= echo MISSING_blastn
#   export PATH=$HOME/bio/ncbi2227/bin:$PATH
# 
#   if [ "X" = "X$trset" ]; then
#     echo "missing env trset=xxxx.tr"; exit -1
#   fi
#   if [ "X" = "X$datad" ]; then
#     echo "missing env datad=/path/to/data"; exit -1
#   fi
# 
#   cd $datad/
#   echo $evigene/prot/tr2aacds.pl -debug -NCPU $ncpu -MAXMEM $maxmem -log -cdna $trset
#   $evigene/prot/tr2aacds.pl -debug -NCPU $ncpu -MAXMEM $maxmem -log -cdna $trset
# 
# =cut

=item blastn fasplit speedup

speedier, same trclass results  
  c: 142 min total, 110 min blastn, blastn-threads
  d:  46 min total,  22 min blastn, blastn-fasplit, 5+x faster
            
catfish1all4.c : using -num_threads 30
#t2ac: BEGIN with cdnaseq= catfish1all4.tr.gz date= Wed Mar 13 09:50:37 PDT 2013
# Class Table for catfish1all4.trclass 
class           okay    drop    okay    drop
althi           5.5     19.5    46419   163572
althi1          0.8     5.8     6947    49022
althia2         0       5.1     0       43011
altmfrag        0.2     0.5     2419    4317
altmfraga2      0       0       291     417
altmid          0.6     0.9     5188    7669
altmida2        0       0       799     809
main            4.3     9.8     36322   82723
maina2          0.2     0.2     1852    2450
noclass         2.1     25.9    18326   217631
noclassa2       0       0       37      191
parthi          0       13.4    0       112502
parthi1         0       2.1     0       17972
parthia2        0       1.9     0       16530
---------------------------------------------
total           14.1    85.8    118600  718816
=============================================
# AA-quality for okay set of catfish1all4.aa.qual (no okalt): all and longest 1000 summary 
okay.top	 n=1000; average=2205; median=1921; min,max=1476,13376; sum=2205718; gaps=5109,5.1
okay.all	 n=56537; average=273; median=125; min,max=40,13376; sum=15437828; gaps=130303,2.3
#t2ac: asmdupfilter_fileset= catfish1all4.okay.tr catfish1all4.okalt.tr catfish1all4.drop.tr catfish1all4.okay.aa catfish1all4.okalt.aa catfish1all4.drop.aa catfish1all4.okay.cds catfish1all4.okalt.cds catfish1all4.drop.cds
#t2ac: tidyup output folders: okayset dropset inputset tmpfiles
#t2ac: DONE at date= Wed Mar 13 12:12:14 PDT 2013
#t2ac: ======================================

catfish1all4.d : using blastn-fasplit(query,30)
#t2ac: BEGIN with cdnaseq= catfish1all4.tr.gz date= Thu Mar 14 12:00:19 PDT 2013
# Class Table for catfish1all4.trclass 
class           okay    drop    okay    drop
althi           5.5     19.5    46420   163576
althi1          0.8     5.8     6947    49019
althia2         0       5.1     0       43012
altmfrag        0.2     0.5     2419    4318
altmfraga2      0       0       291     417
altmid          0.6     0.9     5188    7669
altmida2        0       0       799     809
main            4.3     9.8     36321   82724
maina2          0.2     0.2     1852    2450
noclass         2.1     25.9    18326   217632
noclassa2       0       0       37      191
parthi          0       13.4    0       112498
parthi1         0       2.1     0       17971
parthia2        0       1.9     0       16530
---------------------------------------------
total           14.1    85.8    118600  718816
=============================================
# AA-quality for okay set of catfish1all4.aa.qual (no okalt): all and longest 1000 summary 
okay.top	 n=1000; average=2205; median=1921; min,max=1476,13376; sum=2205718; gaps=5109,5.1
okay.all	 n=56536; average=273; median=125; min,max=40,13376; sum=15437776; gaps=130303,2.3
#t2ac: asmdupfilter_fileset= catfish1all4.okay.tr catfish1all4.okalt.tr catfish1all4.drop.tr catfish1all4.okay.aa catfish1all4.okalt.aa catfish1all4.drop.aa catfish1all4.okay.cds catfish1all4.okalt.cds catfish1all4.drop.cds
#t2ac: tidyup output folders: okayset dropset inputset tmpfiles
#t2ac: DONE at date= Thu Mar 14 12:46:10 PDT 2013
#t2ac: ======================================


=cut

