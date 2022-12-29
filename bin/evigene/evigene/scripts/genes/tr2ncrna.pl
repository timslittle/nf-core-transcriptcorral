#!/usr/bin/perl
# evigene/scripts/genes/tr2ncrna.pl
# parts from altreclass.pl, trclass_resolve_strandmix.pl 

=item about 

  EvidentialGene tr2ncrna 
    Collects putative non-coding RNA transcripts from Evigene tr2aacds results, 
    using replication across assemblies as major validation quality, as:
    
     1.  ncRNAredundant.tr from input.tr (all transcripts) minus okay.mrna (coding transcripts)
     2.  self-align ncRNAredundant.tr, measure replication from near-identicals across assemblies
     3.  classify to locus (main+alt+par and unique/noclass), select subset of representative per locus,
          i.e. keep long+replicated transcripts w/ exon-level variation (alt/par),
          drop near duplicate alternates, short transcripts, partial okay.mrna transcripts

=item usage     

  Usage: tr2ncrna.pl -trset input/name.tr -mrna okayset/name.okay.mrna
   options
  -ncpu 8 : use 8 cpu, cores, for parallel processes
  -log    : write progress to log file
  -min_ncnra=500 : minimum ncRNA transcript size 
  -minalign_ncrna=80 [1-99]: ignore ncrna % alignment below this, as likely paralog
  -reuse_selfblast : test option, 1 or 2 self-align blastn runs on tr subsets
  -updateall  : don't reuse intermediate results
  -debug      : more progress info
    
  Expects Evigene data, from SRA2Genes or tr2aacds, with inputset/dropset transcripts containing ncRNA, 
  and okayset mRNA sequences.  Output to ncrnaset/ with subsets of non-mRNA transcripts and classification
  tables.

=item UPD2222, update 2020.02.22 ..

  algorithm 2, from tests in human/sra2genes_tr19human/tr2aacds_test1908f/try20jaevg4f/evid/ncrna_altparscore.pl

  replace ncrna_replicates: reduce by agreement, not working, ie no good cor of trasm agreement w/ ref ncrna
  with ncrna_classgenes using rank of weighted evidence scores: trlen has most ref cor, 
    and scores of pCDS, codepot, agree
  
  use altparclass to group ncrna by locus, 
  score (im)perfect.dups,frags to drop as per alt.hi1 overabundant subset
  
  output table equivalent to mRNA pubids, with main/uni/alt/altpar classes, drop class of excess and low scores

=item  add pHETERO option as for tr2aacds
  
  need -pHetero=1..9 option to collapse locus/alt calls in presence of heterozygosity
  same as for coding seqs, 
  affects blastn -perc_ident filtering, and altparclass grouping
  $NCBLAST_IDENT/pctident -phet
  remove_mrna_aligned() -phet
  ncrna_align_okcds()
  altparclass2a()
   
=item test code run_evgncrnablself.sh

  1a. remove mrna oids from trset, assume trset=input.tr all transcripts
  1b. remove ~perfect dup+frag to okay.mrna using blastn align to okay.mrna

  
  2. align imperfect notokmrna.tr, calc total align and keep subset w/ poor align to okay.mrna
  2. blagree replacement: use only long-enough notoklong.tr, self-blast to find high-id long aligns, 
    then pick blagree subset replicated over assemblers
  2. FIXME: this way leaves in notokmrna with large overlap to okmrna : need to separate these from ncrna
    step1 removes only contained-in-okmrna subset
    ? use prior step2 + self-agree, add after step3 ? blastn -db uniqrna.tr -query ok.cds -qcov 60-90% ?

  3. altpar classify unique_ncrna.tr
    evigene/scripts/genes/altparclassify.pl -ncpu $ncpu -cds $trname.uniqrna.tr -sizes $trsizes
    
  3b. FIXME: step2 leaves in notokmrna with large overlap to okmrna : need to separate these from ncrna
     step1 removes only contained-in-okmrna subset
     use prior step2 + self-agree, add after step3 ? blastn -db uniqrna.tr -query ok.cds -qcov 60-90% ?

=cut

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/.."); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
use strict;
use Getopt::Long;
use cdna_evigenesub;  
# use evigene_pubsets;
use File::Basename qw(basename dirname fileparse);

use constant VERSION =>  '2020.02.25'; # v2c,d,e..
# '2020.02.10'; # v1a,b ..UPD1912
use constant UPD2222 => 1; # v2..
use constant UPD1912 => 1;
use constant FIXME => 0;
use constant USE_recover_trset => 1; # convert subset.tr <=> subset.trids, works reduces size lots

our $EVIGENES= $FindBin::Bin; $EVIGENES =~ s,scripts/.*,scripts,;
our $EGAPP='tr2ncrna';  
our $EGLOG='ncrna';

my $debug=$ENV{debug}||0; 
my $NCPU=$ENV{ncpu}||1; 
my $tidyup=1;
my $MIN_NCRNA= $ENV{MIN_NCRNA} || 300; # UPD2222 what size opt? 500 too big
my $MINAA_NCRNA= $ENV{MINAA_NCRNA} || 120; # save ncrna aatrans if aasize>min, should be same var as for tr2aacds etc
my $DROPSHORT_NCRNA= $ENV{DROPSHORT_NCRNA} || 0.20; #??
my $PAL_NCRNAMIN= $ENV{PAL_NCRNAMIN} || 80; # reject ncrna align below this 
my $TEST_OKCDS= $ENV{TEST_OKCDS}||0; # needs testing .. should become default
my $MIN_CDSALN= $ENV{MIN_CDSALN}||150; #? for TEST_OKCDS
my $MINBIGAA= $ENV{MINBIGAA} || 300; # remove_bigcdsdrops : 300? 200 ? or 150? 
my $UPDATEALL= 0;
my $REUSE_SELFBLAST= 0; # v2: default or no option..
my $usebself= 1; #v2 default  $ENV{TESTmrnablast} for remove_mrna_aligned() blast method 
my($IDPREncrna,$IDBASEncrna) =("",0); # UPD2222
my($premerge,@submerge)= (0);

# ncrna_altparscore globals
my $MINSCORE=$ENV{minscore} || -999; # need some data on this
my $MINALT= $ENV{minalt}|| 0.50; #??
my $MINDUPALN = $ENV{mindupaln} ||0.80; # ncrna_altparscore option: align(t/r)
my $MINDUPFULL= $ENV{mindupfull}||0.80; # ncrna_altparscore option: size(t/r)

my ($mrna,$trset,$output,$aaqualf,$mapsensetab,$logfile);
my @ARGSAVE=@ARGV;
my %opts= (
  "trset|input=s", \$trset, # require
  "mrna=s", \$mrna, # require, find 
  "logfile:s", \$logfile, # option
  "aaqualf=s", \$aaqualf, # aaqualf required ? w/ codepot; add mrna/cds/aa seq set??
  "MIN_NCRNA=i", \$MIN_NCRNA,  
  "MIN_CDSALN=i", \$MIN_CDSALN,  
  "minalign_ncrna=s", \$PAL_NCRNAMIN,  
  "NCPU=i", \$NCPU, # "MAXMEM=i", \$MAXMEM,  
  "premerge!", \$premerge, "mergesubsets=s", \@submerge, # UPD20mar
  "debug!", \$debug, 
  "UPDATEALL!", \$UPDATEALL, 
  "REUSE_SELFBLAST!", \$REUSE_SELFBLAST, 
  "TEST_OKCDS!", \$TEST_OKCDS, 
  );
# unused   "DROPSHORT_NCRNA=s", \$DROPSHORT_NCRNA,  
# unused  "output:s", \$output, # option NOT YET
  
my $optok= GetOptions(%opts);
my $optlist= ($debug) ? "All opts: ".join(", ",sort keys %opts) : "";
my $DOMERGE= (@submerge > 0)?1:0; #UPD20mar .. inputs not $trset, from submerge[0..n]/ncrnaset/xxx.longnomrna.tr

die "EvidentialGene tr2ncrna, VERSION ",VERSION,"
    Collects putative non-coding RNA transcripts from Evigene tr2aacds results, 
    classifies and ranks to gene locus by sequence qualities (size, non-coding potential,
    replication), removes duplicate and low-quality transcripts, for
    a rough draft ncRNA gene set.
  
Usage: tr2ncrna.pl -input myassembly.tr -mrna okayset/myassembly.okay.mrna 
 opts: 
  -ncpu 8 : use 8 cpu, cores, for parallel processes
  -log    : write progress to log file
  -min_ncnra=$MIN_NCRNA : minimum ncRNA transcript size 
  -minalign_ncrna=$PAL_NCRNAMIN [1-99]: ignore ncrna % alignment below this, as likely paralog
  -reuse_selfblast : test option, 1 or 2 self-align blastn runs on tr subsets
  -updateall  : don't reuse intermediate results
  -debug      : more progress info
  $optlist
  
  Expects Evigene data, from SRA2Genes or tr2aacds, with inputset/dropset transcripts containing ncRNA, 
  and okayset mRNA sequences.  Output to ncrnaset/ with subsets of non-mRNA transcripts and classification
  tables.
\n" unless($optok and $trset);  # and $mrna
# unused?   -dropshort_ncrna=$DROPSHORT_NCRNA [0.01-0.99]: drop ncrna alternate below proportion of longest, as too small

our $DEBUG= $debug;  # set package dbg
$DROPSHORT_NCRNA= $DROPSHORT_NCRNA/100 if($DROPSHORT_NCRNA >=1) ;

openloggit($logfile,basename($trset)); # inputset/name.tr > inputset/ncrna.log is bad
loggit(1, "EvidentialGene tr2ncrna, VERSION",VERSION);
loggit 1, "tr2ncrna  @ARGSAVE \n";

my ($APPblastn,$APPmakeblastdb,$APPfastanrdb);
   $APPfastanrdb= findapp("fastanrdb");  
   $APPblastn= findapp("blastn");  
   $APPmakeblastdb= findapp("makeblastdb");  

my ($okaysetd,$pubsetd); 
our (@ncrnaset, @okayset,@dropset,@inputset,@tmpfiles,@erasefiles,@publicset); # tidyup file sets
my $outputh=undef; # *STDOUT; # unused now

MAIN_tr2ncrna();

#================================================================

=item merge several subasm[123]/ncrnaset/ to merged/ncrnaset/

  my($ntr,$nsub,$mergetrset,$trsizes)= makeSubmergeTr($trname,\@submerge);

  usage?  -merge subasm1/,subasm2/,subasm3/ == paths to evg sub assemblies, w/ okayset/, ncrnaset/, other not needed?
          -premerge == for subasm[123]/ncrnaset/ production, means stop at step long_seqs()/trseq_qual()
             == or -runsteps=stop4|stop.trlongnomrna ?
             
  UPD20mar21: add -merge trlongnomrna[123].tr option for merging results of separate trasm.reduce
   .. modify here to start with collected trlongnomrna.tr inputs, and other needs
   .. trlongnomrna[123].tr should be from subasm[123]/ncrnaset/trlongnomrna.tr
   .. other in: subasm[123]/okayset/okay.{mrna,cds}, 
   .. subasm[123]/aaeval/okayaa_refaa.btall, OR merged/aaeval/okayaa_refaa.btall, and refset/refgenes.aa
   
   FIXME: want also inputset/ aaqualf, aa, cds of subasm matching oids of trlongnomrna subsets,
      for -merge outputs
      
=cut

sub filegz{ my($f)=@_; unless(-f $f){ $f="$f.gz" if(-f "$f.gz"); } return $f; }
sub step_log_or_fini{ my($fini,$msg)=@_; if($fini){ FINISH_tr2ncrna(-1,$msg); } else { loggit( 0, $msg); } }

sub MAIN_tr2ncrna {  # UPD2222 v2 algo
  loggit(0, "BEGIN with input=",$trset,"date=",`date`);
  my $trname=$trset; $trname =~ s/\.gz$//; $trname=~s/.\w+$//; $trname =~ s,^.*/,,;
  
  my @okpubs1a= qw( okayset publicset);
  my @okpubs2b= qw( okayset1st okayset);
  ($okaysetd,$pubsetd)= @okpubs2b; 
  unless( -d $okaysetd and -d $pubsetd ) {  # symlink dir bad here? maybe not, was missing okayset1st
    ($okaysetd,$pubsetd)= @okpubs1a;
    unless( -d $okaysetd and -d $pubsetd ) {
      loggit(LOG_DIE,"Cannot find evigene data directories:  @okpubs2b or @okpubs1a");
    }
  }
  
  # locate data associated w/ trclass
  unless($mrna) { 
    $mrna= filegz("$pubsetd/$trname.okay.mrna");  # mrna.gz ?? fixed remove_mrna_aligned
    }
  unless($aaqualf) { 
    $aaqualf= filegz("inputset/$trname.aa.qual");   # presume trset is all inputset/ or contained in that
    $aaqualf="$trname.aa.qual" unless(-f $aaqualf);
    }
  
  loggit(0,"tr2ncrna( $trset, $mrna)");
  loggit(LOG_DIE,"Cannot find mrna $mrna") unless(-f $mrna);
  # unless(ref $outputh) { $outputh=*STDOUT; }

  my($ntrnew1, $trnomrna, $ntrtotal, $mrnaoids)=(0) x 9;
  my($ntrnew2,$trnomrnab,$trsetagree)=(0) x 9;
  my($ntrnew3,$trlongnomrna,$trlongsizes)=(0) x 9; # DOMERGE sets these..
  
  if($DOMERGE) { # input is trlongnomrna, not trset
    my $nsubset=0;
    ($ntrnew3,$nsubset,$trlongnomrna,$trlongsizes)= makeSubmergeTr($trname,\@submerge);
    step_log_or_fini( ($ntrnew3 < 1),"makeSubmergeTr ntr=$ntrnew3 from $nsubset sets, $trlongnomrna");
    $ntrtotal=$ntrnew1=$ntrnew3; 
    
    # FIXME: replace aaqualf w/ trname_mrg.nclong.aa.qual ?? .. only for trseq_qual(trlong) ?
    (my $trlongnc_aaqual= $trlongnomrna) =~ s/\.\w+$/.aa.qual/; # makename($trlongnomrna,"aa.qual");
    $aaqualf=$trlongnc_aaqual if( -f $trlongnc_aaqual);
    
    # get trsetagree from makeSubmergeTr() trsets.nr ?
    my $trlongnc_agree= "$trlongnomrna.consensus";
    $trsetagree=$trlongnc_agree if( -f $trlongnc_agree);
   
  } else { # full pipeline
  
    ($ntrnew1,$trnomrna,$ntrtotal, $mrnaoids)= remove_mrna_oids($trname,$trset,$mrna); # .gz ok, v2OK
    step_log_or_fini( 0,"remove_mrna_oids kept=$ntrnew1/$ntrtotal, $trnomrna");
    
    # UPD20f15: insert step to remove dropset/perfectdup,perfectfrag where CDS size is large and match:ID is in okayset/
    my($ntrnew1b,$trnomrna1b)= remove_bigcdsdrops($trname, $trnomrna, $mrna, $mrnaoids); # .gz ok, v2OK
    step_log_or_fini( 0,"remove_bigcdsdrops kept=$ntrnew1b/$ntrtotal, $trnomrna1b"); 
    if($ntrnew1b > 0){ $trnomrna= $trnomrna1b; }
    
    ($ntrnew2,$trnomrnab,$trsetagree)= remove_mrna_aligned($trname,$trnomrna,$mrna); # .gz NOW ok yet, v2OK
    step_log_or_fini( 0,"remove_mrna_aligned kept=$ntrnew2/$ntrtotal, $trnomrnab");
    
    ($ntrnew3,$trlongnomrna,$trlongsizes)= long_seqs($trname, $trnomrnab, $MIN_NCRNA); # v2OK
    step_log_or_fini( ($ntrnew3<1), "long_seqs(>=$MIN_NCRNA) kept=$ntrnew3/$ntrtotal, $trlongnomrna");
    
    if($premerge) { 
      # FINISH special, not for premerge? .. no /done okay/
      #   if(USE_recover_trset and $status>=0 and $atstep =~ /done okay/);
      FINISH_tr2ncrna( 1, "done pre-merge of subset $trlongnomrna");
    }
    
  } # unless DOMERGE

  ## ncrna_replicates and altparclass meed trlongsizes.qual file; for new btall method, 
  my($ntrq,$trlongsizef)= trseq_qual($trname, $trlongnomrna, $trlongsizes, $aaqualf); #   v2OK
  
  #** NOTE:   ncrna_selfalign() and altparclass() both do blastn self-align, should merge them
  # but need new parsing for ncrna_replicates : btall blast table, pick out full aligns from that
  # ?? opt for altparclass -use_selfblast $selfblast
  # my $REUSE_SELFBLAST= $ENV{REUSE_SELFBLAST} || 0;

  ##UPD2222 changes here.........
  my($nokcds, $hasokcds) =(0,{}); # not yet tested, dont know value yet
  my($ntrnew4, $trselfblast, $ntrnew5, $truniqrna)=(0) x 9;
  my($ngene,$genetab)=(0,"");
  
if(UPD2222) { # ready

    $REUSE_SELFBLAST= 1; # v2 always
    ($ntrnew5, $truniqrna)= ($ntrnew3,$trlongnomrna); # this will change? truniq = remove drops of classgene
    
    ($ntrnew4, $trselfblast)= ncrna_selfalign($trname, $trlongnomrna, $REUSE_SELFBLAST);
    step_log_or_fini( ($ntrnew4<1),"ncrna_selfalign count=$ntrnew4/$ntrtotal, $trselfblast");

    ($ngene,$genetab)= altparclass2a($trname, $trlongnomrna, $trselfblast, $trlongsizef); # ?? , $hasokcds
    step_log_or_fini(  ($ngene<1), "altparclass loci=$ngene for $ntrnew4 tr, in $genetab");

    #v2: dont need altpar genetab ?? >> needed at make_allevdtab()
    #v2upd: altparclass2a() and makeSelfDupTab() now both read self.btall, do related things, could be merged
    
    my($ntab, $selfduptab, $havedups, $hasduphash, $nagree)
        = ncrna_makeSelfDupTab($trname, $trlongnomrna, $trselfblast, $trsetagree, $trlongsizef);
    step_log_or_fini( 0, "makeSelfDupTab ndups=$havedups, nagree=$nagree");
    # ret: nagree,%blagree, ndup, %hasdup ??

    ## add here? class b/n truniq matching okay.cds, truniq.okcds.tr, and those w/o truniq.ncrna.tr
    $hasokcds= undef;
    (my $okcds=$mrna) =~ s/.mrna/.cds/; # allow cds.gz ?
    if($TEST_OKCDS and -f $okcds) {  # should become default if -f okcds
      ($nokcds, $hasokcds) = ncrna_align_okcds($trname,$trlongnomrna,$okcds);
      step_log_or_fini( 0, "ncrna_align_okcds hasokcds=$nokcds/$ntrnew5");
    } 
    $TEST_OKCDS=0 unless($nokcds>0); # turn off this flag for follow ons
    
    #v2: insert here .. want altpar genetab, blagree/selfdup table, aa.qual + fickett codepot of trset
    # ?? add $hasduphash as allevdtab col ? rather than thru classgenes()
    #    ^^ not good, need to add duptrids as column, classgenes ranks by quals then keeps 1st of dup set
    
    my($nevd,$allevdtab)= make_allevdtab( $trname, $trlongnomrna, $genetab, $trlongsizef,  $selfduptab, $hasokcds );
    step_log_or_fini(  ($nevd<1), "make_allevdtab nevd=$nevd/$ntrnew5 in $allevdtab");

    #v2: final step? ncrna_classgene(@dgenetrs), where @dgenetrs are all trset in gene group, w/ evidence
    # classifier weights evidence for rank by trscore, keeps/drops alt/par of gene group depending on scores
    # output table in publicset/pubids format? w/ class 'main/alt/noclass=uni' and 'drop|cull' qualifier
    # ncrna_classgenes: read evdtab, each tr sorted by genetab dgene ids, with evidence    
    
    my($ngeneok,$ntrok,$ntrdrop,$classtab)= ncrna_classgenes($trname, $trlongnomrna, $genetab, $allevdtab, $hasduphash);
    step_log_or_fini( 0, "ncrna_classgenes ngene=$ngeneok,ok=$ntrok,drop=$ntrdrop/$ntrnew5 in $classtab");
    
    # update ($ntrnew5, $truniqrna) removing drops in classtab
    my($ntrnew5, $truniqrna, $truniqncaa)= ncrna_pubseqset($trname,$trlongnomrna,$classtab);  
    step_log_or_fini( 0, "ncrna_pubseqset kept=$ntrnew5/$ntrtotal, in $truniqrna"); 

    # Maybe insert here ncaablastp( ncrnaset/nc.aa, refaa, okayset.aa.btall) matchup to find uniq ref prots
    # .. or in sra2genes .. need from that refaa, okayset.aa.btall 
    my($aok,$refaa,$okbtall)= ($truniqncaa) ? has_okayset_refblast($trname) : (0); 
    if($aok) {
      my($ntab,$ncaareftab,$ncuniq)= ncrna_aauniq_refblast($trname, $truniqncaa, $refaa,$okbtall);
      ## ntab == -1 on error, needs REFAA, REF_OKAA_BTALL to work
      step_log_or_fini( 0, "ncrna_aauniq_refblast refuniq=$ncuniq/$ntab, in $ncaareftab");
      # output now has  trname.ref_ncaa_okaa.uniqtab , want to pull ncaa/cds/mrna of those to merge w/ pubset?
    }
    
} else {

  my $skip_replicate_filter= $ENV{skip_replicate_filter}||0; # TEST altparclass loci instead, filt by longest/locus
  
  if($skip_replicate_filter) {
    ($ntrnew5, $truniqrna)= ($ntrnew3,$trlongnomrna);
  } else {
    
    ($ntrnew4, $trselfblast)= ncrna_selfalign($trname, $trlongnomrna, $REUSE_SELFBLAST);
    step_log_or_fini( ($ntrnew4<1),"ncrna_selfalign count=$ntrnew4/$ntrtotal, $trselfblast");

    ($ntrnew5, $truniqrna)  = ncrna_replicates($trname,$trselfblast,$trsetagree,$trlongnomrna,$trlongsizes);
    step_log_or_fini( ($ntrnew5<1),"ncrna_replicates kept=$ntrnew5/$ntrtotal, $truniqrna");
  }
  
  ## add here? class b/n truniq matching okay.cds, truniq.okcds.tr, and those w/o truniq.ncrna.tr
  (my $okcds=$mrna) =~ s/.mrna/.cds/;
  if($TEST_OKCDS and -f $okcds) { 
    ($nokcds, $hasokcds) = ncrna_align_okcds($trname,$truniqrna,$okcds);
    step_log_or_fini( 0, "ncrna_align_okcds hasokcds=$nokcds/$ntrnew5");
  }
  
  ($ngene,$genetab)= altparclass($trname, $truniqrna, ($REUSE_SELFBLAST ? $trselfblast :""), $trlongsizef, $hasokcds);
  step_log_or_fini( 0, "altparclass loci=$ngene/$ntrtotal, $genetab");
}
  
  # close($outputh) if($output);
  step_log_or_fini( 0, "tr2ncrna nin=$ntrtotal, nout=$ntrnew5, ngene=$ngene\n");
  
  ## no rewriteTrClass(), but want some ncrna.pubids table output
  FINISH_tr2ncrna(0,'done okay',$trname);
}

#=======================================================
  
# sub MAIN_tr2ncrna_algo1b { ## OLD
#   loggit(0, "BEGIN with input=",$trset,"date=",`date`);
# 
#   my $trname=$trset; $trname =~ s/\.gz$//; $trname=~s/.\w+$//; $trname =~ s,^.*/,,;
#   # my($trname,$trpath,$trsuffix) = fileparse($trset, qr/\.\w*/); # suf-'.trclass' or suf = qr/\.\w*/
#   # chdir($trpath) unless($trpath eq './' or $trpath eq "");
#   
#   my @okpubs1a= qw( okayset publicset);
#   my @okpubs2b= qw( okayset1st okayset);
#   ($okaysetd,$pubsetd)= @okpubs2b; 
#   unless( -d $okaysetd and -d $pubsetd ) {
#     ($okaysetd,$pubsetd)= @okpubs1a;
#     unless( -d $okaysetd and -d $pubsetd ) {
#       loggit(LOG_DIE,"Cannot find evigene data directories: @okpubs1a");
#     }
#   }
#   
#   # locate data associated w/ trclass
#   unless($mrna) { 
#     $mrna="$pubsetd/$trname.okay.mrna";  # mrna.gz ?? fix remove_mrna_aligned
#     }
#   unless($aaqualf) { 
#     # $aaqualf="$okaysetd/$trname.aa.qual"; # which 1st? unless(-f $aaqualf);
#     $aaqualf="inputset/$trname.aa.qual";  # and .gz ??
#     $aaqualf="inputset/$trname.aa.qual.gz"  unless(-f $aaqualf);
#     $aaqualf="$trname.aa.qual" unless(-f $aaqualf);
#     }
#   loggit(LOG_DIE,"Cannot find mrna $mrna") unless(-f $mrna);
# 
# #   my $makeoutpubids= (defined $output) or (not $debug); # always print to file now, unless testing?
# #   if($makeoutpubids and not $output) {
# #     ($output=$pubids) =~ s/\.pubids.*//; $output.=".reorient_pubids";
# #   }
# #   if($output) {
# #     rename($output,"$output.old") if(-f $output);
# #     open($outputh,'>',$output) or die "ERR: writing $output"; 
# #   }
# 
# #? want this? wont have all input trset as many dropped before trclass  
# #  ($trinfo,$aaqualhash)= readTrClass($trclass);  # require trclass here
# #? want this? should include trlen, aalen for all or most input trset  
# #  ($aaqualhash)= getAaQual($aaqualf,$aaqualhash); # cdna_evigenesub;
# #  my($nmapqual, $mapqualh, $alntabh)= ( -f $mapsensetab ) ? readAlignTab($mapsensetab) : (0,0,0);
#   
#   sub step_log_or_fini{ my($fini,$msg)=@_; if($fini){ FINISH_tr2ncrna(-1,$msg); } else { loggit( 0, $msg); } }
#   
#   loggit(0,"tr2ncrna( $trset, $mrna)");
#   # unless(ref $outputh) { $outputh=*STDOUT; }
# 
#   my($ntrnew1,$trnomrna,$ntrtotal, $mrnaoids)= remove_mrna_oids($trname,$trset,$mrna); # .gz ok
#   step_log_or_fini( 0,"remove_mrna_oids kept=$ntrnew1/$ntrtotal, $trnomrna");
#   
#   # UPD20f15: insert step to remove dropset/perfectdup,perfectfrag where CDS size is large and match:ID is in okayset/
#   my($ntrnew1b,$trnomrna1b)= remove_bigcdsdrops($trname, $trnomrna, $mrna, $mrnaoids); # .gz ok
#   step_log_or_fini( 0,"remove_bigcdsdrops kept=$ntrnew1b/$ntrtotal, $trnomrna1b"); 
#   if($ntrnew1b > 0){ $trnomrna= $trnomrna1b; }
#   
#   my($ntrnew2,$trnomrnab,$trsetagree)= remove_mrna_aligned($trname,$trnomrna,$mrna); # .gz NOT ok yet
#   step_log_or_fini( 0,"remove_mrna_aligned kept=$ntrnew2/$ntrtotal, $trnomrnab");
#   
#   my($ntrnew3,$trlongnomrna,$trlongsizes)= long_seqs($trname, $trnomrnab, $MIN_NCRNA);
#   step_log_or_fini( ($ntrnew3<1), "long_seqs(>=$MIN_NCRNA) kept=$ntrnew3/$ntrtotal, $trlongnomrna");
#   # #ncrna: FAILED at step: long_seqs(>=500) kept=/1080380, plJBND.notokmrna.tr
# 
#   ## ncrna_replicates and altparclass meed trlongsizes.qual file; for new btall method, 
#   my($ntrq,$trlongsizef)= trseq_qual($trname, $trlongnomrna, $trlongsizes, $aaqualf);
#   
#   #** NOTE:   ncrna_selfalign() and altparclass() both do blastn self-align, should merge them
#   # but need new parsing for ncrna_replicates : btall blast table, pick out full aligns from that
#   # ?? opt for altparclass -use_selfblast $selfblast
#   # my $REUSE_SELFBLAST= $ENV{REUSE_SELFBLAST} || 0;
# 
#   my $skip_replicate_filter= $ENV{skip_replicate_filter}||0; # TEST altparclass loci instead, filt by longest/locus
#   
#   my($ntrnew4, $trselfblast, $ntrnew5, $truniqrna)=(0) x 9;;
#   if($skip_replicate_filter) {
#     ($ntrnew5, $truniqrna)= ($ntrnew3,$trlongnomrna);
#   } else {
#     
#     ($ntrnew4, $trselfblast)= ncrna_selfalign($trname, $trlongnomrna, $REUSE_SELFBLAST);
#     step_log_or_fini( ($ntrnew4<1),"ncrna_selfalign count=$ntrnew4/$ntrtotal, $trselfblast");
#   
#     ($ntrnew5, $truniqrna)  = ncrna_replicates($trname,$trselfblast,$trsetagree,$trlongnomrna,$trlongsizes);
#     step_log_or_fini( ($ntrnew5<1),"ncrna_replicates kept=$ntrnew5/$ntrtotal, $truniqrna");
#   }
#   
#   ## add here? class b/n truniq matching okay.cds, truniq.okcds.tr, and those w/o truniq.ncrna.tr
#   my($nokcds, $hasokcds) =(0,{});
#   (my $okcds=$mrna) =~ s/.mrna/.cds/;
#   if($TEST_OKCDS and -f $okcds) { 
#     ($nokcds, $hasokcds) = ncrna_align_okcds($trname,$truniqrna,$okcds);
#     step_log_or_fini( 0, "ncrna_align_okcds hasokcds=$nokcds/$ntrnew5");
#   }
#   
#   my($ngene,$genetab)= altparclass($trname, $truniqrna, ($REUSE_SELFBLAST ? $trselfblast :""), $trlongsizef, $hasokcds);
#   step_log_or_fini( 0, "altparclass loci=$ngene/$ntrtotal, $genetab");
#   
#   # close($outputh) if($output);
#   step_log_or_fini( 0, "tr2ncrna nin=$ntrtotal, nout=$ntrnew5, ngene=$ngene\n");
#   
#   ## no rewriteTrClass(), but want some ncrna.pubids table output
#   FINISH_tr2ncrna(0,'done okay');
# }

sub FINISH_tr2ncrna {
  my($status,$atstep,$trname)= @_;
  
  # my $EVG_ncrnaset = 'ncrnaset';
  my $tmpfolder="ncrnaset";
  my $outfolder= $okaysetd; # source files are here, need to merge reorfiles still
  loggit( LOG_WARN,"FAILED at step:",$atstep) if($status<0);
  
  if( $tidyup ) {  #  and $nupdate > 0 
    tidyupFileset($outfolder,@ncrnaset) if(@ncrnaset);  # not used so far
    tidyupFileset($tmpfolder,@tmpfiles) if(@tmpfiles);  # all data here so far ncrnaset/
    my @rmlist;
    if($status<0) { 
      tidyupFileset("ncrnaset_failtmp",@erasefiles) if(@erasefiles);  ## dont erase
    } else {
      foreach my $fn (@erasefiles) { if(-f $fn) { unlink($fn); push @rmlist,$fn; } }  
      if(@rmlist) { my $nrm=@rmlist; my $rml=join" ",@rmlist[0..4]; loggit(0,"tidyup erase: n=$nrm, $rml .."); }
    } 
    
    #UPD2222: cleanup should include turning all name.parts.tr into name.parts.trids list of ids only ..
    if(USE_recover_trset and $status>=0 and $atstep =~ /done okay/) {
      my($tdir,@trs)= getFileset($tmpfolder,'tr$','',$trname);
      for my $tr (@trs) {
        (my $trid= $tr) =~ s/\.\w+$/.trids/;
        ($trid)= faidlist($tr,$trid,'update');          
        unlink($tr) if(-s $trid);  # push @rmlist,$tr;
      }
    }
    
  }
  loggit(0, "DONE at date=",`date`);
  loggit(0, "======================================"); # log may append runs
  exit($status);
}
#------------------------------


sub recover_trset {
  my($trname,$trset,$trsetnew)= @_;

  # also maybe look in subdir ncrnaset/ 
  $trsetnew= "ncrnaset/$trsetnew" if( ! -f $trsetnew and -f "ncrnaset/$trsetnew");
  if( -s $trsetnew ) {
    return( facount($trsetnew), $trsetnew); 
  } else {
    (my $trnewids= $trsetnew) =~ s/\.\w+$/.trids/;
    $trnewids= "ncrnaset/$trnewids" if( ! -f $trnewids and -f "ncrnaset/$trnewids");
    if(-s $trnewids) {
      my($newfa,$ntrnew)= faextract($trset,$trsetnew,$trnewids);
      return($ntrnew, $newfa); 
    }
  }
  return(0,$trsetnew);
}


=item long_seqs algo

  my($ntrnew,$trsetnew,$trsizes)= long_seqs($trname,$trset,$minsize);
  my($ntrnew3,$trlongnomrna,$trlongsizes)= long_seqs($trname, $trnomrnab, $MIN_NCRNA);
  
  env ismrna=1 $evigene/scripts/prot/aaqual.sh $trname.notokmrna.tr
  perl -ne 'BEGIN{ $MINW=$ENV{minw}||500; } if(/^>(\S+)/){ $ok=$ok{$1}; $infa=1; print if $ok; }
  elsif($infa){ print if $ok; } else { ($d,$w)=split; $ok{$d}=1 if($w>=$MINW); } ' \
    $trname.notokmrna.tr.qual $trname.notokmrna.tr > $trname.longnok.tr
  
  trsetlong=$trname.longnok.tr

=cut

sub long_seqs {
  my($trname,$trset,$minsize)= @_;
  my $trsetnew="$trname.longnok.tr"; push @tmpfiles,$trsetnew;
  my $ntrnew= 0;
  my($cmd,$err,$ok,$hin);

if(USE_recover_trset) {  
  unless($UPDATEALL) {
    ($ntrnew,$trsetnew)= recover_trset($trname,$trset,$trsetnew); # handles .trids to .tr recovery
    if($ntrnew > 0) { 
      my $trlongs= fasizes_nogap($trsetnew,"","",1); #  1=only.nt or $fasizes{$id}= join"\t",$nokay,$nt,$ngap; 
      return($ntrnew, $trsetnew, $trlongs); 
    }
  }
} else {  
  if( -s $trsetnew and not $UPDATEALL) {
    my $trlongs= fasizes_nogap($trsetnew,"","",1); #  1=only.nt or $fasizes{$id}= join"\t",$nokay,$nt,$ngap; 
    $ntrnew= scalar(keys %$trlongs);
    return($ntrnew,$trsetnew,$trlongs);
  }
}
  
  # #ncrna: FAILED at step: long_seqs(>=500) kept=/1080380, plJBND.notokmrna.tr

  # my $trsizes= fasizes_nogap($trset); #  $fasizes{$id}= join"\t",$nokay,$nt,$ngap; 
  my $trsizes= fasizes_nogap($trset,"","",1); # only.nt
  my %trlongs=(); for my $d (keys %$trsizes) { 
    my $nt= $trsizes->{$d};  $trlongs{$d}= $nt if($nt >= $minsize); 
    #  my($nok,$nt)= split "\t",$trs;  #?? nt or all trs ?
  } 
  ($trsetnew,$ntrnew)= faextract($trset,$trsetnew,\%trlongs); # grep trsizes >= $MIN_NCRNA
  # NOTE faextract() returns (trset) if no ids in trlongs
  
  return($ntrnew,$trsetnew,\%trlongs);
}

sub makeSubmergeTr {
  my($trname,$submergeh)=@_; 
  
  my($cmd,$err,$ok,$hin);
  my($nsubs,$ntr,$trsizes)=(0,0);
  # my $submergetr="$trname.longnok.tr";
  my $submergetr= $trname . "_mrg.longnok.tr"; #?

#   FIXME here? need also inputset/ aaqualf, aa, cds of subasm matching oids of trlongnomrna subsets,
#     for -merge outputs at pubseqset, refblastp
  
  if( -s $submergetr  and not $UPDATEALL ) { # recover
    $trsizes= fasizes_nogap($submergetr,"","",1); # only.nt
    $nsubs=1; $ntr= scalar(keys %$trsizes);
    return($ntr,$nsubs,$submergetr,$trsizes);
  }
  
  # @$submergeh == global @submerge
  # @submerge ==  @$submergeh from -merge subasm1,subasm2 -merge subasm3,subasm4
  # likely mixed set of paths to subset assemblies w/ ncrnaset/
  
  my @smd= map{ split/[,;\s]+/, $_ } @$submergeh; # @submerge;
  my @subtr;
  for my $smd (@smd) {  
    # my($smnctr)= finddata("$sm/ncrnaset/*.longnok.tr"); #? is .tr.gz or .trids allowed ?
    # $sm= "$sm/ncrnaset" unless($sm=~/ncrna/);
    my $sm=$smd; $sm= "$sm/ncrnaset" if( -d "$sm/ncrnaset");
    my($smdh,$smnctr)= getFileset($sm,'(tr|tr.gz)$',undef,'longnok');
    push @subtr,$smnctr if(-f $smnctr);
  }
  $nsubs= @subtr;
  return($ntr,$nsubs,$submergetr) if($nsubs<1);
  
  my(%trsizes); # ,%stroiddont need stroid? just trsizes{id} : replace %trsizes w/ $trsizes=fasizes(subnrtr)
  open(SF,'>',$submergetr);  push @tmpfiles, $submergetr;    
  for my $str (@subtr) {
    ($ok,$hin)= openRead($str);
    if($ok) { 
      my($id,$len)=(0,0);
      while(<$hin>) { 
        print SF $_; 
        if(/^>(\S+)/) { 
          $trsizes{$id}=$len if($id);  
          $id=$1; $len=0; $ntr++; 
        } else { chomp; $len += length($_); }
      } 
      close($hin); $trsizes{$id}=$len if($id);
      loggit(LOG_DEBUG,"merge $str, nt=$ntr");
    } 
  } close(SF);

  ## ?? Insert here fastanrdb submerge.tr > submerge.nrtr; mv submerge.nrtr submerge.tr;
  ## premerge longnok.tr should be nr reduced. Maybe not enough dups left to bother?
  # sub67mrg.longnok.tr  921051, tr.nr 921002 = 49 shared in sub6,7
  # sub47mrg.longnok.tr 1467462, tr.nr  1467462 = 0 shared

  my $trsizes= \%trsizes;  # use this
  my $subtragree="";
  if(1) { # for pipe data, subtr is from fastanrdb w/ dupids in header, so this will work
    $cmd="$EVIGENES/prot/make_consensus.pl $submergetr";
    $err= runcmd($cmd);
    unless($err){ $subtragree="$submergetr.consensus"; push @tmpfiles, $subtragree if(-f $subtragree); }
  }
  
  if(0) {
    $cmd="$APPfastanrdb -i -f $submergetr > $submergetr.nr";
    $err= runcmd($cmd);
    unless($err) { 
      rename("$submergetr.nr",$submergetr);   # orig trset is tmpfile
      $trsizes= fasizes_nogap($submergetr,"","",1); # only.nt
      $ntr= scalar(keys %$trsizes);
      loggit(LOG_DEBUG,"merged to $submergetr, nr.nt=$ntr");
      
      ## add? make_consensus.pl  .. makec output is tr.consensus
      $cmd="$EVIGENES/prot/make_consensus.pl $submergetr";
      $err= runcmd($cmd);
      $subtragree="$submergetr.consensus"; push @tmpfiles, $subtragree;
    }
  }
  
  ## make also _mrg.longnok.aa,cds,aa.qual .. messy
  my @suba=();
  for my $smd (@smd) {  
    my $sm=$smd; $sm= "$sm/inputset" if( -d "$sm/inputset"); # problem if this isnt 
    my($smdh,$saa,$scds,$saaq);
    ($smdh,$saa)= getFileset($sm,'(aa|aa.gz)$',undef); # NOT ,$trname
    ($smdh,$scds)= getFileset($sm,'(cds|cds.gz)$',$smdh); # ,$trname
    ($smdh,$saaq)= getFileset($sm,'(aa.qual|aa.qual.gz)$',$smdh); # ,$trname
    push @suba, [$saa,$scds,$saaq] if(-f $saa or -f $scds or -f $saaq);
  }
  
  if(@suba > 0) { 
    my %outs; # openwrite _mrg.longnok.aa,cds,aa.qual
    
    sub puto { 
      my($oh,$sin,$trsizes,$istab)=@_; my $n=0; 
      my($ok,$hin)= openRead($sin); return($n) unless($ok);
      while(<$hin>) {  
        if($istab){ my($d)=split; $ok=$trsizes->{$d}||0; $n++ if($ok); } 
        elsif(/^>(\S+)/){ my $d=$1; $ok=$trsizes->{$d}||0; $n++ if($ok); } 
        print $oh $_ if($ok); 
        } 
      close($hin); return($n); 
      }
    
    my($naa,$ncds,$naq)=(0) x 3;
    for my $sf ('aa','cds','aa.qual') {
      (my $onam=$submergetr) =~ s/.\w+$//; $onam.=".$sf";
      open(my $oh,'>',$onam); $outs{$sf}= $oh; push @tmpfiles,$onam; # or erasefiles ??
    }
    for my $suba (@suba) { 
      my($saa,$scds,$saaq)= @$suba;
      $naa += puto($outs{'aa'},$saa,$trsizes,0);  
      $ncds+= puto($outs{'cds'},$scds,$trsizes,0);  
      $naq += puto($outs{'aa.qual'},$saaq,$trsizes,1);  
    } 
    for my $sf (keys %outs) { close( $outs{$sf} ); }
    loggit(LOG_DEBUG,"merge $submergetr, tr=$ntr,aa=$naa,cds=$ncds,aaq=$naq");
  }
  
  return($ntr,$nsubs,$submergetr,$trsizes);
}


#  my($ntrq,$trlongsizef)= trseq_qual($trname, $trlongnomrna, $trlongsizes, $aaqualf);
sub trseq_qual {
  my($trname,$trset, $trsizeh, $aaqualf)= @_;
  my($nts,%aaqual)=(0);
  my $trsizef= "$trset.qual"; 
  if(-s $trsizef and not $UPDATEALL) {
    open(F,$trsizef); while(<F>) { $nts++ if(/^\w/); } close(F); # $nts=`wc -l $trsizef`; chomp($nts);
    return ($nts, $trsizef);
  }
  
  ## BUG fixed .. got all zero $nok? was BAD tr/$pat/$pat/ in fasizes_nogap
  unless(ref $trsizeh) { $trsizeh= fasizes_nogap($trset); } # was BUG HERE: $fasizes{$id}= join"\t",$nokay,$nt,$ngap; 

  if( -f $aaqualf ) {
    my($ok,$hin)= openRead($aaqualf);
    if($ok){ while(<$hin>){ my @v=split; my($id)=$v[0]; if($trsizeh->{$id}){ $aaqual{$id}=\@v;} } close($hin); }
  } 
  open(SF,'>',$trsizef);  push @tmpfiles, $trsizef;  
  foreach my $id (sort keys %$trsizeh) {
    my $ts= $trsizeh->{$id};  my $aq= $aaqual{$id};
    my ($tnok,$tgap,$aaq,$tw,$ofsocp,$ofsicp)=(0) x 9;
    if(ref $aq) { ($tnok,$tgap,$aaq,$tw,$ofsocp,$ofsicp)= @{$aq}[4,2,3,4,5,6]; } # my($id,$aw,$agap,$aaq,$tw,$ofsOrCP,$ofsIfCP)=@v;
    if($ts =~ /\s/){ ($tnok,$tw,$tgap)= split" ",$ts; }
    elsif($ts=~/^\d+$/){ $tnok= $tw= $ts; } #?? fasizes_nogap(,,,1)
    print SF join("\t",$id,$tnok,$tgap,$aaq,$tw,$ofsocp,$ofsicp)."\n"; $nts++;
  }
  close(SF); 
  
  return ($nts,$trsizef);
}


=item remove_mrna_oids algo

  ($ntrnew,$trsetnew)= remove_mrna_oids($trname,$trset,$mrna);
  
  trsetnew=$trname.nomrna.tr
  grep '>' $mrna | perl -ne '($id)=m/>(\S+)/; ($od)= m/oid=([^;,\s]+)/ ? $1 :$id; print "$od\n";' > $mrna.oids
  perl -ne 'if(/^>(\S+)/){ $ok=not $nok{$1}; $infa=1; print if $ok; } elsif($infa){ print if $ok; } elsif(/^\w/){ ($d)=split; $nok{$d}=1; } ' \
   $mrna.oids $trset > $trsetnew

=cut

sub remove_mrna_oids { # STEP1a_
  my($trname,$trset,$mrna)= @_;
  
  my $trsetnew="$trname.nomrna.tr"; push @tmpfiles,$trsetnew;
  my $ntrnew= 0; my $ntrcut=0;
  my ($ok,$hin); 
  
## FIXME: not UPDATEALL here misses IDPREncrna setting from mrna .. recover where?
if(USE_recover_trset) {  
  unless($UPDATEALL) {
    ($ntrnew,$trsetnew)= recover_trset($trname,$trset,$trsetnew); # handles .trids to .tr recovery
    if($ntrnew > 0) { return($ntrnew, $trsetnew, facount($trset)); }
  }
} else {  
  if( -s $trsetnew and not $UPDATEALL) {
    $ntrnew= facount($trsetnew); my $ntotal=facount($trset);
    return($ntrnew, $trsetnew, $ntotal);
  }
}
  
  my(%mrnaoids,%mrnagid,%mrnapre,%mrnaid);  
  ($ok,$hin)= openRead($mrna); return(0) unless($ok);  # .gz ok
  while(<$hin>) { 
    if(/^>(\S+)/){ 
      my $id=$1; my($ods)= m/oid=([^;,\s]+)/;  
      map{ $mrnaoids{$_}= $id; } split(",","$id,$ods");
       
      if($id =~ m/^(\w+\D)(\d+)t(\d+)$/) {
        my($idpre,$gid,$ti)= ($1,$2,$3); 
        $mrnapre{$idpre}++; $mrnagid{$gid}++; $mrnaid{$gid."t$ti"}++;      
      } else {
        $mrnapre{$id}++; # fail evg pubid
      }
      } 
  } close($hin); 

  if(scalar(keys %mrnapre) == 1) { # set pubid vals
    my @gid= sort{ $a <=> $b } keys %mrnagid;
    my $gmax= $gid[-1];
    ($IDPREncrna)= keys %mrnapre;
    $IDBASEncrna= 1000 + 1000*int(0.5 + $gmax/1000);
    loggit(0,"ncRNA public ID: $IDPREncrna,$IDBASEncrna");
  }
  
  #? ($trsetnew,$ntrnew)= faextract($trset,$trsetnew,\%mrnaoids,'notid');
  ($ok,$hin)= openRead($trset); return(0) unless($ok);
  open(TO,'>',$trsetnew) or return(0); 
  $ok=0; 
  while(<$hin>) { 
    if(/^>(\S+)/){ $ok= $mrnaoids{$1}?0:1; 
      if($ok and my($ods)= m/oid=([^;,\s]+)/) {
        $ok=0 if(grep { $mrnaoids{$_} } split(",",$ods)); 
        }
      if($ok){ $ntrnew++; } else { $ntrcut++; }
      } 
    print TO $_ if($ok);
  } close($hin); close(TO); 

  return($ntrnew, $trsetnew, $ntrnew+$ntrcut, \%mrnaoids);
}

=item remove_bigcdsdrops

   * this works well
  
   UPD20f15: insert step to remove dropset/perfectdup,perfectfrag where CDS size is large and match:ID is in okayset/
   my($ntrnew1b,$trnomrna1b)= remove_bigcdsdrops($trname,$trnomrna,$mrna); # .gz ok
   step_log_or_fini( 0,"remove_bigcdsdrops kept=$ntrnew1b/$ntrtotal, $trnomrna1b"); 

  eg. dropset/plHILW.drop.cds.gz
  >Acmpantridba1a_sBn1l1ERR2040842idbaidbtk27Loc10 type=CDS; aalen=3783,92%,complete; clen=12237; strand=+; offs=406-11757; codepot=Code/0.15; 
     evgclass=perfectdup,drop,match:Acmpantridba1a_sBn1l1ERR2040842idbaidbtk27Loc1
  >Acmpantridba1a_sBn1l1ERR2040842idbaidbtk27Loc21 type=CDS; aalen=3563,97%,partial3; clen=10962; strand=-; offs=10689-1; codepot=Code/0.11; 
     evgclass=perfectfrag,drop,match:Acmpantridba1a_sBn1l1ERR2040842idbaidbtk67Loc35

=cut

sub remove_bigcdsdrops {
  my($trname,$trset,$mrna, $mrnaoids)= @_;  
  # reuse %mrnaoids from remove_mrna_oids()
  
  my $trsetnew="$trname.nodropbigcds.tr";  push @tmpfiles,$trsetnew;
  my $ntrnew= 0;
  my($cmd,$err,$ok,$hin);

if(USE_recover_trset) {
  unless($UPDATEALL) {
    ($ntrnew,$trsetnew)= recover_trset($trname,$trset,$trsetnew); # handles .trids to .tr recovery
    if($ntrnew > 0) { return($ntrnew, $trsetnew); }
    }
} else {  
  if( -s $trsetnew and not $UPDATEALL) {
    return(facount($trsetnew),$trsetnew);
  }
}

  # my $dropcds="dropset/$trname.drop.cds";  
  #    $dropcds="dropset/$trname.drop.cds.gz"  unless(-f $dropcds);
  my $dropcds= filegz("dropset/$trname.drop.cds");  
  return(0) unless(-s $dropcds and $mrnaoids and ref($mrnaoids) );

  my($nbigdrop,$ntrcut,%bigdropcds)=(0,0);  
  ($ok,$hin)= openRead($dropcds);  # .gz ok
  if($ok) { while(<$hin>){ 
    if(/^>(\S+)/){ my $oid=$1; 
      my($evc)= m/evgclass=([^;\s]+)/;  
      next unless($evc=~/(perfectdup|perfectfrag)/);
      
      my($aaw,$pcds)= (m/aalen=(\d+),(\d+)/) ? ($1,$2) : (m/aalen=(\d+)/)? ($1,0) : (0,0); 
      next unless($aaw > $MINBIGAA or $pcds > 65); # what %cds?
      next if($aaw < 150 and $pcds < 90); #?? minaa what?
      #x next if($aaw < $MINBIGAA and $pcds < 66); # what %cds?
      #x next if($aaw < 120 and $pcds < 90); #?? minaa what?
      
      # >oid type=CDS; aalen=257,80%,complete; clen=959; strand=+; offs=136-909;
      #? modify this to check %cds, minbigaa small if %cds is high
      #x next unless($aaw >= $MINBIGAA and $evc=~/(perfectdup|perfectfrag)/);
      
      my($mid)= $evc=~m/match:([\.\w]+)/;
      if($mrnaoids->{$mid}) { $bigdropcds{$oid}=$evc; $nbigdrop++; }
      } 
  } close($hin); }

  #? ($trsetnew,$ntrnew)= faextract($trset,$trsetnew,\%bigdropcds,'notid');
  ($ok,$hin)= openRead($trset); 
  if($ok) { 
    open(TO,'>',$trsetnew);  #a push @tmpfiles, $trsetnew;
    $ok=1; while(<$hin>) { 
    if(/^>(\S+)/){ $ok= $bigdropcds{$1}?0:1; 
      if($ok){ $ntrnew++; } else { $ntrcut++; }
      } 
    print TO $_ if($ok);
    } close($hin); close(TO); 
  }

  return($ntrnew, $trsetnew, $ntrnew+$ntrcut);
}


=item remove_mrna_aligned algo

    ($ntrnew,$trsetnew,$trsetagree)= remove_mrna_aligned($trname,$trset,$mrna);

  fastanrdb $trsetnew > $trsetnew.nr
  mv $trsetnew $trsetnew.old; mv $trsetnew.nr $trsetnew;
  $evigene/scripts/prot/make_consensus.pl $trsetnew
  # consensus transcripts n=33652 tabled in pig3c.notokmrna.nrtr.consensus
  trsetagree=$trsetnew.consensus

  makeblastdb -dbtype nucl -in $mrna -logfile /dev/null
  
  # bloptbada="-perc_identity 98 -qcov_hsp_perc 90 -culling_limit 1"
  # bloptbadb="-perc_identity 98 -qcov_hsp_perc 90  -ungapped -xdrop_ungap 4 -dust no"
  # bloptbadc="-perc_identity 98 -qcov_hsp_perc 90 "

  ** UPD/RETEST this, missing very large aligns from drops/perfdup,perffrag .. change to std blastn, btall?
  bloptnu="-perc_identity 98 -qcov_hsp_perc 90 -dust no "
  
  blastn $bloptnu -num_threads $ncpu -db $mrna -query $trsetnew -outfmt 6 -out mrnaperfaln.blastn

  perl -ne 'if(/^>(\S+)/){ $ok=not $nok{$1}; $infa=1; print if $ok; } elsif($infa){print if $ok; } elsif(/^\w/) { ($d)=split; $nok{$d}=1; }' \
    mrnaperfaln.blastn $trsetnew > $trname.notokmrna.tr

=cut

sub remove_mrna_aligned { # STEP1b_ now 1c
  my($trname,$trset,$mrna)= @_; # trset=="$trname.nomrna.tr";
  
  my $trsetnew="$trname.notokmrna.tr";  push @tmpfiles,$trsetnew;
  my $trsetagree="$trset.consensus";
  my $ntrnew= 0;
  my($cmd,$err,$ok,$hin);
  my $pident= 98; # -pHetero fixme
  
if(USE_recover_trset) {
  unless($UPDATEALL) {
    # fixme: where is trsetagree
    ($ntrnew,$trsetnew)= recover_trset($trname,$trset,$trsetnew); # handles .trids to .tr recovery
    if($ntrnew > 0) { return($ntrnew, $trsetnew,$trsetagree); }
    }
} else {  
  if( -s $trsetnew and not $UPDATEALL) {
    return(facount($trsetnew),$trsetnew,$trsetagree);
  }
}
  
  $cmd="$APPfastanrdb -i -f $trset > $trset.nr";
  $err= runcmd($cmd);
  unless($err){ rename("$trset.nr",$trset); } # orig trset is tmpfile

  $cmd="$EVIGENES/prot/make_consensus.pl $trset";
  $err= runcmd($cmd);
  $trsetagree="$trset.consensus"; push @tmpfiles, $trsetagree;
  # trap cmd.warn: consensus transcripts n=13046 tabled in plJBND.nomrna.tr.consensus

  # mrna.gz patch: gunzip -c $mrnagz | makeblastdb -title $mrna -output $mrna
  if($mrna =~ m/.gz$/) {
    my $mrnagz= $mrna; $mrna =~ s/.gz//;
    $cmd="gunzip -c $mrnagz | $APPmakeblastdb -dbtype nucl -title $mrna -out $mrna -logfile /dev/null"; # .gz NOT ok
  } else {
    $cmd="$APPmakeblastdb -dbtype nucl -in $mrna -logfile /dev/null"; # .gz NOT ok
  }
  push @erasefiles, map{ "$mrna.$_" } qw(nsq nin nhr);
  $err= runcmd($cmd); return if($err);

  # ** UPD/RETEST this, missing very large aligns from drops/perfdup,perffrag .. change to std blastn, btall?
  #above: my $usebself= $ENV{TESTmrnablast}||0; # ? TESTmrnablast UPD2222 default
  my $pctalign_trismrna= $PAL_NCRNAMIN;
  my $ralign_trismrna= $pctalign_trismrna/100;
  my $bloptself="-perc_identity $pident -evalue 1e-19 -dust no "; 

  my $bloptnu="-perc_identity $pident  -evalue 1e-19 -qcov_hsp_perc $pctalign_trismrna -dust no "; #? check opts again, maybe no -qcov, -dust yes
  my $blopt= ($usebself)? $bloptself : $bloptnu;
  
  my $blout="$trname.mrnaperf.blastn"; push @tmpfiles, $blout;

  $cmd= "$APPblastn $blopt -num_threads $NCPU -db $mrna -query $trset -outfmt 6 -out $blout";
  $err= runcmd($cmd); return if($err);

  my %ismrna; 
  if($usebself) { # add align hsps ? need tr sizes? ie 80-90% of trsize as per -qcov_hsp_perc
    my $trsizes= fasizes_nogap($trset,"","",1); # only.nt
    my($sal,$ltd,$lmd,$tdsize)=(0) x 9;
    open(B,$blout); 
    while(<B>){ if(/^\w/) { 
      my($td,$md,$pi,$al)=split; 
      if($td ne $ltd or $md ne $lmd){ $sal=0; } 
      if($td ne $ltd) { $tdsize= $trsizes->{$td}||$al; }
      $sal += $al;
      $ismrna{$td}=1 if($sal >= $ralign_trismrna * $tdsize); 
      $ltd=$td; $lmd=$md; 
      } } close(B);
    
  } else { # each align ismrna
    open(B,$blout); while(<B>){ if(/^\w/) { my($d)=split; $ismrna{$d}=1; } } close(B);
  }
  
  #? ($trsetnew,$ntrnew)= faextract($trset,$trsetnew,\%ismrna,'notid');
  ($ok,$hin)= openRead($trset); 
  if($ok) { 
    open(TO,'>',$trsetnew);  # ab: push @tmpfiles, $trsetnew;
    $ok=0; while(<$hin>) { 
    if(/^>(\S+)/){ $ok= $ismrna{$1}?0:1; 
      if($ok and my($ods)= m/oid=([^;,\s]+)/) {
        $ok=0 if(grep { $ismrna{$_} } split(",",$ods)); 
        }
      $ntrnew++ if $ok; 
      } 
    print TO $_ if($ok);
    } close($hin); close(TO); 
  }
  
  return($ntrnew,$trsetnew,$trsetagree);
}


=item ncrna_selfalign algo

  ($ntrnew, $blasttab)= ncrna_selfalign($trname,$trset, $use_altparopt)

  makeblastdb -dbtype nucl -in $trsetlong -logfile /dev/null
  
  ## FIXME odd problem due maybe to qcov_hsp_perc: no hits found for self.blast .. always should have self.tr
  ## switch to evg tr2cds self-blast opts?
  ## selfblopt="-task megablast -ungapped -xdrop_ungap 4 -dust no -perc_identity $CDSBLAST_IDENT -evalue $CDSBLAST_EVALUE "
  ## pident=95, eval=1e-5
  
  # test blastn opts, qcov_hsp_perc? 80 .. 60 ?  pctident? 98 97 96 ?
  # NOTE these opts block any seq with gaps, ie no hits even to self .. solution? dont use -qcov_hsp, use -dust ?
  bloptnu="-perc_identity 98 -qcov_hsp_perc 80 -dust no "
  
  echo blastn $bloptnu -num_threads $ncpu -db $trsetlong -query $trsetlong -outfmt 7 -out $trname.self98.blastn

=cut

sub ncrna_blastopt {
  my($trname, $use_altparopt, $ungapped, $pminalign)= @_; 

  # altparclassify.pl uses this, use same here?
  my $NCBLAST_IDENT=  $ENV{pctident} || 97; # option
  my $NCBLAST_EVALUE= $ENV{evalue} || 1e-9; # option
  
  my $bloptself="-task megablast -dust no -perc_identity $NCBLAST_IDENT -evalue $NCBLAST_EVALUE "; 
  $bloptself .=" -ungapped -xdrop_ungap 4 " if($ungapped);
          #^^ this way assumes will use btall sum of gapped hsps
  
  $pminalign ||= $PAL_NCRNAMIN;
  my $bloptnu="-perc_identity $NCBLAST_IDENT -evalue $NCBLAST_EVALUE -qcov_hsp_perc $pminalign -dust no "; #??
  
  return ($use_altparopt) ? $bloptself : $bloptnu;  
}

sub ncrna_selfalign {
  my($trname,$trset, $use_altparopt)= @_; 

  # altparclassify.pl uses this, use same here?
  my $NCBLAST_IDENT=  $ENV{pctident} || 97; # option
  my $NCBLAST_EVALUE= $ENV{evalue} || 1e-9; # option
  my $selfblopt="-task megablast -ungapped -xdrop_ungap 4 -dust no -perc_identity $NCBLAST_IDENT -evalue $NCBLAST_EVALUE "; 
  my $nuselfblopt="-perc_identity $NCBLAST_IDENT -evalue $NCBLAST_EVALUE -dust no "; # TEST? does it improve agreements

  # my $pminalign= $PAL_NCRNAMIN;
  my $bloptnu="-perc_identity $NCBLAST_IDENT -qcov_hsp_perc $PAL_NCRNAMIN -dust no "; #??
  
  my $blopt= ($use_altparopt) ? $selfblopt : $bloptnu; #??

  # use constant UNGAPPED => 1;
  # my $blopt= ncrna_blastopt($trname, $use_altparopt, UNGAPPED);
  
  # my $trsetnew="$trname.self98"; # 98 => $CDSBLAST_IDENT
  (my $trsetnew= $trset) =~ s/\.\w+$//; $trsetnew .= ".self$NCBLAST_IDENT";
  my $blout="$trsetnew.blastn";  
  my $ntrnew= 0;
  my($cmd,$err,$ok,$hin);
  if(not $UPDATEALL and ! -f $blout and -f "ncrnaset/$blout") { $blout= "ncrnaset/$blout"; }
  if( -s $blout and not $UPDATEALL) {
    ($ntrnew) = `grep -c '# Query' $blout`; chomp($ntrnew); #? need this grep?
    return($ntrnew, $blout);
  }
  
  $cmd="$APPmakeblastdb -dbtype nucl -in $trset -logfile /dev/null";
  $err= runcmd($cmd); return if($err);

  $blout="$trsetnew.blastn"; push @tmpfiles, $blout;
  $cmd= "$APPblastn $blopt -num_threads $NCPU -db $trset -query $trset -outfmt 7 -out $blout";
  $err= runcmd($cmd); return if($err);

  push @erasefiles, map{ "$trset.$_" } qw(nsq nin nhr);
  # for my $suf (qw(nsq nin nhr)) { unlink("trset.$suf"); } # blast idx  
  ($ntrnew) = `grep -c '# Query' $blout`; chomp($ntrnew);
  return ($ntrnew, $blout); #? or parse here
}

=item ncrna_replicates algo

  ($ntrnew, $trsetnew)= ncrna_replicates($trname,$trselfblast,$trsetagree,$trset);

  .. dont need/want  $trname.self98.blagree tmpfile here ?
  .. only need $trname.self98.blagree.trids, maybe as .consensus table as per make_consensus
  
  perl -ne '...' $trsetagree $trname.self98.blastn > $trname.self98.blagree

  grep 'replace=' $trname.self98.blagree | cut -f12 | sed 's/replace=//;' | sort -u > $trname.self98.blagree.dropids
  perl -ne '($d)=@v=split; if(@v==1){ $nok{$d}=1 } else { print "$d\n" unless($nok{$d}); }' \
    $trname.self98.blagree.dropids  $trname.self98.blagree | sort -u > $trname.self98.blagree.trids
  
  perl -ne 'if(/^>(\S+)/){ $ok=$ok{$1}; $infa=1; print if $ok; } elsif($infa){print if $ok; } elsif(/^\w/) { ($d)=split; $ok{$d}=1; }' \
    $trname.self98.blagree.trids $trsetlong > $trname.uniqrna.tr

=cut

sub ncrna_replicates {
  my($trname,$trselfblast,$trsetagree,$trset,$trsizes)= @_; 
  # my($ntrnew5, $truniqrna)  = ncrna_replicates($trname,$trselfblast,$trsetagree,$trlongnomrna,$trlongsizes);
  ## FAILING HERE blagree == 0 .. needs test .. typo if( -f bareword ) not $trselfblast

  my $check_align = ($trsizes and ref $trsizes) ? 1 : 0;#? option

  use constant RepFromBtall => 1;  # probably want this method
  
  (my $btall=$trselfblast) =~ s/\.\w+$/.btall/; # make/use this table for blagree ..
  (my $blagree=$trselfblast) =~ s/\.\w+$/.blagree/;
  my $trsetnew="$trname.uniqrna.tr";
  my $ntrnew= 0;
  my($cmd,$err,$ok,$hin);
  if( -s $trsetnew and not $UPDATEALL) {
    return(facount($trsetnew), $trsetnew);
  }

  our($agreeout,%agree,%blagree,%blreplace,%rsn,@rds,%didid,%did);
  open($agreeout,'>',$blagree); push @tmpfiles,$blagree;
  
  # write to $trname.self98.blagree ? as .consensus table as per make_consensus
  ## BUG trids missing Loc .. consensusof() bug
  ## FIXME block replace= for higher value nr agree{} set
  ## eg: arab6roo3dn4msoapk103loc4405t2 4 rep (soap,velv) ..
  ## replaced at blastn align by arab6roo3dntatrinLocDN33581c0g1t6 no reps.. but scored for 4 of prior.

  sub putagree { 
    my($ltd,$topal)=@_; 
    our($agreeout,%agree,%blagree,%blreplace,%rsn,@rds,%didid,%did);
    my @rsn=sort keys %rsn;  my $nput=0;
    my $nagree= @rsn + $agree{$ltd};
    if($nagree > 1) {
      my $pok=1; my $repok=0; # nu
      if(my $oal= $did{$ltd}) {
        my $od= $didid{$ltd};
        my $agr= $agree{$ltd}||0; my $ogr= $agree{$od} || 0;
        my $oblg= $blagree{$od}||0;
        if($agr == 0 and $ogr == 0) { $agr=$nagree; $ogr=$oblg; }
        if( $topal <= $oal and $agr <= $ogr ) { $pok= -1; }
        elsif($topal > $oal and $agr >= $ogr ) { $pok= 2; $repok=1 unless($oblg>$nagree); }
        elsif($topal > $oal and $agr < $ogr ) { $pok= ($agr>1) ? 2 : 0; } 
      }
      if($pok > 0) {      
        my $rsv= join",",map{ $_.":".$rsn{$_}; } @rsn; 
        $rsv = "equals" unless(@rsn);
        my $repval=0; 
        if($repok and my $od=$didid{$ltd}) { $repval="replace=$od"; $blreplace{$od}=$ltd; }
        $blagree{$ltd}=$nagree; # local dont bother w/ agreeout
        print $agreeout join("\t",$ltd,$nagree,$rsv,"aln=$topal",$repval)."\n"; $nput++;
        $didid{$ltd}=$ltd; $did{$ltd}=$topal; 
        for my $d (@rds) { if($did{$d} < $topal) { $didid{$d}=$ltd; $did{$d}=$topal; }  }
      }
    } 
    return($nput);
  }
  #------------------------------
  
  #   sub putagreeOLD { 
  #     my($ltd,$topal)=@_; 
  #     our($agreeout,%agree,%blagree,%blreplace,%rsn,@rds,%didid,%did);
  #     my @rsn=sort keys %rsn;  my $nput=0;
  #     my $nagree= @rsn + $agree{$ltd};
  #     if($nagree > 1) {
  #       my $pok=1;  # nu
  #       if(my $oal= $did{$ltd}) {
  #         my $od= $didid{$ltd};
  #         my $agr= $agree{$ltd}||0; my $ogr= $agree{$od} || 0;
  #         if( $topal <= $oal and $agr <= $ogr ) { $pok= -1; }
  #         elsif($topal > $oal and $agr >= $ogr ) { $pok= 2; }
  #         elsif($topal > $oal and $agr < $ogr ) { $pok= ($agr>1) ? 2 : 0; } 
  #       }
  #       if($pok > 0) {
  #         my $rsv= join",",map{ $_.":".$rsn{$_}; } @rsn; 
  #         my $repval=0; if(my $od=$didid{$ltd}) { $repval="replace=$od"; $blreplace{$od}=$ltd; }
  #         $blagree{$ltd}=$nagree; # local dont bother w/ agreeout
  #         print $agreeout join("\t",$ltd,$nagree,$rsv,"aln=$topal",$repval)."\n"; $nput++;
  #         $didid{$ltd}=$ltd; $did{$ltd}=$topal; 
  #         for my $d (@rds) { if($did{$d} < $topal) { $didid{$d}=$ltd; $did{$d}=$topal; }  }
  #       }
  #       # if(0) {   #... old
  #       # if(not ($did{$ltd} and $topal <= $did{$ltd}) ) { 
  #       #   my $rsv= join",",map{ $_.":".$rsn{$_}; } @rsn; 
  #       #   my $repval=0; if(my $od=$didid{$ltd}) { $repval="replace=$od"; $blreplace{$od}=$ltd; }
  #       #   $blagree{$ltd}=$nagree; # local dont bother w/ agreeout
  #       #   print $agreeout join("\t",$ltd,$nagree,$rsv,"aln=$topal",$repval)."\n"; $nput++;
  #       # } 
  #       # map{ $didid{$_}=$ltd; $did{$_}=$topal; } ($ltd,@rds); 
  #       # }
  #     } 
  #     return($nput);
  #   }
  #------------------------------
  
  my ($nag)= (0);
  if(-f $trsetagree) {
    ($ok,$hin)= openRead($trsetagree);
    while(<$hin>) { next if(/^\W/);  my($td,$nagr)=split;  $agree{$td}=$nagr; $nag++; }
    close($hin);
    loggit(0,"nagree=$nag in $trsetagree");
  }
  
  ## Maybe replace this w/ makeblastscore btall table parse .. simplifies decisioh,
  ## use sametr ==  min( align/tlen, align/rlen) > 0.90 and min(rlen/tlen, tlen/rlen) > 0.90
 
  if(RepFromBtall and -f $trselfblast) {
    #  $evigene/scripts/makeblastscore3.pl -each -needlen -pminlow 0.5 -sizes bevgc_evg3weed1rd.longnok.tr.qual bevgc_evg3weed1rd.longnok.self97.blastn > bevgc_evg3weed1rd.longnok.self97.btall
    my $trsizef="$trset.qual"; #?? CHECK
    unless(-f $trsizef) {
      if($trsizes and ref $trsizes) {
        open(F,'>', $trsizef); 
        for my $id (sort keys %$trsizes){ print F $id,"\t",$trsizes->{$id},"\n"; } close(F);
      } 
    }
    
    my $cmd= "$EVIGENES/makeblastscore3.pl -each -needlen -pminlow 0.5 -sizes $trsizef $trselfblast > $btall";
    my $err= runcmd($cmd);
    push @tmpfiles, $btall;
    
    my $MINPAL=$ENV{blagree_pal} || 0.80; # smaller? larger?
    my $MINPRW=$ENV{blagree_prw} || 0.70; #? smaller

    my $nokmer= ($ENV{'conkmer'})?0:1; # make_consensus.pl:consensusof() ($aconsensus =~ m/kmer/) ..
    my($ltd,$lrd,$topal,$sumal)= (0) x 9;
    ($ok,$hin)= openRead($btall);
    while(<$hin>) { 
      next if(/^\W/);  
      my($td,$rd,$bs,$ida,$aln,$tw,$rw)=split; # blast btall table OR agree table
      next if($td eq $rd or $tw<1 or $rw<1);
      my $pal= ($tw>$rw)? $aln/$tw : $aln/$rw;
      my $prw= ($tw>$rw)? $rw/$tw  : $tw/$rw;
      ## cancel this skip if $agree{$td} .. but dont add rd rsn
      my $ok = ($pal < $MINPAL or $prw < $MINPRW) ? 0 : 1;
    
      if($td ne $ltd) { $ntrnew += putagree($ltd,$topal) if($ltd); %rsn=@rds=(); } 
      $ltd=$td; $lrd=$rd;
      if($ok) {
        $topal= $aln if($aln > $topal); # $tpal=$pal; $tprw=$prw; 
        foreach my $sd ($td,$rd) {  
          my $st= $sd;
          if($st=~s/[Ll]oc.*$//) { $st =~ s/k\d+$//; } 
          else { $st=substr($sd,0,9); }
          $rsn{$st}++;
        }
        push @rds,$rd;
      }        
    
    } close($hin);
    if(1) { $ntrnew += putagree($ltd,$topal) if($ltd); %rsn=@rds=(); }        
  } # RepFromBtall
 
  if(not RepFromBtall and -f $trselfblast) { # if not, putagree(%agree) ??
    ## this should be from trasm-consensus parser:
    ##  my($s)= (m/^(\w+)k\d+[Ll]oc.*$/)? $1 : substr($_,0,9);
    my $nokmer= ($ENV{'conkmer'})?0:1; # make_consensus.pl:consensusof() ($aconsensus =~ m/kmer/) ..
    
    my($ltd,$lrd,$topal,$sumal)= (0) x 9;
    ($ok,$hin)= openRead($trselfblast);
    while(<$hin>) { 
      next if(/^\W/);  
      my($td,$rd,$pi,$al)=split; next if($td eq $rd);
 
      if($td ne $ltd) { $ntrnew += putagree($ltd,$topal) if($ltd); %rsn=@rds=(); $sumal=$topal=$lrd=0; } 
      $sumal=0 if($rd ne $lrd);
      $sumal += $al; # multiple hsp per td x rd
      
      # need align size / tr size check, unless using blastn -qcov_hsp_pct to enforce long align
      # if this blastn has align split to hsps, need to add hsps .. as with makeblastscore.pl
      my $ok=1;
      if($check_align and not $agree{$td}) { 
        my $tw= $trsizes->{$td}||0;  # bug bad trsizes ?
        my $pal= ($tw>0) ? $sumal/$tw : 0; # use sum align here for multiple hsp per td x rd
        $ok= ($pal*100 >= $PAL_NCRNAMIN); 
        # $ok=1 if($agree{$td});
        # next unless($ok); 
        }
      
      $ltd=$td; $lrd=$rd;
      if($ok) {
        $topal= $sumal if($sumal > $topal);
        # see make_consensus.pl:consensusof()
        foreach my $sd ($td,$rd) {  #?? this changes td, rd ?? shouldnt
          my $st= $sd;
          if($st=~s/[Ll]oc.*$//) { $st =~ s/k\d+$// if($nokmer or $st=~/idba/); } 
          else { $st=substr($sd,0,9); }
          $rsn{$st}++;
        }
        #o my($ts,$rs)= map { my($s)= (m/^(\w+)k\d+[Ll]oc.*$/)? $1 : substr($_,0,9); $rsn{$s}++; $s; } ($td,$rd);  
        push @rds,$rd; #  unless($rd eq $lrd); ## AND ok lrd ..
      }  
        
      #? $ltd=$td; $lrd=$rd; # mangled from above s change??
      } close($hin);
    $ntrnew += putagree($ltd,$topal) if($ltd);
  }  # not RepFromBtall
  
  close($agreeout);
  
  # write list to $trname.self98.blagree.trids ? or just refer to agreeout table
  my @truniqids= sort grep{ not $blreplace{$_} }  keys %blagree;
  $trsetnew="$trname.uniqrna.tr"; push @tmpfiles, $trsetnew;
  ($trsetnew,$ntrnew)= faextract($trset,$trsetnew,\@truniqids) if(@truniqids);
  # NOTE faextract() returns (trset) if no ids in trlongs
  
  return ($ntrnew, $trsetnew);  
}



=item ncrna_classgene / classifylocus

  change output table to match pubids? for merging w/ that:
    pubtrid oid pubgene alti [drop|cull]altclass aaqual pidn/aln notes/tscore/tscorevec
  
=cut

my($didclhdr);

sub ncrna_classgene { # from ncrna_altparscore.pl
  my($outh, $nout, $dgeneid, $dgenesar, $hasduphash, $hasokcds)=@_; #drop: $havedups,
  my($ng,$maxlen,$nok,$ndrop)=(0) x 9; 
  my(%tsc,%tsv,%tscvec,%tdid);
  return($nok,$ndrop) unless(@$dgenesar);

  ## NOW allevd.hdr: IDtrasm	trlen	aalen	pcds	cdpot	agree	dglocus	dgclass
  ## Maybe use @$dgenesar == hash of qualkey=>qualval ? eg fir $TEST_OKCDS:  add  cdsaln
  ##? add col for hasokcds{td} == align to okay.cds
  
  my (@hdr,%ihdr); # IDtrasm >> originalID better
  # pubids: #Public_mRNA_ID	originalID	PublicGeneID	AltNum	Class	AAqual	pIdAln	Notes	Oids

  my @DEFhdr=qw(IDtrasm trlen aalen pcds cdpot	agree	dglocus dgclass); 
  push @DEFhdr, 'cdsaln' if($hasokcds); #$TEST_OKCDS: add hasokcds col here 
  @hdr= @DEFhdr;
  #x if($colhdr and ref($colhdr) =~ /ARRAY/){ @hdr= @$colhdr; } else {  }
  #x for my $i (0..$#hdr){ $ihdr{ $hdr[$i] }= 1+$i; }
  
  #X these are bad .. change to @$dgenesar == hash of qualkey=>qualval
  #X sub valof{ my($tsv,$k)=@_; if(my $i= $ihdr{$k}){ return $$tsv[$i-1] || 0; } else { return 0; } } 
  #X sub setvalof{ my($tsv,$k,$val)=@_; if(my $i= $ihdr{$k}){  $$tsv[$i-1]=$val; } } 
  # my($td,$trlen,$aalen,$pcds,$cdpot,$agree,$dgid,$dgcla,$cdsaln)= map{ valof($tsv,$_) } @DEFhdr;
  
  #  WT may be option? report in out table?
  # my %WT4d=( trlen=> 0.50, pcds => 0.25, cdpot => 0.25, agree => 0.25 ); 
  my %WT5e=( trlen=> 0.70, pcds => 0.15, cdpot => 0.15, agree => 0.25 ); 
  my %WT6f=( trlen=> 0.70, pcds => 0.15, cdpot => 0.15, agree => 0.25, notcds=> 0.25 ); 
  my %WT7g=( trlen=> 0.70, pcds => 0.15, cdpot => 0.15, agree => 0.25, cdsaln=> -0.25 ); 
  my %WT= %WT7g;
  
  # my @dgenes= @$dgenesar;
  for my $tsv (@$dgenesar) {
    my($td,$trlen,$aalen,$pcds,$cdpot,$agree,$dgid,$dgcla)= @$tsv;
    # my($trlen)= valof($tsv,'trlen'); # bad
    $maxlen=$trlen if($trlen>$maxlen); $ng++;
  }
  return(0) unless($ng>0 and $maxlen>0); # Urgh: maxlen==0 
  
  # calc trscore for ranking best > worst
  for my $tsv (@$dgenesar) {
    my($td,$trlen,$aalen,$pcds,$cdpot,$agree,$dgid,$dgcla,$cdsaln)= @$tsv;

    # scale vars to approx 1..100  range, 100 is best score
    my $trlensc= int(0.5 + 100*$trlen/$maxlen);  # 10..100
    my $putr   = 100 - $pcds;  # 10 .. 100
    my $ncpot  = int(100*(1 - $cdpot));  # -50 .. 70 ?
    $agree     = 10*$agree; # 20..60 ?
    
    #?? WT{trlen} may be too low vs others .. WT pcds/cdpot too high?
    ##  saves some rather short alts (10% of longtr) via high score from pcds/cdpot
    my $tscorevec= join",",$trlensc,$putr,$ncpot,$agree;
    my $tscore = $trlensc * $WT{trlen} + $putr * $WT{pcds} + $ncpot * $WT{cdpot} + $agree * $WT{agree};
    if($hasokcds) { 
      $cdsaln ||= 0; 
      my($okid)= ($cdsaln=~s/,(.*)$//) ? $1 : 0;  
      $cdsaln=100 if($cdsaln>100); 
      
      #? if cdsaln > 70..80, this should be dropped as redundant w/ okcds , use -WT or flag drop?
      # OR these are valid mRNA/CDS with long UTR that now bumps them into ncrna class .. relocus by good CDS align?
      # preserve cdsaln,okid
      
      # negative WT.cdsaln now: -0.25
      $tscore +=  $cdsaln * $WT{cdsaln}; $tscorevec.=",$cdsaln";  
      #o my $notcds = int(100 - $cdsaln); # maybe should be 100/0 value .. any scored cdsaln is bad/dupl qual?
      #o $tscore +=  $notcds * $WT{notcds}; $tscorevec.=",$notcds";  
      }
    
    $tsc{$td} = $tscore; $tscvec{$td} = $tscorevec; $tsv{$td}= $tsv; 
  }
  
  my @tdo= sort{ $tsc{$b} <=> $tsc{$a} } keys %tsc; 
  my $nt=@tdo;
  my $td1= $tdo[0];
  my $tsv1= $tsv{$td1}; 
  my $tscvec1= $tscvec{$td1}; 
  my $tsc1 = $tsc{$td1}; 
  my($tlen1)= $tscvec1=~m/^(\d+)/;
  
  #     my($td,$trlen,$aalen,$pcds,$cdpot,$tagree1,$tgid,$dgcla,$cdsaln)= @$tsv;
  my($tagree1)= $$tsv1[5];
  my $tgid= $$tsv1[6]; # must be same for all @tdo
  # my($tagree1)= valof($tsv1,'agree');
  # my($tgid)= valof($tsv1,'dglocus');
  my($mlong,$mpar,$magree)=(0,0,0);
  my $idpre= $IDPREncrna || "EVn"; # EVr EVm EVn .. which?
  
  for my $i (0..$#tdo) {
    my $td=$tdo[$i]; my $tsv=$tsv{$td}; 
    my $tscore= sprintf "%.1f", $tsc{$td}; 
    my $tscorevec= $tscvec{$td};
    
    my($tdX,$trlen,$aalen,$pcds,$cdpot,$tagree,$dgid,$tcla,$cdsaln,$aaqualX)= @$tsv;
    # my $tcla= $$tsv[7]; # DANG was 6 == dgeneID
    #  my $tagree= $$tsv[5];
    # my($tcla)= valof($tsv,'dgclass');
    # my($tagree)= valof($tsv,'agree');
    my($tlen)= $tscorevec=~m/^(\d+)/;
    
    # classify, add drop_redundant using table of nearperfect matches from self.btall
    # ie. perfdupl and perffrag are in drop set as per tr2aacds, 1st of id pair in this tscore order wins
    # maybe add okay/drop column as per trclass, or use 'drop|cull' prefix on class as per pubids 
    
    # mainuni == uni == noclass, reuse noclass?
    # mainpar == new subclass, top-scored ~paralog
    # mainlong == keep-alt as longest, revert to 'alt' class?
    
    my($keep,$drop)= (0,0); # count to decide class prefix
    if($i==0){ $keep++; } 
    my $cla= ($nt==1) ? "mainuni" : ($i==0) ? "main" : "alt";
    if($cla eq "alt" and $tlen == 100 and $tlen1 < 98) { unless($mlong++){ $cla="altlong"; $keep++; } }
    # altag/tagree save may be unwanted .. 
    if($cla eq "alt" and $tagree > 1 and $tagree1 == 0) { unless($magree++){ $cla="altag"; $keep++; } }
    if($cla eq "alt" and $tcla eq "par") { $cla= "altpar";  unless($mpar++) { $keep++; } }
    if($tscore < $MINSCORE) { $drop++; $keep=0; } 
    elsif($cla =~ /^alt/ and $tscore < $MINALT * $tsc1) { $drop++; }

    my($isdup,$dupof)=(0,0);
    if(exists $hasduphash->{$td}) { # $havedups and ..
      for my $d (keys %tdid){ if($isdup= $hasduphash->{$td}{$d}){ $dupof=$d; last; } }
      if($isdup) { $drop++; $cla .= $isdup; } # add value of hasdup to class: ($isdup=~/frag/)?"frag":"dup"; 
    }
    
    # my $okdrop=""; # class prefix: okay = "", drop = "drop", cull = "cull" 
    my($okdrop,$okval)=("",0);
    if($drop > 0 and $keep==0) { $okdrop="drop"; $okval=-1; $ndrop++; }
    else { $okdrop=""; $okval=1; $nok++; }
    $cla=  $okdrop . $cla;
    #o $$tsv[7]= $cla; # replace $tcla # DANG was 6 wrong col.
    # setvalof($tsv,'dgclass',$okdrop . $cla);
    
    my $ialt= 1+$i;
    my $geneid= $idpre . sprintf "%06d", $tgid;
    my $pubid = $geneid . "t$ialt";
    
    unless($didclhdr++) { 
      #chng to pubids? #Public_mRNA_ID	originalID	PublicGeneID	AltNum	Class	AAqual	pIdAln	Notes	Oids
      #O print $outh join("\t","PublicID",@hdr,"trscore","scorevec")."\n";
      #o my @addc= qw(trlen aalen pcds cdpot agree cdsaln); # change this..
      #o my @addc= qw(trlen aalen); # change this..
      #o print $outh join("\t",@pubc,"trscore","scorevec",@addc)."\n";
      
      my @pubc= qw(PublicID	originalID	PublicGeneID	AltNum	Class AAqual TRqual Notes);
      print $outh "#".join("\t",@pubc)."\n";
      
      my $wts= join", ", grep /\w/, map{ ($WT{$_}) ? "$_=".$WT{$_} : "" } @hdr; 
      print $outh "#tscore_weights: $wts\n" if($wts); 
      }
      
    #o my @addc= ($trlen,$aalen,$pcds,$cdpot,$tagree,$cdsaln); # change this..
    #o my @addc= ($trlen,"$aalen,$pcds%"); # change this.. into Notes format = taga:vala,tagb:valb,..
    #o print $outh join("\t",$pubid,@$tsv,$tscore,$tscorevec)."\n"; 
    #o print $outh join("\t", $pubid, $td, $geneid, $ialt, $cla, $tscore,$tscorevec,@addc)."\n"; 
    my $aaqual="$aalen,$pcds%,na"; # add ,complete|partial 
    $aaqual=$aaqualX if($aaqualX and $aaqualX =~ m/,\d+%,/);# preserved val>
    
    my $trqual=$trlen; # or? "$trlen,cp$cdpot";
    my @notes=();
    push @notes, "aaref:$cdsaln" if($cdsaln =~ /,/);
    push @notes, "dupof:$dupof" if($dupof);
    push @notes, "agree:$tagree" if($tagree > 1);
    push @notes, "tscore:$tscore";
    push @notes, "scorevec:$tscorevec" if($debug); 
    my $notes= join",", map{ s=,=/=g; $_ } @notes;  
    
    print $outh join("\t", $pubid, $td, $geneid, $ialt, $cla, $aaqual, $trqual, $notes)."\n";
    
    $tdid{$td}++;#? record okay/drop here for dups? 
    # $tdid{$td}= $okval; # ( $okdrop eq "drop" ) ? -1 : ( $okdrop eq "cull" ) ? -2 : 1;
  }
  return($nok,$ndrop);
}


=item ncrna_classgenes

  v2: final step? ncrna_classgene(@dgenetrs), where @dgenetrs are all trset in gene group, w/ evidence
   classifier weights evidence for rank by trscore, keeps/drops alt/par of gene group depending on scores
   output table in publicset/pubids format? w/ class 'main/alt/noclass=uni' and 'drop|cull' qualifier
   
   ($ngeneok,$ntrok,$ntrdrop)= ncrna_classgenes($ngene,$genetab,$allevdtab) : read evdtab, each tr sorted by genetab dgene ids, with evidence
      my($nok,$ndrop)= ncrna_classgene($outputh, $dgeneid,\@dgenes,$hasduphash);
   
   update ($ntrnew5, $truniqrna) removing drops from classgenes

** FIXME: need to be sure of allevdtab columns, match here using header? row 1
human18ncx.longnok.allevd.tab
IDtrasm	trlen	aalen	pcds	cdpot	agree	dglocus	dgclass
humang010008t1	3070	103	10	0.464	0	000694	main
humang010008t3	2929	103	10	0.464	0	000694	alt
humang010008t4	3025	103	10	0.492	0	000694	alt
humang010008t5	2875	103	10	0.464	0	000694	alt

=cut

sub ncrna_classgenes {
  my($trname, $trset, $genetab, $allevdtab, $hasduphash)= @_;
  # genetab not used here, only allevdtab that should contain genetab data
  
  (my $outclass=$allevdtab) =~ s/.\w+$/.ncrna_class/; # change this .. pubids-like, class?
  my $trsetnew=$outclass; # ??
  my $ntrnew= 0; # == $nloc or $snok ?
  my($cmd,$err,$ok,$hin);
  my($nloc,$nok,$ndrop,$ldglocus,$snok,$sndrop)=(0) x 9;
  
  if(not $UPDATEALL and ! -f $outclass and -f "ncrnaset/$outclass") { $outclass= "ncrnaset/$outclass"; }
  if( -s $outclass and not $UPDATEALL) {
    ($ok,$hin)= openRead($outclass);    
    while(<$hin>){ next if(/^\W/); 
      my($pid,$oid,$gid,$ti,$cla)= split;
      if($cla =~ /^drop/){ $sndrop++; } else { $snok++; }
      $nloc++ if($gid ne $ldglocus);
      $ldglocus= $gid;
      } close($hin);
    return( $nloc, $snok, $sndrop, $outclass);  
  }

  #evd cols: (ID trlen aalen pcds cdpot agree dglocus dgclass);
  my(@hdr,@dgenes);
  $hasduphash = {} unless(ref($hasduphash));
  
  my($nout,$outh)=(0);
  open($outh,'>',$outclass); push @tmpfiles, $outclass;
  
  my $hdr=`head -1 $allevdtab`; chomp($hdr); @hdr= split" ",$hdr; # not via sort
  my $hascdsaln= grep(/cdsaln/, @hdr); #?? pass to classgene()
  
  $cmd="sort -k7,7 -k2,2nr -k1,1 $allevdtab |"; # fixme : dgene now always col 7 ?
  open($hin,$cmd) or loggit(LOG_DIE,"ncrna_classgenes: sort $allevdtab");
  while(<$hin>) { 
    next if(/^\W/);  
    my @v= split; 
    if(/^ID/){ next; } # @hdr=@v;  NOT 1st row due to sort; use to get val cols ;
    my($td,$trlen,$aalen,$pcds,$cdpot,$agree,$dglocus,$dgclass)=@v;
    #old: my($td,$trlen,$pcds,$cdpot,$agree,$dgid,$dgcla)= @tsv= @v[0,1,3,5,6,11,12];
  
    if($dglocus ne $ldglocus) { 
      if(@dgenes) {
      ($nok,$ndrop)= ncrna_classgene($outh, $nout, $ldglocus, \@dgenes, $hasduphash, $hascdsaln); # \@hdr
      $nout++; $nloc++ if($nok); $snok+=$nok; $sndrop+=$ndrop;
      }
      @dgenes=(); 
      }
    push @dgenes, [@v]; $ldglocus=$dglocus; 
  }
  close($hin);
  
  ($nok,$ndrop)= ncrna_classgene($outh, $nout, $ldglocus, \@dgenes, $hasduphash, $hascdsaln); # \@hdr
  $nout++; $nloc++ if($nok); $snok+=$nok; $sndrop+=$ndrop; 
  close($outh);
  
  return($nloc,$snok,$sndrop,$outclass);  
}

sub ncrna_pubseqset {
  my($trname,$trset,$pubidtab)= @_;  
  
  my $trsetnew="$trname.ncrna_pub.fa"; # publicset/ suffix
  # my $trsetnew="$trname.okay.ncrna";   # okayset/ suffix ? which
  (my $trids=$trsetnew) =~ s/.\w+$/.trids/; # cut of col2 of $pubidtab, removing drops. dont really need
  
  ## Also, extract .aa, .cds from inputset/trname.aa,cds.gz unless/if option 
  ## .. if aasize > MINAA, need inputset/*.{aa,cds}[.gz]
  my ($naddseq,%addseq)= (0);
  
  for my $seq (qw(aa cds)) {
    my $ins= filegz( "inputset/$trname.$seq"); 
    # $ins .= ".gz" if( ! -f $ins and -f "$ins.gz");
    
    if($DOMERGE) {  # FIXME: look where? for trname_mrg.nclong.aa,cds  see makeSubmergeTr(); also makes aa,cds
      (my $inset= $trset) =~ s/\.w+$/.$seq/;
      $ins= $inset if( -f $inset);      
    } 
    
    if(-f $ins) { 
      #?? what suffix? trname.ncaa_pub.fa trname.nccds_pub.fa; trname.ncrna_aa.fa trname.ncrna_pub.aa
      (my $outs= $trsetnew) =~ s/ncrna/nc$seq/; 
      $addseq{$seq}= "$ins\t$outs"; $naddseq++;
    }
  }
  
  my $ntrnew= 0; my $ntrcut=0;
  my($cmd,$err,$ok,$hin);
  if(not $UPDATEALL and ! -f $trsetnew and -f "ncrnaset/$trsetnew") { $trsetnew= "ncrnaset/$trsetnew"; }
  if( -s $trsetnew and not $UPDATEALL) {
    return(facount($trsetnew),$trsetnew);
  }

  my(%pubid,%puban,%drops);
  ($ok,$hin)= openRead($pubidtab);  unless($ok){ loggit(1,"pubidtab bad: $pubidtab"); return; }
  # qw(PublicID	originalID	PublicGeneID	AltNum	Class AAqual TRqual Notes);
  while(<$hin>) { 
    next if(/^\W/);
    my($pubid,$oid,$gid,$ti,$cla,$aaq,$trlen,$notes)=split"\t";
    chomp($notes); # notes= aaref: agree: dupid: tscore: .. other
    my($aw,$pcds,$ac)= split",",$aaq;
    
    ## FIXME: need some standard for nctype class, aaq: utrpoor/utrbad useful, package aaqual stat?
    my $nctype= ($aw > 120 and $pcds > 60)?"misc_RNA":"ncRNA"; # other quals to call type? codepot? aaref?
    
    $pubid{$oid}= $pubid; 
    $puban{$oid}= "type=$nctype; aalen=$aaq; clen=$trlen; evgclass=$cla; notes=$notes;";
    $drops{$oid}=$cla if($cla =~ /^drop/);
  } close($hin); 

  # write pubid{oid} = pubid to pubset.trids to match other subset.trids, dont really need but for balance
  open(SF,'>',$trids);  push @tmpfiles, $trids;  
  foreach my $oid (grep { not $drops{$_} } sort keys %pubid) { print SF "$oid\n"; } close(SF); 

  ($ok,$hin)= openRead($trset); 
  unless($ok){ loggit(1,"trset bad: $trset"); return; }
  open(TO,'>',$trsetnew); push @tmpfiles, $trsetnew; # ? @tmpfiles now == ncrnaset/
  $ok=0;
  while(<$hin>) { 
    if(/^>(\S+)/) { 
      my $oid=$1; my $pd= $pubid{$oid}; 
      $ok=($pd and not $drops{$oid})?1:0; 
      unless($ok){ $ntrcut++; next; }
      my $an= $puban{$oid}||"type=ncRNA;"; 
      s/^>.*$/>$pd $an oid=$oid/;
      $ntrnew++; 
      } 
    print TO $_ if($ok);
  } close($hin); close(TO); 
 
  # merge w/ above ncrna output; limit to aasize >= MINAA_NCRNA
  my $ncaaset="";
  for my $seq (sort keys %addseq) {
    my $stype= ($seq eq 'aa')?'protein':$seq;
    my($ins,$outs)= split"\t", $addseq{$seq};
    ($ok,$hin)= openRead($ins); unless($ok){ loggit(1,"trset bad: $ins");  next; }
    open(TO,'>',$outs); push @tmpfiles, $outs;  
    my $nok=0;
    while(<$hin>) { 
      if(/^>(\S+)/) { 
        my $oid=$1; my $pd= $pubid{$oid}; 
        $ok=($pd and not $drops{$oid})?1:0; 
        if($ok){ 
          my($aw)= m/aalen=(\d+)/?$1:0; 
          $ok=0 if($aw>0 and $aw < $MINAA_NCRNA);
        }
        next unless($ok); $nok++;
        my $an= $puban{$oid}||"type=ncRNA;"; 
        $an =~ s/type=\w+/type=$stype/;
        s/^>.*$/>$pd $an oid=$oid/;
        } 
      print TO $_ if($ok);
    } close($hin); close(TO); 
    $ncaaset= $outs if($seq eq 'aa' and $nok > 0);
  }
  
  return($ntrnew, $trsetnew, $ncaaset, $ntrnew+$ntrcut);
}


=item ncrna_aauniq_refblast
  
  this is for tr2ncrna pipe, equiv to sra2genes SCRIPT_ncaablastp
  
  
  this probably fits better into tr2ncrna pipe, than in sra2genes following tr2ncrna.
  
  blastp ncrnaset/ncaa x refaa, check for uniq ref hits not in okayset.aa x refaa (btall)
  eg.  aaset=ncrnaset/plEFMS.ncaa_pub.fa

  Note this uses temp  IDprefix 'EVn' for ncrnaset .. change to that?
gunzip -c ../aaeval/refgenes-plEFMS.okay.btall.gz | cat plEFMSncaa_refplant.btalln - | \
perl -ne 'next unless(/^Tortax/);
($td,$rd,$bs,$idn,$al,$tw,$rw)=@v=split; $sc=$bs;
($rg=$rd)=~s/\.\d+$//;  if($td=~/EVn/){ $osc=$topn{$rg}; if($sc>$osc) {
$topa{$rg}=[@v];  $topn{$rg}=$sc; } } else { if($osc=$topn{$rg}) {
$osp=$topp{$rg}||0; $pok=($sc >= $osc - 1 and $sc > $osp); $ppok=($sc
> $osp); if($pok){ $topp{$rg}=$sc; $toppa{$rg}=[@v]; } elsif($ppok){
$topp{$rg}=$sc; $toppb{$rg}=[@v]; } } } END{ for $rg (sort keys %topn)
{ $topa=$topna=$topa{$rg}; $ntop="ncaa"; $pb="napub";
if($toppa=$toppa{$rg}) { $ntop="pubaa"; $topa=$toppa; }
elsif($toppb=$toppb{$rg}) { $pb=join",", @{$toppb}[0,4,5,6,7]; } 
print join("\t",@$topa,$ntop,$pb)."\n"; } }' \
  > plEFMS_puborncaa_refplant.btall

=item ncaatest plEFMS

# these are only 8 cases (6 ncaa genes for 2 ref parlogs)  w/o better, or nearly same, ref hit in okayset/okay.aa 
grep  ncaa plEFMS_puborncaa_refplant.btallb | grep napub

TortaxEVn040514t1       AT1G01470.1     68.6    49      150     227     151     ncaa    napub
TortaxEVn040514t1       AT2G46140.1     67.8    58      157     227     166     ncaa    napub  *refpar
TortaxEVn053096t1       AT1G32630.1     97.8    42      78      171     623     ncaa    napub
TortaxEVn039970t2       AT1G72600.2     163     78      126     128     208     ncaa    napub
TortaxEVn040609t1       AT2G05310.1     112     49      71      128     125     ncaa    napub  *consv
TortaxEVn040609t1       AT4G13500.1     112     50      71      128     125     ncaa    napub  *refpar
TortaxEVn050066t2       AT4G29540.2     73.6    45      65      152     336     ncaa    napub  *consv
TortaxEVn044638t2       AT5G49800.1     175     99      175     173     242     ncaa    napub  *consv

# .. missed conserved genes
# trclass drops these due to gapbad, parthi flag
TortaxEVn044638t2/Tortaxtrsoap1a_sBn1l1ERR2040863soapk31Loc38490t1        drop    althi   Tortaxtrvelo1a_sBn1l1ERR2040863velvk35Loc21493t1        100/72/.        198,51%,complete-utrpoor-gapbad 0,0,pflag:16
TortaxEVn040609t1/Tortaxtrvelo1a_sBn1l1ERR2040863velvk35Loc10438t3        drop    parthi  Tortaxtridba1a_sBn1l1ERR2040863idbaidbtk27Loc2887       100/49/./Tortaxtridba1a_sBn1l1ERR2040863idbaidbtk27Loc2886      128,56%,complete        0,0,aacons:3,pflag:0
TortaxEVn050066t2/Tortaxtrvelo1a_sBn1l1ERR2040863velvk35Loc10768t2        drop    main    Tortaxtrsoap1a_sBn1l1ERR2040863soapk25Loc478t2  100/27/.        184,72%,partial3-gapbad 0,0,pflag:16
#......
  
=cut


sub ncaa_match_okaa_btall {  # match_ncaa_okaa_btall 
  my($outfile,$ncaabtall,$okaabtall)=@_;
  our(%topa,%topn,%topp,%toppa,%toppb);
  
  sub readbtall { 
    my($src, $btall)=@_;    
    our(%topa,%topn,%topp,%toppa,%toppb);
    my($nok)=(0);
    open(B,$btall) or return(0);
    while(<B>) { 
      next if(/^\W|^Query/);
      my @v=split;  my($td,$rd,$bs,$idn,$al,$tw,$rw)=@v;
      my $sc=$bs; # choice of score: bs,idn,al .. idn may be best but namegenes uses bs
      (my $rg=$rd)=~s/[t\.]\d+$//; # ref_ids may not have gene.transcript id patt
      if($src eq 'ncaa') {
        if($sc > $topn{$rg}) { $topa{$rg}=[@v];  $topn{$rg}=$sc; $nok++; } 
      } elsif(my $osc=$topn{$rg}) {
        my $osp= $topp{$rg}||0; 
        my $ppok=($sc > $osp); $nok++ if($ppok);
        my $pok= ($sc >= $osc - 1 and $ppok); 
        if($pok){ $topp{$rg}=$sc; $toppa{$rg}=[@v]; } 
        elsif($ppok){ $topp{$rg}=$sc; $toppb{$rg}=[@v]; }
      }     
    } close(B);
    return($nok);
  }
  
  my($nok)= readbtall('ncaa',$ncaabtall); # read 1st
  my($aok)= readbtall('okaa',$okaabtall); # only ncaa hits
  return(0) unless($nok); # aok ?
  
  open(O,'>',$outfile); 
  my ($mok,$ncuniq)=(0,0);
  for my $rg (sort keys %topn) { 
    my($topa,$topna,$ntop,$pb,$toppa,$toppb);
    $topa=$topna=$topa{$rg}; $ntop="ncaa"; $pb="napub";
    if($toppa=$toppa{$rg}) { $ntop="pubaa"; $topa=$toppa; }
    elsif($toppb=$toppb{$rg}) { $pb=join",", @{$toppb}[0,4,5,6,7]; } 
    # ^^ refine this to test topa.idaln - toppb.idaln, flag if ncaa >> pubaa
    print O join("\t",@$topa,$ntop,$pb)."\n"; 
    $mok++; $ncuniq++ if($ntop eq "ncaa" and $pb eq "napub");
    # ncaa uniq ref are: ncaa.napub
  } close(O);
  
  return($mok,$ncuniq);
}

sub has_okayset_refblast {
  my($trname)= @_; 

  my $P_REFAA = finddata('REFAA') || finddata('refset/ref*.aa');  
  return (0,'REFAA missing') unless($P_REFAA and -f $P_REFAA);
  
  my $refp= basename($P_REFAA); $refp =~ s/\.\w+$//; 
  my $okbtall="$refp*$trname*.btall"; # SCRIPT_evgblastp form
  my $P_OKAABLAST = finddata($okbtall) || finddata("aaeval/$okbtall");  
  return (0,'OKAABLAST data missing:'.$okbtall) unless($P_OKAABLAST and -f $P_OKAABLAST);
  return( 2, $P_REFAA, $P_OKAABLAST);
}

sub ncrna_aauniq_refblast {
  my($trname, $ncrna_aaset, $refaa, $okbtall)= @_; 

  my $blout="$trname.ref_ncaa.blastp";  
  (my $ncaabtall=$blout ) =~ s/\.\w+$/.btall/; 
  (my $outtab=$ncaabtall) =~ s/\.\w+$/_okaa.difftab/; # = trname.ref_ncaa_okaa.difftab ??
  (my $outuniq=$outtab) =~ s/\.\w+$/.uniqtab/;  # = trname.ref_ncaa_okaa.uniqtab ??
  
  my($cmd,$err,$ok,$hin);
  if(not $UPDATEALL and ! -f $outtab and -f "ncrnaset/$outtab") { $outtab= "ncrnaset/$outtab"; }
  if( -s $outtab and not $UPDATEALL) {
    my($ntrnew) = `wc -l $outtab`; chomp($ntrnew); #? need this grep?
    my($nuniq)= `grep -c 'ncaa.napub' $outtab`; chomp($nuniq);
    return($ntrnew, $outtab, $nuniq);
  }
  
  use constant STEPerr => -1;
  my $P_AASET = $ncrna_aaset; # "ncrnaset/$trname.ncaa_pub.fa"; # FIXME check
  if(-s $P_AASET) { loggit(0,"Query data=$P_AASET"); } 
  else { return (STEPerr,'ncrnaset/ncaa_pub.fa');  }

  my($aok,$P_REFAA,$P_OKAABLAST)= ($refaa and $okbtall) ? (2,$refaa,$okbtall) : (0);  
  ($aok,$P_REFAA,$P_OKAABLAST)= has_okayset_refblast($trname) unless($aok); # or params
  loggit(0,"REFAA data=",$P_REFAA,"OKAABLAST data=",$P_OKAABLAST);
  return (STEPerr,'Missing REFAA blast data') unless($aok);
 
  my($P_BLASTP,$P_NCBIBIN)= findapp('blastp', 1); 
  # my($P_MAKEBLAST)= findapp('makeblastdb', 1); 
  
  unless( -f "$P_REFAA.psq") {
  $cmd="$APPmakeblastdb -dbtype prot -in $P_REFAA -logfile /dev/null"; # P_MAKEBLAST
  $err= runcmd($cmd); return if($err);
  }
  
  my $blopt="-evalue 1e-9 ";  # for blastp refaa
  $cmd= "$P_BLASTP $blopt -num_threads $NCPU -db $P_REFAA -query $P_AASET -outfmt 7 -out $blout";
  $err= runcmd($cmd); return if($err);
  push @tmpfiles, $blout;

  my $trsizef="$P_REFAA.qual,$P_AASET.qual"; # FIXME
  unless(-f "$P_REFAA.qual") { runcmd("$EVIGENES/prot/aaqual.sh $P_REFAA"); }
  unless(-f "$P_AASET.qual") { runcmd("$EVIGENES/prot/aaqual.sh $P_AASET"); push @tmpfiles, "$P_AASET.qual"; }
  
  $cmd= "$EVIGENES/makeblastscore3.pl -each -needlen -pminlow 0.75 -sizes $trsizef $blout > $ncaabtall";
  $err= runcmd($cmd);
  push @tmpfiles, $ncaabtall;

  my($mok,$ncuniq)= ncaa_match_okaa_btall($outtab, $ncaabtall, $P_OKAABLAST);
  push @tmpfiles, $outtab;
  if($mok and $ncuniq > 0) {
  runcmd("grep 'ncaa.napub' $outtab > $outuniq"); # want this final list of likely prots
  push @tmpfiles, $outuniq;
  }
  
  return ($mok, $outtab, $ncuniq); #? or parse here
}


=item make_allevdtab

    #v2: insert here .. want altpar genetab, blagree/selfdup table, aa.qual + fickett codepot of trset
    # .. allevdtab sorted by dgeneid, tr size
    # $evigenes/prot/cdshexprob.pl -onlyfick -test $trset -out $trset.codepot

    ($nevd,$allevdtab)= make_allevdtab( $genetab, $trlongsizef, trset.codepot, $selfduptab )
    
    evd cols: (trlen aalen pcds cdpot agree dglocus dgclass);

=cut

sub make_allevdtab {
  my($trname, $trset, $genetab, $trlongsizef, $selfduptab, $hasokcds)= @_;

  ## FIX: output has all input.tr set, want only trlongsizef set
  
  (my $allevdtab=$trset) =~ s/.\w+$/.allevd.tab/;
  my $ntrnew= 0;
  my($cmd,$err,$ok,$hin);
  if(not $UPDATEALL and ! -f $allevdtab and -f "ncrnaset/$allevdtab") { $allevdtab= "ncrnaset/$allevdtab"; }
  if( -s $allevdtab and not $UPDATEALL) {
    $ntrnew= `wc -l $allevdtab`; chomp($ntrnew);
    return($ntrnew, $allevdtab); #???
  }
  
  my $testokcds= ($TEST_OKCDS and ref($hasokcds));
  # my $trcodepot="$trset.codepot";
  (my $trcodepot=$trset) =~ s/\.gz//; $trcodepot.=".codepot";
  unless(-f $trcodepot) {
    ## FIXed trset.gz was bad here
    $cmd= "$EVIGENES/prot/cdshexprob.pl -onlyfick -test $trset -out $trcodepot";
    $err= runcmd($cmd);
    push @tmpfiles, $trcodepot;
  }
  
  # genetab == dRnaID	OrigID	dGeneID	AltI	APClass	AAqual TRlen	lAltIDs	lParIDs	lChrLocs
  # .. use 2-6 cols here, last 8-10 cols are unused, may change
  # trlongsizef == trid tlen 0 aaqual tlen cpold ..
  # selfduptab == trid agree duptype dupval dupids
  use constant ADDAAQ => 1;
  our( $ld, %val, $outh);
  #orig# my @kcols=qw(trlen aalen pcds cdpot cdfick agree rdmap rdpair ppair prdcov dglocus dgclass);
  my @kcols=qw(trlen aalen pcds cdpot agree dglocus dgclass);
  push @kcols, 'cdsaln' if($testokcds); #$TEST_OKCDS: add hasokcds col here 
  push @kcols, 'aaqual' if(ADDAAQ); # for output pubids
  
  open( $outh, '>', $allevdtab); push @tmpfiles, $allevdtab; # rename($allevdtab,"$allevdtab.old");
  print $outh join("\t","IDtrasm",@kcols)."\n";  
  
  sub putval{ our( $ld, %val, $outh); my $ok=0;
    if($ld and $val{trlen}) { print $outh $ld; for my $k (@kcols) { my $v=$val{$k}||0; print $outh "\t$v"; }; print  $outh "\n"; $ok=1; } 
    %val=(); return $ok; } 
    
  sub fpkm{ my($r,$w)=@_; my $f= $r / ($w/1000); return  sprintf "%.2f",log(1+$f); }  
  
  # check-read idprefix/num from genetab: ($IDPREncrna,$IDBASEncrna) =("",0);
  ($ok,$hin)= openRead($genetab); my $gtop=0;
  while(<$hin>) { next if(/^\W/); 
    my($pubid,$oid,$gid)=split; last if(++$gtop > 9);
    if($gid =~ m/^(\w+\D)(\d+)$/) {
      my($idpre,$gid)= ($1,$2);  # $IDBASEncrna= $gid; # dont need?
      last if($idpre eq $IDPREncrna);
      if($IDPREncrna and $idpre ne $$IDPREncrna) { } # err
      else { $IDPREncrna= $idpre; }
      }
  } close($hin);
    
  $cmd= "cut -f2-6 $genetab | sort - $trlongsizef $trcodepot $selfduptab |";
  open(EV,$cmd) or die "$cmd";
  while(<EV>) {
    next if(/^(\W|name|dRnaID|Ident|OrigID)/);  
    
    my(@v); 
    my($id,$e,$f,$g)=@v=split; 
    #NOT NOW, have right trdata: next unless($trok{$id});
    
    $ntrnew += putval($ld) if($ld and $id ne $ld);
    
    if(/aln=/){ $val{agree}=$e; } # Yes now: trid agree duptype dupval dupids; dupval: aln=aln%,size%
    ## ^? add hasdup evd score or not here?
    ## elsif(@v==5 and $e=~/^\d/ and) { } # duptab @v == 5
        
    elsif($e =~ /^dg\d/){ $e=~s/^dg//; $val{dglocus}=$e; $val{dgclass}= $g; } # @v == 5 for cut -f2-5 dgclass.tab
    elsif($IDPREncrna and $e =~ /^$IDPREncrna\d/){ $e=~s/^$IDPREncrna//; $val{dglocus}=$e; $val{dgclass}= $g; } # @v == 5 for cut -f2-5 dgclass.tab
    ## ^^ FIXMEd not dg0000 if using -idprefix -geneid opts , but try IDPREncrna ?
    
    elsif($g eq "test"){ $val{cdpot}=$f; }  
    elsif($f eq "test"){ $val{cdpot}=$e; } # @v == 4 or 5
    elsif(@v >= 5 and $e=~/^\d/ and $g =~ /,/){  # tr.qual @v >= 5
      my($id,$tw,$gap,$aaw,$twb)=@v; 
      if($tw < 0.5 * $twb) { $tw=$twb; } # is this .aa.qual tab? use twb ?
      my($aw,$pcd,$ac)=split",",$aaw;  # preserve all of aaw = aaqual for output
      $val{aaqual}=$aaw;
      $val{trlen}=$tw; $pcd=~s/%//; $val{pcds}=$pcd; $val{aalen}=$aw; 
      if($testokcds) {
        my ($bokcds,$okid)= split",", ($hasokcds->{$id}||"0"); # val is bases aligned
        # turn to %align of what? okcds? or this trlen?
        my $pokcds= ($bokcds==0) ? 0 : int(100 * $bokcds/$tw);
        $pokcds="$pokcds,$okid" if($okid);
        $val{cdsaln}= $pokcds; 
        }
      }
      
    # elsif(@v==4 and $e=~/^\d/){ # drop this trqual variant
    #   my($id,$tw,$aaw,$cp)=@v; my($aw,$pcd)=split",",$aaw;  
    #   $val{trlen}=$tw; $pcd=~s/%//; $val{pcds}=$pcd;  $val{aalen}=$aw; 
    #   } 
    # elsif(@v==5 and $e=~/^\d/){ # readmap.tab : now now
    #   my($id,$tw,$rdm,$rdp,$unc)=@v; my $cov=$tw - $unc; 
    #   $val{prdcov}= sprintf "%.2f", 100*$cov/$tw; 
    #   $val{ppair}= ($rdm>0)?int(1000*$rdp/$rdm)/10:0; 
    #   $val{rdmap}= fpkm($rdm,$tw); 
    #   $val{rdpair}= fpkm($rdp,$tw); 
    # }  
    
    $ld=$id; 
  } close(EV);
  
  $ntrnew += putval($ld);  close($outh);
  return($ntrnew,$allevdtab);
}

=item ncrna_makeSelfDupTab

  merged table from self.blast with 
    a. classified self duplicates (imperfectdup + imperfectfrag)
    b. classified valid replication (cross assembler duplicates)
    
  mix of ncrna_replicates() and ncrna_altpar:ncrna_makeSelfDupTab()


  my($ntab,$selfduptab, $havedups,$hasduphash, $nagree)
    = ncrna_makeSelfDupTab($trname, $trlongnomrna, $trselfblast, $trsetagree, $trlongsizef);
    
=cut

sub ncrna_makeSelfDupTab {
  my($trname,$trset,$trselfblast,$trsetagree,$trsizef)= @_; 

  #above# use constant RepFromBtall => 1;  #UPD2222 only this method
  #? my $check_align = ($trsizes and ref $trsizes) ? 1 : 0;#? option

  my $MINPAL=$ENV{blagree_pal} || 0.80; # smaller? larger? options..
  my $MINPRW=$ENV{blagree_prw} || 0.70; #? smaller
  
  (my $selfbtall=$trselfblast) =~ s/\.\w+$/.btall/; # make/use this table for blagree ..
  (my $selfdups=$trselfblast) =~ s/\.\w+$/.sduptab/; # CHANGE output name
  # my $newtable=$selfdups; # NO trsetnew here.  selfdup table
  
  # my $trsetnew="$trname.uniqrna.tr"; # NO trsetnew here.  selfdup table
  my $ntrnew= 0;
  my($cmd,$err,$ok,$hin);
  if(not $UPDATEALL and ! -f $selfdups and -f "ncrnaset/$selfdups") { $selfdups= "ncrnaset/$selfdups"; }
  if( -s $selfdups and not $UPDATEALL) {
    $ntrnew= `wc -l $selfdups`; chomp($ntrnew);
    # ?? need valid %hasdup here
    # return($ntab, $selfdups, $nhasdup, \%hasdup, $nblagree); 
    return($ntrnew, $selfdups); #???
  }

  # step A. selfblast > selfbtall
  # upd altparclass2a() now makes selfbtall before this
  #---------------------
  
  ## Maybe replace this w/ makeblastscore selfbtall table parse .. simplifies decisioh,
  ## use sametr ==  min( align/tlen, align/rlen) > 0.90 and min(rlen/tlen, tlen/rlen) > 0.90
  unless( -s $selfbtall) { #? reuse
    
    unless($trsizef and -f $trsizef) {
      $trsizef="$trset.qual"; 
      unless(-f $trsizef) { my $nts;
        ($nts,$trsizef)= trseq_qual($trname, $trset, undef, $trsizef );  # reduced from long.tr sizef
      }
    }
    
    # change -pminlow ?
    loggit(LOG_DIE,"missing selfblast $trselfblast") unless(-s $trselfblast);
    my $cmd= "$EVIGENES/makeblastscore3.pl -each -needlen -pminlow 0.5 -sizes $trsizef $trselfblast > $selfbtall";
    my $err= runcmd($cmd);
    push @tmpfiles, $selfbtall;
  }
  
  
  # step B... read selfbtall, call selfdups and assembler replicates, write table
  #---------------------
  my($ntab,$nhasdup,%hasdup,%valdup)=(0); # selfdup
  my($ltd,$topal,$nag,$nblagree)= (0) x 9;
  our(%agree,%blagree,%blreplace,%rsn,@rds,%didid,%did); # agree

  sub saveagree { 
    my($ltd,$topal)=@_; 
    our(%agree,%blagree,%blreplace,%rsn,@rds,%didid,%did);
    my @rsn=sort keys %rsn;  my $nput=0;
    my $nagree= @rsn + $agree{$ltd};
    if($nagree > 1) {
      my $pok=1; my $repok=0; # nu
      if(my $oal= $did{$ltd}) {
        my $od= $didid{$ltd};
        my $agr= $agree{$ltd}||0; my $ogr= $agree{$od} || 0;
        my $oblg= $blagree{$od}||0;
        if($agr == 0 and $ogr == 0) { $agr=$nagree; $ogr=$oblg; }
        if( $topal <= $oal and $agr <= $ogr ) { $pok= -1; }
        elsif($topal > $oal and $agr >= $ogr ) { $pok= 2; $repok=1 unless($oblg>$nagree); }
        elsif($topal > $oal and $agr < $ogr ) { $pok= ($agr>1) ? 2 : 0; } 
      }
      if($pok > 0) {      
        # my $rsv= join",",map{ $_.":".$rsn{$_}; } @rsn; 
        # $rsv = "equals" unless(@rsn);
        # my $repval=0; 
        if($repok and my $od=$didid{$ltd}) {  
          delete $blagree{$od}; #? but keep agree{} 
          $blreplace{$od}=$ltd; # dont need w/ delete ?
          # $repval="replace=$od"; dont need 
        } 
        ### print $agreeout join("\t",$ltd,$nagree,$rsv,"aln=$topal",$repval)."\n"; 
        $blagree{$ltd}=$nagree; $nput++; # local dont bother w/ agreeout
        $didid{$ltd}=$ltd; $did{$ltd}=$topal; 
        for my $d (@rds) { if($did{$d} < $topal) { $didid{$d}=$ltd; $did{$d}=$topal; }  }
      }
    } 
    return($nput);
  }
  #------------------------------

  if(-f $trsetagree) {
    ($ok,$hin)= openRead($trsetagree);
    while(<$hin>) { next if(/^\W/);  my($td,$nagr)=split;  $agree{$td}=$nagr; $nag++; }
    close($hin);
    loggit(0,"nagree=$nag in $trsetagree");
  }

  ($ok,$hin)= openRead($selfbtall);
  loggit(0, "ncrna_makeSelfDupTab read $selfbtall, ok=$ok");
  while(<$hin>) { 
    next if(/^\W/);
    my($td,$rd,$bs,$ida,$aln,$tw,$rw)= split;
    next unless($tw>0 and $rw>0 and $td ne $rd);  
    # next if($td eq $rd or $hasdup{$rd}{$td}); # no self or circular rd/td and td/rd
 
    ## self-dup code
    my $palt= $aln/$tw; my $prwt= $tw/$rw; # for td isDup/Frag
    if($prwt <= 1.0 and $palt >= $MINDUPALN and not $hasdup{$rd}{$td}) { 
      my $dupflag= ($prwt < $MINDUPFULL)? "frag" : "dup";
      $hasdup{$td}{$rd}= $dupflag; $nhasdup++;
      $valdup{$td}{$rd}= sprintf"%d,%d", int(100*$palt), int(100*$prwt); # append pal,prw ?
    }
   
    ##  agree repl code from above
    if($td ne $ltd) { $nblagree += saveagree($ltd,$topal) if($ltd); %rsn=@rds=(); } 
    $ltd=$td;  
    
    my($pala,$prwa) = ($tw > $rw) ? ($aln/$tw, $rw/$tw) : ($aln/$rw, $tw/$rw); # for agree
    ## cancel this skip if $agree{$td} .. but dont add rd rsn
    my $okagree = ($pala < $MINPAL or $prwa < $MINPRW) ? 0 : 1;
    if($okagree) {
      push @rds,$rd;
      $topal= $aln if($aln > $topal);  
      foreach my $sd ($td,$rd) {  
        my $st= $sd;
        if($st=~s/[Ll]oc.*$//) { $st =~ s/k\d+$//; } 
        else { $st=substr($sd,0,9); } #? cancel if no [Ll]oc ? dont want spurious agree measure; need option
        $rsn{$st}++;
      }
    }
   
  } close($hin);
  $nblagree += saveagree($ltd,$topal) if($ltd); 
  
  my %td=(); map{ $td{$_}=1 } (keys %hasdup, keys %blagree);
  my @td= sort keys %td;
  $ntab= @td; my $dhdr=0;
  loggit(0, "ncrna_makeSelfDupTab n=$ntab to $selfdups");
  open(O,'>',$selfdups) or loggit(LOG_DIE,"write $selfdups");
  push @tmpfiles, $selfdups;
  for my $td (@td) { 
    my @did= sort keys %{$hasdup{$td}};
    my $agree= $blagree{$td} || $agree{$td} || 0; # blagree should contain nrident agree set
    next unless(@did > 0 or $agree > 0);  # BUG above? @did == 0 from hasdup rd,td and td,rd

    my $dupflag=""; my $valdup=0;
    for my $d (@did){ $dupflag= $hasdup{$td}{$d}; $valdup=$valdup{$td}{$d}; last if($dupflag eq "frag"); } 
    my $dupids= join",",@did;
    map{ $_ ||= "0" } ($dupflag,$dupids); # ?data bug
    
    print O join("\t",qw(IDentifier agree duptype dupval dupids))."\n" unless($dhdr++); 
    print O join("\t",$td,$agree,$dupflag,"aln=".$valdup,$dupids)."\n"; 
  } close(O);
  
  return($ntab, $selfdups, $nhasdup, \%hasdup, $nblagree); #? \%blagree,\%hasdup
}

=item ncrna_align_okcds

  ($ncds,$hasokcds,$okcdsblastn)= ncrna_align_okcds($trname,$trset,$okcds);
  
  which way?  -db uniqrna -query okcds or rev?
   makeblastdb -dbtype nucl -in $trname.uniqrna.tr -logfile /dev/null
   bloptnu="-perc_identity 98 -evalue 1e-9 -dust yes "
   blastn $bloptnu -num_threads $ncpu -db $trname.uniqrna.tr -query $okcds -outfmt 7 -out $trname.uniqrna.okcds.blastn

=cut 


sub ncrna_align_okcds {
  my($trname,$trset,$okcds)= @_; # trset=="$trname.uniqrna.tr";
  return (0) unless(-f $okcds);
  my ($ncds,%hasokcds)= (0); 
  my $blout="$trset.okcdsaln.blastn";  
  my $pctident=98; # -pHetero fixme
  
  sub countokcds { 
    my($blout)=@_; 
    open(B,$blout) or return(0); %hasokcds=();
    while(<B>){ if(/^\w/) { my($td,$rd,$pid,$aln)=split; 
      if($aln >= $MIN_CDSALN){ my $oal=$hasokcds{$td}||0; 
        if($aln>$oal){ $hasokcds{$td}="$aln,$rd"; $ncds++; } } 
        } 
    } close(B);
    return scalar(keys %hasokcds);
  }
  
  if(not $UPDATEALL and ! -f $blout and -f "ncrnaset/$blout") { $blout= "ncrnaset/$blout"; }
  if( -s $blout and not $UPDATEALL) { 
    $ncds= countokcds($blout);
    return($ncds,\%hasokcds,$blout);
  }

  ## add here? class b/n truniq matching okay.cds, truniq.okcds.tr, and those w/o truniq.ncrna.tr
  # (my $okcds=$mrna) =~ s/.mrna/.cds/;
  # my($hasokcds)  = ncrna_align_okcds($trname,$truniqrna,$okcds);

  my($cmd,$err,$ok,$hin);
  if($okcds =~ /\.gz$/){ 
    my $okcdsgz=$okcds; $okcds =~ s/.gz$//;
    $cmd="gunzip -c $okcdsgz | $APPmakeblastdb -dbtype nucl -title $okcds -out $okcds -logfile /dev/null";  
  } else {
    $cmd="$APPmakeblastdb -dbtype nucl -in $okcds -logfile /dev/null"; # .gz NOT ok 
  }
  push @erasefiles, map{ "$okcds.$_" } qw(nsq nin nhr);
  $err= runcmd($cmd); return if($err);

  #x my $bloptnu="-perc_identity 98 -qcov_hsp_perc 90 -dust no "; #? check opts again, maybe no -qcov, -dust yes
  my $bloptnu="-perc_identity $pctident -evalue 1e-9 -culling_limit 1";  # ? culling_limit 1 need only 1 hit
  $blout="$trset.okcdsaln.blastn"; push @tmpfiles, $blout;
  $cmd= "$APPblastn $bloptnu -num_threads $NCPU -db $okcds -query $trset -outfmt 6 -out $blout";
  $err= runcmd($cmd); return if($err);
  
  $ncds= countokcds($blout);
  return($ncds,\%hasokcds,$blout);
}


=item altparclass algo

  if [ -f $trname.aa.qual ]; then
    perl -ne '($d)=@v=split; if(@v==1) { $ok{$d}=1; } else { print if($ok{$d}); }' \
      $trname.self98.blagree.trids $trname.aa.qual  > $trname.uniqrna.aa.qual
    trsizes=$trname.uniqrna.aa.qual
  else
    env ismrna=1 $evigene/scripts/prot/aaqual.sh $trname.uniqrna.tr
    trsizes=$trname.uniqrna.tr.qual
  fi
  
  $evigene/scripts/genes/altparclassify.pl -ncpu $ncpu -cds $trname.uniqrna.tr -sizes $trsizes
  
  # 3b. FIXME: step2 leaves in notokmrna with large overlap to okmrna : need to separate these from ncrna
  #    step1 removes only contained-in-okmrna subset
  #    use prior step2 + self-agree, add after step3 ? blastn -db uniqrna.tr -query ok.cds -qcov 60-90% ?
  
   okcds=`echo $mrna | sed 's/.mrna/.cds/;'`
   if [ -f $okcds ]; then 
   makeblastdb -dbtype nucl -in $trname.uniqrna.tr -logfile /dev/null
   bloptnu="-perc_identity 98 -evalue 1e-9 -dust yes "
   blastn $bloptnu -num_threads $ncpu -db $trname.uniqrna.tr -query $okcds -outfmt 7 -out $trname.uniqrna.okcds.blastn
   fi

=cut

sub altparclass {  
  my($trname,$trset,$trselfblast,$trsizef,$hasokcds,$evgapp)= @_;  ## trset == truniqrna; trsizef was aaqualf
  my($ngene,$cmd,$err)= (0);
  
  ## add altparclassify -output $trclass option
  (my $cname= $trset) =~ s/\.\w+$//;
  my $genetab= "$cname.dgclass.tab"; push @tmpfiles, $genetab;  
  if(not $UPDATEALL and ! -f $genetab and -f "ncrnaset/$genetab") { $genetab= "ncrnaset/$genetab"; }
  if(-s $genetab and not $UPDATEALL) {
    ($ngene)= `cut -f5 $genetab | egrep -c '^(main|uni)'`; chomp($ngene);
    return($ngene,$genetab);
  }
  
  # FIXME: trselfblast uses trsizef to keep, but now
  #   trsizef == longnok.tr.qual has more than trset == uniqrna.tr
  my($nts,$tsok)= (0,0);
  if(-s $trsizef and ($trsizef eq "$trset.qual")) {
    my $trids= faidlist($trset,undef,'hash');
    $nts= 0; open(F,$trsizef); while(<F>){ $nts++ if(/^\w/); } close(F);
    $tsok= ($nts == scalar(keys %$trids)) ? 1 : 0;
  }
  unless( $tsok ) {
    ($nts,$trsizef)= trseq_qual($trname, $trset, undef, $trsizef );  # reduced from long.tr sizef
  }
    
  return(0) if($nts < 9);
  
  # use %hasokcds{trid} to reclass as miscRNA vs ncRNA .. miscRNA == pseudogene, missed-mRNA, etc
  
  # altparclassify.pl uses this, use same here? change?
  # my $NCBLAST_IDENT=  $ENV{pctident} || 97; # option .. 95 in altparc
  # my $NCBLAST_EVALUE= $ENV{evalue} || 1e-9; # option .. 1e-5 in altparc
  #NOTE: $REUSE_SELFBLAST reusing $trselfblast which has more ids than trset, added -requiresize
  ## altparclassify bug: TRlen col NOT trlen but aasize*3
  
  $evgapp ||= "$EVIGENES/genes/altparclassify.pl";
  $cmd="$evgapp -ncpu $NCPU -cds $trset -sizes $trsizef";
  $cmd .=" -debug" if($debug);
  $cmd .=" -requiresize -selfblast $trselfblast" if($trselfblast and -f $trselfblast);
  # my($IDPREncrna,$IDBASEncrna) =("",0);
  ##?? BAD for above parsing, expecting 'dg0000' ids only from alrparclass.tab
  $cmd.=" -idprefix $IDPREncrna" if($IDPREncrna); # UPD2222
  $cmd.=" -geneid $IDBASEncrna" if($IDBASEncrna);

  $err= runcmd($cmd);
  # note this makes $trset.aa.qual, trset.self95.blastn
  
  ($ngene)= `cut -f5 $genetab | egrep -c '^(main|uni)'`; chomp($ngene);
  # my ($classcount)= `cut -f5 $genetab | sort | uniq -c `;
  # create uniqrna.main.tr from genetab main|uni ?
  
  return($ngene,$genetab);
}

sub altparclass2a {  #UPD2222
  my($trname,$trset,$trselfblast,$trsizef,$hasokcds)= @_;  ## trset == truniqrna; trsizef was aaqualf
  if(0) {
    return altparclass(@_); #v1 looks okay for v2 ... 
  }

  # FIXME -pHetero reduce altparclassify -altpident 99.0 -pctident 95
  
  my $evgapp="$EVIGENES/genes/altparclassify2b.pl"; # TEST, replace orig if ok
  # orig is mem-pig for largish self-blast, too many internal hashes of all that
  # upd 2b uses self.btall sorted by longest tr, same as made w/ makeselfdup,
  ## added opt: altparclassify2b -output $genetab

  my($ngene,$cmd,$err)= (0);
  (my $cname= $trset) =~ s/\.\w+$//;
  my $genetab= "$cname.dgclass.tab"; push @tmpfiles, $genetab;  
  if(not $UPDATEALL and ! -f $genetab and -f "ncrnaset/$genetab") { $genetab= "ncrnaset/$genetab"; }
  if(-s $genetab and not $UPDATEALL) {
    ($ngene)= `cut -f5 $genetab | egrep -c '^(main|uni)'`; chomp($ngene);
    return($ngene,$genetab);
  }
  
  my($nts,$tsok)= (0,0);
  unless($trsizef and -f $trsizef) {
    $trsizef="$trset.qual"; 
    unless(-f $trsizef) { 
      ($nts,$trsizef)= trseq_qual($trname, $trset, undef, $trsizef );  # reduced from long.tr sizef
      return(0) if($nts < 9);
    }
  }

  loggit(0, "altparclass2a start to $genetab");
  # same as makeSelfDupTab: step A. selfblast > selfbtall
  # trselfblast should be valid file here
  #---------------------
  loggit(LOG_DIE,"missing selfblast $trselfblast") unless(-s $trselfblast);
  (my $selfbtall=$trselfblast) =~ s/\.\w+$/.btall/; # make/use this table for blagree ..
  unless( -s $selfbtall) { 
    my $cmd= "$EVIGENES/makeblastscore3.pl -each -needlen -pminlow 0.5 -sizes $trsizef $trselfblast > $selfbtall";
    my $err= runcmd($cmd);
    push @tmpfiles, $selfbtall;
  }
  
  $cmd="$evgapp -ncpu $NCPU -cds $trset -sizes $trsizef";
  $cmd .=" -debug" if($debug);
  $cmd .=" -btallself $selfbtall"; #  if($selfbtall and -f $selfbtall); # NOT -sortedbtall, presumably
  #o# $cmd .=" -requiresize -selfblast $trselfblast" if($trselfblast and -f $trselfblast);
  $cmd.=" -idprefix $IDPREncrna" if($IDPREncrna); # UPD2222
  $cmd.=" -geneid $IDBASEncrna" if($IDBASEncrna);
  $err= runcmd($cmd);
  if($err) {
    #runcmd logs err: loggit(0,"ERR: $err, ");
  } else {
    ($ngene)= `cut -f5 $genetab | egrep -c '^(main|uni)'`; chomp($ngene);
  }
  return($ngene,$genetab);
}


=item run_evgtr2ncrna.sh

  #! /bin/bash
  ### env trset=inputset/name.tr  mrna=okayset/name.okay.mrna  ncpu=12 datad=path/to/data  qsub -q normal run_evgtr2ncrna.sh
  #PBS -N evg_tr2ncrna
  #PBS -A PutAccountIdHere
  #PBS -l nodes=1:ppn=16,walltime=39:55:00
  #PBS -V
  
  # run_evgtr2ncrna.sh = new evigene/scripts/genes/tr2ncrna.pl
  # inputs: okayset/name.okay.mrna inputset/name.tr
  # v1 opt:  -REUSE_SELFBLAST
  
  # env need: $evigene/scripts/ and blastn on path
  evgapps=$HOME/bio/apps
  evigene=$evgapps/evigene
  export PATH=$evgapps/ncbi/bin:$PATH
  export PATH=$evgapps/exonerate/bin:$PATH
  evgapp=$evigene/scripts/genes/tr2ncrna.pl
  
  if [ "X" = "X$ncpu" ]; then ncpu=8; fi
  if [ "X" = "X$maxmem" ]; then maxmem=64000; fi
  if [ "X" = "X$datad" ]; then echo missing datad=/path/to/data; exit -1; fi
  if [ "X" = "X$mrna" ]; then echo missing mrna='name.mrna'; exit -1; fi
  if [ "X" = "X$trset" ]; then echo missing trset='name.tr'; exit -1; fi
  
  cd $datad/
  echo "#START `date` " 
  
  echo $evgapp -debug -log  -ncpu $ncpu  -mrna $mrna -trset $trset
  $evgapp -debug -log  -ncpu $ncpu  -mrna $mrna -trset $trset
  
  echo "#DONE : `date`"

=cut

__END__

