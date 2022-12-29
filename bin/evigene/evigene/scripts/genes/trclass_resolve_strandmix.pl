#!/usr/bin/env perl
# evigene/genes/trclass_resolve_strandmix.pl
# parts from altreclass.pl : evigene alt-reclassifier

=item about 

  EvidentialGene trclass_resolve_strandmix 
    Examines trclass to reorient mixed-strand alternates
    from evidence of geneval/chr-map.align.tab (stronger) and/or code potential (weaker)
    
  Usage: trclass_resolve_strandmix.pl -trclass mymrna.trclass 
   opts: 
    -pubids myrna.pubids : use this pubids instead of publicset|okayset pubids
    -output mymrna.reorient_pubids : print pubids to this file not STDOUT
    -calccodepot : calculate codepot from mRNA (CDS - UTR)
    -fragcheck   : blast-align new rev CDS to gene set for fragment detection
    -allpubids   : print unchanged + changed pubids
    
   Expects Evigene data, from SRA2Genes or tr2aacds, with trclass, pubids and okayset sequences,
   geneval/chr-map.align.tab is table from SRA2Genes step 9 chr-mapping.
   Output reoriented sequences to reorset/ with updated okayset/ as integrated with tr2aacds4.
   
   Tests on reference transcripts show codepot-only corrects more than it errs (human 10:1, plant 2:1),
   while chr-map.align corrects those with CDS-oriented introns.   Ambiguous cases are tabulated
   and CDS/AA of both ambiguous strands should be further tested with homology evidence.

   The coding gene locus paradigm says that all alternates of one locus should be oriented
   in same direction.  There is a biological subset of CDS where longest ORF calculations 
   return a longer false ORF on reverse strand overlapping a true shorter ORF.  Homology to 
   known reference proteins can correct this false reverse ORF, as can alignment to chromosomes 
   with intron orientation.  False reverse ORFs tend to look more like non-coding amino sequence 
   (lower code potential) but this is a weaker test.
   
   Some paralogs appear as mixed strand alternates of one locus, a problem to keep in mind. 
   Some of the  paralogs in reference gene sets do indeed have reversed ORF orientations,
   which may be supported by experimental evidence, or may be computational reverses.
  
=item original usage

  export PATH=$PATH:$ncbiblastbin # blastn for -fragcheck
  pt=human18ncx
  
  $evigene/scripts/genes/trclass_resolve_strandmix.pl \
   -debug -log -fragcheck -changesonly \
   -trclass $pt.trclass -out $pt.restrand_pubids 

  export DROPSHORT_ANTISENSE=0.01 # can be reset to keep short ref seq

=item draft trclass_resolve_strandmix.pl
  
  a. read trclass, pubids 
  b. collect genes with mixed strands among alt/pars
  c. pick proper strand from largest codepot among those,
     modified by other quals?  aaref homol score, aaqual part/full?, aacons? no
  d. reverse seqs of false strand subset,   
  e. rewrite data: pubids, okayset/prepub seqs, other 2ndary like aa.qual, any tmpfiles?
     .. retain? false strand seqs as .cullrev set?
     
  resolution method: use pubids gene groups, for each group 2+ trs, check trclass for -sense flag,
  if -sense flag, check all aaqual codepot, pick one sense for largest codepot, with caveats? 
  quals aaref* if have it, aasize/complete/aaqual? aacons? no
  then swap revaa for -sense trs, changing okayset seqs
  
  for tr2aacds4_stage2b .. after 1st trclass2pubset ? before blasttrset2exons ?
  *** redo cdsqual, etc for cds and strand changes **

=item update 2020apr

  added dupcds_filter : now removes dup.cds from reor
  added -strandedrna option, and auto-detects from aaqual table
    -strandedrna=no|0 means ignore auto-detect of stranded mRNA, 
    default now is detect and use to break ambiguous orient_r2f/f2r
  added -nomaybe option, pick 1 of ambiguous fwd/rev cds
  
  as tested:
  $evigene/scripts/genes/trclass_resolve_strandmix.pl -debug -nomaybe -strandedrna=1 -calccodepot  -NOfragcheck \
    -log -out -trclass human18ncx.trclass

  options/defaults
    -nomaybe : keep as option? how from tr2aacds?
    -strandedrna : default ok
    -calccodepot : should become default?
    -fragcheck/-nofragcheck : which? adds blastn test useful in reducing excess
    -output ?? is for publicset/name.reor_pubids only, make default .. is this near-obsolete?
    
  $evigene/scripts/genes/trclass_resolve_strandmix.pl -debug -calccodepot -[no]fragcheck -log -trclass $pt.trclass
        
=cut

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/.."); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
use strict;
use Getopt::Long;
use cdna_evigenesub;  
use evigene_pubsets;
# use cdna_proteins;
use File::Basename qw(basename dirname fileparse);

use constant VERSION =>  '2020.04.20'; # dupcds_filter; 
# '2020.01.12'; # UPD1912
use constant UPD1912 => 1;
use constant ANTI1912c => 1; ## antisense checks of asmrna_dupfilter4a.pl

our $EVIGENES= $FindBin::Bin; $EVIGENES =~ s,scripts/.*,scripts,;
our $EGAPP='reorient';  # .. restrand
our $EGLOG='reor';

my($SAME_PI,$SAME_PA,$DIFF_PA)= (100,99,90); # samecds= $pi>=98 and $pa>=97; diffcds= $pa<=80
my $TRUSTPUBIDS= 1; # $ENV{trustpubids}||0; # trust input pubids for class, aaref calls, ie newer than trclass
my $debug=$ENV{debug}||0; 
my $tidyup=1;
my $PUTONLYCHANGES= 1; # $debug; # make this default on, -nochange to get all?
my $putAllPubids= 0; # reverse PUTONLYCHANGES
my $DO_COMPARE2ORIG= $ENV{fragcheck}||0; # runs blastn -db orig.cds -query revnew.cds > perffrag.idtab

# from asmrna_altreclass3c.pl  # DROPSHORT_ANTISENSE below longest fwd aa are cull-able revaa alts
# my($DROPSHORT_AAPART,$DROPSHORT_AAFULL,$DROPSHORT_ANTISENSE, $DROPTINYALT,$DROPSHORT_BADEXONS)
#   = (0.70,0.20,0.60,0.20,0.20); # UPD1908
my $DROPSHORT_ANTISENSE= $ENV{DROPSHORT_ANTISENSE} || 0.20; # some are paralogs; was 0.60 old $DROPSHORT_ANTISENSE;
my $PAL_ANTIMIN= $ENV{PALMIN_ANTISENSE} || 60; # reject antisense align below this, ANTI1912c opt from asmrna_altreclass3c, test w/ refs
my $RNAIsStranded= 0; #UPD20ap
#tests ok# my $noDUPSEQUPD= $ENV{nodupsequpd}||0;  #UPD20ap test need of this .. seems we dont need it
my $NOMAYBE= 0; #UPD20ap
my $AAMIN = $ENV{aamin} || 30; # fixme; was $ENV{MINAA}; tr2aacds4 sets ENV{aamin}, asmdupfilter uses ENV{aamin}

use constant { HNoncode => -1e-4, HCoding => 1e-4 };  # hexamer score, dunno what range to call
my $USE_CODEPOTAB= $ENV{calccodepot}||0; # replace aaqual val, temp test calc w/ evigene prot/cdshexprob.pl

my ($trclass,$pubids,$output,$aaqualf,$mapsensetab,$logfile);
my @ARGSAVE=@ARGV;

my $optok= GetOptions(
  "class|trclass=s", \$trclass, # require
  "pubids=s", \$pubids, # require, find 
  "aaqualf=s", \$aaqualf, # aaqualf required ? w/ codepot; add mrna/cds/aa seq set??
  "output:s", \$output, # option
  "logfile:s", \$logfile, # option
  # "trustpubids!", \$TRUSTPUBIDS, 
  # "revisepublicset!", \$REVISEPUBFILES, 
  "changesonly!", \$PUTONLYCHANGES,  #? make default, yes: -nochanges turns off is bad sense
  "allpubids!", \$putAllPubids,  #? reverse default -changes 
  "calccodepot!", \$USE_CODEPOTAB,   
  "nomaybe", \$NOMAYBE,   # decide on 1 orientation, no maybe both ways; or send maybe reor seq to separate outputs
  "compare2orig|fragcheck!", \$DO_COMPARE2ORIG, 
  "dropshort_antisense=s", \$DROPSHORT_ANTISENSE,  
  "minalign_antisense=s", \$PAL_ANTIMIN,  
  "strandedrna:s", \$RNAIsStranded, # UPD20apr
  "debug!", \$debug, 
  );
## should be opts? DROPSHORT_ANTISENSE PAL_ANTIMIN

die "EvidentialGene trclass_resolve_strandmix VERSION ",VERSION,"
  Examines trclass to reorient mixed-strand alternates
  from evidence of geneval/chr-map.align.tab (stronger) and/or code potential (weaker)
  
Usage: trclass_resolve_strandmix.pl -trclass mymrna.trclass 
 opts: 
  -pubids myrna.pubids : use this pubids instead of publicset|okayset pubids
  -output mymrna.reorient_pubids : print pubids to this file not STDOUT
  -calccodepot : calculate codepot from mRNA (CDS - UTR)
  -fragcheck   : blast-align new rev CDS to gene set for fragment detection
  -allpubids   : print unchanged + changed pubids
  -dropshort_antisense=$DROPSHORT_ANTISENSE [0.01-0.99]: drop antisense alternate below proportion of longest, as too small
  -minalign_antisense=$PAL_ANTIMIN [1-99]: ignore antisense below this % alignment, as likely paralog
  
 Expects Evigene data, from SRA2Genes or tr2aacds, with trclass, pubids and okayset sequences,
 geneval/chr-map.align.tab is table from SRA2Genes step 9 chr-mapping.
 Output reoriented sequences to reorset/ with updated okayset/ as integrated with tr2aacds4.
 Tests on reference transcripts show codepot-only corrects more than it errs (human 10:1, plant 2:1),
 while chr-map.align corrects those with CDS-oriented introns.   Ambiguous cases are tabulated
 and CDS/AA of both ambiguous strands should be further tested with homology evidence.
 \n"
  unless($optok and ($trclass ));  # and $pubids
#   -changesonly : print only changed pubids

our $DEBUG= $debug;  # set package dbg
$PUTONLYCHANGES=0 if($putAllPubids);
$DROPSHORT_ANTISENSE= $DROPSHORT_ANTISENSE/100 if($DROPSHORT_ANTISENSE >=1) ;

openloggit($logfile,$trclass);
loggit(1, "EvidentialGene trclass_resolve_strandmix, VERSION",VERSION);
loggit 1, "trclass_resolve_strandmix  @ARGSAVE \n";
#loggit(1, "ERR: unused arguments:",@ARGV) if(@ARGV>0);

my $APPcdnabest= findevigeneapp("cdna_bestorf.pl"); # allow ENV/path substitutions? ? UPD 2016.07? or cdna_proteins.pm
my ($APPblastn,$APPmakeblastdb);
#   $APPblastn= findapp("blastn"); # defer may not want

#FIXME: use pubid or oid as primary id? trclass has only oids, aaqual also?, pubids has both
my (%trscore,%trv,%trline,%pubid,%newid); # globals now for reclassAlts/putgene
my ($trinfo,$pubinfo,$aaqualhash,$mapinfo, $outvalsh) = (undef) x 9;
my ($okaysetd,$pubsetd); 
our (@okayset,@dropset,@inputset,@tmpfiles,@erasefiles,@publicset); # tidyup file sets

my $outputh=undef; # *STDOUT; # $outh
my $newortable= {}; # FIXME
my $codepotab= {}; # replaces aaqualhash codepot col
my ($APPcdshexprob,$cpotdata)=("","");

MAIN_resolveStrandMix();

=item openRead(STDOUT) bug : fixed now
  still get this unknow err at end of prog .. twice, not seen other progs
  fixed with $outputh=undef instead of =*STDOUT above
  
  Filehandle STDOUT reopened as $hin only for input at /Users/gilbertd/Desktop/genowork/genes2/evigene/scripts/genes/../cdna_evigenesub.pm line 84.
  Filehandle STDOUT reopened as $hin only for input at /Users/gilbertd/Desktop/genowork/genes2/evigene/scripts/genes/../cdna_evigenesub.pm line 84.

=cut
  
sub MAIN_resolveStrandMix {
  loggit(0, "BEGIN with input=",$trclass,"date=",`date`);

  # my $trname= $trclass; $trname=~s/.\w+$//; # .trclass
  my($trname,$trpath,$trsuffix) = fileparse($trclass,qr/\.\w*/); # suf-'.trclass' or suf = qr/\.\w*/
  # trc=`pwd`/iscaplacz.trclass => iscaplacz,/bio/bio-grid/aabugs4/ticks/iticklymec7hetf/,.trclass
  # trc=iscaplacz.trclass => iscaplacz,./,.trclass
  chdir($trpath) unless($trpath eq './' or $trpath eq "");
  
  # need to work with two fileset branches?  stg1=o:okayset,p:publicset stg2=o:okayset1st,p:okayset
  # o: has orig ids, orig rna orient, p: has public ids, mrna orient,
  # .. resolve strands on o: set, using pubid loci;  update p: set or redo p: set
  my @okpubs1a= qw( okayset publicset);
  my @okpubs2b= qw( okayset1st okayset);
  ($okaysetd,$pubsetd)= @okpubs2b; 
  unless( -d $okaysetd and -d $pubsetd ) {
    ($okaysetd,$pubsetd)= @okpubs1a;
    unless( -d $okaysetd and -d $pubsetd ) {
      loggit(LOG_DIE,"Cannot find evigene data directories: @okpubs1a");
    }
  }
  
  # locate data associated w/ trclass
  unless($pubids) { 
    ## default evg file set; okayset/pubids for tr2cds4 n
    # ($pubids=$trname) =~ s,(\w+)$,publicset/$1.pubids,; 
    # $pubids =~ s,publicset,okayset, unless(-f $pubids);
    # $pubids =~ s,okayset/,, unless(-f $pubids);
    #?? fixme here, really want pubids.old if have run all tr2cds4, new pubids from altreclass lacks data
    $pubids="$pubsetd/$trname.pubids";
    if(-f $pubids and -f "$pubids.old") { $pubids= "$pubids.old"; }
    }
  unless($aaqualf) { 
    $aaqualf="$okaysetd/$trname.aa.qual"; # which 1st? unless(-f $aaqualf);
    $aaqualf="inputset/$trname.aa.qual" unless(-f $aaqualf);
    $aaqualf="$trname.aa.qual" unless(-f $aaqualf);
    # unless(-f $aaqualf) {
      # my($okaa,$altaa,$okd)= getOkFileset( $okaysetd,'aa',undef,$trname); # okayset1st or okayset
      # if($altaa) { } # make okall.aa ?
      # ($aaqualf)= makeAaQual($okaa, 1);
    # }
    }
  unless($mapsensetab) { 
    $mapsensetab="geneval/$trname.align.tab";
    $mapsensetab="geneval/$trname.align.tab.gz" unless(-f $mapsensetab);
    $mapsensetab="$trname.align.tab" unless(-f $mapsensetab);
    # ($mapsensetab=$trname) =~ s,(\w+)$,geneval/$1.align.tab,; 
    ## FIXME? what is name of this.align.tab w/ chr name?
    # $mapsensetab =~ s,geneval/,, unless(-f $mapsensetab);
    }

  loggit(LOG_DIE,"Cannot find pubids $pubids") unless(-f $pubids);
  loggit(LOG_DIE,"Cannot find aa.qual $aaqualf") unless(-f $aaqualf);# regenerate from makeAaQual(seqset.aa) ?

  my $makeoutpubids= (defined $output) or (not $debug); # always print to file now, unless testing?
  if($makeoutpubids and not $output) {
    ($output=$pubids) =~ s/\.pubids.*//; $output.=".reorient_pubids";
    #?? got this from where? okayset/human18ncx.old.reorient_pubids ^^ from $pubids= "$pubids.old";
    #?? outname okayset/human18ncx.old.reorient_pubids not good
  }
  
  if($output) {
    rename($output,"$output.old") if(-f $output);
    open($outputh,'>',$output) or die "ERR: writing $output"; 
  }
  loggit 0,"resolveStrandMix( $pubids, $trclass)\n";
  
  ($trinfo,$aaqualhash)= readTrClass($trclass);  # require trclass here
  # elsif($pubids) { ($trinfo,$aaqualhash)= readPubidClass($pubids,$trinfo); }  
  
  # ($aaqualhash)= readAaQual( $trclass||$pubids, $trinfo, $aaqualhash); # required w/ codepot column
  ($aaqualhash)= getAaQual($aaqualf,$aaqualhash); # cdna_evigenesub;

  #UPD20apr:
  # if($rnaIsStranded) ... redo get_bestorf(option=forward) ? use with -reorient
  # NOTE: isStrandedRNA() true if input is merge of subset mRNA, though subset can be unstranded data
  # use isStrandedRNA() test only if flagged? or use unless opt -stranded=no, default use $ismrna test
  # opt -stranded=0|no means dont use ismrna test; -stranded=1 do use ismrna;
  
  if(defined $RNAIsStranded) { $RNAIsStranded ||= 1; } # what?
  my($ismrna,$nfwdstr,$nrevstr,$ntstr)= isStrandedRNA($aaqualhash,$aaqualf);
  #x $RNAIsStranded= ($RNAIsStranded eq "0" or $RNAIsStranded eq "no")? 0 : ($ismrna and $RNAIsStranded) ? 1 : 0;
  $RNAIsStranded= ($RNAIsStranded eq "0" or $RNAIsStranded eq "no")? 0 : ($ismrna) ? 1 : 0;
  loggit 0, "isStrandedRNA=$RNAIsStranded, f:$nfwdstr/r:$nrevstr\n";
  
  if($USE_CODEPOTAB) {
    my $ncpot=0;
    ($ncpot,$codepotab)= calcCodepotOfSeqset($trname,$pubids,$okaysetd,$pubsetd,$aaqualf); # fill in \%codepotab{id}=llcoding
    $USE_CODEPOTAB=0 unless($ncpot>0);
  }
  
  # ($mapinfo)= readMapsenseTab($mapsensetab); # update to add more gmap quals
  # $ENV{MINCOV_ANTISENSE}= 50; # for readAlignTab() w/ goodish cds antisense data
  my($nmapqual, $mapqualh, $alntabh)= ( -f $mapsensetab ) ? readAlignTab($mapsensetab) : (0,0,0);
 
  # also set for restrand_seqs: newortable
  unless(ref $outputh) { $outputh=*STDOUT; }
  
  my ($nin,$nout,$nrenum,$ngdiff,$ngene)
      = resolveStrandMix($outputh, $pubids, $trinfo, $aaqualhash, $mapqualh); 
       
  close($outputh) if($output);
  loggit 0, "resolveStrandMix nin=$nin, nout=$nout, ngdiff=$ngdiff, ngene=$ngene\n";
  
  my($nupdate, $reoridtab, $reorfiles)= restrand_seqs( $trname, $pubids, $output);
  loggit 1, "restrand_seqs n=$nupdate\n" ;
  #  return ($nupd, \%reoridtab, \@outfiles); # = $nupd, \%updidtab, [$fixrna, $fixaa, $fixcds, $fixidfile]

  ## update_mrna_fileset() wont work w/o fileset changes, one set of mrna,aa,cds of fwd+rev seqs
  ## $reoridtab->{$id}=$action;
  # my($upstatus, $upfiles, $uptemp, $upokids)  
  #  = update_mrna_fileset($trpath, $mrnaseq, 'restrand_done', $reoridtab, @$reorfiles) if($nupdate>0);

  if($nupdate > 0) {
    my($nupc,$updclass)= rewriteTrClass( $trclass, $reoridtab);
    if($nupc>0){
      my $sf= (-f "$trclass.orig")?'old':'orig';
      rename($trclass,"$trclass.$sf"); # push @tmpfiles, "$trclass.$sf";
      rename($updclass,$trclass);
      loggit 1, "rewriteTrClass $updclass to $trclass\n" ;
    }
  }
  
  my $tmpfolder="reorset";
  my $outfolder= $okaysetd; # source files are here, need to merge reorfiles still
  if( $tidyup and $nupdate > 0 ) {  
 
    tidyupFileset($outfolder,@$reorfiles) if(@$reorfiles);  
    tidyupFileset($tmpfolder,@tmpfiles) if(@tmpfiles);  
    
    my @rmlist;
    foreach my $fn (@erasefiles) { if(-f $fn){ unlink($fn); push @rmlist,$fn; }}  #($debug) ? loggit(0,"rm $fn") : unlink($fn);
    if(@rmlist) { my $nrm=@rmlist; my $rml=join" ",@rmlist[0..4]; loggit(0,"tidyup erase: n=$nrm, $rml .."); } 
  }
  loggit(0, "DONE at date=",`date`);
  loggit(0, "======================================"); # log may append runs
}

#------------------------------

=item restrand_seqs

  read output table, collec all newor:- ids, cdna/mrna seqs,
      cdna_bestorf -strand opposite of current -cdna reor.cdna -aa out.aa -cds out.cds -mrna out.mrna
        .. how to flag 'opposite of current' ? 
        if input is oriented mrna, -strand reverse
        if input is orig cdna/tr, -strand is -cds strand info (rnaor:+/- in table?)
        
  # see asmrna_trimvec and reuse methods
  my($ntrim, $trimids, @trimfiles)  ## $outf, $outaa, $outcds
      = mrna_trimvec($mrnaseq,$vecscreenf,$genenames);

	my $trimflag='trimvec_done'; 
  ($upstatus, $upfiles, $uptemp, $upokids)  
    = update_mrna_fileset($trpath, $mrnaseq, $trimflag, $trimids, @trimfiles) if(@trimfiles>0); 

  cdna_evigenesub.pm: 
  update_mrna_fileset($trpath, $inmrna, $updatename, $trimids, @trimfiles)
    my($trimmrna, $trimaa, $trimcds, $trimidfile)= @trimfiles; 

  update expects ($trimmrna, $trimaa, $trimcds, $trimidfile)= @trimfiles;  
  and $trimids->{$id} == $action, action not used,ie DROP or OK ignored, all %trimids drawn from trimfiles.seq
  ie. leave DROPs out of trimfiles.seq set if not wanted, or use drop/cull in restrand_pubids ..
  
    $newortable->{$td}=join"\t", $react,"reor=$or",$upaq,$renote,"oid=$od"; # like uvcut.ids: trid, action, info ..

vectrimset/evg4567hetfix.okay.uvcut.ids 
Ixos9fEVm174352t1	OKCUT,	uvcut=Moderate,33,1-33,end5;	aalen=100,25%,complete-utrbad; clen=1184; strand=+; offs=772-1074; organism=Ixodes_scapularis; evgclass=noclass; oid=Ixos4hEVm088514t1,itick1soapk25loc13015t1
Ixos9fEVm210440t2	OKCUT,	uvcut=end3trim,58,0-0,end5;	aalen=81,57%,partial5; clen=431; strand=+; offs=3-248; organism=Ixodes_scapularis; evgclass=althim; oid=Ixos4hEVm032577t1,itick1soapk25loc26037t1utrorf

test ref
sra2genes_tr19human/tr2aacds_test1908f/try1912v4h

pt=human18ncx; 
env DROPSHORT_ANTISENSE=0.01 $evigene/scripts/genes/trclass_resolve_strandmix.pl \
  -debug -changeonly -trclass $pt.trclass -out $pt.restrand_pubidn
  
# trclass_resolve_strandmix( -debug -changeonly -trclass human18ncx.trclass -out human18ncx.restrand_pubidn ), v=2020.01.12
#egr: evigeneapp=cdna_bestorf.pl, path=/Users/gilbertd/Desktop/genowork/genes2/evigene/scripts/genes/../cdna_bestorf.pl
#readTrClass(human18ncx.trclass)= 82479
#egr: getAaQual: naa=152763 in inputset/human18ncx.aa.qual, val1 humang01912t3= 858,66,3,complete,Code/0.2,55-2631:+
#resolveStrandMix( okayset/human18ncx.pubids.old, human18ncx.trclass)
#resolveStrandMix() nin=82479, nout=2077, ngdiff=390, ngene=23282
#egr: CMD= /Users/gilbertd/Desktop/genowork/genes2/evigene/scripts/genes/../cdna_bestorf.pl  -nostop -minaa=30 -noutrorf -codepot -ostrand=fwd -cdna=human18ncx.reor.fwd.tr -outaa=human18ncx.reor.fwd.aa -outcds=human18ncx.reor.fwd.cds
#egr: CMD= /Users/gilbertd/Desktop/genowork/genes2/evigene/scripts/genes/../cdna_bestorf.pl  -nostop -minaa=30 -noutrorf -codepot -ostrand=rev -cdna=human18ncx.reor.rev.tr -outaa=human18ncx.reor.rev.aa -outcds=human18ncx.reor.rev.cds

** all IDutrorf seqs are missing, due to missing from okayset1st/*ok.tr .. FIX? skip?
?? maybe drop all utrorf called as revaa
human18ncx.restrand_pubidn
 cat human18ncx.reor.ids | cut -f2 | sort | uniq -c
 117 OKorient_f2r   # 71 utrorf missing .. 
 498 OKorient_r2f   # all output to reor.fwd.aa
 144 SKIPasis_f2u
  56 SKIPasis_r2r
  51 SKIPasis_r2u

grep -c '>' human18ncx.reor.*.aa 
  human18ncx.reor.fwd.aa:498
  human18ncx.reor.rev.aa:46  # 71 utrorf missing .. thats okay, problem recalc bestorf-strand for utrorfs
  
=cut


sub restrand_seqs {
  my( $trname, $pubids, $neworpubids)= @_;
  my $cdnaseq=""; # input opt? want okayset/name.ok*.tr
  # my($ntrim, $trimids, @trimfiles)= restrand_seqs($output, $pubids, $trname);

  my $outnam= "$trname.reor";     # for reorset/
  my $outoknam= "$trname.okreor"; # for okayset/
  
  # make reor.ids table as per uvcut.ids, and use to generate new reor.{aa,cds,mrna}
  # print newortable to  "$trname.reor.ids", as per uvcut.ids
  # reortab{id} = OKREV  reor=xxx,notes  from reor_pubids
  #  acts: OKREV  DROPREV (ie tiny new aa) PROBLEMREV ?
  
  # my(@fwds,@revs); # id list of cds direction for bestorf
  my($nfwd,$nrev,%fwds,%revs); # id list of cds direction for bestorf
  open(IDTAB,'>',"$outnam.ids");
  for my $td (sort keys %{$newortable}) {
    my $ntab= $newortable->{$td} or next;
    print IDTAB "$td\t$ntab\n";
    ## ntab ==  $react,"reor=$or/$oldor","aalen=$upaq",$renote,"oid=$od"; # like uvcut.ids

    my @tv= split"\t", $ntab;
    my($react,$reor)= @tv[0,1];  # OKf2r OKr2f .. "reor=$or/$oldor",
    my($oid)= ($ntab=~/oid=([^;,\s]+)/) ? $1 : "noid";
    # my($onew,$oold)= split "/", $reor;
    
    #  my $orchange=  $co . '2' . $cn; f2r/r2f/f2f/r2r/ ? any u2r r2u?
    #  my $react= ($upaq =~ /revaa/) ? "orient_$orchange" : "asis_$orchange";
    #  $react=  ( ($react =~ /asis/) ? 'SKIP' : ($upcla =~ /^cull|^drop/) ? 'DROP' : 'OK') . $react;
    #  # orient_f2r orient_r2f orient_f2f orient_r2r .. orient_x2u ? orient_u2x?
    ## act MAYBE = ambiguous, keep both orients, new id 'xxxt1utrrev' for 2nd
    
    if($react =~ /SKIP/){ }
    elsif($react =~ /OK|DROP|MAYBE/) { #  and $react =~ /orient/
      if($react =~ m/r$/){ $revs{$td}= $oid; $nrev++; } #? make these hashes? $rev{pd}=oid; not $revs{$oid}=$td; 
      elsif($react =~ /f$/){ $fwds{$td}= $oid; $nfwd++; } # NOT $fwds{$oid}=$td; 
      # if($react =~ m/r$/){ push @revs, $td; push @revs, $oid if($oid); } #? make these hashes? $rev{pd}=oid
      # elsif($react =~ /f$/){ push @fwds, $td;  push @fwds, $oid if($oid); }
    }
  } close(IDTAB);

  my $TRSUFFIX='tr|cdna|fasta|fna'; # tr2aacd4 forces okayset1st/.tr suffix; NOT .mrna|cds|aa,  
  my($oktr,$alttr,$okd)= getOkFileset( $okaysetd,$TRSUFFIX,undef,$trname); # okayset1st or okayset
  my @origcds= grep/\.cds/, @$okd; 
  
  #? dump also input cds/mixed gene, cmp orig.cds x new.cds for full/substring identicals
  # use pubids to keep genes grouped
  # ugh, add pubids to these so don't get mixed up? or is rev.ids tab enough info?
  # UPD20apr: ** many reor.aa too short for cdna_bestorf, but have reor.tr;  
  # in bestorf_stranded() rewrite allout.tr to have only aa.ids ?
  # UPD20apr: f/r seqs now (mrna,aa,cds); was (tr,aa,cds)
  my($nokf,$fwdseqs)= ($nfwd) ? bestorf_stranded( 'fwd', $outnam, [$oktr,$alttr], \%fwds) :(0);
  my($nokr,$revseqs)= ($nrev) ? bestorf_stranded( 'rev', $outnam, [$oktr,$alttr], \%revs) :(0);

  my($nrecpot,$recodepotab)= (0,{});
  if($USE_CODEPOTAB) {
    my @recds=(); 
    push @recds, $$fwdseqs[2] if($nokf); # grep /.cds/, @$fwdseqs;
    push @recds, $$revseqs[2] if($nokr); 
    my($cpotabf);
    ($nrecpot,$recodepotab,$cpotabf)= calcCodepotOfCDS($trname, @recds) if(@recds);
    push @tmpfiles,$cpotabf if($cpotabf);
  }
  
  #UPD20ap: dup cds filter : large num of reor seqs are duplicates
  # reor_dupcdsfilt() ; dupcds_filter()
  if($nokf+$nokr>0) {
    my($nuni,$nundup,$ndup,$ndupfwd,$nduprev)= 
      dupcds_filter($outnam, \@origcds,\%fwds,$fwdseqs,\%revs,$revseqs) ;
    $nokf -= $ndupfwd;  $nokr -= $nduprev;
  }
  
  my ($nfwdf, $fwdisfragof)=(0);
  my ($nrevf, $revisfragof)=(0);
  if($DO_COMPARE2ORIG) {

    ## FIXME: okayset vs pubids.old, okayset1st not good choice: pubids in okayset are changed from pubids.old
    # my @origcds= getOkFileset('okayset','cds',undef,$trname); # only 1 now? altcds?
    
    ## .. use same seq set as oktr 
    #above# my @origcds= grep/\.cds/, @$okd; 
  
    ($nfwdf, $fwdisfragof)= ($nokf) ? compare_orig2new_cds( 'fwd', $outnam, \@origcds, \%fwds, $fwdseqs) : (0);
    ($nrevf, $revisfragof)= ($nokr) ? compare_orig2new_cds( 'rev', $outnam, \@origcds, \%revs, $revseqs) : (0);
  }

  #?? Filehandle STDOUT reopened as $hin only for input at 
  # evigene/scripts/genes/../cdna_evigenesub.pm line 83.
  # == sub openRead(), called from where? faidlist()? something(undef)?

  # also update newortable IDTAB file w/ new aaquals, fragof info
  my $nupd= $nfwd + $nrev;  
  my $doupdatetab= $nokf+$nokr; # if $nokf | $nokr | $nfwdf | $nrevf ..
  my %reoridtab=(); # == subset of newortable with changes
  if($doupdatetab) {
    my %ups= ( 
      'fwd' => [$nokf, $fwdseqs, $nfwdf, $fwdisfragof, \%fwds],
      'rev' => [$nokr, $revseqs, $nrevf, $revisfragof, \%revs],
      );
    for my $fr (qw(fwd rev)) {
      my($nok,$seqs,$nfrag,$fragof,$ids)= @{ $ups{$fr} };
      next if($nok<1);
      my ($qcds)= grep /\.cds/, @$seqs;
      my ($qhdr,$nq)= ($qcds and -f $qcds) ? faheaders($qcds) : ({},0);
      for my $od (sort keys %$qhdr) { 
        my $hdr= $qhdr->{$od};
        my $pd= $pubid{$od}||$od; # keep global list of oid => pubid
        # aalen=63,41%,complete; clen=468; strand=-; offs=313-122; codepot=Code/0.015; 
        my($aaw,$offs,$or,$cpot)= map{ my($v)= $hdr=~m/$_=([^;\s]+)/?$1:0; $v; } qw(aalen offs strand codepot);
        my $isfrag= ($nfrag==0) ? 0 : $fragof->{$od} || $fragof->{$pd} || 0;
        ##upd  ntab ==  $react,"reor=$or/$oldor","aalen=$upaq",$renote,"oid=$od";  
        my $ntab= $newortable->{$pd} or next;
        $ntab =~ s/aalen=/naalen=$aaw,$offs:$or\toaalen=/;
        
        if($nrecpot) { $cpot= $recodepotab->{$od} || $recodepotab->{$pd} || $cpot; }
        $cpot =~ s,^[NCU]\w+/,,; # val only
        $ntab =~ s/(renote=\S+)/$1,ncpot:$cpot/; # add new cpot, old there as ,cpot:nnn

        ## INSERT here new comparison of orig and restrand seqs, see reverse_gene_trset
        if(1) {
        # my $isCpot2small = ($cpmax - $cprevmax < 0.001)?1:0; # diff too small to use ?
        my $pRevaa2small= $ENV{pRevaa2small} || 0.60; # cancel or keep both if revaa size is too small
        my ($origaaw)=  $aaqualhash->{$od} || $aaqualhash->{$pd} || 0;
        my ($origcpot)= $codepotab->{$od} || $codepotab->{$pd} || 0; #?? no zero vals in these cpotab?
        my $badreor=0;
        my($naaw,$norigaaw)= map{ my($n)= m/^(\d+)/; $n; } ($aaw,$origaaw);
        
        if($origcpot - $cpot > 0.01) { $badreor |= 1; } # ambiguous .. FIXME better test cpot/origcpot ?
        if($naaw < $pRevaa2small * $norigaaw) { 
          $badreor |=2 
            unless($badreor==0 and $cpot - $origcpot > 0.01 and $naaw >= 0.25*$norigaaw); # qualify rev2small w/ >cpot
          } 
        # UPD20apr: RNAIsStranded: dont change ntab
        $badreor=0 if($RNAIsStranded and $ntab =~ m/r2f/);  
        if($badreor) {
           #o: $ntab =~ s/^(OK|DROP)/MAYBE/; # need other action: BOTH? # OKorient_ > SKIPbadreor_ ??
           $ntab =~ s/^(OK)/MAYBE/; # need other action: BOTH? # OKorient_ > SKIPbadreor_ ??
           $ntab =~ s/(renote=\S+)/$1,reorbad:$badreor/;  
           }
        }
        
        if($NOMAYBE) { # other option: write maybe reor seqs to separate file for further tests
          # $ntab =~ s/^MAYBE/DROP/; #?? DROP the reor .. this way also drops orig.
          $ntab =~ s/^MAYBE/SKIP/; #? this way, FIXED: reor seqset.cds,tr now both have reor seq, want none
        }
        
        if($isfrag){ 
          ## FIXME: some of isfrag are paralogs, use pubid gene id to separate, eg. 1200 of 32000 for ixotick 
          # ?? use maplocs, if avail, to confirm diff loci .. not here?
          (my $gid=$pd)=~s/t\d+$//; 
          my @pmd= map{ $pubid{$_}||$_; } split",",$isfrag; 
          my @frthis= grep /^$gid/, @pmd;  my @frpara= grep{ not m/^$gid/ } @pmd; 
          my $isthis= @frthis > 0;  my $fragof= join",", @frthis, @frpara;
          my $tag=($isthis)?"fragof": "fragpar";
          $ntab =~ s/$/\t$tag=$fragof/; #? change fragof= tag unless(isthis) ? fragpar= ?
          $ntab =~ s/^OK/DROP/ if($isthis);
          #   $isfrag=join",", map{ $pubid{$_}||$_ } split",",$isfrag;  ## change isfrag id list to pubids?  
          #   $ntab =~ s/$/\tfragof=$isfrag/; 
          #   $ntab =~ s/^OK/DROP/;
          }
          
        $newortable->{$pd}= $ntab;
        $reoridtab{$pd}= $reoridtab{$od}= $ntab; # both pubid, oid?
      }
    }
    open(IDTAB,'>',"$outnam.ids.upd");
    for my $td (sort keys %{$newortable}) {
      my $ntab= $newortable->{$td} or next;
      print IDTAB "$td\t$ntab\n";
    } close(IDTAB);
    rename("$outnam.ids","$outnam.ids.old"); 
    rename("$outnam.ids.upd","$outnam.ids"); 
    push @tmpfiles,"$outnam.ids","$outnam.ids.old";
  } 
     
  # also merge/replace old okayset seqs w/ new .. in caller?
  
  ## fixme outfiles: update_mrna_fileset($trpath, $inmrna, $updatename, $trimids, @trimfiles)
  ## expects  ($trimmrna, $trimaa, $trimcds, $trimidfile)= @trimfiles;  
  ## FIXME : MAYBE act nees new ID in seq set , IDrevorf only for .aa/.cds, to keep old + new
  
  if($nupd) {
    my @suf= qw( mrna aa cds ); #o qw( tr aa cds ); # @seqs now has mrna,aa,cds; was tr,,
    my @outfiles= map{ "$outoknam.$_"; } @suf;
    my %ups= ( 
      'fwd' => [$nokf, $fwdseqs],
      'rev' => [$nokr, $revseqs],
      );
    for my $fr (qw(fwd rev)) {
      my($nok,$seqs)= @{ $ups{$fr} };
      if($nok) { for my $sf (@suf) {
        my($fi)= grep /\.$sf/, @$seqs; # @seqs now has mrna,aa,cds; was tr,aa,cds
        my($fo)= grep /\.$sf/, @outfiles;
        # read seqs, use reoridtab{id} action to drop/keep seq?
        my($nokw,$ndrop)= writeReorSeqset(\%reoridtab,$fi,$fo,$sf,$fr); 
        push @tmpfiles,$fi;
        }
      }
    }
    system("cp -p $outnam.ids $outoknam.idtab"); # rename or copy? leave outname.ids in reorset/ ?
    push @outfiles, "$outoknam.idtab"; # last
    return ($nupd, \%reoridtab, \@outfiles); # = $nupd, \%updidtab, [$fixrna, $fixaa, $fixcds, $fixidfile]
  }
  
  return (0);
}

=item writereor

>atap16:AT1G14450.1revorf : revorf has reor.aa,reor.cds but not reor.tr
  -- from bestorf.pl
  type=protein; aalen=92,46%,partial5-utrbad; clen=601; strand=+; offs=2-280; codepot=Code/0.0053; 
  -- from orig.tr
  | NADH dehydrogenase (ubiquinone)s | Chr1:4946192-4947616 REVERSE LENGTH=601 | 201606 
  -- from tr2cds
  evgclass=althi1,okay,match:atap16:AT1G14450.2,pct:100/98/-sense; idis=Arath4bEVm025191t3; 
  -- from reor
  reor=orient_r2f,renote=rnaor:-,newor:-,cpot:+0.0044,ncpot:+0.0070;

add reor=  hdr to all seq? drop cpot:,ncpot .. leave to tables

confusing info newor:- when changed to +, cut or change?
  reor=orient_r2f,rnaor:-,newor:-
  
=cut 

sub writeReorSeqset {
  my($reoridtab,$fin,$fout,$suf,$orfr)= @_;
  my($inok,$inh)= openRead($fin); # die "ERR: reading $pubids" unless($ok);
  # $suf now mrna,aa,cds; was tr,aa,cds
  
  #UPD20may ????: mrna needs revcomp() to match reor cds,aa .. test test test
  # .. BUT should not be wrong, as bestorf_stranded() makes new .mrna which is supposed to be correct, output here
  # .. also fix header offs=end-start > offs=start-end ?
  
  open(my $outh,">>",$fout);
  my($ok,$nok,$ndrop)=(0) x 9;
  # $orfr == 'fwd' or 'rev' current/changed orient w/ respect to rna. change renote newor: to that?
  # my $newor= substr($orfr,0,1); $newor=($newor eq 'f')?'+':($newor eq 'r')?'-':'.';
  while(<$inh>) { 
    if(/^>(\S+)/) { my $od=$1; $ok=1;
      if(my $odor= $reoridtab->{$od}) { 
        my($react,$reval,$naaw,$oaaw,$renote,$oid)=split"\t",$odor;
        if($react =~ /^DROP/){ 
          $ok= 0; # drop .reor AND drop .orig via rewriteTrClass()
        } elsif($react =~ /^MAYBE/) {
          $ok= ($suf =~ m/(cds|aa)/) ? 1 : 0; # dont write 2nd mrna/tr seq as same tr for both fwd/rev orf
          my $newid= $od . 'revorf'; # change id : need record? $newids{$od}= $newid;
          s/>$od/>$newid/; 
          
        } elsif($react =~ /^SKIP/) {
          $ok= 0;  # drop .reor BUT keep .orig via rewriteTrClass()
        
        } elsif($react =~ /^OK/) {
          $ok= 1; # change orig to reor
        }
        
        $react =~ s/^[A-Z]+//; s/$/ reor=$react;/; #? add note to all?
        # $renote =~ s/,cpot:.*//; $renote =~ s/renote://;
        # $renote =~ s/newor:./newor:$newor/; #??
        # s/$/ reor=$react,$renote;/; #? drop renote, react has it
      }
    if($ok){ $nok++; } else { $ndrop++; }
    }
    print $outh $_ if($ok);
  }
  close($outh); close($inh);
  return($nok,$ndrop);
}


=item dupcds_filter

  testix3okallreor_dups.info
  ticks/itick19m/pevg4567itickhetf/okayset1st/
  add to trclass_resolve_strandmix.pl
  
  okayset1st uniq ids w/o reor = 1363523 evg4567hetfix.aa.qual
  reorset ids w/ dups = 97436   evg4567hetfix.okreor.tids
  reorset dup ids (in okay &/or dup)= 66845 evg4567hetfix.okreor.dupids
    ^^ id dups, not seq dups
  
  # reduce identical cds (none in okayset source, all from reorset)
  gunzip -c evg4567hetfix.ok{ay,alt}.cds.gz > testix3okallreor.cds
  
  # reorset w/ uniq ids: reor_f2r/r2f + may(be)
  gunzip -c evg4567hetfix.okreor.cds.gz | perl -ne 'if(/^>(\S+)/){ $id=$ido=$1; ($may)= ($id=~s/(revorf)//)?"may":""; ($ro)=m/reor=orient_(\w+)/?$1:"u2u"; $idn=$id."reor_$ro$may";  s/>$ido/>$idn/; } print; ' >> testix3okallreor.cds
  
  fastanrdb -i -f testix3okallreor.cds | grep '^>'  > testix3okallreor.nrcds.hdr
  
  # pick best id of dups: not_reor > reor > reor_maybe , but ID_reor can replace ID_orig
  
  perl -ne 'next unless(/^>/); s/>//; @d=split; next if(@d<2);  @rd=sort grep/reor_/,@d; @md=grep /may$/,@rd; if(@md){ @rd=grep{ not /may$/ } @rd; } @nd=sort grep{ not /reor_/ } @d;  if(@nd) { $id=shift @nd; $dd=join",",@nd,@rd,@md; } else { @rd=(@rd,@md); ($id)=shift @rd; $dd=join",",@rd; } print "$id dup=$dd\n"; ' \
    testix3okallreor.nrcds.hdr > testix3okallreor.nrcds.bestids

=cut

sub dupcds_filter {
  my($outnam, $origcdsf,$fwds,$fwdseqs,$revs,$revseqs)= @_;
  my $alloutcds="$outnam.allokreo.cds";
  my $APPfastanrdb= findapp("fastanrdb");
  my($ok,$inh);
  # loggit 0, "dupcds_filter...\n";
  my @incds=(); #debug list
  open(OC,">$alloutcds");   
  push @erasefiles,"$alloutcds"; 
  for my $ocf (@$origcdsf) {
    ($ok,$inh)= openRead($ocf); next unless($ok); 
    push @incds,$ocf;
    while(<$inh>){ print OC $_; } close($inh);
  }
  
  # my($fcds)= $$fwdseqs[2]; # "$allout.tr","$allout.aa","$allout.cds"
  my($fcds)= grep /\.cds/, @$fwdseqs; 
  my($rcds)= grep /\.cds/, @$revseqs; 
  my(%roids); 
  my %fwc = ( 'r2f' => $fcds, 'f2r' => $rcds );
  for my $fw (qw(r2f f2r)) { # for my $cf ($fcds,$rcds) 
    my $cf= $fwc{$fw} or next;
    ($ok,$inh)= openRead($cf); next unless($ok);  push @incds,$cf;
    while(<$inh>){ 
      if(/^>(\S+)/){
        my $id=$1; my $ido=$id; my($may)= ($id=~s/(revorf)//)?"may":""; 
        my($ro)=$fw; # m/reor=orient_(\w+)/?$1:"u2u"; #<< all u2u, no reor= tag on these data from bestorf
        my $idn= $id."reor_$ro$may";  
        s/>$ido/>$idn/;  $roids{$idn}= $ido;
      }
      print OC $_; 
    } close($inh);
  }
  close(OC);

  # my $cmd="$APPfastanrdb $cdsseq > $cdsnrseq";
  my(%dups,%undup,%uniq); # dont need
  my($nuni,$nundup,$ndup,$ndupfwd,$nduprev)=(0) x 9;
  my $cmd="$APPfastanrdb -i -f $alloutcds |";  
  loggit 0, "dupcds_filter CMD=$cmd\n";
  loggit LOG_DEBUG,"dupcds_filter $alloutcds = @incds";
  $ok= open($inh, $cmd); unless($ok) { }
  while(<$inh>) {
    next unless(/^>/); s/>//; my @d=split; 
    if(@d < 2){ $nuni++; next; } # dont need:  $uniq{$d[0]}=1;
    my($id,@dups);
    my @nd=sort grep{ not m/reor_/ } @d;  # @nd *should be* uniq already
    my @rd=sort grep /reor_/,@d; 
    my @md=grep /may$/,@rd; if(@md){ @rd=grep{ not m/may$/ } @rd; } 
    if(@nd) { ($id)= shift @nd; @dups=(@nd,@rd,@md); } 
    else { @dups=(@rd,@md); ($id)=shift @dups; } 
    
    # FIXME: these IDs are not same as input, fwds,revs, newortab..
    my $ido= $roids{$id} || $id; 
    $undup{$id}= join",",@dups; $nundup++; # print "$id dup=$dd\n"; 
    # here? remove @dd dups from %$fwds,%$revs
    for my $dd (@dups) {
      my $ddo= $roids{$dd} || $dd; next if($dd eq $id); #??next if($ddo eq $ido);
      my $pd= $pubid{$ddo}||$ddo; # keep global list of oid => pubid
      #? $ddo =~ s/revorf//; #?? this too? what is id mixup here?
      my $fw=($fwds->{$pd})?'f':($revs->{$pd})?'r':'u';
      $dups{$dd}= "$id\t$fw"; #? $dups{$ddo}= "$ido\t$fw"; 
      $ndup++;
      if(delete $fwds->{$pd}) { $ndupfwd++; }
      if(delete $revs->{$pd}) { $nduprev++; }
      if(my $ntab= $newortable->{$pd}) {
        ##upd  ntab ==  $react,"reor=$or/$oldor","aalen=$upaq",$renote,"oid=$od";  
        #x $ntab =~ s/^(OK|MAYBE)/DROP/; # ? DROPdup .. should be SKIP not DROP ? DROP means drop orig ? 
        $ntab =~ s/^OK/DROP/; # ? DROPdup .. OK becomes DROP : reor is valid but duplicate of other, drop
        $ntab =~ s/^MAYBE/SKIP/; # ? DROPdup .. MAYBE becomes SKIP : reor ambiguous, keep orig
        $ntab =~ s/$/\tdupof=$id/;  
        # $ntab =~ s/(renote=\S+)/$1,dupof:$id/; # which annot? see fragof=id
        $newortable->{$pd}= $ntab;
      }
    }
  } close($inh); 

  if($ndup>0) {
    my $otab="$outnam.perfdup.ids"; push @tmpfiles, $otab;
    open(FT,'>',$otab); for my $id (sort keys %dups) { print FT "$id\t",$dups{$id},"\n"; } close(FT);
  }

  loggit 0, "dupcds_filter nuni=$nuni, nundup=$nundup, ndup=$ndup, ndupfwd=$ndupfwd, nduprev=$nduprev\n";
  
  #?? here: rewrite reor seq sets, removing ndups
  #?? should also update newortable w/ DROPdup info, above does rewrite IDTAB 
  # .. maybe dont need rewrite seqs, w/ newortable update
  
  # unless($noDUPSEQUPD) 
  # if(0) { # test : dont need this seq update, done in restrand_seqs
  #   my %ups= ( 
  #     'fwd' => [ $ndupfwd, $fwds, $fwdseqs ],
  #     'rev' => [ $nduprev, $revs, $revseqs ],
  #     );
  #   my $oldf= ($debug) ? \@tmpfiles : \@erasefiles; 
  #   for my $fr (qw(fwd rev)) {
  #     my($nok,$idset,$seqs)= @{ $ups{$fr} }; next unless($nok>0);
  #     for my $cf (@$seqs) {
  #       ($ok,$inh)= openRead($cf); next unless($ok);
  #       my $updf= "$cf.upd"; open(my $outh, ">$updf"); my ($nin,$nup)=(0,0);
  #       while(<$inh>){ 
  #         if(/^>(\S+)/){ 
  #           my $od=$1;  
  #           if(1) { my $isdup= $dups{$od}||0; $ok= ! $isdup; }
  #           else { 
  #             my $pd=$pubids{$od}||$od;
  #             $ok= $idset->{$pd} || 0; #?? miss-id bugs? use instead $notok= $dupset->{$od} ?
  #           }
  #           $nin++; $nup++ if($ok); }
  #         print $outh $_ if($ok);
  #       } close($inh); close($outh);
  #       if($nin and $nup){ # bug: nup == 0 due to bad ids?
  #         rename($cf,"$cf.old"); push @$oldf,"$cf.old";
  #         rename($updf,$cf);
  #       }
  #     }
  #   }
  # }  

  return($nuni,$nundup,$ndup,$ndupfwd,$nduprev); #? ,$fwds,$revs dont need to ret \%fwds,\%revs
}

=item compare_orig2new_cds

  my @cmpstats= compare_orig2new_cds( 'rev', $outnam, \@origcds, \%revs, $revseqs) if($nokr);

  need simple test of reor.fwd/rev.cds that are same-strand dups of other valid gene set cds,
  this may work:

  grep 'newor:-' iscaplacz.restrand_pubido | cut -f3 | sort -u > iscaplacz.reor.gids
  
  faextract( okayset/iscaplacz.okay.cds.gz, iscaplacz.reor.orig.cds, iscaplacz.reor.gids)
  
  makeblastdb -dbtype nucl -in iscaplacz.reor.orig.cds
  for pt in qw(fwd rev) do
    blastn -perc_identity 99.9 -qcov_hsp_perc 99  -strand plus  -outfmt 6 \
      -db iscaplacz.reor.orig.cds -query iscaplacz.reor.$pt.cds -out iscaplacz.reor.$pt.origperffrag.blastn
  
    cut -f1 iscaplacz.reor.$pt.origperffrag.blastn | sort -u > iscaplacz.reor.$pt.origperffrag.tids
  
  iscaplacz.reor.orig.cds is subset of okay.cds with reor.gids, from same loci as reor.ids
    9136 iscaplacz.reor.gids, of 55285 okay locus tot : loci w/ newor:-
  136057 iscaplacz.reor.orig.cds for reor.gids, of 290922 okay.cds, ie nearly 1/2 for 1/6 of loci
  
  cut -f1 iscaplacz.reor.fwd.origperffrag.blastn | sort -u | wc -l
      3875 of 11699 fwd appear to be perffrags, check ids
  cut -f1 iscaplacz.reor.rev.origperffrag.blastn | sort -u | wc -l
      5598 of 15259 rev perffrag of orig
      
  grep -c '>' iscaplacz.reor.*.cds
   iscaplacz.reor.orig.cds:136057
   iscaplacz.reor.fwd.cds:11699
   iscaplacz.reor.rev.cds:15259

=cut

sub compare_orig2new_cds {
  my($fwdrev, $outnam, $origcdsf, $idset, $seqset)= @_; 

  $APPblastn=      findapp("blastn",1) unless($APPblastn); 
  $APPmakeblastdb= findapp("makeblastdb",1) unless($APPmakeblastdb);
  if($APPblastn =~ /MISSING/){ loggit LOG_WARN, "compare_orig2new_cds: $APPblastn\n"; return(0); }
  
  my %gids= map{ (my $g=$_) =~ s/t\d+$//; $g => 1; } sort keys %$idset;
  my $allout="$outnam.orig";
  unlink("$allout.cds") if(-f "$allout.cds");  # remove now, recall here.
  my($nxtot)=(0);
  for my $trf (@$origcdsf) {
    next unless($trf and -f $trf); 
    (my $outna=$trf) =~ s/\.\w+$/.orig/;
    my $newtr= "$outna$$.cds";  
    my($tidh)= faidlist($trf, undef, 'ashash');
    # ugh pubids vs oids, may be either in trf, use $pd = $pubid{$od}||$od
    my %idh; 
    for my $od (keys %$tidh) { my $pd= $pubid{$od}||$od; (my $g=$pd) =~ s/t\d+$//; 
      $idh{$od}=$pd if($gids{$g}); }  
    
    my($newtro,$nx)= faextract($trf,$newtr,\%idh,0,'addidval'); #?? idset must be \%hash
    if(-s $newtr) { 
      system("touch $allout.cds"); system("cat $newtr >> $allout.cds"); 
      unlink($newtr); $nxtot+=$nx; }
  }
  return(0) unless($nxtot>0);
  
  my %qisfragof=();
  my ($qcds)= grep /\.cds/, @$seqset; # .gz ?
  if( $qcds and -f $qcds) { 
    my ($cmd,$cerr);
    $cmd="$APPmakeblastdb -dbtype nucl -in $allout.cds";
    $cerr= runcmd($cmd); if($cerr) { loggit LOG_WARN, "$cerr for $cmd\n"; return(0); }
    my $blopt="-perc_identity 99.9 -qcov_hsp_perc 99  -strand plus  -outfmt 6 ";
    if($qcds =~ /\.gz/) {
      $cmd="gunzip -c $qcds | $APPblastn $blopt -db $allout.cds ";
    } else {
      $cmd="$APPblastn $blopt -db $allout.cds -query $qcds ";
    }
    loggit 0, "fragcheck.$fwdrev: CMD=$cmd \n" ; # if($debug);
    open(P,"$cmd | ") or die "err: $cmd";
    while(<P>){ 
      my($td,$rd)=split; next if($td eq $rd);    
      # ($td,$rd)= map{ $pubid{$_} || $_  } ($td,$rd); # want pubids here?
      $qisfragof{$td} .= "$rd,"; 
      } close(P);
    for my $suf (qw(nsq nin nhr)) { unlink("$allout.cds.$suf"); } # blast idx ; or push @erasefiles;
    #no: push @tmpfiles,"$allout.cds"; # 
    }
  unlink("$allout.cds");  # remove now, recall here.
  
  my $nq=scalar keys %qisfragof; #? write to file
  if($nq>0) {
    my $otab="$allout.perffrag.$fwdrev.ids";
    push @tmpfiles, $otab;
    open(FT,'>',$otab); for my $id (sort keys %qisfragof) { print FT "$id\t",$qisfragof{$id},"\n"; } close(FT);
    }
  loggit 0, "fragcheck.$fwdrev: n=$nq \n"; #  if($debug);
  return($nq, \%qisfragof);
}

sub bestorf_stranded {
  my($forrev, $trout, $trset, $idset)= @_;
  
  my(%idh);
  # idset is now hash{pubid} = oid, add rev idh{oid}=pubid;
  for my $pd (keys %$idset){ my $oid=$idset->{$pd}; $idh{$pd}=$oid; $idh{$oid}=$pd; }
  
  # if(ref($idset) =~ /ARRAY/){ %idh= map{ $_ => 1} @$idset; } # now faextract() does this
  # elsif(ref($idset) =~ /HASH/){ %idh= %$idset; }
  # else { warn"ERR: bestorf_stranded idset not list of ids\n"; return 0; }  
  
  # AAMIN now: my $MINAA = $ENV{aamin} || $ENV{MINAA} || 30; # see at top; tr2aacds4 sets ENV{aamin}, asmdupfilter uses
  # UPD20apr: ** many reor.aa too short, need to capture these, not output by cdna_bestorf w/o -minaa=1
  # ? handle minaa here?
  
  my $noutrorf=1; #?
  my $cmdopt  = " -nostop -minaa=$AAMIN";
     $cmdopt .= " -noutrorf" if($noutrorf); # .=" -noturorf" if($noutrorf); ##201609 NO NOT AT END
     $cmdopt .= " -codepot" if(UPD1908); # -codepotential new option, always? 

  #o my $allout="$trout.$forrev";
  my $intr="$trout.$forrev-in.tr";
  my($outaa,$outcds,$outmrna)= map{ "$trout.$forrev.$_" } qw(aa cds mrna);
  my($nxtot)=(0);
  for my $trf (@$trset) {
    next unless($trf and -f $trf); # [oktr, alttr] may have only 1st
    (my $outna=$trf) =~ s/\.\w+$/.$forrev/;
    my $newtr= "$outna$$.tr"; ## trf == okayset1st/human18ncx.okay.tr.gz ..
    my($newtro,$nx)= faextract($trf,$newtr,\%idh,0,'addidval'); #?? idset must be \%hash
    if(-s $newtr) { 
      system("touch $intr"); system("cat $newtr >> $intr"); 
      unlink($newtr); $nxtot+=$nx; }
  }
  return(0) unless($nxtot>0);
  push @erasefiles,$intr; # @tmpfiles,
  
  # upd: -outmrna=$allout.mrna and replace allout.tr w/ this, solves 2 hassles.
  my $runerr= runcmd("$APPcdnabest $cmdopt -ostrand=$forrev -cdna=$intr -outaa=$outaa -outcds=$outcds -outmrna=$outmrna"); 
  #o my $runerr= runcmd("$APPcdnabest $cmdopt -ostrand=$forrev -cdna=$allout.tr -outaa=$allout.aa -outcds=$allout.cds"); 
  my($ntrin,$ntrout,$ntrdrop,$naa)=(0) x 9;
  $ntrin= $nxtot;
  unless($runerr) { 
    # UPD20apr: ** many reor.aa too short for cdna_bestorf, but have reor.tr; need to handle $AAMIN here?  
    # in bestorf_stranded() rewrite allout.tr to have only aa.ids ?
    # Also, return trids for reor.aa too short, flag these as DROP
    my($aaids)= faidlist($outaa, undef, 'ashash');
    $naa= scalar(keys %$aaids);
    # if($naa < 1){ return(0); }
    if($naa < $nxtot) {
      $nxtot= $naa; # if($nx ne $naa) warn ..
   
      #o my($newtro,$nx)= faextract($intr,"$intr.tmp",$aaids); #?? idset must be \%hash
      #o $nxtot= $nx; # if($nx ne $naa) warn ..
      #o rename($newtro,$intr);
      
      ## also reclass no-aa,tooshort as DROP.tr, see also dupcds_filter
      for my $pd (keys %$idset) { 
        my $oid=$idset->{$pd}; # $idh{$pd}=$oid; $idh{$oid}=$pd; 
        my $ok= $aaids->{$oid};
        if(not $ok and my $ntab= $newortable->{$pd}) {
          $ntab =~ s/^OK/DROP/; # OK becomes DROP : reor is valid but aa too short
          $ntab =~ s/^MAYBE/SKIP/; # MAYBE becomes SKIP : reor ambiguous, keep orig
          $ntab =~ s/(renote=\S+)/$1,aatiny:1/; # aamiss? aanone?
          $newortable->{$pd}=$ntab;  
          delete $idset->{$pd}; $ntrdrop++; # remove from fwd/rev id hash
        }        
      }
       
    }
  }
  
  # update wants: $trimmrna, $trimaa, $trimcds, $trimidfile
  return($nxtot, [$outmrna,$outaa,$outcds]); # $ntrdrop,$naa
  #o return($nxtot, ["$allout.tr","$allout.aa","$allout.cds"]); # $ntrdrop,$naa
}


sub resolveStrandMix {
  my($outh, $pubids, $trinfo, $aaqual, $mapinfo )= @_;  

  %trscore=%trv=%trline= (); # globals for putgene()
  my($nin,$nout,$ndrop,$nkeep,$nrenum,$ngdiff,$ngene,$ncols,
    $ncdsrev, $lgd )=(0) x 10;
    
  my($ok,$pubidh)= openRead($pubids); die "ERR: reading $pubids" unless($ok);
  while(<$pubidh>) {

    unless(/^\w/) {  # BUT print header , adding new col names ??
      if(/^#Pub/) { 
        s/$/\tClass\tAAqual\tpIdAln\tNotes/ unless(/\tClass/); 
        $ncols=scalar( split"\t", $_ );
        print $outh $_; 
        }
      next; } 
      
    chomp; 
    my @v= split"\t"; $nin++; # ,$aq,$pia,$cla
    my $nc= @v; $ncols=$nc if($nc>$ncols); 
    my($pd,$od,$gd,$ti)= @v;
    my($aaq,$pia,$cla,$aaref, $chrmap, $clain, $aarefin)= (0) x 9; # trinfo
    my @xcol=();
    if(@v>4) { 
      ($cla,$aaq,$pia,$aaref)= @v[4..7]; # extended pubids, init?
      ($aaref,$chrmap)= clean_aaref($od,$aaref);
      ($clain, $aarefin)= ($cla,$aaref);
      @xcol= @v[ 8..($nc-1) ] if($nc>8);
    } 
    
    unless($gd eq $lgd) { 
      if($lgd) {
        my($nreor,$nsame)= ($ncdsrev) ? restrand_gene($ncdsrev) : (0,1);
        my($aout,$gdiff,$anum)= putgene($outh, $ncdsrev, $nreor); 
        $ngene++; $ngdiff++ if($gdiff); $nrenum+=$anum; $nout+=$aout; 
        }
      $ncdsrev=0; %trscore=%trv=%trline= (); # globals for putgene()
    }

    $pubid{$od}= $pd; # keep global list of oid => pubid
    # $oidof{$pd}= $od; # and rev?
    
    my $trscore= 999999 - $ti; # altnum for (re)ordering gene alts; largest 1st, keep same order as input here
    my $trval= $trinfo->{$od} || $trinfo->{$pd}||"0"; # ==  join"\t",$aaq,$pia,$cl,$aaref; # KEEP 4 cols
    
    my @tv= split"\t", $trval; 
    if(@tv<4) {  # error!!
      $trval= join"\t",$aaq,$pia,$cla,$aaref if($cla);
    } else { 
      # FIXME: pubids col order differs from trline, change trline
      ($aaq,$pia,$cla,$aaref)= @tv; # FIXME: diff class, maybe diff aaref; NOTE diff col order from pubids, fix?
      my $fixval= (@tv > 4 and $cla)?1:0;
      if($clain and $TRUSTPUBIDS) {
        $fixval ||= (($cla ne $clain) or ($aaref ne $aarefin))?1:0;
        $cla= $clain; 
        if($aarefin =~ /aaref:[1-9]/) { $aaref= $aarefin; } #? and aaq
      }
      $trval= join"\t",$aaq,$pia,$cla,$aaref if($fixval); # error, remake $trval?
    }

    my($aw,$ap,$ac)=split",",$aaq;
    my($pi,$pa)=split"/",$pia; 
    my($cdsantisense)= ($pia=~m,\-sense,)?-1:0; # UPD1912: add explicit revaa flag for pubset annot uses
    
    use constant{ kAATINY =>1, kAAUTRBAD =>2, kAAUTRPOOR =>4, kAADUP =>8, kAAGAPS =>16, kNONCODE =>32, }; 
    my($pflag)= ($trval=~m/pflag:(\w+)/)? $1 :0;  # pflag from, see asmrna_dupfilter3c.pl
    my($aacons)= ($trval=~m/aacons:(\w+)/)?$1:0; 
    my($vaaref)= ($aaref=~m/aaref:(\d+)/)?$1:0; # dont have proportion, use raw/bitscore as weight?
    my($idbest)=   $trval =~ m/selfbest:([^,;\s]+)/?$1:0; 
    # NOTE: idbest is oid not pubid .. resolve
      
    if($cdsantisense < 0) {
      # record it, dont resolve till have all altpars of gene; both fwd/rev may be here
      $ncdsrev++;
    } 

    ## change here to PAL_ANTIMIN test
    my $diffcds= ($pa <= $PAL_ANTIMIN)?1:0;
    my $samecds= ($pi >= $SAME_PI and $pa >= $SAME_PA)?1:0;
    # my $diffcds= ($pa<=$DIFF_PA)?1:0; # skip ident score here? alt is diff enough when align score low
    #   $samecds=0 if($diffcds); # shouldnt need..
  
    use constant kAAComplete => 3;     
    my $icla= ($ac =~ /complete/)?kAAComplete: ($ac =~ /partial[35]/)? kAAComplete-1: kAAComplete-2; 
    $icla-- if($ac =~ /gapbad/);
    $icla-- if($cla =~ /frag/); ## check new class tags, 'nc' noncode, maybe? ..

    my $aaqual=  $aaqual->{$od} || $aaqual->{$pd} || $aaq || "";
    # warn "#D aaqual $pd,$od $aaqual\n" if($debug and $nin < 9);
      
    my ($aasize)= ($aaqual =~ m/^(\d+)/)?$1:($aaq =~ m/^(\d+)/)?$1:0; #?? aaqual bad?
    if($aasize>0 and $aasize<$aw) { $aw=$aasize;  } # fixme: update trinfo/trline ??

    my($mapqual,$antisense,$msplit,$mcov,$mexon,$mapor)= (0,0,0,99,2,0); # defaults no mapinfo val
    if(ref $mapinfo) {
		  ## mapqual == "cov=$cov,nexon=$nx,splicex=$nsx,sense=$asense" : keep any more than antisense?
		  ## upd1805: format now as chrmap:99a,99i,1758l,9x,chr23:42516712-42528502:+
		  my $mfo= $mapinfo->{$od} || $mapinfo->{$pd} || ""; # antisense info; which id here? pd or new realt pd or od ?
		  if($mfo) { 
		  $mapqual=$mfo; # was =1; 
		  $mapor    = $mfo =~ m/:([\+\-\.])/?$1:0; # chr-orient
		  $antisense=($mfo =~ /antisense|nonsense/)?-1:0; # nonsense = uncertain antisense
		  $msplit=   ($mfo =~ /split=(\w+)/)?$1:$msplit; 
	    $mcov=     ($mfo =~ /cov=(\d+)/)?$1:$mcov; # NOPATH: paths=-1 or cov=0
		  $mexon=    ($mfo =~ /nexon=(\w+)/)?$1:$mexon; 
		  
		  $mapor = ($antisense) ? '-' : ($mapor eq ".") ? '0' : '+'; # or w/ respect to cds-orient, not chr-orient
      $trval .= ",chrmap:$mfo" unless($trval=~m/chrmap:/); # new chrmap:99a,99i,1758l,9x,chr23:42516712-42528502:+
      $trval .= ",mapor:$mapor"; #?t +/-/0 sense w/respect to cds
		  }
	  }

    my($codepot,$cpval,$cdsor,$aqval)=(0) x 9; # FIXME get codepot number
    if($USE_CODEPOTAB) {
      ($cpval)= $codepotab->{$od} || $codepotab->{$pd} || 0; #?? no zero vals in these cpotab?
      $codepot= ($cpval >= HCoding) ? 'Code' : ($cpval <= HNoncode) ? 'Noncode' : 'Unknown'
          unless($cpval == 0);
      }
    unless($USE_CODEPOTAB and $cpval != 0) { # CODEPOT1607 required here
      ## FIXME aa.qual table now cuts off "/cpval" from Code/cpval .. want to keep it
      if($aaqual=~m/,(Code|Noncode|Unknown)/) { $codepot=$1;
        ($cpval)= ($aaqual =~ m/,(Code|Noncode|Unknown).([^,;\s]+)/)?$2:0;
        }
      }
    if(1) {
      if($pflag & kNONCODE) { $codepot ||= "Noncode"; } #?
      if($aaqual =~ m/,\d+\-\d+:(.)/) { $cdsor=$1; } # cdsor == ,256-12915:+  ,2491-176:-
      my @aq= split",",$aaqual; $aqval=$aq[2] if(@aq>3); # in -3..3 range ? near equiv of $icla ??
      #? add $aqval to trval
      $trval .= ",codepot:$cpval"; # . substr($codepot,0,1);
      $trval .= ",rnaor:$cdsor" if($cdsor); #? debug?
      }

    # $trval= join"\t",$aaq,$pia,$cla,$aaref; # error, remake $trval?
    my $trline=join("\t",$pd,$od,$gd,$ti,$trval); # SHOULD be 8 cols
    
    $trline.= "\t".join("\t",@xcol) if(@xcol>0); #upd1806
    $trline{$pd}= $trline;
    $trscore{$pd}= $trscore; # for reordering gene alts .. keep same as input here
    
    $aaref =~ s/(aacons|pflag).*//; $aaref ||= 0;

    ## %trv for restrand_gene(); oldor == cdsor
    $trv{$pd}= join "\t", $idbest,$cdsantisense,$cdsor,$cpval,$diffcds,$aaqual,$aaref,$aacons,$mapqual,$pflag; #?aaref or vaaref
    #altreclass $trv{$pd}= join "\t", $cla,$aw,$icla,$samecds,$diffcds,$antisense,$aaref,$mapqual,$xchainval,$aacons,$nocode,$pflag;  # UPD19
    
    $lgd=$gd;  
  } close($pubidh);
  
  if($lgd) {
    my($nreor,$nsame)= ($ncdsrev) ? restrand_gene($ncdsrev) : (0,1);
    my($aout,$gdiff,$anum)= putgene($outh, $ncdsrev, $nreor); 
    $ngene++; $ngdiff++ if($gdiff); $nrenum+=$anum; $nout+=$aout; 
    }
  
  return($nin,$nout,$nrenum,$ngdiff,$ngene);
}


sub reverse_gene_trset {
  my( $bestfwdid, $ids, $newor, $oldor, $codepot, $mapsense, $cpmax, $cprevmax)= @_;
  
  # note bestfwdid == codepot max, also want to see longest fwd tr
  #  test reclass as cull the tiny revaa alts, vs mostcpfwd and longestfwd ?
  # @$ids are ordered alt ids of one gene
  my($longfwdid) = grep{ $newor->{$_} > 0 } @$ids;
  $longfwdid ||= $bestfwdid;

  # need mrna seq here to recalc cds/aa, get new quals
  # work here: revise data tables, new cds/aa seq from rev(mrna)
  #  accumulate all restrand ids, dont work on each gene set ..
  
  # internal data changes: 
  # ignore trv{pd} = tmp data values
  #  putgene() needs only trline{} hash to regen pubids table, but need seq set, aaqual, other updates
  # $trline{$pd}= $trline= join("\t",$pd,$od,$gd,$ti,$trval); # SHOULD be 8 cols
  # $trval= join"\t",$aaq,$pia,$cla,$aaref;   << change aaq, pia, aaref-maybe, cla?
  #  aaref includes notes to change: aaref*, aacons*, selfbest, codepot*,  pflag?

  # my $fwdtrline= $trline{$bestfwdid}; 
  # my($pd,$od,$gd,$ti,$aq,$pia,$cla,$notes,@xcol)=split"\t",$fwdtrline; #? use for ref?
  my $topval = $trline{$longfwdid} || $trline{$bestfwdid} || 0; 
  my($mpd,$mod,$mgd,$mti,$maq,$mpia,$mcla,$mnotes)= split"\t", $topval;
  my($mainlen)= ($maq =~ m/(\d+)/)?$1:0;
   
  my $isCpot2small = ($cpmax - $cprevmax < 0.001)?1:0; # diff too small to use ?
  my $pRevaa2small= 0.30; # cancel or keep both if revaa size is too small
  
  my( $nrev, $nfwd, $nzero)= (0) x 9;
  
  for my $td (@$ids) {
    my $or= $newor->{$td}||0; # -1, 0, +1
    my $oldor= $oldor->{$td}||0;  # -1, 0, +1
    my $cpval= $codepot->{$td}||0;
    # my $cmapsense= $mapsense->{$td}||0;
    
    my @trline= split"\t",$trline{$td}; 
    my($pd,$od,$gd,$ti,$aq,$pia,$cla,$notes,@xcol)=@trline;
    
    my $cor=($or<0) ? '-' : ($or>0) ? '+' : 0;
    my $upnote= $notes . ",newor:$cor";
    $upnote .= "1" if($pd eq $bestfwdid);

    # use constant{ kAATINY =>1, kAAUTRBAD =>2, kAAUTRPOOR =>4, kAADUP =>8, kAAGAPS =>16, kNONCODE =>32, }; 
    # my($pflag)= ($notes=~m/pflag:(\w+)/)? $1 :0;  # pflag from, see asmrna_dupfilter3c.pl
    # my($aacons)= ($notes=~m/aacons:(\w+)/)?$1:0; 
    my($vaaref)= ($notes=~m/aaref:(\d+)/)?$1:0; # dont have proportion, use raw/bitscore as weight?
    if($vaaref > 0) {
      $upnote =~ s/aaref:/revaaref:/; # what?
    }

    my($alen)= ($aq =~ m/(\d+)/)?$1:0;
    my $upaq= $aq;  $upaq .= ",revaa" if($or<0); #?? ",reversed";
    
    #? this pia/.sense should refer only to cds pair align sense? not cmapsense? which sets -or also
    #? need separate cdsor and newor?
    my $upia= $pia;  
    my $ors=($or < 0)? '-sense' : '.';   
    if(0) { # turn off for now, not sure want this cds-sense changed
    if($upia =~ m,/.sense,) { $upia =~ s,/.sense,/$ors,; } else {  $upia =~ s,/\.,/$ors,; }
    }
    
    my $upcla= $cla;
    my $ismaybe=0; # not antiok
    if($or < 0) { 
      # $antiok=1;
      if($isCpot2small and $cprevmax - $cpval < 0.001) { $ismaybe=1; } # $ismaye? ambiguous ..
      # if($newaasize < $pRevaa2small * $oldaasize) { $ismaybe= 1; } # alen == oldaasize 
      
      $nrev++;
      my $docull= (not $ismaybe and $alen < $DROPSHORT_ANTISENSE * $mainlen);
      if($docull and $upcla !~ /^cull/) { $upcla = "cull".$upcla; }  
    } elsif($or > 0) { $nfwd++; 
    } else { $nzero++; }
    
    # pubid trline updates: aaq:revaa, upia:-sense, upcla: cullor, upnote: codepot,newor
    @trline[4,5,6,7]= ($upaq, $upia, $upcla, $upnote);
    $trline{$td}= join"\t",@trline;
    
    # record updates, changed newor ..    
   
    # my $or= $newor{$td}||0;
    ## ** need cdsor from aaqual, vs newor, to reorient.  newortable only for cdsor ne newor ?
    ## FIXME: meaning of newor, oldor ambiguous ..
    ## newor < 0 means need to reverse it, should have revaa tag, regless of oldor value
    ## newor > 0 means keep as is?, regardless of rnaor/cdsor/oldor value
    ## oldor ne newor isnt meaningful?
    # ie. two actions here:  +to- or -to+ ; f2r or r2f ; 
    # if($oldor ne $or) # just ne is wrong .. need all of newor=-
    # if($or < 0)  #  need all of newor=- and change to newor=+ ??

# eg. revaa from rna+ to rna- due to codepot rna- > cpot rna+
#	humang0100996637t1 	main	399,76%,complete	100/99/-sense	
#   selfbest:humang0100996637t2,codepot:0.033,rnaor:-,newor:+1	
#	humang0100996637t2 althi1	269,66%,complete,revaa	100/99/-sense	
#   selfbest:humang0100996637t1,codepot:0.017,rnaor:+,newor:-	
    
    if($or < 0 or $oldor ne $or) {  #? want or >= 0, ie all; $or < 0 == revaa
        # stil wrong here.. do need oldor cmp newor
      my $co=($oldor<0) ? 'r' : ($oldor>0) ? 'f' : 'u'; # is this confusing?
      my $cn=($or<0)    ? 'f' : ($or>0)    ? 'r' : 'u'; # OPPOSITE of newor ie change - to +
      
      # my $react= $co . '2' . $cn;
      # $react=  ( ($react =~ /u$/) ? 'SKIP' : ($upcla =~ /^cull|^drop/) ? 'DROP' : 'OK') . $react;
      # my $react= ($upaq =~ /revaa/) ? 'r' : 'f';
      
      ## this seems correct
      $co=($oldor<0) ? 'r' : ($oldor>0) ? 'f' : 'u'; #?? always
      if($or < 0) { $cn=($oldor>0)?'r':'f'; }  
      elsif($or > 0) { $cn=($oldor>0)?'f':'r'; }  
      else { $cn='u'; }
      my $orchange=  $co . '2' . $cn;
      my $react= ($upaq =~ /revaa/) ? "orient_$orchange" : "asis_$orchange";
      my $isskip= ($react =~ /asis/ or $cn eq 'u');

      # FIXME: utrorf problem, cant readily recalc revaa bestorf for these, drop all?
      my $isdrop= (($upcla =~ /^cull|^drop/) or ($od =~ /utrorf/));
      
      #UPD20apr: RNAIsStranded ?? change ismaybe to drop or ok, depending on f2r = drop / r2f = ok
      if($ismaybe and $RNAIsStranded) { 
        ## maybe becomes SKIP not drop ? ie dont change f2r if rnastranded and ambiguous change
        if($co eq 'f' and $cn eq 'r') { $isskip=1; } 
        elsif($co eq 'r' and $cn eq 'f') { } # isok=1 ?
        #x if($co eq 'f' and $cn eq 'r') { $isdrop=1; } elsif($co eq 'r' and $cn eq 'f') { } # isok=1 ?
        $ismaybe=0; }
        
      $react=  ( ($isskip) ? 'SKIP' : ($ismaybe)? 'MAYBE' : ($isdrop) ? 'DROP' : 'OK') . $react;
      
      ## act MAYBE = ambiguous, keep both orients, new id 'xxxt1utrrev' for 2nd
      #?? add BOTHorient_f2r action when have ambiguous evidence ? ie codepot weak score
      # orient_f2r orient_r2f orient_f2f orient_r2r .. orient_x2u ? orient_u2x?
      
      my $renote="renote=";
      $renote .= join",", grep /\S/, map{ my($v)= $upnote =~ m/($_:[^,;]+)/; $v||"";  } qw( mapor rnaor newor);
      $renote .= ",cpot:$cpval"; #?? want codepot here? or notes column?
      # my $aaqual=  $aaqual->{$od} || $aaqual->{$pd} || $upaq; # has cds-offset:or
      #  if($aaqual =~ m/,\d+\-\d+:(.)/) { $cdsor=$1; } # rnaor == cdsor == ,256-12915:+  ,2491-176:-
    
      #?? want oid or pubid as id here?   seqset may have either ; td == pubid always?
      $newortable->{$td}=join"\t", $react,"reor=$or/$oldor","aalen=$upaq",$renote,"oid=$od"; # like uvcut.ids: trid, action, info ..
    }
    
  }
  return($nrev, $nfwd + $nzero);
}

sub restrand_gene {
  my($ncdsrev)= @_;
  ## replaces resolveOneAnti()
  
  my @td= sort{ $trscore{$b} <=> $trscore{$a} or $a cmp $b } keys %trscore;
  # my $t1= shift @td or return 0; # main/longest tr : is this changing from orig main?

  ## call each @tr as fwd or rev strand from max cpval case : NOT always accurate, other revaa quals?
  
  my($ntr,$nrev,$nsame,$hasmap,$cpmax,$cpmaxid,$cprevmax,$cprevmaxid)= (0) x 9;
  my(%pair,%newor,%oldor,%codepot,%mapsense);
  
  foreach my $td (@td) {
    $ntr++;
    
    ## %trv for restrand_gene()
    ## $trv{$pd}=join "\t", $idbest,$cdsantisense,$cdsor,$cpval,$diffcds,$aaqual,$aaref,$aacons,$mapqual,$pflag; #?aaref or vaaref
    my($idbest,$cdsantisense,$cdsor,$cpval,$diffcds,$aaqual,$aaref,$aacons,$mapqual,$pflag)= split"\t", $trv{$td};
    $cdsantisense=1 if($cdsantisense == 0);
    $oldor{$td}= ($cdsor eq '-') ? -1 : ($cdsor eq '+') ? 1 : 0;
    
    # NOTE: idbest is oid not pubid .. resolve
    $idbest= $pubid{$idbest} || $idbest; # keep global list of oid => pubid

    # FIXME here? ** NEED to match this test, to keep from re-reversing low align pairs
    # .. cancel/remove pair{td}{idbest} if pct-align is below PAL_ANTIMIN threshold
    #dupfilt4 does this, so trclass has separated low-align antis as separate loci
    # if(ANTI1912c and $isanti) {  $antiflag= "/." if($pal < $PAL_ANTIMIN); } # turn off, not enough align
    $idbest .= ".diffcds" if($diffcds); #?
    
    $pair{$td}{$idbest}= $cdsantisense;
    $codepot{$td}=$cpval;
    if($cpval > $cpmax){ $cpmax=$cpval; $cpmaxid=$td; }
    if($cpval < $cprevmax){ $cprevmax= $cpval; } # min cp now
    
    if($mapqual) {
		  ## mapqual == "cov=$cov,nexon=$nx,splicex=$nsx,sense=$asense" : keep any more than antisense?
		  ## upd1805: format now as chrmap:99a,99i,1758l,9x,chr23:42516712-42528502:+
		  my ($nx)= ($mapqual =~ /,(\d+)x,/)?$1:0;
      my ($mapsense)= ($mapqual =~ /antisense/) ? -1 
        : ($mapqual =~ /nonsense/) ? 0 # uncertain antisense
        : ($nx>1) ? 1 : 0; # add fwdsense , nonsense ?? or -sense/+sense ?
      $mapsense{$td}=$mapsense;  $hasmap++ if($mapsense);
    }
    # $idbest{$idbest}++;  #? $aaqual{$td}="$aaqual,$aaref,$aacons";
  }    
  
  $newor{$cpmaxid}=1; 
  my($maxor)=(1);

  # rule chr-mapsense takes precedence ? but need to also check cdssense of pair aligns,
  # .. reset newor for mixed gene set some w/ chrmap, some without ? *should* have chrmap for all of gene set
  # for poor/none chrmap should ensure cdssense tracks w/ those having chrmap sense
  if($hasmap) { 
    for my $td (@td) { 
      if( my $or= $mapsense{$td}) { $newor{$td}=$or; $maxor= $or if($td eq $cpmaxid); }
    } 
    for my $td (@td) { 
      unless( $mapsense{$td}) { # set newor via cds-pair link to chr-mapsense set?
        for my $bd (sort keys %{$pair{$td}}) { 
          if(my $nor= $newor{$bd}) { my $oor= $pair{$td}{$bd}; my $or= ($oor < 0) ? -$nor : $nor; $newor{$td}=$or; last; }
        }      
      }
    }
  }
  
  #FIXME codepot scores need to be modified by other quals to determine best orient
  # a. skip codepot if rev aa << smaller than long aa
  # b. USE_CODEPOT calc better choice than current bestorf static cpot table
  # c. ambiguous where diff(cpmax - cprev) < tiny val, ignore cp for fwd/rev then?
  
  for(my $it=0; $it<3; $it++) {
    my $ormis=0;
    for my $td (@td) {
      my $or=0; 
      next if($newor{$td});
      #x if($hasmap and $mapsense{$td}) { $or=$mapsense{$td}; }
      if($td eq $cpmaxid) { $or=$maxor;  } 
      elsif(my $por= $pair{$td}{$cpmaxid}) { 
        $or= $maxor * $por; if($por < 0) { my $tcp= $codepot{$td}; if($tcp>$cprevmax){ $cprevmax=$tcp; $cprevmaxid=$td; } }
        } # need these.. else miss many newor
      elsif(my $por= $pair{$cpmaxid}{$td}) { 
        $or= $maxor * $por; if($por < 0) { my $tcp= $codepot{$td}; if($tcp>$cprevmax){ $cprevmax=$tcp; $cprevmaxid=$td; } }
        } 
      else {
        for my $bd (sort keys %{$pair{$td}}) { 
          if(my $nor= $newor{$bd}) { my $oor= $pair{$td}{$bd}; $or= ($oor < 0) ? -$nor : $nor; last; }
        }
      }
      
      if($or) { $newor{$td}=$or; } else { $ormis++; }
    } 
    last if($ormis < 2);
  } #re-iterate find newor ? resolves 1-2%
  
  if($cpmax - $cprevmax < 1e-4) { } # diff too small to use ?
  
  # my @idrestrand= grep{ $newor{$_} < 0 } @td;
  # my @idokstrand= grep{ not ($newor{$_} < 0) } @td; 
  # my @idunk=  grep{ not $newor{$_} } @td; # are any @td not called in newor ?
  # $nrev= @idrestrand; $nsame= @idokstrand;
  
  ($nrev,$nsame)= reverse_gene_trset( $cpmaxid, \@td, \%newor, \%oldor, \%codepot, \%mapsense, $cpmax, $cprevmax);  
  ## some of restranded will have cds == idok set .. maybe can be dropped?

  return($nrev,$nsame);
}


sub putgene {
  my($outh, $ncdsrev, $nreor)=@_; 
  my($ialt, $changed, $renum)=(0) x 9; # 

  return(0,0,0) if($PUTONLYCHANGES and not $ncdsrev); ## TEST DEBUG
  
  my @tr= sort{ $trscore{$b} <=> $trscore{$a} or $a cmp $b } keys %trscore;
  my $t1= shift @tr or return 0; # main/longest tr : is this changing from orig main?

  foreach my $tr ($t1,@tr) {
    ++$ialt;
    if($tr eq $t1) { # check for noclass drop option

    } else {
    
    }
    
    my $trline= $trline{$tr}; 
    my($pd,$od,$gd,$ti,$aq,$pia,$cla,$notes,@xcol)=split"\t",$trline; 
    # pubid trline updates: aaq:revaa, upia:-sense, upcla: cullor, upnote: codepot,newor
  
    print $outh join("\t",$pd,$od,$gd,$ti,$cla,$aq,$pia,$notes||'.',@xcol)."\n"; 
  }
  
  $changed= $nreor; #?  
  return ($ialt,$changed,$renum); # , scalar(@keeps), scalar(@drops)
}

=item calcCodepotOfSeqset

  calc cds codepot using valid aalong mrna (cds - utr) stats
  a. tabulate loglikecoding,  $EVIGENES/prot/cdshexprob.pl -bestaa -mrna $oktr -aaqualdata $cpotdata
  b. tabulate llcode for cds, $EVIGENES/prot/cdshexprob.pl -aaqualdata $cpotdata -out $okcpot -testcds $okcds $altcds
  
=cut


sub calcCodepotOfSeqset { # $USE_CODEPOTAB 
  my( $trname,$pubids,$okaysetd,$pubsetd,$aaqualf)= @_;
  # ($ncpot,$codepotab)= calcCodepotOfSeqset($trname,$pubids,$okaysetd,$pubsetd,$aaqualf); # fill in \%codepotab{id}=llcoding

  my $REUSEcpot=1; #?? option, check file dates
  my($ncpot,%codepotab)= (0);
  my($okcds,$altcds,$okd)= getOkFileset( $okaysetd,'cds',undef,$trname); # okayset1st or okayset
  return (0) unless($okcds);

  $APPcdshexprob= findevigeneapp("prot/cdshexprob.pl") unless($APPcdshexprob);  # capture stderr to log .. runcmd() patch? 
  
  # need to make sample aaqual.data table from best mrna ? look for pubsetd/mrna, okayset/mrna ..?
  # my($cpotdata); # global now
  my($cmd, $err)= (0) x 9;
  unless($cpotdata and -f $cpotdata){
    ( $cpotdata= $aaqualf) =~ s/\.aa.qual$/.aaqual.data/;  # trname.aa.qual
    ( $cpotdata= $okcds) =~ s/\.\w+$/.aaqual.data/ unless(-f $cpotdata);
  }

# BUG: $oktr == array ..
# sh: -c: line 0: `/oasis/projects/nsf/ind114/ux455375/chrs/evigenes/sra2genes_testdrive/bio/apps/evigene/
#    scripts/prot/cdshexprob.pl -bestaa -mrna ARRAY(0x11681a8) -aaqualdata okayset/plYYPE.okay.aaqual.data'
  
  unless(-f $cpotdata) {
    my($okdir,$oktr,$alttr,$okdtr)=(0) x 9;
    
    ## FIXME may not have evg format mrna here, need APPcdshexprob to use name.tr + name.cds hdrs
    ## BUG 20feb: got wrong okayset/name.okay.tr instead of pubset/name.mrna where pubset/name.pubids is ..
    ## .. But, prior app call is trclass2pubset.pl -onlypub .. but now that does make pubset/name.mrna
    ## bug may be getOkFileset() that gets m/.okay|.okalt/ instead of getFileset()
    ## .. look for oktr.cds if oktr.tr lacks CDSoffs ? cdshexprob.pl does that.. but for file suffix patt bug!
    
    #o ($oktr,$alttr,$okdtr)= getOkFileset( $pubsetd,'mrna',undef,$trname); # try for pubset/mrna 1st?
    ($okdir,$oktr,$alttr,$okdtr)= getFileset( $pubsetd,'mrna',undef,$trname); # try for pubset/mrna 1st?
    unless($oktr) {
      ($okdir,$oktr,$alttr,$okdtr)= getFileset( $okaysetd,'mrna|tr',undef,$trname); # okayset1st 
    }
    # $evigene/scripts/prot/cdshexprob.pl -bestaa=300  -mrna okayset1st/evg4567hetfix.okay.tr -aqual okayset1st/evg4567hetfix.aaqual.data 
    loggit 0, "tabulate loglike-coding from CDS - UTR hexamers of mRNA $oktr\n";
    $cmd="$APPcdshexprob -bestaa -mrna $oktr -aaqualdata $cpotdata";  
    $err= runcmd($cmd); return(0) if($err); #runcmd logs errs
    $REUSEcpot=0;
  }
  
  # return calcCodepotOfCDS($trname, $okcds, $altcds);
  (my $okcpot= $okcds) =~ s/\.\w+$/.codepot/;
  $REUSEcpot= 0 unless(-s $okcpot );
  if( -s $okcpot and $REUSEcpot) {
    my $Mcds= -M $okcds;
    my $Mcpot= -M $okcpot;
    $REUSEcpot= 0 if($Mcds < $Mcpot); # time from now ??
  } 
  unless($REUSEcpot) {
  $cmd="$APPcdshexprob  -aaqualdata $cpotdata -out $okcpot -testcds $okcds $altcds";  
  $err= runcmd($cmd); return(0) if($err);
  }
  
  open(CPT,$okcpot);
  while(<CPT>){ if(/^\w/){ my($id,$val)=split; $codepotab{$id}= $val; $ncpot++; } } close(CPT);
  loggit 0, "calcCodepotOfSeqset($okcds $altcds) n= $ncpot\n";
  return ($ncpot,\%codepotab, $cpotdata);
}

sub calcCodepotOfCDS { # $USE_CODEPOTAB 
  my( $trname,$okcds,@altcds)= @_;
  # my( $trname,$pubids,$okaysetd,$pubsetd,$aaqualf)= @_;

  my($ncpot,%codepotab)= (0);
  # my($okcds,$altcds,$okd)= getOkFileset( $okaysetd,'cds',undef,$trname); # okayset1st or okayset
  return (0) unless($okcds);

  $APPcdshexprob= findevigeneapp("prot/cdshexprob.pl") unless($APPcdshexprob);  # capture stderr to log .. runcmd() patch? 

  my($cmd, $err)= (0) x 9;
  (my $okcpot= $okcds) =~ s/\.\w+$/.codepot/;
  $cmd="$APPcdshexprob  -aaqualdata $cpotdata -out $okcpot -testcds $okcds @altcds";  
  $err= runcmd($cmd); return(0) if($err);
  open(CPT,$okcpot);
  while(<CPT>){ if(/^\w/){ my($id,$val)=split; $codepotab{$id}= $val; $ncpot++; } } close(CPT);
  loggit 0, "calcCodepotOfCDS($okcds @altcds) n= $ncpot\n";
  #  push @tmpfiles,$okcpot;
  return ($ncpot,\%codepotab,$okcpot);
}
  
=item readTrClass/readPubidClass

  * update to read same info from extended .pubids
  pubids.hdr.ext   Public_mRNA_ID originalID PublicGeneID AltNum Class  AAqual pIdAln Notes

  my($ok,$pubidh)= openRead($pubids); die "ERR: reading $pubids" unless($ok);
  while(<$pubidh>) {
    next unless(/^\w/); chomp;  
    my($pd,$od,$gd,$ti,$cla,$aq,$pia,$note)=split"\t"; 
  ...
  
  ## FIXME: *** input aa-size in trclass is wrong, no gap adjust; 
  ##     original trclass/pubids sorting is right from using aa.qual table **

  ## FIXME: add aaref = aablastp ref from flag column 6, if there, as new qual score.
  ## asmrna_dupfilter2.pl -ablast blastp.tall4, adds bitscore,refid, and refbest/refok/refgood flag
  ## empty col6 == 0,0,pflag:0  [pflag == poor if > 0]
  ## aaref col6 == 250,arath:AT4G28760.2,refok,pflag:0
  ## 164.4,arath:AT3G42860.1,pflag:0
  ## 224,arath:AT1G69530.4,aadup:1AB-I11_VO_k30Loc10139t3,pflag:0
  
=cut

sub readPubidClass {
  my($pubids)= @_;
  my($ok,$hin)= openRead($pubids);  die "ERR: reading $pubids" unless($ok);
  my(%oda,%aaqual);  my $ntr=0;
  while(<$hin>) {
    next unless(/^\w/); ## and /\t(okay|maybeok)/ ## classes now: main|noclass|alt|cull|part? ???
    chomp; my @v=split"\t"; 
    ##  Public_mRNA_ID originalID PublicGeneID AltNum Class  AAqual pIdAln Notes; Notes=aaref:xxx,chrmap:xxx,pflag:nnn,feq:xxx,...
    my($pd,$od,$cl,$aaq,$pia,$aaref)=@v[0,1,4,5,6,7];  # same as readTrclass, BUT class, aaref may be updated
    
    my($chrmap);
    ($aaref,$chrmap)= clean_aaref($od,$aaref);
    
    # FIXME: pubids col order differs from trline, change trval to (cla,aaq,pia,aaref)
    my $val=join"\t",$aaq,$pia,$cl,$aaref; # add?
    $oda{$pd}= $oda{$od}= $val; $ntr++;
    $aaqual{$pd}= $aaqual{$od}= $aaq; #? skip readAaQual, but that adds gaps not here
    } close($hin);
    
  loggit 0, "readPubidClass($pubids)= $ntr\n";
  return (\%oda,\%aaqual);
}

sub clean_aaref {
  my($td,$aaref)=@_;
  ## check for chrmap:... other in notes?
  my $chrmap="";
  if($aaref=~/(chrmap|mapq):(\w[;\s]+)/) {
    $chrmap=$1;  
    $chrmap =~ s/,(?:pflag|trscore):\d+.*//; 
    $aaref=~s/[,]?$chrmap//;
    # $mapinfo->{$td}=  $chrmap unless($mapinfo->{$td});
  }
  #? require /aaref:[1-9]../
  #? keep/use pflag qual? 
  $aaref =~ s/^0,0[,]?//; 
  $aaref =~ s/(,pflag:\d+).*/$1/; 
  $aaref =~ s/,aadup:[^,;\s]*//; # ,aadup:idnnn,refbest, ; ok here? keep ref(best|good|ok) flag
  $aaref ||= "0";
  #?? $aaref.=",$chrmap" if($chrmap);
  return($aaref,$chrmap);
}

=item rewriteTrClass

  need to add reor changes to trclass for current usage
  ie updates are in okayset/name.okreor.{seqs,idtab} and idtab == reorhash info to change
  need trclass update so that later consumers have new aaqual, sense, etc info
  ugh.

    # ntab ==  $react,"reor=$or/$oldor","aalen=$upaq",$renote,"oid=$od";  
    $ntab =~ s/aalen=/naalen=$aaw,$offs:$or\toaalen=/;
    if($isfrag){ 
      $ntab =~ s/$/\tfragof=$isfrag/; 
      $ntab =~ s/^OK/DROP/;
      }
    $reoridtab{$pd}= $reoridtab{$od}= $ntab; # both pubid, oid?
  Homsap4aEVm003311t12	OKorient_r2f	reor=-1/-1	naalen=157,23%,complete-utrbad,66-539:+	
    oaalen=197,29%,partial3-utrbad,revaa	renote=rnaor:-,newor:-	oid=humang0389170t2

=cut

sub rewriteTrClass {
  my($trclass, $reoridtab)= @_;
  my($nup,$ntr)=(0,0);
  (my $updclass=$trclass) =~ s/\.gz//; $updclass.=".upreor";
  my($ok,$hin)= openRead($trclass);  loggit(LOG_DIE,"reading $trclass") unless($ok);
  open(UPT,'>',$updclass) or loggit(LOG_DIE,"writing $updclass");
  while(<$hin>) {
    unless(/^\w/) {  next; }# no hdr?
    chomp; my @v=split"\t"; 

    ## trclass cols:
    ##  oid,okay/drop,class,idbestmatch,pIdAln,aaqual,aaref/flags

    my($od,$okdrop,$cl,$idbest,$pia,$aaq,$aaref)= @v;  
    if($reoridtab->{$od} or $reoridtab->{$idbest}) {
      my $asense= ($pia=~/(\-sense)/)?$1:""; #? add -sense if od <> best
      my $odor= $reoridtab->{$od}||"";
      my $beor= $reoridtab->{$idbest}||"";
      my $upskip=0;
      if($odor){ 
        my($react,$reval,$naaw,$oaaw,$renote,$oid)=split"\t",$odor;
        my($ren)= $react=~m/^[A-Z]+(\w+)/?$1:$react; 
        if($react =~ /^DROP/){ 
          $okdrop=~s/^okay|maybeok/drop/; 
        } elsif($react =~ /^SKIP/){ 
          $upskip=1; # dont change UPT @v, not $nup
        } elsif($react =~ /^OK/) { 
        
        } elsif($react =~ /^MAYBE/){ 
          # MAYBE ? adds row: leave current as is, add new IDrevorf
          $okdrop=~s/^okay/maybeok/;  
          $v[1]= $okdrop; $v[6]= $aaref . ",maybeor:$ren";
          print UPT join("\t",@v)."\n"; # put orig
          $od .= 'revorf'; # change id  
          if($cl =~ /^main/){ $cl="althi"; } # cant have main w/ newid .. makes new locus; noclass too?
        }
        unless($upskip) {
          if($asense and not $beor) { $pia=~s/$asense/./; }
          ($aaq)= ($naaw =~ m/aalen=([^;]+)/)?$1:$aaq; $aaq=~s/,\d+\-\d+:.//; # drop offs?
          $aaref .= ",reor:$ren"; # OKorient_r2f .. add old.aaw?
          $nup++;
        }
      }
      if($beor) {
        if($asense and not $odor) { $pia=~s/$asense/./; $nup++; }
      }
      @v=($od,$okdrop,$cl,$idbest,$pia,$aaq,$aaref); #? unless($upskip);
    }
    
    print UPT join("\t",@v)."\n";
  } close(UPT);
  loggit 0, "rewriteTrClass $trclass to $updclass,  nup= $nup\n";  
  return($nup,$updclass);
}

=item readTrClass
  
  basic trclass read, 
  recall that fields idbest,pia,aaref/notes are packed w/ several subfields
  
=cut

sub readTrClass {
  my($trclass)= @_;
  my($ok,$hin)= openRead($trclass);  die "ERR: reading $trclass" unless($ok);
  my(%oda,%aaqual,%idbest,%antis); 
  my $ntr=0;
  while(<$hin>) {
    next unless(/^\w/ and /\t(okay|maybeok)/); ## maybeok !!
    chomp; my @v=split"\t"; 
    ## trclass cols:
    ##  oid,okay/drop,class,idbestmatch,pIdAln,aaqual,aaref/flags
    #o my($od,$cl,$pia,$aaq,$aaref)=@v[0,2,4,5,6];  

    #UPD1912: Fixme for cdsantisense need to record idbestmatch, check/reset -sense flags
    #  if(idbest is dropped), drop -sense flag
    #  if(idbest class/size/qual < this class/size/qual), swap -sense flags
    my($od,$okdrop,$cl,$idbest,$pia,$aaq,$aaref)= @v;  
    
    $idbest{$od}{$idbest}=$pia; #?? use this instead of aaref/notes selfbest: ?
    $antis{$od}="$idbest\t$pia" if($pia=~m/\-sense/);
    
    #dupfilt4 does this, so trclass has separated low-align antis as separate loci
    # ** NEED to match this test, to keep from re-reversing low align pairs
    # if(ANTI1912c and $isanti) {  $antiflag= "/." if($pal < $PAL_ANTIMIN); } # turn off, not enough align
    
    my($chrmap);
    ($aaref,$chrmap)= clean_aaref($od,$aaref);
    # my($piav)= $pia =~ m,^(\d+/\d+),;  #? full form: pi/pa[/optbestmatchid]
    # FIXME: pubids col order differs from trline, change trval to (cla,aaq,pia,aaref)
    
    $aaref .= ",selfbest:$idbest" if($idbest); # bestmatch == selfbest ?
    
    # my $val=join"\t",$aaq,$pia,$cl,$aaref,$idbest; #UPD1912 add idbest?? NO, see below KEEP 4 cols; 
    my $val=join"\t",$aaq,$pia,$cl,$aaref; # keep this way
    $oda{$od}= $val; $ntr++;
    $aaqual{$od}=$aaq; #? skip readAaQual, but that adds gaps not here
    } close($hin);
    
  loggit 0, "readTrClass($trclass)= $ntr\n";  
  return (\%oda,\%aaqual);
}




1;

__END__

# evgmrna2tsa2.pl get_evgtrset(); move to package ***
# sub get_evgtrset { 
#   my($trclass,$cdnaseq,$pubdir)= @_;
#   my($trpath,$trname,$nsra,$sradatah)=("","",0,undef);
# 	my $notokay=0;
#   
#   if($cdnaseq) { 
#     $notokay=1; # dont look in okayset/? look in $pubdir now?
#     $trclass= makename($cdnaseq,".trclass") unless($trclass); 
#   }
# 
#   if($trclass) {
#     my $trpname= makename($trclass,"","trclass"); 
#     if($trpname =~ m,/,) { ($trpath,$trname)= $trpname =~ m,(.*)/([^/]+)$,; } # was BAD?
#     else { $trname=$trpname; }
#     $trpath ||= '.';  
#     
#     my $okpath= ($notokay) ? $trpath :"$trpath/okayset"; # should I check if curdir has okayset files?
# 
#     #??? dont need mRNA made new here, want trset that matches aaqual oids, strands
#     
#     ($cdnaseq)= getmRNA($okpath,$trname,$pubdir,LOG_WARN) if(!$cdnaseq and -d $okpath);
#   }
#   
#   return($cdnaseq,$trpath,$trname,$sradatah);
# }


=item readMapsenseTab replaced

## REPLACEd readMapsenseTab with evigene_pubsets.pm readAlignTab()
##    ($nmapqual, $mapqualh, $alntabh)= readAlignTab($mapqual);

sub readMapsenseTab {
  my($maptab)= @_;
  return {} unless($maptab and -f $maptab);
  my($ok,$hin)= openRead($maptab);  
  unless($ok){ warn "missing maptab $maptab"; return {} }
  my %maps=(); my $nmap=0;
  
  ##old# GenomeID gespan geor AQueryID quspan match qlen cov pid path indels nexon splice aalen offs aamap sense tag
  ##updated 2014: add oid alt ID
  ## GenomeID gespan geor AQueryID quspan match qlen cov pid path indels nexon splice aalen offs aamap sense oid tag
  ## opt use slim map.attr
  ## AQueryID cov pid splice/nexon  GenomeID:gespan:geor  path  oid
  ## antisense flag from gmap align tab is useless for other bad-map quals, eg split, : okay below
  
  my (@hd,@hset); 
  my @hdef=qw(GenomeID gespan geor AQueryID quspan match qlen cov pid path indels nexon splice aalen offs aamap sense oid tag);
  my @Idef= (3,6,7,8,9,11,12,13,16); # @alset index to @hdef
  my @alset=qw(QueryID qlen cov pid path nexon splice aalen sense); # align.tab
  my @maset=qw(QueryID cov pid splice path); # map.attr SKIP FOR NOW ?
  if(1) { # default for no header?
    @hset=@alset; @hd=@hdef;
  }
  while(<$hin>) {
    next unless(/^\w/);
    chomp; my @v=split"\t";   
    my($td,$ql,$cov,$npath,$nx,$nspl,$aw,$sens,$oid,
      $chr,$cbe,$cor,$pid) = (0) x 19;
    
    if(/^[A-Z]\w+ID\t/) { # ^GenomeID|AQueryID header
      if(/\tcov/ and /\tsense/ and /\tsplice/) {
        @hset=@alset;
        @hd=@v; map{ s/AQuery/Query/; } @hd; 
        my %hset=map{ $_=>1 } @hset; for my $h (@hd) { $hset{$h}=2 if($hset{$h}); } # check fields
        my @miss= grep{ $hset{$_} == 1 } @hset;
        if(@miss) { die "ERR: missing Mapinfo table fields: @miss\n"; }
        next;
      } elsif(/\tcov/ and /\tsplice/ and /\tpath/) { # skip this one for now?
        @hset=@maset;
        @hd=@v; map{ s/AQuery/Query/; s,splice/nexon,splice,; } @hd; 
        my %hset=map{ $_=>1 } @hset; for my $h (@hd) { $hset{$h}=2 if($hset{$h}); } # check fields
        my @miss= grep{ $hset{$_} == 1 } @hset;
        if(@miss) { die "ERR: missing Mapinfo table fields: @miss\n"; }
        next;
      }
      
    } elsif(@hd) {
      ($chr,$cbe,$cor)= @v[0,1,2];
      my %v=(); for my $i (0..$#v) { my $h=$hd[$i]; $v{$h}=$v[$i]; }
      ($td,$ql,$cov,$pid,$npath,$nx,$nspl,$aw,$sens)= @v{@alset};  
      #ma: ($td,$cov,$pid,$nspx,$npath)= @v{@maset};  
      $oid=$v{oid}||$td;
    } else {
      ($td,$ql,$cov,$pid,$npath,$nx,$nspl,$aw,$sens)=@v[@Idef]; # == [3,6,7,8,9,11,12,13,16]; # or require this?
      $oid=$v[17]||$td;
    }
    
    ## split flag:   add NOPATH flag? or use cov == 0
    my $msplit= ($npath=~m/(C[123]):/)?$1:0;
    my $mpath = ($msplit or $npath eq "0/0")?0:($npath=~m,^(\d+/\d+),)?$1:($cov==0)?-1:0; # NOPATH=nocov=-1 ?
    my $asense=0;
    my $MINCOV_ANTISENSE= 50; # was 85, depends on mapping qual .. new gmap genes.cdsx should be more valid at low align?
    if($sens<0 and not $msplit) {
      ## new nx == 2.2 9.5 1.1 1.0 for nexon / nvalidsplice; new nspl == 0 for all
      my $nsx = 2; ## bad: ($nx<2)?0:$nspl/$nx; #? bad? trust gmap -sense call? $ nspl == 0 now no good
      my $asc=  ($nx > 3 and $nsx > 1.5 and $cov >= 90) ? 3 # certain?
              : ($nsx >= 1.5 and $cov >= 95) ? 2
              : ($nsx >= 1.8 and $nx > 2 and $cov >= $MINCOV_ANTISENSE) ? 1 : 0;
      $asense= ($asc>0)?"antisense":"nonsense"; # change 0 to maybe-antisense ? nonsense?
    }
     
    # if(MAPQ1805) 
    # chrmap: prefix or not? chrmap:99a,99i,1234l,9x,..
    my $quals= sprintf "%da,%di,%dl,%dx,%s:%s:%s", ($cov,$pid,$ql,$nx,$chr,$cbe,$cor);
    $quals.=",$asense" if($asense); # nonsense ? $quals.=",sense=$asense" if($asense);
    $quals.=",split=$msplit" if($msplit);
    $quals.=",paths=$mpath" if($mpath);
    
    $maps{$td}=  $quals; $nmap++;
    if($oid){ for my $od (split",",$oid) { $maps{$od}= $quals; } }
    } close($hin);
    
  warn "#readMapinfo($maptab)= $nmap\n";
  return ($nmap>0)? \%maps : {}; 
}
=cut

=item readAaQual obsolete

# readAaQual: inputset/*.aa.qual, or use evigenesub:fasize_nogaps(publicset/*.aa_pub.fa.gz,okc='A-WYZa-wyz',gapc='Xx*')

sub readAaQual { 
  my($trclass,$trinfo,$aaqualh)= @_;
  
  # fixme: update prior %aaqual from trclass/pubids, dont make new..
  unless(ref $aaqualh){ my %aaqual=(); $aaqualh= \%aaqual; }
  my $naa=0;
  
  my($name,$path,$suffix) = fileparse($trclass,qr/\.\w*/); # suf-'.trclass' or suf = qr/\.\w*/
  my $aaqualf = "$path/okayset/$name.aa.qual"; 
  $aaqualf = "$path/publicset/$name.aa.qual" unless(-f $aaqualf);
  $aaqualf = "$path/inputset/$name.aa.qual" unless(-f $aaqualf);
  $aaqualf = "$path/$name.aa.qual" unless(-f $aaqualf);

  ## this is a mess, should update aa.qual fields in .trclass, others, to include nnnX gaps and reduce? aasize
  if(-f $aaqualf) {

    ($aaqualh,$naa)= getAaQual($aaqualf); # cdna_evigenesub;
    #egr: getAaQual: naa=39851 in .//ref_arath16ap.aa.qual, UPD200114
    #     val1 atap16:AT4G35090.1= 492,61,3,complete,Code/0.0605,290-1768:+
    # aa.qval ==  "$alen,$pctcds,$acv,$aqual1,$codepot?"
    # acv == numeric of aqual1:complete/partial/.. adds gapbad tag for gaps
    
  }  
  
  warn "#readAaQual($aaqualf)= $naa\n"; # if debug
  return $aaqualh; # scalar(%aaqual)? \%aaqual : undef;
}
=cut


=item resolveOneAnti obsolete here

  now in  altreclass.pl : evigene alt-reclassifier
  -- remove there if this is in pipe
  

sub resolveOneAnti { # ($aid,$trinfo,$apubid); UPD1912, maybe do in reclassAlts calcs
  my($aid,$trinfo,$apubid)= @_;
  my($ab,$upd)= (0,0);  
  $apubid ||= $aid;
  
  my $aval= $trinfo->{$aid} || $trinfo->{$apubid} || "0"; # ==  join"\t",$aaq,$pia,$cl,$aaref; # KEEP 4 cols
  ## my $aval= $trinfo->{$aid}||""; ## $aaq,$pia,$cl,$aaref; # add bid == bestmatchid
  my @aval= split"\t", $aval;
  my($aaw,$apia,$acla,$aref)= @aval; ## NO ,$abestid
  
  my ($bid)= $aref=~m/selfbest:([^,;\s]+)/?$1:0; # $abestid;
  # my($bid)= split"\t", $antis->{$aid};

  my $bval= $trinfo->{$bid} || 0;
  my @bval= split"\t", $bval;
  my($baw,$bpia,$bcla,$bref)= @bval; # NO ,$bbestid

  unless($bval) { # missing/dropped
    $ab += 9;
  } else {
    my($aho)= ($aval=~m/aaref:(\w+)/)?$1:0;
    my($bho)= ($bval=~m/aaref:(\w+)/)?$1:0;
    my($acons)= ($aval=~m/aacons:(\w+)/)?$1:0;
    my($bcons)= ($bval=~m/aacons:(\w+)/)?$1:0;
    
    if($aaw > $baw) { $ab += 1; } elsif($aaw < $baw) { $ab -= 1; }
    if($acla =~ /main|noclass/ and $bcla !~ /main|noclass/) { $ab += 2; }
    elsif($acla !~ /main|noclass/ and $bcla =~ /main|noclass/) { $ab -= 2; }
    if($acons > $bcons) { $ab += 4; } elsif($acons < $bcons) { $ab -= 4; }
    if($aho > $bho) { $ab += 8; } elsif($aho < $bho) { $ab -= 8; }
  }
  
  if($ab > 0) {
    if($apia =~ m/\-sense/) { $apia =~ s,\-sense,.,;  $upd++; }
    if($bpia and $bpia !~ /\-sense/) { $bpia .= "/-sense";  $upd++; }
  } elsif($ab < 0) {
    if($bpia =~ /\-sense/) { $bpia =~ s,\-sense,.,; $upd++; }
    if($apia and $apia !~ /\-sense/) {$apia .= "/-sense";   $upd++; }
  }
  if($upd) {
    if($apia){ $aval[1]= $apia; $trinfo->{$aid}= $aval= join"\t",@aval; }
    if($bpia){ $bval[1]= $bpia; $trinfo->{$bid}= $bval= join"\t",@bval; }
  }
  
  return($upd, $aval); # bval?
} 

=cut

