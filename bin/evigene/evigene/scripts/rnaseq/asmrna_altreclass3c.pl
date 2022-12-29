#!/usr/bin/env perl
# altreclass.pl : evigene alt-reclassifier

## FIXME: *** input aa-size in trclass is wrong, no gap adjust; 
##     original trclass/pubids sorting is right from using aa.qual table **

=item about altreclass.pl

 altreclass.pl : evigene alt-reclassifier
 evigene/scripts/rnaseq/asmrna_altreclass.pl 
 	- follows use of prot/tr2aacds.pl, rnaseq/asmrna_dupfilter, evgmrna2tsa.pl
 	  which create inputs of project/.trclass and project/publicset/.pubids
 	  
  examine evigene mRNA asm publicset and mrna.trclass, to
  identify/remove trival alternates (can be large amount, tends to collect as more trasm added)
  -- drop fragment althi subset: shorter aa, partial/fragment flags
  -- option altrenum: renumber alts so t1 main is longest aa (should have been, not always),
      and t2..tn ordered by alt aasize.  public gene ID is preserved.
  -- option maxaltsame: drop excess complete althi that are same aasize as others 
  
=item usage

  evigene/scripts/rnaseq/asmrna_altreclass.pl -maxalt=19 -altrenum  -trclass kfish2evg367mixx.trclass -out 

 	output: publicset/kfish2evg367mixx.realt_pubids 
  #altreclass: nin/out=499047/499047, ndrop/keep=141000/358047, nrenum=297357, ngenediff=35964/149901

  ## -mapsensetab = gff align.tab, isnt required but -sense == mismapping demotes alts as do utrbad,..
  ## many -sense are gene joins w/ utrbad also, -sense from half mapping to rev-gene
  
  daphmag/rnas/asmrna543all/evg7i/   pt=dmagset56i
  $evigene/scripts/rnaseq/asmrna_altreclass.pl -trclass $pt.trclass -mapsensetab gmap/${pt}pub.align.tab \
    -noclass=70 -altrenum -out -debug > & publicset/$pt.realt.log
  #altreclass: nin/out=218653/218653, ndrop/keep=29076/189577, nrenum=160028, ngenediff=14711/38086

=item use test 1805

  $evigene/scripts/rnaseq/asmrna_altreclass.pl \
   -mapqual ../geneval/pig4321ew-chrpig11c.align.tab \
   -trclass ../pig4321ew.trclass -pubids pig4321ew.pubids.old \
   -altrenum -debug -noclasscut=60 \
   -out pig4321ew.pubidnu.realt > & pig4321ew.pubidnu.realt.log

  # adding -mapqual keeps a few more alts from drop class
  # same ids as orig pubid altreclass, mostly same class counts, but
  pig4321ew.pubidnu.realt | cut -f1,5-
  *? move chrmap before pflag?
  #Public_mRNA_ID	Class	AAqual	pIdAln	Notes
  Susscr4EVm000001t1	main	33067,98%,complete	100/99/.	aaref:60808,human18nc:XP_016860310.1,refbest,pflag:0,chrmap:100a,99i,100233l,304x,chr15:84227558-84499478:-,tscore:63471,
  Susscr4EVm000001t2	althi1	32579,98%,complete	100/99/.	aaref:60604,human18nc:XP_016860310.1,refgood,pflag:0,chrmap:100a,99i,98792l,289x,chr15:84227145-84501223:-,tscore:62881,
  Susscr4EVm000001t3	althi1	32530,98%,complete	100/99/.	aaref:60305,human18nc:NP_596869.4,refbest,pflag:0,chrmap:100a,99i,98645l,290x,chr15:84227145-84501223:-,tscore:62683,
  ...
  Susscr4EVm000002t1	main	8855,95%,complete	100/100/.	aaref:16065,human18nc:XP_006715475.1,refgood,pflag:0,chrmap:99a,99i,27831l,146x,chr1:13715658-14201689:+,paths=1/2,tscore:16888,
  Susscr4EVm000003t1	main	7765,95%,complete	100/100/.	aaref:12703,human18nc:XP_005249372.1,refbest,pflag:0,chrmap:100a,99i,24391l,103x,chr7:28812406-29300219:+,tscore:14117,
  Susscr4EVm000004t1	main	7650,97%,complete	100/100/.	aaref:13629,human18nc:XP_005270753.1,refbest,pflag:0,chrmap:100a,99i,23646l,96x,chr6:95125018-95317348:+,tscore:14465,
  Susscr4EVm000005t1	main	7324,98%,complete	100/99/.	aaref:8734,human18nc:NP_001092093.2,refbest,pflag:0,chrmap:100a,100i,22408l,102x,chr2:51323014-51534239:+,tscore:11691,
  --- vs original pubids ---
  #Public_mRNA_ID	Class	AAqual	pIdAln	Notes
  Susscr4EVm000001t1	main	33067,98%,complete	100/99/.	aaref:60808,human18nc:XP_016860310.1,refbest,pflag:0,tscore:63471,
  Susscr4EVm000001t2	althi1	32579,98%,complete	100/99/.	aaref:60604,human18nc:XP_016860310.1,refgood,pflag:0,tscore:62881,
  ..
  Susscr4EVm000002t1	main	8855,95%,complete	100/100/.	aaref:16065,human18nc:XP_006715475.1,refgood,pflag:0,tscore:16888,
  Susscr4EVm000003t1	main	7763,95%,complete,2X	100/100/.	aaref:12703,human18nc:XP_005249372.1,refbest,pflag:0,tscore:14115,

      
=item FIXME replace bad-main w/ good alt, 2014.05

  - now is not changing any main class?  
  - want to replace poor qual but aalonger main w/ good qual alts 
    poor main = utrorf, utrbad/poor, partial vs 
    good alt  = complete, utrgood, aasize >= 80-90%? of main ?

=item todo

	merge with evigene asmrna_dupfilter (main classifier)
	offer option to rewrite publicset/ : only need to remove drop sequences (to 2nd files)
		.. but may need logic from evgmrna2tsa to do in compatible way (several tables to update).
	
	- add option for more agressive junk removal?
	-- drop incomplete/frag alts; drop noclass/incomplete|short<60aa 
		(largest fraction of noclass are minimal aa size, with few bases diff.  This is bad
		consequence of using %align criteria, need to add minimal base criteria to it.
		e.g. 5% diff for 1000bp = 50 significant, 5% diff for 100bp = 5 ns.
		
	cat  $pt.realt_pubids | grep -v drop | env noclassmin=60 perl -ne\
	'BEGIN{ $NOCLMIN=$ENV{noclassmin}; } $ok=(/\tmain/ or /complete/)?1:0; @v=split; ($aw)=$v[5]=~m/(\d+)/;
	 $ok=0 if(/noclass/ and $aw<$NOCLMIN); print if($ok);' > $pt.hivalue_pubids
	
	for tp in "aa cds mrna"; do {
 		gunzip -c $pt.${tp}_pub.fa.gz | cat $pt.hivalue_pubids - | env idpre=Pita perl -ne\
 		'BEGIN{ $IDPRE=$ENV{idpre}; } if(/^$IDPRE/) { ($d)=split; $ok{$d}=1; } else { 
 		if(/^>(\S+)/) { $d=$1; $ok=$ok{$d}; } print if($ok); }' > $pt.hivalue.$tp
  }		
		
=item notes  

  require sorted by gene at least, also by aa-size? or do that here?
  	-- input pubids should be gene-sorted
  * revised to skip pre-pubidx table creation .. input .trclass + publicset/.pubics
  
 grep okay *.trclass | cat - publicset/kfish2evg367mixx.pubids | perl -ne \
'if(/\tokay/) { @v=split"\t"; ($od,$cl,$pia,$aq)=@v[0,2,4,5]; \
($piav)= $pia =~ m,^(\d+/\d+),; $oda{$od}="$aq\t$piav\t$cl"; } \
elsif(/^Fun/) { ($pd,$od,$td,$ti)=split; $oda=$oda{$od}||0; s/$/\t$oda/; print; }' |\
 sort -k3,3 -k5,5nr -k1,1 > publicset/kfish2evg367mixx.pubidx2

 cat publicset/kfish2evg367mixx.pubidx2 | env altrenum=0 altreclass.pl > reclass.tab
 # redo for .trclass .pubids input, skip pubid2x

=cut

use FindBin;
use lib ("$FindBin::Bin/../"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...

use strict;
use Getopt::Long;
use cdna_evigenesub;  
# use cdna_proteins;
use File::Basename qw(basename dirname fileparse);

use constant VERSION =>  '2019.12.18'; # UPD1912
   # '2019.09.25'; # upd from trset alt exon chains, bestfiltcompare, asmrna_dupfilter3c.pl compat,  
use constant UPD1912 => 1;
use constant UPD1908 => 1;
use constant AACONS => 1;  # aaconsensus: keep tr with consensus of identical aa across asm types
use constant CODEPOT1607 => 1; # add from above codepot work

# use constant VERSION => '2018.06.24'; # upd for evg pubset, revised mrna2tsa
# '2016.07.18'; # update to work on pubid table not trclass, input aaref,maploc quals
#'2014.05.14';  # '2013.08.16';  

#UPD19: review/upd these defaults..
#o my($SAME_PI,$SAME_PA,$DIFF_PA)= (98,97,80); # samecds= $pi>=98 and $pa>=97; diffcds= $pa<=80
my($SAME_PI,$SAME_PA,$DIFF_PA)= (100,99,90); # samecds= $pi>=98 and $pa>=97; diffcds= $pa<=80
my($DROPSHORT_AAPART,$DROPSHORT_AAFULL,$DROPSHORT_ANTISENSE, $DROPTINYALT,$DROPSHORT_BADEXONS)
  = (0.70,0.20,0.60,0.20,0.20); # UPD1908
  # = (0.90,0.20,0.60,0.10); # 201702 upd

#UPD1912: PHETERO should adjust identity values down: $SAME_PI .. others?
my $PHETERO = (UPD1912) ? ($ENV{PHETERO}||0) : 0;

my($ALTGOODasMAIN)=(0.90); # aasize percentof for alt to replace main

my $REVISEPUBFILES=0;
my $arenum=$ENV{altrenum}||0; #? default on ?
my $maxsame=$ENV{maxsame}||999; # for alts complete, same size, drop excess

my $noclasscut= $ENV{noclasscut}||$ENV{aaminnoclass}||0; # minaa size for noclass tiny excess,
## caller trclass2pubset.pl uses $ENV{aaminnoclass}||60; for asmrna_altreclass -noclasscut=60
## NONCODE_CUT modify noclasscut when aa is called Noncode
use constant NOCLASS_NONCODE_CUT => 1; # option? need to improve codepot calc otherwise this cuts valid prots
use constant NOCLASS_NONCODE_CUTAA => 100; #was 140; test off, option? fix codepot calc?

my $NODROPS=$ENV{nodrops}||0; # upd1806

## reclassAlts: aa-poor and antisense demote alt scoring
my $pDemoteAaPoor    = $ENV{scalepooraa} ||0.90; #o0.90; # partial[53],utrbad,utrorf; NOT utrpoor
my $pDemoteAntiSense = $ENV{scaleanti} || 0.95; #o0.90; # antisense map, was 0.75
my $pDemoteSplitmap  = $ENV{scalesplit} || 0.80; #o0.90; # split map , other for <MINCOV ? $cov*score ? not good if other qual high
my $pPromoteAaRef    = $ENV{scaleaaref} || 0.50; # def: $pPromoteAaRef=0.50 assumes bitscore ~ 2 x aa
my $pDemoteNoncode   = $ENV{scalenonecode} ||0.80; #? UPD1908
my $pPromoteAaCons   = $ENV{scaleaacons} || 0.10; # UPD1908 :  trscore * $aacons * pPromote, cons.vals=2,3,4.. add 10% per cons
my $pPromoteGoodExonchain= $ENV{scaleexons} || 0.50; # UPD1908, needs test

my $MINCOV= $ENV{mincov}||75; # chrmap coverage min? 90/80/75/50 ?
my $KEEPMAINatTOP= $ENV{topmain}||0; #? test dont demote orig main too much (t2..t9)
my $TRUSTPUBIDS= $ENV{trustpubids}||0; # trust input pubids for class, aaref calls, ie newer than trclass
my $DID_RESOLVESENSE= $ENV{DID_RESOLVESENSE} || 0; # UPD20j from trclass_resolve_strandmix.pl		

my $debug=0;
my ($trclass,$pubids,$output,$valoutput,$mapsensetab,$exontab);

my $optok= GetOptions(
  "class|trclass=s", \$trclass, # require
  "pubids=s", \$pubids, # option/alt to trclass
  "output:s", \$output, # option
  "valueoutput:s", \$valoutput, # option
  "mapsensetab|mapqual=s", \$mapsensetab, # option, or call opt -mapsensetab ? -gmaptable? or gmap.gff instead?
  "exontable=s", \$exontab, # option, or call opt -mapsensetab ? -gmaptable? or gmap.gff instead?
  #"logfile:s", \$logfile, # option
  "maxsame|maxaltsame=s", \$maxsame, 
  "noclasscut|aaminnoclass=s", \$noclasscut, 
  "altrenum!", \$arenum, 
  "trustpubids!", \$TRUSTPUBIDS, 
  "nodrops", \$NODROPS, 
  "revisepublicset!", \$REVISEPUBFILES, 
  "debug!", \$debug, 
  );

die "EvidentialGene altreclass VERSION ",VERSION,"
  examine Evigene mRNA publicset, mrna.trclass to reclassify trival alternates

Usage: altreclass.pl -trclass mymrna.trclass > mymrna.reclass.pubids
opts: -pubids publicset/mymrna.pubids  -out mymrna.reclass.pubids 
      -mapsensetab=gmap_align.table -[no]altrenum  -maxaltsame=$maxsame\n"
  unless($optok and ($trclass or $pubids));  

sub MAIN_begin {}

unless($pubids) { ($pubids=$trclass) =~ s,(\w+).trclass,publicset/$1.pubids,; } # default evg file set
if(defined $output and not $output) {
	($output=$pubids) =~ s/.pubids//; $output.=".realt_pubids";
}

my $outputh=*STDOUT; # $outh
if($output)	{
  rename($output,"$output.old") if(-f $output);
  open($outputh,'>',$output) or die "ERR: writing $output"; 
}

my $outvalsh= undef;
if(defined $valoutput) {
  unless($valoutput) { $valoutput= $output || "altreclass"; $valoutput=~s/\.\w+$//; $valoutput.=".vals"; }
  open($outvalsh,'>',$valoutput) or die "ERR: writing $valoutput"; 
}

# openloggit($logfile,$trclass);
# loggit(1, "EvidentialGene altreclass.pl VERSION",VERSION);
# loggit(1, "ERR: unused arguments:",@ARGV) if(@ARGV>0);
# loggit(0, "altreclass: in trclass=$trclass, pubids=$pubids, out=$output") if($debug);

warn "#altreclass: in trclass=$trclass, pubids=$pubids, out=$output\n";
warn "#altreclass: opts pDemoteAaPoor=$pDemoteAaPoor, pDemoteAntiSense=$pDemoteAntiSense, pDemoteSplitmap=$pDemoteSplitmap,\n";
warn "#  pPromoteAaRef=$pPromoteAaRef, pPromoteAaCons=$pPromoteAaCons, pDemoteNoncode=$pDemoteNoncode, mincov=$MINCOV, altrenum=$arenum, noclasscut=$noclasscut\n";

my($lgd,%trscore,%trv,%trline,%keepdrop,%newid); # globals now for reclassAlts/putgene
my ($trinfo,$aaqualhash);

#?? read both? maybe publicset/pubids takes precedence now?
if($trclass) { ($trinfo,$aaqualhash)= readTrClass($trclass); }
elsif($pubids) { ($trinfo,$aaqualhash)= readPubidClass($pubids); } # always read if TRUSTPUBIDS ?

#o# my $aaqh= readAaQual($trclass||$pubids, $trinfo); 
#o# $aaqualhash= $aaqh if(ref $aaqh and scalar(%$aaqh));
($aaqualhash)= readAaQual($trclass||$pubids, $trinfo, $aaqualhash); 

my $mapinfo= readMapsenseTab($mapsensetab); # update to add more gmap quals
my $exoninfo= readExonTab($exontab); 

my ($nin,$nout,$ndrop,$nkeep,$nrenum,$ngdiff,$ngene)
  = reclassAlts($outputh, $pubids,$trinfo,$aaqualhash,$mapinfo,$exoninfo);  

close($outputh);

## need logic from evgmrna2tsa make_pubseq(), make_annotab() to do in compatible way (several tables to update) ..
# revisePublicset() if($REVISEPUBFILES);

warn "#altreclass: nin/out=$nin/$nout, ndrop/keep=$ndrop/$nkeep, nrenum=$nrenum, ngenediff=$ngdiff/$ngene\n";

sub MAIN_END {}

#---------------------

=item readMapsenseTab/readGmapQual
  * update to read gmap.attr (brief) and/or gmap.align (wide) tables

	aligntab3.sh  mrna.gmap.gff > mrna.alnsense.tab
	OR add here code of aligntab3.sh ..

kfish2evg367mixx_realt.alnsense.tab
GenomeID        gespan  geor    AQueryID        quspan  match   qlen    cov     pid     path    indels  nexon   splice  aalen   offs    aamap   sense   tag
Scaffold0       5067-10366      -       Funhe2Exy3m032549t1     1-679   678     679     100.0   99.9    0       0       3       6       214,94%,complete        19-663  220     0       gmap
Scaffold0       25351-26047     .       Funhe2Exy3m069279t1     1-697   697     697     100.0   100.0   0       0       1       0       115,49%,complete-utrpoor        295-642 115     0       gmap
Scaffold0       30698-31095     -       Funhe2Exy3m078047t1     1-374   360     374     100.0   95.5    0       1/3     2       0       106,85%,partial3        55-372  45      -1      gmap
     ^^ asense but no valid splice 
Scaffold0       131984-136801   -       Funhe2Exy3m021876t4     1-1261  1254    1261    100.0   99.2    0       0/3     9       18      279,66%,complete        422-1261        270     -1      gmap
     ^^ antisense.valid.ids: hicov, full valid splice/9exons


evg4corn_tgok2oid.map.attr
AQueryID	cov	pid	splice/nexon	GenomeID:gespan:geor	path	oid
Zeamay4gEVm000001t1	100.0	99.9	64/65.62	NC_024464.1:154250538-154294086:+	0	Zeamay4EVm000001t1
Zeamay4gEVm000001t2	100.0	99.3	45/45.45	NC_024464.1:154268669-154294090:+	0	Zeamay4EVm000001t2
Zeamay4gEVm000001t3	52.0	97.9	20/21.18	NC_024464.1:154250975-154268670:+	0	Zeamay4EVm000001t3
Zeamay4gEVm000002t1	57.6	96.9	9/12.8	NC_024465.1:70655982-70669439:+	0	Zeamay4EVm000002t1

------> evg4corn_tgok2pub.align.tab.gz <------
GenomeID	gespan	geor	AQueryID	quspan	match	qlen	cov	pid	path	indels	nexon	splice	aalenoffs	aamap	sense	oid	tag
NC_001666.2	18-1190	.	Zeamay4gEVm016889t1	1-1173	1169	1173	100.0	99.7	0	0	1.1	0	353,90%,complete	41-1102	353	0	Zeamay4EVm014354t1,cornlo2m9slvelvk125Loc2249t1	gmap
NC_001666.2	316-2184	.	Zeamay4gEVm016889t3	1-2007	1731	3298	60.9	91.2	0	144/6	1.1	230,21%,complete-utrbad	1587-2279	278	0	Zeamay4EVm014354t5,cornhi8ms9sgvelvk37Loc3472t1	gmap

=item mapqual_brief

  2018.05 update, new short mapqual statement as per
    evigene_pubsets.pm: sub gene_annot_brief()
    
  annotbrief is brief string of gene attributes,
  used in various evigene tables (eg  trclass2mainalt.pl):
    chrmap:alignpct[a],identpct[i],mrnalength[l],nexons[x],chrlocation
  Reproduces this annotation of trclass2mainalt.pl
  aaref:5767,dapsim:Dapsim1EVm000004t1,chrmap:100a,98i,25555l,33x,scaffold29:324762-356329:+ 
  aaref:496,zfish16nc:NP_938183.2,chrmap:99a,99i,1758l,9x,chr23:42516712-42528502:+
    
=cut

sub readExonTab {
	my($maptab)= @_;
	return {} unless($maptab);
	my($ok,$hin)= openRead($maptab);  die "ERR: reading $maptab" unless($ok);
	my %maps=(); my $nmap=0;

  # exontab from commonexons.pl
  # "\t",$pd,$td,$xn,$icv,$commval,$icstat,$gn
  # pubid,oid,xnidlist, ichain:l14.t1:_ebE___lIhKuFg_Drxxxx	20c:......_lIhKuFg_Drxxx	a42s41u1;sum:a85s71u14/a1606s1592u14;imax:42;icalt:37,icsub:42;iruns:
  
	while(<$hin>) {
		next unless(/^\w/);
		chomp; my @v=split"\t";		
    my($pd,$oid,$xn,$icv,$commval,$icstat,$gn)= @v;

    my($icid)= $icv=~m/^(\w+:[\w\.]+)/?$1:$icv; # split locus parts l14 t1 ?
    $commval =~ s/:.*//; $commval="0".$commval unless($commval =~ /^\d/);
    $icstat =~ s/;.*//;
    #** $gn may be pubid_iclocus for altpars; need to keep, check, change publid locus info??
    ## icstat == "a42s41u1" .. a42c20s22 .. may be useful exon type stat
    ## icid   == "l14.t1" .. "l14.t12" .. useful to separate pars by ic locus?
    my $quals="$commval,$icid,$icstat"; #? use only commval (\d+), also commval == nocomm|noalts

		$maps{$pd}=  $quals; $nmap++;
		$maps{$oid}= $quals if($oid and $pd ne $oid);
		} close($hin);
		
	warn "#readExontab($maptab)= $nmap\n";
	return ($nmap>0)? \%maps : {}; 
}

=item exontab annot eg human

   try9bestevgc/pub3test/trev19human.pubids.realt3h 
   
   * note aacons overweighed aasize for t1,t2 here .. but longest 3 w/o aacons are best ref, now t3,4,5
   Hstra9Evm000001t1 exonch:100%cs,100%cx,9c,icalt:l66.t4,a32c1s31,aacons:3,  main	7425aa
   Hstra9Evm000001t2 exonch:100%cs,85%cx,9c,icalt:l66.t18,a27c1s26,aacons:2   alt	  7326aa
   Hstra9Evm000001t3 exonch:100%cs,100%cx,9c,icalt:l66.t14,a29c1s28,  althim 7478aa bref:humang023499t9
   Hstra9Evm000001t4 exonch:100%cs,100%cx,9c,icalt:l66.t11,a30c1s29     alt  7476aa bref:humang023499t12
   Hstra9Evm000001t5 exonch:100%cs,100%cx,9c,icalt:l66.t12,a30c1s29     alt  7471aa bref:humang023499t16
   Hstra9Evm000001t7 exonch:100%cs,100%cx,9c,icalt:l66.t10,a30c1s29     alt	 7462aa bref:humang023499t19 
   Hstra9Evm000001t20 exonch:89%cs,96%cx,8p,icalt:l66.t2,a33c1s30u2     alt	 7432aa
   Hstra9Evm000001t29 exonch:67%cs,70%cx,6p,icalt:l66.t23,a24c1s23      alt	 5358aa bref:humang023499t31
   Hstra9Evm000001t37 exonch:0%cs,44%cx,nocomm,cx,icsub:l66.t1,a14c1s13,aacons:2  alt	3793aa
   Hstra9Evm000001t50 exonch:22%cs,7%cx,2p,icalt:l66.t47,a3s2u1   cullalt  196aa
   
   Hstra9Evm000011t1	main 4678aa aacons:3,exonch:100%cs,100%cx,20c,icalt:l109.t5,a25c8s17
   Hstra9Evm000011t2	alt	4730aa  aacons:2,exonch:100%cs,100%cx,20p,icalt:l109.t4,a26c8s18  bref:humang023077t29
   Hstra9Evm000011t3	alt	4727aa  aacons:2,exonch:90%cs,100%cx,18p,icalt:l109.t8,a24c8s16  bref:humang023077t28
   Hstra9Evm000011t4	alt	4675aa  aacons:2,exonch:90%cs,100%cx,18p,icalt:l109.t10,a23c8s15 bref:humang023077t26
   Hstra9Evm000011t5	althim 4755aa  exonch:100%cs,100%cx,20c,icdup:l109.t5,a25c8s17  bref:humang023077t9
   Hstra9Evm000011t13	alt	3713aa   exonch:70%cs,71%cx,14p,icalt:l109.t12,a19c8s9u2 bref:humang023077t15
   Hstra9Evm000011t20	alt	2656aa	 exonch:60%cs,67%cx,12p,icsub:l109.t2,a16c8s8
   Hstra9Evm000011t36	cullalt	68aa	aacons:2,pflag:33,exonch:0%cs,0%cx,noalts,cx,ichain:l170718.t1,a2u2

=cut
    

use constant MAPQ1805 => 1;		

## REPLACE readMapsenseTab with evigene_pubsets.pm readAlignTab()
##    ($nmapqual, $mapqualh, $alntabh)= readAlignTab($mapqual);

sub readMapsenseTab {
	my($maptab)= @_;
	return {} unless($maptab);
	my($ok,$hin)= openRead($maptab);  die "ERR: reading $maptab" unless($ok);
	my %maps=(); my $nmap=0;
	
	##old# GenomeID        gespan  geor    AQueryID        quspan  match   qlen    cov     pid     path    indels  nexon   splice  aalen   offs    aamap   sense   tag
  ##updated 2014: add oid alt ID
  ## GenomeID        gespan  geor    AQueryID        quspan  match   qlen    cov     pid     path    indels  nexon   splice  aalen   offs    aamap   sense   oid     tag
  ## opt use slim map.attr
  ## AQueryID	cov	pid	splice/nexon	GenomeID:gespan:geor	path	oid
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
      my $nsx = 2; ## bad: ($nx<2)?0:$nspl/$nx; #? bad? nspl == 0 now no good
      my $asc=  ($nx > 3 and $nsx > 1.5 and $cov >= 90) ? 3 # certain?
              : ($nsx >= 1.5 and $cov >= 95) ? 2
              : ($nsx >= 1.8 and $nx > 2 and $cov >= $MINCOV_ANTISENSE) ? 1 : 0;
      $asense= ($asc>0)?"antisense":"nonsense"; # change 0 to maybe-antisense ? nonsense?
    }
		
		my $quals="";

if(MAPQ1805) {
    # chrmap: prefix or not? chrmap:99a,99i,1234l,9x,..
    my $mapq= sprintf "%da,%di,%dl,%dx,%s:%s:%s", ($cov,$pid,$ql,$nx,$chr,$cbe,$cor);
    $quals = $mapq;
} else {
	  $quals="cov=$cov,nexon=$nx"; # drop: splicex=$nsx, use nexon
}	  
    $quals.=",$asense" if($asense); #$quals.=",sense=$asense" if($asense);
		$quals.=",split=$msplit" if($msplit);
		$quals.=",paths=$mpath" if($mpath);
		
		$maps{$td}=  $quals; $nmap++;
	  if($oid){ for my $od (split",",$oid) { $maps{$od}= $quals; } }
		} close($hin);
		
	warn "#readMapinfo($maptab)= $nmap\n";
	return ($nmap>0)? \%maps : {}; 
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
		# $val.="\t".$mapq if($mapq);
		$oda{$pd}= $oda{$od}= $val; $ntr++;
	  $aaqual{$pd}= $aaqual{$od}= $aaq; #? skip readAaQual, but that adds gaps not here
		#? $pod{$od}=$pd; # vv?
		} close($hin);
		
	warn "#readPubidClass($pubids)= $ntr\n";
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
		
		$idbest{$od}{$idbest}=$pia;
		$antis{$od}="$idbest\t$pia" if($pia=~m/\-sense/);
		
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
	warn "#readTrClass($trclass)= $ntr\n";
	
	## if(UPD1912) { resolveAntis(\%antis,\%idbest,\%oda,\%aaqual); } # maybe not here, like reclassAlts calcs
	
	return (\%oda,\%aaqual);
}


sub resolveOneAnti { # ($aid,$trinfo,$apubid); UPD1912, maybe do in reclassAlts calcs
  my($aid,$trinfo,$apubid)= @_;
  my($ab,$upd)= (0,0);  
  $apubid ||= $aid;
  # return(0) if($DID_RESOLVESENSE);
  
  my $aval= $trinfo->{$aid} || $trinfo->{$apubid} || "0"; # == 	join"\t",$aaq,$pia,$cl,$aaref; # KEEP 4 cols
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



# readAaQual: inputset/*.aa.qual, or use evigenesub:fasize_nogaps(publicset/*.aa_pub.fa.gz,okc='A-WYZa-wyz',gapc='Xx*')
sub readAaQual { 
	my($trclass,$trinfo,$aaqualh)= @_;
	
  # fixme: update prior %aaqual from trclass/pubids, dont make new..
	unless(ref $aaqualh){ my %aaqual=(); $aaqualh= \%aaqual; }
	my $naa=0;
  my $aaqualf = "";
  
  ##old
  my($name,$path,$suffix) = fileparse($trclass,qr/\.\w*/); # suf-'.trclass' or suf = qr/\.\w*/
  my $aaqualf = "$path/inputset/$name.aa.qual"; 
  $aaqualf = "$path/publicset/$name.aa.qual" unless(-f $aaqualf);
  $aaqualf = "$path/$name.aa.qual" unless(-f $aaqualf);

  ## new
  # ($aaqualf=$trclass) =~ s,(\w+).trclass,inputset/$1.aa.qual,; 
  # $aaqualf =~ s,inputset,publicset, unless(-f $aaqualf);
  # $aaqualf =~ s,publicset,okayset, unless(-f $aaqualf);
  # $aaqualf =~ s,okayset/,, unless(-f $aaqualf);

  ## this is a mess, should update aa.qual fields in .trclass, others, to include nnnX gaps and reduce? aasize
  if(-f $aaqualf) {

if(UPD1908) {
    ($aaqualh,$naa)= getAaQual($aaqualf); # cdna_evigenesub;
    # aa.qval ==  "$alen,$pctcds,$acv,$aqual1,$codepot?"
    # acv == numeric of aqual1:complete/partial/.. adds gapbad tag for gaps
} else {
    my($ok,$hin)= openRead($aaqualf);  
    while(<$hin>) {
      next unless(/^\w/);
      my($od,$awnogap,$ngap,$aaq,$clen)=split; $naa++;
      my($aw)= $aaq=~/(\d+)/;
      if($awnogap < $aw) {
        $aaq =~ s/$aw/$awnogap/; $aaq.=",${ngap}X";
          # this isn't effective for publicset/ ids when trinfo has orig ids.. fixme
        if(my $tri= $trinfo->{$od}) { ## change trinfo.aq value?
          $tri =~ s/^\S+/$aaq/; $trinfo->{$od}= $tri;
        }
      }
      $aaqualh->{$od}=$aaq; #? all or just awnogap?
    } close($hin);
}
    
  }  else {
    my $aaseq = "$path/publicset/$name.aa_pub.fa"; 
    #or (my $aaseq=$trclass) =~ s,(\w+).trclass,publicset/$1.aa_pub.fa,; 
    $aaseq .= ".gz" unless(-f $aaseq); 
	  if(-f $aaseq) {
	    $aaqualf= $aaseq; #debug
	    my $fasizes= fasizes_nogap($aaseq,'amino');
	    foreach my $od (sort keys %$fasizes) {
        my($awnogap,$aw,$ngap)=split "\t",$fasizes->{$od}; $naa++;
        if($awnogap < $aw) {
          # this isn't effective for publicset/ ids when trinfo has orig ids.. fixme
          if(my $tri= $trinfo->{$od}) { ## change trinfo.aq value?
            $tri =~ s/^\d+/$awnogap/; 
            $trinfo->{$od}= $tri;
          }
        }
        $aaqualh->{$od}=$awnogap; #? all or just awnogap?
	    }
	  }
	}
	warn "#readAaQual($aaqualf)= $naa\n"; # if debug
	return $aaqualh; # scalar(%aaqual)? \%aaqual : undef;
}

=item revisePublicset

	need logic from evgmrna2tsa to do in compatible way (several tables to update) ..
	revise publicset/$trname.{pubids,mrna_pub.fa,cds_pub.fa,aa_pub.fa,ann.txt} ..
	also revise .trclass okay/drop 
	# my @keepids= sort keys %{$keepdrop{'keep'}}; # new pd, od, old pd ..
	# my @dropids= sort keys %{$keepdrop{'drop'}}; # new pd, od, old pd ..

	## evgmrna2tsa2.pl	
	my($pubmrna,$npm)	= make_pubseq($cdnaseq,'mRNA',\%annothash);
	my($pubaa,$npa) 	= make_pubseq(makename($cdnaseq,'.aa'),'protein',\%annothash);
	my($pubcds,$npc)	= make_pubseq(makename($cdnaseq,'.cds'),'CDS',\%annothash);

=cut

sub revisePublicset 
{
	warn "#revisePublicset NOT READY\n";  

# 	my($pubd,@ft)= getFileset('publicset','mrna_pub.fa|cds_pub.fa|aa_pub.fa');   
# 	foreach my $inseq (@ft) {
# 		my($drop,$nin,$fkeep,$fchange,$fdrop)=(0) x 9;
# 		my($ok,$inh)= openRead($inseq);
# 		(my $outfa=$inseq) =~ s/_pub/_pubrealt/;  #realt_pubids
# 		my $outh; $ok= open($outh,'>',$outfa) if($ok);
# 		while(<$inh>) { 
# 			if(/^>(\S+)/) {  my $d=$1; # old pubid
# 				$drop= ($keepdrop{'drop'}{$d}) ? 1 : 0;
# 				my $nd= $newid{$d}; if($nd and $nd ne $d and not $drop) { s/>$d/>$nd/; $fchange++;}
# 				$nin++; $fkeep++ unless($drop); $fdrop++ if($drop);
# 			} 
# 			print $outh $_ unless($drop);
# 		}
#   	close($inh); close($outh);
# 		#loggit(0,"revised $inseq to $outfa: nin=$nin, nkeep=$fkeep, nchangeid=$fchange"); 
# 	}
	
# 	my($annotab, $ntra1, $ncdsa1) 
# 		= make_annotab($cdnaseq,$trclass,$skiptrrun); # add main/alt pub ids, other geneinfo 
# 	loggit(0,"revised publicset: ",$pubmrna,$pubaa,$pubcds,$annotab); 

}


sub reclassAlts {
	my($outh, $pubids,$trinfo,$aaqual,$mapinfo,$exoninfo)= @_;  
	
  %trscore=%trv=%trline= (); # globals for putgene()
	my($nin,$nout,$ndrop,$nkeep,$nrenum,$ngdiff,$ngene,$ncols)=(0) x 10;
	my($ok,$pubidh)= openRead($pubids); die "ERR: reading $pubids" unless($ok);
	
	while(<$pubidh>) {
		#pubids.hdr:orig  Public_mRNA_ID originalID PublicGeneID AltNum
		#pubids.hdr:ext   Public_mRNA_ID originalID PublicGeneID AltNum Class  AAqual pIdAln Notes
		#1805upd add Oids col9
		#Public_mRNA_ID originalID      PublicGeneID    AltNum  Class   AAqual  pIdAln  Notes   Oids
		
		unless(/^\w/) {  # BUT print header , adding new col names
		  if(/^#Pub/) { 
		    s/$/\tClass\tAAqual\tpIdAln\tNotes/ unless(/\tClass/); 
		    $ncols=scalar( split"\t", $_ );
		    print $outh $_; }
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
			if($lgd) { my($aout,$gdiff,$anum,$akeep,$adrop)= putgene($outh,$outvalsh); 
				$ngene++; $ngdiff++ if($gdiff); $nrenum+=$anum; $nout+=$aout; $ndrop+=$adrop; $nkeep+=$akeep; }
			%trscore=%trv=%trline= (); # globals for putgene()
		}
		
		my $trval= $trinfo->{$od} || $trinfo->{$pd}||"0"; # == 	join"\t",$aaq,$pia,$cl,$aaref; # KEEP 4 cols
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
        if($aarefin =~ /aaref:[1-9]/) { 
          # AcmpanEVm000001t1 pubido: aaref:5143,AT3G02260.2.refbest  trclass: 0,0,aacons:2,pflag:0; 
          # pubidnu: aaref:5143,AT3G02260.2.refbest,,exonch:.. no aacons
          my($aacons)= ($aaref=~m/aacons:(\w+)/)?$1:0;  #UPD20mar
          if($aacons and not ($aarefin=~/aacons:/)) { $aarefin .=",aacons:$aacons"; }
          $aaref= $aarefin;  #? and aaq UPD20Mar: lost aacons: not in pubids.old, in trclass
          }
		  }
			$trval= join"\t",$aaq,$pia,$cla,$aaref if($fixval); # error, remake $trval?
	  }

		my($cdsantisense)=0;
		if(UPD1912) {
		  ($cdsantisense)= ($pia=~m,\-sense,)?-1:0; # UPD1912: add explicit revaa flag for pubset annot uses
		  # compare w/ chrmap antisense, one way is true, other false ?
      ## FIXME: cdsantisense showing up on single tr/locus and on largest aa of alts, versus
      ##   short alts w/ fwd sense .. should swap these, turn off for 1 tr/locus
      if($cdsantisense < 0 and not $DID_RESOLVESENSE) {
        my($didra,$trvalupd)= resolveOneAnti($od,$trinfo,$pd); # ($aid,$trinfo,$apubid); UPD1912, maybe do in reclassAlts calcs
        if($didra){
          $trval= $trvalupd;
          ($aaq,$pia,$cla,$aaref)= split"\t",$trval;  
          ($cdsantisense)= ($pia=~m,\-sense,)?-1:0; # UPD1912: add explicit revaa flag for pubset annot uses
        }
      }
    }
    
		my($aw,$ap,$ac)=split",",$aaq;
		my($pi,$pa)=split"/",$pia; 
    ## UPD19 use pflag bits, asmrna_dupfilter3c.pl		
    use constant{ kAATINY =>1, kAAUTRBAD =>2, kAAUTRPOOR =>4, kAADUP =>8, kAAGAPS =>16, kNONCODE =>32, }; 
    my($pflag)= ($trval=~m/pflag:(\w+)/)? $1 :0;  # pflag from, see asmrna_dupfilter3c.pl
    
    my($aacons)=0;
    if(AACONS) { ($aacons)= ($trval=~m/aacons:(\w+)/)?$1:0; } # UPD1908
    
    my($hasmap,$antisense,$msplit,$mcov,$mexon)= (0,0,0,99,2); # defaults no mapinfo val
    if(ref $mapinfo) {
		  ## mapinfo == "cov=$cov,nexon=$nx,splicex=$nsx,sense=$asense" : keep any more than antisense?
		  ## upd1805: format now as chrmap:99a,99i,1758l,9x,chr23:42516712-42528502:+
		  my $mfo= $mapinfo->{$od} || $mapinfo->{$pd} || ""; # antisense info; which id here? pd or new realt pd or od ?
		  if($mfo) { 
		  $hasmap=$mfo; # was =1; 
		  $antisense=($mfo =~ /antisense/)?-1:$antisense; # or -1:0 or "antisense" ?
		  #or $antisense=($mfo =~ /antisense|nonsense/)?-1:0; # nonsense = uncertain antisense
		  $msplit=   ($mfo =~ /split=(\w+)/)?$1:$msplit; 
	    $mcov=     ($mfo =~ /cov=(\d+)/)?$1:$mcov; # NOPATH: paths=-1 or cov=0
		  $mexon=    ($mfo =~ /nexon=(\w+)/)?$1:$mexon; 
		  }
	  }
	  
	  ## commonexon qual: modify to give qual as pct of best alt (longest common exon chain),
	  ## then use pct-xchain to demote those w/ low/0 score, or promote those w/ good chain
	  ## exonch val now: 100%cs,100%cx,20c:...,cx:nx/nxt
	  my($hasxchain, $xchainval)=(0,0);
	  if(ref $exoninfo){
	    my $xfo= $exoninfo->{$od} || $exoninfo->{$pd} || "";
	    if($xfo) {  
	      $hasxchain= $xfo; $xchainval= $xfo =~ m/^(\d+)/?$1:0; 
        $trval .=",exonch:$hasxchain"; #  
	      }
	  }
	  ## xchainval:  add goodness increment, nothing for poor, 
	  ## eg:   trscore += trscore * 0.50 * [100/100, 90/100, 50/100, 0/100 ] 
	  # if($hasxchain){  $trscore += int(0.5 + $trscore * $pPromoteGoodExonchain * $xchainval ); }
	  
	  if($antisense<0) { #? move to mapqual? or keep w/ aaqual?
	    if($cdsantisense<0) { } #UPD1912 what? both tests likely mean bad orf call .. discard this one?
	  	$trval=join("\t","$aaq,antisense",$pia,$cla,$aaref); # more sensible to add to aaqual field
	  } elsif($cdsantisense < 0 and not $hasmap) {	  	
	  	#NO. cdsanti isnt aa-qual, but comparison to other tr quality.  leave as pia/-sense
	  	#NO $trval=join("\t","$aaq,antisense",$pia,$cla,$aaref); ##UPD1912 : likewise add to aaqual field
	  }
	  
    my $mapqual= "mapok"; # need as default for no map info?
    if($hasmap) {
if(MAPQ1805) {
      $mapqual= $hasmap;
      $trval .=",chrmap:$mapqual" unless($trval=~m/chrmap:/); # new chrmap:99a,99i,1758l,9x,chr23:42516712-42528502:+
} else {    
      $mapqual= ($msplit or $mcov<$MINCOV)?"mappoor":"mapok"; # or mexon<2 ?
      $mapqual.=",${mcov}a,${mexon}x"; 
      $mapqual.=",Split$msplit" if($msplit);
      $trval .=",mapq:$mapqual"; # trinfo may have chrmap:vals
}      
    }
	  
		
		my $samecds= ($pi>=$SAME_PI and $pa>=$SAME_PA)?1:0;
		my $diffcds= ($pa<=$DIFF_PA)?1:0; # skip ident score here? alt is diff enough when align score low
		   $samecds=0 if($diffcds); # shouldnt need..
	
	  use constant kAAComplete => 3;	   
	  # bug here icla bad w/ Noncode, drops aaw=139 + Noncode see below
		my $icla= ($ac =~ /complete/)?kAAComplete: ($ac =~ /partial[35]/)? kAAComplete-1: kAAComplete-2; 
		$icla-- if($ac =~ /gapbad/);
		$icla-- if($cla =~ /frag/); ## check new class tags, 'nc' noncode, maybe? ..

		#old.dont demote poor# $icla-- if($ac =~ /utr(bad|poor)/ or $od =~ /utrorf/); ## 201405: utrbad/poor/orf class demote
	  if(UPD1908) {
	    ; #skip: $icla-- if($od =~ /utrorf/); ## UPD1908: drop utrorf only, utrbad may not be bad.. utrorf may be good also
	  } else {
	    $icla-- if($ac =~ /utrbad/ or $od =~ /utrorf/); ## 201405: utrbad/poor/orf class demote
    }
    
    ## FIXME: *** $aw aa-size from trclass is wrong, no gap adjust; 
    ##     original trclass/pubids sorting is right from using aa.qual table **
    ## FIXME2: add $aaref : must-keep if $aaref =~ /ref(best|good|ok)/, maybe-keep if $aaref=~ bitscore,refid
    ##     add to trv{pd} BUT, dont change trscore sort order for must-keep. ??

    my $aaqual=  $aaqual->{$od} || $aaqual->{$pd} || $aaq || "";
    my ($aasize)= ($aaqual =~ m/^(\d+)/)?$1:($aaq =~ m/^(\d+)/)?$1:0; #?? aaqual bad?
    if($aasize>0 and $aasize<$aw) { $aw=$aasize;  } # fixme: update trinfo/trline ??

    #UPD1908: $aaqual from getAaQual() is "$alen,$pctcds,$acv,$aqual1,$codepot"
    # $aqval acv == numeric of aqual1:complete/partial/.. adds gapbad tag for gaps
    my($codepot,$cpval,$cdsor,$aqval)=(0) x 9; # FIXME get codepot number
    if(UPD1908) { # == CODEPOT1607
      # ALSO new class tags: althincmaybe  'nc' = noncode/maybe, should have pflag & kNONCODE
      if($aaqual=~m/,(Code|Noncode|Unknown)/) { $codepot=$1;
        ($cpval)= ($aaqual =~ m/,(Code|Noncode|Unknown).([^,;\s]+)/)?$2:0;
        }
      if($pflag & kNONCODE) { $codepot ||= "Noncode"; }
      if($aaqual =~ m/,\d+\-\d+:(.)/) { $cdsor=$1; } # cdsor == ,256-12915:+  ,2491-176:-
      my @aq= split",",$aaqual; $aqval=$aq[2] if(@aq>3); # in -3..3 range ? near equiv of $icla ??

      if(0) { # update? per trclass_resolve_strandmix.pl
      $trval .= ",codepot:$cpval" unless( $trval =~ /codepot:/ ); # . substr($codepot,0,1);
      $trval .= ",rnaor:$cdsor" if($cdsor and  $trval !~ /rnaor:/ ); #? debug?
      }
    }
    
    my $trscore= $aasize || (99999 - $ti); #?? 99999 - ti bad default? fail if no aasize?

    ##  antisense often from partial align 1st intron mismap; need opts? cancel if utrbad? same cause?
    if(UPD1912) {
    if($antisense < 0 or $cdsantisense < 0) {
      # not sure should use cdsantisense, unless matches chr antisense, either fwd/ref cds may be true
      my $pdemoteanti= # def=0.95
          ($antisense < 0 and $cdsantisense < 0) ? $pDemoteAntiSense * $pDemoteAntiSense #? =0.902, or other
        : ($antisense < 0) ? $pDemoteAntiSense
        : 1; #?? ($cdsantisense < 0) ? 0.99 : 1;
      $trscore = int(0.5 + $trscore * $pdemoteanti); 
      }
    } else {
		$trscore = int(0.5 + $trscore * $pDemoteAntiSense) if($antisense<0); 
		}
		
		if($msplit) { $trscore = int(0.5 + $trscore * $pDemoteSplitmap); }
		elsif($mcov<$MINCOV) { $trscore = int(0.5 + $trscore * 0.98); } 	#ignore?# 
		
		## $trscore= $trscore * (97+$icla)/100 if($icla < kAAComplete); #?? this can change main/alt sort order.. want?
		## for icla==0, paapoor = 0.9 * 0.9 * 0.9 * 0.9 == 0.65, ok
		if($icla < kAAComplete) {
		  my $paa=1; for(my $i=$icla; $i < kAAComplete; $i++) { $paa= $pDemoteAaPoor * $paa; }
		  $trscore = int(0.5 + $trscore * $paa);
		}
    
		## but see above note: dont change trscore ? use $pPromoteAaRef == 0 for no change ?
		my($vaaref)= ($aaref=~m/aaref:(\d+)/)?$1:0; # dont have proportion, use raw/bitscore as weight?
		if($vaaref>=9) { $trscore += int(0.5 + $pPromoteAaRef * $vaaref); } # def: $pPromoteAaRef=0.50 assumes bitscore

		if(AACONS and $aacons>0) { $trscore += int( 0.5 + $trscore * $pPromoteAaCons * $aacons); } # with pPromoteAaCons ~ 0.10 * 2,3,4
		
    my $nocode=0;
    if(UPD1908 and $codepot) { # == CODEPOT1607
      $nocode=($codepot =~ /^Noncode/i)?1:0;
      ## aasize > 139 : use SMALLMAXAA param $smallclassmax= $ENV{'SMALLMAXAA'} : 100;  
      my $iscode=( $vaaref >= 9  or $aasize > 139 or ($aasize > 99 and $icla >= kAAComplete)); # block
      if(UPD1912) { $iscode=0 if($cdsantisense < 0 and $vaaref < 9); }
 		  $trscore = int(0.5 + $trscore * $pDemoteNoncode) if($nocode and not $iscode);
    }

	  if(UPD1908 and $hasxchain){ $trscore += int(0.5 + $trscore * $pPromoteGoodExonchain * $xchainval ); }

		my $trline=join("\t",$pd,$od,$gd,$ti,$trval); # SHOULD be 8 cols for pubids, trval has 4 col
		$trline.= ",tscore:$trscore" if($debug); #debug
		$trline.= "\t".join("\t",@xcol) if(@xcol>0); #upd1806
		$trline{$pd}= $trline;

		$trscore{$pd}= $trscore; ## add icla to score? or $aw * (97+$icla)/100
		
		#UPD1912? if($cdsantisense < 0 and not $hasmap) : should cdsantisense be set in trval? will drop some where revaa is true value
		# note: aaref now has other fields appended .. drop? aacons:2,pflag:32
		$aaref =~ s/(aacons|pflag).*//;  $aaref ||= 0;
		$trv{$pd}=join "\t", $cla,$aw,$icla,$samecds,$diffcds,$antisense,$aaref,$mapqual,$xchainval,$aacons,$nocode,$pflag;  # UPD19
		
		$lgd=$gd;  
	} close($pubidh);
	
	if($lgd) { my($aout,$gdiff,$anum,$akeep,$adrop)= putgene($outh,$outvalsh); 
		$ngene++; $ngdiff++ if($gdiff); $nrenum+=$anum; $nout+=$aout; $ndrop+=$adrop; $nkeep+=$akeep; }
	## add summary counts:  ndrop, nrenum, ngenechanged, nin, nout
	return($nin,$nout,$ndrop,$nkeep,$nrenum,$ngdiff,$ngene);
}

my $didtrvhdr=0;
sub putgene {
  my($outh,$outvalh)=@_; 
  my($changed,$renum)=(0,0); # 
  $outvalh ||= 0;
  
  #? use orig t1 to start, change only if isbad?
  # option to not demote orig t1 too much, nor drop? eg to move split mains from top spot but keep otherwise
  my @tr= sort{ $trscore{$b} <=> $trscore{$a} or $a cmp $b } keys %trscore;
  my $t1= shift @tr or return 0; # main/longest tr : is this changing from orig main?
  if($KEEPMAINatTOP and not($t1 =~ /t1$/) ) {
    my $i1=-1; for(my $i=0;$i<$#tr;$i++) { if($tr[$i]=~m/t1$/){ $i1=$i; last; } }
    if($i1>=9) { my $t1=$tr[$i1]; splice(@tr,$i1,1); @tr=($t1,@tr); } #too high, shuffle down?
  }
  
	my @TRVCOLS=qw(class aasize aaqual samecds diffcds asense aaref mapqual exonchain aacons noncode pflag);
	#	$trv{$pd}=join "\t", $cla,$aw,$icla,$samecds,$diffcds,$antisense,$aaref,$mapqual,$xchainval,$aacons,$nocode,$pflag;  # UPD19
	
  my($t1cla,$t1aw,$t1icl,$ti1same,$t1diff,$t1anti,$t1aaref,$t1mapq,$t1exonc,$t1aacons,$t1noncode)= split"\t", $trv{$t1};
  ## if t1 is antisense, need to change below drops, add antisense -trscore ??
  $t1aw= $trscore{$t1} if($t1anti<0); ## reduce t1aw, use trscore?
  my $t1cull= ($t1cla =~ /^cull/)?1:0; # other cla?
  
  my (@drops,@keeps); # culls are in @keeps, but cull prefix to class
  my $ialt=1; 
  $t1aw=1 if($t1aw<1); ## $t1aw ||= 1;  #?? bug div by zero below
  my ($ltaw)=($t1aw);

	my $AABITS_MIN_NOCL = 69; # or what?
	my $AABITS_MIN_ALT  = 99; # or what?
	my $altreplacemain=""; my $mainisbad=0;
  my $aaconsmax= $t1aacons;
  
  if(CODEPOT1607) { ## UPD19 add this evd
  # $tcode=0 if($aw > 139 or ($aw > 99 and not $ispoor)); #?? which? many utrorf/utrbad for aa>99
  # $ispoor |= kNONCODE if($tcode =~ /^Noncode/); #?
	}
	
  foreach my $tr ($t1,@tr) {
    my ($drop,$keep,$cull)=(0,0,0); 
    $t1aw=1 if($t1aw<1);  # BUG below: division by zero .. cant be 
    
    if($tr eq $t1) { # check for noclass drop option
    	## NO: t1icl == number not trclass
    	$cull= $t1cull; # ($t1cla =~ /^cull/)?1:0; # other cla?
    	
      if($noclasscut>0 and $t1cla =~ /noclass/) {  # check for noclass drop option
        # UPD1908 add noncode val to this
        # FIXME: this drops many valid plant ref due to poor codepot vals
        # BUT it also cuts large subset of ~useless tinycds things, ie 50,000 for tick test set
        # need more reliable codepot calc, cptable from this-trasm longest cds (code) - utr (nocode)
        # ? use this to increase noclasscut val? $pDemoteNoncode   = $ENV{scalenonecode} ||0.80;  

        my $noclminaa= $noclasscut;
        if(NOCLASS_NONCODE_CUT and $t1noncode>0) { 
          # increase aacut by 20% if pdemote = 0.80, ie from 60 aa to 72 aa? more? 50% better
          #? $noclminaa= $noclminaa + (1-$pDemoteNoncode) * $noclminaa; 
          $noclminaa= NOCLASS_NONCODE_CUTAA if($noclminaa < NOCLASS_NONCODE_CUTAA); 
        }
        if($t1aacons>1) { $noclminaa -= 10 * ($t1aacons-1); }
        $drop=($t1aw < $noclminaa)?1:0;
        if($t1aaref=~/ref(best|good|ok)/) { $keep=1; $cull=$drop=0; } ##?? PROBLEM FAILING TEST 
        elsif($t1aaref=~/^(\d+)/) { my $bits=$1; $drop=0 if($bits>$AABITS_MIN_NOCL); } # what?
        
      } else { # should be main, but maybe not ..
      
        $keep=1; #UPD1908 : drop some main+alt loci when all alts are poor qual
        # eg. add $locusisbad, use similar qual to noclass drop, evd: aasize, aaref, noncode, aacons, exonc, mapq, 
        
        $cull=0 if($t1aaref=~/ref(best|good|ok)/);
        #FIXME: check main qual, reclass mainbad <> altgood
        ## mainisbad for Split also
        $mainisbad= ($t1icl < 3 or not ($t1mapq=~/mapok/))?1:0;
        $mainisbad=1 if($cull);
        $mainisbad=0 if($t1aaref=~/refbest/);
      }
      
    } else { # if($tr ne $t1)
      my ($tcla,$taw,$ticl,$tisame,$tidiff,$tianti,$tiaaref,$timapq,$tiexonc,$tiaacons,$tinoncode)= split"\t", $trv{$tr};
   		$taw= $trscore{$tr} if($tianti<0); ## reduce taw/use trscore?
      my $paw= $taw / $t1aw;  # BUG: division by zero .. cant be

      # FIXME1806: for t1culls, also cull alts, if paw <= 1, else reclass as new t1..
    	$cull= ($tcla =~ /^cull/ or ($t1cull and $paw <= 1))?1:0; # drop? other ecla?
      $cull=0 if($tiaaref=~/ref(best|good|ok)/);
      if(AACONS and $tiaacons>=$aaconsmax) { $keep=1; $cull=0; $aaconsmax=$tiaacons;  } # maybe.. or not, should check vs t1/other alt aacons
      # my($ecla,$eclb)= $t1cla=~/^(\w*)(main|noclass|alt|part|noncode)/ ? ($1,$2) : ("",$cla);
      
      # 1702 problem: very short alts, aw=50aa vs t1aa=5000 aa,  pi=100, pa=50, should drop , bypassed here
      # ticl == t1icl, but tidiff == 1, tisame == 0 : same/diff cds class at fault for tiny alts
      # drop=1 if($paw < 0.20? 0.10? );
      ## ?? change: ($ticl>=$t1icl and $tisame and $paw < $DROPSHORT_AAFULL) to ($paw < $DROPSHORT_AAFULL)
      ## or add ($paw < $DROPTINYALT); DROPTINYALT = 0.1nnn == frag

      $drop=( ($paw < $DROPTINYALT) # 0.10d
          or ($ticl< $t1icl and $paw<$DROPSHORT_AAPART and not $tidiff) # 0.90d
      	  or ($ticl>=$t1icl and $tisame and $paw < $DROPSHORT_AAFULL)) ? 1 :0;  # 0.20d
      ## ** ^^ need to check if dropping all slightly short partials is good idea; 
      ## check classify value of these:  ($ticl< $t1icl) ; $tidiff; $tisame; $tianti; $tiexonc
      
      $drop=1 if($tianti < 0 and $paw < $DROPSHORT_ANTISENSE); ## antisense check .. this is chranti, not cdsanti
			
			if(UPD1908) { $drop=1 if($tiexonc < 25 and $paw < $DROPSHORT_BADEXONS); }
			
			# check this altreplace classing
      $altreplacemain=$tr if( $mainisbad and not $altreplacemain 
        and not ($drop or $cull)
        and ($timapq=~/mapok/) and $ticl >= $t1icl and $paw >= $ALTGOODasMAIN);

			## aablastp ref score keeps some drops, doesn't change alt order
			## HOWEVER, may have several alts same aaref score.. modify maxsame check? add trhoscore?
			if($tiaaref) { ## and $drop < NOT this
			 if($tiaaref=~/ref(best|good|ok)/) { $keep=1; $cull=$drop=0; } # cull cancel ??
			 elsif($tiaaref=~/^(\d+)/) { my $bits=$1; $drop=0 if($bits>$AABITS_MIN_ALT); }
      }
      
      if($ialt > $maxsame and not $drop and not $keep) { $drop=1 if($taw > $ltaw-12); } # $taw == $ltaw or 
      ## ^ this classing need checking. may be dropping valid alts
      
      $ltaw=$taw; # unless drop ??
    }
    
    $cull++ if($drop and $NODROPS and not $keep); # cull instead of drop? try it
    $drop=0 if($keep or $NODROPS); # fix what?
    ## new col for aaref? missing often; append to Notes col
    # FIXME: pubids col order differs from trline, change trline
    my $trline=$trline{$tr}; 
    if($drop) {
      push @drops, $trline;
    } else { 
      if($altreplacemain and $mainisbad) { 
        my $tr1l= shift @keeps; unshift @keeps, $trline; $trline=$tr1l; $mainisbad=0; $altreplacemain=0; 
      }
      
      # two kinds of culls: b. added here, a. from input pubid table (from keepdrop class table)
      # b.cull may be canceled above, not a.cull ?
      
      if($cull) { # change evgclass in trline?
        my @trl=split"\t",$trline; 
        # my($pd,$od,$gd,$ti,$aq,$pia,$cla,$notes,@xcol)= @trl; # FIXME: changing col order
        unless($trl[6] =~ /^cull/){ $trl[6]="cull".$trl[6]; $trline= join"\t",@trl; }
      } 
    	push @keeps, $trline;
    	$ialt++; # need to count good alts here for maxsame test. reset for output
    }
    
    if(UPD1908 and $outvalh) {
      my $kept=($drop||$cull)?0:"keep";
      print $outvalh join("\t","PublicID", @TRVCOLS,"keep","trscore")."\n" unless($didtrvhdr++); # $trv{$tr},"\t",$trscore{$tr},"\n";
      print $outvalh join("\t",$tr,$trv{$tr},$kept,$trscore{$tr}),"\n";
    }
        
  }
	
	# maybe recheck @keeps here for excess same alts
	# save drop/keep pd,od,oldpd in hash for other uses, ie rewrite publicset files

	## UPDATE unless($arenum) put newID in notes column, and maybe put ialt into ti column?
	## FIXME: change classes: 1st class to 'main', change 'dropmain' to 'dropalt', main>alt if not 1st
	##  .. and preserve oldclass:$cla in Notes 
  ## BAD Output :  ^,newid:PitaTv1R000040t115  n=110585; CHOMP bug ??  on trline?? 
	
	$ialt=0; # reset for output.
  foreach my $trline (@keeps) {
    my($pd,$od,$gd,$ti,$aq,$pia,$cla,$notes,@xcol)=split"\t",$trline; 
    
    # my($cullc)= $cla=~/^(\w+)(?:main|alt|noclass|part)/ ? $1 : "";
    my($ecla,$eclb)= $cla=~/^(\w*)(main|noclass|alt|part|noncode)/ ? ($1,$2) : ("",$cla);
    $ecla ||="";

   	$ialt++; my $oldid=$pd; 
   	# FIXME1806: 'aaref:0,0,' is new bug from where?		
   	if($notes) { if($notes =~ /^(\d+)/){ my $av=$1; $notes="aaref:$notes" if($av>0); } $notes.=","; }
   	else { $notes=""; } # not "0"
   	
   	if($ialt == 1 and $eclb=~/alt/) { $eclb="main"; } # FIXME
   	elsif($ialt > 1 and $eclb=~/main/) { $eclb="althim"; } # FIXME
   	
   	if($ti != $ialt) { my $newd=$gd.'t'.$ialt; 
  		if($arenum) { $notes.="oldid:$pd,"; $pd=$newd; $ti=$ialt; $renum++; $changed++; }
   		else { $notes.="newid:$newd,";  $ti=$ialt; } # change ti or not?
   	}
   	
    #bad: $cla=$cullc.$cla;
    # FIXME: pubids col order differs from trline, change trline
    print $outh join("\t",$pd,$od,$gd,$ti,$ecla.$eclb,$aq,$pia,$notes||'.',@xcol)."\n"; 
    my @ids=($pd,$od); if($oldid ne $pd) { push @ids, $oldid; $newid{$oldid}=$pd; } $keepdrop{'keep'}{@ids}=$trline; #?
  }
  
  foreach my $trline (@drops) {
    my($pd,$od,$gd,$ti,$aq,$pia,$cla,$notes,@xcol)=split"\t",$trline; 
    my($ecla,$eclb)= $cla=~/^(\w*)(main|noclass|alt|part|noncode)/ ? ($1,$2) : ("",$cla);

    $ialt++; my $oldid=$pd;
   	# FIXME1806: 'aaref:0,0,' is new bug from where?		
   	if($notes) { if($notes =~ /^(\d+)/){ my $av=$1; $notes="aaref:$notes" if($av>0); } $notes.=","; }
   	else { $notes=""; } # not "0"
   	if($ti != $ialt) { my $newd=$gd.'t'.$ialt; 
  		if($arenum) { $notes.="oldid:$pd,"; $pd=$newd;  $ti=$ialt; $renum++; }
   		else { $notes.="newid:$newd,";  $ti=$ialt; } # change ti or not?
   	}
   	
    print $outh join("\t",$pd,$od,$gd,$ti,'drop'.$eclb,$aq,$pia,$notes||'.',@xcol)."\n"; 
    my @ids=($pd,$od); if($oldid ne $pd) { push @ids, $oldid; $newid{$oldid}=$pd; } $keepdrop{'drop'}{@ids}=$trline; #?
    $changed++;
  }
  return ($ialt,$changed,$renum, scalar(@keeps), scalar(@drops));
}

__END__

=item value table correlations

vtabf="try9weedevgc.reclass3h.valref"
wtabf="try9humanevgc.reclass3h.valref"
wtab=read.delim(wtabf,header=T)
wtabc= wtab[,vtabcols]; rownames(wtabc)= wtab[,1]
wtabd= wtabc; wtabd$notref= 2 - as.numeric(wtabd$notref); colnames(wtabd)[10]= "isref"
wcor= cor(wtabd)

wtrain=sample(1:125469,50000)
wlda= lda(notref ~ ., wtabc, prior = c(1,1)/2, subset = wtrain)
wpred=  predict(wlda, wtabc[-wtrain, ])
table( wtabc[wpred$class == "notref","notref"])

=== weed data ===
> head(vtab)
         PublicID  class aasize aaqual samecds diffcds asense    aaref mapqual exonchain aacons noncode pflag keep
1  AtEv3g000001t3 althi1   5393      3       1       0      0 aacons:3   mapok       100      3       0     0    1
2  AtEv3g000001t1   main   5438      3       1       0      0  pflag:0   mapok       100      0       0     0    1
3  AtEv3g000001t2 althi1   5431      3       1       0      0  pflag:0   mapok       100      0       0     0    1
4  AtEv3g000001t4 althi1   5360      3       1       0      0  pflag:0   mapok       100      0       0     0    1
5  AtEv3g000001t5 althi1   5231      3       1       0      0  pflag:0   mapok       100      0       0     0    1
6 AtEv3g000001t13 althi1   4803      3       1       0      0 aacons:2   mapok        90      2       0     0    1
  trscore notref
1  357561  isref
2  277338 notref
3  276981 notref
4  273360 notref
5  266781 notref
6  265144 notref


head(vtabd)
                aasize aaqual samecds diffcds exonchain aacons noncode pflag trscore isref
AtEv3g000001t3    5393      3       1       0       100      3       0     0  357561     1
AtEv3g000001t1    5438      3       1       0       100      0       0     0  277338     0
AtEv3g000001t2    5431      3       1       0       100      0       0     0  276981     0
AtEv3g000001t4    5360      3       1       0       100      0       0     0  273360     0
AtEv3g000001t5    5231      3       1       0       100      0       0     0  266781     0
AtEv3g000001t13   4803      3       1       0        90      2       0     0  265144     0

> head(vcor,20)
           aasize  aaqual samecds diffcds exonchain  aacons noncode   pflag trscore   isref
aasize     1.0000  0.0923  0.1297 -0.1771    0.3802 -0.0255 -0.2002 -0.2052   0.843  0.1138
aaqual     0.0923  1.0000 -0.0843  0.0216    0.1883  0.2557  0.0632  0.0654   0.181  0.0964
samecds    0.1297 -0.0843  1.0000 -0.5343    0.0926  0.1391 -0.1061 -0.1099   0.152  0.0838
diffcds   -0.1771  0.0216 -0.5343  1.0000   -0.1633 -0.0903  0.2103  0.2157  -0.194 -0.0725
exonchain  0.3802  0.1883  0.0926 -0.1633    1.0000  0.1684 -0.1818 -0.1844   0.665  0.2738 *
aacons    -0.0255  0.2557  0.1391 -0.0903    0.1684  1.0000  0.0318  0.0379   0.179  0.2418 *
noncode   -0.2002  0.0632 -0.1061  0.2103   -0.1818  0.0318  1.0000  0.9893  -0.150 -0.0463
pflag     -0.2052  0.0654 -0.1099  0.2157   -0.1844  0.0379  0.9893  1.0000  -0.154 -0.0468
trscore    0.8433  0.1806  0.1524 -0.1937    0.6645  0.1792 -0.1499 -0.1538   1.000  0.2436 +
isref      0.1138  0.0964  0.0838 -0.0725    0.2738  0.2418 -0.0463 -0.0468   0.244  1.0000

====  human data ===

> head(wtabd)
                   aasize aaqual samecds diffcds exonchain aacons noncode pflag trscore isref
Hstra9Evm000001t15   7425      3       1       0       100      3       0     0  492303     0
Hstra9Evm000001t20   7326      3       1       0       100      2       0     0  448341     0
Hstra9Evm000001t1    7478      3       1       0       100      0       0     0  381378     1
Hstra9Evm000001t2    7476      3       1       0       100      0       0     0  381276     1
Hstra9Evm000001t3    7471      3       1       0       100      0       0     0  381021     1
Hstra9Evm000001t4    7468      3       1       0       100      0       0     0  380868     0
> wcor= cor(wtabd)
> head(wcor,20)
          aasize  aaqual samecds diffcds exonchain aacons noncode  pflag trscore  isref
aasize     1.000  0.1447  0.1899  -0.234     0.308 -0.073  -0.106 -0.118   0.875  0.088
aaqual     0.145  1.0000 -0.0045  -0.090     0.215  0.067  -0.041 -0.041   0.184  0.168
samecds    0.190 -0.0045  1.0000  -0.630     0.093  0.136  -0.047 -0.056   0.182  0.136
diffcds   -0.234 -0.0904 -0.6296   1.000    -0.191 -0.073   0.090  0.103  -0.225 -0.125
exonchain  0.308  0.2149  0.0930  -0.191     1.000  0.071  -0.105 -0.114   0.561  0.266  *
aacons    -0.073  0.0669  0.1365  -0.073     0.071  1.000   0.084  0.095   0.064  0.335  *
noncode   -0.106 -0.0412 -0.0473   0.090    -0.105  0.084   1.000  0.986  -0.087 -0.049
pflag     -0.118 -0.0412 -0.0558   0.103    -0.114  0.095   0.986  1.000  -0.097 -0.052
trscore    0.875  0.1839  0.1821  -0.225     0.561  0.064  -0.087 -0.097   1.000  0.188
isref      0.088  0.1676  0.1361  -0.125     0.266  0.335  -0.049 -0.052   0.188  1.000
 

====
print(vlda)
Call:
lda(notref ~ ., data = vtabc, prior = c(1, 1)/2, subset = vtrain)

Prior probabilities of groups:
 isref notref 
   0.5    0.5 

Group means:
       aasize aaqual samecds diffcds exonchain aacons noncode pflag trscore
isref     435   2.93   0.640   0.119      83.2  1.967 0.00872 0.316   23380
notref    321   2.79   0.492   0.228      46.5  0.809 0.04331 1.455   10232

Coefficients of linear discriminants:
                 LD1
aasize     0.0014370
aaqual     0.0097449
samecds   -0.1704066
diffcds   -0.0436147
exonchain -0.0121088
aacons    -0.3914457
noncode   -0.5874602
pflag      0.0269834
trscore   -0.0000459

vpred=  predict(vlda, vtabc[-vtrain, ])

* Not very good discrim ..
> table( vtabc[vpred$class == "isref","notref"])
 isref notref 
  6270  72681 
> table( vtabc[vpred$class == "notref","notref"])
 isref notref 
 18299 207637 

==== human ====

print(wlda)
Call:
lda(notref ~ ., data = wtabc, prior = c(1, 1)/2, subset = wtrain)

Prior probabilities of groups:
 isref notref 
   0.5    0.5 

Group means:
       aasize aaqual samecds diffcds exonchain aacons noncode pflag trscore
isref     570    3.0    0.62    0.18        77   1.57 0.00045 0.067   26989
notref    436    2.8    0.42    0.35        48   0.49 0.01912 0.707   14056

Coefficients of linear discriminants:
                LD1
aasize     0.000061
aaqual    -0.603972
samecds   -0.354711
diffcds    0.006556
exonchain -0.013521
aacons    -0.728780
noncode   -2.025442
pflag      0.091224
trscore   -0.000004

> table( wtabc[wpred$class == "isref","notref"])
 isref notref 
  4407  28111 
> table( wtabc[wpred$class == "notref","notref"])
 isref notref 
 12492  80459 

=== attr means ===

rbind(isref=colMeans(wtabd[wtabd$isref==1,]),notref=colMeans( wtabd[wtabd$isref==0,]))

human  aasize aaqual samecds diffcds exonchain aacons  noncode  pflag trscore isref
isref     559   2.99   0.619   0.179      76.8  1.551 0.000473 0.0633   26373     1
notref    437   2.79   0.420   0.352      48.2  0.493 0.018771 0.6998   14048     0

plant  aasize aaqual samecds diffcds exonchain aacons noncode pflag trscore isref
isref     436   2.93   0.648   0.117      83.2  1.981 0.00948 0.343   23448     1
notref    321   2.78   0.494   0.226      46.6  0.806 0.04294 1.440   10273     0

=== attr quantiles ===
tab=wtabc
rbind( isref=quantile( tab[tab$notref=="isref",va], probs=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1) ), notref=quantile( tab[tab$notref=="notref",va], probs=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1) ))

va=aasize
human  0% 10% 20% 30% 40% 50% 60% 70% 80%  90% 100%
isref  44 177 241 302 364 432 513 614 763 1052 7478
notref 18 101 144 192 244 305 378 474 611  883 7468

plant  0% 10% 20% 30% 40% 50% 60% 70% 80% 90% 100%
isref  31 152 211 267 319 371 431 501 601 787 5393
notref 19  91 128 165 205 250 303 369 461 628 5438

exonchain
> quantile( wtabc[wtabc$notref=="isref","exonchain"], probs=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1) )
human   0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
isref    0   20   50   70   83   92  100  100  100  100  100 
notref   0    0    0   20   35   50   63   75   88  100  100 

plant   0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
isref    0   35   69   83   93  100  100  100  100  100  100 
notref   0    0    0   20   33   50   60   71   83  100  100 

aacons
> quantile( wtabc[wtabc$notref=="isref","aacons"], probs=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1) )
human   0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
isref    0    0    0    0    0    2    2    3    3    3    3 
notref   0    0    0    0    0    0    0    0    2    2    3 

plant   0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
isref    0    0    0    0    2    2    3    3    4    4    4 
notref   0    0    0    0    0    0    0    2    2    3    4 

va=trscore
human  0%  10%  20%   30%   40%   50%   60%   70%   80%   90%   100%
isref  42 3131 8159 12251 16116 20298 25047 30855 38964 54478 381378
notref 17  146  527  2186  4463  7245 10811 15403 22032 34805 492303

plant  0%  10%  20%   30%   40%   50%   60%   70%   80%   90%   100%
isref  25 3921 8619 12342 15759 19437 23670 28356 34988 46563 362508
notref 16  137  340  1607  3568  5858  8560 12098 17064 25947 327675

=cut

=item work

evigene/scripts/rnaseq/asmrna_altreclass.pl -maxalt=19 -altrenum  -trclass kfish2evg367mixx.trclass -out 
#altreclass: nin/out=499047/499047, ndrop/keep=141000/358047, nrenum=297357, ngenediff=35964/149901

## adding gmapsense..
pt=kfish2evg367mixx
$evigene/scripts/rnaseq/asmrna_altreclass.pl -trclass $pt.trclass -mapsensetab publicset/$pt.gmapsense.tab \
  -out $pt.asenrealt2.idtab -maxalt=19 -altrenum -debug
#altreclass: nin/out=499047/499047, ndrop/keep=110027/389020, nrenum=299194, ngenediff=34890/149901
  ^^ more are kept here, why?
cat kfish2evg367mixx.asenrealt2.idtab | grep antisens | wc -l =    7552
	drop  = 4502; keep = 3050
	dropmain = 340 << these need to be checked w/ blastp for bad antisense calls.
	Funhe2Exx3m006258t29    Fungr1EG3m002911t1      Funhe2Exx3m006258       29      dropmain        747,78%,complete,antisense      99/100  oldid:Funhe2Exx3m006258t1
	Funhe2Exx3m010025t7     Fungr1EG3m004981t1      Funhe2Exx3m010025       7       dropmain        573,76%,complete,antisense      99/100  oldid:Funhe2Exx3m010025t1
		^^ antisense yes, but revaa/tr may be best model, should recalc fwd aa,cds, not drop
		
cat kfish2evg367mixx.asenrealt2.idtab | grep antisens | grep -v drop | head
Funhe2Exx3m001262t1     Funhe2Emap3m000720t1    Funhe2Exx3m001262       1       main    1502,89%,complete,antisense     99/99   .
Funhe2Exx3m006254t1     Funhe2Emap3m004077t1    Funhe2Exx3m006254       1       main    748,21%,complete-utrbad,antisense       99/100  .
Funhe2Exx3m006448t3     Fungr1EG3m003035t1      Funhe2Exx3m006448       3       main    734,62%,complete,antisense      98/97   oldid:Funhe2Exx3m006448t1
Funhe2Exx3m006448t4     Funhe2Eq7m051374t1      Funhe2Exx3m006448       4       altmid  700,63%,complete,antisense      98/97   .
		^^ 448t3,4 : was t1,t2; old t3,4 became new t1,t2 due to antisense.
Funhe2Exx3m007789t3     Funhe2Emap3m005117t1    Funhe2Exx3m007789       3       main    662,51%,complete-utrpoor,antisense      99/100  oldid:Funhe2Exx3m007789t1

grep Funhe2Exx3m006448t kfish2evg367mixx.asenrealt2.idtab
** FIXME: main reset here, need change class output..
## human:UniRef50_Q9H0H3 Ectoderm-neural cortex protein 2, 589aa matches new 448t2, t1 may be bogus partial5.
Funhe2Exx3m006448t1     Funhe2E6bm006723t2      Funhe2Exx3m006448       1       altmid  654,84%,partial5        98/91/Funhe2Eq7m051374t1        oldid:Funhe2Exx3m006448t3
Funhe2Exx3m006448t2     Funhe2E6bm006723t3      Funhe2Exx3m006448       2       althi1  589,82%,complete        98/100/Funhe2Eq7m051374t1       .
Funhe2Exx3m006448t3     Fungr1EG3m003035t1      Funhe2Exx3m006448       3       main    734,62%,complete,antisense      98/97   oldid:Funhe2Exx3m006448t1
Funhe2Exx3m006448t4     Funhe2Eq7m051374t1      Funhe2Exx3m006448       4       altmid  700,63%,complete,antisense      98/97   .


eg. keep partial?  t25 part3 here is bigger than t6+, 99/92 says 8% differs from old main
.. map view says likely good form, almost complete
Funhe2Exx3m001853t25    Funhe2Eq7m065601t1      Funhe2Exx3m001853       25      1264,95%,partial3       99/92   dropalthi1      oldid:Funhe2Exx3m001853t10
..
Funhe2Exx3m001853t1     Funhe2E6bm001704t3      Funhe2Exx3m001853       1       1513,95%,complete       99/87   althi   oldid:Funhe2Exx3m001853t3
Funhe2Exx3m001853t2     Funhe2E6bm001704t6      Funhe2Exx3m001853       2       1286,89%,complete       99/90   althi   oldid:Funhe2Exx3m001853t6
Funhe2Exx3m001853t3     Funhe2E6bm001704t37     Funhe2Exx3m001853       3       1279,90%,complete       99/90   althi1  oldid:Funhe2Exx3m001853t48
Funhe2Exx3m001853t4     Funhe2E6bm001704t40     Funhe2Exx3m001853       4       1277,91%,complete       99/91   althi1  oldid:Funhe2Exx3m001853t41
Funhe2Exx3m001853t5     Funhe2E6bm001704t35     Funhe2Exx3m001853       5       1270,89%,complete       99/92   althi1  oldid:Funhe2Exx3m001853t42
Funhe2Exx3m001853t6     Funhe2E6bm001704t14     Funhe2Exx3m001853       6       1116,83%,complete       99/94   althi1  oldid:Funhe2Exx3m001853t25
Funhe2Exx3m001853t7     Funhe2E6bm001704t28     Funhe2Exx3m001853       7       1116,90%,complete       99/94   althi1  oldid:Funhe2Exx3m001853t37
Funhe2Exx3m001853t8     Funhe2E6bm001704t38     Funhe2Exx3m001853       8       1116,87%,complete       99/94   althi1  oldid:Funhe2Exx3m001853t45
Funhe2Exx3m001853t9     Funhe2E6bm001704t29     Funhe2Exx3m001853       9       1108,93%,complete       99/94   althi1  oldid:Funhe2Exx3m001853t31
Funhe2Exx3m001853t10    Funhe2E6bm001704t25     Funhe2Exx3m001853       10      997,79%,complete        99/94   althi1  oldid:Funhe2Exx3m001853t26
...
Funhe2Exx3m001853t20    Funhe2E6bm001704t1      Funhe2Exx3m001853       20      1292,99%,partial        99/91   dropmaina2      oldid:Funhe2Exx3m001853t1
 > case of dropped main .. is this right? partial but 2nd longest. prob same as Funhe2Exx3m001853t6 = 1286aa complete + 6aa; maybe bad orf call


/bio/bio-grid/kfish2/rnas/kf2evgr/trevg367mixx
poor alt removal addition to trclassing 
may need to do after making pubids ?
separate altpoor classifier?

#......
gmapz/altpoordrop.info

## redo / add to trclass to remove many poor/partial alts,
   use pubids + trclass info: aaqual, class, maybe also p-id/aln score?
   sort per gene by aasize, qual?
   per gene remove (to other file) cases where alt.aasize << main.aasize, alt class=frag, alt.aaqual = partial
   especially for common case of % id/aln ~ 99/100 (ie subset alt)
 
# redo, add $pia to pubidx 
grep okay *.trclass | cat - publicset/kfish2evg367mixx.pubids | perl -ne\
'if(/\tokay/) { @v=split"\t"; ($od,$cl,$pia,$aq)=@v[0,2,4,5]; \
($piav)= $pia =~ m,^(\d+/\d+),; $oda{$od}="$aq\t$piav\t$cl"; } \
elsif(/^Fun/) { ($pd,$od,$td,$ti)=split; $oda=$oda{$od}||0; s/$/\t$oda/; print; }' |\
sort -k3,3 -k5,5nr -k1,1 > publicset/kfish2evg367mixx.pubidx2

# alt-reclassifier: require sorted by gene at least, also by aa-size? or do that here?
# .. keep orig pubid or renum alts? need option
cat publicset/kfish2evg367mixx.pubidx2 | env altrenum=0 perl altreclass.pl

=cut
