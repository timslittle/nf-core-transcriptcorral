#!/usr/bin/perl
# evgs/genes/trclass2pubset.pl 

=item trclass2pubset.pl

 revision of parts of  evgmrna2tsa2.pl
 following genes/altbest2pubset.pl (pubset from annot GFF)
 roughly same as MAIN_pubsetonly() of evgmrna2tsa2
 
 updated parts:
   trclass2maintab() =  public class/id table from trclass, okayset of tr2aacds
   
=item usage (maybe)

 set pt=pig4321ew; 
 $evigene/scripts/genes/trclass2pubset.pl -debug -outbase pig4ewc \
  -keepdrop $pt.keepdrop.tab -keepoldids $pt.keepids -names $pt.names \
  -mrna tmpfiles/$pt.mrna -trclass $pt.trclass

 .. this is fairly messy set of opts, need all?
 .. expect/require? inputs of both -mrna (with okayset oids) and -trclass table
    name.mrna should also have name.aa, cds seq files associated
    -mrna could be from okayset/, or from vec trimset result
 .. want output of all publicset/, pubids, ann.txt
 .. -outbase should be option, default to publicset/$trname
 .. fix -keepdrop <> -keepoldids ?
     -keepoldids is unusual option, but want for now
     -keepdrop table may be common input, for culled and dropped okayset

=item UPD20apr merge usage
  
  trclass2pubset.pl -rna ixos9fok1st.okreor.mrna,ixos9fmn.uvcut.mrna,ixos9fmn.uvcut.ncrna,ixos9mpub4f.mrna,ixos9ncode.rna \
    -pubids xxx.pubids -keepdrop xxx.keepdrop -names xxx.names -log -debug

  -rna list is ordered precedence, first-in rna IDs are only output, ie. okreor.mrna replaces > uvcut.mrna > pub.mrna
  -keepdrop table for this needs work, OK/DROP/SKIP acts now, 
    .. is this per-seq file, or for all seqset?
    .. SKIP means ignore ID per-seq file, DROP ID means drop from all seqset (?? check)
  -pubids table in instead of -trclass, reusing orig pubids ; -keepids not needed here
  for rna/cds/aa input files, .gz allowed but leave off -rna list
  
=item UPD1911 altreclass

  add exontab chain analysis from genes/blasttrset2exons2.pl
  ** need updated rnaseq/asmrna_altreclass3c.pl
  asmrna_altreclass3c.pl -debug -nodrops -noclasscut=60 -altrenum \
    -trclass $trclass -pubids $pubids -exontable $exontab -out $realtids
 
  my $APPaltreclass= findevigeneapp("rnaseq/asmrna_altreclass.pl",1); # upd 201405 ;  NODIE

=cut

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/.."); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
use strict;
use warnings;
use Getopt::Long;
use cdna_evigenesub; 
use evigene_pubsets; # now has some of below subs
use evgpubsetsum; # debug; move up

use constant VERSION => '2020.03.15'; # '2019.11.05'; 
# '2018.06.18'; #  from evgmrna2tsa2.pl & altbest2pubset
use constant UPD1911 => 1;
use constant UPD20mar => 1;

our $EVIGENES= "$FindBin::Bin/.."; # ="$FindBin::Bin";  # WRONG place ; this is in evgs/genes/ subdir ..
our $DEBUG= $ENV{debug}|| 0;
our $EGAPP='class2pub';  # or okay2pub ?
our $IDPREFIX= $ENV{idprefix} || 'EVGm'; # *! no default, else -idprefix param replaced by ORGANISM !*
our $ORGANISM= $ENV{ORGANISM} || ""; # no default
my  $NEWIDpre= $ENV{idtag}||'n'; #??
my  $SRC= $ENV{source}||"evm"; # fixme: for gff output, source.class col2 src_(main/alt/other);
our $SHOWDROPS=1; #default on? only for mainalt.tab outputs
our $SIZESORT=1; 
our $FAHDRisValid=1; # export from pubset.pm for annots from mrna hdr
our ($pubid_format,$altid_format,$GBPROID);
our $pubidnum_start=0;
our $preserveOldIds=undef; # change to ok pattern, IDPREOK
our $SORT_PUBSET= 0;
our $DATE;  
our $RNATYPES; # our $RNATYPES='mRNA|ncRNA'; # in evigene_pubsets.pm; transcribed_RNA mrnatypes exontypes 
our ($KEEPDROPX,%keepdropxtab); # evigene_pubsets.pm
my $KEEPOKPUBID= 0; # set depending on preserveOldIds vs trclass ids
my $USEOKCULL= 0; # UPD20mar: related to KEEPOKPUBID
my $ONLY_PUBIDS= 0; my $CLASSbySOURCE= 0;
my $CULLOTHERS= 0; # how to treat unclassed input genes, unused now

# main opt:
my $AAMIN_NOCLASS=$ENV{aaminnoclass}||60; # asmrna_altreclass -noclasscut=$AAMIN_NOCLASS; drop noclass (needs user opt), but rescued by aaref
my $NOALTDROPS=0; # turn on if have keepdrop table? add user opt? nodrop turns the to culls

my ($trclass,$output,$logfile,$keepdropin,$trexontab);  
my ($trpath,$trname, $sradatah, %settings);
my ($oname, @genes, @seqset, $genenames,$pubidin);
my $doREALT= ($ENV{norealt})?0:1;

#UPD20mar: evgdirs: trimset=>vecset, add ncrnaset, tmpsets: okayset1st, prepubset, reorset, publicset_oldor
# sra2genes: inputset tmpfiles dropset okayset ncrnaset publicset submitset \
#       aaeval geneval okayset1st reorset vecset prepubset publicset_orold 

our @evgdirs = qw(okayset dropset inputset trimset tmpfiles erasefiles publicset submitset);
our (@okayset,@dropset,@inputset,@tmpfiles,@erasefiles,@publicset,@submitset); # tidyup file sets
our (%pubids,%pubidinfo,%puballoids); ## cdna_evigenesub/evigene_pubsets globals 
my $tidyup=1; my $NCPU=1;

my @ARGSAVE=@ARGV;
my $optok= GetOptions(
  "trclass|class|input=s", \$trclass, # many gff inputs : use perl input pipe while(<>) ?
  "mrna|cdna|rna|sequences=s", \@seqset, # many inputs? 
  "oname|outbase|output=s", \$oname,  # FIXME: idpre option  overwritten by spppref
  "idprefix=s", \$IDPREFIX,  
  "tagidnew=s", \$NEWIDpre,   
  "source=s", \$SRC,   
  "names=s", \$genenames,   
  "exontab=s", \$trexontab,   # UPD1911
  "pubids|pubidin=s", \$pubidin,   # UPDapr; pubidin maybe== preserveOldIds
  "keepoldids|preserveOldIds=s", \$preserveOldIds,  # allow -idpre=XXX -keepold w/o 2nd param
  "keepdropin=s", \$keepdropin,   # other opt key? reclassin?
  "logfile:s", \$logfile,
  "dropshow!", \$SHOWDROPS,  
  "sizesort!", \$SIZESORT,  # default:on ; may be required to get proper ids w/ preserveOld
  "pubsortseq!", \$SORT_PUBSET,  # default:off?
  "CULLOTHERS!", \$CULLOTHERS,  # or otherclass=cull|skip|remainder name?
  "NOALTDROPS!", \$NOALTDROPS,  #  -noaltdrops == culls instead
  "realt!", \$doREALT,  # UPD1911-norealt instead of ENV{norealt}
  "onlypubids|ONLY_PUBIDS!", \$ONLY_PUBIDS,  #CLASSbySOURCE also : same opt?
  #gff: "classbysource!", \$CLASSbySOURCE,  # default:off?
  # "dryrun|n!", \$dryrun, 
  "NCPU=s", \$NCPU, # not used here, leave in opts   
  "tidyup!", \$tidyup, # default on ?
  "debug!", \$DEBUG, 
  );

#? push @genes, grep(/gff/, @ARGV); # remainder gff ? trclass here?

die "EvidentialGene trclass2pubset -trclass myspp.trclass [-idprefix $IDPREFIX ] 
  makes tables of public ids and main-alt linkage, from results of tr2aacds
  opts: -idprefix Thecc1EG  -mrna myspp.mrna  -names mrna.names 
     -keepdrop keep_drop_ids.table -preserveOldIds=old.pubids
     -nosizesort -[no]pubsortseq -debug
  version ", VERSION, "\n"
  unless($optok and ($trclass or $pubidin)); # (@genes or @seqset) 

$NEWIDpre="" unless($NEWIDpre=~/^[A-Za-z]/); # allow same old id?

my $IDPREOK=$IDPREFIX;
# if(defined $preserveOldIds) {
#   if($preserveOldIds and $preserveOldIds=~/^[a-zA-Z]/) { $IDPREOK= $preserveOldIds; }
#   else { $preserveOldIds=$IDPREOK=$IDPREFIX; }
# } else {
#   $preserveOldIds=0;
# }  

# data globals
my($mainindex,$ntr)= (0) x 9;
my(%annotes, %main, %mainsize,%alt,%oclass,%altsize,%altdrops,%didaltdrops,%balt,%drop);  # globals of trclass2maintab
#unused here: my(%aaqual, %piad, %newmain, %altdefer,%maindrops);

my($nkeepdrop, $keepdroph, )= (0);
my($nmapqual, $mapqualh, $alntabh)= (0); #? globals from map align.tab, if exists
# %dropids, %dropid not used yet; change to culldropkeep table: id, actclass, actinfo, oid?, other..
# orig dropids,dropid is  short list of gff.source tags, e.g. 'cull' in source 'zf17evgcull'

openloggit($logfile,$trclass||$pubidin);  
loggit(1, "EvidentialGene trclass2pubset VERSION",VERSION);
loggit(1, "$0  @ARGSAVE");

MAIN();

sub MAIN
{
  ## default output to publicset/oname
  ## maybe default input seqs from okayset/ ; see mrna2tsa:get_evgtrset()
  @seqset= map{ split(/[,;\| ]/,$_) } @seqset;
  unless($oname){ $oname=$trclass || $pubidin || $seqset[0]; $oname=~s/\.\w+$//; }

  $DATE=`date '+%Y%m%d'`; chomp($DATE);
  loggit(0, "BEGIN $0  input=",$oname,"date=",$DATE);
	#? do_settings("restore",$trclass);  
	
  my($upokgenes, $upokids, $upstatus)= (0) x 9; 
  my($pubids,$maltab)= map{ "$oname.$_" } qw(pubids mainalt.tab);
  $pubids= $pubidin if($pubidin and -f $pubidin); # UPD20apr
  my($aaseq,$cdsseq,$addseq)=("","",undef); # set if paired w/ mrna
  my $generef= (@genes > 0)? \@genes : undef;
  my $mrnaseq= (@seqset > 0)? $seqset[0] : ""; # many?? all input .mrna .cds and/or .aa
  
  #UPD20Mar: get_evgtrset() should fill in @seqset instead, for extra mrna,aa,cds .. for makePublicset()
  # USEOKCULL for get_evgtrset() ?? $addseq
  ($mrnaseq,$trpath,$trname,$aaseq,$cdsseq,$addseq) #UPD20j get_evgtrset add aa,cds
    = get_evgtrset($trclass,$mrnaseq,"publicset"); # cdnaseq == mrnaseq here from getmRNA

  loggit(0, "get_evgtrset=",$mrnaseq,$trpath,$trname);  
  loggit(LOG_DIE, "Missing mrna",$mrnaseq) unless(($mrnaseq and -s $mrnaseq) or $ONLY_PUBIDS);
  
  ($nkeepdrop,$keepdroph)= ($keepdropin) ? readKeepDrop($keepdropin) : (0); # fill  global %keepdrop
  $NOALTDROPS=1 if($nkeepdrop);

  if($genenames and -f $genenames) { 
    #UPD1912 load here names for annot in makePubIdtabFromTrclass(), via gene_annot_brief()
    #cdna_evigenesub.pm sets globals: %genenames %genedbxref %genenamepct %namedgenes %cddnames
    my($nnamed,$namin)= parse_genenames($genenames,NAME_NOEDITS);
    loggit(0, "names $genenames n=$nnamed\n"); 
  }

  ## FIXME19: pubset wants trname.align.tab , gmap now makes trname-chrname.align.tab ..

  # (my $mapqual=$trclass)=~ s/\.\w+$/.align.tab/;
  # if(! -f $mapqual and -d "geneval") { $mapqual="geneval/$mapqual"; }
  my($mok,$mapqual)= getMapqual($trname);
  if($mok) {
    ($nmapqual, $mapqualh, $alntabh)= readAlignTab($mapqual);
  }
  
  # Separate steps here, call w/ input pubid tab to make pubgff, pubseq sets
  unless( -s $pubids ) {
    ($upokids)= makePubIdtabFromTrclass($pubids,$maltab,$trclass); # read attr from mrnaseq?
	  $upstatus++ if($upokids>0);
	  #see makePubIdtab:  altreclass_block($trclass,$pubids); # reclassifies alts in pubids table, drops some
  } 
  
  if($ONLY_PUBIDS) { # want this opt?
    loggit(0, "done -onlypubids makePubIdtab n=",$upokids);
    #x exit(0); # go to tidyup below
  } else {
    ## change this mrnaseq to array= grep 'mrna|cdna|xxx' @seqset
    #UPD20Mar: makePublicset() should use @seqset instead, for extra mrna,aa,cds
    ($upokgenes)= makePublicset($pubids, $oname, $mrnaseq,$aaseq, $cdsseq, $addseq, $generef); 
    $upstatus++ if($upokgenes>0);
  }
  
	#? do_settings("log|save",$trclass||$mrnaseq,); # or after last call?  ("log|restore|save");
  if( $tidyup and $upstatus > 0 ) {  
    tidyupFileset("publicset",@publicset);  
  	tidyupFileset("submitset",@submitset);  #? after or before tsa_tbl2asn? need path in fileset
    tidyupFileset("tmpfiles",@tmpfiles);  
    my @rmlist;
    foreach my $fn (@erasefiles) { if(-f $fn) { unlink($fn); push @rmlist,$fn; } }  
    if(@rmlist) { my $nrm=@rmlist; my $rml=join " ", splice(@rmlist,0,5); loggit(0,"tidyup erase: n=$nrm, $rml .."); } 
  }

  loggit(0, "DONE at date=",`date`);
  loggit(0, "======================================"); # log may append runs
}
  
# ---------------------------

## readKeepDrop moved to evigene_pubsets.pm
# my($nkeepdrop,$keepdrop_idhash)= readKeepDrop($keepdropin) if($keepdropin); # fill  global %keepdrop
use constant { kdDROP => -999, kdCULL => -1, kdOK => 1, kdMUSTKEEP => 2, kdOther => 0 };  

sub getMapqual {
  my($trname)= @_;
  ## FIXME19: pubset wants trname.align.tab , gmap now makes trname-chrname.align.tab ..
  # (my $mapqual=$trclass)=~ s/\.\w+$/.align.tab/;
  my $mapqual = $settings{'mapqual'} || finddata("$trname*.align.tab") || finddata("geneval/$trname*.align.tab") ; 
  my $ok= ($mapqual and -s $mapqual)?1:0;
  return ($ok,$mapqual);
}

# sub makePubIdtabFromGeneGFF
# add sub makePubIdtabFromTrclass

sub makePubIdtabFromTrclass
{
  my($pubids,$maltab,$trclass)= @_; # ,$genes
  
  # FIXME: 18.03/04 got 300 dup pubids in 300k now, mixup recent..
  
	# ($pubid_format,$altid_format)= make_IDPREFIX(); # default abbrev of $organism now
    #^^ wait for input gff IDs as existing id?
  loggit(0, "makePubIdtabFromTrclass($pubids,$maltab)");

  my($ok, $ndone, $nmiss, $nkeepids)= (0) x 9; 
  my($upstatus,$upfiles,$uptemp,$upokids) = (0) x 9;  

  # genenames/old pubids need trclass oid > old.pubid map
  # fill cdna_evigenesub global %pubids w/ keepids? 
  # ** ^^ problem later? global takes precedence over this pubid calc
  # -- need 2nd pubids{oid} for genenames, other intermediate annots w/ pre-publicset, post-trasm ids
  # read_pubids() sets %pubidinfo w/ $alloids : use that?
  # my($oid,$gid,$alti,$class,$aqual,$pia,$notes,$alloids)= split"\t", $pubidinfo{$pubid};

  # annotab2tblinfo?? pubidin vs global %pubids may be problem, pubid from this sub used as final? public id
  # .. need here only for genenames annot via gene_annot_brief() > parse_evgheader() ; fix those
  #  cdna_evigenesub: global %pubids, only parse_evgheader() uses that
  #  evigene_pubsets: gene_annot_brief() uses parse_evgheader,  
  #    our (%pubids,%pubidinfo,%puboidc,%puballoids); ## cdna_evigenesub globals %pubids,%pubidinfo
  #    putPubAnnot(),make_annot_from(),make_annotab(),make_pubseq() use %pubids
  #    read_pubids() and read_annotab2() replace %pubids
  
  if($preserveOldIds and  -f $preserveOldIds and open(F, $preserveOldIds)) {
    while(<F>){ next if(/^\W/);  
      my($pd,$oid,$gid,$alti,$tclass)=split;
      # next if($oid !~ /^\w/ or ($tclass and $tclass =~ /^(drop|cull)/)); #UPD20mar ?? problems skipping cull
      next if($oid !~ /^\w/); #UPD20mar ?? problems skipping cull
      #BAD# $pubids{$oid}=$pd; 
      $puballoids{$oid}{$pd}++; # use this?
      $nkeepids++;
      } close(F);
  }

  # add trclass reading parts of evgmrna2tsa2.pl MAIN_pubsetonly() 
  #	$upstatus=0; 	
  #  ($upstatus, $upfiles, $uptemp, $upokids)  
  #    = update_mrna_fileset($trpath, $cdnaseq, $trimflag, $trimids, @trimfiles); #? leave out here?
  
  ($upokids)= readTrclass($trclass,"publicset"); # was trclass2maintab()
  $KEEPOKPUBID= ($nkeepids >= 0.66 * $upokids) ? 1 : 0; # preserve pubids thru altreclass if enough
  
  reclassGenes(); # link alts to mains via %main,%alt,%balt id hashes
  
  ($pubid_format,$altid_format,$GBPROID)= make_IDPREFIX() unless($pubid_format); # default abbrev of $organism now
  
  $nmiss=$upokids;
  my($outpubidh,$outh);
  $ok= open(OPD,'>',$pubids); $outpubidh=*OPD;
  if($ok) { $ok= open(OMA,'>',$maltab);  $outh=*OMA; }
  if($ok) {
    ($ndone,$nmiss)= putPubidTab($outh,$outpubidh); # $mainindex,$ntr
    close($outh); close($outpubidh);
    push @publicset, $pubids, $maltab;
  } else {
    loggit(1, "#ERR: cant putPubidTab($pubids)");
  }
  loggit(0, "done makePubIdtab($pubids,$maltab) nin=$upokids, npub=$ndone, nmiss=$nmiss");
  
	## OPTION, dont want as is when have keepdrop cull table, BUT do want its added mapqual
	if($doREALT){ altreclass_block($trclass,$pubids); } # reclassifies alts in pubids table, drops some

  return($ndone,$upokids);
}

=item altreclass or not

  no altreclass: $ENV{norealt} = 1 
  altreclass w/o dropalt:
    -noclasscut= 1
    -maxaltsame= 99999
    -nodrops : added
    >> not yet ENV/opts
    DROPTINYALT = 0
    DROPSHORT_AAPART = 0
    DROPSHORT_AAFULL = 0
    DROPSHORT_ANTISENSE = 0
    
cat publicset/pigevg4wc.pubids.old | cut -f5 | sed 's/a2//; s/hi1/hi/; s/midfrag/frag/; s/cull[12]/cull/; s/cullalt.*
/cullalt/; s/maybe//; ' | sort | uniq -c | head -40

before altreclass
376671 althi
  4264 altmid
  3840 altfrag
 32299 cullalt
 33160 cullmain
 61469 cullnoclass
     3 cullparthi
 30668 main
 11970 noclass
  3515 parthi

after altreclass : is cutting cull tag
357756 althi
  2984 altfrag
 13687 althim       | reclass mains
  4371 altmid
   971 dropaltfrag  | want
 36619 dropalthi    | these 
   373 dropaltmid   | altdrops ?
 64259 main    : but has culls
 73321 noclass : but ann.txt has cull tags
  3518 parthi

=cut

=item UPD1911 altreclass

  add exontab chain analysis from /genes/blasttrset2exons2.pl
  ** need updated rnaseq/asmrna_altreclass3c.pl
  my $APPaltreclass= findevigeneapp("rnaseq/asmrna_altreclass.pl",1); # upd 201405 ;  NODIE
  
if [ -s $trnrname.exontab ]; then
  $evigenes/rnaseq/asmrna_altreclass3c.pl -debug -nodrops -noclasscut=60 -altrenum \
  -trclass $trname.trclass -pubids publicset/$trname.pubids  \
  -exontable $trnrname.exontab  -out publicset/$trname.pubids.realt \
  > publicset/$trname.altreclass3c.log 2>&1

  if [ -s  publicset/$trname.pubids.realt ]; then
    mv publicset/$trname.pubids publicset/$trname.pubids.old
    mv publicset/$trname.pubids.realt publicset/$trname.pubids;  
  fi
 

=cut

sub altreclass_block {
  my($trclass,$pubids)= @_;

  ## FIX 201405: insert here asmrna_altreclass.pl, .realt_pubids to become new .pubids ..
  my $APPaltreclass= findevigeneapp("rnaseq/asmrna_altreclass.pl",1); # upd 201405 ;  NODIE
  if($APPaltreclass =~ /MISSING/) { $APPaltreclass=""; }
  
  if($doREALT and $APPaltreclass and -x $APPaltreclass) { # was not $ENV{norealt} 
    # ADD OPTIONS: -noclass=$MINAA == drop noclass short things if no other good qualities
    #  -MAXALTSAME=n == drop excessive number of althi class that appear same by aa size/aa cluster, 
    #  asmrna_altreclass to be run or not?   other opts?
    
    my $realtids="$pubids.realt";
    my $aopts="-debug -noclasscut=$AAMIN_NOCLASS";
    $aopts.=" -nodrops" if($NOALTDROPS);

    # class calls: pubids now override trclass class, from culls, etc. 
    # need flag to altreclass    
    $aopts.= " -trustpubids"; # or $ENV{trustpubids}=1; # 
    
    ## aopts:  "valueoutput:s", \$valoutput, # option to dump table of classing scores
    $aopts.=" -valueoutput" if($DEBUG > 1);
    
    (my $trname=$trclass)=~ s/\.trclass.*//;
    my($mok,$mapqual)= getMapqual($trname);
    if($mok) { $aopts.=" -mapqual $mapqual"; }
    
    if(UPD1911) {
      # ** need updated rnaseq/asmrna_altreclass3c.pl >> rename asmrna_altreclass4v.pl ?
      # asmrna_altreclass3c.pl -debug -nodrops -noclasscut=60 -altrenum \
      #   -exontable $exontab -trclass $trclass -pubids $pubids -out $realtids
      if($trexontab and -s $trexontab) { 
        # NOT HERE? let caller set properl altreclass.pl .. BUT need be ok for -exontable
        # .. maybe test date of APPaltreclass >= 2019;  grep -x UPD1908 $app ?  grep readExonTab $app ?   
        #   perl  rnaseq/asmrna_altreclass3c.pl -BOBexontab xxx
        #     Unknown option: bobexontab  << use this test?
        my $updaltreclass= findevigeneapp("rnaseq/asmrna_altreclass3c.pl",1);  
        $APPaltreclass=$updaltreclass if($updaltreclass !~ /MISSING/ and -x $updaltreclass); # NOT best way,
              
        # $trnrname.exontab == tmpfiles/trname_nrcds.exontab, from blasttrset2exons2.pl
        $aopts.=" -exontable $trexontab"; 
      }
    }
    
    ##UPD20mar: -altrenum may not be wanted, to preserve okayset/pubids = publicset/pubids when have annots w/ okayset pubids **
    # use constant KEEPOKPUBID => 0; # * need option, use only after making okayset/pubids, for publicset/pubids
    $aopts .= ($KEEPOKPUBID) ? " -noaltrenum" : " -altrenum";
    
    my $cmd="$APPaltreclass $aopts -trclass $trclass -pubids $pubids -out $realtids > $realtids.log 2>&1";
    my $runerr= runcmd($cmd);
    
    if(-s $realtids) {
      rename($pubids,"$pubids.old"); rename($realtids,$pubids); 
      push @publicset, "$pubids.old", "$realtids.log"; # or tmpfiles ?
    } else {
      loggit(LOG_WARN,"ERR: failed update $realtids by $APPaltreclass");  # 2014 add notice...
      push @publicset, $realtids, "$realtids.log"; # or tmpfiles ?
    }
    my $realtstat=`tail -n1 $realtids.log`; chomp($realtstat); $realtstat ||= "$APPaltreclass updated $pubids";
    loggit(0, $realtstat, "log=$realtids.log"); # log tail has stats 

  } else {
    loggit(LOG_WARN,"PLEASE ALSO RUN publicset fixup:\n",  # 2014 add notice...
      "evigene/scripts/rnaseq/asmrna_altreclass.pl -trclass $trclass -altrenum -out -debug");
  } 

}



sub makePublicset {
  my($pubids, $outname, $cdnaseq,$aaseq,$cdsseq, $addseqh, $genegffset)= @_;
  # new($pubids, $oname, $mrnaseq,$aaseq,$cdsseq, $genegff);
  # from evgmrna2tsa sub MAIN_pubsetonly()
  # $cdnaseq => \@mrnaseq ?
  
  my $skiptrrun=0;
  #global opt now# my $genenames=""; #?? look for
  our @publicset; #global in pubsets.pm
	my $pubdir="publicset";

  loggit(0,"makePublicset($pubids,$outname)"); # LOG_NOTE,LOG_DEBUG // if($DEBUG);

  my(@mrna,@aa,@cds);
  #UPD20mar: makePublicset()  @seqset may have same as cdnaseq,aaseq,cdsseq
if(UPD20mar) {
  #?? allow .seq.gz here
  my @addseq= (ref $addseqh) ? @$addseqh : ();
  @mrna= grep /\.(mrna|ncrna|rna|cdna|tr)$/, ($cdnaseq,@seqset,@addseq);  
  @aa  = grep /\.(aa|pep)$/, ($aaseq,@seqset,@addseq);  
  @cds = grep /\.(cds)$/,  ($cdsseq,@seqset,@addseq);

} else {  
  ## change this cdnaseq/mrnaseq to array= grep 'mrna|cdna|xxx' @seqset
  @mrna= grep /\.(mrna|ncrna|rna|cdna|tr)$/, @seqset;
  @aa  = grep /\.(aa|pep)$/, ($aaseq,@seqset); # if empty but @mrna try s/.mrna/.aa/ ?
  @cds = grep /\.(cds)$/, ($cdsseq,@seqset);
  if(@mrna) { $cdnaseq= $mrna[0]; } elsif($cdnaseq) {  @mrna=($cdnaseq); } 
}
  
  # BUG fixed, was $f =~ s/\.\w+/.aa/; == bad for name.okay.mrna 
  if(@mrna and not @aa) { @aa= map{ (my $f=$_) =~ s/\.\w+$/.aa/; $f; } @mrna; }
  if(@mrna and not @cds) { @cds= map{ (my $f=$_) =~ s/\.\w+$/.cds/; $f; } @mrna; }
  
  if($genenames and -f $genenames) { }
  elsif( -f "$outname.names") { $genenames= "$outname.names"; } #? need -option
  else { my $pn= makename($pubids,'.names'); $genenames=$pn if(-f $pn); }
  
	my($npubid, $pubidh)= read_pubids($pubids, $cdnaseq); # if($pubids); 
	  # NOT cdnaseq here it isnt publicset/pubid file?? but make_annotab uses cdnaseq annots/IDs
    # return($nred, \%pubids);

	my($annotab, $ntra1, $ncdsa1)=(0) x 9;	
  use constant OLDANN => 0;
  if(OLDANN) {
    # ** FIXME: make_annotab input  @mrna  AND/OR @genegff annot input
    # UPD: make_annotab() uses pubset.pm global pubids hash info : ADD gff mRNA annots for names, other
    # uses pack global %pubids from read_pubids
	  ($annotab, $ntra1, $ncdsa1) 
		  = make_annotab($cdnaseq, $genenames, $skiptrrun, $outname); # add main/alt pub ids, other geneinfo 
  
  } else {
    $annotab =  makename($outname,".ann.txt"); 
    if(not -f $annotab and -d $pubdir) {
       my($pubd,$ft)= getFileset($pubdir,'.ann.txt');  
       $annotab=$ft if($ft);
    }
    if(ref($genegffset)) { ## not in trclass2pubset
      # DUPID fix: may need keepdrop{ oid } handling here to skip dropped dups by oid..
      ($ntra1, $ncdsa1)= make_annot_from("gff", $annotab, $genegffset,  $genenames);
    } elsif(@mrna) {
      ($ntra1, $ncdsa1)= make_annot_from("mrna", $annotab, \@mrna,  $genenames);
    }
  }
  #NO, make_annot does: push @publicset, $annotab; # if($ntra1);
   
  $FAHDRisValid=1; # export from pubset.pm for annots from mrna hdr
  my($nann,$annothash)= read_annotab2($annotab);
  $FAHDRisValid=2 if($npubid and $nann > 0); ## upd1809; FLAG 2 == sensible merge only missing vals

  my(@sumtables)= geneset_summary($pubidh, $annothash);
  if(@sumtables) { 
    my $gs = makename($outname,".genesum.txt");  
    open(TO,'>',$gs); for my $tab (@sumtables) { print TO $tab,"\n"; } close(TO); 
    push @publicset, $gs; 
    }

  # may already have .aa, .cds, .mrna from prior steps. is this ok?
  # add input -seqset params ??  our $RNATYPES='mRNA|ncRNA';

	my($pubmrna,$npm,$minfo) = make_pubseq(\@mrna,'mRNA',$annothash, "$outname.mrna");
	my($pubaa,$npa,$ainfo) 	= make_pubseq(\@aa,'protein',$annothash, "$outname.aa"); # makename($cdnaseq,'.aa')
	my($pubcds,$npc,$cinfo)	= make_pubseq(\@cds,'CDS',$annothash, "$outname.cds"); # makename($cdnaseq,'.cds')
  
  #make_pub does: push @publicset, $pubmrna,$pubaa,$pubcds; # if($npm);
  #UPD: now make_pub puts *.checktab to tmpfiles if no err
  
      # DUPID fix: may need keepdrop{ oid } handling here to skip dropped dups by oid..
  my($pubgff,$ngenes,$ginfo)= make_pubgff($genegffset,$SRC,$annothash, $outname);
  
  my $nout= $npm || $npa || $npc || 1;
  loggit(0,"publicset: ",$pubmrna,$minfo,$pubaa,$ainfo,$pubcds,$cinfo,$pubgff,$ginfo,$annotab); 
  return($nout);
}



sub reclassGenes {
  #?? keepdrop/cull handle here?
  my @amain= grep { not $main{$_} } sort keys %alt; # dropmain here now only for SHOWDROPS !
  foreach my $am (@amain) { 
    my $md= $balt{$am} || $am;  
  	if(!$main{$md} and $drop{$md}) { my $md1= $balt{$md}||""; if($md1 and $main{$md1}) { $md=$md1; } }
  	
    if($main{$md}) { my @at=keys %{$alt{$am}}; map{ $alt{$md}{$_}=$alt{$am}{$_}; $balt{$_}=$md; } @at; } 
    elsif($md) { 
     $main{$am}="NOMAIN";  # FIXME: get rid of these by finding alt main
     } 
  }
  
  foreach my $td (keys %balt) {
    my $md= $balt{$td} || $td; 
    #?need: unless($md) { (undef,undef,$md)= evigene_idparts($td); }
    $main{$md}="NOMAIN" unless($main{$md});# FIXME: get rid of these by finding alt main
  }
}


# evgmrna2tsa2.pl get_evgtrset(); move to package?
sub get_evgtrset { 
  my($trclass,$cdnaseq,$pubdir)= @_;
  my($trpath,$trname,$aaseq,$cdsseq)=("") x 9; # my $sradata=undef; my $nsra=0;
	my $notokay=0;
  my @addseq=();
  
  if($cdnaseq) { 
    $notokay=1; # dont look in okayset/? look in $pubdir now?
    $trclass= makename($cdnaseq,".trclass") unless($trclass); 
  }

  if($trclass) {
    my $trpname= makename($trclass,"","trclass"); 
    if($trpname =~ m,/,) { ($trpath,$trname)= $trpname =~ m,(.*)/([^/]+)$,; } # was BAD?
    else { $trname=$trpname; }
    $trpath ||= '.';  $trname=~s/\.gz//;
    
    my $okpath= ($notokay) ? $trpath :"$trpath/okayset"; # should I check if curdir has okayset files?

    #OBSOLETE# my($mrnaOfUtrorf,$nutrorf)= getmRNA_utrorf($okpath,$trname);

    $USEOKCULL=1; #?? if( -f "$okpath/$trname.cull.mrna"); #??
    
    my $mrnaopts="addreor,"; # UPD20ja was ADDutrorf=1, now opts= noutrorf, addreor|addrestrand  for okayset/na.okreor.seqs
    $mrnaopts .= "addcull," if($USEOKCULL); # okayset => publicset reclass
    
    #oo ($cdnaseq)= getmRNA($okpath,$trname,$pubdir,LOG_WARN,$mrnaopts) if(!$cdnaseq and -d $okpath);
    #ob ($cdnaseq,$aaseq,$cdsseq)= getmRNA($okpath,$trname,$pubdir,LOG_WARN,$mrnaopts) if(!$cdnaseq and -d $okpath);
    my @ret_mrna_aa_cds= (-d $okpath) ? getmRNA($okpath,$trname,$pubdir,LOG_WARN,$mrnaopts) : ();

#UPD20may: bug got okreor, need okay.mrna -- fix in getmRNA()
#egr: EvidentialGene trclass2pubset VERSION 2020.03.15
#egr: get_evgtrset= ./okayset/plDZQM.okreor.mrna . plDZQM <<< == cdnaseq from splice(@ret_mrna_aa_cds,0,3); 

    if(@ret_mrna_aa_cds > 0) { 
      ($cdnaseq,$aaseq,$cdsseq)= splice(@ret_mrna_aa_cds,0,3); 
      @addseq= @ret_mrna_aa_cds; # culls,etc are any extras; need all 6 w/ addcull, BUT dont dup 0..2
    }
  }
  
  return($cdnaseq,$trpath,$trname, $aaseq,$cdsseq, \@addseq); #UPD20j add aa,cds; drop ,$sradatah; \@seqset ??
}


sub readTrclass # from sub trclass2maintab
{
  my($trclass,$pubdir, $okids)=@_;
  my $ntr=0;  my $nerr=0;
  my $mainindex= $pubidnum_start;
  my $hasokids= ($okids and ref($okids) and scalar(%$okids))?1:0;
  
  #... not here?? no output written in this sub
  my $maintab = makename($trclass,".mainalt.tab","trclass");  # > $pt.mainalt.tab
  my $pubidtab= makename($trclass,".pubids","trclass");   
  # if(not -f $pubidtab and $pubdir and -d $pubdir) {
  #   my($pubd,$ft);
  #   ($pubd,$ft)= getFileset($pubdir,'pubids',$pubd);  $pubidtab=$ft if($ft);  
  #   ($pubd,$ft)= getFileset($pubdir,'mainalt.tab',$pubd);  $maintab=$ft if($ft);  
  # }
  # return($maintab,$pubidtab,$mainindex,$ntr) if( -s $maintab and -s $pubidtab);# or dryrun ..

  #globals: (%annotes, %main,%mainsize,%alt,%oclass,%altsize,%altdrops,%didaltdrops,%balt,%drop);  # globals of trclass2maintab
  #unused here: %newmain, %altdefer, %maindrops [unused],  %aaqual, %piad [in annots]    from mrna2tsa
  
  my($ok,$inh,$outh,$outpubidh);
  ($ok,$inh)= openRead($trclass);
  unless($ok) { loggit(1,"ERR: parse $trclass TO $maintab"); return; }

  while(<$inh>) {
    next unless(/^\w/); chomp;
    my($td,$ok,$cl,$md,$piad,$aq,$notes)=split "\t";
    unless($cl and $md and $aq) { $nerr++; loggit(1,"ERR: trclass inline:$_"); next; } # what?? report?

		## maybeok needs handling, drop? or keep/mark/cull ??
		## should revise asm dupfliter to avoid these maybes, includes 'refbest', others some good/bad
		## now all are from exoneq flag, all 'althi1|part|frag' classes; 5702 refbest(dups?) of 13059
		## ** maybe keep all maybes ; found some refbest/good being dropped here.. but maybes are large subset
		## should filter by hoscore if avail
		
		## UPD1807 keepdropin should override trclass drop/okay : not working below ***
    # UPD1911: new noncode-maybe tag .. change to noncode here? mainnc althinc noclassnc; also noncode == pflag & 32
    # if($cl=~m/nc$/) { } 
		
		my $dropit=0;
		my $clorig=$cl;
		if($ok eq 'maybeok') { 
		   if($notes=~/refbest|refgood/) { $ok="okay"; $cl.="maybe"; } 
		   else { $ok="okay"; $cl.="maybe"; } #??
		}
		
    if($ok ne 'okay') { $cl=$ok.$cl; $drop{$td}=$cl; $dropit=1;  } #  OPTION: include drops?
    elsif($hasokids and not $okids->{$td}) { $dropit=1; $cl='dropid'.$cl;  $drop{$td}=$cl; } 
      #^^ replace okids with keepdroph 
 
    my($kdxact, $kdxid)=(0,0);
    ## if($KEEPDROPX) { ($kdxact, $kdxid)= check_keepdrop('seqid',"$td,$alloids"); }
    ## if($kdxact =~ /cull|drop/) { $evgclass=$kdxact.$evgclass unless($evgclass =~ /^(cull|drop)/); }
     
    if($keepdropin and (my $kdv= $keepdroph->{$td})) {      
      #? my $oid= $tblinfo->{oid}; # for gff
      # for my $od (split",",$oid){ if(my $kdvo= $keepdroph->{$od}){ $kdv=$kdvo; last; } }
      
      my($kdi,$kact)= split" ",$kdv; #no need kdi?
      if($kdi < 0) {   # drop or cull, now also kdSKIP == -3 .. not for trclass? but for merges of same ID
        if($kdi <= kdDROP){ 
          $dropit=1; $cl='dropid'.$cl;  $drop{$td}=$cl;
          $kdxact="drop"; $kdxid=$td;
          # $mclass=cDROP; $ok=0; $ndrop++; #  $kact=~/drop/ else drop anything else?
          }
        elsif($kdi == kdCULL){ 
          # $cl= ($isalt or $mainid)?"cullalt":"cull"; # variants?
          # $cl .= "cull"; # suffix? prefix
          $cl = "cull$cl"; # suffix? prefix
          $kdxact="cull"; $kdxid=$td;
          # $mclass= cCULL; 
        } 
      } elsif($kdi>0) { 
        if($dropit) { $dropit=0; delete $drop{$td}; $cl=$clorig; $ok="okay"; }
        # $mclass= cMAIN unless($mclass == cALT || $mclass == cNEWLOC);  
      }  
      
    }
      
    ## all piad entries should have pd and sense/asense fields, placeholder where needed '.' ?
    #o.my($pi,$pa,$asense,$pd)=split"/",$piad; # asense before pd ** NEED To revise asm dupfilter to clarify
    my($pi,$pa,$asense,$pd)= split"/",$piad; # asense before pd ** NEED To revise asm dupfilter to clarify
    map{ $_ ||=""; } ($pi,$pa,$asense,$pd);
    # piad == 100/100/self1 << self1 misused here. means mainid == td ?
    
    if($asense =~ /sense/) { $pd="" unless($pd =~ /^\w/); }
    elsif($asense =~ /^\w/) { $pd=$asense; $asense=""; }
    # if($pd eq "self1") { $md=$td; } # is this right? only for noclass, where td == md already above
    $md=$pd if($pd =~ /^\w/ and not $pd=~/self/);
   	$cl=~s/a2$//;  #? dont need 

		## FIXME here?  new altmap class from eqgene hides main class now; all mains are reclassed altmap, 1 should be left as main.
 		if($dropit and not $SHOWDROPS) {
 			next unless($cl =~ /main|noclass/); # this is enough, skip dropalts
 		}
 		
 		$ntr++;
 		my $aasize= ($aq =~ m/(\d+).(\d*)/) ? "$1.$2" : 0; # size.pCDS for best sort
 		
 		my $mapq= ($nmapqual) ? $mapqualh->{$td} : "";
 		# my $maptab= ($nmapqual) ? $alntabh->{$td} : "";
    $notes.= ",chrmap:$mapq" if($mapq);
 		
 		my($gannot, $tblinfo)= gene_annot_brief($td,"mRNA",$notes);
 		
		$annotes{$td}{aaqual}= $aq;
		$annotes{$td}{piad}= $piad;
	  $annotes{$td}{notes}= $gannot;
		$annotes{$td}{oid}= $tblinfo->{oid};
  	$mainsize{$td}= $altsize{$td}= $aasize; #? dont need both?  aasize{$id}= $aasize;
 		
    if($cl =~ /main|noclass/)  # include cullmain|cullnoclass 
    { 
    	$main{$td}=$cl; $balt{$td}=$td;
    } else { 
    	$alt{$md}{$td}= $cl; $balt{$td}=$md; 
    	# NOMAIN fix: always add dropmain here 
    }  
  }

	if($SHOWDROPS) {
		# drops from dropset/*.aa headers for perfect_dups, perfect_frags info not in trclass ..
		my $ndr=0;
  	my($dset,$droptr)= getFileset("$trpath/dropset",'drop.tr|drop.aa|drop.cdna');  # dangit need fixed tr/cdna/fa suffix here
  	my($ok,$hin)= ($droptr) ? openRead($droptr) : (0,0);
  	if($ok) { 

  	  ## FIXME.160911: MISSING evgclass=,drop, listings in mainalt.tab .. should include all dropset/*.aa hdr ids
  		# FIXME2: spurious drops getting to pubid table via this $md == main of drop td, but may also be a drop!
  		# below via @amain from alt{md} != this main?
  		## new problem: md here may well not be in mainlist or linked there .. will then leave
  		## these dropped main/alt out of mainalt.tab .. maybe right, these are items of no value?
  		
  		while(<$hin>) { if(/^>(\S+)/) {  
  		  my $td=$1;
  			my($cl,$ok1,$md)= m/evgclass=(\w+),(\w+),match:([^\s;,]+)/;
  			# others:  evgclass=noclass,drop;  and no evgclass=
  		 	if($cl and $md and not $balt{$td}) { 
    			$ok1="drop"; # ensure no bad cases
  			  my $mmd= $balt{$md}||$md; # locate new main
  		 	  $drop{$td}="$ok1.$cl,$md"; 
  		 	  $altdrops{$mmd}{$td}= $ok1.$cl; $ndr++; 
  		  }
  		  # BUG160911: altdrops{md} miss when md reclass to not-main;  not: $balt{$td}=$md; 
  		} 
  	} close($hin); }
		loggit(0,"readTrclass: dropset has $ndr alts"); 
	}
	
  return($ntr); 
  # reclassGenes(); # called after this reader     
}



=item gannot pubid annots

  * moved to evigene_pubsets gene_annot_brief
  reproduce this annotation (trclass2mainalt.pl) in notes hash
    aaref:5767,dapsim:Dapsim1EVm000004t1,chrmap:100a,98i,25555l,33x,29sc:324762-356329:+,pflag:0
    
  my @ANKEY_MRNA= qw(aalen  cov pid nexon clen namealn Dbxref scoresum);
  leave out special annots of altbest, alttr, altid/mainid/newlocusid

  gannot moved to evigene_pubsets gene_annot_brief($id,@mrnarow)

=cut

sub putPubidTab {
  my($outh, $outpubidh)= @_; # globals?
 
  #globals: (%annotes, %main,%mainsize,%alt,%oclass,%altsize,%altdrops,%didaltdrops,%balt,%drop);  # globals of trclass2maintab
 
  ## headers 
  #mainalt: originalID     MainClass  Alternates
  #pubid  : Public_mRNA_ID originalID      PublicGeneID    AltNum  Class   AAqual  pIdAln  Notes

  print $outh '#'.join("\t",qw(originalID MainClass Alternates))."\n";
  print $outpubidh '#'.join("\t",qw(Public_mRNA_ID originalID PublicGeneID AltNum Class AAqual pIdAln Notes Oids))."\n"
    if($outpubidh);
  #upd1802: pubid add Oids col
  
  my %doneid=();
  my @mainlist;
  # %oclass adds problems == main or alt
  
  #* BUG bad %main id set?? %balt has all; this doesnt fix missing doneid
  # problem in %alt{md} **
  foreach my $ad (keys %balt) {
    unless( $drop{$ad} ) {
    my $md= $balt{$ad} || $ad;  
    unless( $main{$md} ) { $main{$md}="NOMAIN"; } # NOMAIN ?
    unless( $md eq $ad or $alt{$md}{$ad} ) { 
      my $omc= $oclass{$ad} || "altmiss"; # $cull= $omc if($omc); #?
      $alt{$md}{$ad}= $omc; } # culls here ? also mains only for ad ne md?
      }
  }

  if($SIZESORT) { # or NOT IDSORT ?
  @mainlist= sort{ ($mainsize{$b}||0) <=> ($mainsize{$a}||0) or $a cmp $b } keys %main;
  } else { # IDsort, only special cases, like merge old/update, keeping old order for most
  @mainlist= sort{  $a cmp $b } keys %main;
  }

  #upd1708: preserveOldIds option ..  keep input gene ids as best feasible, including altnum?
  #  need all alts per main, check evg gene ids for 1 x many, also need check across loci for 2+ new loci w/ old geneids
  #  where cant preserve geneid for all of locus, do what? new geneid, or modify old?
  #  UPD: want/need main = t1, cant preserve old alt nums if alt classing changes ..

  my ($nnewids,$newidh)= ($preserveOldIds) ? preserveOldIds(\@mainlist, \%alt, \%drop, \%altsize) : (0, undef);
      #,$newpubidh ret
      
  foreach my $md (@mainlist) { 

    my @ad= sort{ $alt{$md}{$a} cmp $alt{$md}{$b}
      or $altsize{$b} <=> $altsize{$a} or $a cmp $b } keys %{$alt{$md}}; 

    my $culls={}; #? ($CULLXEQ) ? cullExonEq($md,\@ad,\%alt,\%notes) : {}; #?? here
    my $ad= join",",map{ "$_/".$alt{$md}{$_} } @ad; 
    my $mc= $main{$md}; 
    ## oclass cullalt mix up main class
    # if($mc=~/other/){ $mc=$oclass{$md}; }
    
    # FIXME: change $mc eq "noclass" to "main" if alts exist, and vversa
    # FIXME2: preserve $mc eq "cull" .. use %culls{id} ?
    my($cull)= ($mc =~ m/^cull/) ? "cull":""; # this is problem? replace w/ omc=oclass{md}
    my $omc= $oclass{$md} || ""; # $cull= $omc if($omc); #?
    
    # FIXME cullalt in @ad dont count to main/noclass
    my $altcount= @ad;
    $altcount= scalar( grep{ my $c= $oclass{$_}||""; $c !~ m/cull/ } @ad );
    
    if($mc =~ /^NOMAIN/) { 
      # do below
    } elsif($altcount>0 and $mc !~ /^main/) {
      $mc =~ s/^\w+/main/;
      $main{$md}= $mc= $cull.$mc; 
    } elsif( $altcount==0 and $mc !~ /^noclass/) {
      $mc =~ s/^\w+/noclass/;
      $main{$md}= $mc= $cull.$mc; 
    }
    
    if($SHOWDROPS) {
    	my @add= sort{$altdrops{$md}{$a} cmp $altdrops{$md}{$b}} keys %{$altdrops{$md}};  
    	map{ $didaltdrops{$_}=1 } @add; # later dump not didaltdrops
    	my $add=join",",map{ "$_/".$altdrops{$md}{$_} } @add; 
    	$ad .= ",$add" if ($add);
    }
    
    my $mcout=$mc;
    if($omc) { $mcout=$omc.$mc if($omc=~/^(cull|drop)/); } #? both?
    print $outh join("\t",$md,$mcout,$ad)."\n";  # mainalt.tab
    
    if($outpubidh) { # should be required ??
      my $ialt= 0; my $needmain=0;
      my($cla,$aaq,$pida,$nots,$oids);
      $cla= $mcout;  # cla=$main{$td}=$cl;  ; changed above?
      # ** $cla=  $main{$md}||"main"; # need main == cull now
      # $aaq= $aaqual{$md}||"noaa";
      # $pida=$piad{$md}||0;  
      # $nots=$notes{$md}||"nonote";  
      $aaq= $annotes{$md}{aaqual}||"noaa";
      $pida=$annotes{$md}{piad}||0;  
      $nots=$annotes{$md}{notes}||"nonote";  
      $oids=$annotes{$md}{oid}||"noid";  

      if($mc eq "NOMAIN") { $cla=$cull . (($altcount>0)?"main":"noclass"); } ## needs to change, to main? to noclass?
      elsif($mc =~ /^alt/) { } # is this were nomain show up? or @ad?

      my @sad = @ad;      
      if(1) {  # $SIZESORT .. want to use altbest sort order, from tNN altnum
        @sad= sort{ $altsize{$b} <=> $altsize{$a} or $a cmp $b } @ad;      
      }
      
      if($drop{$md} or $doneid{$md}) { $needmain=1; }
      else {
      	$mainindex++; $needmain=0; # BUG: move below drop{}
      	my($pubmrnaid,$pubgeneid,$pubti)= make_pubid($md, $mainindex, ++$ialt, $newidh);
      	print $outpubidh join("\t",$pubmrnaid,$md,$pubgeneid,$pubti,$cla,$aaq,$pida,$nots,$oids)."\n";  #n
      	$ntr++; $doneid{$md}++;
      	}
      
      foreach my $ad (@sad) {
        unless($drop{$ad} or $doneid{$ad}) {
        $cla=  $alt{$md}{$ad}||"nocl"; 
        $aaq= $annotes{$ad}{aaqual}||"noaa";
        $pida=$annotes{$ad}{piad}||0;  
        $nots=$annotes{$ad}{notes}||"nonote";  
        $oids=$annotes{$ad}{oid}||"noid";  
        # $cull= $culls->{$ad}||""; # $cla == "cullalt..";
      	if($needmain) { 
      	  $mainindex++; $needmain=0; 
      	  if($cla=~/^alt/){ $cla=(@sad>1)?"main":"noclass"; } 
      	}  
        my($altmrnaid,$altgeneid,$alti)= make_pubid($ad, $mainindex, ++$ialt, $newidh);
        print $outpubidh join("\t",$altmrnaid,$ad,$altgeneid,$alti,$cla,$aaq,$pida,$nots,$oids)."\n"; 
        $ntr++; $doneid{$ad}++;
        }
      }
    }
  }
  
  # double check not-done .. missing what? LOTS .. missing from %main ??
  my($notdone)=(0);
  foreach my $ad (keys %balt) {
    unless( $drop{$ad} or $doneid{$ad} ) {
    my $md= $balt{$ad} || $ad; 
    $altdrops{$md}{$ad}++;  $notdone++;
    }
  }
  
  if($SHOWDROPS) { # UPD.160911
    my @mains= sort keys %altdrops;
    for my $md (@mains) {
      my @misdrop= grep{ not $didaltdrops{$_} } keys %{$altdrops{$md}};
      if(@misdrop) {
        @misdrop= sort{$altdrops{$md}{$a} cmp $altdrops{$md}{$b}} @misdrop;  
        map{ $didaltdrops{$_}=1 } @misdrop; # later dump not didaltdrops
        my $misdrop= join",",map{ "$_/".$altdrops{$md}{$_} } @misdrop; 
        my $mc= $main{$md}||"NOMAINd";  # fixme: oclass{$md}
        print $outh join("\t",$md,$mc,$misdrop)."\n";  # mainalt.tab
      }
    }
  }
  return($ntr,$notdone);  
  # close($inh); close($outh);
  #? return($maintab,$pubidtab,$mainindex,$ntr);  # return main,alt id hashes ....
}


sub preserveOldIds {  # for trclass/okayseq set with input -preserve oldpubid.tab 
  my($mainlist, $altsOfMain, $drop, $altsize)= @_;
  my(%gids, %gnums, %gdone, %newids, %newpubids, $gprefix);
  my $nids=0;
  #above# 
  my $NEWIDpre=''; # 'n';
  
  ## prefix, gnum problem: have mixed gprefices .. diff gnum set for each
  ## FIXME  gnum only for $gprefat =~ /$IDPREOK/
  # change for evgmrna2tsa, read table of pubid,oid from old gene set to preserve
  # @$mainlist is of oids, not old pubids;

  my $GNEXTNUM=0;
  my (%pod,%idparts);
  #UPD20mar: have read preserveOldIds already for global %pubids{oid}= pubid .. do this same time?
  if( -f $preserveOldIds and open(F, $preserveOldIds)) {
    while(<F>){ next if(/^\W/);
      #o: my($pd,$oid)=split; next unless($oid =~ /\w/);
      my($pd,$oid,$gid,$alti,$tclass)=split; # likely is pubids table, maybe not
      # next if($oid !~ /^\w/ or ($tclass and $tclass =~ /^(drop|cull)/)); #UPD20mar ?? problems skipping cull
      next if($oid !~ /^\w/); #UPD20mar ?? problems skipping cull
      $pod{$oid}=$pd; 
      my($gd,$gpre,$gnum,$ti)=(0,0,0,0);
      if($pd =~ m/^(\w+[A-Za-su-z])(\d\d+)t(\d+)$/) { # basic evg id form
        ($gpre,$gnum,$ti)=($1,$2,$3);
        $gd= $gpre.$gnum;
        $gnums{$gnum}++;
        $idparts{$pd}= [$gpre,$gnum,$ti]; # unless($idparts{$pd}); # dups? shouldnt be
        $gprefix= $gpre unless($gprefix);
        # $gids{$gnum}{$ti}= $pd;
        # $newids{$oid}= $pd;
      }      
    } close(F);
  }
  return unless(%pod);
  
  # foreach my $md (@$mainlist) {
  #   my @okd = grep{ not($drop->{$_}) } ($md, keys %{$altsOfMain->{$md}} );
  #   for my $id (@okd) {
  #     my($gpre,$gnum,$gd,$ti,$isplit)= evigene_idparts($id);
  #     next unless($gnum and $gpre =~  m/$IDPREOK/);
  #     $gnums{$gnum}++; # $gids{$gd}{$ti}= $id;
  #     $gprefix= $gpre unless($gprefix);
  #     }
  # }
  
  my($glast)= sort{ $b <=> $a } keys %gnums;
  $GNEXTNUM= 9 + $glast;
  
  my $idformat= $NEWIDpre . $gprefix . '%06d'; #?  use pubid_format of make_IDPREFIX  ?
  my(%havepubg);

  foreach my $md (@$mainlist) {
  
    my @ad= sort{ $a cmp $b } keys %{$altsOfMain->{$md}}; 
    my @okd = grep{ not($drop->{$_}) } ($md,@ad);
  
    my ($mnum,$nd,$timax)=(0,0,0);
    for my $oid (@okd) {
      my $pubid= $pod{$oid} or next;
      # check all alts for diff old pubgene id?
      my($gpre,$gnum,$ti)= @{$idparts{$pubid}};
      $mnum= $gnum unless($mnum);
      
      $timax = $ti if($ti > $timax); #UPD20mar
      #UPD20mar: Dont change orig ti ? but need to record use, fit new alts b/n, after ?
      # .. but not changing pubid,ti unless newpubids has it
      # use constant KEEPIDALT=>1;
      # if(KEEPIDALT) {  
      #   if($ti > $timax) { $timax = $ti; } 
      # } else { # old      
      #   if($ti > $timax) { $timax = $ti; } else { $ti= ++$timax; }
      # }  
    
      if($gnum ne $mnum or $newpubids{$pubid} ){
        while(1) {
          $pubid= $gprefix . $mnum . sprintf( $altid_format, $ti); # FIXME ti clash
          if( $newpubids{$pubid} ) { $ti= ++$timax; } else { last; }
        }
      }
      
      $newids{$oid}= $pubid; $nd++; $nids++;
      $newpubids{$pubid}= $oid;
    }
    
    if($nd == 0 and @okd > 0) {
      $mnum= ++$GNEXTNUM;
    }
    if($nd < @okd) { # finish, some or all new per locus / mnum 
      for my $oid (@okd) {
        next if($newids{$oid});
        my($pubid,$ti);
        do {
          $ti= ++$timax; # keep list unused ti and reuse?
          $pubid= $gprefix . $mnum . sprintf( $altid_format, $ti);
        } while( $newpubids{$pubid} );
        $newids{$oid}= $pubid; $nd++; $nids++;
        $newpubids{$pubid}= $oid;
      }
    }
  } # mainlist
  
  return($nids,\%newids,\%newpubids);
} # sub preserveOldIds

sub preserveOldIds_FOR_GFF {
  my($mainlist, $altsOfMain, $drop, $altsize)= @_;
  my(%gids, %gnums, %gdone, %newids, %newpubids, $gprefix);
  my $nids=0;
  #above# my $NEWIDpre='n';
  
  ## prefix, gnum problem: have mixed gprefices .. diff gnum set for each
  ## FIXME  gnum only for $gprefat =~ /$IDPREOK/

  my $GNEXTNUM=0;
  foreach my $md (@$mainlist) {
    my @okd = grep{ not($drop->{$_}) } ($md, keys %{$altsOfMain->{$md}} );
    for my $id (@okd) {
      #o.my($gpre,$gnum,$gd,$ti,$isplit)= evigene_idparts($id);
      my($gpre,$gnum,$ti,$isplit,$gdup)= evgIdParts($id);
      next unless($gnum and $gpre =~  m/$IDPREOK/);
      $gnums{$gnum}++; # $gids{$gd}{$ti}= $id;
      $gprefix= $gpre unless($gprefix);
      }
  }
  my($glast)= sort{ $b <=> $a } keys %gnums;
  $GNEXTNUM= 9 + $glast;
  
  my $idformat= $NEWIDpre . $gprefix . '%06d'; #?  use pubid_format of make_IDPREFIX  ?
  my(%havepubg);

  foreach my $md (@$mainlist) {
  
    # my @ad= sort{ $altsOfMain->{$md}{$a} cmp $altsOfMain->{$md}{$b} # class sort
    #  or $altsize->{$b} <=> $altsize->{$a} or $a cmp $b } keys %{$altsOfMain->{$md}}; 
    my @ad= sort{ $a cmp $b } keys %{$altsOfMain->{$md}}; 
    my @okd = grep{ not($drop->{$_}) } ($md,@ad);
  
    %gnums=(); %gids=();
    my $gnumfirst=0; my $gprefat=0;
    for my $id (@okd) {
      #o.my($gpre,$gnum,$gd,$ti,$isplit)= evigene_idparts($id);
      my($gpre,$gnum,$ti,$isplit,$gdup)= evgIdParts($id);
      next unless($gnum); # for _G2,n.. new loci
      $gnums{$gnum}++; $gids{$gnum}{$ti}= $id;
      $gnumfirst=$gnum unless($gnumfirst or $gdone{$gnum}); 
      $gprefat= $gpre unless($gprefat); #?? problems? yes
    }
    
    my @gnums= sort{ $gnums{$b}<=>$gnums{$a}  or $a <=> $b } keys %gnums;
    unless($gnumfirst) { # pick most or first in (main) <<
      for(my $i=0; $i<=$#gnums; $i++) {
        unless($gdone{$gnums[$i]}) { $gnumfirst=$gnums[$i]; last; }
      }
    }
    # usage assumption: -keepids=IdPreA -idpre IdPreA, so inval gprefat will go to ++GNEXTNUM 
    unless($gprefat and $gprefat =~ /$IDPREOK/){ $gprefat=$IDPREFIX;  $gnumfirst=0; }
    $idformat= $NEWIDpre . $gprefat . '%06d'; # maybe new for each main id??
    my($pubgn,$haveit)=(0,0);
    do { 
      unless($gnumfirst) { $gnumfirst= ++$GNEXTNUM; } #??
      $pubgn = sprintf( $idformat, $gnumfirst);
      $haveit= $havepubg{$pubgn} || 0;
      ## fixme dup gene ids now, check with t1
      unless($haveit) { $haveit=1 if($newpubids{$pubgn . 't1'}); }
      $gnumfirst=0 if($haveit);
    } while ($haveit);
    $havepubg{$pubgn}=1; 

    # try2
    my %tidone=(); my $timax= @okd; #? not same using orig ti > @okd
    for my $id (@okd) {
      #o.my($gpre,$gnum,$gd,$ti,$isplit)= evigene_idparts($id);
      my($gpre,$gnum,$ti,$isplit,$gdup)= evgIdParts($id);
      my $tinew=0;
      if($gnum eq $gnumfirst and not $tidone{$ti}) { $tinew=$ti; }
      if($tinew==0) { do { $tinew++; } while( $tidone{$tinew} ); } # can put lowqual alts at ti top
      $tidone{$tinew}++; 
      $gdone{$gnumfirst}++;
      #above# my $pubgn = sprintf( $idformat, $gnumfirst); 
      my $pubti = sprintf( $altid_format, $tinew); # same ti or new? cant use same ti w/o checks
      my $pubid=  $pubgn . $pubti;
      $newids{$id}= $pubid;
      $newpubids{$pubid}= $id;
      $nids++;
    }
  } # mainlist
  
  return($nids,\%newids,\%newpubids);
} # sub preserveOldIds




1;

__END__
