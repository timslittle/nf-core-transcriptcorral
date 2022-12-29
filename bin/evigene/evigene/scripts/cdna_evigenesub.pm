# cdna_evigenesub.pm

# package cdna_evigenesub;
package main;

use strict;
use warnings;
use FindBin;
use File::Basename qw(basename dirname fileparse);

use constant TRAA2CDS_2018 => 1; # prot/traa2cds.pl usage updated for utrorf handling
use constant UPD1908 => 1; 
use constant UPD1912 => 1;
use constant UPD20UORF => 1; # 2020mar26: update utrorf annots thruout evigene

#? require Exporter;
# our @ISA = qw (Exporter);
# our @EXPORT = qw (xxxx);

use vars qw ( $EVIGENES $EGAPP $GDB_PREFIX $DEBUG $dryrun
  $MIN_NAMEIDENT $MIN_IDLIKE $EVGLOGH %EVG_SEQSUFFIX
  %genenames %genedbxref %genenamepct %namedgenes %cddnames $HAVE_gnames
  %pubids %pubidinfo %puballoids $APPtraa2cds
  $AAQUALF $AAQUALH $BAD_GAPS
  );

## add globals for getAaQual here? should be in main caller?
use constant { kAAQUAL_MAX => 3, kAAQUAL_MIN => -3, kAAQUAL_NONE => 0, }; #  aaqualscore() range
our ($AAQUALF,$AAQUALH) = ("",undef,undef);
#NOT YET# $AASIZEH 
#  $AAQUALH->{$id}="$alen,$pctcds,$acv,$aqual1";  $acv == numeric score of aqual1
#  our %AAQUALS = (); %AASIZES = (); # global hash, keep in sync w/ file name

our $DEBUG=0;
our $dryrun=0; ## $DRYRUN ?
our $EVIGENES="$FindBin::Bin"; #??
our $EGAPP='mrna2tsa'; # FIXME?
our $EGLOG='egr';
our $GDB_PREFIX='gnl|Evigene|';  #see below;  use IDPREFIX? No, ID has
our $APPtraa2cds= undef; #findevigeneapp("prot/traa2cds.pl"); # move to cdna_evigenesub for get_mRNA

## name stuff 
our $MIN_NAMEIDENT = 35;  # min similar% for protein evid align/equivalence, for naming; JCVI uses 35%  
our $MIN_IDLIKE   = 15;  # low for now; 15..20 seems right;add Note with loqualname
our $BAD_GAPS= 25;  # % gaps in AA

# UPD1911 required seq file suffices, see makename(), getOkFileset() getFileset() getmRNA()
# makename.BAD  $insuf ||= 'aa|blast|cdna|mrna|cds|tr|trclass|tbl|fasta|faa|fsa|fa';  # need insuf: tr|fasta|fa
our %EVG_SEQSUFFIX = ( 'tr' => 'tr|cdna', 'aa' => 'aa|pep', 'cds' => 'cds', 
                       'mrna' => 'mrna', 'ncrna' => 'ncrna'); 
                      

use constant { LOG_NOTE => 0, LOG_WARN => 1, LOG_DIE => -1, LOG_DEBUG => 2, };
our $EVGLOGH= undef; # renamed EVGLOGH from logh; package local? now exported

sub loggit{ 
	# my $dowarn=shift; my $s= join(' ',@_); # dang warn @_ empty join for 1st call here, from where ??
	my($dowarn,@msg)= @_; return unless($dowarn or @msg);
	return if($dowarn== LOG_DEBUG and not $DEBUG);
  #my $s= join(' ',@msg); ## Use of uninitialized value $msg[3] in join
  my $s= join ' ', map{ defined($_) ? $_ : '.' } @msg;
  chomp($s); $s="FATAL $s" if($dowarn == LOG_DIE);
  if($EVGLOGH){ print $EVGLOGH "#$EGLOG: $s\n"; } elsif($dowarn>0||$DEBUG){ warn "#$EGLOG: $s\n"; }
  if($dowarn == LOG_DIE) { die "#$EGLOG: $s\n" ; }
}

sub openloggit {
  my($logfile,$trname)= @_;
  if(not $logfile and defined $logfile) { # use output name
    $logfile= $trname || $EGLOG;
    $logfile= makename($logfile,".$EGAPP.log");  # need program suffix??
  }
  if($logfile) { 
    open($EVGLOGH, '>>', $logfile) or die $logfile; 
    ## EVGLOGH should have immediate write, $| = 1 ??, like STDERR
    my $lastio = select(STDOUT);
    select($EVGLOGH); $| = 1;
    select($lastio); 
  } 
}

sub openRead { # add to cdna_evigenesub.pm
  my($fna)= @_; my($ok,$hin)= (0,undef);
  if($fna and -f $fna) { 
    $ok= ($fna =~ /\.gz$/) ? open($hin,"gunzip -c $fna|") : open($hin,$fna); }
  loggit(1,"ERR: openRead $fna") unless($ok);
  return ($ok,$hin);
}

## note these are in cdna_protein also; need more package local privacy.
sub _min1 { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max1 { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }



sub parse_genenames
{
  my($genenames, $noEdits)= @_;
  $noEdits ||=0; 
  # returns in globals: (%genenames,%genenamepct,%genedbxref,%namedgenes,%cddnames) 
  
  if($HAVE_gnames and $genenames and $HAVE_gnames eq $genenames) { # calling here multiple times, diff subs
    my $ngot= scalar(keys %genenames);
    my($gwc)= `wc -l $genenames`; my($nin) = $gwc =~ m/(\d+)/;
    return($ngot,$nin) if($ngot > 1 and $ngot >= 0.5 * $nin);
  }
  
  my($ngot,$nin)=(0,0);
  %genenames=%genenamepct=%genedbxref=%namedgenes=%cddnames=();
  return($ngot,$nin) unless($genenames and -f $genenames);
  
  ## FIXME2: ** use uniq names vs ERR: too short to keep valid tr, e.g. 
  # er2g: ERR: too short:183 #LitvaEG0018688t4     oid=litovavel2k35Loc15824t1     len=183 name=CDD: clpS, ATP-dependent Clp protease ada..-like    
  # grep  'CDD: clpS,' *.names = 1 only = litovavel2k35Loc15824t1
  # FIXME: need better reader; 2+ rows/id; pick best .. format may change..
  # names.tab ==  id, name, pctalign, refid, repid  : now
  #  trid1  C-ets-2 protein, putative       89%,103/116,197 RefID:UniRef50_E0VFI2   RepID:E0VFI2_PEDHC
  #  trid2  DBH-like monooxygenase protein 1        73%,445/613,516 RefID:UniRef50_Q6UVY6   RepID:MOXD1_HUMAN

  #UPD1912: score refbest/good from all ref blast scores , as per trclass
  my(%aablast, %aablastref, %bscore);
  my($ok,$inh)= openRead($genenames); 
  unless($ok) { loggit(1,"ERR: parse_genenames reading $genenames"); return; }
  
  while(<$inh>) { 
    next unless(/^\w/ and /\t/);
    chomp; $nin++;
    my($id,$name,$pctalign,$refid,$repid)=split"\t"; # may have only id, name

    my $xtra=""; 
    ($name,$xtra)=split";",$name,2 unless($noEdits); #??? xtra may be valid some have ;
    $name =~ s/\s+$//;

    ## BUG in data: missing pctalign fields ; dont know why.
    ## output of evigene/scripts/prot/namegenes.pl 
    ## whitefly1vel5k45Loc9888t2       Synaptojanin-1-like protein     RefID:UniRef50_B4E1Z3   UniProt:B4E1Z3_HUMAN
    if($refid and not defined $repid and $pctalign and $pctalign =~ /^\D/) {
    	$repid=$refid; $refid=$pctalign; $pctalign="";
    }
    $pctalign||=""; $refid||=""; $repid||=""; # missing ok
    
    # FIXME: 2 names/id maybe: CDD: xxx and gene xxx; keep both in ann.txt ? and pctalign?
    ## pctalign == 100%,450/450,446 : pct,naln/nref,ntrg
    # pctalign format expected: 72%,3270/4555,3282 ;  may be '72' or '72%' only
    
    # my ($pcta,$aln)= $pctalign =~ m=^(\d+)%,(\d+)\b=; # partial match ok
    # unless($pcta) { ($pcta)= $pctalign=~m/(\d+)/; $aln= $pcta; }
    
    my($pcta,$aln)= ($pctalign =~ m=^(\d+)%,(\d+)\b=) ? ($1,$2)    
                  : ($pctalign=~m/(\d+)/) ? ($1,$1) : (0,0);
    my $haspct= ($pcta > 0)?1:0; # zero is missing value; $pctalign =~/^\d/ and

    if(!$noEdits and $haspct and $pcta < $MIN_NAMEIDENT) { # use MIN_IDLIKE : name-like ? got some '0%' align
      ## bad: Uncharacterized protein-like ; Nuclease HARBI1, putative-like; Protein-like
      ## *should* leave these name changes to nameclean()
      if($pcta >= $MIN_IDLIKE) { 
       unless($name =~ /\blike|^Uncharacterized/) {
        $name =~ s/, putative//; 
        unless( $name =~ s/\s+protein$/-like protein/ ) { $name .= '-like'; } ## fixme: 'xxxx protein -like'
        }
      } else { next; } ## caller should decide? should we preserve for ann.txt table ? as Unchar ?
    }

        # DBXREF_RECODE fix for ncbi's  dbx: restrictions : evgmrna2tsa2.pl:putTblFsa()
        # FIXME999: more problems w/ gene.names table having odd/local DBprefix:ID
        #   .. fix where? should have here list of valid NCBI/ISxxx db prefixes. from where?
    
    ## fixme: CDD:206692,cd04107,RefID:UniRef50_Q9NX57,UniProt:RAB20_HUMAN, 
    ## drop  RefID:; drop? cd04107
    $refid =~ s/^RefID://;  $repid =~ s/^RepID://;  ## RepID: also 
    ## ?? try here add right DbPrefix: ? Uniprot/Uniref easy, others a mess.
    map { if(/:/) { } # asis
    	elsif(/^UniRef/i or /^[A-Z0-9]+_[A-Z][A-Z]+$/) { $_="TrEMBL:$_"; } 
    	elsif(/^ENS\w+\d\d+$/) { $_="ENSEMBL:$_"; }
      } ($refid,$repid);
    
    #old# $genedbxref{$id} .= "$refid," if($refid);
    #old# $genedbxref{$id} .= "$repid," if($repid and not $refid =~ /^CDD:/);
    
    $namedgenes{$name} .= "$id,"; #? if($pctalign >= $MIN_NAMEIDENT); # for uniq name retention
    
    ## FIXME: keep CDD names for .ann.txt, maybe .tbl submit as 2nd note
    ## UPD20mar: add BUSCO/ORTHODB/ORTHOMCL/ORTHOxxx as cddnames?
    my $iscdb= ($name =~ /(CDD|BUSCO|ORTHO\w*):/ or $refid =~ /^(CDD|BUSCO|ORTHO\w*):/)?1:0;
    $cddnames{$id}= $name if($iscdb and not $cddnames{$id});
    #o: $cddnames{$id}= $name if($name =~ /CDD:/ and not $cddnames{$id});

    if(UPD1912 and $haspct and $refid) {
      unless($aablast{$id} and $bscore{$id} >= $aln) {
        $aablast{$id}="$aln,$refid";  $bscore{$id}= $aln; }
      unless( $aablastref{$refid} and $bscore{$refid} > $aln ) {
        $aablastref{$refid}= "$aln,$id";  $bscore{$refid}= $aln;
        $aablastref{$id}= "$aln,$refid";   # revhash?
      }  
    }
            
    ## FIXME: 1st in dont replace.. not just for CDD: ?
    unless($genenames{$id}) { #was and $name =~ /CDD:/ 
      $pctalign ||= 0; $refid ||= 0;
      $genedbxref{$id} .= "$refid," if($refid);
      $genedbxref{$id} .= "$repid," if($repid and not $refid =~ /^CDD:/); #? CDD
      $genenames{$id}= $name;  $ngot++;
      $genenamepct{$id}= $pctalign;
       # repid Yes or No? this is by default RefID=UniRef50_xxxx and RepID=UniProt:xxxx_HUMAN
       
    } else { # UPD1912, append dbxrefs
      if($refid and not $genedbxref{$id} =~ m/$refid/){
        $genedbxref{$id} .= "$refid,";
      }
    }
  } close($inh);

  if(UPD1912) {
    my($nrefqual)= genename_refquals(\%aablast, \%aablastref);
  }
  
  $HAVE_gnames= $genenames if($ngot > 1 and $ngot >= 0.5 * $nin);  
  return($ngot,$nin);
}

sub genename_refquals {
  my( $aablast, $aablastref)= @_;
  my $MINBLASTSCORE= 60; # aablast only? bitscore always?
  my($nrefqual,$nid)=(0);
  
  for my $tid (sort keys %genedbxref) {
    $nid++;
    # my @rid= split ",", $refids;

    my $tbits= $$aablast{$tid} || "0,0";
    my($tbscore,$tbref)=split",",$tbits;

    my($tbrscore,$tbrefbest)= ($tbref) ? (split",", $$aablastref{$tbref}) : (0,0); 
    my $rbits2= $$aablastref{$tid} || "0,0";
    my($tbrscore2,$tbref2)= split",", $rbits2; # this is revhash tid => score,ref refbest

    unless( $tbscore == 0 or $tbits=~/^0,0/) {
      my $risbest= (($tbrefbest eq $tid) or ($tbscore >= $tbrscore))?1:0;
      my $risgood= ($tbrscore > $MINBLASTSCORE and $tbscore >= 0.90 * $tbrscore)?1:0;
      my $risok=   ($tbref2 and $tbrscore2 >= $MINBLASTSCORE)?1:0;
      if($risok and $risgood) { $risgood=0 if($tbrscore2 > $tbrscore); } 
      my $rfl= ($risbest)?"best" : ($risgood)?"good" : ($risok)?"ok" :"";
      if($rfl) {
        $tbits .= ",ref$rfl";
        my $refids= $genedbxref{$tid}||"";  
        $refids =~ s/$tbref\b/$tbref.ref$rfl/; #?? this is new place for qual?
        $genedbxref{$tid}= $refids;  $nrefqual++;
      }
    }
  }
  return($nrefqual,$nid);
}

sub parse_evgheader
{
  my($oid,$hdr,$trlen,$seqoid)= @_;
  $seqoid ||= $oid; # for split/dup gff gene IDs: Id_C1,2 Id_G2,3,.. names have seqoid
  $trlen  ||= 0;
    ## this becomes param or not?
    ## oid maybe pubid, check hdr for others
  #o# my $pubid= $pubids{$oid} || $oid; # is it ERR if no pubid{oid} ?
  my @xpubids=();  #UPD20mar for intermediate annot ids
  if(my $pd=$pubids{$oid}) { push @xpubids, $pd; }
  if(exists $puballoids{$oid}){ push @xpubids, sort keys %{$puballoids{$oid}}; }

  my $pubid= $xpubids[0] || $oid;
  my $protid= $pubid;
  $protid=~s/t(\d+)$/p$1/; # genbank requires diff protid from mrnaid
  $protid= $GDB_PREFIX.$protid if($GDB_PREFIX);
  #  protein_id      gnl|CacaoGD|Thecc1EG016762p1

#     if(my $lotagpre= $settings{'LOCUSTAG'}) {
#       my($pubidnum)= $pubid =~ m/$IDPREFIX(\d+)/;
#       $locustag= "$lotagpre$pubidnum" if($pubidnum);
#     }

  my %tblinfo= (pubid => $pubid, oid => $oid,  protid => $protid, locustag => 0,
      aaqual => "", trlen => $trlen, cdsoff => "", cdsor => 1, 
      name => "", namepct => 0, dbxref => '', cdd => ''); ## not yet "na" 
  # tblinfo.Specials:  Selcstop
  # tblinfo.addkeys:   locus/location/maploc, mapqual
  
  use constant CHECK_NAMEOIDS => 1; # yes keep
  my $nameoid= $seqoid;
  unless( $genenames{$nameoid} ) { #CHECK_NAMEOIDS and 
    my $nok=1; my @alloids=();
    if($hdr =~ m/\boid=([^\s;]+)/) {  @alloids=split",",$1; }
    for my $d (@alloids,@xpubids){ if( $genenames{$d} ) { $nameoid=$d; $nok=0; last; }  }
    # if($nok and $pubid ne $nameoid and $genenames{$pubid}) { $nameoid= $pubid; }
  }  

  if( $genenames{$nameoid} ) {
    $tblinfo{'name'}= $genenames{$nameoid};
    $tblinfo{'nameoid'}= $nameoid; 
    my $nap= $genenamepct{$nameoid} || 0;
    if($nap) {
      $tblinfo{'namepct'}= $nap;
      my ($nal)= $nap =~ m/,(\d+)/ ? $1 : $nap =~ m/^(\d+)/ ? $1 : 0; 
      $tblinfo{'namealn'}= $nal ; # UPD1912 want aa align, add from genenames??
      }
    $tblinfo{'nameref'}= $tblinfo{'dbxref'}=  $genedbxref{$nameoid}; # ||"na" # should this be 'nameref' instead?
    $tblinfo{'cdd'}=  $cddnames{$nameoid}; # ||"na"
  }
        
  my($cdsb,$cdse,$aafull)=(0,0,0);
  if($hdr =~ m/\b(?:offs|cdsoff)=([\d-]+)/) { my $cdsoff=$1; $tblinfo{'cdsoff'}= $cdsoff; 
    ($cdsb,$cdse)= split/[-]/,$cdsoff;  } # do in putseq
  if($hdr =~ m/\b(?:aaqual|aalen)=([^\s;]+)/) { my $aq=$1; $tblinfo{'aaqual'}= $aq; 
    ($aafull)= $aq =~ m/(complete|partial\w*)/; }
  if($hdr =~ m/\bclen=(\d+)/) { my $ln=$1; unless($trlen){ $trlen=$ln; $tblinfo{'trlen'}= $trlen; } } #  and not $trlen ; skip? 
  if($hdr =~ m/\bstrand=([^\s;]+)/) { $tblinfo{'cdsor'}= $1; } # expect all '+' from traa2mrna
  if($hdr =~ m/\bSelcstop=([^\s;]+)/) { $tblinfo{'Selcstop'}= $1; }  # Selcstop update 14.12.30
  if($hdr =~ m/\boid=([^\s;]+)/) { my $od=$1; 
    if($od eq $oid){$od="";} elsif($od=~/$oid/){} elsif($oid ne $pubid){$od="$oid,$od";} #?
    $tblinfo{'oid'}= $od if($od); } ## oid param maybe pubid, check others
  #^ add (?:Target|trg)= as oid alternate of gff
  if($hdr =~ m/\b(?:Target|trg)=([^\s;]+)/) { my $trg=$1; my $od=$tblinfo{'oid'}||"";
    unless($od=~/$trg/){ $od.="," if($od); $od.= $trg; $tblinfo{'oid'}= $od; } }
  if($hdr =~ m/\b(?:[Nn]ame|[Pp]roduct|[Pp]roduct[Nn]ame)=([^=\n;]+)/) { # Product_Name= ?
     my $na=$1; $tblinfo{'name'}= $na unless($tblinfo{'name'}); }
  if($hdr =~ m/\b[Nn]amepct=([^\s;]+)/) { my $nap=$1; 
    $tblinfo{'namepct'}= $nap unless($tblinfo{'namepct'}); }
  elsif($hdr =~ m/\b[Nn]amealn=([^\s;]+)/) { my $naa=$1;
    ## namepct =>  namealn=58p,17578/30278,17622; .. keep align stats?
    my($nap)= $naa=~m/^(\d+)/?$1:0;  my($nal)= $naa=~m/,(\d+)/?$1:0;
    $tblinfo{'namepct'}= $nap if($nap and not $tblinfo{'namepct'});
    $tblinfo{'namealn'}= $nal if($nal);
    }
  if($hdr =~ m/\bevgclass=([^\s;,]+)/) { $tblinfo{'evgclass'}= $1; } # upd1803
  if($hdr =~ m/\btype=(\w+)/) { $tblinfo{'seqtype'}= $1; } # upd1803; type| seqtype| moltype ?
    
    #?? merge hdr dbxref and genenames dbxref?  
    ## problems w/ new/old dbxref near same, drop old, treat genenames as accurate
    # .. check dup uniprot _species now..
  if($hdr =~ m/\b(?:[Dd]bxref|db_xref)=([^\s;]+)/) { my $dbx=$1; 
    if(my $tdx=$tblinfo{'dbxref'}) { 
      for my $dx (split",",$dbx) { 
        my($d,$x)=split":",$dx; 
        if($x and not $tdx=~/$x/) {
          my $ok=1;
          if($d=~/SwissProt|TrEMBL|UniProt/) { 
            my($sp)= $x=~m/(_\w+)$/; $ok=0 if($sp and $tdx =~ /$sp/);
          }
        $tdx.=",$dx" if($ok);  
        }
      }
      $tblinfo{'dbxref'}= $tdx;
    } else { $tblinfo{'dbxref'}= $dbx; }
  }

  map{ $tblinfo{$_} ||= "na" } qw(aaqual dbxref cdd);
  return \%tblinfo;
}


sub getOkFileset
{
  my($okpath,$SUFFIX,$okfiles,$trname)= @_;
  
  #FIXME.1712: added trname opt, may have many projects in okayset/ 
  # see also getmRNA() .. too complex to use here, elsewhere
  #    ($trf)= grep /$trname\.(mrna$|mrna.gz$)/, (@pubd, @okd); # drop okd here?
  #UPD1911:    $TRSUFFIX='tr|cdna|fasta|fna'; # bad** use instead NOT .mrna|cds|aa|log
  #UPD1911.was: $SUFFIX='tr|cdna|fasta|fna' unless($SUFFIX); # is this enough? NOT .mrna  
  # see also EVG_SEQSUFFIX
  
  my($oktr,$alttr)=("","");
  my @okd=(); 
  if($okfiles and ref($okfiles) and @$okfiles > 0) { @okd= @$okfiles; }
  elsif(-d $okpath) { opendir(D,$okpath); @okd= map{ chomp; "$okpath/$_" } readdir(D);  closedir(D); }
  if($trname) { @okd= grep /$trname/, @okd; } # fix for different subprojects in okayset/ ..
      
  if($SUFFIX) {
    ($oktr) = grep /.okay\.($SUFFIX)/, @okd;  
    ($alttr)= grep /.okalt\.($SUFFIX)/, @okd; 
  } 
  unless($SUFFIX and $oktr) { #UPD1911
    ($oktr) = grep{ not m/\.(aa|cds|mrna|log)/} grep /\.okay\./, @okd; 
    ($oktr) = grep{ not m/\.(aa|cds|mrna|log)/} grep /\.okalt\./, @okd; 
  }
  return($oktr,$alttr,\@okd); ## change to? (\@okd,$oktr,$alttr)  
}

sub getFileset
{
  my($okpath,$SUFFIX,$okfiles,$trname)= @_;
  #UPD1911: SUFFIX requirement is bug
  $SUFFIX='tr|cdna|fasta|fna' unless($SUFFIX); # is this enough? NOT .mrna  
  #FIXME.1712: added trname opt, may have many projects in okayset/ 
  my (@okd,@files); 
  if($okfiles and ref($okfiles) and @$okfiles > 0) { @okd= @$okfiles; }
  elsif(-d $okpath) { opendir(D,$okpath); @okd= map{ chomp; "$okpath/$_" } readdir(D);  closedir(D); }
  if($trname) { @okd= grep /$trname/, @okd; } # fix for different subprojects in okayset/ ..

  if($SUFFIX) {
    @files = grep /\.($SUFFIX)/, @okd;
  } else {
    @files= @okd; #UPD1911 ???
  }

  return(\@okd, @files); 
}

sub tidyupFileset { 
	my($tod,@td)= @_;  
	return unless($tod and @td);
	my @tdlist;
	mkdir($tod) unless(-d $tod);
	foreach my $fn (@td) { 
  	if($fn and not ($fn =~ m,$tod/,)) { 
    #   if(-f $fn ) {  # old only -f
    #     (my $tfn=$fn) =~ s,^[\w\.]+/,,;   ## ASSUMES tod is subdir in current path !!
    #     rename($fn,"$tod/$tfn");  push @tdlist, "$tod/$tfn";
    #   }
     
      #? (my $tfn=$fn) =~ s,^[\w\.]+/(.+),$1,;   ## ASSUMES tod is subdir in current path !!
      my $tfn= basename($fn); #^old way bad?
      
      my $todfn= "$tod/$tfn"; #? check exists? replace?
      if(-f $fn ) { 
        rename($fn, $todfn); push @tdlist, $todfn;
      } elsif( -d $fn ) { # move subdirs ok this way?
        rename($fn, $todfn); push @tdlist, $todfn;
      }
  	}
  } 
	my $n=@tdlist; my $n1= _min1($n-1,4); 
	loggit(0,"tidy: n=",$n, @tdlist[0..$n1]); 
}


=item AaQual evigene attribute
  
  AaQual is a transcript attribute extensively used by Evigene.
  Value is a tuple: "#aa-length,#coding-percent,Completeness"
  where aa-length is count of aa residues, including gaps (usually)
  coding-percent is %(CDS-length/mRNA-length)
  Completeness is controlled vocabulary: complete|partial3|partial5|partial 
    (partial=missing 5' and 3' ends, partial5=missing 5', ..)
    with other appended: -utrbad|-utrpoor|-gapbad|..
  It is calculated from proteins of mRNA transcripts following ORF translation.
  
  Evigene ORF sequences (.aa and .cds) and size table (.aa.qual) have this and
  other  ORF translation values, 
    offs=CDS-offset (b-e) in mRNA/cDNA, (e-b) for revcomp
    strand=+|- in cDNA/mRNA,
    clen=cDNA/mRNA length, 
    aalen=AaQual tuple or simple aa-length,

=item AaQual score
   
   This is integer value of Completeness vocabulary, with "complete" only as highest value.
   "partial", "utr" and "gap" attributes reduce score.
   Current range is -3..+3 (kAAQUAL_MIN..kAAQUAL_MAX) 
   
   A single numeric comparison of transcript Aa would include aa-size, coding% and completeness,
   for instance for selecting or sorting transcripts / proteins.
   Perhaps aascore = aa-length * codingpct/100 * (aaqual - kAAQUAL_MIN) / (kAAQUAL_MAX - kAAQUAL_MIN)
   
=cut

#above# use constant { kAAQUAL_MAX => 3, kAAQUAL_MIN => -3, kAAQUAL_NONE => 0, };
sub aaqualscore
{
  my($mqual)= @_;  $mqual ||="missing";
  my $mqv= kAAQUAL_NONE; 
  if($mqual =~ /complete/) { $mqv = kAAQUAL_MAX; } 
  elsif($mqual =~ /partial[35]/) { $mqv = kAAQUAL_MAX - 1; }
  if($mqual =~ /utrbad/) { $mqv = ($mqv >= kAAQUAL_MAX-1) ? kAAQUAL_NONE - 1 : $mqv - 2; } 
  elsif($mqual =~ /utrpoor/) { $mqv -= 1; }
  if($mqual =~ /gapbad/) { $mqv -= 1; } # or -2?
  return $mqv; # range is now -3..+3, was -2 .. +2
}


## replacing evigene/scripts/prot/aaqual.sh
sub makeAaQual {
  my($aaseq,$aaqualSetHash)= @_;
  my $doff=1; # $ENV{doff}; 
  my $doid=1; # no oids here? drop this? option?  $ENV{doid};
  our $makeAaQual_SETHASH= $aaqualSetHash||0;
  
  # minor update bug: this always makes new, should instead reuse old aasize..
  my $aasize= makename($aaseq,".aa.qual"); 
  if( -s $aasize ) {
    my $aqhash= getAaQual($aasize) if($makeAaQual_SETHASH); # yes or not?
    return($aasize);
  }
  
  my ($ok,$inh)= openRead($aaseq);  
  $ok= open(AAQ,'>',$aasize) if($ok);
  return unless($ok);
  
  # use constant  makeAaQual_SETHASH => 1;
  my($id,$aat,$aag,$al,$cl,$naa)= (0) x 9;
  my(@amo);
  
  sub puta { 
    my($id,$aat,$aag,$al,$cl,@amore)= @_; our $makeAaQual_SETHASH;
    $al=($aat+$aag).",na" if($al eq "na"); 
    print AAQ join("\t",$id,$aat,$aag,$al,$cl,@amore)."\n"; 
    ## option here: set global hashes
    if($makeAaQual_SETHASH) {
       my($aww,$pctcds,$aqual1)=split",",$al; 
       $pctcds =~ s/\%//;      
       $aqual1 .= "-gapbad" if($aag>0 and (100*$aag/($aat+$aag) > $BAD_GAPS)); # add qual flag for gap/(alen+gap) > MAXGAP:  
       my $acv= aaqualscore($aqual1);
       #NOT YET# $AASIZEH->{$id}=$aat; 
       #o $AAQUALH->{$id}="$aat,$pctcds,$acv,$aqual1"; 
       $AAQUALH->{$id}= join",", $aat,$pctcds,$acv,$aqual1,@amore; #UPD22a fix need @amore
     }
    return 1;
  }

  if($makeAaQual_SETHASH) {  $AAQUALF=$aasize; $AAQUALH={}; } # $AASIZEH={}; 
  while(<$inh>) {
    if(/^>(\S+)/) { my $td=$1; 
      $naa+= puta($id,$aat,$aag,$al,$cl,@amo) if($id);  
      $id=$td; $aat=$aag=0;  @amo=();
      ($al)=m/aalen=([^;\s]+)/; $al||="na";
      ($cl)=m/clen=(\d+)/; $cl||=0;
      if($doff){ 
        my($or) = (m/strand=(.)/)?$1:"."; 
        my($ofs)= (m/offs=([\d-]+)/)? "$1:$or" : 0; 
        push @amo, $ofs;  #o $cl.= ($ofs)?"\t$ofs:$or":"\t0";  
        } 
      if($doid){ my $oid;
        unless(($oid)=m/oid=([^;\s]+)/) { ($oid)=m/gene[=:]([^;\s]+)/; } 
        $oid||="noid"; push @amo,$oid; #o $cl.="\t$oid"; 
      } 
      if(UPD1908) { # read new attr: codepot=Noncode/Code/.. 
        if( my($cp)= m/codepot=([^;\s]+)/ ) { unshift(@amo,$cp); } # cds.qual column after clen, before ofs,
        #O if( my($cp)= m/codepot=(\w+)/ ) { unshift(@amo,$cp); } # cds.qual column after clen, before ofs,
      }
    } else { s/\*$//; $aat += tr/A-WYZa-wyz/A-WYZa-wyz/; $aag += tr/Xx\*/Xx\*/; }
  } 
  $naa+= puta($id,$aat,$aag,$al,$cl,@amo);  
  close(AAQ); close($inh);
  loggit(0, "makeAaQual: naa=$naa IN $aasize\n") if($DEBUG);
  ## optionally set global hashes: $AAQUALF,$AAQUALH,$AASIZEH
  return($aasize);
}

=item isStrandedRNA
  
  cds-orient of longest 1000 aa tells if transcripts are from stranded RNA.
  $ismrna= isStrandedRNA($aaqualHash [,$aaqualFile]);
  ($ismrna,$nfwd,$nrev,$nt)= isStrandedRNA($aaqualHash);
  
=cut
 
sub isStrandedRNA {  #UPD20apr
  my($inaaqualh,$aaqualf,)= @_; # one param as hashref or filename?
  my($aaqualh);
  if($inaaqualh and ref($inaaqualh) =~ /HASH/) {
    $aaqualh= $inaaqualh;  
  } elsif($AAQUALH and ref($AAQUALH) =~ /HASH/) {
    $aaqualh= $AAQUALH;  
  } elsif($aaqualf and -f $aaqualf) {
    $aaqualh= getAaQual($aaqualf);
  }
  return 0 unless($aaqualh and ref($aaqualh) =~ /HASH/); # -1 for err?
  
  # cds-orient of longest 1000 aa tells if transcripts are from stranded RNA.
  use constant { ksMinTest => 100, psMinFwd => 0.90 };
  my %nor;  
  my @aaord= sort{ my($aw)= split",",$aaqualh->{$a}; my($bw)= split",",$aaqualh->{$b};
                  $bw <=> $aw or $a cmp $b } keys %$aaqualh;
  for my $id (splice(@aaord,0,1000)) {
    my @v= split",",$aaqualh->{$id}; next unless(@v > 4);
    my($alen,$pctcds,$acv,$aqual1,$cpot,$ofs)=@v; 
    # my($alen,$pctcds,$acv,$aqual1,$cpot,$trw,$ofs)=@v; #UPD22a is trw==cl there??? need fuckin offset 
    #? not here, but want other aa1k stats? ave(alen,pctcds,complete,acv,cpot)
    $ofs=$cpot if(@v==5 or $cpot=~m/\d\-\d/);
    unless($ofs =~ m/(\d+)\-(\d+)/) { my($offf)= grep{ m/(\d+)\-(\d+)/ } @v; $ofs= $offf if($offf); }#UPD22a
    if($ofs =~ m/(\d+)\-(\d+)/) { 
      my($cb,$ce)=($1,$2); 
      my($or)= ($ofs=~m/$ce:(.)/)?$1:0;
      if($cb > $ce or $or eq 'r'){ $or='-'; } 
      elsif($or eq '0' or $or eq 'f'){ $or='+'; }
      $nor{$or}++; 
      }
  }   
  my $nfwd= $nor{'+'} || 0; my $nrev= $nor{'-'} || 0; my $nt=$nfwd+$nrev;
  my $ismrna= ($nt >= ksMinTest and $nfwd > $nrev and $nfwd >= psMinFwd * $nt)?1:0;
  return (wantarray) ? ($ismrna,$nfwd,$nrev,$nt) : $ismrna;
}


# getAaQual returns $aaqual->{$id}="$alen,$pctcds,$acv,$aqual1"; see also asmdupfilter:readSizes()
sub getAaQual {
  my($aaqualf,$inaaqualh)= @_;
  my($naa,$ntr,$nerr,$ok,$inh)=(0) x 10;
  
  # our ($AAQUALF,$AAQUALH,$AASIZEH) = ("",undef,undef);
  if($AAQUALF and $AAQUALF eq $aaqualf and ref($AAQUALH)) { 
    $naa= scalar(keys %$AAQUALH);
    if($DEBUG) { my($id1)= (keys %$AAQUALH)[0];
      loggit(0, "getAaQual: naa=$naa in $AAQUALF, val1 $id1=",$AAQUALH->{$id1});
      }
    return (wantarray) ? ($AAQUALH,$naa) : $AAQUALH; # %aasize also? NOT YET ,$AASIZEH
  }
  
  my %aaqual = (); # use global hash, keep in sync w/ file name
  if($inaaqualh and ref($inaaqualh)) {
    $AAQUALH= $inaaqualh; # incase read aaqualf fails
    %aaqual= %$inaaqualh; # copy contents ? ids not in file?
  }
  
  #NOT YET# my %aasize =(); # use global hash, keep in sync w/ file name
  if($aaqualf) {  
    ## fix for aacount gaps: id,size,gaps : NOT NOW, aa.qual: id,size-gaps,gaps,..
    ## drop faCount? require aa.qual here?

    ($ok,$inh)= openRead($aaqualf); # $ok= open($inh,$aaqualf);
    unless($ok) { loggit(1,"ERR: getAaQual $aaqualf"); return 0; } # die "FAIL: read $aaqualf ..."; 

    while(<$inh>) { 
      next if(/^\W/ or /^total/); 
      my($id,$alen,$gap,$aqual,$tlen,$codepot,$ofs,$oid)=(0) x 9;
      my @v=split; # aa.qual cols; gap is removed from alen
      ($id,$alen,$gap,$aqual,$tlen)=@v;
      if(UPD1908 and @v > 5) {
        # read new attr: codepot=Noncode/Code/..   cds.qual column after clen, before ofs,
        my $nv= @v; my $o=5; 
        if($o<$nv and $v[$o] =~ /^(No|Co|Un)/) { $codepot=$v[$o]; $o++;  $codepot=~s/,.*//; } # codepot likely not there
        if($o<$nv and $v[$o] =~ /^\d/){ $ofs=$v[$o]; $o++; unless($ofs =~ m/:/) { my($b,$e)=split/\-/,$ofs; $ofs.=":" . (($b>$e)?'-':'+'); } } # likely
        if($o<$nv and $v[$o] =~ /^\w/){ $oid=$v[$o]; $o++; } # maybe
      }
      unless($alen =~ /^\d/) { $nerr++; next; }
      
      if($aqual) { 
        $aqual .= "-gapbad" if($gap>0 and (100*$gap/($alen+$gap) > $BAD_GAPS)); # add qual flag for gap/(alen+gap) > MAXGAP:  
        my $acv= aaqualscore($aqual);
        my($aww,$pctcds,$aqual1)=split",",$aqual;  
        $pctcds =~ s/\%//;  
        my $qval="$alen,$pctcds,$acv,$aqual1";
        $qval.=",$codepot" if($codepot); # should be Coding/Noncode/Unknown
        $qval.=",$ofs" if($ofs);  
        $aaqual{$id}=$qval; 
        #?? $aqv{$id}= $alen * $pctcds * $acv; # want something like this, better, for 1 num score
        #? want# if($tlen =~ /^\d/) { $trsize{$id}= $tlen; $ntr++; }
      } else {
        my($pctcds,$aqual1)=(99,"na"); my $acv= aaqualscore($aqual1);
        $aaqual{$id}="$alen,$pctcds,$acv,$aqual1"; 
      }
      #NOT YET# $aasize{$id}=$alen; #? need both aaqual and aasize hashes?
      $naa++; 
    } close($inh); 
    
    if($naa) { $AAQUALF=$aaqualf; $AAQUALH= \%aaqual; } ## dont need:  $AASIZEH= \%aasize;
  }
  
  if($DEBUG) { my($id1)= (keys %$AAQUALH)[0];
    loggit(0, "getAaQual: naa=$naa in $AAQUALF, val1 $id1=",$AAQUALH->{$id1});
    }
  ## AAQUALH isnt reset unless $naa..
  return (wantarray) ? ($AAQUALH,$naa) : $AAQUALH; # which ?
}



sub getmRNA   ## move to cdnasubs.pm ?
{
  my($okpath,$trname,$pubdir,$dieOrWarn,$options)= @_; # opts was ,$ADDutrorf
  my($cdnaseq)=(""); # == mrna (oriented), not cdna.parts.tr
  my @ret_mrna_aa_cds=(); # for return (wantarray) ? (mrna,aa,cds) : mrna
  $dieOrWarn ||= LOG_DIE; # FIXME: dont die here, for some uses at least..  
  $options   ||= "";
  
	use constant ALSOMAKE_AACDS => 1;
  my $ADDutrorf= ($options =~ /noutr/)?0:1; # what? always check for okayset/*.utrorf.mrna ?
  # UPD20j: $ADDreor = 1; # for okayset/name.okreor.{aa,cds,tr} restranded set of trclass_resolve_strandmix.pl
  my $ADDreor = ($options =~ /reor|restrand/)?1:0;
  my $ADDcull = ($options =~ /addcull|okcull/)?1:0; # UPD20mar
    
  # UPD1911: TRSUFFIX requirement is BUG .. do away w/ that, or change evg pipes to create files w/ that
  # maybe change TRSUFFIX to not(aa|cds|mrna|log|xxx) ? see %EVG_SEQSUFFIX
  my $TRSUFFIX='tr|cdna|fasta|fna'; # is this enough? NOT .mrna|cds|aa,  
  ## ?? add .fa ?? BUT conflict in now publicset/*.{mrna,cds,aa}_pub.fa **

  my($oktr,$alttr,$okd)= getOkFileset($okpath,$TRSUFFIX,undef,$trname);
  my @okd= @$okd;
  my @pubd=();
  if($pubdir and -d $pubdir) { 
    # UPD1911: TRSUFFIX needed for pubdir, has other data
    my($pubd)= getFileset($pubdir, $TRSUFFIX, undef, $trname); 
    @pubd= @$pubd; 
    }
 
  ## another bad call: #egr: get_evgtrset= publicset/locust1all5asm.p4.mrna . locust1all5asm
  ##   instead of publicset/locust1all5asm.mrna .. my test version mistake..
  ## Ugh! bad call: #er2g: get_evgtrset= ./okayset/pogonus1all3.mrna0.ann.txt . pogonus1all3
  ## need grep /\.mrna$|\.mrna.gz$|.mrna.fa$/ ??? 
  
  # UPD1807: allow other okayset/okXXX parts? or add to okay/okalt files?
  
  # UPD20ja .. ugh if 2+ okayset/.mrna, then which? dont want .cull
  # UPD20may .. ugh2, got okreor.mrna ? instead of okay.mrna + okalt.mrna  ADDreor bug
  # UPD20may: bug got okreor, need okay.mrna -- fix in getmRNA()
  #egr: EvidentialGene trclass2pubset VERSION 2020.03.15
  #egr: get_evgtrset= ./okayset/plDZQM.okreor.mrna . plDZQM <<< == cdnaseq from splice(@ret_mrna_aa_cds,0,3); 
  #.. this bug is in tr2aacds4 stage2 calls, after reorient (which places .okreor.mrna,aa into okayset/ )
  #..  at this stage1-2 transit, okayset/ contains name.okay,okalt.tr and name.okreor.mrna 
  #..  publicset/ exists in transit, renamed publicset_orold before 2nd trclass2pubset makes new publicset > okayset2nd
  #..  prior preference here for .mrna suf picks only okayset/okreor.mrna
  #..  but above getOkFileset(okayset,suf='tr|cdna|fasta|fna') gets the orig okay,okalt.tr set
  
  my $trf=""; my @trf=();
  if(@pubd) {
    # *should be* trname.mrna_pub.fa(.gz), avoiding cull/reor/ncrna/others
    @trf= grep { m/$trname\.mrna_pub\.fa/} (@pubd);   # publicset/ final .mrna_pub
    if(@trf==1) { $trf=$trf[0]; }
    else {
      @trf=  grep { m/$trname.*\.(mrna|mrna.gz)$/ and not m/cull/ } (@pubd); # ?? is this useful? below is okayset stg2 mrna
      if(@trf==1) { $trf=$trf[0]; }
    }  
  }
  if(not $trf and @okd) { # tr2aacds4 stage2 okayset/trname.okay.mrna
    @trf= grep { m/$trname.okay\.(mrna$|mrna.gz$)/} (@okd);  # this should be okayset stage2 main mrna seq
    if(@trf==1) { $trf=$trf[0]; }
  }
  
  # ?? this way safely picks okayset/name.okay.mrna of stg2 ?? 
  # .. maybe not, $TRSUFFIX='tr|cdna|fasta|fna' excludes .mrna
  # .. safer to skip this choice, fall into 'Make mRNA $cdnatmp from okayset transcripts'
  if(not $trf and $oktr and not $alttr) { 
    #? $trf=$oktr; #? is this okayset1st/trname.okay.tr ? is that wanted here?
  }
  
  ## not this way, mistakes of precedence oksubset.mrna > okay,okalt.tr
  # if( 0 and not $trf) { 
  #   @trf=  grep {  m/$trname.*\.(mrna$|mrna.gz$)/ and not m/cull/ } (@okd); # drop okd here?  UPD1912 tr2aacds4 okayset/*.okay.mrna is desired
  #   if(@trf==1 and $trf[0] =~ m/\.okay/) { $trf=$trf[0]; } 
  #   # elsif(grep/okalt/,@trf) { } # SKIP this trf, need to (re)make pubset / okaystage2
  #   # else also remake? other name tags ?
  #   #??no unless($trf) { ($trf)= grep /\.mrna_pub\.fa|\.mrna$|\.mrna.gz$/, (@pubd); } # , @okd want this or not? YES UPD1912
  # }

  #o ($trf)= grep /$trname.*\.(mrna$|mrna.gz$)/, (@pubd, @okd); # drop okd here?  UPD1912 tr2aacds4 okayset/*.okay.mrna is desired
  #o20ja ($trf)=  grep {  /$trname.*\.(mrna$|mrna.gz$)/ and not m/cull/ } (@pubd, @okd); # drop okd here?  UPD1912 tr2aacds4 okayset/*.okay.mrna is desired
  #o20ja unless($trf) { ($trf)= grep /\.mrna_pub\.fa|\.mrna$|\.mrna.gz$/, (@pubd, @okd); } # , @okd want this or not? YES UPD1912

  if($trf and -s $trf) { 
    $cdnaseq= $trf; 
    @ret_mrna_aa_cds= ($cdnaseq); 
    foreach my $suf (qw(aa cds)) { (my $okf= $trf) =~ s/\.mrna/.$suf/;  push @ret_mrna_aa_cds, $okf; } #? if( -f $okf);
    
    ## need okay.mrna,aa,cds first in @ret
    if($ADDcull) {  # check here for okayset2nd/trname.cull.seq ? add to @ret_mrna_aa_cds
      my($trc)= grep { m/$trname.*\.(mrna$|mrna.gz$)/ and  m/cull/ } (@pubd, @okd);
      if($trc) { push @ret_mrna_aa_cds, $trc;
        foreach my $suf (qw(aa cds)) { (my $okf= $trc) =~ s/\.mrna/.$suf/;  push @ret_mrna_aa_cds, $okf; } #? if( -f $okf);
      }
    }
     
     #UPD20ja: missing .aa,.cds partial .mrna in pubdir from new stg2 okayset/name.okay.{mrna,aa,cds}
     #  if(ALSOMAKE_AACDS and $pubdir) .. need to make pubdir/ copy over okayset/name.okay.seq to pubdir/name.seq ??
     #  mimic publicset/name.{mrna,aa,cds} of  below ALSOMAKE_AACDS  
     # BUT shouldnt need this .. caller does same from trf.mrna > .aa,.cds .. bug is in trclass2pubset.pl makePublicset
    if(0 and ALSOMAKE_AACDS and $pubdir and not ($trf =~ m,$pubdir,) ) {
      mkdir($pubdir) unless(-d $pubdir);      
  		foreach my $suf (qw(mrna aa cds)) {
  	    (my $okf= $trf) =~ s/\.mrna/.$suf/; #"$okpath/$trname.$pt.$sf";
        my $pubf= "$pubdir/$trname.$suf"; 
        #?? $cdnaseq=$pubf if($suf eq 'mrna'); # switch dirs for links ??
        symlink("../$okf", $pubf) if( -f $okf);
  		}
    }
     
  } else { 
    my $okall=0;
    my $cdnatmp= ($pubdir) ? $pubdir : $okpath;
    $cdnatmp .= "/$trname.mrna";  
    mkdir($pubdir) if($pubdir and not -d $pubdir);

    ## FIXME: utrorf : made okayset/*.utrorf.{mrna,aa,cds} ; merge into update_mrna_fileset() or getmRNA/okayset ??

    loggit(0,"Make mRNA $cdnatmp from okayset transcripts");
    my($okaa) = grep /.okay\.aa$|.okay\.aa.gz$/, @okd;  
    my($altaa)= grep /.okalt\.aa$|.okalt\.aa.gz$/, @okd; 
    # #.. problem here .aa.gz not .aa; should makename() do optionally?
    my @okreor= grep /.okreor\.(aa|cds|mrna|tr)/, @okd; #UPD20may: should now be okreor.{aa,cds,mrna}[.gz]
    
    ## FIXME: fail here if missing oktr ?
		$APPtraa2cds= findevigeneapp("prot/traa2cds.pl") unless($APPtraa2cds); #  
		my $tropts="-trout "; $tropts .= "-nomiss " if($ADDutrorf);
		
=item traa2cds new usage

  .. detect version of traa2cds? traa2cds.pl -help : VERSION >= 2018
  .. need all tr inputs to resolve utrorfs:  
  prot/traa2cds.pl -mrnaout -aa $pt.okay.aa -aa $pt.okalt.aa -cdna $pt.okay.cdna -cdna $pt.okalt.cdna -out $pt.okall.mrna
  .. OR can make mrna of okay,okalt separately using one -aa, but both -cdna (only aa-valid IDs written)

=cut
 		
 		our %headerupd= (); # for UPD20UORF
    sub headup { 
      my($id,$h)=@_; our(%headerupd);
      if(exists $headerupd{$id}) { 
        for my $k (sort keys %{$headerupd{$id}}) {
          my $v= $headerupd{$id}{$k}||0;
          unless($h =~ s/\b$k=[^;\s]+/$k=$v/){ $h =~ s/$/ $k=$v;/; }
        }
      } 
      return $h; 
      } 
        

if(TRAA2CDS_2018) {
    $ADDutrorf=0;
    my $cmd= "$APPtraa2cds -mrnaout -out $cdnatmp";
    my $okm=($oktr and -f $oktr and -f $okaa);
    my $oka=($alttr and -f $alttr and -f $altaa);
    
    if($ADDreor and @okreor) {
      ## UPD20may: shouldnt need to recalc okreor.mrna,  BUT need okreor.mrna at top of output cdnatmp
      ## .. now only way w/o rewrite traa2cds is add okreor.mrna as if .tr/.cdna
      my($reaa)= grep /\.(aa|aa.gz)$/, @okreor;
      my($retr)= grep /\.(mrna|mrna.gz)$/, @okreor;
      unless($retr) { ($retr)= grep /\.(tr|tr.gz)$/, @okreor; }
      $cmd.=" -ignoredupids -aa $reaa -cdna $retr" if($reaa and  $retr);
    ##... UPD20mayreor mrna <> cds bug due to old okreor.tr lookup, now okreor.mrna
    #  # FIXME these ids take precedence over ok/alt set in merge to mrna set
    #  # need APPtraa2cds to handle dup id seqs properly UPD20jan
    #  my($reaa)= grep /.okreor\.aa$|.okreor\.aa.gz$/, @okd; 
    #  my($retr)= grep /.okreor\.tr$|.okreor\.tr.gz$/, @okd; # UGH, this is now .mrna not .tr !!!
    #  $cmd.=" -ignoredupids -aa $reaa -cdna $retr" if($reaa and  $retr);
    }
    
    $cmd.=" -aa $okaa -cdna $oktr" if($okm);
    $cmd.=" -aa $altaa -cdna $alttr" if($oka);
    my $err= ($okm or $oka)? runcmd($cmd) : -2;
    $okall +=2  unless($err);
    
    if(UPD20UORF) {
      # look in mrna hdr for uorfspan=2680-3115/3115, update aa,cds hdr w/ new offs,clen
      # >idutrorf type=mRNA; aalen=104,72%,complete; clen=436;  strand=+; offs=10-324;uorfspan=2680-3115/3115;
      my($hmrna,$nmrna)= faheaders($cdnatmp);
      for my $id (keys %$hmrna) {
        my $hd=$hmrna->{$id} or next;
        if($hd =~ m/uorfspan=([^;\s]+)/) { my $uosp=$1;
          $headerupd{$id}{uorfspan}= $uosp;
          for my $k (qw(clen offs strand)) {
            my($upv)= $hd =~ m/$k=([^;\s]+)/?$1:0;
            $headerupd{$id}{$k}=$upv if($upv); 
          }  
        }
      }
    }  
    
} else {
    if($oktr and -f $oktr and -f $okaa) { 
      my $err= runcmd("$APPtraa2cds $tropts -cdna $oktr -aa $okaa -out stdout >> $cdnatmp");
      $okall++ unless($err);
    }
    if($alttr and -f $alttr and -f $altaa) { 
      my $err= runcmd("$APPtraa2cds $tropts -cdna $alttr -aa $altaa  -out stdout >> $cdnatmp");
      $okall++ unless($err);
    }
    
    if($ADDutrorf and $okall > 0) {
 			my($okin) = grep /.utrorf.mrna$|.utrorf.mrna.gz$/, @okd;  
   		if($okin) { runcmd("cat $okin >> $cdnatmp"); loggit(0,"add $okin to $cdnatmp"); } # err check? loggit?
   		else { $ADDutrorf=0; } # dont do .aa,cds
    }
}
    
    $cdnaseq= $cdnatmp if(-s $cdnatmp);
    @ret_mrna_aa_cds= ($cdnaseq); 
    my $seqok= ($okall>1 and -s $cdnaseq)?1:0;
    unless($seqok) {
      map{ $_ ||="missing"; } ($oktr,$okaa,$alttr,$altaa);
      loggit($dieOrWarn,"make mRNA $cdnatmp: missing files oktr=$oktr,okaa=$okaa,alttr=$alttr,altaa=$altaa,
    	from (oktr,alttr)=getOkFileset(path=$okpath,suffix=$TRSUFFIX)") ; 
    	return ($cdnaseq); # may exist .. what then?
    }
    
    # FIXmaybe: also make pubdir/.aa,.cds along with .mrna ? see hassle in update_mrna_fileset 
    # FIXME2: ok*.aa and utrorf.aa have dup entries, use utrorf if there.
    if(ALSOMAKE_AACDS and $pubdir) {
  		my($ok,$hin,$hout,$fout,$okin,$altin,$utrin,$reorin); 
  		my %ids=();
      # foreach my $suf (qw(aa cds))  
  		foreach my $suf (".aa",".cds") {
				$fout= makename($cdnatmp,$suf);  
				push @ret_mrna_aa_cds, $fout;
				($okin) = grep /.okay$suf$|.okay$suf.gz$/, @okd;  
				($altin)= grep /.okalt$suf$|.okalt$suf.gz$/, @okd; 
				($utrin)= ($ADDutrorf) ? grep(/.utrorf$suf$|.utrorf$suf.gz$/, @okd) : (); 
         # ADDreor ids take precedence over ok/alt set in merge to mrna set
				($reorin)= ($ADDreor) ? grep /.okreor$suf$|.okreor$suf.gz$/, @okd : (); 
				
        #UPD20mar: opt to add cull.seq : replacing okalt?
        # to handle this okayset4v/name.{okay,cull}.seq > prepubset/name.seq > refilter_newevd > publicset/name.seq
        # NOT here, this is for okayset1st/ with okay,okalt seqs
        #if($ADDcull and not $altin) { ($altin)= grep /.cull$suf$|.cull$suf.gz$/, @okd;  }
        
        # UPD20UORF : need to alter utrorf header info now for aa,cds to match utr.mRNA: new cdsoff=, mrna clen=
        # .. traa2cds -mrnaout updates mrna header info
        # orfbug_Auschi.aa,cds
        # >AuschiNNNutrorf type=protein; aalen=104,72%,complete; uorfcut=2680-3115;uorfoffs=10-324:+; clen=3115; strand=+; offs=2689-3003;
        # orfbug_Auschi.traa2mrna
        # >AuschiNNNutrorf type=mRNA; aalen=104,72%,complete; clen=436;  strand=+; offs=10-324; AuschiEVm020094t1.utrorfbug.true 

				if(($okin or $altin) and not -f $fout) { #UPD20mar: was (okin and altin)
					$ok= open($hout,'>',$fout);	my(%did,$did);
					if($reorin) { ($ok,$hin)= openRead($reorin); 
					  # FIXME?? Check: did{IDrevorf} needs also did{ID} ?  dont want both aa forms
						while(<$hin>){ if(/^>(\S+)/) { $did=($did{$1}++)?1:0; $_= headup($1,$_); } print $hout $_ unless($did); } close($hin); }
					if($utrin) { ($ok,$hin)= openRead($utrin); 
						while(<$hin>){ if(/^>(\S+)/) { $did=($did{$1}++)?1:0; $_= headup($1,$_); } print $hout $_ unless($did); } close($hin); }
					if($okin) { ($ok,$hin)= openRead($okin);  
						while(<$hin>){ if(/^>(\S+)/) { $did=($did{$1}++)?1:0; $_= headup($1,$_); } print $hout $_ unless($did); } close($hin); }
					if($altin) { ($ok,$hin)= openRead($altin); 
						while(<$hin>){ if(/^>(\S+)/) { $did=($did{$1}++)?1:0; $_= headup($1,$_); } print $hout $_ unless($did); } close($hin); }
					close($hout);
					map{ $ids{$_}{$suf}=$did{$_} } keys %did;
					}
				}
				
  	  ## should check for ID agreement among output .mrna, .aa, .cds  
			my $nidok=0;
			($ok,$hin)= openRead($cdnaseq); while(<$hin>){ if(/^>(\S+)/) { $ids{$1}{'.mrna'}++; } } close($hin);
			foreach my $id (sort keys %ids) {
				my @suf= sort keys %{$ids{$id}}; 
				if(@suf==3) { $nidok++; } 
				else { loggit(1,"ERR:getmRNA-misspart $id:",@suf); } # FIXME: okreor/revorf seen here
			}	
			loggit(0,"getmRNA $pubdir/mrna,aa,cds nid=",$nidok);
    }
  }
  
  return (wantarray) ? @ret_mrna_aa_cds : ($cdnaseq);  # , $aaseq, $cdsseq    
}


## from asmrna_trimvec : generalize mrna_update_fileset ..
## see tr2aacds.pl:asmdupclass_fileset
sub update_mrna_fileset
{
  my($trpath, $inmrna, $updatename, $trimids, @trimfiles)= @_; 
  my $upstatus=0;
  # outputs now should go to pubdir, but use inmrna path.
  my($trimmrna, $trimaa, $trimcds, $trimidfile)= @trimfiles; # hash instead? for below %fset
	map{ $_||="" } ($trimmrna, $trimaa, $trimcds, $trimidfile);
	my %okids=(); ## return hash of valid ids * AND add oids if found, for cross checks
	
	my $ntrim=0;
	if($trimids and ref($trimids)) {
		$ntrim= scalar(keys %$trimids);

	} elsif($trimidfile and -f $trimidfile) { # hash of ids in trimfiles, regenerate from trimidfile
   	$trimids={};  
   	my($ok,$hin)= openRead($trimidfile); 
   	while(<$hin>) { if(/^\w+/) { 
   		my($id,$action,@info)=split; # parse more of trimidfile? expect: id, action=DROP/OKCUT/PROBLEM, trimnotes, aanewhdr 
   		$trimids->{$id}=$action; $ntrim++; } 
   	} close($hin);

	} else { # ERROR unless -f flagtrimvec
 		#below# loggit(1, "ERR update_fileset  empty trim ids= $trimidfile"); 
	}
	
## ........ **#@&@!*% this file name wrangling is a big waste of time ..........
## ........ use fewer naming choices ???  no okdir for pubdir set
  ##upd1807: work on mrna_pub.fa => aa_pub.fa, cds_pub.fa names
  
	my $flagtrimvec= makename($inmrna,'.'.$updatename); 
	my($outaa,$outcds);
	if($inmrna =~ m/mrna_\w+\./) { # publicset/mrna_(pub|cull|xxx).fa
	  ($outaa = $inmrna) =~ s/mrna_/aa_/;
	  ($outcds = $inmrna) =~ s/mrna_/cds_/;
	} else {
    $outaa = makename($inmrna,".aa");  
    $outcds= makename($inmrna,".cds"); 
  }
  
  my($pubdir);
  #badlocal#(my $mrnapath= $inmrna) =~ s,/[^/]+$,,; ## THIS MAY BE WRONG NOW .. publicset/ vs okayset/
  ## ^^ FIXME localdir ./name.mrna == name.mrna 
  (my $mrnapath= $inmrna) =~ s,[^/]+$,,; 
  unless($mrnapath) { $mrnapath="./"; } else { $mrnapath =~ s,/$,,; }
  ($pubdir)= getFileset($mrnapath);
  
  ## drop okalt here? see above getmRNA/ALSOMAKE_AACDS, expect/require pubdir/aa,cds?
  my $aapatt= basename($outaa);
  my($inaaseq) = grep /$aapatt$|$aapatt.gz$/, (@$pubdir);  
  if($inaaseq) { $outaa=$inaaseq;  } else { $inaaseq=""; } # fail? ignore miss
  
  my $cdspatt=basename($outcds);
  my($incdsseq) = grep /$cdspatt$|$cdspatt.gz$/, (@$pubdir);  
  if($incdsseq) { $outcds=$incdsseq; } else { $incdsseq=""; } # fail?
  
  ## FIXME:for -novectrim set updatename == 
  if( $updatename =~ /SKIP/ or -f $flagtrimvec) {
  	my($ok,$hin)= openRead($inmrna); 
  	while(<$hin>) { if(/^>(\S+)/) { my $oid= (m/\boid=([^,;\s]+)/)?$1:1;  $okids{$1}=$oid; } } 
  	close($hin);
  	return ( 0, $inmrna, $outaa, $outcds, \%okids);   	# what if outaa,outcds missing?
  }
  if($ntrim < 1 or ! -s $trimmrna) { # fixme for -novectrim
		loggit(1, "ERR update_fileset  empty trim files= $trimmrna, $trimidfile"); 
		return (-1, $inmrna, $outaa, $outcds, \%okids);
  	}

  my $upmrna  = makename($inmrna,".mrna_upd"); 
  my $upaaseq = makename($inmrna,".aa_upd"); # failed empty outaa from above
  my $upcdsseq= makename($inmrna,".cds_upd");  # failed empty outcds from above
  
  $outaa =~ s/.gz$//; $outcds =~ s/.gz$//; # output fupname
  (my $outmrna= $inmrna)=~ s/.gz$//;
  my %fset= (
            # $nup,$nsame,$fin,$ftrim,$fup,$fupname
    mrna => [ 0, 0, $inmrna, $trimmrna, $upmrna, $outmrna ],
    aa   => [ 0, 0, $inaaseq, $trimaa, $upaaseq, $outaa],  #?? fix here for .okay.aa + .okalt.aa ?
    cds  => [ 0, 0, $incdsseq, $trimcds, $upcdsseq, $outcds],
  );
  
  foreach my $suf (sort keys %fset) { # is inmrna always .mrna ?
    my($ok,$hin,$hup,%keptids);
    my($nup,$nsame,$fin,$ftrim,$fup,$fupname)= @{$fset{$suf}};

      # is one only allowed here? -s ftrim but no fin?
    if(-s $fin and -s $ftrim) {
      %keptids=();   
      $ok= open($hup,'>',$fup); # unless($ok)...
      ($ok,$hin)= openRead($fin);  
      $ok=0; while(<$hin>) { 
      	if(/^>(\S+)/) { my $d=$1; my $oid= (m/\boid=([^,;\s]+)/)?$1:1; 
      	$ok=($trimids->{$d})?0:1; if($ok){ $keptids{$d}=$oid; $nsame++; } } 
      	print $hup $_ if($ok); }
      close($hin);
      
    	## pull trimset/uvcut.$suf, check? for trimids{id} and/or collect above kept ids
      ($ok,$hin)= openRead($ftrim); 
      $ok=0; while(<$hin>) { 
      	if(/^>(\S+)/) { 
      	  my $d=$1; my $oid= (m/\boid=([^,;\s]+)/)?$1:1; 
      	  $ok=($trimids->{$d} and not $keptids{$d})?1:0; 
      		if($ok){ $keptids{$d}=$oid; $nup++; } }
      	print $hup $_ if($ok); } 
      close($hin);
      close($hup);
      $fset{$suf}->[0]= $nup; $fset{$suf}->[1]= $nsame;  # $fset{$suf}= [$nup,$nsame,$fin,$fin2,$ftrim,$fup,$fupname];
      $upstatus++ if($nup>0);
    	%okids= %keptids if($suf eq 'mrna');
    } else {
      ## error, warn/loggit, skip?
      loggit(1, "ERR update_fileset.$suf empty $fin or $ftrim"); 
    } 
  }
  ## end merge loop

  my (@outfiles,@tmpfiles);
  if($upstatus == 3)  { # $upstatus == 3 or > 0?
    ## rename input files to input.old, better: input.untrim .notrim? .pretrim? .old?
    ## rename newmerge files to input
    foreach my $suf (sort keys %fset) { 
      my($nup,$nsame,$fin,$ftrim,$fup,$fupname)= @{$fset{$suf}};  
      if(-s $fupname) { rename($fupname,"$fupname.untrim"); push @tmpfiles, "$fupname.untrim"; }
      rename($fup,$fupname); push @outfiles, $fupname;
      ## push @tmpfiles, $ftrim; # from @trimfiles. dont call tmpfile? 
      loggit(0, "update_fileset.$suf upd=$nup, same=$nsame, $fin + $ftrim > $fupname"); 
    }
    runcmd("touch $flagtrimvec"); ## touch flag-file that new input has uvtrim results ..
  } else {
    loggit(1, "ERR update_fileset missing status=$upstatus/3");  # list fupnames?
    foreach my $suf (sort keys %fset) { # sort: aa,cds,mrna
      my($nup,$nsame,$fin,$ftrim,$fup,$fupname)= @{$fset{$suf}};  
      push @outfiles, $fup;
      loggit(1, "ERR update_fileset.$suf upd=$nup, same=$nsame, $fin + $ftrim >$fup"); 
    }
  }
  
  return ($upstatus, \@outfiles, \@tmpfiles, \%okids); # $upaaseq, $upcdsseq,$upmrna, 
}


sub facount {
  my $fa=shift; my $n=0;
  my($ok,$hin)= openRead($fa);
  # my $ok= ($fa =~ /\.gz$/) ? open(F,"gunzip -c $fa|") : open(F,$fa); # openRead()
  if($ok) { while(<$hin>) { $n++ if(/^>/); } close($hin); }
  return $n;  
}

sub fasize {  
  my $fa=shift; my $n=0;
  my($ok,$hin)= openRead($fa);
  # my $ok= ($fa =~ /\.gz$/) ? open(F,"gunzip -c $fa|") : open(F,$fa);# openRead()
  if($ok) { while(<$hin>) { $n += (length($_) - 1) unless(/^>/); } close($hin); }
  return $n; 
}

sub faidlist 
{
  my ($fa,$faids,$options)=@_; # other opt to return hash not file.
  my ($n,$noupdate,$ashash)=(0) x 9;  
  $options||=""; 
  $noupdate=($options=~/noupdate/)?1:0;
  $ashash=($options=~/hash/)?1:0;
  if($ashash) {
    my %faids=(); $faids= \%faids;
    my ($ok,$hin)= openRead($fa); 
    if($ok) { while(<$hin>) { if(/^>(\S+)/){ my $id=$1; $n++; $faids{$id}=1; } } close($hin); }
  } else {
    $faids= makename($fa,".ids") unless($faids); 
    return $faids if($noupdate and -s $faids);
    my ($ok,$hin)= openRead($fa);  return 0 unless($ok);
    rename($faids,"$faids.old") if(-s $faids); #? OR option to return faids as is if exists...
    $ok= open(my $hout,'>',$faids);   
    if($ok) { while(<$hin>) { if(/^>(\S+)/){ my $id=$1; $n++; print $hout "$id\n"; } } }
    close($hout); close($hin); 
  }
  return $faids;  
}

sub faheaders {  # use with parse_evgheader($id,$hdr)..
  my($fa,$headhash)=@_; 
  my $n=0; $headhash={} unless(ref $headhash); 
  my($ok,$hin)= openRead($fa);
  if($ok){ while(<$hin>){ if(/^>(\S+)/){ my $id=$1; $headhash->{$id}= $_; $n++; }} close($hin); }
  return ($headhash,$n); 
}


sub faextract #  in cdna_evigenesub.pm
{
  my ($fa,$newfa,$faids,$noupdate,$options)=@_; 
  my ($nout,$hids,%ids)=(0);  
  $noupdate||=0; $options||="";
  $newfa= makename($fa,".extract") unless($newfa); 
  return $newfa if($noupdate and -s $newfa);
  my ($ok,$hin)= openRead($fa); return 0 unless($ok);
  my $addidval= ($options=~/addidval/)?1:0;
  my $NOTid= ($options=~/notid/)?1:0;
  
  ## faids may be ref() HASH, ARRAY or filename
  if(ref($faids) =~ /HASH/) {
    %ids= %$faids;
  } elsif(ref($faids) =~ /ARRAY/) {
    %ids= map{ $_ => 1 } @$faids;
  } else {
    $faids= makename($fa,".ids") unless($faids); 
    ($ok,$hids)= openRead($faids); return 0 unless($ok); 
    while(<$hids>) { if(/^\w/) { my($id)=split; $ids{$id}=1; } } close($hids);
  }
  unless(scalar(%ids)) { 
    loggit(LOG_WARN,"ERR: faextract no IDs for $fa TO $newfa"); return($fa); 
    }
  rename($newfa,"$newfa.old") if(-s $newfa); #? OR option to return faids as is if exists...
  $ok= open(my $hout,'>',$newfa);   
  if($ok){ $ok=0; while(<$hin>) { 
    if(/^>(\S+)/){ my $id=$1; $ok= $ids{$id}||0; $ok= not $ok if($NOTid); 
      if($ok){  $nout++; s/$/ idis=$ok;/ if($addidval); } } 
    print $hout $_ if($ok); } 
  }
  close($hout); close($hin); 
  return (wantarray) ? ($newfa,$nout) : $newfa;  
}

sub fadupids #  in cdna_evigenesub.pm
{ 
  my $fa=shift; my($ndup)=(0,0); my %ids; 
  my($ok,$hin)= openRead($fa);
  if($ok) { while(<$hin>) { if(/^>(\S+)/) { my $id=$1; my $ni= ++$ids{$id}; 
    if($ni>1) { $ndup++; loggit(1,"ERR: dup id:$id in $fa"); } 
    } } close($hin); }
  # die "ERR: $ndup duplicate ids in $fa\n" if($ndup>0); # leave to caller ?
  return (wantarray) ? ($ndup,\%ids) : $ndup;
}

sub fasizes_nogap {  
  my($fa,$okayc,$gapc,$onlysize)= @_;
  $onlysize ||= 0;
  my $cokn=1; $okayc ||= 'ACGTagct'; $gapc ||= 'Nn'; 
  if($okayc =~ /^aa|amino|prot/) { $okayc='A-WYZa-wyz'; $gapc='Xx\*'; $cokn=2; }
  my ($nokay,$ngap,$nt,$id)=(0,0,0,0); # NOT total, per record
  my %fasizes=();
  my($ok,$hin)= openRead($fa);
  if($ok) { while(<$hin>) { 
    if(/^>(\S+)/) {
      if($id) { $fasizes{$id}= ($onlysize)? $nt : join("\t",$nokay,$nt,$ngap); }
      $id=$1; ($nokay,$ngap,$nt)=(0,0,0);
    } else {
      chomp; s/\*$//;
      $nt += length($_);
      if($cokn == 1) { $nokay += tr/ACGTagct/ACGTagct/; $ngap  += tr/Nn/Nn/;}
      elsif($cokn == 2) { $nokay += tr/A-WYZa-wyz/A-WYZa-wyz/; $ngap  += tr/Xx*/Xx*/; }
      else { $nokay += $nt; }
      #BAD $nokay += tr/$okayc/$okayc/; # BADD no tr/$vars/ !! tr/A-WYZa-wyz/A-WYZa-wyz/; ## aa chars; na chars=/ACGTagct/ gaps=/Nn/
      #BAD $ngap  += tr/$gapc/$gapc/;  # tr/Xx\*/Xx\*/;
    }
  } close($hin); }
  if($id) { $fasizes{$id}= ($onlysize) ? $nt : join("\t",$nokay,$nt,$ngap); }
  return \%fasizes; # ($nokay,$nt,$ngap);
}


sub fasplit { #  this splits by size; need alternate split by count: facount/ncpu
  my($fa, $spldir, $npart, $splsize, $fasuf)=@_;
  $fasuf ||= "fa";
  my @splist= (); my ($ok,$hin)= (0,undef);
  my($atsize,$atfile)=(0,0);
  my $fabase= basename($fa);
  # if($fa =~ /\.gz$/) { $ok= open(F,"gunzip -c $fa|"); } else { $ok= open(F,$fa); }
  ($ok,$hin)= openRead($fa);
  unless($ok) { loggit(LOG_DIE,"ERR: fasplit $fa $spldir $npart $splsize"); return @splist; }
  mkdir($spldir) unless(-d $spldir);
  while(<$hin>) {
    if($atfile==0 or ($atsize > $splsize && /^>/)) {
      $atfile++; if($atfile>$npart) { } # finishup???
      close(SPL) if($atfile>1);
      my $spname= $spldir . makename($fabase,".split$atfile.$fasuf",$fasuf.'|tbl|fasta|fsa|fa'); ##  choose suffix
      $ok= open(SPL,'>',$spname);  $atsize= 0;
      unless($ok) { loggit(LOG_DIE,"ERR: fasplit $atfile $spname"); return @splist; }
      push @splist, $spname;
      }
    print SPL;
    $atsize+= (length($_) - 1) unless(/^>/);
  } 
  close(SPL); close($hin);
  return @splist;
}

sub fasplitcount { # alternate split by count: facount/ncpu
  my($fa, $spldir, $npart, $splcount, $fasuf)=@_;
  $fasuf ||= "fa";
  my @splist= ();
  my($ok,$hin, $irec,$atfile)=(0) x 10; 
  my $fabase= basename($fa); #? is this bad
  # if($fa =~ /\.gz$/) { $ok= open(F,"gunzip -c $fa|"); } else { $ok= open(F,$fa); }
  ($ok,$hin)= openRead($fa);
  unless($ok) { loggit(LOG_DIE,"ERR: fasplit $fa $spldir $npart $splcount"); return @splist; }
  mkdir($spldir) unless(-d $spldir);
  
  while(<$hin>) {
    if($atfile==0 or ($irec >= $splcount && /^>/)) {
      $atfile++;  $irec=0;
      close(SPL) if($atfile>1);
      my $spname= $spldir . makename($fabase,".split$atfile.$fasuf",$fasuf.'|tbl|fasta|fsa|fa'); ##  choose suffix
      $ok= open(SPL,'>',$spname); 
      unless($ok) { loggit(LOG_DIE,"ERR: fasplit $atfile $spname"); return @splist; }
      push @splist, $spname;
      }
    print SPL;
    $irec++ if(/^>/);
    ## $atsize+= length($_) unless(/^>/);
  } 
  close(SPL); close($hin);
  return @splist;
}


##upd1804 fasort
sub fasort {
  my($infa,$outf)=@_;
  our $LOGT=($dryrun||$DEBUG)? LOG_WARN : LOG_DIE; 

  sub fasort_put {
    my($sout,$id,$hd,$fa)=@_;
    my $sid=$id;
    my($ni,$ti)= ($id=~m/(\d+)t(\d+)/)?($1,$2):(0,0); #evigene id form
    if($ni) {
      (my $gp=$id)=~s/${ni}t${ti}.*//;
      $sid="$gp\t$ni\t$ti"; 
    } else { # alphasort id
      $sid="$id\t0\t0";
    }
    map{ chomp($_); s/\n/\t/g; } ($hd,$fa);
    print $sout "$sid\t$hd\t$fa\n";
  }

  sub fasort_sort {
   my($tmpsrt,$outfa)=@_; our($LOGT);
   my $ns=0;
   open(OUT,">",$outfa) or loggit( $LOGT, "write $outfa"); 
   open(SI,"sort -T ./ -k1,1 -k2,2n -k3,3n $tmpsrt |") or loggit( $LOGT, "fail sort $tmpsrt");  
   while(<SI>) {
     my ($sa,$sb,$sc,$val)= split"\t",$_,4; 
     $val=~s/\t/\n/g;  chomp($val);
     print OUT $val,"\n"; $ns++
   }
   close(SI); close(OUT);
   return($ns);
  }
  
  my $tmpsrt= makename($infa,"$$.so"); 
  $outf= makename($infa,".sorted") unless($outf); 
  my($nt,$id,$hd,$fa)=("") x 9;  
  my($ok,$hin)= openRead($infa);
  unless($ok) { loggit( LOG_WARN, "cant sort $infa"); return -1; } 
  open(SO,">",$tmpsrt) or loggit( $LOGT, "write $tmpsrt");  
  while(<$hin>) {
    if(/^>(\S+)/){ my $d=$1; 
      fasort_put(*SO,$id,$hd,$fa) if($fa); 
      $id=$d; $hd=$_; $fa=""; }
    else { $fa.=$_; }
  } close($hin);
  fasort_put(*SO,$id,$hd,$fa) if($hd);
  close(SO); 
  
  my($ns)= fasort_sort($tmpsrt,$outf);
  unlink($tmpsrt);
  return($ns);
}

## in cdna_evigenesub.pm
sub findapp
{
  my($aname, $nodie)=@_;
  my $app=""; $nodie ||= 0;
  $app=$ENV{$aname} if(not $app and $ENV{$aname});  
  $app=$ENV{uc($aname)} if(not $app and $ENV{uc($aname)});  
  $app=`which $aname` unless($app); 
  chomp($app);
  ## #tr2aacds: app=blastn, path=no blastn in 
  my $dol=0; if(not $app or $app =~ /^no $aname/) { 
    $app="echo MISSING_$aname"; 
    $dol=($dryrun||$DEBUG||$nodie)? LOG_WARN : LOG_DIE; 
    }
  loggit( $dol, "app=$aname, path=$app");
  return (wantarray) ? ($app,dirname($app)) : $app;
}

sub finddata # new 1712
{
  my($aname, $options)=@_;
  #o# my($aname, $nodie)=@_;
  my $dfile="";  $aname||="finddata.missing"; $options||="";
  my $nodie=1; # for now;  
  my $gzok= ($options=~/gz|gunzip/)?1:0; 
  # add other suffixes like .nsq ? $sufok=($options =~ /suf=(\w+)/)?$1:0;
  $dfile=$ENV{$aname} if(not $dfile and $ENV{$aname});  
  $dfile=$ENV{uc($aname)} if(not $dfile and $ENV{uc($aname)});  
  #NO# $dfile=`which $aname` unless($dfile); 
  ## check file exists:  $isfile= -f $dfile
  ## BUT may be ncbi db prefix: UniVec == UniVec.{nsq,..} or directory
  ## also check file wildcards: *?
  my $dexists=0;
  if($dfile and -e $dfile){ $dexists=1; }  
  elsif($dfile and $gzok and -e "$dfile.gz"){ $dfile="$dfile.gz"; $dexists=1; }  
  elsif(-f $aname){ $dfile=$aname; $dexists=1; }  
  elsif($gzok and -f "$aname.gz"){ $dfile="$aname.gz"; $dexists=1; }  
  elsif($aname =~ m/[\*\?]/) {
    my $dls=`/bin/ls $aname 2>/dev/null`; chomp($dls);  # eat sys err==ls: geneval/human18ncx*.align.tab: No such file or directory
    if($dls and -e $dls) { $dfile=$dls; $dexists=1; }
    elsif($gzok) {
      $dls=`/bin/ls $aname.gz 2>/dev/null`; chomp($dls); 
      if($dls and -e $dls) { $dfile=$dls; $dexists=1; }
    }
  } 
  my $dol=0; if(not $dfile or $dfile =~ /^no $aname/) { 
    #NO# $dfile="echo MISSING_DATA_$aname"; 
    $dol=($dryrun||$DEBUG||$nodie)? LOG_WARN : LOG_DIE; 
    $dfile=""; # dont return 'no aname'
    }
  loggit( $dol, "data=$aname, exists=$dexists, path=$dfile");
  return (wantarray) ? ($dfile,dirname($dfile)) : $dfile;
}

## in cdna_evigenesub.pm
sub findevigeneapp
{
  my($aname, $nodie)=@_;
  my $app=""; 
  my $ename= basename($aname,'.pl');
  $app= $ENV{$ename} if($ENV{$ename});  
  $app= $aname if(not $app); $nodie ||= 0;
  # my $EVIGENES="$FindBin::Bin/.."; # ok?
  # my $APPcdnabest= findevigeneapp("$EVIGENES/cdna_bestorf.pl"); # allow ENV/path substitutions?
  $app="$EVIGENES/$aname" unless(-x $app);
  my $dol=0; 
  unless( -x $app) { 
    $app="echo MISSING_$aname"; 
    $dol=($dryrun||$DEBUG||$nodie)? LOG_WARN : LOG_DIE; 
    }
  loggit( $dol, "evigeneapp=$aname, path=$app");
  return($app);
}

sub cat_splitset
{
  my($tofile, $splitset)=@_;
  my $nin= (ref($splitset))? @$splitset : 0; 
  loggit( ($dryrun) ? 1 : 0,"CMD=","cat_splitset to $tofile from n=$nin,",$$splitset[0],".. ");  
  return 0 if($dryrun); return -1 if($nin<1);
  my (@qfail,@ofail);
  if(-s $tofile) { } # save?  
  my $nok= 0;
  my $ok= open(O,'>', $tofile);
  if($ok) { for my $i (0..$nin-1) {
    my $outf= $splitset->[$i];
    my $isend= 0;
    if($outf and -s $outf) {
      my $lastl="";      
      $ok= open(I,$outf); while(<I>) { $lastl=$_; print O $_; } close(I);
      $isend=$ok; #not blast# $isend++ if( $lastl =~ m/^# BLAST processed/); 
    }
    if($isend>0) { $nok++; } else { push @ofail, $outf; }
  }
  close(O); }
  my $runerr= $nin - $nok;
  if($runerr) { loggit(1,"ERR=$runerr ","cat_splitset to $tofile"); } 
  return ($runerr, $nok, \@ofail);
}

sub runcmd
{
  my @cmd= @_;
  ## fail if $cmd[0] =~ /MISSING_/
  my $isdry=($dryrun or $cmd[0] =~ m/^echo\b/)?1:0;
  loggit( ($isdry) ? 1 : 0,"CMD=",@cmd);  
  my $err= ($isdry) ? 0 : system(@cmd);  ## ?? add option run: @result=`$cmd`;
  if($err) { loggit(1,"ERR=$err ",$cmd[0]); } # ..
  return $err;
}

sub forkcmd
{
  my @cmd= @_;
  loggit( ($dryrun) ? 1 : 0,"forkCMD=",@cmd);  
  unless($dryrun) {
    my $pid= fork();
    if($pid) { # parent
      return $pid;
    } else { # child
      my $err= system(@cmd);
      exit($err);
    }
  }
  # if($err) { loggit(1,"ERR=$err ",$cmd[0]); } # ..
}

sub makename
{
  my($infile,$osuf,$insuf)=@_;
  #BAD?#  $insuf ||= 'aa|blast|cdna|mrna|cds|tr|trclass|tbl|fasta|faa|fsa|fa';  ## fixme need insuf: tr|fasta|fa
  ## in: dmagset56ri.tr.split4.fa : OUT dmagset56ri_split/dmagset56ri.aa  LOST split4; all splits to one.aa
  unless($infile) { warn "cdna_evigenesub:makename MISSING infile"; $infile="Noname$$"; } # or die / buggy soft
  my $outfile= $infile; $outfile =~ s/\.gz$//; # bad for infile empty/undef .. do what?

  #BADold# $outfile =~ s,\.($insuf)[^\/\s]*$,,; ## ADD \. cut only final suffix !!!
  ## or instead use: smallest at of $at=rindex($outfile,'.'.$insuf[i]); $outfile=substr($outfile,0,$at)
  #o2# $insuf ||= 'aa|blast|cdna|mrna|cds|fasta|fa';  ## fixme need insuf: tr|fasta|fa
  #o2# $outfile =~ s,\.($insuf)[^\.\/\s]*$,,; ## FIXMEd, or s,\.($insuf)\w*$,,; or s,\.($insuf)$,,;
  
  if($insuf) { $outfile =~ s,\.($insuf)[^\.\/\s]*$,,; } # insuf can have \.: 'aa.qual|aa.size|...'
  else { $outfile =~ s,\.\w*$,,;  } ## default insuf maybe chop any \.\w+$ instead of guess list?
  
  $outfile.= $osuf if($osuf); 
  $outfile.= "_out" if($outfile eq $infile);
  return $outfile;
}

=item evigene_config

	GetOptions(.. "config=s", \$configfile, "cadd=s", \@configadd,);
	$configh= evigene_config($configfile, \@configadd); # always even if $config null

	$DBID= $configh->{general}->{dbid} || "MyDBID"; 
	$LOCUSTAG= $configh->{general}->{locus_tag} || $DBID; 
	$IDPrefix= $configh->{pubopt}->{publicid} || "Evigene"; 
	
 	from  evigene/scripts/evigene2genbanktbl.pl
	but drop global special hashes: %evidence; %public_options; %geneset; %programs; 
	change to: return \%config # replace other config hashes w/ this 2-level: config{part}{key} = val

=cut

sub evigene_config {
  my($cfile, $addoptions)= @_;
  my %config=(); my $ctype=0;

  use constant{ kGENERAL => 'general', kEVFILE => 'evidence', kEVPROG => 'programs', kPUBOPT => 'pubopt',  };
  	#old: kEVOPT => 'evoption', kANOPT => 'anoption',  kEVGENES => 'geneset',
  	
  if($cfile) { #  and -f $cfile
    open(F,$cfile) or die "ERROR reading config: $cfile";
  
    my @CONFIG= <F>;
    push @CONFIG, @$addoptions if(ref $addoptions);
    
    my ($lastkey, $lastval);
    foreach (@CONFIG)
    {  # key => value
      s/^\s+//;  s/\#.*$//; s/\s+$//; # end of line comments 

    ## need now to handle continuation lines, end with \
    ## or prefix with "+" << use this

      my($key,$val);
      if(/^[+\.]/ and $lastkey) { $key= $lastkey; s/^.//; $val= $lastval.$_; }
      elsif($lastval =~ m,\\$, and $lastkey) { $key= $lastkey; $lastval=~s,\\$,,; $val= $lastval.$_; }
      else { ($key,$val)= split(/[=\s]+/, $_, 2); }
      
      next unless($key =~ /^\w/); 
      $val =~ s/\s*$//;  $val =~ s/\\n/\n/g; $val =~ s/\\t/\t/g; # dequote tabs,newlines
      # allow for val == regex in match/m or subs/s;  
      #  names:
      #   cutdbx = s/\((InterPro|TAIR):[\w\.-]+\)//g ?
      #   isgeneid = m/^(Os|At|AT)\d{1,2}g\d{3,}/

# revise to parse '^section:' into separate hash
      if($key =~ s/:$//) { $ctype=$key; }
      
      if($key =~ /^evidence/) { $ctype= kEVFILE; } # old/fixed set
      elsif($key =~ /^pubopt/) { $ctype= kPUBOPT; }
      elsif($key =~ /^program/) { $ctype= kEVPROG; }
      elsif($key =~ /^end$/) { $ctype= 0; }
#       elsif($key =~ /^evoption/) { $ctype= kEVOPT; }
#       elsif($key =~ /^anoption/) { $ctype= kANOPT; }
#       elsif($key =~ /^geneset/) { $ctype= kEVGENES; }
#.. drop special config hashes...
#       elsif($ctype eq kEVFILE ) { $evidence{$key}= $val; } # /gff/ was bad for geneset
#       elsif($ctype == 0 and $val =~ /\.gff/) { $evidence{$key}= $val; } 
# #       elsif($ctype eq kEVOPT) { $evaluate_options{$key}= $val; } 
# #       elsif($ctype eq kANOPT ) { $annotate_options{$key}= $val; } 
# #      elsif($ctype eq kEVGENES) { $geneset{$key}= $val; }
#       elsif($ctype eq kPUBOPT ) { $public_options{$key}= $val; } 
#       elsif($ctype eq kEVPROG) { $programs{$key}= $val; }

      # generic keys: name date genome .. other?
      if($key =~ /\w/ and $val =~ /\w/) { 
        my $ogroup= $ctype || "general";
        $config{$ogroup}->{$key}= $val; # this one #which?  $config{$ogroup}{$key}= $val;
      }
      
      # also for now : overlap, other progs
      ($lastkey, $lastval)=($key, $val);
    } close(F);
  }
  return \%config;
}


1;
__END__
