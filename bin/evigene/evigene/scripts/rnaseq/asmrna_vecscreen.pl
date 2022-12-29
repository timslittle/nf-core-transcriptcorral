#!/usr/bin/env perl
# evigene/scripts/rnaseq/asmrna_vecscreen.pl

=item about asmrna_vecscreen

  update for asmrna_trimvec, 2020apr
  
  update vecscreen() portion of evigene/scripts/rnaseq/asmrna_trimvec.pl
  a. replace NCBI vecscreen with blastn + parsing (vecscreen binary no longer widely available,
    same results using blastn -db UniVec +ncbi-vecscreen-options | vecscreen_parse_blastn
  b. add blastn -db contam_in_euks.fa 
  c. optional blastn -db other_contam.fa, eg. -db foreign_genes.fa

  Note 2020mar: UniVec 2017 latest at ncbi, and so is contam dbs at ftp://ftp.ncbi.nih.gov/pub/kitts/

=item UPD20apr example 4d

  these calls *should* work
  
  pt=evg4567hetfix; 
  
  A.... screen mrna,ncrna on 2 contam dbs
  $evigene/scripts/rnaseq/asmrna_vecscreen.pl -merge=$pt.allvec.vtab -log -debug -NCPU 3 \
    -mrna okayset/$pt.okay.mrna.gz -ncrna ncrnaset/$pt.ncrna_pub.fa  \
    -vectordb $contamf/evg_uvncbi.fa -contamdb $contamf/contam_in_euks.fa -rrnadb no.rrna.fa  
  
  # -merge=$pt => output to $pt.allvec.vtab .. change to fullname -merge=$pt.allvec.vtab ?

  B.... update mrna,ncrna seq sets w/ trimmed seqs, or -defer make uvcut set only
  
  $evigene/scripts/rnaseq/asmrna_trimvec.pl  -deferupdate -MINSIZE 180 -NCPU 2 -log \
    -names $pt.names -vectors $pt.allvec.vtab -mrna okayset/$pt.okay.mrna.gz -ncrna ncrnaset/$pt.ncrna_pub.fa
  
  #* add -ncbipolicy, not default, MINSIZE=200, ENDGAP2 = 50
  # test -nodeferupdate, wont handle ncrna update yet. is mrna update ok still?
  $evigene/scripts/rnaseq/asmrna_trimvec.pl  -nodeferupdate -ncbipolicy -MINSIZE 180 -NCPU 2 -log \
    -names ixos9fevgok.names  -vectors vectrimset/$pt.allvec.vtab \
    -mrna okayset/$pt.okay.mrna.gz -ncrna ncrnaset/$pt.ncrna_pub.fa
  
=cut

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/.."); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname);
use cdna_proteins;
use cdna_evigenesub;  

use constant VERSION =>  '2020.04.04'; # 4v

our $EVIGENES= $FindBin::Bin; $EVIGENES =~ s,scripts/.*,scripts,;
our $EGAPP='egvecscreen';  
our $EGLOG='vecs';
our $DEBUG= $ENV{debug}||0;  # set package dbg

my $NCPU=$ENV{ncpu}||1; 
my $tidyup=1;
my $MINSIZE_NCBI=200; # NCBI
my $MINSIZE= $ENV{min} || 150; # lower to check cuts w/ valid homology?
my $MAXGAP=$ENV{maxgap}|| 15; # NCBI decides, changes..

use constant { tyVECTOR => 'vector', tyCONTAM => 'contam', tyRRNA => 'rrna', tyFOREIGN => 'foreign', };
my $UniVecDB= $ENV{UniVec} || "UniVec_Core";  # other?
my $contamDB= $ENV{contamdb} || "contam_in_euks";  
my $foreignDB=$ENV{foreigndb} || "";
my $rrnaDB=$ENV{rrna} || "";

my $VECSUF= "vector.tab"; #fixed in  vecscreen: makename($cdnaseq,".vector.tab");
my $outsuf= $ENV{outsuf}|| "uvcut"; # was uvcut.mrna; 

our($dryrun);
my (@mrna,$logfile,$UPDATEALL,$merge,$nameout);
my @ARGSAVE=@ARGV;
my %opts= (
  "mrna|ncrna|cdna|input=s", \@mrna, # require, maybe several files
  "logfile:s", \$logfile, # option
  "vectordb|UniVecDB=s", \$UniVecDB,  
  "contamDB=s", \$contamDB, "foreignDB=s", \$foreignDB, "rrnaDB=s", \$rrnaDB,
  #"MIN_NCRNA=i", \$MIN_NCRNA,  
  #"MIN_CDSALN=i", \$MIN_CDSALN,  
  "NCPU=i", \$NCPU, # "MAXMEM=i", \$MAXMEM,  
  "UPDATEALL!", \$UPDATEALL, 
  # "merge!", \$merge, "nameout=s", $nameout, # only for merge .. use -merge=nameout ?
  "merge:s", \$merge,  
  "debug!", \$DEBUG,  "dryrun|n!", \$dryrun, 
  );

my $optok= GetOptions(%opts);
my $optlist= ($DEBUG) ? "All opts: ".join(", ",sort keys %opts) : "";

die "EvidentialGene asmrna_vecscreen VERSION ",VERSION,"
mRNA vector screen/removal and gap trimming
Usage: asmrna_vecscreen.pl -mrna mrna.fa -ncrna ncrna.fa  -NCPU $NCPU -debug
  screening databases: -vectordb $UniVecDB -contamdb $contamDB -rrna $rrnaDB -foreign $foreignDB 
" unless($optok and (@mrna));     

my ($APPblastn,$APPmakeblastdb);
my ($pubsetd,$mrna,$trname); 
our( @tmpfiles, @vecset, @erasefiles);
#  my $tmpfolder="vectrimset"; # @tmpfiles go here
#  my $outfolder= $pubsetd; # @vecset go here ??

MAIN_vecscreen();

#================================================================

sub MAIN_vecscreen {

  $mrna= $mrna[0]; #must have, maybe in subdir
  $trname=$mrna; $trname =~ s/\.gz$//; $trname=~s/.\w+$//; $trname =~ s,^.*/,,;
  
  openloggit($logfile,$trname); # ||$trclass
  loggit(0, "EvidentialGene asmrna_vecscreen.pl VERSION",VERSION);
  loggit(0, "asmrna_vecscreen.pl",@ARGSAVE);
  loggit(0, "BEGIN with input=",@mrna,"date=",`date`);
  $APPblastn= findapp("blastn");  #die if missing? 
  $APPmakeblastdb= findapp("makeblastdb");  

  if(1) { #?? dont need this, assume @mrna is proper file list
    $pubsetd="publicset";
    $pubsetd="okayset" unless( -d $pubsetd ); 
    # loggit(LOG_WARN,"Cannot find evigene data directories:  $pubsetd") unless( -d $pubsetd );
    # locate data associated w/ trclass  
    $mrna= filegz("$pubsetd/$trname.mrna_pub.fa") unless(-f $mrna);   
    $mrna= filegz("$pubsetd/$trname.okay.mrna") unless(-f $mrna); 
    $mrna= filegz("$pubsetd/$trname.mrna") unless(-f $mrna);   
    
    loggit(LOG_DIE,"Cannot find mrna $mrna") unless(-f $mrna);
    loggit(LOG_DIE,"Cannot find blastn: $APPblastn") if($APPblastn =~ /MISSING/);
    
    ##? relocate other rnas ?? ugh, assume @mrna is correct
    # if(@mrna>1 and $mrna ne $mrna[0]) { 
    #   my $si= index($mrna,$trname);
    #   my $pre= substr($mrna,0,$si+length($trname));
    #   for my $j (0..$#mrna) { 
    #     my $si= index($mrna[$j],$trname);
    #     my $suf= substr($m,$si+length($trname));
    #     my $m= $pre.$suf; $mrna[$j]=$m if(-f $m);
    #     }
    #   }  
    
    $mrna[0]= $mrna;
  }
  
  # FIXME: output names should be $trname.vecNNNN.vtab or $mrna.vecNNN.vtab
  #  for consistancy w/ other evg trname set data.
  
  my($nvout,$vtab)= screen_vectors($UniVecDB,\@mrna);
  my($ncout,$ctab)= screen_contams($contamDB,\@mrna);
  my($nrout,$rtab)= screen_rrna($rrnaDB,\@mrna);
  my($nfout,$ftab)= screen_foreign($foreignDB,\@mrna);
  
  if(defined $merge and $nvout+$ncout+$nrout+$nfout > 0) {
    #o my $oname="vall-$trname.vtab"; # what name?
    #? my $oname="$trname-vall.vtab"; # what name?
    #? my $oname= $vtab || $ctab || $rtab || $vtab; $oname =~ s/\-\w+/-vall/;
    #* Only here use $nameout == $merge name, else @mrna will overwrite blastn/vtabs
    my $oname;
    if($merge and $merge =~ m/.vtab$/) { $oname=$merge; } # full table name
    else { $oname= outname( $merge || $mrna, 'allvec', 'vtab'); }
    
    my $cmd="sort -k1,1 -k2,2n -k3,3nr -k4,4r -k5,5";
    $cmd.=" $vtab" if($nvout);
    $cmd.=" $ctab" if($ncout);
    $cmd.=" $rtab" if($nrout);
    $cmd.=" $ftab" if($nfout);
    $cmd.=" > $oname";
    my $err= runcmd($cmd);
    push @tmpfiles,$oname;  
  }
  
  FINISH_vecscreen(0,'done okay',$trname);
}

sub blast2vectab {
  my($bldb,$vtype,@bltab)= @_;
  # blast output to simpler vector table,
  # qname-evg_contam.blastn > qname-evg_contam.vtab
  # assume for now   $outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" ;

  my($nout,$otab,$ltd)= (0,0,0);
  return($nout,$otab) unless(@bltab and -s $bltab[0]);

  ($otab=$bltab[0]) =~ s/\.\w+$/.vtab/;
  # maybe fixme here: bltab = bob.mrna.blastn + bob.ncrna.blastn > bob.vtab ?
  if( -s $otab and not $UPDATEALL and (-M $otab <= -M $bltab[0])) {  # -M $bltab newer than $otab
    $nout=`wc -l $otab`; chomp($nout);  return($nout,$otab); 
    } 

  # $VDB ||= ""; # prefix for vector db
  # my $vn= $bldb; $vn =~ s/\.\w+$//; $vn =~ s,^.*/,,;
  my $VDB= $vtype.":"; # "v$vn:";
  my $MIND=50; # bp distance for diff hits
  my $TERMD=40; # bp from ends for opts
  
  # FIXME: need option variants for vectors, contams, foreign 
  my $Smin=$ENV{Smin}||"24,30";  # Strong, End/Inner raw-score min
  my $Mmin=$ENV{Mmin}||"19,25";  # Moderate, End/Inner raw-score min
  my($Smint,$Smini)= split",",$Smin; 
  my($Mmint,$Mmini)=split",",$Mmin;  
  if($vtype eq tyCONTAM) {
    # mins: ($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200) where $3 == $pi, $4 == $aln
  } elsif($vtype eq tyVECTOR) {
  
  } elsif($vtype eq tyRRNA) {
  
  } elsif($vtype eq tyFOREIGN) {
  
  }
  
  our @ovr=();
  sub nover { 
    my($tb,$te)=@_; my $ov=0; our @ovr;
    for my $be (@ovr) { 
      if($tb < $$be[1] and $te > $$be[0]) {
        return 0 if($tb >= $$be[0] and $te <= $$be[1]); 
        my $db=$$be[0] - $tb; my $de=$te - $$be[1];
        return 0 unless($db >= $MIND or $de >= $MIND);  
      } } 
    return 1; 
  } 

  open(my $outh,'>',$otab) or return -1; # loggit
  for my $bltab (@bltab) { 
    # open(my $inh, $bltab) or next; # loggit warn
    my($oki,$inh)= openRead($bltab); next unless($oki);
    while(<$inh>) {
      next if(/^\W/); 
      my @v=split; 
      my($td,$vd,$pi,$aln,$mm,$nd,$tb,$te,$vb,$ve,$evl,$bs,$tlen,$vlen)=@v; #outfmt 6 + tlen,vlen
      # next if($OK and $ok{$td});
      my($smin,$mmin,$sclass)=(0,0,"Weak");
      my $sc=int($bs/2); # ~raw-score = bitscore/2
      $tlen= $te unless($tlen); # error .. or use sizes table? rna hdr size?

#FIXME: extend vecspan to start/end if (a)  TERMD-near, (b) have weak match at end (both/either?) 
use constant DOXTND => 1;

      if($vtype eq tyCONTAM) {
        ($smin,$mmin)=(700,50); #? raw score =~ align bases
        #  $smin=($sc >= 700 or $sc >= 0.75 * $tlen)?$sc : 
        $sclass=($sc >= 700 or $sc >= 0.75 * $tlen)?"Strong":"Moderate"; # add Weak?
        next unless($pi>=98.0 && $aln>=50)||($pi>=94.0 && $aln>=100)||($pi>=90.0 && $aln>=200);
      } elsif($vtype eq tyVECTOR) {
        ($smin,$mmin)=($tb < $TERMD or $te > $tlen - $TERMD)? ($Smint,$Mmint): ($Smini,$Mmini); 
        next if($sc < $mmin); 
        $sclass=($sc >= $smin)?"Strong":"Moderate"; # add Weak?
      } elsif($vtype eq tyRRNA) {
        ($smin,$mmin)=(100,100); 
        next unless($aln>=100);
        $sclass=($sc >= $smin)?"Strong":"Moderate"; # add Weak?
      } elsif($vtype eq tyFOREIGN) {
        ($smin,$mmin)=(200,120); #? raw score =~ align bases
        next unless($pi>=98.0 && $aln>=120);
        $sclass=($sc >= $smin)?"Strong":"Moderate"; # add Weak?
     }
      
      @ovr=() if($td ne $ltd);
      if( nover($tb,$te) ) { 
        if(DOXTND and ($tb < $TERMD or $te > $tlen - $TERMD)) { 
          $tb=1 if($tb < $TERMD);
          $te=$tlen if($te > $tlen - $TERMD); # flag it?
        }
        (my $vdp=$vd)=~s/^gi.\d+.//; 
        $vdp=~s/^gnl\|(\w+)\|/$1:/; $vdp=~s/(\w+)\|/$1:/; 
        push @ovr,[$tb,$te];
        #above# $sclass=($sc >= $smin)?"Strong":"Moderate"; # add Weak?
        print $outh join("\t",$td,$tb,"$te/$tlen",$sclass."_match","$VDB$vdp")."\n"; $nout++;
        } 
      $ltd=$td;
    } close($inh); 
  } close($outh);
  
  push @tmpfiles,$otab; # or @vecset ?
  return($nout,$otab);
}

sub hasdb { return($_[0] and not ($_[0] =~ /^(no|0)/)); }

sub outname { 
  my($rna,$dbname,$suf)=@_; 
  my $oname= basename($rna); 
  $oname =~ s/\.gz//;
  $oname="$oname-$dbname.$suf";
  return $oname; 
}


sub makeblastdb {
  my($vdb)= @_;
  my($cmd,$err);
  my $bldb= finddata("$vdb*nsq") || finddata("$vdb*");
  if($bldb =~ /\.nsq$/){ $bldb =~ s/.nsq//; }
  elsif(-f $bldb) {
    if($vdb =~ /.gz/) { 
      (my $vdbnam= $vdb) =~ s/.gz//;
      $cmd="gunzip -c $vdb | $APPmakeblastdb -dbtype nucl -title $vdbnam -out $vdbnam -logfile $vdbnam.mblog"; 
      $bldb= $vdbnam;
    } else {
      $cmd="$APPmakeblastdb -dbtype nucl -in $vdb -logfile $vdb.mblog"; 
      $bldb= $vdb;
    }
    $err= runcmd($cmd);  loggit(LOG_DIE,$cmd) if($err);
  } else {
    loggit(LOG_WARN,"Cannot find data: $vdb"); return 0;
  }
  return($bldb);
}

sub blast1rna {
  my($rna,$blopt,$bldb,$outfmt)=@_;
  # ($oname,$err,$blcmd)= blast1rna($rna,$blopt,$bldb,$outfmt);
  my($cmd,$err);
  $outfmt ||= "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" ;

  my $dbn= $bldb; $dbn =~ s/\.\w+$//; $dbn =~ s,^.*/,,; $dbn="v".$dbn;
  # my $oname= basename($rna); 
  # $oname =~ s/\.\w+$//; #  dont cut suffix here, eg. bob.mrna, bob.ncrna will replace blastn
  $cmd= "$APPblastn $blopt -num_threads $NCPU -db $bldb -outfmt '$outfmt'";
  if($rna =~ m/\.gz$/){ $cmd= "gunzip -c $rna | $cmd"; } #  $oname =~ s/\.gz//;  -out $oname
  else { $cmd="$cmd -query $rna"; } #  -out $oname
  
  # FIXME here: change output names to same trname-xxxx.sss as other evg data
  # $oname="$oname-$dbn.blastn"; #old: $oname="$dbn-$oname.blastn"; 
  my $oname= outname($rna, $dbn, 'blastn');
  $cmd .=" -out $oname";
  if( -s $oname and not $UPDATEALL) { $err=0; } 
  else { $err= runcmd($cmd); loggit(LOG_DIE,$cmd) if($err); } #? return if($err);
  return($oname,$err,$cmd);
}

sub blastrnas {
  my($vdb,$vtype,$mrnas,$blopt,$outfmt)=@_;

  # if($vtype eq tyCONTAM) {
  # } elsif($vtype eq tyVECTOR) {
  # } elsif($vtype eq tyRRNA) {
  # } elsif($vtype eq tyFOREIGN) {
  # }

  my @blout=();
  my $bldb= makeblastdb($vdb) or return;
  for my $rna (@$mrnas) {
    my($oname,$blerr,$blcmd)= blast1rna($rna,$blopt,$bldb,$outfmt);
    if(-f $oname) { push @blout,$oname;  push @tmpfiles,$oname; }
  }
  return(@blout);
}

sub screen_vtype {
  my($vtype, $vdb, $mrnas, $blopt, $outfmt)= @_;

  # ??which ofmt .. same for all?
  # $outfmt="6 qseqid sseqid pident length qstart qend qlen sstart send slen score" ;
  $outfmt ||= "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" ;

  return (0) unless( hasdb($vdb));
  my @blout= blastrnas($vdb, $vtype, $mrnas, $blopt, $outfmt);
  my($nout,$vecout)= blast2vectab($vdb, $vtype, @blout);
  loggit(0, "screen $vtype n=$nout in $vecout from $vdb");
  return($nout,$vecout);
}

=item screen_vectors

$nbin/blastn -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true \
 -evalue 700 -searchsp 1750000000000  \
 -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen score"  \
 -db evg_contamuvncbi.fa  -query $queryfa -out qname.evg_contam.blastn

=cut

sub screen_vectors {
  my($vdb, $mrnas)= @_;
  my $blopt="-task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000";

  return screen_vtype(tyVECTOR,$vdb, $mrnas, $blopt);

  #   my $outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" ;
  #   return (0) unless( hasdb($vdb));
  #   my @blout= blastrnas($vdb, tyVECTOR, $mrnas,$blopt,$outfmt);
  #   my($nout,$vecout)= blast2vectab($vdb, tyVECTOR, @blout);
  #   loggit(0, "vectors $vdb n=$nout, $vecout");
  #   return($nout,$vecout);
}

=item screen_contams

#  this way for these longish contams, vs vecscreen contam opts
$nbin/blastn -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 \
 -dust yes -evalue 0.0001 -perc_identity 90.0 \
 -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen score"  \
 -db contam_in_euks.fa -query $queryfa | \
 awk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' \
 > qname.contam_in_euks.blastn

=cut

sub screen_contams {
  my($vdb, $mrnas)= @_;
  my $blopt="-task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1"
    ." -dust yes -evalue 0.0001 -perc_identity 90.0";

  return screen_vtype(tyCONTAM,$vdb, $mrnas, $blopt);
  
  #   return (0) unless( hasdb($vdb));
  #   
  #   # ??which ofmt
  #   # $outfmt="6 qseqid sseqid pident length qstart qend qlen sstart send slen score" ;
  #   my $outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" ;
  #   
  #   # my $bldb= makeblastdb($vdb);
  #   # for my $rna (@$mrnas) {
  #   #   my($oname,$blerr,$blcmd)= blast1rna($rna,$blopt,$bldb,$outfmt);
  #   #   push @blout,$oname; ## unless($blerr);
  #   # }
  # 
  #   my @blout= blastrnas($vdb, tyCONTAM, $mrnas,$blopt,$outfmt);
  #   my($nout,$vecout)= blast2vectab($vdb, tyCONTAM, @blout);
  #   loggit(0, "contam $vdb n=$nout, $vecout");
  #   return($nout,$vecout);
}

sub screen_foreign {
  my($vdb, $mrnas)= @_;
  
  #?? what opts? like contam?  bump perident to high value: 98 ?
  my $blopt="-task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1"
    ." -dust yes -evalue 0.0001 -perc_identity 98.0";
  
  return screen_vtype(tyFOREIGN,$vdb, $mrnas, $blopt);
}

=item screen_rrna  

$nbin/blastn -task megablast -template_length 18 -template_type coding -window_size 120 -word_size 12 -xdrop_gap 20 \
 -no_greedy -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 1E-9 -gapextend 2 -gapopen 4 \
 -penalty -4 -perc_identity 95 -reward 3 -soft_masking true \
 -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen score"  \
 -db rrna -query $queryfa | awk '$4>=100' > qname.rrna.blastn

=cut

sub screen_rrna {
  my($vdb, $mrnas)= @_;

  my $blopt="-task megablast -template_length 18 -template_type coding -window_size 120 -word_size 12 -xdrop_gap 20"
    ." -no_greedy -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 1E-9 -gapextend 2 -gapopen 4"
    ." -penalty -4 -perc_identity 95 -reward 3 -soft_masking true";
    
  return screen_vtype(tyRRNA,$vdb, $mrnas, $blopt);
}


sub filegz{ my($f)=@_; unless(-f $f){ $f="$f.gz" if(-f "$f.gz"); } return $f; }
# sub step_log_or_fini{ my($fini,$msg)=@_; if($fini){ FINISH_vecscreen(-1,$msg); } else { loggit( 0, $msg); } }

sub FINISH_vecscreen {
  my($status,$atstep,$trname)= @_;
  
  my $tmpfolder="vectrimset";
  my $outfolder= $pubsetd; # source files are here, need to merge reorfiles still
  loggit( LOG_WARN,"FAILED at step:",$atstep) if($status<0);
  
  if( $tidyup ) {  #  and $nupdate > 0 
    tidyupFileset($outfolder,@vecset) if(@vecset);  # not used so far
    tidyupFileset($tmpfolder,@tmpfiles) if(@tmpfiles);  # all data here so far ncrnaset/
    my @rmlist;
    if($status<0) { 
      tidyupFileset("vecset_failtmp",@erasefiles) if(@erasefiles);  ## dont erase
    } else {
      foreach my $fn (@erasefiles) { if(-f $fn) { unlink($fn); push @rmlist,$fn; } }  
      if(@rmlist) { my $nrm=@rmlist; my $rml=join" ",@rmlist[0..4]; loggit(LOG_DEBUG,"tidyup erase: n=$nrm, $rml .."); }
    } 
    
  }
  loggit(0, "DONE at date=",`date`);
  loggit(0, "======================================"); # log may append runs
  exit($status);
}

__END__

=item NCBI vecscreen

VecScreen Match Categories
Vector contamination usually occurs at the beginning or end of a sequence; therefore, different criteria are applied for 
terminal and internal matches. VecScreen considers a match to be terminal if it starts within 25 bases of the beginning of 
the query sequence or stops within 25 bases of the end of the sequence. Matches that start or stop within 25 bases of another
match are also treated like terminal matches. Matches are categorized according to the expected frequency of an alignment 
with the same score occurring between random sequences.

Strong Match to Vector
(Expect 1 random match in 1,000,000 queries of length 350 kb.)
    Terminal match with Score <E2><89><A5> 24.
    Internal match with Score <E2><89><A5> 30.
Moderate Match to Vector
(Expect 1 random match in 1,000 queries of length 350 kb.)
    Terminal match with Score 19 to 23.
    Internal match with Score 25 to 29.
Weak Match to Vector
(Expect 1 random match in 40 queries of length 350 kb.)
    Terminal match with Score 16 to 18.
    Internal match with Score 23 to 24.
Segment of Suspect Origin
    Any segment of fewer than 50 bases between two vector matches or between a match and an end. 

=item adaptors_for_screening_euks

ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/adaptors_for_screening_euks.fa

qseq=geneset.mrna
$nbin/vecscreen -text_output -outfmt 0 -db $contamf/adaptors_for_screening_euks.fa \
 -query $qseq -out $qseq.vecscreen.tmp > & $qseq.vecs.log
env ncbicxx=1 vecscreen=$nbin/vecscreen $evigene/scripts/rnaseq/asmrna_trimvec.pl \
 -parsevec $qseq.vecscreen.tmp

=item contam_in_euks

ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_euks.fa.gz

qseq=geneset.fa
$nbin/blastn -query $qseq  -db $contamf/contam_in_euks.fa \
 -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 \
 -dust yes -evalue 0.0001 -perc_identity 90.0 \
 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | \
 awk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' \
 > $qseq.contam_in_euks.blastn

Contaminant matches from (1)  are merged if they are from the same class of
sequence (VECTOR, E.coli, IS, PHG) and they overlap or are separated by 50 bases
or less.

a. If the total coverage of contaminant matches from (1) is >75% of the sequence
length then flag the sequence as a contaminant to be excluded.

b. If the contaminant is classed as VECTOR, E.coli, IS:.*, PERM:.* or PHG:* and
the contaminant location is within 100 bases of the the start or end of the
sequence (or gap is the sequence is not contiguous), or within 100 bases of
another contaminant match that is at an end, flag the contaminant span for
trimming.

If the contaminant is one of the above, and the match is longer than 700 bases
flag the contaminant span for trimming. .. Other matches may be false alarms.
Treat them as suspect spans and reBLAST the hit span plus 10 Kbp of flanking
sequence on each side against nr, HTGS, related and unrelated chromosomes (as
described below).

Flag all adaptor spans for trimming.

=item rrna

ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/rrna.gz

qseq=geneset.fa

$nbin/blastn -query $qseq -db $contamf/rrna \
 -task megablast -template_length 18 -template_type coding -window_size 120 -word_size 12 -xdrop_gap 20 \
 -no_greedy -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 1E-9 -gapextend 2 -gapopen 4 \
 -penalty -4 -perc_identity 95 -reward 3 -soft_masking true -outfmt 6 | \
  awk '$4>=100' > $qseq.rrna.blastn

Ribosomal RNA genes are the cause of many false positives because the include some segments that align to distantly
related organisms. Segments that match rRNA genes are identified so that such segments are not reported as being foreign.


=cut
