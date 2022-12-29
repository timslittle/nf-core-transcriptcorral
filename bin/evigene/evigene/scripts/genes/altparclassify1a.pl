#!/usr/bin/env perl
# altparclassify.pl
# cut from  evigene altpars/reclassaltpar.pl

use FindBin;
use lib ("$FindBin::Bin"); # assume evigene/scripts/genes; layout: main, prot/, rnaseq/, ...

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname);
#? use cdna_evigenesub;  

# our $EVIGENES= $ENV{evigenes} || "$FindBin::Bin";  # fixme for other loc
our $EVIGENES= $FindBin::Bin; $EVIGENES =~ s,scripts/.*,scripts,;

our $EGAPP='altparclass';  
our $EGLOG='apc';
our $dryrun=0; ## $DRYRUN ?
our $DEBUG= $ENV{debug}|| 0;

my $CDSBLAST_IDENT=  $ENV{pctident} || 95; # option
my $CDSBLAST_EVALUE= $ENV{evalue} || 1e-5; # option
my $NCPU=$ENV{ncpu}||1; 
my $MAXMEM=$ENV{maxmem}||8000; # in Mb not used
my $MINA=  $ENV{mina}||90; 
my $ALTPI= $ENV{altpi} || 99.0; #hi qual genes use 99.990; # or use mism <= 1/2 , ndel <= 0/1
my $PDIGITS = $ENV{digits}||3;
my $IDPRE= $ENV{idprefix} || "dg"; #?
my $GID  = $ENV{geneidnum} || 0; # allow user opt, ie to add to evg pubids

my($cdsseq,$logfile,$runname,$cdsselfblast);
my($genesizes, $asCDSsize, $idset,$requiresize)= (0,1, 0, 0);

my $optok= GetOptions(
  "cdsseq=s", \$cdsseq,   
  "sizes|qual=s", \$genesizes,   
  "selfblast=s", \$cdsselfblast, # optional input self.blast table
  "requiresize|needsize!", \$requiresize, 
  # "runname=s", \$runname,   
  "logfile:s", \$logfile,
  "pctident=s", \$CDSBLAST_IDENT,  
  "minalign=s", \$MINA,  
  "altpident=s", \$ALTPI,  
  "idprefix=s", \$IDPRE,  
  "geneidnum|GID=i", \$GID,  
  "NCPU=i", \$NCPU, # "MAXMEM=i", \$MAXMEM,  
  "dryrun|n!", \$dryrun, 
  "debug!", \$DEBUG, 
 ); 

$cdsseq= shift(@ARGV) unless($cdsseq);
die "# usage: altparclassify  $cdsseq\n opts: -ncpu 8 -log" unless($optok and $cdsseq);

my($nokids, $idokset)=(0,0);
my($nsizes,$aasizeh,$trsizeh,$cdspanh,$aaqualh)= (0,0,0,0);
my $ALTAL= ($MINA>59) ? $MINA : 90; # ENV{altal}

sub MAIN_STUB {}

  # openloggit($logfile,$runname);  
  # loggit(1, "EvidentialGene altparclassify.pl"); ##  "(-help for info), VERSION",VERSION);

  (my $cname= $cdsseq) =~ s/\.\w+$//;
  my $genetab= "$cname.dgclass.tab"; # "refarath16ap.selfaln.tab";
  #FIXME: ids w/ \W bad here??    map{ s/\W+/_/g; } ($rd,$td); # cleanids .. need to fix all cds input ids: size, etc

  ## make cds/aa.qual unless exists.. evigene/scripts/prot/cdsqual.sh
  ##  or package sub makeAaQual ?  only good for evg headers
  $genesizes ||= "$cname.aa.qual";
  unless(-s $genesizes) {
    my $cmd="$EVIGENES/prot/cdsqual.sh $cdsseq";
    my($err)= runcmd($cmd);
  }
  ($nsizes,$aasizeh,$trsizeh,$cdspanh,$aaqualh) = ($genesizes) ? readSizes($genesizes) : (0,{},{});
  
  if($requiresize and $nsizes) { # this forces selfalign() to use only genesizes ids
    $idokset= $aasizeh; $nokids= $nsizes; 
  }
  
  if($cdsselfblast and -s $cdsselfblast) {
    warn "# selfalign blast input $cdsselfblast\n" if $DEBUG; 
  } else {  
    (my $blerr, $cdsselfblast)= blastself($cdsseq, $cname, $NCPU);
  }
  
  my($nseqs,$ndups)= selfalign($cdsselfblast);# , $oselfseq
  unless($nseqs){ die "#ERR: no input data, selfblast=$cdsselfblast\n"; }
  
  selfclass();
  dgenetable($genetab);

#=======================

sub blastself {
  my($cdsseq,$cname,$ncpu)= @_;
  
  my $isgz= ($cdsseq =~ /\.gz/)?1:0;
  if($isgz) {
    die "# FIXME: please gunzip input $cdsseq for this";
  }
  my $cdsdb= $cname || $cdsseq;
  # my $cdsname=`basename $cdsseq .cds | sed 's/.gz//; s/.cds_pub.fa//; s/\.fa.*$//'`;
  # (my $cdsdb= $cdsseq) =~ s/\.\w+$//;
  my $blout= "$cname.self$CDSBLAST_IDENT.blastn";
  my $selfblopt="-task megablast -ungapped -xdrop_ungap 4 -dust no -perc_identity $CDSBLAST_IDENT -evalue $CDSBLAST_EVALUE "; 
  my $outfmt='7'; # or? '7 std qseq'
  my($cmd, $err, $dropdb)= (0) x 9;
  
  if( -f $blout) { return(0,$blout); }
  elsif( -f "$blout.gz") { return(0,"$blout.gz");  }
  
  unless( -f "$cdsdb.nsq" ) {
  $cmd="makeblastdb -dbtype nucl -in $cdsseq -out $cdsdb";  
  $err= runcmd($cmd); $dropdb=1 unless($err);
  # my $pid= forkcmd($cmd);  
  }

  $cmd="blastn $selfblopt -db $cdsdb  -query $cdsseq -num_threads $ncpu -outfmt \"$outfmt\" -out $blout";
  $err= runcmd($cmd);
  if($dropdb){ unlink("$cdsdb.nsq"); unlink("$cdsdb.nhr"); unlink("$cdsdb.nin"); }
  return($err,$blout);
}


=item selfalign/chralign

  NOTE: run_refinealtpar.sh and this chralign extract refchrset.fa from evigene.gff of chr.blastn,
    may want update method for refchrset.fa, see 
    run_refinealtpar.sh gets refself.fa in same way as here, but separate SKIP_PARTS handling
    
  upd for both self/chr align? , SKIP_PARTS for putfa
  if [ $SKIP_PARTS = 1 ]; then
    perl -ne 'if(/^>/) { $ok=(m/alid=self/)?1:0; } print if $ok;' $refselfset > $refselfset.noparts
    mv $refselfset $refselfset.allparts
    mv $refselfset.noparts $refselfset
  fi

=cut


sub min{ ($_[1] < $_[0])? $_[1] : $_[0]; }
sub pctof{ my($nu,$de,$pt)=@_; my $dp= 10 ** ($pt||$PDIGITS);  return ($de>0)? int(0.5 + $dp*100*$nu/$de) / $dp : 0; }
sub round{ my($nu,$pt)=@_; my $dp= 10 ** ($pt||$PDIGITS);  return int(0.5 + $dp*$nu) / $dp; }

my(%aself,%aln,%alnsum,%dgid,%apclass);
#NOTUSEDHERE: my(%chrloc,%chrexon,%chrpar);

sub selfalign {
  my($inblast,$outfa)=@_;
  my($gok,$ndup,$nseq)=(0) x 9;
  #unused: ,%shdr,%seqh
  my($rpi,$ral,$rial,$rbs,$rnx)= (0) x 9; # rd.sums
  my $lrd=0;
  
  warn "# selfalign $inblast\n" if $DEBUG; # loggit
  if($inblast=~/\.gz/) {
    open(IN,"gunzip -c $inblast |") or die $inblast; 
  } else {
    open(IN,$inblast) or die $inblast;  
  }
  
  while(<IN>) {
    next if(/^\W/); chomp; my @v=split"\t"; 
    my($rd,$td,$pi,$al,$mism,$ndel,$rb,$re,$tb,$te,$evl,$bs,$rsq)=@v;
    # map{ s/\W+/_/g; } ($rd,$td); # cleanids ??
    next if($al < $MINA);
    next if($nokids and not ($idokset->{$rd} and $idokset->{$td})); #? td also?
    
    my($ad,$bd); my $or="f"; if($te < $tb){ ($tb,$te)=($te,$tb); $or="r"; } 
    if($rd eq $td){ 
      ($ad,$bd)=($rd,"self:$tb:$te:$or");  $nseq++;
      $aself{$rd}= [$al,$pi,$tb,$te,$or,$mism,$ndel,$bs]; 
      #u $shdr{$ad}=" alid=$bd;aln=$al,$pi%;";
    } else {
      $ad="$rd:$rb:$re"; $bd="$td:$tb:$te:$or"; $nseq++;
      push @{$aln{$rd}{$td}}, [$al,$pi,$rb,$re,$tb,$te,$or,$mism,$ndel,$bs]; #add mism,ndel?
      #u $shdr{$ad}=" alid=$bd;aln=$al,$pi%;"; 
    }
    
    #o $pi= int(0.5+$pi); $rsq||="NNN"; map{ s/:/_/g; } ($ad,$bd);
    #o if($shdr{$ad}){ $shdr{$ad} .= " alid=$bd;aln=$al,$pi%;"; $ndup++; } 
    #o else{ $seqh{$ad}=$rsq; $shdr{$ad}=" alid=$bd;aln=$al,$pi%;"; $nseq++ }
    $lrd= $rd;
    }
  close(IN);
  
  # alignseqout('selfalign',$SKIP_PARTS,$outfa,\%shdr, \%seqh);
  warn "# selfalign nseq=$nseq, ndups=$ndup\n" if $DEBUG; 
  return($nseq,$ndup);
}


sub selfclass {

  #top# my $ALTAL= 90;
  #top# my $ALTPI= 99.990; # or use mism <= 1/2 , ndel <= 0/1
  #top# my $GID=0; # allow user opt, ie to add to evg pubids
  
  my @rd= sort{ $aself{$b}->[0] <=> $aself{$a}->[0] or $a cmp $b} keys %aself; # long 1st
  my $nrd= @rd;
  warn "#selfclass ntr=$nrd\n" if $DEBUG;

  # calc alnsum from aln{rd,td} here..
  # UPD20fe: add class quals: althi1 == alt.redundant, part/frag classes
  # and evd/qual scores, as per asm dupfilter & altreclass,
  
  for my $rd (@rd) {
    my @td= sort keys %{$aln{$rd}}; 
    for my $td (@td) {
      my($rpi,$ral,$rial,$rbs,$rnx)=(0) x 9;
      my $apclass="";
      my @ax= @{$aln{$rd}{$td}};
      for my $ax (@ax) {
        my($al,$pi,$rb,$re,$tb,$te,$or,$mism,$ndel,$bs)= @$ax;
        my $ial= $al - $mism - $ndel;
        $rnx++; $ral += $al; $rbs+= $bs; $rpi+= $pi; $rial += $ial;
        # $ralhi += $al if($pi >= $ALTPI);
        if($al >= $ALTAL and $pi >= $ALTPI) {
          # if($apclass eq "par"){ $apclass="alt" if($ralhi >= 2*$ALTAL); } else {  $apclass="alt"; }
          $apclass="alt";
        } elsif($al >= $ALTAL) {
          $apclass="par" unless($apclass);
        }
      }
      $rpi=  ($rnx>0)? round( $rpi / $rnx) : $rpi; ## OR pctof($rpi/100,$rnx||1); 
      $alnsum{$rd}{$td}= [$rbs,$rial,$ral,$rpi,$rnx,$apclass||"nap"]; 
    }
    
    ## NOT USED HERE: %chrpar
    # if($chrpar{$rd}) { # add chrmap paralogs 
    #   my($rlen,$rpi,$rsb,$rse,$ror,$rmism,$rndel,$rbs)= @{$aself{$rd}};
    #   my @cpar=sort{ $chrpar{$rd}{$b}<=>$chrpar{$rd}{$a} } keys %{$chrpar{$rd}};
    #   for my $rdp (@cpar) {
    #     my $pval= $chrloc{$rdp}||"0,0,0,nochr"; ## "$rd:$rb:$re:$ro,$al,$pi%";
    #     my($ploc,$pal,$ppi)=split",",$pval;
    #     $aself{$rdp}= [$rlen,$ppi,$rsb,$rse,$ror,$rmism,$rndel,$rbs];# bogus vals for cpar
    #     $alnsum{$rd}{$rdp}= [$pal,$pal,$pal,$ppi,0,"parc"];
    #     #? apclass{gid??}{rdp}= "parc"; alnsum{rdp} ??
    #   }
    # }

  }

  for my $rd (@rd) {
    my($rlen,$rpi,$rsb,$rse,$ror,$rmism,$rndel,$rbs)= @{$aself{$rd}};
    # aself, alnsum include self-only == uni, aln has altpars
    my @apd= keys %{$alnsum{$rd}}; # Was keys %{$aln{$rd}};
    my $nap= @apd; my($gid,$ti)=(0,0);
    unless($gid= $dgid{$rd}) { $dgid{$rd}= $gid=  ++$GID; }
    unless(@apd) {
      #? $apclass{$rd}{self}= [ $gid, "uni", $ti, $nap,1,$rlen,$rpi,$rlen ];
      $ti=1; $apclass{$gid}{$rd}= [ "uni", $ti ];
    } else {
      #? swap rd for longest of rd,@apd here? so it is main of set
      #.. but rd has longest self-align, should be longest or eq long, or have gaps/repeats
      @apd= sort{ $alnsum{$rd}{$b}->[0] <=> $alnsum{$rd}{$a}->[0] or $a cmp $b } @apd; #? sort bitscore.sum or ialign
      $ti=1; $apclass{$gid}{$rd}= [ "main", $ti ]  unless($apclass{$gid}{$rd});
      for my $td (@apd) {
        ++$ti;
        my ($apclass,$amax,$pmax,$asum,$nx)=(0) x 9;
        my ($rbs,$rial,$ral,$rpi,$rnx,$rapclass)= @{$alnsum{$rd}{$td}};
    
        unless($dgid{$td}){ $dgid{$td}= $gid; }
        #? $apclass{$rd}{$td}= [ $gid, $apclass||"nap", $ti, $nap, $nx, $amax,$pmax,$asum ]; # $rbs,$ral,$rial,$rpi,$rnx
        $apclass{$gid}{$td}= [ $rapclass || "nap", $ti ] unless($apclass{$gid}{$td});
      } 
    }
  }
  warn "# selfclass dgenes=$GID for ntr=$nrd\n" if $DEBUG; 
}

=item dgenetable

  draft genes table, similar to evg pubids, row per transcript seq
  dg_trid oid(source_trid) dg_gid paralti apclass aaqual  laltids  lparids  lchrlocs

=cut

sub dgenetable {
  my($otab)= @_;
  
  warn "#dgenetable $otab\n" if $DEBUG;
  my $outh= *STDOUT;
  if($otab) {
    rename($otab,"$otab.old") if(-f $otab);
    open(OUT,">$otab"); $outh= *OUT;
  }
  
  my $DGHDR=0; # my($ndg)=scalar(%dgid);
  for my $rd (sort{ $dgid{$a}<=>$dgid{$b} } keys %dgid) {
    my $gid= $dgid{$rd};
    next unless($aself{$rd}); #?? err?
    my($rlen,$rpi,$rsb,$rse,$ror,$rmism,$rndel,$rbs)= @{$aself{$rd}};
    my @apd= sort{ $alnsum{$rd}{$b}->[0] <=> $alnsum{$rd}{$a}->[0] or $a cmp $b }keys %{$alnsum{$rd}};   
    my ($apclass,$ti)= @{ $apclass{$gid}{$rd} || ["nap",0] };

    my $gids= $IDPRE . sprintf "%06d",$gid;
    my $trid= $gids."t$ti";    
    my($aaqual,$trlen)=(0,0);
    if($nsizes and ref($aaqualh)) {
      $aaqual=$aaqualh->{$rd}||0; $trlen= $trsizeh->{$rd}||0;
    } 
    unless($aaqual and $trlen) {
      $aaqual="$rlen,1%,ac"; $trlen=$rlen; #this is self-align len
    }
    
    my ($laltids, $lparids, $lchrlocs)= ("noalt","nopar","noloc");
    my @altd= grep{ $apclass{$gid}{$_}->[0] =~ m/^(alt|main)/ } @apd; # eq "alt" 
    my @pard= grep{ $apclass{$gid}{$_}->[0] =~ m/^par/ } @apd; # par,parc
    $laltids= join",",@altd if(@altd);
    $lparids= join",",@pard if(@pard);
    ## NOT HERE: if(my $cloc=$chrloc{$rd}) { $lchrlocs= $cloc; }  # "$rd:$rb:$re:$or,$al,$pi%;"
    
    #UPD2020: add TRlen col : BAD trlen here, is aalen*3
    print $outh join("\t",qw(dRnaID OrigID  dGeneID AltI  APClass AAqual TRlen lAltIDs lParIDs lChrLocs))."\n" unless($DGHDR++);
    print $outh join("\t",$trid,$rd,$gids,$ti,$apclass,$aaqual,$trlen,$laltids,$lparids,$lchrlocs)."\n";
    
    # my $PRINTPAIRS= $ENV{dopairs}||0;
    # if($PRINTPAIRS) { # another table?
    # for my $td (@apd) {
    #   my($rbs,$rial,$ral,$rpi,$rnx,$rapclass)= @{$alnsum{$rd}{$td}};
    #   print $outh join("\t",$trid,"$rd,$td",$gids,0,$rapclass,$rbs,$rial,$ral,$rnx)."\n";
    # }
    # }
    
  }
  close(OUT) if($otab);
}

sub readSizes {
  my(@inf)= @_; 
  my $nt=0;
  my (%alen,%trlen,%cdspan,%aaqual); %alen=(); 
  my $hasgap=0; # ($AAGAP)?1:0;#  && $iscount .. my $iscount=1; 
  my $hasspan=0; # test for it.. $CDSSPAN; #  collect %trlen,%cdspan ? ** TEST sizes input for this?
  my $testspan=1; #(defined $CDSSPAN)?0:1; #? always test, ignore $CDSSPAN ?
  $asCDSsize= 0; #??? not set anywhere?
  
  foreach my $aaf (grep/\w/, map{ split(/[,\s]+/,$_) }@inf) {
    if(open(F,$aaf)) { my $n=0;
      while(<F>){ 
        next if(/^\W/); chomp; 
        my($id,$aw,@ac)=split"\t"; 
        #x if($hasgap) { if(@ac>1 or @ac==0){ $hasgap=0; } else { $aw -= $ac[0]; }} # OLD data/code, drop?
        ## dang new cds.qual has Code/Noncode col before cdspan col ..
        $alen{$id} = $aw; $n++; 
        $trlen{$id}= (@ac>=3 and $ac[2]=~/^\d/) ? $ac[2] : 3*$aw; # default, not ($asCDSsize)? 3*$aw : $aw;
        #^v -- want correct trlen, care less here of cdsspan 
        if($hasspan or $testspan) { 
          if(@ac>=3 and $ac[2]=~/^\d/) { 
	          my($gp,$aq,$tw,$csp,$cspx)=@ac; # csp == span OR Code/Noncode col ..
		        $csp||=""; $cspx||="";
		        if($tw) { $hasspan=1; $testspan=0; } else { $tw=3*$aw; }
            #old# $tw= $aw unless($tw); #? or 3*aw if aa > cds size?
            
	          my $isutrorf=(m/utrorf/)?1:0; # key in $oid may be missing
	          my $cspan= ($csp=~/^\d/) ? $csp : ($cspx=~/^\d/) ? $cspx : 0;
	          if($testspan) {
	            if($cspan) { $hasspan=1; $testspan=0; }
	            else { if(++$testspan>9) { $hasspan=0; $testspan=0; } }
	          }
	          if($asCDSsize) {
	            my $cw= $aw*3; $cw+=3 if($aq=~m/complete|partial5/);
	            if($cspan){ my($b,$e)= $cspan=~m/(\d+).(\d+)/; 
	              if($b and $e){ ($b,$e)=($e,$b) if($b>$e); $cw= 1 + $e - $b; } 
	              }
	            #NOT: $tw=$cw;
	            $tw= $cw if($cw > $tw);
	          }
            $aaqual{$id}= $aq; 
            $trlen{$id}= $tw;
            $cdspan{$id}= $cspan; # ($csp=~/^\d/) ? $csp : ($cspx=~/^\d/) ? $cspx : 0; # Code-col=3?  
            } 
          else { if(++$testspan >19 ) { $hasspan=0; $testspan=0; } } 
        } 
        
      } close(F); 
      
      $nt+=$n; warn  "# readSizes n=$n from $aaf\n" if $DEBUG;
    } else {
      warn "# cant read sizes from $aaf\n" ;# if $DEBUG
    }
  }
  return($nt,\%alen,\%trlen,\%cdspan,\%aaqual); # change to 1 hash w/ fields?
}


# from cdna_evigenesub.pm
sub runcmd
{
  my @cmd= @_;
  my $isdry=0; # ($dryrun or $cmd[0] =~ m/^echo\b/)?1:0;
  warn("#$EGLOG CMD=@cmd\n") if($DEBUG); # loggit( ($isdry) ? 1 : 0,"CMD=",@cmd); 
  my $err= ($isdry) ? 0 : system(@cmd);  
  if($err) { warn("# $EGLOG ERR=$err\n") } ## loggit(1,"ERR=$err ",$cmd[0]); } # ..
  return $err;
}

__END__

# sub alignseqout {
#   my($oname,$SKIP_PARTS,$outfa,$rshdr,$rseqh)=@_;
#   if($outfa) {
#   warn "#$oname faout=$outfa\n" if $DEBUG;
#   rename($outfa,"$outfa.old") if(-f $outfa);
#   open(OUT,">$outfa");
#   for my $ad (sort keys %$rshdr){ 
#     my $hdr=$rshdr->{$ad}; my $rsq=$rseqh->{$ad}; 
#     my $ok=($SKIP_PARTS and not $hdr=~m/alid=self/) ? 0 : 1;
#     print OUT ">$ad $hdr\n$rsq\n" if($ok); 
#     }
#   close(OUT);
#   }
# }

#================
