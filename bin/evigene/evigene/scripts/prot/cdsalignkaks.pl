#!/usr/bin/env perl
# cdsalignkaks.pl
# from  evigene run_cdsaln_kaks.sh, blastcds2axt.pl,
# expand this to add more cds align stats? basic refcds x testcds qual stats, as per aa homol

use FindBin;
use lib ("$FindBin::Bin"); # assume evigene/scripts/genes; layout: main, prot/, rnaseq/, ...

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname);
#? use cdna_evigenesub;  

# our $EVIGENES= $ENV{evigenes} || "$FindBin::Bin";  # fixme for other loc
our $EVIGENES= $FindBin::Bin; $EVIGENES =~ s,scripts/.*,scripts,;

our $EGAPP='cdsalnkaks';  
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

# blastcds2axt.pl
my $OVSLOP = $ENV{ovslop}||39; # what? not great, should use blast out order, keep all but major overlap
my $pOVSLOP= $ENV{povslop}||0.49; # use this for now
my $SKIPALT= $ENV{skipalt}||0; # skip query alts w/ ID pat > =~ /t1$/
# FIXME: upd skipalt skips if both query/targ gids are >t1, misses some valid matches (cross species)
#   change SKIPALT to (a) skip only query.t2+, a/o (b) keep best query.t1 x targ.tany
# already have this: didq# my $ONEQUERY=1; # skip 2ndary aligns, eg. db = all-alt cds db
my $ONEONLY=$ENV{oneonly}||0; #  self-cds can have  2+ paralogs, keep all aligns >= MINDUP*maxaln
my $MINDUP=$ENV{mindup}||0.50;
my(%didq,%didqr);

# FIXME findapp("blastn")
# FIXME findapp("kaks_calculator")
# my $kkbin="~/bio/apps/other/kaks_calculator/bin";
# my $kakscalc="$kkbin/KaKs_Calculator";

my($cdsseq,$refcdsseq,$logfile,$runname);
my($genesizes, $asCDSsize, $idset)= (0,1, 0);
my($nofreq,$nosigtab)=(0,0); #out opts: -[no]freq -[no]sigtab

my $optok= GetOptions(
  "cdsseq|testcds=s", \$cdsseq,   
  "refcds=s", \$refcdsseq,   
  "oneonly!", \$ONEONLY, 
  "skipalt!", \$SKIPALT, 
  "name|runname=s", \$runname,   
  #unu# "sizes|qual=s", \$genesizes,   
  #unu# "logfile:s", \$logfile,
  "pctident=s", \$CDSBLAST_IDENT,  
  "minalign=s", \$MINA,  
  #unu# "altpident=s", \$ALTPI,  
  "NCPU=i", \$NCPU, # "MAXMEM=i", \$MAXMEM,  
  #unu# "dryrun|n!", \$dryrun, 
  "debug!", \$DEBUG, 
  "nofreq!", \$nofreq, "nosigtab!", \$nosigtab, 
);

# two inputs: refcdsseq=dpt=subject.cds cdsseq=qpt=query.cds; exit -1; 
$cdsseq= shift(@ARGV) unless($cdsseq);
unless($refcdsseq){ $refcdsseq= shift(@ARGV); } # option self?
if($refcdsseq =~ /^(self|\-)/){ $refcdsseq=$cdsseq; }

die "# usage: cdsalnkaks testcds refcds \n opts: -[no]oneonly -[no]skipalt\n" 
  unless($optok and $cdsseq);

use constant kDIEonERR => 666;
my($errbl,$blastpath)= findapp("blastn",kDIEonERR);
my($errkk,$kkpath)= findapp("KaKs_Calculator",kDIEonERR);

my($nokids, $idokset)=(0,0);
my($nsizes,$aasizeh,$trsizeh,$cdspanh,$aaqualh)= (0,0,0,0);
# my $ALTAL= ($MINA>59) ? $MINA : 90; # ENV{altal}

sub MAIN_STUB {}

  # openloggit($logfile,$runname);  
  # loggit(1, "EvidentialGene cdsalignkaks.pl"); ##  "(-help for info), VERSION",VERSION);

  (my $cname= $cdsseq) =~ s/\.\w+$//;
  (my $rname= $refcdsseq) =~ s/\.\w+$//;
  my $pairname= $runname || "$cname-$rname.cdsaln"; # "refarath16ap.selfaln.tab";
  
  ## make cds/aa.qual unless exists.. evigene/scripts/prot/cdsqual.sh
  # $genesizes ||= "$cname.aa.qual";
  # unless(-s $genesizes) {
  #  my $cmd="$EVIGENES/prot/cdsqual.sh $cdsseq";
  #  my($err)= runcmd($cmd);
  #}
  #
  #? ($nsizes,$aasizeh,$trsizeh,$cdspanh,$aaqualh) = readSizes($genesizes) if($genesizes);
    
  my($blerr, $cdsblast)= blastcoding($cdsseq, $refcdsseq, $pairname, $NCPU);
  if($blerr){ die "#ERR: blastcoding($cdsseq, $refcdsseq)\n"; }
  
  my $cdsaxtc= $pairname.".axtc";
  my($erraxt)= blastcds2axt($cdsblast,$cdsaxtc); #UPD20apr: fixed for ncbi blastn fmt change
  
  my $cdskaks= $pairname.".kaks";
  my($errkk)= kkcalc($cdsaxtc,$cdskaks);
  if($errkk){ die "#ERR: kkcalc($cdsaxtc, $cdskaks)\n"; }

  my($errout)= kktables($cdskaks,$nosigtab,$nofreq); # out sigkaks.tab, freqks.tab .. other?

#=======================

sub kkcalc {
  my($cdsaxtc,$cdskaks)= @_;
  my $cmd="KaKs_Calculator -m MYN -i $cdsaxtc -o $cdskaks > $cdskaks.klog 2>&1";
  my($err)= runcmd($cmd);
  return($err);
}

=item kktables

cat $pt.kaks | perl -pe \
'($sd)=@v=split"\t"; if($sd=~/\slen=\d/){ ($td)= $sd=~/^([\w\.-]+)/; ($lw)= $sd=~m/len=(\d+)/; 
($qd)= $sd=~m,/([\w\.-]+),; ($aw)= $sd=~m/aln=(\d+)/; $tda="$td,$aw/$lw,$qd"; }
elsif($v[1]=~/^len=\d/){ ($td,$awl,$qd)=@v[0,1,2]; ($lw,$aw)= $awl=~m/len=(\d+).aln=(\d+)/;
$qd=~s,^[\d:]+/,,; $qd=~s/:\d.*$//; $tda="$td,$aw/$lw,$qd"; splice(@v,1,2); }
if($lw and $qd) { $v[0]=$tda; } else { $v[0]= "0$sd"; } $_=join"\t",@v; ' \
 | cut -f1,3,4,5,6 | grep -v 'NA' > $pt.sigkaks.tab

echo "#DONE $trset : `date`"

#... freq distr for kaks.tab
# set pt=tsaGDQRpinealb-tsaGBYRok.one

env kamax=9990.2 all=1 perl -ne 'BEGIN{ $KAMAX=$ENV{kamax}||999;
$TOP=$ENV{top}||0; $ALL=$ENV{all}||0; } ($d)=@v=split; if(@v==1) {
$nok++; $ok{$d}=1; } else { ($davdb,$ka,$ks,$kas,$pv)=@v; next if($ka
> $KAMAX or $ks=~/NaN/); ($da,$al,$db)=split",",$davdb; if(($ALL or
not $nok or ($ok{$da} and $ok{$db}))){  $rks=int($ks*10)/10;
$fks{$rks}++; $sks+=$ks; push @ks,$ks;  $nt++; } } END{ @fks=sort{$a
<=> $b} keys %fks; @sks=sort{$a <=> $b}@ks; $mdks=$sks[int($nt/2)];
$aks=sprintf "%.3f",$sks/$nt; print "nt=$nt; aveks=$aks,
medianks=$mdks \n",join("\t",qw(Ks pKs nKs))."\n"; for $ks (@fks) {
$c=$fks{$ks}; $p=int(0.5 + 1000*$c/$nt)/10; print "$ks\t$p\t$c\n"; } }' \
 $pt.sigkaks.tab 

=cut

sub kktables {
  my($cdskaks,$nosigtab,$nofreq)= @_;
  my($dosig,$dofreq)=(! $nosigtab, ! $nofreq);
  #my($dosig,$dofreq)=(1,1); $flags||="";
  #$dosig=0 if($flags=~/nosig/);
  #$dofreq=0 if($flags=~/nofreq/);
  
  my @ISK= map{ $_ - 1} (1,3,4,5,6); # == ($davdb,$ka,$ks,$kas,$pv)
  (my $tabsigkaks= $cdskaks) =~ s/\.\w+$/.sigkaks.tab/; 
  (my $tabfreq   = $cdskaks) =~ s/\.\w+$/.freqks.tab/; # from kaks or sigkaks ?? option?
  my $KAMAX=$ENV{kamax}||9990.2; #??
  my $fALL=$ENV{ksall}||0; 
  
  open(F,$cdskaks) or die "reading $cdskaks";
  if($dosig) { open(OSK,'>',$tabsigkaks) or die "writing $tabsigkaks"; }
  
  my($sd,$td,$qd,$lw,$aw,$tda,$awl,$qd)=(0) x 9;
  my($sks,$nt,@ks,%fks,)= (0,0); # dofreq
  
  while(<F>) {
    chomp; my @v=split"\t"; 
    ($sd)= @v;
    if($sd=~/\slen=\d/) { 
      ($td)= $sd=~/^([\w\.-]+)/; ($lw)= $sd=~m/len=(\d+)/; 
      ($qd)= $sd=~m,/([\w\.-]+),; ($aw)= $sd=~m/aln=(\d+)/; 
      $tda="$td,$aw/$lw,$qd";
    } elsif($v[1]=~/^len=\d/) { 
      ($td,$awl,$qd)=@v[0,1,2]; ($lw,$aw)= $awl=~m/len=(\d+).aln=(\d+)/;
      $qd=~s,^[\d:]+/,,; $qd=~s/:\d.*$//; 
      $tda="$td,$aw/$lw,$qd"; splice(@v,1,2); 
    }
    if($lw and $qd) { $v[0]=$tda; } else { $v[0]= "0$sd"; } 
    my($hasna)= grep /^(NA|NaN)$/i, @v[@ISK]; #? just NA or both? dang macos has 'nan' instead of NaN
    if($dosig) { print OSK join("\t",@v[@ISK])."\n" unless($hasna); }
    
    if($dofreq) {
      my($davdb,$ka,$ks,$kas,$pv)= @v[@ISK]; 
      unless($ka > $KAMAX or $ks=~/NaN/i) {
        my($da,$al,$db)=split",",$davdb; 
        if(1){  ## if($fALL or not $nok or ($ok{$da} and $ok{$db})) == id subset
          my $rks=int($ks*10)/10; 
          $fks{$rks}++; $sks+=$ks; push @ks,$ks;  $nt++; 
        }
      }
    }
    
  } close(F); 
  close(OSK) if($dosig);

  if($dofreq) {
    open(OFQ,'>',$tabfreq) or die "writing $tabfreq"; 
    my @fks=sort{$a <=> $b} keys %fks; 
    my @sks=sort{$a <=> $b} @ks; 
    my $mdks=$sks[int($nt/2)];
    my $aks=sprintf "%.3f",$sks/$nt; 
    print OFQ "# Ks frequency for nt=$nt, aveks=$aks, medianks=$mdks \n";
    print OFQ join("\t",qw(Ks pKs nKs))."\n"; 
    for my $ks (@fks) {
      my $c=$fks{$ks}; my $p=int(0.5 + 1000*$c/$nt)/10; 
      print OFQ "$ks\t$p\t$c\n"; 
      } 
    close(OFQ);
  }    
  
}

sub blastcoding {
  my($cdsseq,$refcds,$cname,$ncpu)= @_;
  
  my $isgz= ($cdsseq =~ /\.gz/)?1:0;
  if($isgz) { die "# FIXME: please gunzip input $cdsseq for this"; }
  my $refisgz= ($refcds =~ /\.gz/)?1:0;
  if($refisgz) { die "# FIXME: please gunzip input $refcds for this"; }

  my $cdsdb= $cname || $cdsseq;
  my $blout= "$cname.dcblast1n"; # cname == "$cname-$rname.cdsaln"
  
  if( -f $blout) { return(0,$blout); }
  elsif( -f "$blout.gz") { return(0,"$blout.gz");  }
  
  my $cdblopt="-task dc-megablast -template_type coding -template_length 18 -evalue 1e-9";
  my $outfmt=''; # default pairalign blastout ; NOT '7 std qseq'
  my($cmd, $err, $dropdb)= (0) x 9;
  
  unless( -f "$cdsdb.nsq" ) {
    $cmd="makeblastdb -dbtype nucl -in $cdsseq -out $cdsdb";  
    $err= runcmd($cmd); $dropdb=1 unless($err);
  }

  unless($err) {
    $cmd="blastn $cdblopt -db $cdsdb -query $refcds -out $blout"; 
    $err= runcmd($cmd);
  }
  
  if($dropdb){ unlink("$cdsdb.nsq"); unlink("$cdsdb.nhr"); unlink("$cdsdb.nin"); }
  return($err,$blout);
}



sub _min{ ($_[1] < $_[0])? $_[1] : $_[0]; }
sub _max{ ($_[1] > $_[0])? $_[1] : $_[0]; }
sub pctof{ my($nu,$de,$pt)=@_; my $dp= 10 ** ($pt||$PDIGITS);  return ($de>0)? int(0.5 + $dp*100*$nu/$de) / $dp : 0; }
sub round{ my($nu,$pt)=@_; my $dp= 10 ** ($pt||$PDIGITS);  return int(0.5 + $dp*$nu) / $dp; }
sub cleanid { my($d)=@_; $d=~s/^gi\|\d+\|//; $d=~s/\|$//; $d=~s/^(\w+)\|/$1:/; $d=~s/[^\w\.]/_/g; $d; } 

# from evigene/scripts/prot/blastcds2axt.pl
my( %al, %bal, %bb, %ee, );

sub blastcds2axt { 
  my($blastin,$axtout)= @_;
  
  if(-s $axtout) { return(0); } # or replace??
  my $ok;
  if($blastin =~ /.gz$/){ $ok= open(BLIN,"gunzip -c $blastin |"); }
  else { $ok= open(BLIN,$blastin); } 
  die "reading $blastin" unless($ok);
  open(AXOUT,'>',$axtout) or die "writing $axtout";
  my $outh= *AXOUT;
  my($qd,$qlen,$rd,$ina,$st,$bs,$err)=(0) x 9;
  while(<BLIN>) {
    
    if(/^Query=\s*(\S+)/){ $qd=cleanid($1); $qlen=0; } 
    elsif(/^Length=(\d+)/){ $qlen=$1 unless($qlen); } 
    elsif(/^Lambda/) { if(%al) { puta($rd); putb($outh,$qd,$rd,$qlen); } $rd=$ina=0; } 
     # fix mult align here? putb() only first rd? ignore splits?
     #UPD20Apr: dang ncbi blast format change, blastn fmt0 refid 
     # ncbi201[78] was: m/^> (\S+)/; now ncbi2019: m/^>(\S+)/
    elsif(/^>\s*(\S+)/) { if(%al) { puta($rd); putb($outh,$qd,$rd,$qlen); } $rd=cleanid($1); $ina=1; }
    elsif($ina) { 
      if(/^(Query|Sbjct)/) { $st=$1;
        my($sc,$bi,$al,$be)=split; 
        $bb{$st}=$bi if($bb{$st}<1); $ee{$st}=$be;
        $al{$st}.=$al; 
      } elsif(/^\s*Score =\s+(\d+)/) { 
        puta($rd) if(%al); $bs=1; $ina++; 
      } 
    } 
  } close(AXOUT); close(BLIN);
  
  return($err);
}

sub puta { 
  my($rd)= @_;
  return unless(%al);
  my($aq,$qb,$qe, $as,$rb,$re)= map{ ($al{$_},$bb{$_},$ee{$_}) } qw(Query Sbjct); 
  my $or=0; 
  if($qb>$qe) { $or++; ($qb,$qe)=($qe,$qb);}
  if($rb>$re) { $or++; ($rb,$re)=($re,$rb);} 
  $or=($or==1)?"-":"+";
  $bal{$qb}= [$qb,$qe,"$rd:$rb:$re:$or",$aq,$as]
    unless($bal{$qb});  
  %al=%bb=%ee=(); 
} 

sub isover {
  my($qb,$qe,$spans)=@_;
  for my $xbe (@$spans) { 
    my($xb,$xe)=@$xbe;
    if($qb<$xe and $qe>$xb) {
      my $qw=$qe - $qb; my $ov= _min($xe,$qe) - _max($xb,$qb);
      return 1 if($ov/$qw > $pOVSLOP);
      # return 1 unless($qb+$OVSLOP > $xe or $qe-$OVSLOP < $xb);
    }
  }
  return 0;
}

sub putb {
  my($outh, $qd,$rd,$qlen)=@_;
  my($lb,$le,$bb,$ml,$nskip)=(0) x 9;
  my($hd,$saq,$sas)=("") x 9;
  my @alnspan;
  $qlen||=0; 
  if($qd eq $rd) { %bal=(); return 0; } # SKIPSELF always, for self-cdsblast paralog tests
  # FIXME: upd skipalt skips if both query/targ gids are >t1, misses some valid matches (cross species)
  #   change SKIPALT to (a) skip only query.t2+, a/o (b) keep best query.t1 x targ.tany
  ## SKIPALT bug2: 1st blastin not best align (need all hsp, not best/partial hsp).
    # use: $didq{$qd} = $oaln
  my($qg,$qi)= ($qd=~m/(\w+)t(\d+)$/) ? ($1,$2) : ($qd,1); 
  my($rg,$ri)= ($rd=~m/(\w+)t(\d+)$/) ? ($1,$2) : ($rd,1); 
  if(1) { #($SKIPALT) .. SKIPSELF-ALT also
    if($qg eq $rg or ($SKIPALT and ($qi>1))) { %bal=(); return 0; } 
    #ob if($qg eq $rg or ($SKIPALT and ($qi>1 or $didqr{$qg.$rg}))) { %bal=(); return 0; } 
    #o if($qg eq $rg or ($SKIPALT and ($qi>1 or $ri>1))) { %bal=(); return 0; } 
  }
  ## below now
  # if(my $oaln=$didq{$qd}) { if($ONEONLY or $aln < $MINDUP*$oaln) { %bal=(); return 0; } }

  for my $i (sort{$a<=>$b}keys %bal) { 
    my $bx=$bal{$i}; 
    my($qb,$qe,$rloc,$aq,$as)=@$bx; my $d=0; 
    if(isover($qb,$qe,\@alnspan)) { $nskip++; next; }
    if($le>0 and $qb<=$le) { $d=1+$le-$qb; } 
    elsif($le>0 and $qb>$le+1) { $d=$qb - (1+$le); $d= 3 - ($d % 3); }
    if($d > 0) { $as=substr($as,$d); $aq=substr($aq,$d); $qb+=$d; }
    $sas.=$as; $saq.=$aq; $hd.="$qb:$qe/$rloc,"; 
    $bb=$qb if($bb==0);
    ($lb,$le)=($qb,$qe); push @alnspan, [$qb,$qe]; 
  } 

  # trim cds.query '-' align crap; after size checks?
  if($saq =~ m/\-/) {
    my(@nq,@ns); my @saq=split"",$saq; my @sas=split"",$sas;
    for(my $i=0; $i<=$#saq; $i++) { if($saq[$i] ne '-') { push @nq,$saq[$i]; push @ns,$sas[$i]; } }
    $saq=join"",@nq; $sas=join"",@ns;
    # while( my $i= rindex($saq,'-') >= 0 ) { #?? bad
    #  $saq=substr($saq,0,$i).substr($saq,$i+1); 
    #  $sas=substr($sas,0,$i).substr($sas,$i+1); 
    #}
  }
  
  $ml= ($bb-1) % 3;
  if( $ml > 0) { 
    my $d= 3-$ml; map{ $_= substr($_,$d); } ($saq,$sas); $bb+=$d;  # trim to codon start
  }
  my($lq,$ls)= map{length($_)} ($saq,$sas);
  if($lq ne $ls) {
    if($ls < $lq) { $lq=$ls; $saq=substr($saq,0,$lq); } 
    elsif($ls > $lq) { $sas=substr($sas,0,$lq); } 
  } 
  unless( ($ml= $lq % 3) == 0) { $lq-=$ml; map{ $_= substr($_,0,$lq); } ($saq,$sas); } 
  $hd=~s/,$//;
  my $aln= $saq =~ tr/A-Za-z/A-Za-z/;
  my $paln= ($qlen<1)?0:int(100*$aln/$qlen);
  my $oaln=$didq{$qd}||0;
  %bal=(); 
    # FIXME: self-cds can have  2+ paralogs, didq no good... need did-qd-rd?
  if($oaln) { return 0 if($ONEONLY or $aln < $MINDUP*$oaln); }
  print $outh "$qd\tlen=$qlen,aln=$aln,$paln%\t$hd\n";
  print $outh "$saq\n$sas\n\n";
  $didq{$qd}= $aln if($aln > $oaln);  #old: $didq{$qd}++;
  $didqr{$qg.$rg}= $aln;
  ## DANG missing saq or sas in large subset..
}  



sub readSizes {
  my(@inf)= @_; 
  my $nt=0;
  my (%alen,%trlen,%cdspan,%aaqual); %alen=(); 
  my $hasgap=0; # ($AAGAP)?1:0;#  && $iscount .. my $iscount=1; 
  my $hasspan=0; # test for it.. $CDSSPAN; #  collect %trlen,%cdspan ? ** TEST sizes input for this?
  my $testspan=1; #(defined $CDSSPAN)?0:1; #? always test, ignore $CDSSPAN ?
  
  foreach my $aaf (grep/\w/, map{ split(/[,\s]+/,$_) }@inf) {
    if(open(F,$aaf)) { my $n=0;
      while(<F>){ 
        next if(/^\W/); chomp; my($id,$aw,@ac)=split"\t"; 
        if($hasgap) { if(@ac>1 or @ac==0){ $hasgap=0; } else { $aw -= $ac[0]; }} 
        ## dang new cds.qual has Code/Noncode col before cdspan col ..
        $trlen{$id}= ($asCDSsize)? 3*$aw : $aw;

        if($hasspan or $testspan) { 
          if(@ac>=3 and $ac[2]=~/^\d/) { 
	          my($gp,$aq,$tw,$csp,$cspx)=@ac; # csp == span OR Code/Noncode col ..
		        $csp||=""; $cspx||="";
            $tw= $aw unless($tw); #? or 3*aw if aa > cds size?
	          my $isutrorf=(m/utrorf/)?1:0; # key in $oid may be missing
	          my $cspan= ($csp=~/^\d/) ? $csp : ($cspx=~/^\d/) ? $cspx : 0;
	          if($testspan) {
	            if($cspan) { $hasspan=1; $testspan=0; }
	            else { if(++$testspan>9) { $hasspan=0; $testspan=0; } }
	          }
	          if($asCDSsize){
	            my $cw= $aw*3; $cw+=3 if($aq=~m/complete|partial5/);
	            if($cspan){ my($b,$e)= $cspan=~m/(\d+).(\d+)/; 
	              if($b and $e){ ($b,$e)=($e,$b) if($b>$e); $cw= 1 + $e - $b; } 
	              }
	            $tw=$cw;
	          }
            $aaqual{$id}= $aq; $trlen{$id}=$tw;
            $cdspan{$id}= $cspan; # ($csp=~/^\d/) ? $csp : ($cspx=~/^\d/) ? $cspx : 0; # Code-col=3?  
            } 
          else { if(++$testspan>9) { $hasspan=0; $testspan=0; } } 
        } 
        $alen{$id}=$aw; $n++; 
      } close(F); 
      
      $nt+=$n; warn  "# readSizes n=$n from $aaf\n" if $DEBUG;
    } else {
      warn "# cant read sizes from $aaf\n" ;# if $DEBUG
    }
  }
  return($nt,\%alen,\%trlen,\%cdspan,\%aaqual); # change to 1 hash w/ fields?
}


# from cdna_evigenesub.pm
#  my($err,$kkpath)= findapp("KaKs_Calculator",kDIEonERR);
# pn=KaKs_Calculator;  pp=`which $pn`; pp=`basename $pp`; if [ "$pp" != "$pn" ]; then echo missing path to $pn; exit -1; fi
sub findapp {
  my($appname,$dieonerr)= @_;
  my($pn,$pp,$ppb,$err)= (0) x 9;
  $pn=$appname;  
  $pp=`which $pn`; chomp($pp); # or use sys status $? != 0
  $ppb=`basename $pp`;  chomp($ppb);
  if($ppb ne $pn){ $err="ERR: missing path to $pn"; die $err if($dieonerr); } # dieonerr == kDIEonERR   
  return($err,$pp);
}

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

