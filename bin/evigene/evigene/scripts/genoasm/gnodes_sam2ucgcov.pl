#!/usr/bin/perl
# gnodes_sam2ucgcov.pl
# special case of cov depth calc for uniq conserved genes, busco etc;
# C.UCG = Lread * Nread / Width_gene
# with filters to reject genes w/ duplicate reads, and reads with low map qual

use strict;
use Getopt::Long; # replace ENV w/ -opts

my $debug= $ENV{debug}||1;
my $MINE = $ENV{mine}|| 999; # NM:i:\d+ err, use 2..9, data dependent? dont want?
my $RDLEN=$ENV{rdlen}||150; # drop for measure sam col[9] = read
my $MINRDLEN=$ENV{minrdlen}||50; #? mainly for variable len w/ tiny read frags in some SRR
my $pMAXDUP= 0.05; #?? need tests
my $DOHETZ=0; #UPD19mar, test heterozygosity measure/base/gene, default?
my $HOMIN=0.80; my $HOMHET=0.98 - $HOMIN; # DOHETZ options? 

use constant { SAMFLAG_nomap => 0x04,SAMFLAG_rev => 0x10, SAMFLAG_2ndary => 0x100 }; # 0x100 == 256; # supplemental => 0x800;  

my ($crclassf,$outtab,$inbam,$dosumtab)=("") x 9; 
my @inbam; # opt to allow many -bam inputs, -bam one.bam two.bam three.bam ... or .sam ;  

my $optok= GetOptions( 
  'bam=s', \@inbam, # bam or sam or STDIN ? for samtools view call, add opt many @inbam
  'output|covtab=s',\$outtab, 
  'crclassf|idclassf=s',\$crclassf, # alt table chr => class for savereadids
  'maxerr|minident=s',\$MINE, #?? reuse for NM:i:(err) filter
  'maxdup|mindupident=s', \$pMAXDUP,  #? reuse MIN_DUPIDENT for dup reads/total reads filter
  'minrdlen|minreadlen=i',\$MINRDLEN, # add
  'rdlen|readlen=i',\$RDLEN, #  drop
  'heterozygosity!', \$DOHETZ, 
  'HOMIN=s', \$HOMIN, # for DOHETZ
  'summary|sumtab!', \$dosumtab, 
  'debug!', \$debug, 
  );

if(@ARGV) { push @inbam, grep(/\.(bam|sam)$/, @ARGV); }
$inbam= shift @inbam if(@inbam and not $inbam);

# my $optinok= (($inbam and -s $inbam) or (defined $MERGE and ($outtab or @ARGV)) );
unless($optok) {
  die "usage: gnodes_sam2ucgcov.pl -bam uniq_conserved_genes_readmap.bam -output uniq_conserved_genes.covtab
    or: gnodes_sam2ucgcov.pl < uniq_conserved_genes_readmap.sam > uniq_conserved_genes.covtab
    opts: -idclass genes.idclass : gene id table w/ UCG|busco class flag
    -maxdup=$pMAXDUP portion for uniq class, -maxerr=$MINE : max read map err, 
  ";
}

sub MAIN_stub {}

if($dosumtab){
  ucgcov_sumtab(); 
  exit;
}  

#? allow many -bam inputs, -bam one.bam two.bam three.bam ... or .sam ; does samtools view handle input list? no
my $inh = *STDIN; if($inbam) { open($inh,"samtools view -h $inbam |") or die "FAIL samtools view $inbam"; }
my $outh= *STDOUT; if($outtab) { rename($outtab,"$outtab.old") if(-s $outtab); open($outh, '>', $outtab) or die "writing $outtab"; }

my $topline= "#gnodes_sam2ucgcov options: maxdup=$pMAXDUP, maxerr=$MINE, idclass=$crclassf\n";
print $outh $topline; warn $topline if($debug and $outtab);

$HOMHET=0.98 - $HOMIN; 
# read map stats: n_mapok == nmapid? sam2cov uses n_mapok for all maps > nmapid
#  nr == n_mapok

my($ngene,  $nmapid, $stw, $scm, $srdlen, $srdlen_nomap, $lrid)= (0) x 9; # $nreadid, $nr
my($n_readid, $n_nomap, $n_mapok, $n_mapbad, $n_mapnotucg, $n_rdtooshort, # other: $n_dupbad, $n_partb,  
   $n_intron, $n_insert, $n_delete, $n_softclip, $n_mismatch)= (0) x 19; # globals?

my(%tw, %nrd, %nrdlen, %loc, %locdup, %cmd);
my $lrdlen= 0; # $RDLEN; #upd.18mar drop RDLEN usage, get from input.sam; add MINRDLEN

# assume inbam has all genes, want only UCG/busco subset, need idclassh for that
my($nidclass,$idclassh,$idclasslist)= ($crclassf) ? read_idclass($crclassf,1) : (0); # cds,te class by id, $idclass->{id} = class

foreach my $itbam ($inbam,@inbam) {
  if($itbam ne $inbam) { close($inh);  open($inh,"samtools view -h $itbam |") or next; }
  warn "#  samtools view $itbam\n" if($debug);
  
  while(<$inh>) {
    if(/^\@/) { 
      if(m/SN:/){ 
        my($td)=m/SN:(\S+)/; my($tw)=m/LN:(\d+)/; 
        #FIXME here? busco gene filter? maybe record all SN/LN, separate UCG later
        if($tw and isUCGene($td)){ $ngene++; $stw += $tw;  $tw{$td}=$tw; } # FOXME: ngene++ wrong for @inbam>1, use scalar(%tw)
      }
      next; } 
    elsif(/^\W/){ next; } 
  
    my @v=split; 
    my($srid,$fl,$td,$tb,$qs,$cig,$rseq)= @v[0,1,2,3,4,5,9];

    # UPDmar18: LNtot fix, variable rdlen: add sum_rdlen_notok and n_notok, need sum_rdlen_total= srdlen_ok+srdlen_no & n_ok+n_no for proper total LN
    # srdlen_no,n_no need to include nomap, and notUCGene for flag < 2ndary
    # add opt to skip all rdlen < $MINRDLEN =~ 50..100, ie tiny frag reads in some SRR sets
    # now: $srdlen_tot= $srdlen+$srdlen_nomap; $n_rdtot= $n_mapok + $n_nomap + $n_mapnotucg + $n_mapbad; #  n_rdtot should == n_readid **
    
    my $rdlen=length($rseq);  $rdlen=$lrdlen if($rdlen<3); # empty seq for some, use last rdlen?
    if($rdlen<$MINRDLEN and $fl < SAMFLAG_2ndary){ $n_rdtooshort++; next; }  #BUT check flags < SAMFLAG_2ndary  

    my $nextid= ($srid ne $lrid);
    $n_readid++ if($nextid);  $lrid=$srid; # assumes orig map/readid order
    
    #o: next if($fl & 0x04 or not isUCGene($td)); 
    #UPD: record all readid mapt to CDS genes, then filt by isUCGene
    if($fl & SAMFLAG_nomap) { $n_nomap++; $srdlen_nomap+=$rdlen; next; }
    elsif( not isUCGene($td) ) {
       if($nextid and $fl < SAMFLAG_2ndary) { $n_mapnotucg++; $srdlen_nomap+=$rdlen;  } #? $snotucg_rdlen += $rdlen; 
       next;
    }
   
    my($nmi)= m/NM:i:(\d+)/; # in @v[11..]
    my($zd,$zb,$zig)= (m/SA:Z:([\w\:\.\-]+),(\d+),.,(\w+)/) ? ($1,$2,$3) : (0,0,0);
    # id pat fix for prefix:id >> SA:Z:([\w\:\-]+), ?? or use td chars

use constant UPD_Cigarz => 1;  # finally debugged?

    if($fl >= SAMFLAG_2ndary or $nmi > $MINE) { 
      if($fl & SAMFLAG_2ndary) { 
if(UPD_Cigarz) {      
        my($alen,$lenc,$cend,$nlen)= addCigarz("D",$fl,$td,$tb,$cig, $zd,$zb,$zig) unless($nmi > $MINE); # skip err 2nd
        #not for 2ndary/suppl: if($nmi>0){ $alen -= $nmi; $n_mismatch += $nmi;  } 
} else {        
        addb("D",$fl,$td,$tb,$cig); # addCigar()
        addz($td,"D",$fl,$zd,$zb,$zig) if($zig); 
}        
      } elsif($fl < SAMFLAG_2ndary and $nmi>0) {
        $n_mismatch += $nmi;  $locdup{$td}{'Uerr'} += $nmi;  # count nmi err        
        $n_mapbad++; $srdlen_nomap+=$rdlen; #?? this way
      }

    } else { # good map ..
if(UPD_Cigarz) {      
      #UPD19mar: add $rseq
      my($alen,$lenc,$cend,$nlen)= addCigarz("U",$fl,$td,$tb,$cig, $zd,$zb,$zig, (($DOHETZ)?$rseq:""));
      if($nmi>0){ $n_mismatch += $nmi;  $locdup{$td}{'Uerr'} += $nmi; $alen -= $nmi; }  
} else {        
      addb("U",$fl,$td,$tb,$cig); # addCigar()
      addz($td,"U",$fl,$zd,$zb,$zig) if($zig);
}
      $nmapid++ if($nextid); # should be same as $n_mapok, 2nd above, be sure
      $n_mapok++; $srdlen += $rdlen; # if($nextid) ??
      $nrd{$td}++; $nrdlen{$td} += $rdlen;  #<< valuable stats
    }

    $lrdlen= $rdlen;
  } # itbam
} # for itbam (inbam,@inbam)

# END

## fixmed: add median rct from $nrdlen{$td}/$tw{td}
## FIXMEd: use only "uniq" set for summary stats .. need to calc after read loop
##   my $uok= ($sln<9)?"zero":($sld/$sln > $pMAXDUP)?"dupl":"uniq";

use constant { sZERO => 9, cZERO => 2 };

# fix @inbam $ngene, and $stw= sum(tw)
$ngene= scalar(keys %tw);
$stw=0; for my $td (sort keys %tw){ $stw += $tw{$td}; }

my $didhdr=0;
# totals including dup genes
my $rct=($stw<sZERO)?0:sprintf "%.1f", $srdlen / $stw; #was $n_mapok * $RDLEN / $stw ;
my $cmt=($stw<sZERO)?0:sprintf "%.1f", $scm / $stw;  
my $ardlen= ($n_mapok<sZERO)?0:sprintf "%.0f", $srdlen / $n_mapok; 

#UPD18mar new LN totals: readlen*n_readid
my $srdlen_tot= $srdlen+$srdlen_nomap; 
my $nrdlen_tot= $n_mapok + $n_nomap + $n_mapnotucg + $n_mapbad; #  nrdlen_tot should == n_readid **

my($urdlen,$utw,$ucm,$ucmw,$ssucmw)=(0) x 9; my(@crl,@cmw);
for my $td (sort keys %tw) { 
  my $w=$tw{$td} or next; 
  my $sln= $locdup{$td}{'N'}||0;  
  my $sld= $locdup{$td}{'D'}||0;
  my $isuniq= ($sln<sZERO)?0:($sld/$sln > $pMAXDUP)? 0 : 1;
  next unless($isuniq);
  
  my $cmd= $cmd{$td}; my $rclen=$nrdlen{$td}; # any zero/missing?  
  # dont int() these, lost precision.. for output?  
  my $crl= $rclen/$w; my $cmw= $cmd/$w; # was int()  
  $urdlen += $rclen; $ucm += $cmd; $utw += $w; 
  $ucmw += $cmw; $ssucmw += $cmw*$cmw; # use ucmw for ave, not ucm/utw
  push @crl,$crl; push @cmw,$cmw;  
}

@crl= sort{ $b <=> $a } @crl; @cmw= sort{ $b <=> $a } @cmw;
my $nc= @cmw; my $nch= int($nc/2);
my $mclr= $crl[$nch]; my $mcmw= $cmw[$nch];

#my $ucmwAve= ($nc<sZERO)? 0 : sprintf "%.1f", $ucmw/$nc;
my($ucmwMedn,$ucmwAve,$ucmwSD,$ucmwSE)=($cmw[$nch],0,0,0);
if($nc > sZERO) {
  $ucmwAve=$ucmw/$nc; 
  my $var=($ssucmw - $ucmwAve*$ucmwAve)/($nc-1); $ucmwSD=sqrt($var); 
  $ucmwSE=$ucmwSD/sqrt($nc);
}

# problem of extremes in  ucm/utw, use instead meanof(@cmw) as better stat
my $urct=($utw<sZERO)?0:sprintf "%.1f", $urdlen / $utw; #was $n_mapok * $RDLEN / $stw ;
my $ucmt= $ucmwAve; # old: ($utw<sZERO)?0:sprintf "%.1f", $ucm / $utw;  

my $nmapcdsid= $nmapid+$n_mapnotucg;
my $pmapucg= ($n_readid<sZERO)?0:sprintf"%.1f", 100*$nmapid/$n_readid; # pmr
my $pmapcds= ($n_readid<sZERO)?0:sprintf"%.1f", 100*$nmapcdsid/$n_readid;

# UPD add map err stats, maybe mod Gsize est
# $n_readid, $n_nomap, $n_mapok,  $n_mapnotucg,
#  base errs:  $n_insert, $n_delete,  $n_mismatch
#  relative to either $srdlen (nucg reads * rdlen) or $scm = mapt bases
my $nberr = $n_mismatch + $n_insert + $n_delete; # mid or sid; see also new locdup{td}{Uerr}
my $pberr = ($scm<sZERO)?0:sprintf "%.2f", 100*$nberr / $scm;  

# FIXME : med: mcmw or ave: ucmt ? see $ucmwAve update, best stat? 
# with clean data, med =~ ave, w/ unclean data median is more reliable est of true value
use constant C_MEDIAN => 1; # use median not ave, tested both, ave influenced by extremes

my $Cmap= (C_MEDIAN) ? $ucmwMedn : $ucmwAve; # was: $mcmw : $ucmt;
my $Genosize_mb= ($Cmap<cZERO) ? 0 :  sprintf"%.1f",  $ardlen * $n_readid / $Cmap / 1_000_000; 
my $CDS_mb= ($Cmap<cZERO) ? 0 :  sprintf"%.1f",  $ardlen * $nmapcdsid / $Cmap / 1_000_000; 

# GS = LN/C, L=ardlen N=n_readid,  C=C.Map ave or med?
# Nrmap=$n_mapok,  == Nrid.map=$nmapid .. should be same due to dup filter.. 
# ** Add cds-reads/tot-reads calc for %CDS/genome .. count all mapt cds-reads, record all input cds-spans? tw{}
# C.LN/W=$mclr,$urct >> C.LN/W=$mclr

my $GSEeqn= sprintf "L*N/C= %d * %d / %.1f",$ardlen,$n_readid,$Cmap;

#UPD18mar: other .. maybe print only when rdlen is variable or debug
# my $skipgsalt= ( $srdlen_tot == $ardlen*$n_readid) ? 1 : 0; # ardlen= $srdlen / $n_mapok
my $addGSEalt="";
my $dogsalt= (($n_mapok>=sZERO) and ($srdlen_tot == $srdlen * $n_readid / $n_mapok)) ? 0 : 1;  
if($dogsalt){
my $Genosize_alt= ($Cmap<cZERO) ? 0 :  sprintf"%.1f", $srdlen_tot / $Cmap / 1_000_000; 
my $ardlen_alt  = ($nrdlen_tot<1)? 0 : sprintf"%.1f", $srdlen_tot / $nrdlen_tot;
my $GSEeqn_alt  =  sprintf "LN/C= %d / %.1f for Lr: $ardlen_alt, Nr: %d = %d+ %d+ %d+ %d [ok,nomap,noucg,bad]", 
  $srdlen_tot, $Cmap,
  $nrdlen_tot, $n_mapok, $n_nomap, $n_mapnotucg, $n_mapbad; # srdlen_tot sum of all readlens counted, *should* == ave(rdlen) * n_readid
  $addGSEalt=  "# Gsize_alt  = $Genosize_alt Mb of $GSEeqn_alt\n";
}

my $CmapAveSE= sprintf "%.1f, %.1f +/-%.2f (mdn,ave,sem)",$ucmwMedn,$ucmwAve,$ucmwSE; # replace C.Map/W=$mcmw, $ucmt
$mclr=int($mclr);
my $topinfo= 
"# All  Gene Cov Depth n=$ngene, C.LN/W=$rct, C.Map/W=$cmt ave, for W.genes=$stw, LN.reads=$srdlen\n" .
"# Uniq Gene Cov Depth n=$nc, C.LN/W=$mclr,  C.Map/W=$CmapAveSE for W.genes=$utw, LN.reads=$urdlen\n" .
"# Genome_size= $Genosize_mb Mb of $GSEeqn, CDS_size= $CDS_mb Mb,\n" .
$addGSEalt .
"#   for Nr.ucgcds=$nmapid ($pmapucg%), Nr.allcds=$nmapcdsid ($pmapcds%), Nr.total=$n_readid, Lrdlen=$ardlen, MapErr=$nberr ($pberr%), Nshort=$n_rdtooshort\n"; 
#UPD18mar: add new LN.reads=srdlen_tot, others above
#my $srdlen_tot= $srdlen+$srdlen_nomap; 
#my $nrdlen_tot= $n_mapok + $n_nomap + $n_mapnotucg + $n_mapbad; #  nrdlen_tot should == n_readid **

print $outh $topinfo; warn $topinfo if($debug and $outtab);

for my $td (sort keys %tw) { 
  my $w=$tw{$td}; my $rc=$nrd{$td}; my $cmd=$cmd{$td}; 
  my $rclen=$nrdlen{$td};
  # dont int() for calcs, only output
  my $cmw= $cmd/$w; my $crl= $rclen/$w; # was int($rc*$RDLEN/$w); 

  my($sld,$slu,$sln,$slerr,$slun,$neq,$nnz, $shet, $shetho)=(0) x 19;
  # ($sld,$slu,$sln,$slerr)= map{ $locdup{$td}{$_} || 0 } qw( D U N Uerr);
  $slerr= $locdup{$td}{'Uerr'}||0;
  #>> change this per gene loc{td}[i] to locdup{td}{(D U N)} ?? need neq, nnz ?
  for(my $i=0; $i<$w; $i++) { 
    my($ld,$lu,$ln)= map{ $loc{$td}[$i]{$_}||0 } qw(D U N); 
    if($DOHETZ){
      my($ba,$bc,$bg,$bt)= map{ $loc{$td}[$i]{$_}||0 } qw(A C G T); 
      my $bn=$ba+$bc+$bg+$bt; # should == $ln
      my ($b1h,$b2h)=sort{ $b <=> $a} ($ba,$bc,$bg,$bt);
      # UPDhetz: add cds 1,2,3 position stats ($i % 3)
      if($bn > sZERO) {
       #above: my $HOMIN=0.90; my $HOMHET=0.98 - $HOMIN; # DOHETZ options? 
       my $ishet= ($b1h < $bn*$HOMIN and $b2h >= $bn*$HOMHET)?1:0; # does this hetz stat require b2h be high, > 10%?
       $shet++ if($ishet); $shetho++;
      }
    }
    $sld+=$ld; $slu+=$lu; $sln+=$ln; $nnz++ if($ln>0);
    if($ln>0 and $lu == $ln){ $slun+=$lu; $neq++; }
  }
  my $ceq=($neq<sZERO)?0:sprintf"%.1f", $slun/$neq;
  my $cnz=($nnz<sZERO)?0:sprintf"%.1f", $sln/$nnz;
  my $uok= ($sln<sZERO)?"zero":($sld/$sln > $pMAXDUP)?"dupl":"uniq";
  my $merr=($slu<sZERO)?0:sprintf"%.1f",  100*$slerr/$slu;
  my $phet=($shetho<sZERO)?0:sprintf"%.0f,%d/%d",  100*$shet/$shetho,$shet,$shetho;
  
  # fixme better output: keep: td,w,rc,crl,cmw, add:uok, maybe: cnz, ceq, maybnot: sld,u,n,un
  # Gene_ID Glen C.LN C.M Uniq C.nz C.eq Nread Rdlen S:dup,unq,nt,equ
  print $outh join("\t",qw(Gene_ID Glen Nread Rdlen C.LN C.M Uniq Merr C.nz S.Dup,Unq,Nt,Equ Hetz))."\n" 
    unless($didhdr++); 
  print $outh join("\t",$td,$w,$rc,$rclen, int($crl), int($cmw), $uok, $merr,
    "$cnz,$nnz", "$sld,$slu,$sln,$slun", $phet)."\n"; # off: ,"$ceq,$neq"
}  

my $botinfo= join", ", grep /\w/, map{ my $bn=`basename $_ .bam`; chomp($bn); ($bn)?$bn:""; } ($inbam,@inbam);
if($botinfo){
  $botinfo="# inbam: $botinfo\n";
  print $outh $botinfo;  
}

#--------------------------------------

sub isUCGene {
  my($id)= @_;
  if($nidclass) { # allow also for all-are-ucg option
    my $idclass= $idclassh->{$id} || "";
    return($idclass =~ m/UCG|busco/i)?1:0; # 
  } else {
    return 1; # no choice
  }  
}

sub addz { 
  my($od,$bt,$fl,$td,$tb,$cig)=@_; $bt="D" if($od ne $td); 
  addb($bt,$fl,$td,$tb,$cig); # addCigar() 
}

sub addCigarz {
  my($bt,$fl,$td,$tb,$cigar, $zd,$zb,$zig, $rseq)=@_; 
  my($alen,$lenc,$cend,$nlen,$softclip)= addCigar($bt,$fl,$td,$tb,$cigar, $rseq);
  if($zig) {
    $bt="D" if($zd ne $td); 
    my($zalen,$zlenc,$zcend,$znlen)= addCigar($bt,$fl, $zd,$zb,$zig, $rseq);
    $alen+= $zalen; $nlen+=$znlen; $lenc+=$zlenc; #?
    $cend=$zcend if($zd eq $td and $zcend > $cend);
  }
  return($alen,$lenc,$cend,$nlen,$softclip);
}

sub addb { ## old
  my($bt,$fl,$td,$tb,$cigar)=@_; 
  my($ci,$cm)=($tb,0);  

  while($cigar=~m/(\d+)([A-Z])/g){ 
    my($w,$t)=($1,$2); 
    unless($t eq "S" or $t eq "H"){ 
      if($t eq "M"){ $cm+=$w; for(my $i=0;$i<$w;$i++){ $loc{$td}[$ci+$i]{$bt}++; } }
      for(my $i=0;$i<$w;$i++){ $loc{$td}[$ci+$i]{"N"}++; } 
      $locdup{$td}{'D'} += $w if($t eq "M" and $bt ne "U"); #? or change to 'U' uniq count?
      $locdup{$td}{'U'} += $w if($t eq "M" and $bt eq "U"); #?  
      $locdup{$td}{'N'} += $w;
      }
    $ci+=$w; 
  } 
  if($bt eq "U"){ $cmd{$td}+=$cm; $scm+=$cm; } 
}

sub addCigar { #upd, also add err count/td : D+I+Nmi 
  my($tun,$fl,$td,$tb,$cigar,$rseq)=@_; 
  my($alen,$lenc,$cend,$nlen,$indel,$softclip)= (0) x 9;
  # cend == ci, alen == cm, bi == w # $cm+=$w; 
  
  # $cend= $tb; # my($ci,$cm)=($tb,0);  
  $tb--; $cend=$tb; # FIXME: tb/cend is 1-base, need 0-base for [cend+i] below

  #UPD19mar: add seq base count to $loc{$td}[$cend+$i]{A,C,G,T} to calc heterozygosity/base/gene
  my $dobases= ($DOHETZ and $rseq and $rseq=~/[ACGT]/);
  ## FIXME: $fl & SAMFLAG_rev => reverse @rseq ????
  my @rseq= ($dobases)?split("",$rseq):(); # this is read, 
  #DEBUG.off# if($dobases and $fl & SAMFLAG_rev){ @rseq= reverse @rseq; } #?? is this right
  # index $ir=$cend - $tb; $atbase= $rseq[$cend - $tb + $i]
  
  # ??replace w/ smocig() but need loc iter
  #  my($alen,$lenc,$cend,$softclip)= smokeCigar($cig,$tb);
  
  while($cigar =~ m/(\d+)([A-Z])/g) { 
    my($bi,$bt)=($1,$2); # my($w,$t)=($1,$2); 
 
    unless($bt eq 'H' or $bt eq 'S' or $bt eq 'I') { #? include both I/D indels for N? not I
      $nlen+=$bi; for(my $i=0;$i<$bi;$i++){ $loc{$td}[$cend+$i]{"N"}++; } 
    }
    if($bt eq 'H') { 
      $lenc += $bi; $bi=0; 
    } elsif($bt eq 'S') {
      $softclip += $bi; $lenc += $bi; $bi=0; 
    } elsif($bt eq 'M') {
      for(my $i=0;$i<$bi;$i++){ $loc{$td}[$cend+$i]{$tun}++; }
      if($dobases){ my $ri=$cend - $tb; for(my $i=0;$i<$bi;$i++){ 
        my $rb=$rseq[$ri + $i]; $loc{$td}[$cend+$i]{$rb}++ if($rb=~m/[ACGT]/); } 
      }
      $alen+= $bi; $lenc += $bi;  $cend += $bi;  
    } elsif($bt eq 'N') { 
      $n_intron++;  $cend += $bi; # cend not changed for HISP; 
    } elsif($bt eq 'I') {
      $n_insert++; $lenc += $bi; $bi=0; $indel++;
    } elsif($bt eq 'D') {  # P also?
      $n_delete++; $cend += $bi; $indel++;
    } elsif($bt eq 'P') {  # what is P now?
      $cend += $bi;
    } else {
      # unknown here, what? record?
    }
  }
  $locdup{$td}{'D'} += $alen if($tun ne "U"); # alen only for M
  $locdup{$td}{'U'} += $alen if($tun eq "U"); #? both
  $locdup{$td}{'N'} += $nlen; # cend - cbeg rather than alen
  #x $locdup{$td}{'Err'} += $indel if($tun eq "U"); # only record for Uniqs ?, ie dont want D/paralog read err to skew Uniq err
  if($tun eq "U"){ 
    $cmd{$td}+=$alen; $scm+=$alen; 
    $locdup{$td}{'Uerr'} += $indel; # or 'ME' for maperr, or Uerr to match U align
  } 
  
  #/////old////  
  # while($cigar =~ m/(\d+)([A-Z])/g) {  my($w,$t)=($1,$2); 
  #   unless($t eq "S" or $t eq "H"){ 
  #     if($t eq "M"){ $cm+=$w; for(my $i=0;$i<$w;$i++){ $loc{$td}[$ci+$i]{$bt}++; } }
  #     for(my $i=0;$i<$w;$i++){ $loc{$td}[$ci+$i]{"N"}++; } 
  #     $locdup{$td}{'D'} += $w if($t eq "M" and $bt ne "U"); #? or change to 'U' uniq count?
  #     # $locdup{$td}{'U'} += $w if($t eq "M" and $bt eq "U"); #?  
  #     $locdup{$td}{'N'} += $w;
  #     }
  #   $ci+=$w; 
  # } 
  
  return($alen,$lenc,$cend,$nlen,$softclip);#?
}

=item ucgcov_sumtab

  .. brief row summary from ucgcovtabs + metadata
  cols now: Species_SRAset G_fcyto G_asm G_ucg  C.M.mdn,ave,sem  N_ucg W_ucg Nr_ucg,% Nr_cds,% Nr_tot Len_rd SRA_Info 
  
  ( for gf in *20gnodes; do { head -7  $gf/*.metad $gf/*.ucgcovtab $gf/*ucg.covtab; } done ) | ./gnodes_ucgcov_sumtab.pl
  see daphnia/gnodes_ucgcov_sumtab.pl

=cut
  
sub ucgcov_sumtab { 

  ## drop these special cases
  #? dromel=SRR6366285 is that bact.contam set?
  #? arath=SRR10178322 is that insect.contam set?
  my $skiprdset="arath=ERR3624574,arath=SRR10178322,pig=SRR6263260,fig=DRR187753,dromel=SRR10512945,dromel=SRR6366285";
  my $skiprdpat=join"|", map{ my($k,$v)=split"=",$_; $v; }split",",$skiprdset;
  
  my @hd=qw(Species_SRAset G_fcyto G_asm G_ucg  C.M.mdn,ave,sem  N_ucg W_ucg Nr_ucg,% Nr_cds,% Nr_tot Len_rd SRA_Info );
  print "Table UC1. UCG coverage depth summary\n\n";
  print join("\t",@hd)."\n";

  my ($ingn,$gnset,$gnfile,$gse,$gsf,$nskip)=("0") x 9; 
  my ($idcla,$srrid,$lenrd,$nrucg,$nrcds,$nrtot,$sppid,$flocy,$chrasm,$cmu,$nucg,$wucg,$lnucg)= ("0") x 19;
  my ($cmd,$cma,$cse,$maperr,$lucg,$lcds,$prucg,$prcds,$pmaperr)=("0") x 19;      
  my (%meta);
  while(<>) { 
    my $inl=$_;
    
      # ==> apis20gnodes/apismel20asm.metad <==
      # ==> apis20gnodes/apismel14evg3cds_SRR9108936a.ucgcovtab <==
      # ==> aweed20gnodes/arath20asm.metad <==
      # ==> aweed20gnodes/arath18tcds_ERR3624574a.ucgcovtab <==
    if(m,^==\> (\w+)/(\S+),) { $gnset=$1; $gnfile=$2; }
  
    ## spp.metad
    # asmid=chrpig11c
    # asmtotal=2501 Mb
    # flowcyto=2924-3139 Mb
    # species=Sus_scrofa
    if(/^asmid=/){ %meta=(); }
    if(/^(species|flowcyto|asmtotal|asmname|asmid)\s*=\s*(.+)$/) { my($mt,$mv)=($1,$2); $mv=~s/#.*$//; $meta{$mt}=$mv; }
    
    ## gnodes_ucgcov summaries
    if(/([SDE]RR\d+[ab12_]*)/) { ($srrid)= $1; $srrid=~s/_1/a/; $srrid=~s/_2/b/;  $srrid=~s/_$//; }
    if(/^#gnodes_sam2ucgcov/) {
      ($idcla)=m/idclass=(\S+)/?$1:"noidc"; $ingn=1;
      ($lenrd,$nrtot,$sppid,$flocy,$chrasm,$cmu,$nucg,$wucg)= ("0") x 19;
      ($sppid)= $idcla=~m/^(\w+)/;
      
    } elsif(/# Uniq Gene Cov Depth/) {
      ($nucg, $wucg, $lnucg)= map{ my($v)= $inl =~ m/\b$_=\s*(\S+)/?$1:0; $v; } qw(n W.genes LN.reads );
      ($cmu)= $inl =~ m,C.Map/W=(.+).mdn.ave, ? $1 : "na";
      ($cmd,$cma,$cse)= split(/[,\s]+/,$cmu);
      
    } elsif(/# Genome_size=/) {
      ($gse)= $inl =~ m,Genome_size=\s*(\S+), ? $1 : "na";
      ($gsf)= $inl =~ m=LN/C.\s*([^,]+)= ? $1 : "na"; #? print or not, printd parts
      
    } elsif(/Nr.total=/) { # last,print row
      $lucg= (m/(Nrd.ucgmap)/)?$1:"Nr.ucgcds";
      $lcds= (m/(Nrd.allcds)/)?$1:"Nr.allcds";
      ($nrucg,$nrcds,$nrtot,$lenrd,$maperr)= 
        map{  my($v)= $inl =~ m/\b$_=\s*(\S+)/?$1:0; $v;  } 
        ($lucg,$lcds,"Nr.total","Lrdlen","MapErr");
      ($prucg,$prcds,$pmaperr)= 
        map{ my($v)= $inl =~ m/\b$_=\s*\S+\s+.([\d\.]+%)/?$1:0; $v;  } 
        ($lucg,$lcds,"MapErr");
        
      my $gse_new= ($cmd<1)?0: sprintf "%.1f", $lenrd * ($nrtot/1_000_000) / $cmd;
      
      my $flowcyto= $meta{flowcyto}||0;
      my $asmtotal= $meta{asmtotal}||0;
      my $asmname= $meta{asmname}||$meta{asmid}||0;
      my $species= $meta{species}||0;
      #xx if($species and $sppid) { $sa= sppabbr($species); $sppid=$sa.$sppid unless($sppid=~m/$sa/i); } 
      my $xinfo= join",", grep( /\w/, $asmname,$srrid,$species);
      map{ s/\s*Mb\b//i } ($flowcyto, $asmtotal, $gse_new);
      map{ s/\s+$//; s/\s*,$//; } ($cmu, $nucg, $wucg, $nrtot, $lenrd);
      
      my $sppidc=$sppid; map{ s/(_t1cds|_cds1t|_cds|t1cds|1cds|ucgcds|t1m_cds)$//; } ($sppidc);
      if($srrid =~ m/$skiprdpat/){ $nskip++; }
      else {
      printf "%-15s\t",$sppidc;
      print join("\t", $flowcyto, $asmtotal, $gse_new, $cmu, $nucg, $wucg, "$nrucg,$prucg", "$nrcds,$prcds", $nrtot, $lenrd, $xinfo)."\n";
      }
      $ingn=0; $srrid="0"; 
      # %meta=(); #<< FIXME need meta w/ spp tag, read all, use $meta{asmid}{@vals} ??
    }
  }
}

sub sppabbr{
  my($spp)=@_;
  my($gn,$sn)= split/[_\s]/,$spp,2;
  my $sa=$spp;
  if($sn){ $sa=substr($gn,0,3) . substr($sn,0,3); }
  else { $sa=substr($spp,0,6); } return ($sa);
}


sub read_idclass {
  my($crclassf,$allclasses,$crlenh)=@_; 
  my($nid,$ncl,$ok,$inh)=(0,0,0);
  my %crclass=(); my @crclass=();
  $allclasses||=0;

  use constant kMAXIDCLASS => 9; # idclass limit, using 3-4 now
  my $CRTPAT=''; # no default, see sam2covtab
  
  if($crclassf and ($ok,$inh)= openRead($crclassf) and $ok) { 
    while(<$inh>){ next if(/^\W/); 
      my($cr,$crclass,@clx)=split;  # may have more columns .. keep all ie CDS,BUSCO 
      if($allclasses and @clx){ $crclass=join(" ",$crclass,@clx); }
      $crclass{$cr}= $crclass || 0;  
    } close($inh); 
    $nid= scalar(keys %crclass);
  }
  #.. insert other way from sam hdr ids and CRTPAT?
  # if($nid==0 and $CRTPAT and ref($crlenh)){
  #   for my $cr (sort keys %$crlenh) {
  #     my($crt)= ($cr =~ m/($CRTPAT)/)? $1 : 'UNK'; 
  #     $crclass{$cr}= $crt;
  #   }
  #   $nid= scalar(keys %crclass);
  # }
  
  unless($nid>0) { @crclass=('UNK'); }
  else {
    my %crcl=(); # make list of values = classes
    map{ $crcl{$crclass{$_}}++; } keys %crclass;  # or values %crclass > not uniq
    @crclass= sort keys %crcl; # sort by count?
  }
  
  $ncl= @crclass; 
  if($ncl > kMAXIDCLASS) { # problem, cancel..
    warn "#ERR too many idclasses n=$ncl from nid=$nid of $crclassf  \n";# or sam ids x CRTPAT='$CRTPAT'
    $nid=0; %crclass=(); @crclass=();
  }
  warn "# read nid=$nid, nclass=$ncl from $crclassf\n" if($debug);
  return($nid,\%crclass,\@crclass);
}  

sub openRead {
  my($fn)=@_; my($ok,$inh)=(0); 
  return(0,undef) unless($fn and -f $fn);
  if($fn =~ /\.gz/) { $ok= open($inh,"gunzip -c $fn |"); } else { $ok= open($inh,$fn); }
  return($ok,$inh);
}

__END__

=item test outputs
# See gnodes_sam2ucgcov.sum.txt for much of these test outputs
#----------------------------------

=item try daphplx

perl gnodes_sam2ucgcov.pl -debug  -idclass dplx20gnodep/dplx20cdste.idclass -bam dplx20gnodep/daphplx17evgt1m_cds_SRR3090572_1_bwa.bam  -out dplx20gnodep/daphplx17evgt1m_cds_SRR3090572a.try3ucgcovtab

head -15 dplx20gnodep/daphplx17evgt1m_cds_SRR3090572a.try4ucgcovtab
#gnodes_sam2ucgcov options: maxdup=0.05, maxerr=999, idclass=dplx20gnodep/dplx20cdste.idclass
# All  Gene Cov Depth n=853, C.LN/W=54.8, C.Map/W=36.9 ave, for W.genes=1230513, LN.reads=67390049
# Uniq Gene Cov Depth n=725, C.LN/W=53,53.9  C.Map/W=36,36.6 (mdn,ave) for W.genes=1001661, LN.reads=54007474
# Genome_size= 210.2 Mb of LN/C= 228 * 33737295 / 36.6, CDS_size= 104.2 Mb,
#   for Nrd.ucgmap=296136 (0.9%), Nrd.allcds=16733408 (49.6%), Nr.total=33737295, Lrdlen=228

Gene_ID                 Glen    Nread   Rdlen   C.LN    C.M     Uniq    C.nz            S.Dup,Unq,Nt,Equ
Daplx7pEVm000149t1      8235    1628    369268  44      34      dupl    79.3,8223       357780,283700,652239,3914
Daplx7pEVm000164t1      8016    1864    422365  52      46      uniq    49.2,7996       15868,375309,393045,329730
Daplx7pEVm000244t1      6726    1847    418950  62      45      dupl    64.4,6723       118626,307083,433235,99886
Daplx7pEVm000387t1      6006    1176    266604  44      37      uniq    37.5,5996       487,223588,224600,216063
Daplx7pEVm000546t1      5430    1239    282165  51      42      uniq    42.9,5420       0,229860,232451,215102
Daplx7pEVm000553t1      5115    1010    229694  44      33      uniq    35.1,5112       3197,171337,179669,145550
Daplx7pEVm000573t1      5307    1256    284393  53      39      uniq    40.2,5300       2000,209830,212833,190814
Daplx7pEVm000640t1      4884    1099    250476  51      38      uniq    39.9,4882       4422,188309,194596,168642
Daplx7pEVm000650t1      5229    1156    261203  49      38      dupl    53.2,5228       72987,199709,277930,108123

head  dplx20gnodep/daphplx17evgt1m_cds_SRR3090572a.try3ucgcovtab
#gnodes_sam2ucgcov options: maxerr=999, maxdup=0.05, idclass=dplx20gnodep/dplx20cdste.idclass
# All  Gene Cov Depth n=853, C.LN/W=54.8, C.Map/W=36.9 ave, for Width.genes=1230513, LN.reads=67390049
# Uniq Gene Cov Depth n=725, C.LN/W=53,53.9  C.Map/W=36,36.6 (mdn,ave) for Width.genes=1001661, LN.reads=54007474
# Genome_size=210.2 Mb of LN/C= 228 * 33737295 / 36.6, for Nrid.map=296136 (0.9%), Nrid.in=33737295, Rdlen=228

Gene_ID                 Glen    Nread   Rdlen   C.LN    C.M     Uniq    C.nz           S.Dup,Unq,Nt,Equ
Daplx7pEVm000149t1      8235    1628    369268  44      34      dupl    79.3,8223       357780,283700,652239,3914
Daplx7pEVm000164t1      8016    1864    422365  52      46      uniq    49.2,7996       15868,375309,393045,329730
Daplx7pEVm000244t1      6726    1847    418950  62      45      dupl    64.4,6723       118626,307083,433235,99886
Daplx7pEVm000387t1      6006    1176    266604  44      37      uniq    37.5,5996       487,223588,224600,216063
Daplx7pEVm000546t1      5430    1239    282165  51      42      uniq    42.9,5420       0,229860,232451,215102

daphplx17evgt1m_cds_SRR3090572a.try3ucgcovtab  stats for nt=726
Item    Median  Mean    SEM     Nitem   StDev   Sum
  Glen   1152   1380.77 60.60   726     1632.83 1002438
 Nread   278    327.23  14.17   726     381.93  237571
  C.LN    53    55.39   2.12    726     57.01   40216
   C.M    36    35.89   1.36    726     36.57   26056
  C.nz    36    36.57   1.38    726     37.19   26546

----

perl gnodes_sam2ucgcov.pl -debug  -idclass dplx20gnodep/dplx20cdste.idclass -bam dplx20gnodep/daphplx17evgt1m_cds_SRR3090572_1_bwa.bam  -out dplx20gnodep/daphplx17evgt1m_cds_SRR3090572a.ucgcovtab
#gnodes_sam2ucgcov options: maxerr=999, maxdup=0.05, idclass=dplx20gnodep/dplx20cdste.idclass
# read nid=34298, nclass=4 from dplx20gnodep/dplx20cdste.idclass
#Uniq Gene Cov Depth n=725, C.LN/W=53,53.9 med/ave, C.Map/W=36,36.6 med/ave for Wgenespan=1001661, Nrdlen=54007474

head dplx20gnodep/daphplx17evgt1m_cds_SRR3090572a.ucgcovtab
#gnodes_sam2ucgcov options: maxerr=999, maxdup=0.05, idclass=dplx20gnodep/dplx20cdste.idclass
#Uniq Gene Cov Depth n=725, C.LN/W=53,53.9 med/ave, C.Map/W=36,36.6 med/ave for Wgenespan=1001661, Nrdlen=54007474
#All Gene Cov Depth n=853, C.LN/W=54.8 ave, C.Map/W=36.9 ave for Wgenespan=1230513, Nrdlen=67390049, Nread=296136, Rdlen=228
Gene_ID                 Glen    Nread   Rdlen   C.LN    C.M     Uniq    C.nz    S.Dup,Unq,Nt,Equ
Daplx7pEVm000149t1      8235    1628    369268  44      34      dupl    79.3,8223       357780,283700,652239,3914
Daplx7pEVm000164t1      8016    1864    422365  52      46      uniq    49.2,7996       15868,375309,393045,329730
Daplx7pEVm000244t1      6726    1847    418950  62      45      dupl    64.4,6723       118626,307083,433235,99886
Daplx7pEVm000387t1      6006    1176    266604  44      37      uniq    37.5,5996       487,223588,224600,216063
Daplx7pEVm000546t1      5430    1239    282165  51      42      uniq    42.9,5420       0,229860,232451,215102
Daplx7pEVm000553t1      5115    1010    229694  44      33      uniq    35.1,5112       3197,171337,179669,145550

grep -v dupl dplx20gnodep/daphplx17evgt1m_cds_SRR3090572a.ucgcovtab | grep -v '^#' | sed 's/^Gene_/#Chr/;' | env icols=1,2,3,4,5,7 median=1 sum=1 ./cds_meanvar.pl
 stats for nt=726
  Item    Median  Mean    SEM     Nitem   StDev   Sum
    Glen   1152   1380.77 60.60   726     1632.83 1002438
   Nread   278    327.23  14.17   726     381.93  237571
    C.LN    53    55.39   2.12    726     57.01   40216
     C.M    36    35.89   1.36    726     36.57   26056
    C.nz    36    36.57   1.38    726     37.19   26546
  
  ==> dplx20gnodep/daphplx16mlupr1a_SRR3090572_sum.txt <==
  Source=Daphplx16ml, KUlow=36, KUhigh=36, FlowcytSize=215-264 Mb Formula_LN/C=185.9 Mb (daphplx_gasm16ml)
  ==> dplx20gnodep/daphplx18paupr1a_SRR3090572_sum.txt <==
  Source=Daphplx19ml, KUlow=31, KUhigh=36, FlowcytSize=215-264 Mb Formula_LN/C=190.3-221.0 Mb (daphpulex_pa42v2)
  ==> dplx20gnodep/daphplx20ma1bupr1a_SRR3090572_sum.txt <==
  Source=Daphplx20ma1b, KUlow=32, KUhigh=35, FlowcytSize=215-264 Mb Formula_LN/C=211.6-231.5 Mb (daphplx20maca1b_findeg)

  ==> dplx20gnodep/daphplx16mlupr1a_SRR3090572_sum.txt <==
    Genome Size Est= 201.5 Mb (Nread), 185.9 Mb (Maprd), for readset SRR3090572,
     for Size=LN/C, Cov=36,36, N_reads=33737296, N_maprd=31129065,92.3%, L_readlen=215
    Uniq Conserved Gene Cover: median=36, sem=1.140, n=1076
  
  ==> dplx20gnodep/daphplx18paupr1a_SRR3090572_sum.txt <==
    Genome Size Est= 203.4-236.2 Mb (Nread), 190.3-221.0 Mb (Maprd), for readset SRR3090572,
     for Size=LN/C, Cov=31,36, N_reads=33737296, N_maprd=31578447,93.6%, L_readlen=217
    Uniq Conserved Gene Cover: median=36, sem=1.188, n=931
  
  ==> dplx20gnodep/daphplx20ma1bupr1a_SRR3090572_sum.txt <==
    Genome Size Est= 213.0-233.0 Mb (Nread), 211.6-231.5 Mb (Maprd), for readset SRR3090572,
     for Size=LN/C, Cov=32,35, N_reads=33737296, N_maprd=33515655,99.3%, L_readlen=221
    Uniq Conserved Gene Cover: median=35, sem=1.062, n=1137


=item dropse20 
  .. rename Width > W, or Wid, short for equation
  ** CDS_size= 32.2 Mb close to dropse20cur3a_SRR11813283_sum.txt 33 Obs mb, 34-35 Est.mb
  
head -15 dropse20gnodes/dropse20t1cds_SRR11813283a.try4ucgcovtab
#gnodes_sam2ucgcov options: maxdup=0.05, maxerr=999, idclass=dropse20gnodes/dropse20cdste.idclass
# All  Gene Cov Depth n=1048, C.LN/W=116.8, C.Map/W=94.4 ave, for Width.genes=1475412, LN.reads=172296300
# Uniq Gene Cov Depth n=1020, C.LN/W=116,116.1  C.Map/W=93,94.2 (mdn,ave) for Width.genes=1413654, LN.reads=164150100
# Genome_size= 168.1 Mb of LN/C= 150 * 105560783 / 94.2, CDS_size= 32.2 Mb,
#   for Nrd.ucgmap=1148642 (1.1%), Nrd.allcds=20241345 (19.2%), Nr.total=105560783, Lrdlen=150 

Gene_ID                 Glen    Nread   Rdlen   C.LN    C.M     Uniq    C.nz    S.Dup,Unq,Nt,Equ
dropseg117183176t1      618     668     100200  162     99      uniq    97.3,617        0,60027,60029,59966
dropseg117183718t1      2094    1729    259350  123     111     uniq    111.5,2093      0,233370,233374,232978
dropseg117184007t1      1356    1204    180600  133     117     uniq    117.7,1355      0,159457,159458,159331
dropseg117184773t1      702     613     91950   130     100     uniq    99.9,701        0,70037,70037,70037
dropseg13036350t1       7896    5928    889200  112     107     uniq    108.0,7885      0,851810,851810,851810
dropseg26531963t1       3261    2302    345300  105     92      dupl    107.9,3260      55106,296290,351609,260942
dropseg26531986t1       2361    2291    343650  145     107     uniq    107.4,2359      0,253326,253328,253309
dropseg26532704t1       681     689     103350  151     104     uniq    101.3,680       0,68861,68861,68861
dropseg26533243t1       1260    1511    226650  179     141     dupl    204.2,1258      78667,178261,256928,62134


perl gnodes_sam2ucgcov.pl -debug  -idclass dropse20gnodes/dropse20cdste.idclass -bam dropse20gnodes/dropse20t1cds_SRR11813283_1_bwa.bam -out dropse20gnodes/dropse20t1cds_SRR11813283a.try3ucgcovtab
#gnodes_sam2ucgcov options: maxerr=999, maxdup=0.05, idclass=dropse20gnodes/dropse20cdste.idclass
# read nid=15225, nclass=4 from dropse20gnodes/dropse20cdste.idclass

#Uniq Gene Cov Depth n=1020, C.LN/W=116,116.1 med/ave, C.Map/W=93,94.2 med/ave for Wgenespan=1413654, Nrdlen=164150100
#All Gene Cov Depth n=1048, C.LN/W=116.8 ave, C.Map/W=94.4 ave, Genome_size(LN/C)=168 Mb 
  GSE Accurate, add nums? (150*105.560783/94.2)
# for Wgenespan=1475412, Nrdlen=172296300, Nrmap=1148642, Nrid.map=1148642, Nrid.in=105560783 (1.1%), Rdlen=150


=item try from full cds.bam .. diff from test busco.sam

perl gnodes_sam2ucgcov.pl -debug -maxerr=5 -maxdup=0.05 -idclass dropse20gnodes/dropse20cdste.idclass -bam dropse20gnodes/dropse20t1cds_SRR11813283_1_bwa.bam -out dropse20gnodes/dropse20t1cds_SRR11813283a.ucgcovtab
# read nid=15225, nclass=4 from dropse20gnodes/dropse20cdste.idclass

head dropse20gnodes/dropse20t1cds_SRR11813283a.ucgcovtab
#gnodes_sam2ucgcov options: maxerr=5, maxdup=0.05, idclass=dropse20gnodes/dropse20cdste.idclass
#Uniq Gene Cov Depth n=1019, C.LN/W=115,115.2 med/ave, C.Map/W=92,93.4 med/ave for Wgenespan=1413321, Nrdlen=162768750
#All Gene Cov Depth n=1048, C.LN/W=115.7 ave, C.Map/W=93.5 ave for Wgenespan=1475412, Nrdlen=170756850, Nread=1138379, Rdlen=150
Gene_ID         Glen    Nread   Rdlen   C.LN    C.M     Uniq    C.nz            S.Dup,Unq,Nt,Equ
g117183176t1    618     666     99900   161     99      uniq    97.0,617        0,59860,59862,59799
g117183718t1    2094    1721    258150  123     110     uniq    110.9,2093      0,232179,232183,231791
g117184007t1    1356    1196    179400  132     116     uniq    116.8,1355      0,158312,158313,158186
g117184773t1    702     611     91650   130     99      uniq    99.5,701        0,69737,69737,69737
g13036350t1     7896    5885    882750  111     107     uniq    107.2,7885      0,845614,845614,845614
g26531963t1     3261    2285    342750  105     91      dupl    107.1,3260      55106,293870,349189,258672

perl gnodes_sam2ucgcov.pl -debug -maxerr=99 -maxdup=0.05 -idclass dropse20gnodes/dropse20cdste.idclass -bam dropse20gnodes/dropse20t1cds_SRR11813283_1_bwa.bam -out dropse20gnodes/dropse20t1cds_SRR11813283a.err99ucgcovtab
>> high maxerr raises C by 1 here, closer to chrmap C.UCG

# read nid=15225, nclass=4 from dropse20gnodes/dropse20cdste.idclass
head dropse20gnodes/dropse20t1cds_SRR11813283a.err99ucgcovtab
#gnodes_sam2ucgcov options: maxerr=99, maxdup=0.05, idclass=dropse20gnodes/dropse20cdste.idclass
#Uniq Gene Cov Depth n=1020, C.LN/W=116,116.1 med/ave, C.Map/W=93,94.2 med/ave for Wgenespan=1413654, Nrdlen=164150100
#All Gene Cov Depth n=1048, C.LN/W=116.8 ave, C.Map/W=94.4 ave for Wgenespan=1475412, Nrdlen=172296300, Nread=1148642, Rdlen=150
Gene_ID                Glen    Nread   Rdlen   C.LN    C.M     Uniq    C.nz    S.Dup,Unq,Nt,Equ
dropse20uc:g117183176t1 618     668     100200  162     99      uniq    97.3,617        0,60027,60029,59966
dropse20uc:g117183718t1 2094    1729    259350  123     111     uniq    111.5,2093      0,233370,233374,232978
dropse20uc:g117184007t1 1356    1204    180600  133     117     uniq    117.7,1355      0,159457,159458,159331
dropse20uc:g117184773t1 702     613     91950   130     100     uniq    99.9,701        0,70037,70037,70037
dropse20uc:g13036350t1  7896    5928    889200  112     107     uniq    108.0,7885      0,851810,851810,851810
dropse20uc:g26531963t1  3261    2302    345300  105     92      dupl    107.9,3260      55106,296290,351609,260942

grep -v dupl dropse20gnodes/dropse20t1cds_SRR11813283a.err99ucgcovtab | grep -v '^#' | sed 's/^Gene_/#Chr/;' | env icols=1,2,3,4,5,7 median=1 sum=1 ./cds_meanvar.pl
 stats for nt=1020
Item    Median  Mean    SEM     Nitem   StDev   Sum
  Glen   1137   1385.94 53.83   1020    1719.34 1413654
 Nread   880    1072.88 41.47   1020    1324.48 1094334
  C.LN   116    123.43  5.10    1020    162.74  125896
   C.M    93    95.00   3.12    1020    99.54   96895
  C.nz    93    95.17   3.13    1020    99.91   97069
  
head -1  dropse20gnodes/*sum.txt
  ==> dropse20gnodes/dropse20chrs_cov7b_sum.txt <==
  Source=Dropse20uc, KUlow=94, KUhigh=95, FlowcytSize=161-180 Mb Formula_LN/C=174 Mb (dropse20chrs)
  ==> dropse20gnodes/dropse20cur3a_SRR11813283_sum.txt <==
  Source=Dropse20uc, KUlow=94, KUhigh=95, FlowcytSize=161-180 Mb Formula_LN/C=174 Mb (dropse20chrs)
  ==> dropse20gnodes/dropse20t1cds_SRR11813283_sum.txt <==
  Source=dropse20t1cds, KUlow=87, KUhigh=88, FlowcytSize=0 Formula_LN/C=0 (dropse20t1cds)

=item try3 cucumber

perl gnodes_sam2ucgcov.pl -debug  -idclass cucum20gnodes/cucum19nc_cds1t.idclass -bam cucum20gnodes/cucum19nc_cds1t_SRR11300859_1_bwa.bam -out cucum20gnodes/cucum19nc_cds1t_SRR11300859a.try3ucgcovtab 

#gnodes_sam2ucgcov options: maxerr=999, maxdup=0.05, idclass=cucum20gnodes/cucum19nc_cds1t.idclass

 head cucum20gnodes/cucum19nc_cds1t_SRR11300859a.try3ucgcovtab
#gnodes_sam2ucgcov options: maxerr=999, maxdup=0.05, idclass=cucum20gnodes/cucum19nc_cds1t.idclass
# All  Gene Cov Depth n=1364, C.LN/W=39.3, C.Map/W=27.0 ave, for Width.genes=2147307, LN.reads=84483000
# Uniq Gene Cov Depth n=1305, C.LN/W=38,38.3  C.Map/W=25,26.1 (mdn,ave) for Width.genes=2046972, LN.reads=78304950
# Genome_size=422.9 Mb of LN/C= 150 * 73590832 / 26.1, for Nrid.map=563220 (0.8%), Nrid.in=73590832, Rdlen=150

Gene_ID Glen    Nread   Rdlen   C.LN    C.M     Uniq    C.nz    S.Dup,Unq,Nt,Equ
cucum19nc_g101202736t1  591     102     15300   25      20      uniq    21.0,589        0,12349,12349,12349
cucum19nc_g101202737t1  1329    375     56250   42      31      uniq    32.3,1308       0,42299,42299,42299
cucum19nc_g101202749t1  1899    368     55200   29      27      uniq    27.9,1896       0,52860,52860,52860
cucum19nc_g101202761t1  507     102     15300   30      20      uniq    19.6,503        0,9840,9840,9840
cucum19nc_g101202769t1  2085    695     104250  50      33      uniq    33.6,2070       0,69188,69526,66991

cucum19nc_cds1t_SRR11300859a.try3ucgcovtab stats for nt=1304
Item    Median  Mean    SEM     Nitem   StDev   Sum
  Glen   1374   1568.91 51.55   1304    1861.44 2045865
 Nread   350    400.04  14.31   1304    516.63  521648
  C.LN    38    39.09   1.19    1304    43.01   50970
   C.M    25    25.34   0.75    1304    27.00   33047
  C.nz    25    25.82   0.76    1304    27.49   33665

cuc chrasm
==> cucum19cgir9bs_SRR11300859_sum.txt <==
Source=Cucumber19CGI, KUlow=25, KUhigh=25, FlowcytSize=372-509 Mb Formula_LN/C=427.1 Mb (cucum19cgi_chr)
  Total span=226.6 Mb, covspan=225.3, gaps=0 Mb for assembly cucum19cgi_chr
  Genome Size Est= 441.5 Mb (Nread), 427.1 Mb (Maprd), for readset SRR11300859,
   for Size=LN/C, Cov=25,25, N_reads=73590833, N_maprd=71188878,96.7%, L_readlen=150
  Uniq Conserved Gene Cover: median=25, sem=0.663, n=1427
  
==> cucum20pccr9bs_SRR11300859_sum.txt <==
Source=Cucumber20PCC, KUlow=25, KUhigh=25, FlowcytSize=372-509 Mb Formula_LN/C=408.1 Mb (cucum20pcc_asm)
  Total span=341.8 Mb, covspan=333.6, gaps=0 Mb for assembly cucum20pcc_asm
  Genome Size Est= 441.5 Mb (Nread), 408.1 Mb (Maprd), for readset SRR11300859,
   for Size=LN/C, Cov=25,25, N_reads=73590833, N_maprd=68011549,92.4%, L_readlen=150
  Uniq Conserved Gene Cover: median=25, sem=0.655, n=1449

=item try4 aweed

head -15 aweed20gnodes/arath18tcdsbusco_ERR4586299a.try4ucgcovtab
#gnodes_sam2ucgcov options: maxdup=0.05, maxerr=999, idclass=aweed20gnodes/arath18tair1cds.idclass
# All  Gene Cov Depth n=1397, C.LN/W=73.8, C.Map/W=51.7 ave, for W.genes=2168564, LN.reads=160021650
# Uniq Gene Cov Depth n=1318, C.LN/W=75,73.9  C.Map/W=51,51.7 (mdn,ave) for W.genes=2026028, LN.reads=149770650
# Genome_size= 156.5 Mb of LN/C= 150 * 53925997 / 51.7, CDS_size= 62.1 Mb,
#   for Nrd.ucgmap=1066811 (2.0%), Nrd.allcds=21410853 (39.7%), Nr.total=53925997, Lrdlen=150

Gene_ID         Glen    Nread   Rdlen   C.LN    C.M     Uniq    C.nz    S.Dup,Unq,Nt,Equ
AT1G01140t1     1344    954     143100  106     53      dupl    54.5,1343       4065,69034,73199,60968
AT1G01220t1     2601    1244    186600  71      58      uniq    59.0,2599       0,153055,153394,145998
AT1G01290t1     813     357     53550   65      54      uniq    54.4,811        0,44081,44081,44081
AT1G01350t1     1032    397     59550   57      47      dupl    95.2,1025       49394,47913,97578,2157
AT1G01760t1     1260    639     95850   76      48      uniq    49.2,1254       0,61555,61661,60679
AT1G01770t1     1899    983     147450  77      47      uniq    47.7,1893       0,90374,90375,90316
AT1G01880t1     1797    751     112650  62      46      uniq    47.1,1795       0,84129,84580,82932
AT1G01930t1     1743    827     124050  71      52      uniq    51.2,1741       0,89159,89159,89159
AT1G01970t1     1230    455     68250   55      48      uniq    48.4,1228       0,59490,59491,59438


perl gnodes_sam2ucgcov.pl -debug -idclass aweed20gnodes/arath18tair1cds.idclass -bam aweed20gnodes/arath18tair1cds_ERR4586299_1_bwa.bam -out aweed20gnodes/arath18tcdsbusco_ERR4586299a.try3ucgcovtab
head aweed20gnodes/arath18tcdsbusco_ERR4586299a.try3ucgcovtab

#gnodes_sam2ucgcov options: maxerr=999, maxdup=0.05, idclass=aweed20gnodes/arath18tair1cds.idclass
# All  Gene Cov Depth n=1397, C.LN/W=73.8, C.Map/W=51.7 ave, for Width.genes=2168564, LN.reads=160021650
# Uniq Gene Cov Depth n=1318, C.LN/W=75,73.9  C.Map/W=51,51.7 (mdn,ave) for Width.genes=2026028, LN.reads=149770650
# Genome_size=156.5 Mb of LN/C= 150 * 53925997 / 51.7, for Nrid.map=1066811 (2.0%), Nrid.in=53925997, Rdlen=150

Gene_ID         Glen    Nread   Rdlen   C.LN    C.M     Uniq    C.nz    S.Dup,Unq,Nt,Equ
AT1G01140t1     1344    954     143100  106     53      dupl    54.5,1343       4065,69034,73199,60968
AT1G01220t1     2601    1244    186600  71      58      uniq    59.0,2599       0,153055,153394,145998
AT1G01290t1     813     357     53550   65      54      uniq    54.4,811        0,44081,44081,44081
AT1G01350t1     1032    397     59550   57      47      dupl    95.2,1025       49394,47913,97578,2157
AT1G01760t1     1260    639     95850   76      48      uniq    49.2,1254       0,61555,61661,60679

arath18tcdsbusco_ERR4586299a.try3ucgcovtab  stats for nt=1318
Item    Median  Mean    SEM     Nitem   StDev   Sum
  Glen   1377   1537.20 50.13   1318    1820.10 2026028
 Nread   682    757.57  24.00   1318    871.32  998471
  C.LN    75    76.57   2.16    1318    78.34   100913
   C.M    51    51.02   1.41    1318    51.18   67244
  C.nz    51    51.14   1.41    1318    51.31   67400

  
perl gnodes_sam2ucgcov.pl -debug -maxerr=5 -maxdup=0.05 -idclass aweed20gnodes/arath18tair1cds.idclass -bam aweed20gnodes/arath18tair1cds_ERR4586299_1_bwa.bam -out aweed20gnodes/arath18tcdsbusco_ERR4586299a.ucgcovtab
# read nid=27444, nclass=2 from aweed20gnodes/arath18tair1cds.idclass

head  aweed20gnodes/arath18tcdsbusco_ERR4586299a.ucgcovtab
#gnodes_sam2ucgcov options: maxerr=5, maxdup=0.05, idclass=aweed20gnodes/arath18tair1cds.idclass
#Uniq Gene Cov Depth n=1318, C.LN/W=74,72.7 med/ave, C.Map/W=50,50.6 med/ave for Wgenespan=2026028, Nrdlen=147240900
#All Gene Cov Depth n=1397, C.LN/W=72.3 ave, C.Map/W=50.4 ave for Wgenespan=2168564, Nrdlen=156769650, Nread=1045131, Rdlen=150
Gene_ID         Glen    Nread   Rdlen   C.LN    C.M     Uniq    C.nz    S.Dup,Unq,Nt,Equ
AT1G01140t1     1344    952     142800  106     52      dupl    54.4,1343       4065,68858,73018,60860
AT1G01220t1     2601    1134    170100  65      53      uniq    53.5,2599       0,139060,139076,137606
AT1G01290t1     813     357     53550   65      54      uniq    54.4,811        0,44081,44081,44081
AT1G01350t1     1032    395     59250   57      46      dupl    94.9,1025       49394,47626,97291,2157
AT1G01760t1     1260    599     89850   71      45      uniq    45.8,1254       0,57402,57403,57340
AT1G01770t1     1899    979     146850  77      47      uniq    47.5,1893       0,89898,89899,89840

arath18tcdsbusco_ERR4586299a.ucgtab stats for nt=1318
Item    Median  Mean    SEM     Nitem   StDev   Sum
  Glen   1377   1537.20 50.13   1318    1820.10 2026028
 Nread   664    744.77  23.62   1318    857.62  981606
  C.LN    74    75.37   2.13    1318    77.24   99340
   C.M    50    50.02   1.38    1318    50.19   65929
  C.nz    50    50.19   1.39    1318    50.36   66155

head -2  aweed20gnodes/a*_sum.txt
  ==> aweed20gnodes/arath18cds1a_sum.txt <==
  Source=arath18cds1a, KUlow=90, KUhigh=92, FlowcytSize=0 Formula_LN/C=0 (arath18cds1a)
  ==> aweed20gnodes/arath18tair1a_sum.txt <==
  Source=Arath18TAIR, KUlow=52, KUhigh=52, FlowcytSize=157-166 Mb Formula_LN/C=0 (arath18tair_chr)
  ==> aweed20gnodes/arath20maxr1a_sum.txt <==
  Source=Arath20Max, KUlow=52, KUhigh=52, FlowcytSize=157-166 Mb Formula_LN/C=0 (arath20max_chr)

=item try4

  .. maxerr=1,3,5 make largish diff in C.UCG (83 .. 94 .. 99) ; any way to resolve read map errs?
  .. should apply same maxerr filter to chrasm read map
  
gunzip -c dropse20gnodes/dropse20cdsbusco_SRR11813283_2_bwa.sam.gz | $gw/daphnia/gnodes_sam2ucgcov.pl -debug -maxerr 1  | less

#gnodes_sam2ucgcov options: maxerr=1, maxdup=0.05, idclass=
#Uniq Gene Cov Depth n=1017, C.LN/W=83,83.0 med/ave, C.Map/W=65,65.4 med/ave for Wgenespan=1403262, Nrdlen=116439900
#All Gene Cov Depth n=1048, C.LN/W=83.2 ave, C.Map/W=65.3 ave for Wgenespan=1475412, Nrdlen=122758950, Nread=818393, Rdlen=150
Gene_ID      Glen    Nread   Rdlen   C.LN    C.M     Uniq    C.nz    S.Dup,Unq,Nt,Equ
g117183176t1 618     515     77250   125     74      uniq    73.4,616        0,45208,45208,45208
g117183718t1 2094    1254    188100  89      79      uniq    79.3,2093       0,166039,166040,165945
g117184007t1 1356    856     128400  94      82      uniq    82.9,1350       0,111859,111860,111752
g117184773t1 702     432     64800   92      71      uniq    71.5,700        0,50036,50036,50036
g13036350t1  7896    4298    644700  81      77      uniq    77.3,7885       0,609472,609472,609472
g26531963t1  3261    1547    232050  71      60      dupl    75.2,3259       49556,195295,244998,173187
g26531986t1  2361    1797    269550  114     83      uniq    83.0,2360       0,195843,195845,195690
g26532704t1  681     512     76800   112     73      uniq    69.6,675        0,46952,46952,46952
g26533243t1  1260    1053    157950  125     95      dupl    155.3,1257      75430,119821,195251,39919
g26533276t1  780     340     51000   65      55      uniq    55.3,778        0,43052,43052,43052
g26533533t1  1209    725     108750  89      75      uniq    75.5,1208       0,91217,91217,91217

=item try3

  .. maxerr=5 vs maxerr=3 makes largish diff in C.UCG (94 .. 99) ; try maxerr=1 ?

gunzip -c dropse20gnodes/dropse20cdsbusco_SRR11813283_2_bwa.sam.gz | $gw/daphnia/gnodes_sam2ucgcov.pl -debug  | less

#gnodes_sam2ucgcov options: maxerr=5, maxdup=0.05, idclass=
#Uniq Gene Cov Depth n=1019, C.LN/W=99,99.7 med/ave, C.Map/W=79,79.2 med/ave for Wgenespan=1409835, Nrdlen=140592900
#All Gene Cov Depth n=1048, C.LN/W=100.2 ave, C.Map/W=79.3 ave for Wgenespan=1475412, Nrdlen=147842100, Nread=985614, Rdlen=150
Gene_ID      Glen    Nread   Rdlen   C.LN    C.M     Uniq    C.nz            S.Dup,Unq,Nt,Equ
g117183176t1 618     577     86550   140     84      uniq    83.6,616        0,51486,51498,51315
g117183718t1 2094    1479    221850  105     94      uniq    94.1,2093       0,196926,196927,196826
g117184007t1 1356    1009    151350  111     97      uniq    97.8,1350       0,132045,132046,131925
g117184773t1 702     505     75750   107     84      uniq    84.8,700        0,59331,59331,59331
g13036350t1  7896    5137    770550  97      92      uniq    92.4,7885       0,728676,728677,728553
g26531963t1  3261    1975    296250  90      77      dupl    91.7,3259       49556,249118,298829,218659
g26531986t1  2361    2074    311100  131     96      uniq    96.3,2360       0,227269,227275,227076
g26532704t1  681     612     91800   134     89      uniq    85.8,675        0,57883,57883,57883
g26533243t1  1260    1337    200550  159     123     dupl    183.4,1257      75430,155127,230560,56230
g26533276t1  780     432     64800   83      71      uniq    71.3,778        0,55469,55470,55397
g26533533t1  1209    877     131550  108     92      uniq    92.1,1208       0,111316,111316,111316
g26533670t1  654     528     79200   121     96      uniq    95.4,652        0,62213,62216,62036

=item try2

  .. maybe drop C.Map as spurious, unclear result
  
gunzip -c dropse20gnodes/dropse20cdsbusco_SRR11813283_2_bwa.sam.gz | \
  $gw/daphnia/gnodes_sam2ucgcov.pl -debug -minident 3 -maxdup 0.05  

All=1048, C.LN/W=94.7 ave, C.Map/W=74.7 ave for Wgenespan=1475412, Nread=931339, Nrdlen=139700850
Uniq=1018, C.LN/W=94,94.3 med/ave, C.Map/W=74,74.7 med/ave for Wgenespan=1406430, Nrdlen=132685650
Gene_ID      Glen    Nread   Rdlen   C.LN    C.M     Uniq    C.nz            Sdup,unq,tot,u=t
g117183176t1 618     554     83100   134     80      uniq    79.2,616        0,48763,48763,48763
g117183718t1 2094    1411    211650  101     89      uniq    89.7,2093       0,187733,187734,187635
g117184007t1 1356    955     143250  105     92      uniq    92.6,1350       0,125041,125042,124925
g117184773t1 702     484     72600   103     80      uniq    81.0,700        0,56705,56705,56705
g13036350t1  7896    4865    729750  92      87      uniq    87.6,7885       0,690681,690681,690681
g26531963t1  3261    1847    277050  84      71      dupl    86.1,3259       49556,230923,280631,201738
g26531986t1  2361    1986    297900  126     92      uniq    92.0,2360       0,217168,217172,216981
g26532704t1  681     574     86100   126     82      uniq    79.3,675        0,53539,53539,53539
g26533243t1  1260    1216    182400  144     112     dupl    172.9,1257      75430,141900,217330,51080
g26533276t1  780     404     60600   77      66      uniq    66.4,778        0,51629,51629,51629
g26533533t1  1209    823     123450  102     86      uniq    86.4,1208       0,104362,104362,104362

=item try1

  .. dont need both C.nz, C.eq for uniq set, differ only for dup set, drop C.eq col
  
gunzip -c dropse20gnodes/dropse20cdsbusco_SRR11813283_2_bwa.sam.gz | \
 $gw/daphnia/gnodes_sam2ucgcov.pl -debug -maxdup 0.05 | less

Items=1048, C.LN/W=109,109.8 med/ave, C.Map/W=mcmw,87.8 med/ave for Wgenespan=1475412, Nread=1079901, Nrdlen=161985150

Gene_ID      Glen    Nread   Rdlen   C.LN    C.M     Uniq    C.nz            C.eq            Sdup,unq,tot,u=t
g117183176t1 618     611     91650   148     90      uniq    89.5,616        89.8,612        0,55116,55128,54945
g117183718t1 2094    1605    240750  114     102     uniq    102.3,2093      102.2,2073      0,213995,214032,211761
g117184007t1 1356    1120    168000  123     108     uniq    109.1,1350      109.6,1304      0,147263,147322,142864
g117184773t1 702     548     82200   117     92      uniq    93.1,700        93.2,680        0,65123,65143,63384
g13036350t1  7896    5683    852450  107     102     uniq    102.2,7885      102.3,7681      0,805780,806030,786006
g26531963t1  3261    2210    331500  101     87      dupl    101.8,3259      84.2,2949       49556,282139,331850,248447
g26531986t1  2361    2232    334800  141     104     uniq    104.7,2360      104.7,2352      0,247164,247174,246274
g26532704t1  681     648     97200   142     95      uniq    92.4,675        90.7,641        0,62333,62367,58142
g26533243t1  1260    1462    219300  174     136     dupl    196.7,1257      190.3,337       75430,171805,247239,64123

=cut
