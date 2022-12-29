#!/usr/bin/env perl
# cdshexprob.pl
# test code for scoring coding potential via cds hexamers
# update for mrna w/ cds offset info, calc code/nocode=utr likelyhood
# llHexcode{hex} = log((1+code)/(1+noncode)) as per cpat

use strict;
use Getopt::Long;

use constant MIN_MRNA_SAMPLE => 200; # 100 too low ?  
use constant ADD_STOPC => 1;
## FIXME compute HNoncode, HCoding from sample vals, ie HNoncode < 05% or 01% llscores if( llval < 0), HCoding ?> 0
use constant { HNoncode => -1e-4, HCoding => 1e-4 };  # hexamer score, dunno what range to call

my($NSTAT,$NSAMP,$MAXTEST)=(1000, $ENV{nsample}||5, $ENV{maxtest}||0 );
my $USE_BEST_MRNA=$ENV{bestaa}||0;
my $BEST_MINAA= undef; # set below ($USE_BEST_MRNA > 99) ? $USE_BEST_MRNA : 300;

my($ncds,$nhex,$nsamp,$avehex,$logave)=(0) x 9;
my $output= ""; # undef; # means stdout, or always file?
my($hexcode,$hexnoncode,%samp)= ({}, {});
use constant DO_CODONS => 1; # now only tabulate in aaqual.data
my($ncodon,$codon)= (0,{});
my $debug= $ENV{debug}||0;
my $use_fickett= 1; # option
my $only_fickett= 0; # skip hexcode/aaqualdata

my( $nmrna, $n_hexcode, $n_hexnoncode, 
    $llhexcode, $n_llhexcode, $llhex_isCoding, $llhex_isNoncode)=(0) x 9;
my($incds,$testcds,$aaqualdata); #= @ARGV; pre-getopts
my @testcds; # option list of test.cds files

my $optok= GetOptions(
  "incds|mrna|ref=s", \$incds, # require
  "testcds=s", \@testcds, # require, find  .. allow @testcds array ?
  "aaqualdata=s", \$aaqualdata, # aaqualf required ? w/ codepot; add mrna/cds/aa seq set??
  "output:s", \$output, # option
  "bestaa:i", \$BEST_MINAA, 
  "maxtest=i", \$MAXTEST,  
  "nsample=i", \$NSAMP,    
  "fickett!", \$use_fickett, 
  "onlyfickett!", \$only_fickett, 
  "debug!", \$debug, 
  );

if(@ARGV) { 
  $incds= shift @ARGV if(@ARGV and not ($incds or $aaqualdata)); #? drop? @testcds only?
  $testcds= shift @ARGV if(@ARGV and not $testcds);
  $aaqualdata= shift @ARGV if(@ARGV and not $aaqualdata); #?? reorder? mrna, qualdata, testcds1 testcds2 ..
  push @testcds,  @ARGV if( @ARGV);
}
$testcds= shift @testcds if(@testcds and not $testcds);

die "Evigene measure coding potential via CDS hexamers, draft version
usage: cdshexprob.pl validsample.cds testset.cds  OR 
  cdshexprob.pl -mrna valid.mrna | -incds valid.cds -test test.cds -aaqualdata valid.aaqual.data
   opts: -output codepot table, -bestaa, use valid.mrna with best aa quality,
         -maxtest=100 limit -test output
  -- calculate coding likelihood from CDS, UTR in mRNA, requires Evigene header type=mRNA; offs=50-500;
  cdshexprob.pl validsample_evg.mrna testset.cds  validsample_evg.aaqual.data
  -- use coding likelihood table in validsample_evg.aaqual.data
  cdshexprob.pl validsample_evg.aaqual.data othertestset.cds
\n" unless( $optok and ($incds or ($testcds and ($aaqualdata or $only_fickett)) ) ); #  and $testcds : option?

if(defined $BEST_MINAA){ $USE_BEST_MRNA ||= 1; }
$BEST_MINAA= ($BEST_MINAA > 69)? $BEST_MINAA : ($USE_BEST_MRNA > 99) ? $USE_BEST_MRNA : 300;
$use_fickett=2 if ($only_fickett); $only_fickett=0 unless($use_fickett);

if($incds =~ /aaqual.data$/) { $aaqualdata= $incds; $incds=""; }
elsif($aaqualdata){ $aaqualdata .= ".aaqual.data" unless($aaqualdata =~ /aaqual.data$/); }
elsif($incds =~ /\.mrna/) { ($aaqualdata=$incds) =~ s/.mrna.*/.aaqual.data/; }

my $OUTH=undef;
if(lc($output) eq 'stdout' or $output eq '-') { $OUTH= *STDOUT; }
elsif(defined $output) {
  unless($output) { ($output= $testcds) =~ s/\.\w+$//; $output ||="cdshexprob"; $output.=".codepot"; }
  rename($output,"$output.old") if(-f $output);
  open($OUTH, '>', $output);
} else {
  $OUTH= *STDOUT; # write to $output
}

if($aaqualdata and -s $aaqualdata) {
  ($n_llhexcode, $llhexcode)= read_loglikecoding($aaqualdata);
} 
if($aaqualdata and not $incds and not $n_llhexcode) {
  die "missing loglikecoding table in aaqualdata: $aaqualdata\n";
}

($ncds,$nmrna,$nhex)= ($incds and not $n_llhexcode) ? readValidCds($incds) : (0,0,0); 
  # with incds == mrna,  make_mrna_loglikecoding(); print_loglikecoding(aaqualdata);

# valid sample stats, only from readValidCds() ......
if($ncds) { 
  my($nsum,$hsum)=(0) x 3; 
  for my $id (sort keys %samp) {  my $cds=$samp{$id};
    my($phex,$llc)= iscode_hexamer($id,$cds); $hsum +=$phex; $nsum++; 
  } 
  $avehex= ($nsum>0) ? $hsum/$nsum : 0; 
  $logave= log(1+$avehex);
  
  my $ns=0; my $icds=0;
  for my $id (sort keys %samp) { 
    if(++$icds % 50 == 33 and $ns++ < $NSAMP) { my $cds=$samp{$id};
    my($phex,$llc)= iscode_hexamer($id,$cds); 
    my($fkval,$fkcode) = ($use_fickett) ? iscode_fickett($id, $cds):(0,'u');
    printhex($id,$phex,$llc,'sample',$fkval); 
    }
  } 
}

for my $tcds ($testcds, @testcds) {
  measure_testcds($tcds,$MAXTEST) if($tcds);
} 

close($OUTH) if($output);

#-----------

sub measure_testcds {
  my($testcds, $maxtest)= @_;
  my($id,$cds,$ntest)=(0,0,0); 
  warn "measure_testcds: $testcds, output to ", ($output||"STDOUT"), "\n" if $debug;
  
  if($testcds =~ /.gz/){ open(TCDS,"gunzip -c $testcds |") or die "measure_testcds: failed $testcds"; } 
  else { open(TCDS,$testcds) or die "measure_testcds: failed $testcds";}
  
  while(<TCDS>) {
    if(/^>(\S+)/){ 
      if($cds) { 
        my($phex,$llc)= ($only_fickett) ? (0,0) : iscode_hexamer($id,$cds); 
        my($fkval,$fkcode) = ($use_fickett) ? iscode_fickett($id, $cds):(0,'u');
        # $fkval .= substr($fkcode,0,1); #?
        printhex($id,$phex,$llc,'test',$fkval); 
        }
      $id=$1; $cds=""; $ntest++; last if($maxtest and $ntest > $maxtest); }
    else { chomp; $cds.=$_; }
  } close(TCDS);
  if($cds) { 
    my($phex,$llc)= ($only_fickett) ? (0,0) : iscode_hexamer($id,$cds); 
    my($fkval,$fkcode) = ($use_fickett) ? iscode_fickett($id, $cds):(0,'u');
    # $fkval .= substr($fkcode,0,1); #?
    printhex($id,$phex,$llc,'test',$fkval); 
    }
}

sub readValidCds {
  my($incds)= @_;
  my($id, $cds, $ismrna, $offs, $aaq, $clen,)= (0) x 9;
  $ncds=$nmrna=$nhex= 0; # globs
  warn "readValidCds: $incds\n" if($debug);
  
  ## FIXME may not have evg format mrna here, need APPcdshexprob to use name.tr + name.cds hdrs
  my($ncdsinfo, %cdsinfo)=(0,0);
  if( $incds !~ /\.(mrna|cds|aa)/ ) { # may be unoriented transcripts tr/cdna/fasta w/ evg.cds
    (my $cdsf=$incds) =~ s/\.\w+$/.cds/;
    warn "readValidCds: '$incds' is not mRNA? checking $cdsf\n" if($debug);
    if( -f $cdsf and open(my $hin, $cdsf)) { 
      while(<$hin>) { if(/^>(\S+)/) { my $id=$1; 
        ($aaq)= (m/aalen=([^;\s]+)/)?$1:0; ($clen)= (m/clen=(\d+)/)?$1:0; ($offs)= (m/offs=([^;\s]+)/)?$1:0; 
        if($offs and $aaq) { $cdsinfo{$id}= "$aaq\t$clen\t$offs"; $ncdsinfo++; }
        } 
      } close($hin); 
    }
  }
  
  my $ok;
  if($incds =~ /.gz/){ $ok= open(INCDS,"gunzip -c $incds |"); }
  else { $ok=open(INCDS,$incds); }
  die "readValidCds fail: $incds" unless($ok);; # evg.openRead($incds) to use seq.gz
  while(<INCDS>) {
    if(/^>(\S+)/){ 
    
      if($cds and $id){ 
        if($ismrna and $offs) {
          my($nxc,$nxu,$cdsof)= count_mrna_hexamer($id,$cds,$offs,$clen,$aaq); 
          if($nxc){ $nhex += $nxc; $nmrna++; $cds= $cdsof; }
        } else {
          $nhex += count_cds_hexamer($id,$cds); 
        }
        $ncds++; 
        if($ncds % 5 == 3 and $nsamp++ < $NSTAT) { $samp{$id}=$cds; }
      } 
      
      $id=$1; $cds=""; 
      $ismrna=$offs=$aaq=$clen=0; # strand also? expect only +strand for mRNA
      if(/type=mRNA/){ $ismrna=1; 
        ($aaq)= (m/aalen=([^;\s]+)/)?$1:0; ($clen)= (m/clen=(\d+)/)?$1:0; ($offs)= (m/offs=([^;\s]+)/)?$1:0; 
      } elsif($ncdsinfo) {
        ($aaq,$clen,$offs)= split"\t",$cdsinfo{$id}; $ismrna=1 if($offs);
      }
      
    } else { chomp; $cds.=$_; }
  } close(INCDS);
  
  if($cds){ 
    if($ismrna and $offs) {
      my($nxc,$nxu,$cdsof)= count_mrna_hexamer($id,$cds,$offs,$clen,$aaq); 
      if($nxc){ $nhex += $nxc; $nmrna++; $cds= $cdsof; }
    } else {
      $nhex += count_cds_hexamer($id,$cds); 
    }
   $ncds++; 
  }
  warn "readValidCds done: ncds=$ncds,nmrna=$nmrna,nhex=$nhex\n" if($debug);

  if($nmrna >= MIN_MRNA_SAMPLE) { 
    my $llmakelog;
    ($n_llhexcode,$llmakelog)= make_mrna_loglikecoding(); 
    print_loglikecoding($aaqualdata, $n_llhexcode, $llhexcode,$llmakelog) 
      if($n_llhexcode and $aaqualdata); 
    # globals ( $nmrna, $n_hexcode, $n_hexnoncode, $llhexcode, $n_llhexcode)=(0) x 9;
  }
    
  return($ncds,$nmrna,$nhex);
}


my($printhdr,$fmtHex)=(0,0,0);

sub printhex { 
  my($id,$phex,$llc,$tag,$cdficket)=@_;
  my $dhex= log(1+$phex) - $logave;
  $cdficket ||= 0;
  
  unless($printhdr) {
    my(@f,@h,@l);
    my $idw= length($id); $idw= 5 + 5*int($idw/5); $idw=15 if($idw<15);
    @f=("%-$idw"."s"); @h=("%-$idw"."s"); @l=('Ident');
    if($only_fickett) { }
    elsif($n_llhexcode and not $ncds) { push @f, "%+.4f";  push @h, "%6s"; push @l, 'llCode'; }
    elsif($n_llhexcode) { push @f, "%+.4f","%+.4f","%+.4f";  push @h, "%4s","%4s", "%6s"; push @l, qw(pHex dHex llCode); }
    else { push @f, "%+.4f","%+.4f";  push @h, "%4s","%4s"; push @l, qw(pHex dHex); }
    if($use_fickett) { push @f, "%.3f"; push @h, "%s"; push @l, 'Fickett'; }
    push @f, "%s"; push @h, "%s"; push @l, 'Type';
    $fmtHex = join("\t", @f) . "\n";
    my $fmtHdr = join("\t", @h) . "\n";
    printf $OUTH $fmtHdr, @l; $printhdr=1;
  }

  $cdficket= $tag unless($use_fickett);
  if($only_fickett) {  
  printf $OUTH $fmtHex, $id, $cdficket, $tag;
  } elsif($n_llhexcode and not $ncds) { 
  printf $OUTH $fmtHex, $id, $llc, $cdficket, $tag;
  } elsif($n_llhexcode) { #? drop dHex
  printf $OUTH $fmtHex, $id, $phex, $dhex, $llc, $cdficket, $tag;
  } else {
  printf $OUTH $fmtHex, $id, $phex, $dhex, $cdficket, $tag;
  }
  
  # if(0) {    
  #   if($n_llhexcode and not $ncds) { 
  #   printf $OUTH "%-20s\t%6s\t%s\t%s\n", qw(Ident llCode cdFicket Type) unless($printhdr++);
  #   printf $OUTH "%-20s\t%+.4f\t%.3f\t%s\n", $id, $llc,$cdficket, $tag;
  #   } elsif($n_llhexcode) { #? drop dHex
  #   printf $OUTH "%-20s\t%4s\t%4s\t%6s\t%s\t%s\n", qw(Ident pHex dHex llCode cdFicket Type) unless($printhdr++);
  #   printf $OUTH "%-20s\t%.4f\t%+.4f\t%+.4f\t%.3f\t%s\n", $id, $phex, $dhex, $llc, $cdficket, $tag;
  #   } else {
  #   printf $OUTH "%-20s\t%4s\t%4s\t%s\t%s\n", qw(Ident pHex dHex cdFicket Type) unless($printhdr++);
  #   printf $OUTH "%-20s\t%.4f\t%+.4f\t%.3f\t%s\n", $id, $phex, $dhex, $cdficket, $tag;
  #   }
  # }  
}

sub iscode_hexamer {
  my($id, $cds)=@_;
  my $n=length($cds); $cds=uc($cds);
  my($val,$nx,$nc, $llval)=(0) x 9;
  for(my $i=0; $i <= $n-6; $i+=3) {
    my $hex=substr($cds,$i,6); $nc++;
    if(my $h= $hexcode->{$hex}) { $val += $h; $nx++; } # $h= $h/$nhex;  below
    if($n_llhexcode) { my $h= $llhexcode->{$hex}||0; $llval += $h; }	    
    }
  # my $ftype= ($val <= HNoncode) ? "Noncode" : ($val >= HCoding)? "Code" : "Unknown";
  #o my $aval= 100 * $val / $nc;
  my $aval= ($nhex < 1 or $nc < 1) ? 0 : 100 * $val / ($nhex*$nc);
  return($aval,$llval,$nx,$nc);
}

sub count_hexamerh {
  my( $hexhash, $id, $seq, $iscds)=@_;
  my $n=length($seq);  $seq=uc($seq); $iscds||=0;
  my($nx)=(0);
  for(my $i=0; $i <= $n-6; $i+=3) {
    my $hex= substr($seq,$i,6);
    $hexhash->{$hex}++; $nx++;
    if(DO_CODONS and $iscds) { 
      my $cod=substr($hex,0,3); unless($cod =~ m/N/) { $codon->{$cod}++; $ncodon++; } 
      }
    }
  return $nx;
}

sub count_cds_hexamer {
  my($id, $cds)=@_;
  return count_hexamerh($hexcode,$id,$cds,1);
  # my $n=length($cds); $cds=uc($cds);
  # my($val,$nx)=(0,0);
  # for(my $i=0; $i<$n-6; $i+=3) {
  #   my $hex=substr($cds,$i,6);
  #   $hexcode->{$hex}++; $nx++;
  #   if(DO_CODONS) { my $cod=substr($hex,0,3); $codon->{$cod}++; $ncodon++; }
  #   }
  # return $nx;
}

sub count_mrna_hexamer {
  my($id, $mrna, $offs, $clen, $aaq)=@_;  
  $offs||="";
  my($nxcds, $nxutr)=(0,0);
  my($cb,$ce)= $offs =~ m/(\d+).(\d+)/;
  return(0) unless($cb and $ce > 0);
  
  my $or=0; if($cb > $ce){ ($cb,$ce,$or)= ($ce,$cb,-1); }
  if($or<0) { } # fixme? expect only mrna +or
  
  my $cw= 1+$ce-$cb;
  return (0) if($USE_BEST_MRNA and $cw/3 >= $BEST_MINAA and $aaq =~ /complete/ and $aaq !~ /utrbad/);
  if(ADD_STOPC and $aaq and $aaq=~/complete|partial5/ and $ce+3 <= $clen) { $cw += 3; } # count stop codons
  my $cds= substr($mrna, $cb-1, $cw);
  ($nxcds)= count_hexamerh($hexcode,$id,$cds,1);

  # mrna_loglikecoding nmrna=60949, nxcode=52590452, nxutr=37020163, hexs=4108,nocds=107,noutr=8
  # this way drops 1/2 of utronly NNN codes from useless spacer, gives nocds ~ noutr
  if($cb > 6) {  $nxutr += count_hexamerh($hexnoncode,$id,substr($mrna, 0, $cb-1)); }
  if($ce < $clen - 6) { $nxutr += count_hexamerh($hexnoncode,$id,substr($mrna, $ce)); }

  #   my $utr="";
  #   # my $us = ($cb-1) % 3; $us=3 if($us==0); my $usp= substr('NNNN',0,$us);
  #   # $utr= substr($mrna, 0, $cb-1) . $usp . substr($mrna, $ce);
  #   if($cb > 6) { 
  #     my $us = ($cb-1) % 3; $us=3 if($us==0); #? drop excess NNN ?
  #     my $usp= substr('NNNN',0,$us);
  #     $utr .= substr($mrna, 0, $cb-1) . $usp; 
  #     }
  #   if($ce < $clen - 6) { $utr .= substr($mrna, $ce); }
  #   ($nxutr)= count_hexamerh($hexnoncode,$id,$utr);
  
  
  $n_hexcode += $nxcds; # global counts
  $n_hexnoncode += $nxutr;
  return($nxcds, $nxutr, $cds);
}

sub print_loglikecoding {
  my($aaqualdata, $n_llhexcode, $llhexcode, $llmakelog)= @_;
  my($nout)=(0);
  if(-f $aaqualdata){ rename($aaqualdata,"$aaqualdata.old"); }
  warn "print_loglikecoding n=$n_llhexcode of $incds\n";
  open(LLO,'>',$aaqualdata) or die "print_loglikecoding $aaqualdata";
  print LLO "#loglikecoding n=$n_llhexcode of $incds\n";
  print LLO "#$llmakelog\n";

  if($n_hexcode > 0 and $n_hexnoncode > 0) {
    print LLO join("\t","#Hexcode","loglikeCDS","pCDS","pUTR")."\n";
    for my $h (sort keys %$llhexcode) { 
      my $llc= $llhexcode->{$h};
      my $hc= ($hexcode->{$h}||0) / $n_hexcode;      
      my $hn= ($hexnoncode->{$h}||0) / $n_hexnoncode; 
      printf LLO "%s\t%.4g\t%.3g\t%.3g\n",$h,$llc, $hc, $hn; $nout++;
    }
  } else {    
    print LLO "#Hexcode\tloglikeCDS\n";
    for my $h (sort keys %$llhexcode) { printf LLO "%s\t%.4g\n",$h,$llhexcode->{$h}; $nout++; }
  }
  print LLO "#Hexcode_END\n";
  
  if(DO_CODONS and $ncodon > 0) {
    my @codon= sort keys %$codon; my $nc= @codon;
    print LLO "\n#Codon_Table ncodon=$nc\n";
    print LLO "#Codon\tFreq\tCount\n";
    for my $cod (@codon) {
      #cut in counter: next if($cod =~ m/N/); #  dont print gapped codon
      my $nc= $codon->{$cod}; 
      printf LLO "%s\t%.4g\t%d\n", $cod, $nc/$ncodon, $nc; 
    }
    print LLO "#Codon_END\n";
  }
  
  close(LLO); 
  return($nout);
}



sub read_loglikecoding {
  my($aaqualdata)= @_;
  open(LLO,$aaqualdata) or die "read_loglikecoding $aaqualdata";
  my($nhex,$inv,%llhexcode)= (0,0);
  while(<LLO>) {
    if(/^#Hexcode\tloglikeCDS/) { $inv=1; }
    elsif(/^#Hexcode_END/) { $inv=0; last; }
    elsif($inv and /^\w/){ my($h,$v)=split; $llhexcode{$h}=$v if($v=~/\d/); $nhex++; }
  }
  warn "read_loglikecoding n=$nhex of $aaqualdata\n";
  return($nhex,\%llhexcode);
}

sub make_mrna_loglikecoding {
  # llHexcode{hex} = log((1+code)/(1+noncode)) as per cpat
  # need frag counts per hex hash for freq
  my @hcode= sort keys %$hexcode;
  my @hnonc= sort keys %$hexnoncode;
  use constant HMIN => 300; # expect about 4,000 hexcodes for each
  $llhexcode= {}; # clear global
  my($nhex, $nhcode, $nhutr)=(0) x 9;
  if(@hcode > HMIN and @hnonc > HMIN) {
    my %h= map{ $_,1 } (@hcode,@hnonc);
    my @h=sort keys %h; $nhex= @h;
    for my $h (@h) {
      my $hc= ($hexcode->{$h}||0); $nhcode++ if($hc>0);         
      my $hn= ($hexnoncode->{$h}||0); $nhutr++ if($hn>0);  
      my $llc= log(1 + $hc / $n_hexcode) - log(1 + $hn / $n_hexnoncode); 
      $llhexcode->{$h}= $llc;
    }
  my $ncz= $nhex - $nhcode; my $nuz= $nhex - $nhutr;
  my $llmakelog="make_llcds nmrna=$nmrna, nxcode=$n_hexcode, nxutr=$n_hexnoncode, hexs=$nhex,nocds=$ncz,noutr=$nuz";
  warn "$llmakelog\n";
  return($nhex,$llmakelog);
  }
  return(0);
}

## UPD 20.02.20 Fickett code score

sub _mina { my @v=@_; my $m=shift @v; for(@v){ $m=$_ if($_<$m); }; return $m; }
sub _maxa { my @v=@_; my $m=shift @v; for(@v){ $m=$_ if($_>$m); }; return $m; }

use constant { FNoncode => 0.74, FCoding => 0.95 };  
our (%position_prob, %position_weight, @position_para, %content_prob, %content_weight, @content_para);

# Fickett TESTCODE data; NAR 10(17) 5303-531
BEGIN{ # dangit perl..
our %position_prob =(
'A'=>[0.94,0.68,0.84,0.93,0.58,0.68,0.45,0.34,0.20,0.22],
'C'=>[0.80,0.70,0.70,0.81,0.66,0.48,0.51,0.33,0.30,0.23],
'G'=>[0.90,0.88,0.74,0.64,0.53,0.48,0.27,0.16,0.08,0.08],
'T'=>[0.97,0.97,0.91,0.68,0.69,0.44,0.54,0.20,0.09,0.09]
);
our %position_weight=('A'=>0.26,'C'=>0.18,'G'=>0.31,'T'=>0.33);
our @position_para  = (1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0);

our %content_prob=(
'A'=>[0.28,0.49,0.44,0.55,0.62,0.49,0.67,0.65,0.81,0.21],
'C'=>[0.82,0.64,0.51,0.64,0.59,0.59,0.43,0.44,0.39,0.31],
'G'=>[0.40,0.54,0.47,0.64,0.64,0.73,0.41,0.41,0.33,0.29],
'T'=>[0.28,0.24,0.39,0.40,0.55,0.75,0.56,0.69,0.51,0.58]
);
our %content_weight=('A'=>0.11,'C'=>0.12,'G'=>0.15,'T'=>0.14);
our @content_para  =(0.33,0.31,0.29,0.27,0.25,0.23,0.21,0.17,0);
}

sub fickett_position_prob {
  my($value, $base, $debug)= @_;
  return 0 if($value<0);
  for my $i (0..$#position_para) {
    if ($value >= $position_para[$i]) { return $position_prob{$base}->[$i] * $position_weight{$base}; }
  }
  return 0;
}

sub fickett_comp_prob {
  my($value, $base)=@_;
  return 0 if($value<0);
  for my $i (0..$#content_para) {
     if ($value >= $content_para[$i]) { return $content_prob{$base}->[$i] * $content_weight{$base}; }
  }
  return 0;
}


sub iscode_fickett {
  my($id,$cds)=@_;
  $cds= uc($cds); # ensure, fix data bug below
  my @cds= split"",$cds;
  my %ACGT=(A=>1,C=>2,G=>3,T=>4);
  my @kACGT=qw(A C G T); # (0=>'A',1=>'C',2=>'G',3=>'T');

  my $ndebug=0;
  my $ncd=0;
  my($c,$f,$i,$j,$k,$v);
  my(@scores);
  for($i=0; $i<@cds; $i++) {
    if( $k=$ACGT{$cds[$i]} ) { $scores[$k-1][ $i % 3 ]++; $ncd++; }
  }
  my ($tsum,@bsum);
  for $k (0,1,2,3) {
    for $j (0,1,2) { my $s= $scores[$k][$j]; $bsum[$k]+=$s; $tsum+=$s; }
  }
  return(0,0,0) unless($tsum); # data BUG tsum == 0 from lc(cds)
  my @comp= map{ $_/$tsum } @bsum;
  my @posi;
  for $k (0,1,2,3) {
    my $vmax= _maxa(@{$scores[$k]});
    my $vmin= _mina(@{$scores[$k]});
    $posi[$k]= $vmax / (1+$vmin);
  }  

  my $fscore=0;
  for $k (0,1,2,3) {  
    my $base=$kACGT[$k];
    my ($pk,$ck);
    my $tdebug=0; # ($debug and $ndebug < 9)?1:0;
    my $pp= fickett_position_prob(($pk=$posi[$k]), $base, $tdebug); 
    my $cp= fickett_comp_prob(($ck=$comp[$k]), $base, $tdebug); 
    $fscore += $pp + $cp;
  }
  
  my $ftype= ($fscore <= FNoncode) ? "Noncode" : ($fscore >= FCoding)? "Code" : "Unknown";  
  return($fscore,$ftype,$ncd);
}

__END__

=item codon table

  env maxtest=10 $evigene/scripts/prot/cdshexprob.pl arathv4j.oklong.mrna arathv4j.oklong.cds  arathv4j_longcds.aaqual.data | head -30 
  make_llcds nmrna=24159, nxcode=14108968, nxutr=3880219, hexs=4109,nocds=106,noutr=13
  print_loglikecoding n=4109 of arathv4j.oklong.mrna
  
  Ident               	pHex	dHex	llCode	Type
  Arath4bEVm000518t1  	0.0445	-0.0000	+0.1511	sample
  Arath4bEVm001092t1  	0.0467	+0.0021	+0.0990	sample
  Arath4bEVm001707t1  	0.0471	+0.0025	+0.0772	sample
  Arath4bEVm002464t1  	0.0431	-0.0014	+0.0553	sample
  Arath4bEVm003015t1  	0.0439	-0.0006	+0.0674	sample
  
  Arath4bEVm000611t5  	0.0460	+0.0014	+0.1394	test
  Arath4bEVm000611t4  	0.0496	+0.0049	+0.1729	test
  Arath4bEVm000611t1  	0.0483	+0.0036	+0.2832	test
  Arath4bEVm000101t1  	0.0482	+0.0035	+0.2556	test
  Arath4bEVm004121t1  	0.0501	+0.0053	+0.0959	test
  Arath4bEVm006596t1  	0.0388	-0.0055	+0.0457	test
  Arath4bEVm006548t1  	0.0378	-0.0064	+0.0391	test
  Arath4bEVm006182t1  	0.0375	-0.0068	+0.0340	test
  Arath4bEVm000894t1  	0.0446	+0.0000	+0.1236	test
  Arath4bEVm005121t1  	0.0481	+0.0034	+0.0642	test

  head arathv4j_longcds.aaqual.data
  #loglikecoding n=4109 of arathv4j.oklong.mrna
  #make_llcds nmrna=24159, nxcode=14108968, nxutr=3880219, hexs=4109,nocds=106,noutr=13
  #Hexcode	loglikeCDS	pCDS	pUTR
  AAAAAA	-0.003281	0.000533	0.00382
  AAAAAC	-0.0005982	0.000533	0.00113
  AAAAAG	-0.0004546	0.000793	0.00125
  AAAAAT	-0.0008231	0.000506	0.00133
  AAAACA	-0.0005557	0.00059	0.00115
  AAAACC	-0.0002657	0.000441	0.000707
  AAAACG	-0.0001027	0.00022	0.000323
  
  tail -n75 arathv4j_longcds.aaqual.data | head
  
  #Codon_Table ncodon=71
  #Codon	Freq	Count
  AAA	0.03158	445609
  AAC	0.02005	282863
  AAG	0.03139	442893
  AAT	0.02421	341531

From Arath ref mrna, commonest
  GAT	0.03799	536055
  GAA	0.0363	512158
  GAG	0.03179	448536
  AAA	0.03158	445609
  AAG	0.03139	442893
  GTT	0.0272	383796
  GCT	0.02711	382497
  TCT	0.02603	367251
  CTT	0.02511	354235
  AAT	0.02421	341531
  ATG	0.02389	337014  # start
  GGA	0.02328	328518
  TTT	0.02294	323702

Note stop codons are not in CDS offset, but some are in data
  TGA	4.465e-06	63
  TAA	2.41e-06	34
  TAG	2.41e-06	34
ADD_STOPC counts for nmrna=24159
  TGA	0.000677	9567
  TAA	0.0005461	7717
  TAG	0.0003152	4454

Arath ref mrna has some odd nucleic bases: K S W M Y
  GAK	1.418e-07	2
  GKT	1.418e-07	2
  KTG	1.418e-07	2
  STA	1.418e-07	2
  AAW	7.088e-08	1
  GCM	7.088e-08	1
  GYG	7.088e-08	1

=cut

=item ref tests plant

plant  $weed/tr2aacds_test1908f/aaeval/arathv4j.oknoncode_cds.llcode.tab
  Ident               	pHex	dHex	llCode	Type
  Arath4bEVm000518t1  	0.0446	-0.0001	+0.1566	sample
  Arath4bEVm001092t1  	0.0468	+0.0020	+0.1045	sample
  Arath4bEVm001707t1  	0.0473	+0.0024	+0.0827	sample
  Arath4bEVm002464t1  	0.0433	-0.0014	+0.0593	sample
  Arath4bEVm003015t1  	0.0440	-0.0006	+0.0715	sample
  Arath4bEVm025809t4  	0.0405	-0.0040	+0.0053	test
  Arath4bEVm026542t1  	0.0446	-0.0001	+0.0004	test
  Arath4bEVm023916t2  	0.0344	-0.0099	-0.0014	test
  Arath4bEVm023327t2  	0.0416	-0.0030	-0.0026	test
  1310 of 3033, 43%, have neg llval for those called Noncode by old table
  
ref long cds neg values,
pick HNoncode cut-off somewhere in 32 neg llCode vals, of 11450, bottom 0.27 % of llvals
cut -f1,2,4,5 $weed/tr2aacds_test1908f/aaeval/arathv4j.oklong_cds.llcode.tab | grep -v sample | sort -k3,3n | grep -c '       -'
  32 of 11448, 0.27% long cds have neg llval
    Ident               	pHex	dHex	llCode	Type
  Arath4bEVm000136t1  	0.0465	-0.3693	test
  Arath4bEVm002806t1  	0.0384	-0.1438	test
  Arath4bEVm006642t1  	0.0375	-0.0893	test
    ..
  Arath4bEVm006117t1  	0.0399	-0.0224	test
  Arath4bEVm003594t1  	0.0394	-0.0211	test
  Arath4bEVm006610t1  	0.0336	-0.0181	test  <
  Arath4bEVm005548t3  	0.0391	-0.0113	test  <
  Arath4bEVm002054t1  	0.0422	-0.0111	test  <
  Arath4bEVm004720t1  	0.0375	-0.0108	test  <  ?
  Arath4bEVm003334t1  	0.0380	-0.0100	test  <  ? cut-off
  Arath4bEVm002164t1  	0.0399	-0.0089	test


=item ref tests human

human sra2genes_tr19human/tr2aacds_test1908f/try1912v4i/
  grep -c '>' try1912v4i.oklong.mrna try1912v4i.oknoncode.cds
  try1912v4i.oklong.mrna:60949  # aamin=300
  try1912v4i.oknoncode.cds:2438
  
  $evigene/scripts/prot/cdshexprob.pl try1912v4i.oklong.mrna try1912v4i.oknoncode.cds | head -10
  mrna_loglikecoding nmrna=60949, nxcode=52529503/4542 miss, nxutr=37118262/6 miss, llhexcode=8442
  or
  mrna_loglikecoding nmrna=60949, nxcode=52529503, nxutr=37118262, hexs=8442,nocds=4542,noutr=6
  
  Ident               	pHex	dHex	llCode	Type
  Homsap4aEVm000728t1 	0.0498	+0.0033	+0.2360	sample
  Homsap4aEVm001668t1 	0.0414	-0.0047	+0.1098	sample
  Homsap4aEVm002605t1 	0.0508	+0.0043	+0.2042	sample
  Homsap4aEVm003500t2 	0.0473	+0.0009	+0.1440	sample
  Homsap4aEVm004268t1 	0.0441	-0.0022	+0.0880	sample
  Homsap4aEVm020342t2 	0.0393	-0.0067	-0.0439	test
  Homsap4aEVm017640t2 	0.0325	-0.0133	-0.0097	test
  Homsap4aEVm015398t5 	0.0309	-0.0149	-0.0110	test
  Homsap4aEVm022415t2 	0.0306	-0.0151	-0.0019	test

cat try1912v4i.oknoncode_cds.llcode.tab | grep -v sample  | cut -f1,4 | sort -k2,2n | grep -c '  -'
  848 of 2438, 35%, ref called Noncode by old codepot table are still -llval

ref long cds neg vals

  cat try1912v4i.oklong_cds.llcode.tab | grep test | cut -f1,4 | sort -k2,2n | grep -c '   -'
  372 of 60949 = 0.61 %

  cat try1912v4i.oklong_cds.llcode.tab | grep test | cut -f1,2,4 | sort -k3,3n | head
  Ident               	pHex	  llCode
  Homsap4aEVm006313t1 	0.0209	-0.0588
  Homsap4aEVm003805t5 	0.0325	-0.0378
  Homsap4aEVm008804t1 	0.0335	-0.0350
  Homsap4aEVm004874t2 	0.0335	-0.0337
  Homsap4aEVm000531t1 	0.0380	-0.0299
  Homsap4aEVm004811t2 	0.0339	-0.0283
  Homsap4aEVm004874t1 	0.0340	-0.0280
  Homsap4aEVm003805t4 	0.0341	-0.0277
  Homsap4aEVm012440t1 	0.0311	-0.0274
  Homsap4aEVm008002t1 	0.0344	-0.0267
    after -100
  Homsap4aEVm007807t1 	0.0347	-0.0102
  Homsap4aEVm007659t3 	0.0351	-0.0101
  Homsap4aEVm011635t2 	0.0320	-0.0101
  Homsap4aEVm004591t4 	0.0373	-0.0100
  Homsap4aEVm007408t15	0.0410	-0.0100
  Homsap4aEVm013626t1 	0.0348	-0.0100
  Homsap4aEVm004940t1 	0.0385	-0.0099
  Homsap4aEVm007408t14	0.0410	-0.0099
  Homsap4aEVm007807t2 	0.0350	-0.0097
  Homsap4aEVm007807t3 	0.0349	-0.0097

=cut
