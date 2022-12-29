#!/usr/bin/perl
# evigene/scripts/genoasm/chrfasta2agp.pl from abyss-fatoagp.pl
# Convert a FASTA file of scaffolds to a FASTA file of contigs and an AGP file.
# Written by Shaun Jackman <sjackman@bcgsc.ca>.
# dgg mods: input.fa wrapped to lines, MINgap opt to split at (eg. 10 NNN)
# update at evigene/scripts/genoasm/chrfasta2agp.pl

use strict;
use warnings;
use Getopt::Long ;

my $MINGAP= 10; ## require this many NNN to split, may want 1 default?
  my $MAXGAP= 0; # only for gapends ? limit max gap sizes
my $FAWIDTH= 100; # orig: 60; 
my $GAPENDW= 200; # gapends contig spans
my $MINSCLEN= 200; # of scaffolds
my $MINCTLEN= 50; # of split contigs
my $ENDTRIM=1; # opt? silently trim N+ at ends of scaf
my $debug=1;
my($chrasm,$agp,$contigs,$ninsc,$noksc,$ncout,$sseq,$sgap)= (0) x 10;
my $gapends= undef; # ** defined == 0 means do gapends !

my $optok= &GetOptions (
  "chrfasta|input=s" => \$chrasm, # input chr dir OR all.fasta
  "output|agpout=s" => \$agp, # output(s), STDOUT is option for agp, not for ctg.fasta
  "contigs|fasta=s" => \$contigs, # output(s), STDOUT is option for agp, not for ctg.fasta
  "gapends:i"=>\$gapends, 
  "MINGAP=i"=>\$MINGAP,  "MAXGAP=i"=>\$MAXGAP,   
  "MINSCLEN=i"=>\$MINSCLEN,   
  "MINCTLEN=i"=>\$MINCTLEN,   
  "FAWIDTH=i"=>\$FAWIDTH,   
  "debug!"=>\$debug, 
  );
die "usage: chrfasta2agp.pl -agp asm_contigs.agp -contigs asm_contigs.fa -in assembly.fa 
  OR:  chrfasta2agp.pl < assembly.fasta > asm_contigs.agp 
  -mingap=$MINGAP -minsclen=$MINSCLEN -minctlen=$MINCTLEN -fawidth=$FAWIDTH
  -gapends=$GAPENDW -contigs output is instead gap with surrounding contig seq
   Split genome assembly fasta to no-gap contigs, output as AGP and optional fasta seq\n"
  unless($optok);

my $inh;
if(@ARGV and not $chrasm){ $chrasm= shift(@ARGV); }
if($chrasm) { 
  if($chrasm =~ /\.gz/){ open($inh,"gunzip -c $chrasm |") or die "reading $chrasm\n"; }
  else { open($inh,$chrasm) or die "reading $chrasm\n"; }
} else { $inh=*STDIN; } # fixme perl <> 

my $OUTH=*STDOUT;
if($agp) { open $OUTH, ">$agp" or die "error: $agp: $!\n" if $agp; }
if($contigs) { open FASTA, ">$contigs" or die "error: $contigs: $!\n"; }

use constant USE_MINGAP => 1;
my $DOGAPENDS=0; # also require $contigs
if(defined $gapends){ $DOGAPENDS=1; $gapends=$GAPENDW unless($gapends>0); }

$MINGAP||=1; my $GAPP= "N" x $MINGAP;
my ($scafid,$scafseq)=("","");
$agp||=""; $contigs||="";

while (<$inh>) {
  chomp;
  if(/^>/){ 
    putsc($scafid,$scafseq) if($scafseq);
	  ($scafid, undef) = split ' ', $_, 2;
	  substr $scafid, 0, 1, '';
    $scafseq="";
  } elsif(/^\w/) { 
	  $scafseq .= $_;
  }
}
putsc($scafid,$scafseq) if($scafseq);

my($mbs,$mbg)= map{ int($_/100_000)/10 } ($sseq,$sgap); 
warn "#fa2agp nin=$ninsc, nok=$noksc, nout=$ncout, seq=$mbs.mb, gaps=$mbg.mb, to $agp, $contigs\n" if $debug;

sub putfa {
  my($ctgid, $ctgseq, $scafid, $xb, $xe)=@_;
  # my($xb,$xe)=($x+1,$x+$len);
  my $len= length $ctgseq;
  $ctgseq =~ s/(.{$FAWIDTH})/$1\n/g; 
  chomp($ctgseq) if(substr($ctgseq, -1, 1) eq "\n");
  print FASTA '>', $ctgid, " len=$len; orig=$scafid:$xb-$xe;\n", $ctgseq, "\n";
}

sub putsc {
  my($scafid,$scafseq)=@_;
  my @ctgseqs;
  $ninsc++;
	# mask scaftigs shorter than -S threshold with "N"s
if(USE_MINGAP){ 
	@ctgseqs = split /($GAPP+)/i, $scafseq;
} else {
	@ctgseqs = split /([Nn]+)/, $scafseq;
}
	foreach my $ctgseq (@ctgseqs) {
		next if $ctgseq=~m/^[nN]/;
		if (length($ctgseq) < $MINCTLEN) {
			$ctgseq = "N" x length($ctgseq);
		}
	}

	# rejoin and split to merge adjacent stretches of "N"s
	$scafseq = join '', @ctgseqs;
	return(0) unless $scafseq =~ /[^nN]/;

	# trim leading/trailing "N"s that may result  from masking short contigs (-S option)
	if($ENDTRIM) {
	$scafseq =~ s/^[nN]+//g;
	$scafseq =~ s/[nN]+$//g;
  }
  
	# skip scaffold if length less than -s threshold
	my $scaflen = $scafseq =~ tr/ACGTacgt//;
	return(0) if $scaflen < $MINSCLEN;
  $noksc++;
  
if(USE_MINGAP){ 
  @ctgseqs = split /($GAPP+)/i, $scafseq;
} else {
  @ctgseqs = split /([Nn]+)/, $scafseq;
}

	my ($i,$x,$lctg,$lgapw,$lgaps,$lx,$lctgx) = (0) x 9; 
	for my $ctgseq (@ctgseqs) {
		my $len = length $ctgseq;
    if ($len == 0){ $i++; next; }
		#agp1: object object_beg object_end part_number
		print $OUTH $scafid, "\t", $x + 1, "\t", $x + $len, "\t", $i + 1, "\t";
		if ($ctgseq =~ /^[nN]/) {
			#agp2: component_type gap_length gap_type linkage
			print $OUTH "N\t", $len, "\tscaffold\tyes\tpaired-ends\n";
			$sgap+=$len; $lgapw=$len; $lgaps= $ctgseq; 
		} else {
			my $iid= 1 + int($i/2);
			my $ctgid = 'c'.$iid.'_'.$scafid ;
			#agp2: component_type component_id component_beg component_end orientation
			print $OUTH "W\t", $ctgid, "\t1\t", $len, "\t+\n"; 
			$ncout++; $sseq+=$len;
			
      if($contigs){
        if($DOGAPENDS) { 
          my $llen= length $lctg; 
          use constant MINCTEND => 9;
          my $okend= (($MAXGAP>0 and $lgapw > $MAXGAP) or $lgapw < $MINGAP or $llen < MINCTEND or $len < MINCTEND) ? 0 : 1;
          if($okend) {
          my $aw= ($gapends > $llen) ? 0 : $llen - $gapends;
          my $ctga= substr($lctg,$aw);
          my $bw= ($gapends > $len)? $len : $gapends;
          my $ctgb= substr($ctgseq,0,$bw);
          my $gapid = 'gap'.($iid - 1).'_'.$scafid ;
          putfa( $gapid, $ctga . $lgaps . $ctgb, $scafid, $lctgx+$aw, $x+$bw);
          }
          
        } else {
          putfa( $ctgid, $ctgseq, $scafid, $x+1,$x+$len); 
        }
        $lctg= $ctgseq; $lctgx= $x;
			}
		}
		$i++; $lx= $x; $x += $len; 
	}
	return 1;
}
