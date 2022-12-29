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
my $FAWIDTH= 100; # orig: 60; 
my $MINSCLEN= 200; # of scaffolds
my $MINCTLEN= 50; # of split contigs
my $debug=1;
my($chrasm,$agp,$contigs,$ninsc,$noksc,$ncout,$sseq,$sgap)= (0) x 10;

my $optok= &GetOptions (
  "chrfasta|input=s" => \$chrasm, # input chr dir OR all.fasta
  "output|agpout=s" => \$agp, # output(s), STDOUT is option for agp, not for ctg.fasta
  "contigs|fasta=s" => \$contigs, # output(s), STDOUT is option for agp, not for ctg.fasta
  "MINGAP=i"=>\$MINGAP,   
  "MINSCLEN=i"=>\$MINSCLEN,   
  "MINCTLEN=i"=>\$MINCTLEN,   
  "FAWIDTH=i"=>\$FAWIDTH,   
  "debug!"=>\$debug, 
  );
die "usage: chrfasta2agp.pl -agp asm_contigs.agp -contigs asm_contigs.fa -in assembly.fa 
  OR:  chrfasta2agp.pl < assembly.fasta > asm_contigs.agp 
  -mingap=$MINGAP -minsclen=$MINSCLEN -minctlen=$MINCTLEN -fawidth=$FAWIDTH
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

sub putsc {
  my($scafid,$scafseq)=@_;
  my @ctgseqs;
  $ninsc++;
	# mask scaftigs shorter than -S threshold with "N"s
use constant USE_MINGAP => 1;
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
	$scafseq =~ s/^[nN]+//g;
	$scafseq =~ s/[nN]+$//g;

	# skip scaffold if length less than -s threshold
	my $scaflen = $scafseq =~ tr/ACGTacgt//;
	return(0) if $scaflen < $MINSCLEN;
  $noksc++;
  
if(USE_MINGAP){ 
  @ctgseqs = split /($GAPP+)/i, $scafseq;
} else {
  @ctgseqs = split /([Nn]+)/, $scafseq;
}

	my ($i,$x) = (0,0); 
	for my $ctgseq (@ctgseqs) {
		my $len = length $ctgseq;
    if ($len == 0){ $i++; next; }
		#agp1: object object_beg object_end part_number
		print $OUTH $scafid, "\t", $x + 1, "\t", $x + $len, "\t", $i + 1, "\t";
		if ($ctgseq =~ /^[nN]/) {
			#agp2: component_type gap_length gap_type linkage
			print $OUTH "N\t", $len, "\tscaffold\tyes\tpaired-ends\n";
			$sgap+=$len;
		} else {
			my $iid= 1 + int($i/2);
			my $ctgid = 'c'.$iid.'_'.$scafid ;
			#agp2: component_type component_id component_beg component_end orientation
			print $OUTH "W\t", $ctgid, "\t1\t", $len, "\t+\n"; 
			$ncout++; $sseq+=$len;
      if($contigs){
 			  $ctgseq =~ s/(.{$FAWIDTH})/$1\n/g; my($xb,$xe)=($x+1,$x+$len);
        chomp($ctgseq) if(substr($ctgseq, -1, 1) eq "\n");
			  print FASTA '>', $ctgid, " len=$len; orig=$scafid:$xb-$xe;\n", $ctgseq, "\n";
			}
		}
		$i++; $x += $len; 
	}
	return 1;
}
