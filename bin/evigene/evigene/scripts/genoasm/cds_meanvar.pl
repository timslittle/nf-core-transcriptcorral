#!/usr/bin/perl

# loctab meanvar
# cat dmag19skasm.cds_rdcovk27cb.loctab | \
# egrep '^#scaf|^LG5	[0-9]*	100	100	100	0' | \
# env median=1 ntave=0 icols=6,7,8 ./cds_meanvar.pl

my $DOMED= (defined $ENV{median}) ? $ENV{median} : 1;
my $ZEROS= $ENV{zero}||0;
my $doTRIM= $ENV{trim}||0;
my $showSUM= $ENV{sum}||0;
my $showRANGE= $ENV{range}||undef;  #? make default, ENV{var} replaces showSUM,range?
if($showSUM or $ENV{var}) { $showRANGE=0; } elsif( not defined $showRANGE){ $showRANGE=1; }
  $DOMED=1 if ($showRANGE);
my $title=$ENV{title}||"";
my $colhdr=$ENV{colhdr} || $ENV{header} || "scaf|ChrID|GeneID";
   # opt/default 1st row = colhdr, eg Rstats header=T

my $rotab= $ENV{rotab}||0; # rotate/transpose table, R: t(tab), cols>rows

# Item 6..11: k27cov kcn1 kcn2x rdcovt rdcovm rdcovu
my(@IVAR,@lvar);
BEGIN {
  if(my $ic= $ENV{icols}){ 
     @ic=split(",",$ic);  @iv=();
     for $i (@ic){ ($j,$k)=split"-",$i; push @iv,$j; while($j<$k){ push @iv, ++$j; } }
     @IVAR= @iv; # split(",",$ic);  # or icol-1 ?
  } else {  # new general default: all cols 1-n, skip 0 ?

  }
  # } else { @IVAR=(6,7,8); } # FIXME need opts for diff tables
  # @obIVAR=(6,7,10,11); # skip kcn2x rdcovt  
  # @oIVAR=(6,7,8,9,10,11); @oldIVAR=(6,7,8,9); 
} 

$tophdr= ($colhdr eq 1 or $colhdr eq "T");
while(<>) {
  my @v=split; 
  unless(@IVAR) {  @IVAR=( 1..$#v ); }
  if($tophdr) { @lvar=@v[@IVAR]; $tophdr=0; }
  elsif(/^#/){ @lvar=@v[@IVAR] if(/^#($colhdr)/); if(/^#title\s*(.+)$/){ $title=$1; }} 
  else{ sum( @v[@IVAR]);} 
}

END{ meanvar(); } 

sub sum{ my @s=@_; 
 my $nz=0; map{ $nz++ if($_); } @s; return unless($nz); # or $ZEROS ??
 $nt++; for $i (0..$#s){ $s=$s[$i]; if($s>0 or $ZEROS) {
   $sum[$i]+=$s; $nit[$i]++; $svar[$i]+=$s*$s; 
   if($DOMED){ push @{$sval[$i]},$s; }
 }
} 
return $nz; } 

sub trim {
  my($i,$pCut)=@_;
  if($pCut <= 0.0001){ $pCut = 0.05; }
  elsif($pCut == 1) { $pCut = 0.05; }
  elsif($pCut > 1){ $pCut = $pCut/100; }

  my($nscut)=(0);
  my @ss=sort{$a <=> $b} @{$sval[$i]};
  my $ns=@ss;
  my $nend= int($ns * $pCut);
  splice(@ss,-$nend); splice(@ss,0,$nend);
  $nscut= $ns - @ss; # my $icut= $ns - @ss; $nscut=$icut if($icut>$nscut);
  $sum[$i]= $svar[$i]= $nit[$i]= 0;
  for my $s (@ss) { $sum[$i]+= $s;  $nit[$i]++; $svar[$i]+=$s*$s; }
  $sval[$i]= \@ss;
  return($nscut,$pCut);
}

sub meanvar { 
my $tis="$title stats for nt=$nt"; 
$tis.=" trim=$doTRIM" if($doTRIM);
print "$tis\n";

my $lvar=($showRANGE)?"Range":($showSUM)?"Sum":"Var";

## drop this
# if($ENV{ntave}) { 
# print join("\t",qw(Item__ Mean SEM StDev),$lvar)."\n"; 
# for my $i (0..$#sum){ $lvar=$lvar[$i]||"it$i"; $s=$sum[$i]; 
#   $ave=$s/$nt; $var=($svar[$i] - $ave*$ave)/($nt-1); $sd=sqrt($var); $se=$sd/sqrt($nt); 
#   $var=$s if($showSUM);
#   printf "%6s\t%.2f\t%.2f\t%.2f\t%.2f\n",$lvar,$ave,$se,$sd,$var; } 
# }

@hdr=qw(Mean SEM Nitem StDev); push(@hdr,$lvar); unshift(@hdr,"Median") if($DOMED);
$fmt="%.2f\t%.2f\t%d\t%.2f"; 
$fmt.=($showRANGE)?"\t%s":($showSUM)?"\t%d":"\t%.1f"; 
  #below  $fmt="% 4d\t".$fmt if($DOMED); #<< not d sometimes, but %.1f
my $fmt0= $fmt;

if($rotab) {
 $srow= join("\t","Item  ",@hdr); @rows=($srow);
} else {
 print join("\t","Item  ",@hdr)."\n"; 
}
for $i (0..$#sum){ 
  my($ntrim)= ($doTRIM) ? trim($i,$doTRIM) : 0;  
  $lvar=$lvar[$i]||"it$i"; $s=$sum[$i]; $nit=$nit[$i]; ($ave,$var,$sd,$se)=(0) x 4;
  if($nit>0){ $ave=$s/$nit; if($nit>1){ 
    #WRONG: $var=($svar[$i] - $ave*$ave)/($nit-1); 
    if($ENV{dvar}) { $var=0; for my $si ( @{$sval[$i]} ){ $v=$si - $ave; $var+= $v*$v; } $var=$var/($nit-1); } 
    else { $var=($svar[$i] - $nit*$ave*$ave)/($nit-1); }
    $sd=sqrt($var); $se=$sd/sqrt($nit); 
    } 
  }
  @mv=($ave,$se,$nit,$sd); 
  if($DOMED) { my @ss=sort{$a <=> $b} @{$sval[$i]}; $med= $ss[ int($nit/2) ];  unshift(@mv,$med); 
    $fmt=($med =~ /\.\d/)? "%4.2f\t".$fmt0 : "% 4d\t".$fmt0;
    if($showRANGE) { push @mv, "$ss[0]-$ss[-1]"; }
    }
  unless($showRANGE) { push @mv,(($showSUM)?$s:$var); }
  if($rotab) { $srow=sprintf("%6s\t$fmt",$lvar,@mv); push @rows, $srow; }
  else { printf "%6s\t$fmt\n",$lvar,@mv; }
} 

if($rotab) {
  # print "#debug: @rows \n";
  @rc=(); $nc= @rows; $nr= scalar(split"\t",$rows[0]);
  for($i=0; $i<$nc; $i++){ @c= split"\t",$rows[$i]; $nr=@c; for ($j=0; $j<$nr; $j++){ $rc[$j][$i]=$c[$j]; } }
  for($j=0; $j<$nr; $j++) { $rv=$rc[$j]; print join("\t",@$rv); print "\n"; }
  # for ($j=0; $j<$nr; $j++) { for( $i=0; $i<$nc; $i++){ print "\t" if($i>0); print $rc[$j][$i]; } print "\n"; }
} 

}

