#!/usr/bin/perl
# evigene/genoasm/chragp2fasta.pl from chragp2contigs.pl

use strict;
use warnings;
use Getopt::Long ;

my $GAPSIZE= 100; ## need option, unknown_gap_size .. use chr.agp N sizes ??
my $FAWIDTH= 60; 
my $MINLEN= 200;
my $IDPREFIX= 'Chr'; 
#o# my $SCAFTAG= 'SG'; # SG{number} is always part of scafid, or see UG ungrouped set?
my $DIEonERR= 1; # debug opt
my($chrdir,$chrasm,$agpdir,$contigs,$unplaced,$placedscafs,$debug,$addhead)= (0) x 10;
my $GAPfromAGP= $ENV{agpgap}||0;
my $CTID_ONLY=$ENV{ctidonly}||0;

my $optok= &GetOptions (
  "contigs|input|chrdir=s" => \$chrdir, # input chr dir OR all.fasta
  "agpdir=s" => \$agpdir, # input agp dir OR all.agp
  "addhead=s" => \$addhead, # add header info
  "output=s" => \$chrasm, # output chrasm fasta, is STDOUT option?
   # ^ also maybe unplaced.agp > unplaced.fasta
   #x "contigs|ctgdir=s" => \$contigs, # output contigs folder
  "IDPREFIX=s"=>\$IDPREFIX,   
  "MINLEN=i"=>\$MINLEN,   
  "GAPSIZE=i"=>\$GAPSIZE,   
  "unplaced!"=>\$unplaced, # agp == unplaced_scaffolds.agp, cid == scafid input.fa, dont use agp locs
  "placedscafs!"=>\$placedscafs, #variant of unplaced_scaffolds, output scafs in agp=chromosomes.agp
  "gapfromagp|gapagp!"=>\$GAPfromAGP, 
  "CTID_ONLY!"=>\$CTID_ONLY, 
  "debug!"=>\$debug,   
  );

# for cfa in seq/mm_ref_GRCm38.p3_*.fa.gz; do { .. }
$chrdir= shift @ARGV unless($chrdir); # or list of fa files?
$agpdir= shift @ARGV unless($agpdir); # or list of agp files?

unless($chrasm) {
   $chrasm= ($placedscafs)? "placed_scafs.fa" : ($unplaced) ? "unplaced_scafs.fa" : "chromosomes.fa"; 
}
  
# if($contigs) { #? use like chrdir, fasta input 
#   warn "# later .. see chragp2contigs.pl\n"; exit -1;
#   # mkdir($contigs) unless(-e $contigs);
# }

my (@agp,@chr,%agp,%agpc,%chrfa,@chrord,%did);
my (%chrgap);
my ($tlen,$tgap,$nout,$ncr,$nskip)= (0) x 9; # globals now
my $outh = undef;

if( -d $chrdir) {
  opendir(D,$chrdir); @chr= sort grep /\.fa/, map{ chomp; $_; } readdir(D);  closedir(D); 
} elsif( -f $chrdir ) {
  @chr=($chrdir); $chrdir=".";
}
if( -d $agpdir) {
  opendir(D,$agpdir); @agp= sort grep /\.agp/, map{ chomp; $_; } readdir(D);  closedir(D); 
} elsif( -f $agpdir ) {
  @agp=($agpdir); $agpdir=".";
}

=item agp eg

  ==> lgmap5u_dmag20sk4abysseal/CHRR_genome.agp <==
  DC10    19414159        19414258        10382   U       100     scaffold        yes     map
  DC10    19414259        19431703        10383   W       dmag20sk4kaby_sc4787036 1       17445   +
  DC10    19431704        19431803        10384   U       100     scaffold        yes     map
  DC10    19431804        19432130        10385   W       dmag20sk4kaby_sc2128334 1       327     +
  
  ==> lgmap5u_dmag20sk4abysseal/CHRR_scaffolds.agp <==
  dmag20sk4kaby_sc4799818 1       6540    1       W       ctg502135       1       6540    +
  dmag20sk4kaby_sc4799818 6541    7203    2       N       663     fragment        yes     
  dmag20sk4kaby_sc4799818 7204    9910    3       W       ctg502136       1       2707    +
  dmag20sk4kaby_sc4799818 9911    10110   4       N       200     fragment        yes     
  dmag20sk4kaby_sc4799818 10111   13378   5       W       ctg502137       1       3268    +
  
  ==> lgmap5u_dmag20sk4abysseal/CHRR_unplaced_scaffolds.agp <==
  dmag20sk4kaby_sc4799702 2813    8805    3       W       ctg501767       1       5993    +
  dmag20sk4kaby_sc4799787 1       2099    1       W       ctg502012       1       2099    +
  dmag20sk4kaby_sc4799787 2100    2299    2       N       200     fragment        yes
  dmag20sk4kaby_sc4799787 2300    5691    3       W       ctg502013       1       3392    +

=cut

 
for my $agpin (@agp) { # only 1 likely ..
  my($pt,$chrin,$ino,$ok,$nid,$id,$fa,$fahdr,$infa,$nout)=(0) x 19;
  ($pt= $agpin)=~s/\.agp.*//; 
  ($chrin)= grep/$pt\.fa/,@chr; # fixme for not same name?
  unless($chrin) { if(@agp==1 and @chr==1) { $chrin=$chr[0]; } }
  
  # warn "# extract contigs: $agpdir/$agpin, $chrdir/$chrin to $contigs/$pt.ctg.fa\n" if($debug);
  warn "# write chr.fasta from $agpdir/$agpin, $chrdir/$chrin to $chrasm\n" if($debug);
  die "ERR:missing $agpdir/$agpin" unless(-f "$agpdir/$agpin");
  die "ERR:missing $chrdir/$chrin" unless(-f "$chrdir/$chrin");
    
  $ino= ($agpin=~/\.gz/) ? "gunzip -c $agpdir/$agpin|" : "$agpdir/$agpin";
  $ok= open(AGP,$ino) or die "open $ino";
  %agp= %agpc=();  @chrord=();
  
  while(<AGP>) {

    if(/^\w\S+\t\d/) { 
      my @v= split; 
      my($cid,$cb,$ce,$ci,$FN,$sid,$sb,$se,$so)=@v;      
      #add?  $chragp{$cid}{$ci}=[@v];
      
      if($FN =~ m/^[NU]/) {
        # dmag20sk6tmaca83v_scf100463	17569	17700	20	N	132	scaffold	yes	paired-ends
        if($GAPfromAGP){ 
          # my $gapr=[$cid,$cb,$ce,$ci,$sid,$sb,$se,$so];
          $chrgap{$cid}{$ci}= $sid; # sid == gap size, save all agpr?
        }
      }
      next if($FN =~ m/^[NU]/ or not($sid=~/\w/ and $se>0));
      
      # check $FN eq "W" ?? maybe keep NU for gap sizes; 
      # unplaced_scaf.agp gap entries need keeping or else write full scaf
      unless($agp{$cid}){ push @chrord, $cid; }
      my $agpr=[$cid,$cb,$ce,$ci,$sid,$sb,$se,$so];
      push @{$agp{$cid}},$agpr; #? drop
      push @{$agpc{$sid}},$agpr;#  sid == input contig/scaf id, use this
    }
  } close(AGP); 

  $ino= ($chrin=~/\.gz/) ? "gunzip -c $chrdir/$chrin |" : "$chrdir/$chrin";
  $ok= open(CHR,$ino) or die "open $ino";
  
  if($chrasm =~ /^stdout|^\-/i) { $outh= *STDOUT; }
  else {
    rename($chrasm, "$chrasm.old") if( -s $chrasm );
    $ok= open(OUT,">$chrasm"); $outh= *OUT;
  }
  
  while(<CHR>) { # not CHR.fasta here, but contigs.fa
    if(/^>/) { 
      $nout += puts($id,$nid,$fa,$fahdr,$outh) if($id and $fa); 
      ($id)= m/^>(\S+)/;
      ($fahdr= $_) =~ s/^>$id//; chomp($fahdr);
      $id=~s/gi.\d+.ref.//; $id=~s/\|$//; 
      $nid=$id;
      # if( my($cn)=m/chromosome (\w+)/) { $id="chr$cn"; } # not here
      $infa=1;  $fa=""; 
    } elsif($infa) { chomp; $fa.=$_; } 
  } close(CHR); 
  
  $nout += puts($id,$nid,$fa,$fahdr,$outh) if($fa); 
  if($unplaced||$placedscafs) { $ncr= $nout; }
  else {
    ($ncr,$tlen,$tgap)= putchrfa($outh);  
  }
  close($outh) if(ref($outh));
  
  my ($mb,$mbg)= map{ int($_/100_000)/10 } ($tlen,$tlen-$tgap);
  $ncr .= " ($nskip too short)" if($nskip);
  warn "# chr.fasta: ncontig=$nout to nchr=$ncr, size_Mb=$mb,$mbg(nogap) in $chrasm\n" if($debug);
}


#--------

sub putfa {
  my($outh, $chrid, $chrfa, $fahdr, $nctg)= @_;
  my $slen= length($chrfa); 
  my($ngap)= $chrfa =~ tr/Nn/Nn/; # count gaps?
  $tlen+= $slen; $tgap+= $ngap; # globals now
  $chrfa =~ s/(.{60})/$1\n/g; # width opt
  $chrfa.="\n" unless($chrfa =~ m/\n$/);
  
  if($fahdr) {
    my($olen)= $fahdr =~ m/(?:len|length|size)\W(\d+)/?$1:0;
    if($olen) { 
      if($slen eq $olen) { $fahdr =~ s/(?:len|length|size)\W(\d+)[;]?//; } 
      else {  $fahdr =~ s/(len|length|size)(\W$olen)/orig$1$2/; }
    }
    $fahdr=~s/\s*$/;/;
  } else { $fahdr=($nctg)?"contigs=$nctg;":""; }
  if($addhead){ $fahdr .= " $addhead;"; }
 
  print $outh ">$chrid length=$slen; gaps=$ngap; $fahdr\n",$chrfa; 
  return 1;
}

sub revc{ my $s=reverse($_[0]); $s=~tr/ACGTacgtRYrySWswKMkmBDHVbdhv/TGCAtgcaYRyrWSwsMKmkVHDBvhdb/; $s;}

sub putchrfa {
  my($outh,$fahdr)= @_; # , $chrord ?
  my $NNN= ('N') x $GAPSIZE;
  $fahdr||= "";
  my $nout=0; # my ($tlen,$tgap,$nout)= (0) x 9; # globals now
  for my $cid (@chrord) {
    # @$agp{$cid} == [$cid,$cb,$ce,$ci,$sid,$sb,$se,$so] list
    # $chrfa{$cid}{$ci}= $s;
    my @it= sort{ $a <=> $b } keys %{ $chrfa{$cid} };
    my $cna= ($cid =~ /^\d+$/) ? sprintf("%02d",$cid) : $cid;
    my $chrid= ($IDPREFIX) ? $IDPREFIX.$cna : $cna;
    my $chrfa="";
    for my $ci (@it) {
      my $s= $chrfa{$cid}{$ci} or next; # err?
      if($GAPfromAGP) {
        if($ci>1) {
        my $gapi= $ci - 1;
        ## FIXME: NOt GAPSIZE, but zero if chrfa == F fill replacing N gap !!!!!
        ## maybe need full chragp here? $agp= $chragp{$cid}{$ci - 1}; if($$agp[4] eq 'N') ..
        my $gapw= $chrgap{$cid}{$gapi} || 0; # $GAPSIZE; # sid == gap size, save all agpr?
        if($gapw > 0) { my $nnn = ('N') x $gapw; $chrfa .= $nnn; }
        }
      } else {
        if($chrfa) { $chrfa .= $NNN; }
      }
      $chrfa .= $s;
    }
    
    my $nctg= @it;
    $nout += putfa($outh, $chrid,$chrfa,$fahdr,$nctg); 
  }
  return($nout,$tlen,$tgap);
}

sub puts { 
  my($id,$nid,$fa,$fahdr,$outh)=@_; my $nout=0;
  # $id,$nid,$fa,$fahdr,$outh
  # NOW nid (or id) is contig id for contig.fa, use agpc
  # use  @{$agpc{$sid}},[$cid,$cb,$ce,$ci,$sid,$sb,$se,$so];#  sid == input contig/scaf id, use this
  
  my $agp;  
  if($placedscafs){
    $agp= $agpc{$id}||$agp{$id}||$agp{$nid}; # placed scaf is 2nd agp id, chr is 1st
  } elsif($unplaced) {
    $agp= $agp{$id}||$agp{$nid}||$agpc{$id}; # unplaced scaf is main/1st agp id
  } else {
    $agp= $agpc{$id}||$agpc{$nid}||$agp{$id};  # id is in scaffold-id hash
  }
  return 0 unless($agp); # likely many contig ids not in chr agp

  my $falen= length($fa);
  if($unplaced or $placedscafs) {
    if($unplaced and $falen < $MINLEN) { $nskip++; return 0; } # add MINLEN option for unplaced
    if($placedscafs) { # add chr:loc to fahdr
      my $ag= $agp->[0];
      my($cid,$cb,$ce,$ci,$sid,$sb,$se,$so)=@$ag;  # test sid == nid
      if($ce and $cb){
        my $cloc="chrloc=$cid:$cb-$ce/this:$sb-$se:$so; ";
        $fahdr= $cloc.$fahdr; 
        }
    }
    my $scafid= ($IDPREFIX =~ /^Chr/i or not $IDPREFIX) ? $id : $IDPREFIX.$id; # use IDPRE or not?
    return putfa($outh, $scafid, $fa, $fahdr); 
  }
  
  for my $ag (@$agp) { 
    my($cid,$cb,$ce,$ci,$sid,$sb,$se,$so)=@$ag;  # test sid == nid
    # use ci or cb as chrfa hash index, both should be ok, ci should be 1..n iterator

    my $serr=""; my $skipit=0;
    (my $sidc= $sid) =~ s/^c\d+_//;
    if($sid eq $id) { } #ok
    elsif($CTID_ONLY) { $skipit=1; next; } # $sid ne $id
    # elsif($sidc eq $id) { } #ok, see BUG2 below, use sb,se = cb,ce
    elsif($sidc ne $id) { $skipit=1; $serr.="sid:$sid ne fa:$id,";
      #? warn "#err.puts: $id,$cid.$ci = $serr\n"; # too many warns
      next;
      } # err
    
    # if($cid eq $id) { } # ok ?
    # elsif($sidc ne $id) { $serr.="sid:$sid ne fa:$id,";} # err
    # ^^ FIXME agp sid adds c1_ to id
    # BUG is for F fill sid != cid, but call puts(cid,cidfa) AND puts(sid,sidfa)
    # .. both have same @agp ?   $agpc{$sid} vs  $agp{$cid}
    
    ## BUG2: when id == cid, cN_cid, need to use cb,ce NOT sb,se for cN contig
    if($id eq $cid and $sid ne $sidc) { ($sb,$se)= ($cb,$ce); } #? drop sid ne sidc

    if(!$skipit and $se > $falen) { $serr.="se:$se>$falen,"; $se=$falen; } # err  
    if(!$skipit and $sb < 1) { $serr.="sb:$sb<1,"; $sb=1; } #err;
    
      # substr outside of string at chragp2fasta.pl line 254, <CHR> line 533. ## WARN
    my $s= ($skipit)? "" : ($sb==1 and $se == $falen) ? $fa : substr($fa,$sb-1,1+$se-$sb);
    unless($s) { $serr.="fa.empty:$sb-1;1+$se-$sb,"; }
    elsif( my $sa= $chrfa{$cid}{$ci} ) { 
      # {$cid}{$ci} should be uniq in agp
      $sa=revc($sa) if($so eq "-"); 
      if($sa ne $s) {
      my $ls=length $s; my $lsa=length $sa;
      $serr.= "duploc $cid.$ci $ls <> $lsa,";  # dup loc? err; *** FIXME bad bug? or replacing same
      $s=""; # dont replace ??
      }
    }
    
    if($serr){ warn "#err.puts: $id,$cid.$ci = $serr\n"; }
    if($s) {  
    $s=revc($s) if($so eq "-"); 
    $chrfa{$cid}{$ci}= $s; # duploc err?
    $did{$sid}++; $nout++;  
    }
    
    #----- for ctg out
    ## not here, need to hash all ctg 
    # my $s=substr($fa,$cb-1,1+$ce-$cb); $s=revc($s) if($so eq "-"); 
    # my $slen=length($s);
    # if($did{$sid}){ for my $sp (qw(a b c d e f g h i j)){ unless($did{$sid.$sp}){ $sid.=$sp; last; } } }
    # $s=~s/(.{60})/$1\n/g; 
    # print OUT ">$sid ofs=$sb-$se:$so; chrloc=$id:$cb-$ce; chrid=$nid; ci=$ci; length=$slen\n";
    # print OUT "$s\n"; $did{$sid}++; $nout++;
  } 
  
  return $nout;
}

=item results

pt=lgmap5u_dmag14bgi2vtop5k; 
./chragp2fasta.pl -debug -idpre dmagbgi5u_ -agp $pt/CHRR_genome.agp  -contigs  genome/dmag14bgi2vtop5k.fa.gz -out $pt/CHRR_chrs.fasta
# write chr.fasta from ./lgmap5u_dmag14bgi2vtop5k/CHRR_genome.agp, ./genome/dmag14bgi2vtop5k.fa.gz to lgmap5u_dmag14bgi2vtop5k/CHRR_chrs.fasta
# chr.fasta: ncontig=2796 to nchr=10, size_ Mb=168.4,153.3(nogap) in lgmap5u_dmag14bgi2vtop5k/CHRR_chrs.fasta

pt=dmag20skspad9k; 
./chragp2fasta.pl -debug -idpre dmagsk9spad_ -agp lgmap5u_$pt/CHRR_genome.agp  -contigs  genome/$pt.fa.gz -out lgmap5u_$pt/CHRR_chrs.fasta
# write chr.fasta from ./lgmap5u_dmag20skspad9k/CHRR_genome.agp, ./genome/dmag20skspad9k.fa.gz to lgmap5u_dmag20skspad9k/CHRR_chrs.fasta
# chr.fasta: ncontig=57060 to nchr=10, size_ Mb=183,147.7(nogap) in lgmap5u_dmag20skspad9k/CHRR_chrs.fasta

pt=dmag20sk4abysseal
./chragp2fasta.pl -debug -idpre dmagsk4abys_ -agp lgmap5u_$pt/CHRR_genome.agp  -contigs  genome/$pt.fa.gz -out lgmap5u_$pt/CHRR_chrs.fasta
# write chr.fasta from ./lgmap5u_dmag20sk4abysseal/CHRR_genome.agp, ./genome/dmag20sk4abysseal.fa.gz to lgmap5u_dmag20sk4abysseal/CHRR_chrs.fasta
# chr.fasta: ncontig=60042 to nchr=10, size_ Mb=203.3,181.9(nogap) in lgmap5u_dmag20sk4abysseal/CHRR_chrs.fasta

pt=dmag20ugspad4d_nobac
./chragp2fasta.pl -debug -idpre dmagug4spad_ -agp lgmap5u_$pt/CHRR_genome.agp  -contigs  genome/$pt.fa.gz -out lgmap5u_$pt/CHRR_chrs.fasta
# write chr.fasta from ./lgmap5u_dmag20ugspad4d_nobac/CHRR_genome.agp, ./genome/dmag20ugspad4d_nobac.fa.gz to lgmap5u_dmag20ugspad4d_nobac/CHRR_chrs.fasta
# chr.fasta: ncontig=14030 to nchr=10, size_ Mb=106.3,95.5(nogap) in lgmap5u_dmag20ugspad4d_nobac/CHRR_chrs.fasta

#... unplaced : ** need bacteria/foreign dna screen

pt=dmag20skspad9k; idpt=dmagsk9spad_;  

./chragp2fasta.pl -debug -idpre 0  -unplaced -agp lgmap5u_$pt/CHRR_unplaced_scaffolds.agp  -contigs  genome/$pt.fa.gz -out lgmap5u_$pt/${idpt}unscaf.fa
# write chr.fasta from ./lgmap5u_dmag20skspad9k/CHRR_unplaced_scaffolds.agp, ./genome/dmag20skspad9k.fa.gz to lgmap5u_dmag20skspad9k/dmagsk9spad_unscaf.fa
# chr.fasta: ncontig=7280 to nchr=7280, size_Mb=6,5.9(nogap) in lgmap5u_dmag20skspad9k/dmagsk9spad_unscaf.fa

head lgmap5u_dmag20skspad9k/dmagsk9spad_unscaf.fa
>dmag20spad9k_sc002102 length=28963; gaps=10;  ; cov=25.956131;
GGCCCGCAGCCGAAGGCGAGGACACGGGCCACCCCGGAGGGCAAAGCGAACCCTTGACCC

./chragp2fasta.pl -debug -idpre $idpt -agp lgmap5u_$pt/CHRR_genome.agp  -contigs  genome/$pt.fa.gz -out lgmap5u_$pt/${idpt}chrs.fa
# write chr.fasta from ./lgmap5u_dmag20skspad9k/CHRR_genome.agp, ./genome/dmag20skspad9k.fa.gz to lgmap5u_dmag20skspad9k/dmagsk9spad_chr.fa
# chr.fasta: ncontig=57060 to nchr=10, size_Mb=183,147.7(nogap) in lgmap5u_dmag20skspad9k/dmagsk9spad_chr.fa

head lgmap5u_dmag20skspad9k/dmagsk9spad_chr.fa
>dmagsk9spad_DC01 length=22972158; gaps=4156396; contigs=7029;
CGTCGACAGTTAATTGCAACATGGACGGAATGGCCACTGTTCCATGATTAGCAGAACCTA
CAGGACCTCCGGATGCGACCAGCTCATCATTCTTGTGAATGTTGACGACTAAATCGGTGC

pt=dmag14bgi2vtop5k; idpt=dmagbgi5u_;  
./chragp2fasta.pl -debug -idpre $idpt  -unplaced -agp lgmap5u_$pt/CHRR_unplaced_scaffolds.agp  -contigs  genome/$pt.fa.gz -out lgmap5u_$pt/${idpt}unscaf.fa
# write chr.fasta from ./lgmap5u_dmag14bgi2vtop5k/CHRR_unplaced_scaffolds.agp, ./genome/dmag14bgi2vtop5k.fa.gz to lgmap5u_dmag14bgi2vtop5k/dmagbgi5u_unscaf.fa
# chr.fasta: ncontig=2168 to nchr=2168, size_Mb=12.2,7.2(nogap) in lgmap5u_dmag14bgi2vtop5k/dmagbgi5u_unscaf.fa

pt=dmag20sk4abysseal; idpt=dmagsk4abys_;  
./chragp2fasta.pl -debug -idpre $idpt  -unplaced -agp lgmap5u_$pt/CHRR_unplaced_scaffolds.agp  -contigs  genome/$pt.fa.gz -out lgmap5u_$pt/${idpt}unscaf.fa
# write chr.fasta from ./lgmap5u_dmag20sk4abysseal/CHRR_unplaced_scaffolds.agp, ./genome/dmag20sk4abysseal.fa.gz to lgmap5u_dmag20sk4abysseal/dmagsk4abys_unscaf.fa
# chr.fasta: ncontig=435287 to nchr=435287 (7 too short), size_Mb=153.2,152.4(nogap) in lgmap5u_dmag20sk4abysseal/dmagsk4abys_unscaf.fa
#  nchr=435287 (7 too short) << too many at MIN=200 ? maybe cut at MINW=500

./chragp2fasta.pl -minlen=500 .. large diff in short-scaf count and total Mb 
# chr.fasta: ncontig=43336 to nchr=43336 (391958 too short), size_Mb=38.7,37.8(nogap) in lgmap5u_dmag20sk4abysseal/dmagsk4abys_unscaf.fa

pt=dmag20ugspad4d_nobac; idpt=dmagug4spad_;  
./chragp2fasta.pl -minlen=500  -debug -idpre 0  -unplaced -agp lgmap5u_$pt/CHRR_unplaced_scaffolds.agp  -contigs  genome/$pt.fa.gz -out lgmap5u_$pt/${idpt}unscaf.fa 
# write chr.fasta from ./lgmap5u_dmag20ugspad4d_nobac/CHRR_unplaced_scaffolds.agp, ./genome/dmag20ugspad4d_nobac.fa.gz to lgmap5u_dmag20ugspad4d_nobac/dmagug4spad_unscaf.fa
# chr.fasta: ncontig=24517 to nchr=24517, size_Mb=70.3,68.6(nogap) in lgmap5u_dmag20ugspad4d_nobac/dmagug4spad_unscaf.fa


=cut

