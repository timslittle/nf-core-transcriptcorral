#!/usr/bin/perl
# simple_sam2covtab.pl

=item 

[gilbertd@h2 dplx20gnodeq]$ 
input data :
 ls {daphplx_gasm16ml,daphpulex_pa42v2,dplx20maca4pkr_dc}_SRR13333791_b8_bwa.bam
  daphplx_gasm16ml_SRR13333791_b8_bwa.bam  daphpulex_pa42v2_SRR13333791_b8_bwa.bam@  
  dplx20maca4pkr_dc_SRR13333791_b8_bwa.bam

compare to
  daphplx_gasm16ml_chr-SRR13333792uc3c.al90covtab ..

env bam=${pt}_SRR13333791_b8_bwa.bam out=${pt}_chr-SRR13333792bwa1a.al90covtab \
   mina=90 ridpat=SRR13333792 sbatch run_simpsam2covtab.sh


=cut

my $BN=100; my $BN1=$BN - 1;
my $MINAL=$ENV{mina}||1;
my $RIDP=$ENV{ridpat}||""; # SRR13333792 for comparison to dchr-SRR13333792uc3c.blastn
@bam= grep/\.bam/, @ARGV;

for $bam (@bam) {
  open($insam,"samtools view $bam |") or die "samtools view $bam";
  while(<$insam>) {
    my($id,$fl,$cr,$cb,$cx,$cig,@samx)=split;
    if($fl & 0x4){ $n_nomap++; next; } 
    if($RIDP){ next unless($id =~ m/$RIDP/); }
    my($nmi)= (m/NM:i:(\d+)/)?$1:0;    
    my($alen,$lenc,$cend)= addCigar($cr,$cb,$cig);  
    #NOTE: addCigar cend was 0-based, now cb+cend
    # if($nmi>0){ $alen -= $nmi; $n_mismatch += $nmi; } 
    if($alen < $MINAL){  
      $nlo++; 
    } else {
      $nok++; $taln+= $alen; $trdn++;
      my($ib,$ie)= map{ int($_/$BN) } ($cb,$cend);
      for(my $i=$ib; $i <= $ie; $i++) { 
        $ab=$i*$BN; $ab=$cb if($ab<$cb);
        $ae=$i*$BN + $BN1; $ae=$cend if($cend<$ae);
        $abin{$cr}[$i] += 1 + $ae - $ab;
      }
      $amax{$cr}= $ib if($ib>$amax{$cr});
      $aln{$cr}+=$alen; $rdn{$cr}++;
    }
  } close($insam);
}

putCovtab();

sub putCovtab {
  print "#".join("\t",qw(ChrID Pos aCovT))."\n"; 
  for $cr (sort keys %amax) {
    $maxib= $amax{$cr}; $aln=$aln{$cr}; $rdn=$rdn{$cr};
    $maxb= $BN*$maxib;
    print "# chr=$cr, maxb=$maxb, aln=$aln, readn=$rdn\n";
    for($ib=0; $ib <= $maxib; $ib++){ 
      $act=$abin{$cr}[$ib]||0; $cb=$BN*(1+$ib); 
      print join("\t",$cr,$cb,$act)."\n"; 
    } 
  }
}

sub addCigar { 
  my( $cr,$cb,$cigar)=@_;
  my($alen,$lenc,$cend,$nlen,$indel,$softclip)= (0) x 9;
 
  while($cigar =~ m/(\d+)([A-Z])/g) { 
    my($bi,$bt)=($1,$2);
    if($bt eq 'H') { 
      $lenc += $bi; $bi=0; 
    } elsif($bt eq 'S') {
      $softclip += $bi; $lenc += $bi; $bi=0; 
    } elsif($bt eq 'M') {
      # for(my $i=0;$i<$bi;$i++){ $thisdepth[$cend+$i]=1; } #  $aligndepth->[$cend+$i]++; 
      $alen+= $bi; $lenc += $bi;  $cend += $bi;  
    } elsif($bt eq 'N') { 
      $n_intron++;  $cend += $bi; # cend not changed for HISP; 
    } elsif($bt eq 'I') {
      $n_insert += $bi; #UPD7f change for base count eq n_mismatch += nmi
      $lenc += $bi; $bi=0; $indel++;
    } elsif($bt eq 'D') {  # P also?
      $n_delete += $bi; #UPD7f change for base count
      $cend += $bi; $indel++;
    } elsif($bt eq 'P') {  # what is P now?
      $cend += $bi;
    } else {
      # unknown here, what? record?
    }
  }
  
  return($alen,$lenc,$cb+$cend,$softclip);#?
}

