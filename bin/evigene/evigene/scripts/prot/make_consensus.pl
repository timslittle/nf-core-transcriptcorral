#!/usr/bin/env perl
# make_consensus.pl cut from tr2aacds2c.pl
# in: mygenes.nr.cds of fastanrdb;  out: cdsnrseq.consensus

use strict;

use constant AACONS => 1;  # aaconsensus: keep tr with consensus of identical aa across asm types
use constant UPD1912 => 1;
 
my ($USE_AACONS, $aconsensus, $nagree, $agreeset, )= (AACONS, "", 0, {}); # AACONS
my @agrees= grep /agree$|consadd$|aacons$/, @ARGV; # UPD1912
my ($cdsnrseq)= grep { not m/agree$|consadd$|aacons$/ } @ARGV; # all? $ARGV[0];
die "usage: make_consensus.pl mygenes.nr.cds : fastanrdb of mygenes.cds, optional mygenes.agree \n"
  unless($cdsnrseq and (-f $cdsnrseq or $cdsnrseq =~ /stdin|^-$/));
  
# UPD1912 
# my($maind,@agrees)= getFileset('./','agree$|consadd$|aacons$');   # consensus
for my $agf (@agrees) { my($nadd,$agreesetr)= add_consensus_idset($agf, $agreeset); $nagree+= $nadd; }

my @res= make_consensus($cdsnrseq);  

#-------------

sub add_consensus_idset
{
  my($agreefile, $agreelist)=@_; 
  my($ncons,$consfile)=(0,"none");
  $agreelist={} unless(ref $agreelist);

  if($agreefile and -f $agreefile) {
    open(F,$agreefile);  
    while(<F>) { 
      next unless(/^\w\S+\s+\d/);  
      my($id,$agscore,$agsrc,$agids)= split;
      if($agscore > 1) {
        $agreelist->{$id} += $agscore; $ncons++;
        if($agids) { map{ $agreelist->{$_} += $agscore; } split(",",$agids); }
      }
    } close(F);
  }  
  warn("add_consensus_idset n=$ncons from $agreefile\n"); # if $debug; 
  return($ncons,$agreelist);
}

sub make_consensus # ($cdsnrseq, $aaqual)
{
  my($cdsnrseq, $aaqual)=@_;
  my ( %better, %aconsensus ); # %aq,
  my ($nbetter,$nrec,$nrtot,$naddcons)=(0) x 9;
  
  #n our $aaqualh= getAaQual($aaqual); # ret: ($qualh,$sizeh,$naa) : $qualh # may be global hash, read once?
    
  # my($ok,$hin)= openRead($cdsnrseq); # open(F,$cdsnrseq); 
  my($ok,$hin)=(0,0);
  $cdsnrseq||="stdin";
  if($cdsnrseq =~ /stdin|^-$/) { $ok=1; $hin=*STDIN; }
  else { $ok= open(F,$cdsnrseq); $hin=*F; }

  while(<$hin>) { if(/^>/) { 
    s/>//; my @dp=split; $nrtot++;
    if(@dp>1) {
      my $nrflag= 'nrdup';
      
      if($USE_AACONS) {
      my($cscore,$cmain,$ctop,$cset)= consensusof('nrdupcds',@dp); 
      if($nagree){ map{ $cscore += $agreeset->{$_}||0; } @dp;  } #UPD1912
      # * add consensus flag in cdsnrseq header, dont change aaqual better hdr order tho
      if($cscore) {
        $nrflag= "nrdup,agree=$cscore"; #?
        my $cids= join",",@dp; # YES, need these as nofrag rewrites nr.cds headers
        $aconsensus{ $dp[0] }= join"\t",$cscore,$cset,$cids;  
        }
      }
      
      $nrec++; 
      #n $nbetter += nrcheckaaqual($nrflag, $aaqualh, \%better, @dp);
    } elsif($USE_AACONS and $nagree) { #UPD1912
        my $cids= $dp[0];
        my $cscore = $agreeset->{$cids}||0;   # FIXME: cset == prefix of cid .. from consensusof()??
        if($cscore) { $aconsensus{ $dp[0] }= join"\t",$cscore,$cids,$cids; $naddcons++; }
      }
    } 
  } close($hin);

  # AACONS write table trset.consensus for asmrna_dupfilter3c.pl
  if(AACONS and $USE_AACONS) {
    my $maxcon= int($nrtot * 0.40);  # NOT $nrec; dunno how many to allow
    my @cids= sort keys %aconsensus;
    my $ncons= @cids; my $nconstest= $ncons;
    $nbetter= $ncons; # ret this for cons only
    if($naddcons) { #UPD1912 maxcon for addcons/agreeset
      $nconstest= $ncons - $naddcons;
    }
    
    # my $cname= makename($cdsnrseq,".consensus"); #? global? need below for asmrna_dupfilter3c.pl
    my $cname= $cdsnrseq . ".consensus";
    
    if($nconstest > $maxcon) { 
      my $pcons= ($nrtot>0)?int(100*$ncons/$nrtot):0;
      warn("too many cds have consensus to be usable: $ncons of $nrtot ($pcons%); skipping consensus..");  
      #o loggit(1,"too many cds have consensus to be usable: $ncons of $nrtot ($pcons%); skipping consensus..");  

    } else {
      #* put list of all nr ids into trset.consensus? col3, Yes, need these later
      if(-s $cname) { rename($cname, "$cname.old"); }
      my $ok= open(B,">$cname");
      for my $cid (@cids) { print B "$cid\t",$aconsensus{$cid},"\n"; }
      close(B);
      #o loggit(0,"consensus transcripts n=$ncons tabled in $cname"); 
      warn("consensus transcripts n=$ncons tabled in $cname\n"); 
    }
  } 
  
  # if($nbetter>0) { # rewrite $cdsnrseq headers 
  #   my $ok= open(B,">$cdsnrseq.best");
  #   if($ok) { 
  #     ($ok,$hin)= openRead($cdsnrseq); 
  #     while(<$hin>) { if(/^>(\S+)/) { if(my $hdr= $better{$1}) { s/>.*$/>$hdr/; } } print B $_; } 
  #     close(B); close($hin); 
  #     my $cmd="mv $cdsnrseq $cdsnrseq.old; mv $cdsnrseq.best $cdsnrseq";
  #     system($cmd); #? runcmd($cmd);
  #     push @tmpfiles, "$cdsnrseq.old";
  #   }
  # }
   
  return($nbetter,$nrec); # what?
}


sub consensusof { 
  my($nrtest, $dmain, @d2)= @_;
  my($cscore,$cmain,$ctop,$cset)= (0) x 9;
  $cmain= $dmain;
if($USE_AACONS) {

  # assume for now trasm IDs are from evg trformat.pl with common patterns:
  #  subsetprefix | asmmethod | kKMER | [Ll]oc NNNN tIII where locusNNNN alternate tIII are chopped off
  # FIXME: need to use $aconsensus param as options for how to parse ids
  # eg: acons='(trin|velv|idba|soap|pacb)' is source pattern
  # or: acons='(trin|idba|pacb|(?:velv|soap)\w*k\d+)' for kmer setting
  # or: acons='^(\w+)' .. etc
  
  my $userparse= ($aconsensus =~ m/\w/)?1:0;
  my $nokmer= (($aconsensus =~ m/kmer/) or $ENV{'conkmer'})?0:1; #UPD19
  my @dsrc=();
  for my $d ($dmain,@d2) {
    my($s,$sdef)=("",$d);
    if(1) { # default parse
      unless($sdef=~s/[Ll]oc.*$//) { # need checks, if no LocNNN chomp numbers
        $sdef=~s/utrorf//; $sdef=~s/_//g;
        while($sdef=~s/(\d)\D\d+$/$1/) { ; } # AT2G1234.1 > AT2G1234; XM_1234.2 > XM_1234; TrinDN1234_c0_g2_i1 TrinDN1234c1g2t1  > TrinDN; 
        $sdef=~s/\d+$//; 
        }
      $sdef=~s/k\d+$// if($nokmer or $sdef=~/idba/); # fixme for asmbler tags
      $s= $sdef;
    } 
    if($userparse) { ($s)= ($d=~m/($aconsensus)/)?$1:$sdef; }
    push @dsrc, $s if($s=~/\w/);
  }
    
  my %dsrc=();  map{ $dsrc{$_}++ } @dsrc;
  @dsrc= sort{ $dsrc{$b}<=>$dsrc{$a} or $a cmp $b } keys %dsrc; # dont care for sort by count?
  if(@dsrc>1) {
    $cscore= @dsrc; # number of sources, for 1st test
    $ctop= $dsrc[0]; 
    $cset=join",", sort @dsrc; # alpha sort for display
    $cmain= $dmain; # for now
  }
}  
  return($cscore,$cmain,$ctop,$cset);
}


