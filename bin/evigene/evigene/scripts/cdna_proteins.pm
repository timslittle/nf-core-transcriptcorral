# cdna_proteins.pm

# package cdna_proteins;
package main;

use strict;
use warnings;

#? require Exporter;
# our @ISA = qw (Exporter);
# our @EXPORT = qw (translate_sequence get_protein reverse_complement);

use vars qw ($DEBUG $MINAA $MINEXON $MINGOOD $USEGOODLEN 
      $AA_cdna_GT_genome $ORF_FULLvPART $KEEPSAMECDS $NoStopCodon $InnerStopToX
      $pCDSbad $pCDSpoor $MINUTRORF  $USESelenocysteine $TRIMCDSENDGAP $LIKELY_HEXCODE
      );

use constant cdna_proteins_VERSION  => '20200326'; #UPD20UORF; '20191115'; 
use constant UPD1908 => 1;
use constant UPD20UORF => 1; # 2020mar26: update utrorf annots thruout evigene
use constant TESTPOLYA => 0; # UPD20ja; uncommon, 27/48000 in arath16ap, ~3/1000 in ixotick trasm ; not useful

# UPD1908: keep $ORF_FULLvPART = 0.85 default, tested with UPD1908, gives good balance now 
#   for aacomplete for ref genes (eg arabid) vs long aapartial for incomplete trasms (eg. 1000Plants)
# 20190829 updates for qual report, cdshexcode, 
# .. maybe lower default ORF_FULLvPART from 0.85 to 0.50? or lower? to match "standard" practice reference gene data sets
# '20141231'; # Selc annots added
# '20140317';  # '20131124' 20130818 0228'20120721'; 

our $DEBUG=0;
our $MINAA= 30;  # used 
our $MINEXON= 60; # for cut
our $MINGOOD= 0.75; # filter out prots w/ fewer good aminos
# our $USE_CDSEXONS = 0;
our $USEGOODLEN=1;
our $AA_cdna_GT_genome= 1.50;  # for prot(cdna) > prot(genome) test; option?
# FIXME: adjust when to take partial vs complete, eg partial5 often is a few aa longer than complete M
our $ORF_FULLvPART = $ENV{ORF_FULLvPART} || 0.85; # full == 0.0001 for always complete/full
our $KEEPSAMECDS= 0; # prefer keep same CDS exons but can extend/shorten protein bounds
our $NoStopCodon=0; 
our $InnerStopToX=0;  # 3 values?  0=leave *, 1= tr/*/X/, -1 or 2=skip orf
    ## need global params for utrbad/poor
our $pCDSbad = $ENV{pcdsbad}  ||30; # adjust DOWN *?
our $pCDSpoor= $ENV{pcdspoor} ||60;
our $MINUTR= $ENV{minutr}||300; # ~ 300b "fixed" average utr sizes, maybe too low, 
#our $BAD_GAPS= $ENV{aagapmax} || $ENV{BAD_GAPS} || 15;  # % gaps in AA

our $MINUTRORF= $ENV{minutrorf}||300; # was 300; #?
our $TRIMCDSENDGAP= 1; #201503 add; some bugs here
our $LIKELY_HEXCODE= UPD1908; #20190829 add; likely_coding_hexamer() stats

# our @stop_codons = qw(TAA TAG TGA); # allow changes, esp TGA => SelC
# parts from PASA/PasaLib/Nuc_translater.pm 
use vars qw ($currentCode %codon_table @stop_codons %amino2codon_table %llHexcode);
use constant backtrans_AmbiguousNucs => 0;  # default off
our (%amino2codon_ambig,%amino2codon_fixed); # see backtranslate_protein & BEGIN

# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }


sub bestorf_test
{
  my($bestorf,$nextorf) = @_;
  my ($bestprot)= orfParts($bestorf);
  my ($nextprot)= orfParts($nextorf);
  my $MinCDS= 3*$MINAA;

       # FIXME bestorf_test: option too high -ratiocdna 1.25; at least should replace cdna XXXX gaps with cdnain perfect
       # test for near-sameprot but for XXX gaps? 
  
  my $oki= ( $nextprot =~ /\w/ 
    && ($nextorf->{goodlen} >= $MinCDS)  # dang, goodlen,leng are cdslen not aalen
    && ($nextorf->{goodlen}/$nextorf->{length} >= $MINGOOD)) ? 1 : 0;
    
  if($oki) {
    if( $bestprot ) {
      ## problem here for rev=0,1 missing huge partial rev for short fwd .. but get huge if just rev
      my $arat = $nextorf->{goodlen} /  $bestorf->{goodlen}; # ok,> MINAA
      my $adiff= $nextorf->{goodlen} -  $bestorf->{goodlen};  

      if( $arat > 1.0 and  $arat < $AA_cdna_GT_genome and $bestprot =~ m/XX/) { ##  $bestorf->{goodlen}/$bestorf->{length} < 0.99
        (my $bp= $bestprot) =~ s/X/./g; 
        # if(length($bp) > length($nextprot) { } # chomp some?
        return(1, $nextprot, $nextorf) if($nextprot =~ m/$bp/);
      }
      
      if( $arat > $AA_cdna_GT_genome
        or ( $adiff >  0 and ($nextorf->{complete}==3 or $bestorf->{complete} <3) ) 
        or ( $adiff >= 0 and $nextorf->{complete}==3 and $bestorf->{complete} <3 ) )
       { 
        return(1, $nextprot, $nextorf);
       }
    } else {
      return(1, $nextprot, $nextorf);
    }
  }
  return( 0, $bestprot, $bestorf);
}

sub utrorf_test
{
  my($bestorf, $allorfs, $cdnasize) = @_;

  my ($utrorf,$utrosize)= getUtrOrf($bestorf, $allorfs, $cdnasize);  
  if ($utrorf) {
    my ($orfok) = bestorf_test(undef,$utrorf);
    return($utrorf,$utrosize) if($orfok);  
  }
  return(0,0);  
}

sub revorf_report
{
  my($bestorf, $allorfs, $cdnasize, $issorted) = @_;

  # UPD1908: add qual stat of bestorf vs allorfs, ie how much longer is best?
  # turn this into allorf_report ?? revorf special case
  my $ALLORFINFO= UPD1908; # option?
  my($raalen,$rinfo)=(0,"");
	my($rall,$secbestorf,$longerorf)=("",0,0);
	
	if($LIKELY_HEXCODE) {
	  my($codeval, $codetype)= likely_coding_hexamer("any", $bestorf);
	  $rinfo .=";" if ($rinfo);  $rinfo .= "codepot=$codetype/$codeval";
	}
	
  my ($revorf,$revosize)= getRevOrf($bestorf, $allorfs, $cdnasize, $issorted);  
	if($ALLORFINFO) {
	  ($rall,$secbestorf,$longerorf)= reportAllOrf($bestorf, $allorfs, $cdnasize, $issorted);
	  if($secbestorf eq $revorf) { $rall =~ s/aa2nd=/aarev=/; }
	  $rinfo .=";" if ($rinfo);  $rinfo .= $rall;
	}
  
  if ($revorf and $revorf ne $secbestorf and $revorf ne $longerorf) {
    my($rrept, $revaalen)= reportOrfQual("aarev", $bestorf, $revorf, $cdnasize);
    $raalen= $revaalen; ## $revosize; # this is cds-size; want **aasize**
    $rinfo .=";" if ($rinfo);  $rinfo .= $rrept;
	}
  return($raalen,$rinfo,$revorf);  
}

sub reportOrfQual {
  my($tag,$bestorf,$revorf,$cdnasize)= @_;
  return("$tag=0%,none") unless($revorf and ref($revorf));
 
  my($aalen,$pcds,$compl,$orflen,$fahead) = proteindoc($revorf, $cdnasize); 
   ## note: fahead= "aalen=99,90%,complete; clen=350; strand=-; offs=9-309;";
   ## aarev=55%,aa99,90%,complete,s-,o9-309;  << use this?  or aarev=99,...,9-309:-; 
  my $bestsize= $bestorf->{goodlen} || 1; # orflen == cds-length
  my $revsize= $revorf->{goodlen} || 1; # orflen == cds-length
  my $plong= int(0.5 + 100*$revsize / $bestsize);
  
  $fahead =~ s/aalen=(\d+),\d+%/${1}aa/;  # =~ s/aalen=/aa/; 
  $fahead =~ s/clen=\d+;//;  $fahead =~ s/strand=/s/; $fahead =~ s/offs=/o/; 
  $fahead =~ s/=/:/g; $fahead =~ s/ //g; $fahead =~ s/;\s*$//;  $fahead =~ s/;/,/g; 
  if(UPD20UORF) {
    $fahead =~ s/uorfcut=/uorfc=/; $fahead =~ s/uorfoff=/uorfo=/; #??
  }
	if($LIKELY_HEXCODE) {
	  my($cv, $ctyp)= likely_coding_hexamer("any", $revorf);
	  $fahead .= ",cp$cv"; #NO .substr($ctyp,0,1);
	}
  return ("$tag=$plong%,$fahead", $aalen);  
}

sub reportAllOrf # report2ndBest 
{
  my ($bestorf, $orfs, $cdnalen, $issorted)= @_;
  my $MinCDS= 3*$MINAA; $issorted||=0;
  my ($secbest, $secsize, $nall, $ngood, $gotbest, $longerorf, $longsize)=(0) x 9;
  return("") unless(ref($bestorf) and $bestorf->{goodlen} >0);
  my $bestsize= $bestorf->{goodlen};
  
	foreach my $orf (@$orfs) {  # expect long sorted, or check all? 
		my $ogood= $orf->{goodlen};
	  $nall++;  $ngood++ if($ogood>=$MinCDS);
	  if($orf eq $bestorf) { $gotbest=1; next; }
		# my ($oorient,$ogood,$osize)= ($orf->{orient},$orf->{goodlen},$orf->{length});
		unless($gotbest or $longerorf) { $longerorf=$orf; $longsize= $ogood; }
	  if($gotbest and not $secbest) { $secbest= $orf; $secsize= $ogood; }
   }  
  # my $plong= int(0.5 + 100*$secsize / $bestsize); # use same as for revorf
  # my $secrept= "aa2nd=$plong%,${secsize}aa2,${ngood}ng/${nall}n";
  my($secrept)= reportOrfQual("aa2nd", $bestorf, $secbest, $cdnalen); # secbest == 0 ok
  $secrept.=",ng${ngood}"; # $secrept.=",${ngood}ng/${nall}n";
  if($longerorf) {
    my($lorept)= reportOrfQual("aalonger", $bestorf, $longerorf, $cdnalen);
    $secrept.= ";$lorept"; 
  }
  return($secrept, $secbest, $longerorf);    
}

sub proteindoc
{
  my($orf, $cdnalen, $cdnarev, $forFASTA) = @_;
  my($aalen, $compl, $pcds, $istop, $Selc) = (0) x 10;
  my $uorfcut=""; #UPD20UORF
  
  my( $orfprot, $prostart, $proend, $orflen, $orient) = orfParts($orf); # [qw(protein start stop length orient cdnalen complete)]
  $cdnarev ||= $orient; # fill in blank; use always? shouldnt need to pass cdnarev as param.
    # Ooops, orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
  if($orfprot) {
    $pcds  = ($cdnalen>0 && $orflen>0) ? int(100*$orflen/$cdnalen) : 0;
    my $urev= ($prostart>$proend)?1:0;
    my $u1len= ($urev) ? $cdnalen - $prostart : $prostart - 1; 
    my $u2len= ($urev) ? $proend - 1 : $cdnalen - $proend;

    if(UPD20UORF) {
      if(my $ucut=  $orf->{'uorfcut'}) { # 'uorfcut=123-456' should be key annot to others?
        my $uo= $orf->{'uorfoff'}||"0";
        $uorfcut=" uorfcut=$ucut;uorfoff=$uo;"; # $ucut/$cdnalen ??
        # modify pcds, u1len, u2len, ? others for uorf offs,span
        my($ucb,$uce)= split /\W/,$ucut;
        my($uob,$uoe,$uoo)= split /\W/,$uo; ($uoo)= substr($uo,-1,1); # +/-
        unless($uce > 0 and $uoe > 0) { 
          # bad vals??
        } else {
        my $uclen= ($ucb > $uce)? 1 + $ucb - $uce : 1 + $uce - $ucb;
        $pcds  = ($uclen>0 && $orflen>0) ? int(100*$orflen/$uclen) : 0;
        $urev= ($uob > $uoe)?1:0;
        $u1len= ($urev) ? $uclen - $uob : $uob - 1; 
        $u2len= ($urev) ? $uoe - 1 : $uclen - $uoe;
        }
      }
    }
      
    $aalen= length($orfprot); # $bestorf->{length}; # this is cds-len
    if(substr($orfprot,-1,1) eq '*') { $aalen--; if($NoStopCodon) { $orfprot =~ s/\*$//; } }
    $istop= $orf->{innerstop} || 0; # add 201403:
    $compl= $orf->{complete}; #? trust? or check ^M..*$
    $compl= ($compl==3)?"complete":($compl==2)?"partial5":($compl==1)?"partial3":"partial";
    
    ##? not bad if partial? if u1len or u2len == 0
    ## need global params for utrbad/poor
    if($cdnalen - $orflen <= $MINUTR) { } # ignore pcds if utr small
    elsif($pcds < $pCDSbad or $u1len > $orflen or $u2len > $orflen) { $compl.="-utrbad"; }  
    elsif($pcds < $pCDSpoor) { $compl.="-utrpoor";  } #?? maybe change to use prostart OR protend3 > 35%? 40% ?
    ##? add istop flag to compl ??
if(TESTPOLYA) {  # TEST UPD20ja .. uncommon, 27/48000 in arath16ap, ~3/1000 in ixotick trasm           
    if( my $pa= $orf->{polya} ) { $compl.=",polya$pa"; } 
}
    
    # 2014.12 add Selc flags 
    # Funhe2Exx11m009903t6 aalen=556,80%,complete,selcstop; Selcstop=index; .. Name=Selenoprotein N 
    if($USESelenocysteine and index($orfprot,'u')>0) { 
      $compl.=",selcstop"; 
      my $isel= $orf->{Selc}||1; $Selc="Selcstop=$isel"; # want mRNA/CDS index * 1-origin,may be list: 123,456,888
      ## NCBI needs mRNA position, as range (3b codon), .. do any have >1 u, yes
      ## CDS  /transl_except=(pos:1002..1004,aa:Sec);
      # >Funhe2Exx11m009903t6 aalen=556,80%,complete; Selcstop=1; clen=2067; strand=+; offs=128-1798; pubid=Funhe2EKm010847t1; oid=Funhe2EKm010847t1,Funhe2Emap3m010942t3; organism=Fundulus_heteroclitus; type=protein; isoform=1; Name=Selenoprotein N (100%P); genegroup=FISH11G_G7797; Dbxref=TrEMBL:UniRef50_Q9NZV5,TrEMBL:SELN_HUMAN,; tblerr=E:InternalStop,E:MisMatchAA,
      # GALDDQSC>>u<<GSGRTLRETVLESSPVLALLNQSFVSSWSLVRELENMQADEENPALSEKAR
      ## 2+ u cases: Funhe2Exx11m049881t1/Funhe2EKm029327t1 SelM, aalen=135,99%,partial; Selcstop=1; tho 2nd 'uu' may be real stop.
      ## Funhe2Exx11m129535t7/Funhe2EKm029571t1 aalen=321,92%,complete; Selcstop=1; Name=SelP; 
    }
    
    if($forFASTA) { $orfprot =~ s/(.{60})/$1\n/g; } #? only for fasta output
  } else {
    $orfprot="X"; # ? make dummy orfprot?
  }
  
  ## fixme: add goodlen or gaps count to aadoc    
  # my $aagap= $orf->{length} - $orf->{goodlen};
  # my $aagap= $orfprot =~ tr/X/X/;
  
  my $fahead= "aalen=$aalen,$pcds%,$compl; clen=$cdnalen; strand=$cdnarev; offs=$prostart-$proend;$uorfcut";
  $fahead .= " $Selc;" if($Selc);
  $fahead .= " innerstop=$istop;" if($istop>0);
  
  # UPD1908: add code potential attr, in proteindoc(), may be in revinfo..
  if(0 and $LIKELY_HEXCODE and $orfprot ne "X") {
    # stick codeval into $orf->{codepot} so dont need recalc?
    my($codeval, $codetype)= likely_coding_hexamer("any", $orf);
    $fahead .= " codepot=$codetype/$codeval;";
  }
          
  
  if(my $orflags= $orf->{flags}) { $fahead .= " orflags=$orflags;"; }
  return($aalen,$pcds,$compl,$orflen,$fahead,$orfprot); # add codepot as sep item?
}

sub proteinqual # short version proteindoc
{
  my($orf, $cdnalen) = @_;
  $cdnalen||=0;
  my( $orfprot, $prostart, $proend, $orflen, $orient, $ocdnalen, $complete, $istop) 
    = orfParts($orf, [qw(protein start stop length orient cdnalen complete innerstop)]);
 
  $cdnalen= $ocdnalen unless($cdnalen > $ocdnalen);

  my $pcds = ($cdnalen>0 && $orflen>0) ? int(100*$orflen/$cdnalen) : 0;
  my $urev= ($prostart>$proend)?1:0; # or $orient eq '-' ??
  my $u1len= ($urev) ? $cdnalen - $prostart : $prostart - 1; 
  my $u2len= ($urev) ? $proend - 1 : $cdnalen - $proend;
    # $istop= $orf->{innerstop} || 0; # add 201403:
  my $compl= ($complete==3)?"complete":($complete==2)?"partial5":($complete==1)?"partial3":"partial";
  if($cdnalen - $orflen <= $MINUTR) { } # ignore pcds if utr small
  elsif($pcds < $pCDSbad or $u1len > $orflen or $u2len > $orflen) { $compl.="-utrbad"; }  
  elsif($pcds < $pCDSpoor) { $compl.="-utrpoor";  } #?? maybe change to use prostart OR protend3 > 35%? 40% ?
  my $cdsoff="$prostart-$proend"; $cdsoff.=":$orient" if($orient);
  
  my $Selc=0;
  my $aalen= int($orflen/3); # length($orfprot); 
  if($orfprot) {
    $aalen= length($orfprot); # $bestorf->{length}; # this is cds-len
    if(substr($orfprot,-1,1) eq '*') { $aalen--; if($NoStopCodon) { $orfprot =~ s/\*$//; } }
    # Funhe2Exx11m009903t6 aalen=556,80%,complete,selcstop; Selcstop=index; .. Name=Selenoprotein N 
    if($USESelenocysteine and index($orfprot,'u')>0) { 
      $compl.=",selcstop"; 
      my $isel= $orf->{Selc}||1; $Selc="Selcstop=$isel"; # want mRNA/CDS index * 1-origin,may be list: 123,456,888
      ## NCBI needs mRNA position, as range (3b codon), .. do any have >1 u, yes
      ## CDS  /transl_except=(pos:1002..1004,aa:Sec);
    }
  }

  $compl.=",innerstop$istop" if($istop>0);
  if(wantarray) {
    return($aalen,$pcds,$compl,$orflen,$cdnalen,$cdsoff,$Selc);# ,$orfprot
  } else {
    my $fahead= "aalen=$aalen,$pcds%,$compl;cxlen=$orflen/$cdnalen;cdsoff=$cdsoff;";
    $fahead .= "$Selc;" if($Selc);
    return $fahead;
  }
}


sub getBestProt2
{
  my($ptype, $cdna, $exongff, $oldStart_b, $oldStart_e)= @_;
  $oldStart_b||=0; $oldStart_e||=0;
  
  # fix this to return longest full prot, and longest partial (if longer) .. test which is best.
  # FIX2: add test intron overlaps : retained = intron inside exon; err = intron rev at splice/inside
  # FIXME3: option to mark/return 2ndary orf(s) in aberrant long-utr transcripts, that dont overlap 1st orf
    
  my $strands= ($ptype =~ /(both|fwd|for|rev)/i)? $1 : "fwd";
  my ($longorf,$longfull,$orfs) = getAllOrfs($cdna, $strands, $ptype);  
  
  ## 201402 option: ($longorf)= getOneOrf($cstart,$cend,$seq,$strand); # return (lor,for,@orfs)?
  #  where strands=null/fwd/rev/both, and strands=rev means do revcomp on cdna
  # Ooops, orf->{start,stop} are reversed for strand=rev; start>stop
  
  if(ref($longorf)) {
    my $lookmore= ($ptype =~ /long/)?0:1; 
    my($lostart, $loend, $cdnalen, $locomplete, $loflags) = orfParts($longorf, [qw( start stop cdnalen complete flags)]);
    $lookmore=0 if($loflags && $loflags=~/bestspan/);# 201402 orthobest; fix2: return not bestspan..
    
    if($KEEPSAMECDS and $oldStart_b > 0) { 
      # may not be right yet; try this?  { $_->{start} < $oldStart_e && $_->{stop} > $oldStart_b }
      my($samestartorf);
      if($ORF_FULLvPART <= 0.8) {
        ## oldStart_e is end point of orf.start NOT orf.stop
      #BUG:($samestartorf) = grep { $_->{complete} == 3 and $_->{start} >= $oldStart_b and  $_->{stop} <= $oldStart_e } @$orfs;
      ($samestartorf) = grep { $_->{complete} == 3 and $_->{start} >= $oldStart_b and  $_->{start} <= $oldStart_e } @$orfs;
      } else {
      ($samestartorf) = grep { $_->{start} >= $oldStart_b and  $_->{start} <= $oldStart_e } @$orfs;
      #BUG:($samestartorf) = grep { $_->{start} >= $oldStart_b and  $_->{stop} <= $oldStart_e } @$orfs;
     }
      if(ref $samestartorf) { $longorf= $samestartorf; $lookmore=0; } #NOT# else { return (undef); } # not found == no change, here
    } 
    
    $lookmore=0 if($locomplete >= 3 or not ref($longfull) or $longorf eq $longfull);
    if($lookmore) { #  and $locomplete < 3 and ref($longfull) 
      # my $keylen=($USEGOODLEN)?"goodlen":"length";
      my $lsize= $longorf->{'goodlen'}; #was $keylen
      my $fsize= $longfull->{'goodlen'}; #was $keylen

if(UPD1908) { ## fix of fix of fix...
        # UPD1908: adjust here for partial at end of transcript .. keep partial if $ORF_FULLvPART < 0.80? and end-part < ~4 bp
        # problem w/ lots of 5'coding precedes known AAstart, up to start of mRNA == uORF/oORF biology
      use constant PARTatEND => 4;  # 6 is too big, 5?
      my $useful=0; # $longorf= $longfull if($useful); 
      ($lostart,$loend)= ($loend,$lostart) if($lostart > $loend); # not?
      if($lostart < PARTatEND or ($cdnalen - $loend) < PARTatEND) {
        if($fsize >= $ORF_FULLvPART * $lsize) {
          unless(ref($exongff)) {
            $useful=1;            
          } else {
            my ($cdslong, $attrL, $pcodeL, $maxutrL)= getCDSgff2( $exongff, $longorf);
            my ($cdsfull, $attrF, $pcodeF, $maxutrF)= getCDSgff2( $exongff, $longfull);
            $useful=1 if( $maxutrF < 3 or ($pcodeF >= $ORF_FULLvPART * $pcodeL));
          }
        }
      } else { # no partials unless they abut ends of tr
        $useful=1 unless($fsize < 0.20 * $lsize); # not for huge diff ?
      }
    $longorf= $longfull if($useful);
} else {
    # FIXME: adjust when to take partial vs complete, eg partial5 often is a few aa longer than complete M
    if($fsize >= $ORF_FULLvPART * $lsize) {
      if(ref($exongff)) {
        my ($cdslong, $attrL, $pcodeL, $maxutrL)= getCDSgff2( $exongff, $longorf);
        my ($cdsfull, $attrF, $pcodeF, $maxutrF)= getCDSgff2( $exongff, $longfull);
        $longorf= $longfull if( $maxutrF < 3 or ($pcodeF >= $ORF_FULLvPART * $pcodeL));
      } else {
        $longorf= $longfull;
      }
    }
}
    }
  }

  return($longorf, $orfs); # return allorfs now
  # old# return($orfprot, $prostart5, $proend3, $longorf, $orfs); ##, $utrorf 
}

sub getBestProt # old 
{
  # my($ptype, $cdna, $exongff, $oldStart_b, $oldStart_e)= @_;
  my ($longorf, $orfs)= getBestProt2( @_ );
  my ($orfprot,$prostart5,$proend3)= orfParts($longorf);
    # my $getutrorf=1; # always test for utr protein? moved out of here
    # my ($utrorf,$utrosize)= getUtrOrf($longorf, $orfs, length($cdna)); # ($getutrorf) ? xxx() : (0,0);
    # ^^ remove from here, add to utrorf_test
  return($orfprot, $prostart5, $proend3, $longorf); ##, $utrorf 
}

sub getUtrOrf
{
  my ($longorf, $orfs, $cdnalen)= @_;
  my ($utrorf,$utrosize)=(0,0);
  my $lsize= $longorf->{length};
  # my $lgood= $longorf->{goodlen};
  return($utrorf,$utrosize) unless(($cdnalen - $lsize >= $MINUTRORF));  # test even if lsize/cdna > 60% ?  
  # my $dotest= (($cdnalen - $lsize > $MINUTRORF) || ((100*$lsize/$cdnalen) < $pCDSpoor));
  #  ^^ this is bad test; some large tr with pCDS > pCDSpoor are hiding other orf, usually in 3+utr exons on 1 end
  if(1) { 
    # Ooops, orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
    my($lb,$le)= ($longorf->{start},$longorf->{stop});  ($lb,$le)=($le,$lb) if($lb>$le);
    foreach my $orf (@$orfs) {  # orfs can be rev of longorf, start/stop are fwd tho
      my($ob,$oe,$ogood,$osize)= ($orf->{start},$orf->{stop},$orf->{goodlen},$orf->{length},);
      ($ob,$oe)=($oe,$ob) if($ob>$oe);
      if(($ob > $le or $oe < $lb) and ($ogood>=$MINUTRORF) and $ogood>$utrosize) { # size or $osize > 0.5*$lsize ??
        $utrorf= $orf; $utrosize= $ogood;        
        }
      }
    }
  
  # 201402 fixme: do splitutrorf mrnaseq here, see also introncut mrnaseq ..
  # longorf and utrorf should both have ->{mrnaseq} with proper portions cut..
  # this also changes quality utrbad, etc. and output of mrnaseq
    
  return($utrorf,$utrosize);    
}

sub getRevOrf # or revorf_test  # find/report on longest/best orf in reverse of longorf
{
  my ($longorf, $orfs, $cdnalen, $issorted)= @_;
  my $MinCDS= 3*$MINAA; $issorted||=0;
  my ($revorf,$revosize)=(0,0);
  my $lorient= $longorf->{orient};
  # my $lsize= $longorf->{length}; # goodlen ?
	##  orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
	## my($lb,$le)= ($longorf->{start},$longorf->{stop});  ($lb,$le)=($le,$lb) if($lb>$le);
	foreach my $orf (@$orfs) {  # expect long sorted, or check all? 
		my($oorient,$ogood,$osize)= ($orf->{orient},$orf->{goodlen},$orf->{length});
		# $ob,$oe, = $orf->{start},$orf->{stop}, # ($ob,$oe,$oorient)=($oe,$ob,'-') if($ob>$oe);
		if($oorient ne $lorient and $ogood>=$MinCDS and $ogood > $revosize) {  
			$revorf= $orf; $revosize= $ogood; 
			last if $issorted;
			}
   }  
  return($revorf,$revosize);    
}

# replace old sub getCDSgff($exons,$orfprot,$prostart5,$proend3) w/ getCDSgff2($exons,$orf)
sub getCDSgff2 # cdna_proteins.pm
{
  my($exons,$orf)= @_;

  my ($orfprot,$prostart5,$proend3,$orflen0)= orfParts($orf);
  my $cdnalen= $orf->{cdnalen}||0;
#   my($exons,$orfprot,$prostart5,$proend3,$cdnalen) = @_;
#   # FIXME: need trlen= length(cdnain) for -cdna, and/or use gmap qlen= tag
#   # FIXMEs: exon Split= needs copy to CDS .. losing it, here???
    
  my ($cds5,$cds3,$cds3last,$strand)=(0,0,0,'.');
  my @cds= ();
  my @utr= ();
  ## for phase; need reverse @exons
  my ($cdna1,$inc5,$inc3,$nt_length, $nu5, $nu3)= (0) x 10;
  
  $cdna1= 0; # was 1; # offby1 at end?
  $nt_length= 0; # $prostart5 % 3; #??
  # my $KEEPAN='Split|err|gaps'; #? |Target|trg |gapfill|gapfix ?
  my $DROPAN='Target|trg|gapfill|gapfix|splice';
  
  # ** FIXME 2011Dec : stopcodon split intron >> CDS ends w/o final 1,2 bases ** WRONG
  # .. must make next exon(if exists) part of CDS stop
  
  #?? rev bug here? got bad CDS for good prot/cdna after completeCDSb, revgene
  my $addat="";
  if($exons->[0]->[6] eq "-") {
    my @xbeg= @{$exons->[0]}[3,4,6]; # b,e,o
    my @xend= @{$exons->[-1]}[3,4,6];
    if($xbeg[0] < $xend[0]) {
      $addat.=",badrev"; 
      my @xrev= reverse @$exons; $exons= \@xrev; 
      }
  }
  
  foreach my $exon (@$exons) {
    my ($ref,$src,$xtyp,$xend5, $xend3,$xv,$xo,$xph,$xattr,$gid)= @{$exon};
    
    my $cdsattr=$xattr; 
    $cdsattr=~s/;($DROPAN)=[^;\n]+//g; #** Target|trg has spaces **
    $cdsattr=~s/Parent=[^;\s]+[;]?//; #? or leave on should be same as $gid
    $cdsattr="" unless($cdsattr=~/\w+=/);

    $strand= $xo;
    ($xend5, $xend3)= ($xend3,$xend5) if($xo eq "-"); #patch rev?
    my $xd= abs($xend3 - $xend5); # ?? +1 for width
    
    $cdna1++; # add 1 here, not end loop 
    my $cdna2= $cdna1 + $xd; # ?? +1 for width
    # ** offby1 here ?? YES, need <=, >= to get full CDS stop,start split by intron; see below cdna1=cdna2+1
    #OLD.if($cdna1 < $proend3 and $cdna2 > $prostart5) 
    if($cdna1 <= $proend3 and $cdna2 >= $prostart5) 
    { # overlap
                
      my $d5= ($cdna1 >= $prostart5) ? 0 : $prostart5 - $cdna1; # pos
      my $c5= ($xend5 > $xend3) ? $xend5 - $d5 : $xend5 + $d5;    				  
      
      my $d3= ($cdna2 <= $proend3) ? 0 : $proend3 - $cdna2; # neg
      my $c3= ($xend5 > $xend3) ? $xend3 - $d3 : $xend3 + $d3;
  
      my $elength = 1 + abs($c3 - $c5);
      $nt_length  += $elength;
      $inc3        = $nt_length % 3;
      $inc5        = ($elength - $inc3) % 3; # only care about this one
      # $frame       = ($c5 + $inc5) % 3;
      if ($inc5 == -1) { $inc5 = 2; }
      my $phase= $inc5; # is this right?
      
      my($cb,$ce)= ($c5 > $c3) ? ($c3,$c5): ($c5,$c3); #? rev patch
      # exon Split=, other xattr need copy to CDS HERE ***
      my $rloc= [$ref,$src,"CDS",$cb,$ce,".",$xo,$phase,"Parent=$gid;$cdsattr",$gid]; 
      push @cds, $rloc;
     
      if($cdna1 <= $prostart5) { $cds5=$c5; }
      if($cdna2 >= $proend3) { $cds3=$c3; } else { $cds3last= $c3; } # cdna2 < proend3 here
 
    } elsif(1) { # $addutr .. not used?
      my $d5= ($cdna1 >= $proend3) ? 0 : $proend3 - $cdna1; # pos
      my $u5= ($xend5 > $xend3) ? $xend5 - $d5 : $xend5 + $d5;    				  
      my $d3= ($cdna2 <= $prostart5) ? 0 : $prostart5 - $cdna2; # neg
      my $u3= ($xend5 > $xend3) ? $xend3 - $d3 : $xend3 + $d3;

      my($ub,$ue)= ($u5 > $u3) ? ($u3,$u5): ($u5,$u3); 
      my $up= ($cdna1 < $prostart5) ? "five" : ($cdna2 > $proend3) ? "three" : "odd";
      if($cdna1 < $prostart5) { $nu5++; } elsif($cdna2 > $proend3) { $nu3++; }
      my $rloc= [$ref,$src, $up."_prime_utr",$ub,$ue,".",$xo,0,"Parent=$gid;$cdsattr",$gid]; 
      push @utr, $rloc;
    }
    
    ##$cdna1= $cdna2+1;  # is this off-by-1 now? yes, dont +1 here, do above
    $cdna1= $cdna2;  
  }       
  $cds3||=$cds3last; # partial3
  
  use constant forFASTA => 0;
  my $trlen= ($cdnalen>$cdna1) ? $cdnalen : $cdna1;  
  my($aalen,$pcds,$compl,$orflen)  
     = proteindoc($orf,$trlen,$strand,forFASTA);
  # my($aalen,$pcds,$compl,$orflen,$ocdnalen,$ocdsoffs,$oSelc) 
  #   = proteinqual($orf,$trlen);

  my $mattr="cxlen=$orflen/$trlen;aalen=$aalen,$pcds%,$compl";
  $mattr.= ";protein=$orfprot" if($orfprot);
  $mattr.= ";cdsoff=$prostart5-$proend3"; #? as per ;utroff=$ustart-$uend
  $mattr.= ";cdsspan=$cds5-$cds3$addat"; #?add for other uses? rev if needed? or not?
  $mattr.= ";utrx=$nu5,$nu3" if($nu5 > 2 or $nu3 > 2); # ;utrx=$u5,$u3
  ## mattr keys: cxlen,aalen,protein,utrx
  
  #? resort @cds by loc, not reversed. : let caller do
  # @cds = sort _sortgene @cds;

  ## return also: $clen, $trlen or $utrlen or $ap, $nu5+$nu3, 
  ## ($cdslong, $attrL, $pcodeL, $maxutrL)
  return (\@cds, $mattr, $pcds, _max($nu5,$nu3), \@utr); 
}

  
sub getAllOrfs { # added strands to return both, or rev
  my ($input_sequence, $strands,$flags) = @_;
  $strands ||= "fwd";
  $flags ||="";
  # add fwd,rev orfs here ? gff-caller wants only 1 strand at a time, need option 
  
  return undef unless ($input_sequence or length ($input_sequence) >= 3) ;
  $input_sequence = uc ($input_sequence); # was lc() change all to uc()
  
#  ## fixme: this screws seq indices: start,stop; instead of chomp, restrict @starts,@stops
#   if($flags =~ /dropnnn|chomp/i) { $input_sequence =~ s/^N+//; $input_sequence =~ s/N+$//;  }

  my $seqlen=length($input_sequence);
  
  # 201402 fix: orthobest input cds span for bestorf..
  # .. innerstop problem: retained introns looks like cause,
  # .. can add ref/tr hsp list, which gives some info on introns in un-aligned spans..
  # .. then need messy orf-calling using tr-hsp list with extensions up to find intron
  # .. or some other intron splice site seek to cut it out of trasm for mRNA/CDS/aa
  # .. try scan for splice sites fwd:AG/GT (rev:AC/CT )  ?
  
  my $CheckCutIntrons= 1; # debug
  
  my @bestspan=(); my @besthsps=(); # add @besthsps,@refhsps for innerstop,introncuts ?
  my $bestspanorf=undef; my $checkspanframe=0; my $hasHSPs=0;
  my $bestspanflag="";
  if($flags =~ /bestspan=([\w,\:\+\-]+)/) { # \d to \w to get all parts?
    $bestspanflag=$1; # my $bspan=$1;  ## now "bestspan=$tralnspan:$cla:$rid:$refaaspan";
    $checkspanframe= ($flags=~/frameck/)?1:0;
    # FIXME: spans here are 1origin, beststart/getOneOrf want 0origin ??
    my($cbe,$co)=split":", $bestspanflag; ## 2341-1364:-
    ## allow hsps here?  22-333,444-555:+ ?
    if($cbe=~/,\d/) {
      my($mcb,$mce,$lb,$le)=(-1,0,0,0); 
      for my $ch (split",",$cbe) { 
        my($cb,$ce)= split '-',$ch; next unless($ce>0);
        if($co eq '-') { ($cb,$ce)=($seqlen-$ce,$seqlen-$cb); } 
        else { $cb--; $ce--; } # 0origin fix, BUT not rev !!
        if($le>0 and $co eq '-' and $ce>=$lb) { $ce=$le; $besthsps[-2]=$cb if($cb<$lb); }
        elsif($le>0 and $cb<=$le) { $cb=$lb; $besthsps[-1]=$ce if($ce>$le); }
        else { push @besthsps, $cb,$ce; }
        ($lb,$le)=($cb,$ce);
        $mcb=$cb if($mcb==-1 or $cb<$mcb); $mce=$ce if($mce<$ce);
      }
      $hasHSPs=(@besthsps>2)? scalar(@besthsps) : 0;
      ## push @besthsps,$co;       
      @bestspan=($mcb,$mce,$co);
    } else {
    my($cb,$ce)= split /[\-]/, $cbe; 
    if($cb>$ce) { $co='-'; ($cb,$ce)=($ce,$cb); } # doublerev, or (cb,ce)=(seqlen-cb,seqlen-ce);
    if($co eq '-') { ($cb,$ce)=($seqlen-$ce,$seqlen-$cb); } 
    else { $cb--; $ce--; } # 0origin fix !!! BUT not rev !!
    @bestspan=($cb,$ce,$co); # but need to extend to @start,@stop !
    }
  }
  
  my (@starts, @stops, @orfs);
  my $workseq= $input_sequence;
  unless($strands =~ /^(r|\-)/) {  # forward_strand_only();
    @stops  = identify_putative_stops($workseq);
    @starts = identify_putative_starts($workseq,\@stops);
    push @orfs, get_orfs (\@starts, \@stops, $workseq, '+');
    if(@bestspan and $bestspan[2] eq '+') { 
      my ($bb,$be,$bo)=@bestspan; 
      
      # BUG: got shorter aa from bestspan, but same start .. found innerstop here, but not w/o span
      # .. seems to be problem of extended aln, from possibly alt-exon in same trasm: t8453-8854/r2703-2836 adds stops
      # def  =socatfishv1k39loc150483t1 aalen=2617,74%,complete; clen=10510; strand=+; offs=213-8066; 
      # bspan=socatfishv1k39loc150483t1 aalen=2593,74%,complete; clen=10510; strand=+; offs=213-7994;
      #   orflags=bestspan:213-4487,4905-6608,8453-8854:+:clmiss:ref:Funhe2EKm027765t1:1-1314,1399-1972,2703-2836,incut:20-stops:1751-bp; 
      ## aln: Funhe2EKm027765t1  socatfishv1k39loc150483t1  213-8854:+      over:213-8066:+
      
      # FIXME: need to find best frame, bb can be offby 1,2.. maybe not, blastn may be on-frame
      ## fix for stopatend: (substr($orfprot,-1) eq '*')  
      my $aa0= getOneOrf($bb,$be,$workseq,$bo); 
      my $nstop= $aa0->{innerstop}; 
      my $aa0len= $aa0->{goodlen};
      
      # # .. try scan for splice sites fwd:AG/GT (rev:AC/CT )  ?
      if($nstop>0 and $CheckCutIntrons and $hasHSPs ) {  
        my @icut=();
        for(my $i=1; $i<$hasHSPs-1; $i+=2) {
          my($b1,$e1,$b2,$e2)= @besthsps[$i-1,$i,$i+1,$i+2];
          my $ospan= getOneOrf($b1,$e2,$workseq,$bo);  # check if bad segment!
          my $instop= $ospan->{innerstop};
          next unless($instop);

          my $es= index($workseq,'AG',$e1); 
          my $bs= rindex($workseq,'GT',$b2);
          if($es>$e1 and $bs>$es) { 
            push @icut,$es,$bs+2;  
          } else {
            push @icut,$e1+1,$b2;  # yes, this doubles n unstopped cases
          }
        }
        if(@icut) {
          my $bat=0; my $mrna=""; my $ncut= $#icut; my $spancut=0;
          my $cb=$besthsps[0]; my $ce=$besthsps[-1]; # needs adjust..
          for(my $i=0; $i<$ncut; $i+=2) { 
            my($sb,$se)= @icut[$i,$i+1]; 
            $mrna .= substr($workseq,$bat,$sb-$bat);
            $bat = $se; my $cut=$se-$sb; $ce -= $cut; $spancut+=$cut;
            }
          $mrna .= substr($workseq,$bat);  
          my $aac= getOneOrf($cb,$ce,$mrna,$bo); 
          my $nsc= $aac->{innerstop};  
          my $aa1len= $aac->{goodlen};  
          # DEBUG out here..
          warn "#DBG CutIntrons: aa0w=$aa0len, aa1w=$aa1len, stop0=$nstop, stop1=$nsc, hsps=@besthsps, icut=@icut, flags=$flags\n"
            if($DEBUG);
          ## this works for ~50/500 cases in catfish1evg8.. enough to use? any other fixes?
          ## using full cut e1..b2 ups that to 120/500 fixed no innerstop, usable
          ## .. possibly test for stops in each mrnasegment add, scan over istop? ie shift se=>bat above?
          ## .. some of crap is not full intron insert, but partial .. test cut e1..b2 w/o AG/GT scan?
          
          if($nsc < 1) {
            $bestspanflag.=",incut:$nstop-stops:$spancut-bp";
            ($bb,$be)= ($cb,$ce); $nstop=$nsc; 
            $workseq= $mrna; $seqlen=length($mrna);
            @stops  = identify_putative_stops($workseq);
            @starts = identify_putative_starts($workseq,\@stops);
          }
        }
      }   
      
      if($nstop>0 and $checkspanframe) {
        my $aa1= getOneOrf($bb+1,$be,$workseq,$bo);  my $ns1= $aa1->{innerstop};  
        my $aa2= getOneOrf($bb+2,$be,$workseq,$bo);  my $ns2= $aa2->{innerstop};  
        my $f=0; 
        if($ns1 < $nstop) { $f=1; $nstop=$ns1; } 
        if($ns2 < $nstop) { $f=2; }
        $bb += $f;
      }
      
      $bb=beststart($bb,0,\@starts,\@stops);
      $be=beststop($be,$bb,$seqlen,\@stops);
      $bestspanorf = getOneOrf($bb,$be,$workseq,$bo);   @bestspan= ($bb,$be,$bo);
      $bestspanorf->{flags}= "bestspan:$bestspanflag";
      # unshift @orfs,$bestspanorf;  # check all orfs for this tho, no dupl.
      }
    }
    
  if($strands =~ /^(b|r|\-)/) { # rev|r|both
    my $revseq= revcomp($input_sequence); 
    $workseq= $revseq;
    @stops  = identify_putative_stops($workseq);
    @starts = identify_putative_starts($workseq,\@stops);
    push @orfs, get_orfs (\@starts, \@stops, $workseq, '-');
    if(@bestspan and $bestspan[2] eq '-') { 
      my ($bb,$be,$bo)=@bestspan; 
      # FIXME: need to find best frame, bb can be offby 1,2..
      my $aa0= getOneOrf($bb,$be,$workseq,$bo);
      my $nstop= $aa0->{innerstop}; # {protein} =~ tr/*/*/;
      if($nstop>0 and $checkspanframe) {
        my $aa1= getOneOrf($bb+1,$be,$workseq,$bo);  my $ns1= $aa1->{innerstop}; #{protein} =~ tr/*/*/;
        my $aa2= getOneOrf($bb+2,$be,$workseq,$bo);  my $ns2= $aa2->{innerstop}; #{protein} =~ tr/*/*/;
        my $f=0; if($ns1 < $nstop) { $f=1; $nstop=$ns1; } 
        if($ns2 < $nstop) { $f=2; }
        $bb += $f;
      }

      $bb= beststart($bb,0,\@starts,\@stops);
      $be= beststop($be,$bb,$seqlen,\@stops);
      $bestspanorf = getOneOrf($bb,$be,$workseq,$bo);  @bestspan= ($bb,$be,$bo);
      $bestspanorf->{flags}= "bestspan:$bestspanflag";
      # unshift @orfs,$bestspanorf; # unless grep{$_ eq $bestorf}@orfs; # check all orfs for this tho, no dupl.
      }
  }

  
  if (@orfs or $bestspanorf) {  # get both longest complete, longest partial
    ## set in order of decreasing length
    my $keylen=($USEGOODLEN)?"goodlen":"length";
    @orfs = sort {$b->{$keylen} <=> $a->{$keylen} or $b->{complete} <=> $a->{complete}} @orfs;
    
    # FIXME: for bestspan.. but cancel if innerstop > 0 or > 1 .. esp. scrambled bestspan's w/ ## stops
    # FIXME2: cancel if bestspan goodlen << longest goodlen, happening w/ odd hsps, maybe alt exons; need ortho vals
    if($bestspanorf) {
      my $nstop= $bestspanorf->{innerstop}; 
      if($nstop > 0 and $InnerStopToX<0) { 
         # skip but add info to orf0
        my $ncds= $bestspanorf->{length};
        if(@orfs) { $orfs[0]->{flags}= "cancelspan:innerstop.$nstop,cdsw.$ncds,$bestspanflag"; } #??
      } else { 
        if($nstop > 0 and $InnerStopToX==1) { substr( $bestspanorf->{protein}, 0,-1) =~ tr/*/X/; }
        unshift @orfs,$bestspanorf; 
        # $longest = $bestspanorf;
      }
    }
    
    # FIXME: option to mark/return 2ndary orf(s) in aberrant long-utr transcripts, that dont overlap 1st orf
    # FIXME 1505, partial 1aa > partial3 bug from utrorf.mrna cut just past partial3 cds.
    # ?? look for 2nd longest, complete > longest.complete ?

    my $longest = $orfs[0];   
    if($longest->{innerstop}) { # 201504 fix: drop innerstops unless desired
      unless($InnerStopToX > 0) { ($longest) = grep { $_->{innerstop} == 0 } @orfs; }
      #?? if($InnerStopToX <= 0) { ($longest) = grep { $_->{innerstop} == 0 } @orfs; }
      # else { } # convert to X
    }
    my($longfull) = grep { $_->{complete} == 3 && $_->{innerstop} == 0 } @orfs; # add  $_->{innerstop} == 0
    
    ## upd1712: cancel this for $flag =~ /long/ 
    ## UPD1908: $ORF_FULLvPART can be too low for this change, use constant ~ 0.90 to pick shorter part5,3
    ## what is utrorf partial53 problem from 2015.05 that this fixes?
    if($longest->{complete} == 0 and not ($flags =~ m/long/)) { # utrorf bug 1505
      my $pLONGPART= 0.90; ## (UPD1908) ? 0.90 : $ORF_FULLvPART; # dont use $ORF_FULLvPART at all? 
      my($loc)= grep { $_->{complete} > 0 &&  $_->{innerstop} == 0 } @orfs; 
      if(not $loc) { } elsif($longfull and $longfull eq $loc) { } # avoid undefs
      elsif($loc->{goodlen} > $pLONGPART * $longest->{goodlen}) { $longest=$loc; }
    }

    return ($longest, $longfull, \@orfs);
  } else {
    return undef;
  }
}

=item newOrf

  my $orf= newOrf( sequence => 'acgtcgat', protein => "MABCD", complete =>  3, start => 111, stop => 999, orient => '+', length => 1999 );

  see get_orfs() makes these from data
        my $orf = { sequence => $orfSeq, protein => $protein,
          start=>$start, stop=>$stop, orient=>$direction,
          length=>$orflen, goodlen=>$goodlen, cdnalen=>$seq_length,
          complete=>$isfull, innerstop=>$innerstop,
          };

   name? orfOf(xxx) or  newOrf(xxx)
   
=cut 

my @ORFPARTS= qw(sequence protein start stop orient length goodlen cdnalen complete innerstop);
# add Selc ? other sometimes keys?

sub newOrf { 
  my(%parts)= @_; 
  my %orf=(); unless(%parts){ %parts= (); }
  my $ORP=join '|', @ORFPARTS;
  for my $k (@ORFPARTS){ $orf{$k}= $parts{$k}||0; }
  for my $k (grep{ not m/$ORP/ } keys %parts) { $orf{$k}= $parts{$k}; } #? add all parts not std?

  # check data: start,stop required, stop < start  when -orient 
  # minimal: newOrf( protein=>"xxxx", start=>99,  stop=>999, cdnalen=>1999);
  #  NoStopCodon problems, *should* only chop for output
  #   if(substr($orfprot,-1,1) eq '*') { $aalen--; if($NoStopCodon) { $orfprot =~ s/\*$//; } }

  my $len= ($orf{stop} < $orf{start})? 1 + $orf{start} - $orf{stop} : 1 + $orf{stop} - $orf{start}; 
  $orf{length}= $len if($orf{length} < $len);
  $len= $orf{length};
  $orf{goodlen} ||= $len;
  #? $orf{cdnalen}= $len if($orf{cdnalen} < $len); #? or 0 for missing?
  $orf{orient} ||= '.'; $orf{protein} ||= '';  $orf{sequence} ||= '';
  
  if(my $protein= $orf{protein}) {
    my $isfull= $orf{complete}; # either comp|2 or protein* here
    if(substr($protein,0,1) eq 'M') { $isfull |= 1; }
    if(substr($protein,-1,1) eq '*') { $isfull |= 2; } elsif($NoStopCodon) { } #? need to trust {complete}
    my $innerstop= $protein =~ tr/*/*/;   
    $innerstop-- if($innerstop>0 and ($isfull & 2));  
    $orf{complete}= $isfull if($isfull > $orf{complete});
    $orf{innerstop}= $innerstop if($innerstop>0 and $orf{innerstop}<1);
    }
  
  return(\%orf);
}

sub orfParts
{
  my($orf,$parts) = @_;
  # $parts= [ qw(protein start stop length orient) ] unless(ref $parts); 
  if(ref $orf) {
    if(ref $parts) { return @$orf{ @$parts }; }
    return($orf->{protein}, $orf->{start}, $orf->{stop}, $orf->{length}, $orf->{orient});  
  } else {
    return("",0,0,0,0);
  }
}


sub beststart {
  my( $instart,$minpos,$starts_ref,$stops_ref) = @_;
  my $inframe = $instart % 3; 
  my $spat=$instart;
  ## fixme need to stop at stops before starts
  my @startrev= reverse grep { $_ <=  $instart } @{$starts_ref};
  my @stoprev = reverse grep { $_ <=  $instart } @{$stops_ref};
  my $stopat=0;
  foreach my $ss (@stoprev) { if($ss % 3 == $inframe) { $stopat= $ss; last; } }
  foreach my $sp (@startrev) { if($sp<$stopat) { last; } elsif($sp % 3 == $inframe) { $spat= $sp; last; } }
  # $spat=$stopat+3 if($spat < $stopat);
  return $spat;
#  
#  #......
#   foreach my $sp (@{$starts_ref}) {
#     last if($sp > $instart); next if($sp<$minpos);
#     $spat= $sp if($sp % 3 == $inframe);
#   }
#  return $spat;
}

sub beststop {
  my( $inend,$instart,$maxpos,$stops_ref) = @_;
  my $inframe = $instart % 3; 
  my $spat=$inend;
  foreach my $sp (@{$stops_ref}) {
    last if($sp > $maxpos); next if($sp < $inend - 3);
    if($sp % 3 == $inframe) { $spat= $sp; last; }
  }
  return $spat;
}

sub getOneOrf {
  my( $start_pos, $stop_pos, $mrnaseq, $direction, $flags) = @_;
  
  #* Need option to extend start,stop; ** Need to look for start codon, stop codon if start/stop not on these.
  #* Assume seq == revcomp(seq) if dir eq '-' ??
  
  my $seq_length = length ($mrnaseq);
  my $sp_frame = $start_pos % 3;
  
  ## FIXME: ?? bad for dir=-
#  my $isfull5= ($start_pos+3 <= $seq_length and substr($mrnaseq, $start_pos, 3) eq "ATG")?1:0;
#  
#  #?? or use: foreach my $start_pos (@{$starts_ref})
#  ## use also: foreach my $stop_pos (@{$stops_ref}) 
#   if(!$isfull5 and $start_pos > 2 and $flags =~ /full|seekstart/i ) { #??
#     for (my $b=$start_pos - 3; $b > 2 && !$isfull5; $b -= 3) {
#       if(substr($mrnaseq, $b, 3) eq "ATG") { $isfull5=1; $start_pos=$b; last; }
#       }
#   }
#  
#   my $start_pos_frame = $start_pos % 3;
#   if($isfull5) { $last_stop_full{$start_pos_frame} = $stop_pos; }
#   else { $last_stop_part{$start_pos_frame} = $stop_pos; }

  my $stopplus3= _min( $stop_pos+3, $seq_length); #dgg patch; getting overruns.
  my $orflen = ($stopplus3 - $start_pos);
  my $stopoff = $orflen % 3;
  if($stopoff>0) {  $orflen -= $stopoff; $stopplus3 -= $stopoff; } # do before pull seq..
  my ($start_pos_adj, $stop_pos_adj) = ( ($start_pos+1), $stopplus3);
  my ($start, $stop) = ($direction eq '-') 
    ? (revcomp_coord($start_pos_adj, $seq_length), revcomp_coord($stop_pos_adj, $seq_length))
    : ($start_pos_adj, $stop_pos_adj);
  
  my $orfSeq =  substr ($mrnaseq, $start_pos, $orflen); #include the stop codon too.
  my $protein= translate_sequence($orfSeq, 1); ## from Nuc_translator
  $orflen= length($orfSeq);

  my $isfull= 0;
  $isfull |= 1 if(substr($protein,0,1) eq 'M');  
  $isfull |= 2 if(substr($protein,-1,1) eq '*');  
  my $innerstop= $protein =~ tr/*/*/;   
  $innerstop-- if($innerstop>0 and ($isfull & 2));  
  my $nxxx= $innerstop + $protein =~ tr/X/X/;
  my $goodlen= $orflen - 3*$nxxx; # goodlen to avoid XXXXXXX* crap
  
  if($DEBUG>1){  
    my $aalen = length($protein); my $u=index($protein,'u');
    print STDERR "#DBG: orf1,olen=$orflen,alen=$aalen,at=$start-$stop,selc=$u,aa1=",
      substr($protein,0,9),"..",substr($protein,-3,3),"\n";
  }
  
  # orf maybe add: mRNAseq => $mrnaseq, as output help, e.g. for rev and intron-cut seqs
  my $orf = { sequence => $orfSeq, protein => $protein, 
    mrna => $mrnaseq, #? only for getOneOrf 
    start=>$start, stop=>$stop, orient=>$direction,
    length=>$orflen, goodlen=>$goodlen, cdnalen=>$seq_length,
    complete=>$isfull, innerstop=>$innerstop,
    };

  # USESelenocysteine add Selc position in orf for exception reporting, NCBI /transl_except=(pos:1002..1004,aa:Sec); in mRNA span
  if($USESelenocysteine and (my $iu=index($protein,'u')) > 0) {
    my $Selc=""; ## do any have >1 u ?? yes
    while($iu >= 0) { 
      my $uoff=$start_pos_adj + 3*$iu;
      $uoff= revcomp_coord($uoff, $seq_length) if($direction eq '-');
      $Selc.= "$uoff,"; $iu= index($protein,'u',$iu+1);  
      }
    $Selc=~s/,$//; $orf->{Selc}=$Selc;
  }
    
  return $orf;
  # return (wantarray) ? ($orf, $orf, [$orf]) : $orf;
}




sub get_orfs {
  my( $starts_ref, $stops_ref, $mrnaseq, $direction) = @_;
  
  # $direction here only affects start,stop relative to forward cdna seq .. use it.
  # .. input seq is already rev(seq) for -dir
  # Ooops, orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
  # .. fix some places, leave that way others?
    
  # unless ($starts_ref && $stops_ref && $mrnaseq && $direction) { warn "get_orfs: params not appropriate"; }  
  # want only max orf, complete + partial
  # FIXME: option to mark/return 2ndary orf(s) in aberrant long-utr transcripts, that dont overlap 1st orf
   
	my %last_stop_part = ( 0=>-1, 1=>-1,  2=>-1); # position of last chosen stop codon in spec reading frame.
	my %last_stop_full = ( 0=>-1, 1=>-1,  2=>-1); # position of last chosen stop codon in spec reading frame.
  my @orfs;
  my $seq_length = length ($mrnaseq);
  my $norf=0;
  
  ##  $mrnaseq .= "###" if($DEBUG); # why're we getting bad prot at partial end? isnt prot but cds range needs adjust for -2,-1
  
  foreach my $start_pos (@{$starts_ref}) {
		my $start_pos_frame = $start_pos % 3;
		# includes partial5: starts at 0,1,2 unless real start
		# ** PROBLEM? -- for internal mis-assemblies, should
		#   test partial-starts also following each stop? dont require ATG start inside?
		
		my $isfull5= ($start_pos+3 <= $seq_length and substr($mrnaseq, $start_pos, 3) eq "ATG")?1:0;
		
		foreach my $stop_pos (@{$stops_ref}) {
		  if ( ($stop_pos > $start_pos) && #end3 > end5
				 ( ($stop_pos - $start_pos) % 3 == 0) && #must be in-frame
				   # ($start_pos > $last_delete_pos{$start_pos_frame}) # dgg: bad for partial5 + full
				 ($isfull5 ? $start_pos > $last_stop_full{$start_pos_frame} : $start_pos > $last_stop_part{$start_pos_frame} )
				 ) #only count each stop once.
			{
		    # includes partial3: stops at end- 0,1,2 unless real stop
				
				#dgg.no# $last_delete_pos{$start_pos_frame} = $stop_pos;
			  my $oldlast_stop_full= $last_stop_full{$start_pos_frame}; # TRIMCDS fixup
				if($isfull5) { $last_stop_full{$start_pos_frame} = $stop_pos; }
			  else { $last_stop_part{$start_pos_frame} = $stop_pos; }
			
			  my $stopplus3= _min( $stop_pos+3, $seq_length); #dgg patch; getting overruns.
				my $orflen = ($stopplus3 - $start_pos);
				my $stopoff = $orflen % 3;
				if($stopoff>0) {  $orflen -= $stopoff; $stopplus3 -= $stopoff; } # do before pull seq..
			  
				my ($start_pos_adj, $stop_pos_adj) = ( ($start_pos+1), $stopplus3);
				# my ($start, $stop) = ($direction eq '+') ? ($start_pos_adj, $stop_pos_adj) 
				#	: (&revcomp_coord($start_pos_adj, $seq_length), &revcomp_coord($stop_pos_adj, $seq_length));
        my ($start, $stop) = ($direction eq '-') 
          ? (revcomp_coord($start_pos_adj, $seq_length), revcomp_coord($stop_pos_adj, $seq_length))
          : ($start_pos_adj, $stop_pos_adj);

        ## reuse sub getOneOrf() above here ?				
	      my $orfSeq =  substr ($mrnaseq, $start_pos, $orflen); #include the stop codon too.
				
        # FIXME here? gap trim NNN/XXX at start/stop of coding seq, adjusting start/stop
        # .. this is causing bug for MxxxxxABCD cases.. > innerstops despite block; phase problem; make i%3 == 0 
	      # * should do this before/in getStarts() otherwise trim can miss M a few aa after gap.

			  if($TRIMCDSENDGAP) {
			    my $at= index($orfSeq,'NNN');
			    if($at >= 0 and $at < 4) {
			      my $i=$at+2; while($i<$orflen and substr($orfSeq,$i,1) eq 'N') { $i++; }
 			      $i++ while($i % 3 > 0); # 1505 patch solves innerstop bug
			      $orfSeq= substr($orfSeq,$i); 
			      $orflen= length($orfSeq);
			      $start += $i; $start_pos+= $i;
			## more fixups:
			      my $newfull5= (substr($orfSeq, 0, 3) eq "ATG")?1:0;
                              if($isfull5 and not $newfull5) { 
			         $last_stop_full{$start_pos_frame}= $oldlast_stop_full; ##not delete $last_stop_full{$start_pos_frame};  
                                 $last_stop_part{$start_pos_frame} = $stop_pos; $isfull5=$newfull5;
				}

			    }
			    # $at= rindex($orfSeq,'NNN'); # problems w/ stop codon change ..
			    # if($at >= $orflen - 3) { }
			  }
			  	
				my $protein= translate_sequence($orfSeq, 1); ## from Nuc_translator
			  $orflen= length($orfSeq);
  
				my $isfull= 0;
        $isfull |= 1 if(substr($protein,0,1) eq 'M'); # $protein =~ /^M/
        $isfull |= 2 if(substr($protein,-1,1) eq '*'); # $protein =~ /\*$/
 	      my $innerstop= $protein =~ tr/*/*/;   
        $innerstop-- if($innerstop>0 and ($isfull & 2));  
        ## flag problem if XXX count > non-X count
        my $nxxx= $innerstop + $protein =~ tr/X/X/;
        my $goodlen= $orflen - 3*$nxxx;
        $norf++;
        # opt to drop *stopcodon, here or where, 

             
        my $orf = { sequence => $orfSeq, protein => $protein,
          start=>$start, stop=>$stop, orient=>$direction,
          length=>$orflen, goodlen=>$goodlen, cdnalen=>$seq_length,
          complete=>$isfull, innerstop=>$innerstop,
          };

if(TESTPOLYA) {          
        sub haspolya {
          my($mrna,$stop_pos,$seq_length)=@_;
          my $iaaa= index($mrna,'AAAAAAAAAAAAAAAAAAAAA',$stop_pos); # 21 A = 7 K; 'AAAAAAAAAAAAAAAAAA' 18 A = 9 K
          return ($iaaa<0) ? 0 : 1+$iaaa;
        }   
        if( my $pa= haspolya($mrnaseq, $stop_pos, $seq_length) ){ $orf->{polya}=$pa; }
}
        
        # add Selc position in orf for exception reporting, NCBI /transl_except=(pos:1002..1004,aa:Sec); in mRNA span
        if($USESelenocysteine and (my $iu=index($protein,'u')) > 0) {
          my $Selc=""; ## do any have >1 u ?? yes
          while($iu >= 0) { 
            my $uoff=$start_pos_adj + 3*$iu;
            $uoff= revcomp_coord($uoff, $seq_length) if($direction eq '-');
            $Selc.= "$uoff,"; $iu= index($protein,'u',$iu+1);  
            }
          $Selc=~s/,$//; $orf->{Selc}=$Selc;
        }
          
				push (@orfs, $orf);
				
        if($DEBUG>1) { 
          print STDERR "#\n" if($norf==1);  
  				my $aalen = length($protein); my $u=index($protein,'u'); my $or=$direction;
          my $cpd=""; if(UPD1908) {
            my($cpv, $cptyp)= likely_coding_hexamer("orf$norf", $orf);
            $cpd=",cpot=$cpv".substr($cptyp,0,1);
          }
          print STDERR "#DBG: orf$norf,olen=$orflen,alen=$aalen,at=$start-$stop:$or,selc=$u$cpd,aa1=",
            substr($protein,0,9),"..",substr($protein,-3,3),"\n";
        }
        
				last;   
			}
		}
  }
  return (@orfs);
}



sub identify_putative_starts {
  my ( $seq, $stops_aref) = @_;
  my %starts;
  my %stops;
  foreach my $stop (@$stops_aref) {
		$stops{$stop} = 1;
    }
	
  if(1) {  # ($self->{ALLOW_5PRIME_PARTIALS} || $self->{ALLOW_NON_MET_STARTS}) 
    my $i=0;
    if($seq =~ /^N/) {  # FIXME: skip leading NNN
      my $n= length($seq);
      $i++ while( $i<$n && substr($seq,$i,1) eq 'N');
    }
		$starts{$i} = 1 unless $stops{$i};
		++$i; $starts{$i} = 1 unless $stops{$i};
		++$i; $starts{$i} = 1 unless $stops{$i};
    }
    
  if (1) {  # ! $self->{ALLOW_NON_MET_STARTS}  #Look for ATG start codons.
		my $start_pos = index ($seq, "ATG");
		# cases of MxxxxxABCD gaps are problem for a few.. cancel if ATGnnnnnnn ?
		while ($start_pos != -1) {
			$starts{$start_pos} = 1;
			$start_pos = index ($seq, "ATG", ($start_pos + 1));
		  }
    } 

  my @starts = sort {$a<=>$b} keys %starts;
  return (@starts);
}


sub identify_putative_stops {
  my ($seq) = @_;
  my %stops;
  
  if(1){  # $self->{ALLOW_3PRIME_PARTIALS}
	## count terminal 3 nts as possible ORF terminators.
	my $e = length ($seq);
  if($seq =~ /N$/) { # FIXME: skip trailing NNN
    $e-- while( $e>1 && substr($seq,$e-1,1) eq 'N');
  }
	$stops{$e} = 1;
	$e--; $stops{$e} = 1;
	$e--; $stops{$e} = 1;
  }
  
  # my @stop_codons = @{$self->{stop_codons}};
  # my @stop_codons = &Nuc_translator::get_stop_codons(); # live call, depends on current genetic code.
  #global# my @stop_codons = qw(TAA TAG TGA);
## add frameshift detect option here?  only if remaining seq >> nnn, only if %cds/utr falls below pCDSpoor ?
## need indel max to test: offby -2,-1,1,2 only to shift $i
  
  foreach my $stop_codon (@stop_codons) {
    my $stop_pos = index ($seq, $stop_codon);
    while ($stop_pos != -1) {
      $stops{$stop_pos} = 1;
      $stop_pos = index ($seq, $stop_codon, ($stop_pos + 1)); #include the stop codon too.
      }
  }
  my @stops = sort {$a<=>$b} keys %stops;
  return (@stops);
}


sub revcomp {
    my ($seq) = @_;
    my $reversed_seq = reverse ($seq);
    $reversed_seq =~ tr/ACGTacgtyrkm/TGCAtgcarymk/;
    return ($reversed_seq);
}


sub revcomp_coord {
    my ($coord, $seq_length) = @_;
    return ($seq_length - $coord + 1);
}

sub backtranslate_protein {
  my ($sequence, $useAmbiguousNucs) = @_;
  $sequence = uc ($sequence); # redundant ??
  my $seq_length = length ($sequence);
  my $cds_sequence="";
  if(defined $useAmbiguousNucs) {
    if( $useAmbiguousNucs ) { %amino2codon_table= %amino2codon_ambig ;  }
    else { %amino2codon_table= %amino2codon_fixed; } 
  }
  
  for (my $i = 0; $i < $seq_length; $i++) {
      my $codon;
      my $amino = substr($sequence, $i, 1);  # deal with non-alpha
      if (exists($amino2codon_table{$amino})) {
        $codon = $amino2codon_table{$amino};
      } elsif($amino =~ /[A-Z]/) {
        $codon = 'NNN'; # fixme
      } else {
        $codon= ''; # eat it?
      }
      $cds_sequence .= $codon;
  }
  return($cds_sequence);
}
      

sub translate_sequence {
  my ($sequence, $frame) = @_;
    
  $sequence = uc ($sequence); # redundant now
	$sequence =~ tr/U/T/;
  my $seq_length = length ($sequence);
  unless ($frame >= 1 and $frame <= 6) { 
		warn "Frame $frame is not allowed. Only between 1 and 6"; # die?
		return -1; # 
	}
	
	if ($frame > 3) {
		# on reverse strand. Revcomp the sequence and reset the frame
		$sequence = revcomp($sequence);
		if ($frame == 4) {
			$frame = 1;
		}
		elsif ($frame == 5) {
			$frame = 2;
		}
		elsif ($frame == 6) {
			$frame = 3;
		}
	}
	
  # $sequence =~ tr/T/U/; # dont need this; change codon_table
  my $start_point = $frame - 1;
  my $protein_sequence="";
  for (my $i = $start_point; $i < $seq_length; $i+=3) {
      my $codon = substr($sequence, $i, 3); # problem here for i>seq_length-3 ?? or in caller getting +1,+2 > true len

## add frameshift detect option here?  if codon = stop_codons ; need it above other places
## need indel max to test: offby -2,-1,1,2 only to shift $i

      my $amino_acid;
      if (exists($codon_table{$codon})) {
          $amino_acid = $codon_table{$codon};
      } else {
          if (length($codon) == 3) {
              $amino_acid = 'X';
          } else {
              $amino_acid = "";
          }
      }
      $protein_sequence .= $amino_acid;
  }
  return($protein_sequence);
}

sub useSelenocysteine
{
  my ($turnon)= @_;
  $turnon=1 unless(defined $turnon);
  my $lastu= $USESelenocysteine; 
  unless($turnon) { # ?? dont reset global $USESelenocysteine .. BUT above uses,
    @stop_codons = qw(TAA TAG TGA); # or push @stops, 'TGA' unless(grep/TGA/,@stops);
    $codon_table{'TGA'} = '*';
    $currentCode = "universal"; $USESelenocysteine=0;
  } else {
    @stop_codons = grep !/TGA/, @stop_codons; #  qw(TAA TAG); # not TGA
    $codon_table{'TGA'} = 'u';
    $currentCode = "universalSelC"; $USESelenocysteine=1;
  }
  return $lastu;
}

our $USE_TGAw=0; 
sub useTGAw
{
  my ($turnon)= @_;
  $turnon=1 unless(defined $turnon);
  unless($turnon) {
    @stop_codons = qw(TAA TAG TGA); 
    $codon_table{'TGA'} = '*';
    $currentCode = "universal"; $USE_TGAw=0;
  } else {
    @stop_codons = qw(TAA TAG); # not TGA
    $codon_table{'TGA'} = 'w';
    $currentCode = "universal_TGAw"; $USE_TGAw=1;
  }
  return $turnon;
}

#==== UPD 20190829 llHexcode ====

=item about likely_coding_hexamer

  hexcode table is of 4096 codon-pair hexamers, as log likelihood of coding/noncoding,  
  from max lcn of four ref species {human,mouse,zfish,fly}hexamer.tsv
 
  likely_coding($cds) is sum of hexcode ll vals in cds, coding >0, noncode <0, with uncertain cutoffs
 
  perl -ne '($hex,$cd,$nc)=split; next if(/^hex/); $lf=log( (1+$cd)/(1+$nc)); 
  if($ov=$hexv{$hex}) { $hexlo{$hex}= ($lf>$ov)?$ov:$lf; $lf=$ov if($ov>$lf); } $hexv{$hex}=$lf; 
  END { for $hex (sort keys %hexv) { $lf=$hexv{$hex}; $lmin=$hexlo{$hex}||$lf;  
  printf "%s\t%.2g\t%.2g\n", $hex,$lf,$lmin; } } ' \
     $gs/codingpot/cpatf/dat/*_Hexamer.tsv
     
  2014 codingpot/cpatf/dat/Human_Hexamer.tsv
  2014 codingpot/cpatf/dat/Mouse_Hexamer.tsv
  2014 codingpot/cpatf/dat/fly_Hexamer.tsv
  2014 codingpot/cpatf/dat/zebrafish_Hexamer.tsv

=cut

# above in vars: our %llHexcode;
use constant { HNoncode => -1e-4, HCoding => 1e-4 };  # hexamer score, dunno what range to call

# usage: ($codeval, $codetype)= likely_coding_hexamer($id, $orfORcds);
sub likely_coding_hexamer {
  my($id, $orfORcds)=@_;
  my($cds,$isref)=(0,0);
  if(ref($orfORcds)) { 
    if(my $val= $orfORcds->{codepot}) { 
      my $ftype= ($val <= HNoncode) ? "Noncode" : ($val >= HCoding)? "Code" : "Unknown";
      $val= sprintf("%.2g",$val);
      return ($val,$ftype,99);
    } 
    $cds= $orfORcds->{sequence}; $isref=1;
  }
  elsif($orfORcds=~/[ACGTacgt]/) { $cds=$orfORcds; } 
  else { return(0); }
  my $n=length($cds); $cds=uc($cds); # dont assume
  my($val,$nx)=(0,0);
  for(my $i=0; $i < $n-6; $i+=3) {
    my $hex=  substr($cds,$i,6); # assume uc(cds)
    if(my $h= $llHexcode{$hex}) { $val+=$h; $nx++; }  # loglike value of coding+1 /noncode -1,    
    }

  my $ftype= (($val <= HNoncode) ? "Noncode" : ($val >= HCoding)? "Code" : "Unknown");
  #o $val= sprintf("%.4g",$val); # chomp excess -0.000265000000000001, 4.59999999999692e-06
  $orfORcds->{codepot}= $val if($isref);
  $val= sprintf("%.2g",$val); # chomp excess -0.000265000000000001, 4.59999999999692e-06
  return($val,$ftype,$nx); # drop nx
}

#============= protein code DATA ====================

BEGIN {
  $currentCode = "universal"; 
  @stop_codons = qw(TAA TAG TGA);
  # stops: TAG  TGA  TAA; Note: TGA also == Selenocysteine valid non-stop
  ## add/allow N in silent positions; need other codon table w/ extended AA vals for ambiguous
  %codon_table = (    
  TTT => 'F', TTC => 'F', TTA => 'L', TTG => 'L',
  CTT => 'L', CTC => 'L', CTA => 'L', CTG => 'L',  CTN => 'L',
  ATT => 'I', ATC => 'I', ATA => 'I', ATG => 'M',
  GTT => 'V', GTC => 'V', GTA => 'V', GTG => 'V',  GTN => 'V',
  TCT => 'S', TCC => 'S', TCA => 'S', TCG => 'S',  TCN => 'S',
  CCT => 'P', CCC => 'P', CCA => 'P', CCG => 'P',  CCN => 'P',
  ACT => 'T', ACC => 'T', ACA => 'T', ACG => 'T',  ACN => 'T',
  GCT => 'A', GCC => 'A', GCA => 'A', GCG => 'A',  GCN => 'A',
  TAT => 'Y', TAC => 'Y', TAA => '*', TAG => '*',
  CAT => 'H', CAC => 'H', CAA => 'Q', CAG => 'Q',
  AAT => 'N', AAC => 'N', AAA => 'K', AAG => 'K',
  GAT => 'D', GAC => 'D', GAA => 'E', GAG => 'E',
  TGT => 'C', TGC => 'C', 
    TGA => '*',  # alternate UGA => 'U' Sec / SeC Se-Cys
  TGG => 'W',
  CGT => 'R', CGC => 'R', CGA => 'R', CGG => 'R',  CGN => 'R',
  AGT => 'S', AGC => 'S', AGA => 'R', AGG => 'R',
  GGT => 'G', GGC => 'G', GGA => 'G', GGG => 'G',  GGN => 'G', 
  );
  

  #DEBUG# $USESelenocysteine=1; #  # fixme, need to define it off? or expect caller to?
  $USESelenocysteine=0 unless(defined $USESelenocysteine); # fixme, need to define it off? or expect caller to?
  useSelenocysteine($USESelenocysteine) if($USESelenocysteine); # no?
  # useTGAw(0); # DEBUG
  
  # should use freq/codon or other choice of best codon per aa.
  # our %amino2codon_ambig,%amino2codon_fixed
  
  # use constant backtrans_P2C => 1;  
  # backtrans pep2codon from pal2nal.pl
  # ambiguous dna codes: R=AG, Y=TC, 
  %amino2codon_ambig = (
    "F" => "TTy", # (T|C|Y)
    "L" => "YTn", # "(CT.)|(TTR)", # A|G|R
    "I" => "ATn", # (T|C|Y|A)
    "M" => "ATG",
    "V" => "GTn",
    "S" => "TCn", # "(TC.)|(AGY)", # (T|C|Y)
    "P" => "CCn",
    "T" => "ACn",
    "A" => "GCn",
    "Y" => "TAy", # (T|C|Y)
    "*" => "Trn", # "(TAR)|(TGA)", # (A|G|R)
    #"_" => "(TA(A|G|R))|(TGA)",
    "H" => "CAy", # (T|C|Y)
    "Q" => "CAr", # (A|G|R)
    "N" => "AAy", # (T|C|Y)
    "K" => "AAr", # (A|G|R)
    "D" => "GAy", # (T|C|Y)
    "E" => "GAr", # (A|G|R)
    "C" => "TGy", # (T|C|Y)
    "W" => "TGG",
    "R" => "CGn", # "(CG.)|(AGR)", # (A|G|R)
    "G" => "GGn",
    "X" => "nnn",
    );

  foreach my $codon (sort keys %codon_table) {
    my $aa= $codon_table{$codon};
    $amino2codon_fixed{$aa}= $codon unless($amino2codon_fixed{$aa}); # multiple codons/aa, what to do?
  }
 
if(backtrans_AmbiguousNucs) { # default ; now opt to backtranslate_protein
  %amino2codon_table= %amino2codon_ambig;  
} else {
  %amino2codon_table= %amino2codon_fixed;  
}  

# UPD 201908
# too maxy: maximum codon-pair stat log((1+code)/(1+noncode)) of cpat/dat/{Human,Mouse,zebrafish,fly}_hexamer.tsv
# average of 4 ref species: maybe right..
  %llHexcode= ( # %llHexcode_ave
AAAAAA => -0.0017, AAAAAC => -0.00033, AAAAAG => 7.3e-05, AAAAAT => -0.00072, AAAACA => -0.00057, AAAACC => -2.1e-05, 
AAAACG => -0.00012, AAAACT => -0.00028, AAAAGA => -0.00038, AAAAGC => -0.00015, AAAAGG => -0.00019, AAAAGT => -0.00031, 
AAAATA => -0.00083, AAAATC => -9.9e-05, AAAATG => -0.00029, AAAATT => -0.00047, AAACAA => -0.00053, AAACAC => -0.00018, 
AAACAG => 0.00012, AAACAT => -0.00033, AAACCA => -4.3e-05, AAACCC => 9.9e-05, AAACCG => 1.5e-05, AAACCT => 2.5e-05, 
AAACGA => -5.4e-05, AAACGC => 5.7e-05, AAACGG => 2.5e-05, AAACGT => -6.6e-05, AAACTA => -0.00019, AAACTC => 6.1e-05, 
AAACTG => 0.00022, AAACTT => -0.00015, AAAGAA => 0.0001, AAAGAC => 0.00033, AAAGAG => 0.00055, AAAGAT => 0.00024, 
AAAGCA => -6.8e-05, AAAGCC => 0.00026, AAAGCG => -1.2e-05, AAAGCT => 9.1e-05, AAAGGA => 5.6e-05, AAAGGC => 0.00013, 
AAAGGG => -6.1e-05, AAAGGT => 1.6e-05, AAAGTA => -0.00018, AAAGTC => 3.7e-05, AAAGTG => 0.00014, AAAGTT => -0.00016, 
AAATAA => -0.0011, AAATAC => -7.5e-05, AAATAG => -0.00036, AAATAT => -0.00058, AAATCA => -0.00027, AAATCC => 2.3e-05, 
AAATCG => -5.9e-05, AAATCT => -0.00014, AAATGA => -0.00063, AAATGC => -0.00019, AAATGG => -0.00017, AAATGT => -0.00049, 
AAATTA => -0.00058, AAATTC => -0.00011, AAATTG => -0.00024, AAATTT => -0.00048, AACAAA => -0.00027, AACAAC => 0.00045, 
AACAAG => 0.00055, AACAAT => 7e-05, AACACA => -7.9e-05, AACACC => 0.00025, AACACG => 0.00011, AACACT => -1.9e-05, 
AACAGA => -0.00015, AACAGC => 0.00036, AACAGG => -1.7e-05, AACAGT => 4.9e-05, AACATA => -0.00012, AACATC => 0.00051, 
AACATG => 0.0003, AACATT => -7.4e-05, AACCAA => -0.00011, AACCAC => 0.00012, AACCAG => 0.00043, AACCAT => -3.7e-05, 
AACCCA => 3.5e-05, AACCCC => 0.00019, AACCCG => 6.4e-05, AACCCT => 7.4e-05, AACCGA => 3.6e-05, AACCGC => 0.00019, 
AACCGG => 0.00011, AACCGT => 4.1e-05, AACCTA => -1.3e-05, AACCTC => 0.00022, AACCTG => 0.00059, AACCTT => -4.5e-06, 
AACGAA => 6e-05, AACGAC => 0.00025, AACGAG => 0.00044, AACGAT => 0.00019, AACGCA => -9.1e-07, AACGCC => 0.00026, 
AACGCG => 2.6e-05, AACGCT => 6.2e-05, AACGGA => 0.00016, AACGGC => 0.00026, AACGGG => 9.7e-05, AACGGT => 6.7e-05, 
AACGTA => -3.6e-05, AACGTC => 0.00011, AACGTG => 0.00025, AACGTT => -3.9e-05, AACTAA => -0.00034, AACTAC => 0.00032, 
AACTAG => -0.00016, AACTAT => 4.1e-05, AACTCA => -5.8e-05, AACTCC => 0.00021, AACTCG => 0.00013, AACTCT => 2.8e-06, 
AACTGA => -0.00039, AACTGC => 9.5e-05, AACTGG => 9.3e-05, AACTGT => -9e-05, AACTTA => -0.00019, AACTTC => 0.00033, 
AACTTG => 1.2e-05, AACTTT => -4.7e-05, AAGAAA => 0.00037, AAGAAC => 0.00056, AAGAAG => 0.0015, AAGAAT => 0.00014, 
AAGACA => 1.9e-05, AAGACC => 0.00042, AAGACG => 0.00021, AAGACT => 0.0001, AAGAGA => -2.9e-05, AAGAGC => 0.00033, 
AAGAGG => 0.00015, AAGAGT => 8.2e-05, AAGATA => -6e-05, AAGATC => 0.00063, AAGATG => 0.00044, AAGATT => 0.00016, 
AAGCAA => -3.9e-05, AAGCAC => 0.00022, AAGCAG => 0.0007, AAGCAT => -8.8e-06, AAGCCA => 0.00014, AAGCCC => 0.00036, 
AAGCCG => 0.0002, AAGCCT => 0.00015, AAGCGA => 0.0001, AAGCGC => 0.00031, AAGCGG => 0.00023, AAGCGT => 0.00012, 
AAGCTA => 2.7e-05, AAGCTC => 0.00026, AAGCTG => 0.00093, AAGCTT => 5.8e-05, AAGGAA => 0.00049, AAGGAC => 0.00071, 
AAGGAG => 0.0014, AAGGAT => 0.00058, AAGGCA => 0.00014, AAGGCC => 0.00059, AAGGCG => 0.00021, AAGGCT => 0.00032, 
AAGGGA => 0.0001, AAGGGC => 0.0004, AAGGGG => 3.3e-05, AAGGGT => 0.00016, AAGGTA => -4.9e-06, AAGGTC => 0.00027, 
AAGGTG => 0.00065, AAGGTT => 8.6e-05, AAGTAA => -0.00035, AAGTAC => 0.0004, AAGTAG => -0.00019, AAGTAT => 3.9e-05, 
AAGTCA => -3.6e-05, AAGTCC => 0.00025, AAGTCG => 0.00014, AAGTCT => 7.1e-05, AAGTGA => -0.00035, AAGTGC => 9.3e-05, 
AAGTGG => 5.2e-05, AAGTGT => -7.9e-05, AAGTTA => -0.00015, AAGTTC => 0.00043, AAGTTG => 5.4e-05, AAGTTT => 2.7e-05, 
AATAAA => -0.0009, AATAAC => -2.7e-05, AATAAG => 2.2e-05, AATAAT => -0.0004, AATACA => -0.00027, AATACC => 1.3e-05, 
AATACG => -5.4e-06, AATACT => -0.00015, AATAGA => -0.00018, AATAGC => -2.3e-05, AATAGG => -0.00011, AATAGT => -0.00011, 
AATATA => -0.00042, AATATC => -2.7e-05, AATATG => -8.1e-05, AATATT => -0.00047, AATCAA => -0.00019, AATCAC => -4.3e-05, 
AATCAG => 0.00013, AATCAT => -0.00019, AATCCA => 1.1e-05, AATCCC => 0.00012, AATCCG => 8.8e-05, AATCCT => 4.9e-05, 
AATCGA => -3.7e-05, AATCGC => 7.8e-05, AATCGG => 1.8e-05, AATCGT => -1.2e-05, AATCTA => -9.1e-05, AATCTC => 1.7e-05, 
AATCTG => 0.00022, AATCTT => -9e-05, AATGAA => 0.00013, AATGAC => 0.00029, AATGAG => 0.00046, AATGAT => 0.00012, 
AATGCA => 1.9e-05, AATGCC => 0.00041, AATGCG => 7.2e-05, AATGCT => 6.8e-05, AATGGA => 0.00023, AATGGC => 0.00035, 
AATGGG => 8.7e-05, AATGGT => 8.5e-05, AATGTA => -0.00023, AATGTC => 6.6e-05, AATGTG => 0.00024, AATGTT => -0.00024, 
AATTAA => -0.0007, AATTAC => -7.7e-05, AATTAG => -0.00028, AATTAT => -0.0004, AATTCA => -0.00021, AATTCC => -3.8e-05, 
AATTCG => -1e-05, AATTCT => -0.00018, AATTGA => -0.00037, AATTGC => -0.00013, AATTGG => -0.00014, AATTGT => -0.00034, 
AATTTA => -0.0005, AATTTC => -0.00019, AATTTG => -0.00021, AATTTT => -0.00063, ACAAAA => -0.00054, ACAAAC => -0.0002, 
ACAAAG => -7.5e-05, ACAAAT => -0.00029, ACAACA => -0.00022, ACAACC => -9.5e-06, ACAACG => 5.5e-06, ACAACT => -0.0001, 
ACAAGA => -0.00019, ACAAGC => -7.9e-05, ACAAGG => -0.00011, ACAAGT => -9.7e-05, ACAATA => -0.00023, ACAATC => -6.1e-05, 
ACAATG => -6.3e-05, ACAATT => -0.0002, ACACAA => -0.00025, ACACAC => -0.00065, ACACAG => -9.6e-06, ACACAT => -0.00024, 
ACACCA => 2e-05, ACACCC => 5.6e-05, ACACCG => 5.4e-05, ACACCT => -2.4e-05, ACACGA => -3.6e-05, ACACGC => -2e-05, 
ACACGG => -1.2e-05, ACACGT => -4.5e-05, ACACTA => -9.6e-05, ACACTC => -3.7e-05, ACACTG => 6.6e-05, ACACTT => -0.00018, 
ACAGAA => 3.8e-05, ACAGAC => 0.00017, ACAGAG => 0.00031, ACAGAT => 0.00012, ACAGCA => -5.9e-05, ACAGCC => 0.00015, 
ACAGCG => 1.2e-05, ACAGCT => 3.5e-05, ACAGGA => 5.6e-05, ACAGGC => 5.3e-05, ACAGGG => -5.3e-06, ACAGGT => -3.2e-06, 
ACAGTA => -8.9e-05, ACAGTC => 1.3e-05, ACAGTG => 0.00014, ACAGTT => -9.7e-05, ACATAA => -0.00036, ACATAC => -9.6e-05, 
ACATAG => -0.00017, ACATAT => -0.00026, ACATCA => -0.00012, ACATCC => 4.8e-06, ACATCG => 2.3e-05, ACATCT => -5.2e-05, 
ACATGA => -0.0003, ACATGC => -0.00013, ACATGG => -0.0001, ACATGT => -0.0002, ACATTA => -0.00022, ACATTC => -8.3e-05, 
ACATTG => -0.00012, ACATTT => -0.00046, ACCAAA => 0.00012, ACCAAC => 0.00036, ACCAAG => 0.00062, ACCAAT => 0.00019, 
ACCACA => 0.00013, ACCACC => 0.00045, ACCACG => 0.00014, ACCACT => 0.00013, ACCAGA => -9.1e-05, ACCAGC => 0.00031, 
ACCAGG => -5.3e-06, ACCAGT => 0.00013, ACCATA => 6.4e-06, ACCATC => 0.00055, ACCATG => 0.00031, ACCATT => 0.00018, 
ACCCAA => -5.7e-05, ACCCAC => 4.5e-05, ACCCAG => 0.00024, ACCCAT => -2.1e-05, ACCCCA => 1.9e-05, ACCCCC => 2.8e-05, 
ACCCCG => 4.9e-05, ACCCCT => 4e-05, ACCCGA => 1.8e-05, ACCCGC => 9.6e-05, ACCCGG => 5.8e-05, ACCCGT => 2.3e-05, 
ACCCTA => -1.9e-05, ACCCTC => 0.00011, ACCCTG => 0.00039, ACCCTT => -3.5e-05, ACCGAA => 7.4e-05, ACCGAC => 0.00016, 
ACCGAG => 0.00028, ACCGAT => 0.00014, ACCGCA => 2.2e-05, ACCGCC => 0.00021, ACCGCG => 8e-06, ACCGCT => 3.7e-05, 
ACCGGA => 7.9e-05, ACCGGC => 0.00016, ACCGGG => 4.2e-05, ACCGGT => 5.1e-05, ACCGTA => 2.1e-05, ACCGTC => 0.00014, 
ACCGTG => 0.00024, ACCGTT => 2.8e-05, ACCTAA => -0.00017, ACCTAC => 0.0003, ACCTAG => -0.00011, ACCTAT => 9e-05, 
ACCTCA => -1.5e-06, ACCTCC => 0.00014, ACCTCG => 0.0001, ACCTCT => -1.8e-05, ACCTGA => -0.00026, ACCTGC => 0.00014, 
ACCTGG => 1.2e-05, ACCTGT => -2.4e-07, ACCTTA => -7.5e-05, ACCTTC => 0.00038, ACCTTG => 3.3e-05, ACCTTT => 3.4e-06, 
ACGAAA => -9.9e-05, ACGAAC => 1.4e-05, ACGAAG => 3.5e-05, ACGAAT => -2.7e-05, ACGACA => -1.7e-05, ACGACC => 1.9e-05, 
ACGACG => 3e-05, ACGACT => -1.4e-05, ACGAGA => -6e-05, ACGAGC => 1.5e-05, ACGAGG => -4.4e-05, ACGAGT => -3.1e-05, 
ACGATA => -3e-05, ACGATC => 1.8e-05, ACGATG => 2.5e-05, ACGATT => -3.5e-05, ACGCAA => -2.3e-05, ACGCAC => 1.6e-05, 
ACGCAG => 0.00012, ACGCAT => -2.8e-05, ACGCCA => 7.3e-05, ACGCCC => 0.00019, ACGCCG => 9.9e-05, ACGCCT => 2.9e-05, 
ACGCGA => -2.6e-06, ACGCGC => 4.1e-05, ACGCGG => 2.2e-05, ACGCGT => -1.1e-06, ACGCTA => 1.2e-05, ACGCTC => 0.0001, 
ACGCTG => 0.00032, ACGCTT => -9.7e-07, ACGGAA => 7.1e-05, ACGGAC => 0.00018, ACGGAG => 0.00033, ACGGAT => 0.00014, 
ACGGCA => 4.7e-05, ACGGCC => 0.00022, ACGGCG => 9.7e-05, ACGGCT => 8.7e-05, ACGGGA => 7.1e-05, ACGGGC => 0.00018, 
ACGGGG => 2.4e-05, ACGGGT => 4.6e-05, ACGGTA => -1.4e-05, ACGGTC => 7.4e-05, ACGGTG => 0.00025, ACGGTT => 2.8e-05, 
ACGTAA => -9.1e-05, ACGTAC => 6.1e-05, ACGTAG => -5.8e-05, ACGTAT => -1.5e-05, ACGTCA => -1.9e-05, ACGTCC => 5.3e-05, 
ACGTCG => 3e-05, ACGTCT => -7.2e-06, ACGTGA => -0.00011, ACGTGC => 2.1e-06, ACGTGG => -2e-05, ACGTGT => -6.1e-05, 
ACGTTA => -5.2e-05, ACGTTC => 4.4e-05, ACGTTG => -2.9e-06, ACGTTT => -5.7e-05, ACTAAA => -0.00014, ACTAAC => -4e-05, 
ACTAAG => 1.7e-06, ACTAAT => -0.0001, ACTACA => -8.3e-05, ACTACC => 1.2e-05, ACTACG => -4.8e-07, ACTACT => -6.1e-05, 
ACTAGA => -0.0001, ACTAGC => -1.9e-05, ACTAGG => -6.9e-05, ACTAGT => -4.3e-05, ACTATA => -0.00013, ACTATC => 4.3e-06, 
ACTATG => -3.4e-05, ACTATT => -0.0001, ACTCAA => -0.0001, ACTCAC => -4.5e-05, ACTCAG => 0.0001, ACTCAT => -7.8e-05, 
ACTCCA => 2.9e-05, ACTCCC => 2.4e-05, ACTCCG => 6e-05, ACTCCT => 4e-05, ACTCGA => -1.1e-05, ACTCGC => 1e-05, 
ACTCGG => 1.5e-05, ACTCGT => -5.8e-06, ACTCTA => -5.1e-05, ACTCTC => -3.3e-05, ACTCTG => 0.00015, ACTCTT => -9.6e-05, 
ACTGAA => 2.9e-05, ACTGAC => 9.6e-05, ACTGAG => 0.00017, ACTGAT => 8.8e-05, ACTGCA => -6.4e-05, ACTGCC => 0.00015, 
ACTGCG => -3.2e-06, ACTGCT => 3.5e-05, ACTGGA => 0.00013, ACTGGC => 0.00011, ACTGGG => 4.6e-06, ACTGGT => 4.8e-05, 
ACTGTA => -0.00013, ACTGTC => 5.2e-05, ACTGTG => 0.00017, ACTGTT => -8.3e-05, ACTTAA => -0.00037, ACTTAC => -7.1e-06, 
ACTTAG => -0.00017, ACTTAT => -0.00013, ACTTCA => -0.00013, ACTTCC => -4.9e-05, ACTTCG => -1.5e-05, ACTTCT => -0.00012, 
ACTTGA => -0.0003, ACTTGC => -8.9e-05, ACTTGG => -0.00016, ACTTGT => -0.00018, ACTTTA => -0.00024, ACTTTC => -6.6e-05, 
ACTTTG => -0.00014, ACTTTT => -0.00041, AGAAAA => -0.00058, AGAAAC => -0.00018, AGAAAG => -0.00018, AGAAAT => -0.00033, 
AGAACA => -0.00024, AGAACC => -6.4e-05, AGAACG => -6.1e-05, AGAACT => -0.00015, AGAAGA => -0.00033, AGAAGC => -0.00018, 
AGAAGG => -0.00022, AGAAGT => -0.00017, AGAATA => -0.00023, AGAATC => -8.8e-05, AGAATG => -0.00015, AGAATT => -0.00021, 
AGACAA => -0.00017, AGACAC => -0.00013, AGACAG => -5.8e-05, AGACAT => -0.00015, AGACCA => -0.00011, AGACCC => -2.9e-05, 
AGACCG => -2.4e-05, AGACCT => -5.3e-05, AGACGA => -3.9e-05, AGACGC => -8.3e-06, AGACGG => -3.6e-05, AGACGT => -3.7e-05, 
AGACTA => -8.1e-05, AGACTC => -6.8e-05, AGACTG => -5.6e-05, AGACTT => -0.00016, AGAGAA => -0.00014, AGAGAC => 3.8e-05, 
AGAGAG => -7.5e-05, AGAGAT => 1.1e-06, AGAGCA => -0.00017, AGAGCC => -2.2e-05, AGAGCG => -8.4e-05, AGAGCT => -0.0001, 
AGAGGA => -0.00013, AGAGGC => -9.6e-05, AGAGGG => -0.00018, AGAGGT => -9.2e-05, AGAGTA => -0.00011, AGAGTC => -6.2e-05, 
AGAGTG => -5e-05, AGAGTT => -0.00014, AGATAA => -0.0003, AGATAC => -3.4e-05, AGATAG => -0.00016, AGATAT => -0.00018, 
AGATCA => -0.00015, AGATCC => -7.5e-05, AGATCG => -4.9e-05, AGATCT => -0.00013, AGATGA => -0.00036, AGATGC => -0.00019, 
AGATGG => -0.00025, AGATGT => -0.00024, AGATTA => -0.00018, AGATTC => -0.00012, AGATTG => -0.00015, AGATTT => -0.00032, 
AGCAAA => 2.1e-05, AGCAAC => 0.0003, AGCAAG => 0.0004, AGCAAT => 0.00012, AGCACA => 4.3e-05, AGCACC => 0.00033, 
AGCACG => 7.8e-05, AGCACT => 4.7e-05, AGCAGA => -0.00019, AGCAGC => 0.00052, AGCAGG => -0.00011, AGCAGT => 0.00018, 
AGCATA => -9e-05, AGCATC => 0.00031, AGCATG => 0.0002, AGCATT => -1.9e-05, AGCCAA => -7.7e-05, AGCCAC => 3.2e-05, 
AGCCAG => 0.00031, AGCCAT => -0.0001, AGCCCA => -7e-06, AGCCCC => 0.00012, AGCCCG => 2.3e-05, AGCCCT => 2.7e-05, 
AGCCGA => 1.3e-05, AGCCGC => 0.00012, AGCCGG => 7.1e-05, AGCCGT => 8.5e-06, AGCCTA => -4.1e-05, AGCCTC => 7.2e-05, 
AGCCTG => 0.00032, AGCCTT => -5.6e-05, AGCGAA => 7.6e-05, AGCGAC => 0.00026, AGCGAG => 0.00033, AGCGAT => 0.00017, 
AGCGCA => -3.3e-06, AGCGCC => 0.00019, AGCGCG => -2e-06, AGCGCT => 2.3e-05, AGCGGA => 6.9e-05, AGCGGC => 0.00023, 
AGCGGG => 4.1e-05, AGCGGT => 4.9e-05, AGCGTA => -1.5e-05, AGCGTC => 6.7e-05, AGCGTG => 0.00015, AGCGTT => 3.3e-06, 
AGCTAA => -0.00024, AGCTAC => 0.00019, AGCTAG => -0.00016, AGCTAT => 3.2e-05, AGCTCA => -4.7e-05, AGCTCC => 0.00024, 
AGCTCG => 7.1e-05, AGCTCT => 7.1e-06, AGCTGA => -0.00039, AGCTGC => -3.6e-05, AGCTGG => -0.00012, AGCTGT => -0.00015, 
AGCTTA => -0.00011, AGCTTC => 0.00021, AGCTTG => -3.8e-06, AGCTTT => -6.5e-05, AGGAAA => -0.00016, AGGAAC => 5.5e-05, 
AGGAAG => 9.6e-05, AGGAAT => -8.1e-05, AGGACA => -0.00017, AGGACC => -1.7e-05, AGGACG => -2.4e-05, AGGACT => -0.0001, 
AGGAGA => -0.00022, AGGAGC => -9.3e-05, AGGAGG => -0.00023, AGGAGT => -0.00011, AGGATA => -8.8e-05, AGGATC => 2.1e-05, 
AGGATG => -4.5e-05, AGGATT => -9.2e-05, AGGCAA => -0.00013, AGGCAC => -5.5e-05, AGGCAG => -6.4e-05, AGGCAT => -0.00012, 
AGGCCA => -0.00012, AGGCCC => -3.9e-05, AGGCCG => -2.1e-05, AGGCCT => -8.4e-05, AGGCGA => -3.5e-05, AGGCGC => 1.4e-05, 
AGGCGG => -3.3e-05, AGGCGT => -3.7e-05, AGGCTA => -8.7e-05, AGGCTC => -5.5e-05, AGGCTG => -5e-05, AGGCTT => -0.00014, 
AGGGAA => -5.7e-05, AGGGAC => 6.2e-05, AGGGAG => 8.8e-05, AGGGAT => 3.8e-05, AGGGCA => -0.00012, AGGGCC => 8.1e-06, 
AGGGCG => -1.7e-05, AGGGCT => -8.4e-05, AGGGGA => -0.00014, AGGGGC => -4e-05, AGGGGG => -0.00013, AGGGGT => -7.4e-05, 
AGGGTA => -7e-05, AGGGTC => -3.1e-05, AGGGTG => -2.4e-05, AGGGTT => -0.00012, AGGTAA => -0.00018, AGGTAC => -9.5e-06, 
AGGTAG => -0.00016, AGGTAT => -6.1e-05, AGGTCA => -0.00014, AGGTCC => -4.6e-05, AGGTCG => -2.7e-05, AGGTCT => -9.6e-05, 
AGGTGA => -0.00025, AGGTGC => -0.00011, AGGTGG => -0.00019, AGGTGT => -0.00018, AGGTTA => -0.0001, AGGTTC => -3.3e-05, 
AGGTTG => -0.00011, AGGTTT => -0.00014, AGTAAA => -0.00015, AGTAAC => 7.3e-06, AGTAAG => -1.4e-05, AGTAAT => -8.3e-05, 
AGTACA => -8.8e-05, AGTACC => 8e-06, AGTACG => -1.3e-05, AGTACT => -7.4e-05, AGTAGA => -0.00013, AGTAGC => 1.7e-06, 
AGTAGG => -9.5e-05, AGTAGT => -4.5e-05, AGTATA => -0.00014, AGTATC => -2.5e-05, AGTATG => -4.1e-05, AGTATT => -0.00018, 
AGTCAA => -0.00011, AGTCAC => -6.1e-05, AGTCAG => 8.1e-05, AGTCAT => -0.0001, AGTCCA => 3.7e-06, AGTCCC => 7.3e-05, 
AGTCCG => 6.1e-05, AGTCCT => 3.8e-05, AGTCGA => 5.2e-07, AGTCGC => 6e-05, AGTCGG => 1.5e-06, AGTCGT => 1.4e-05, 
AGTCTA => -6.8e-05, AGTCTC => -5.7e-05, AGTCTG => 0.0001, AGTCTT => -0.0001, AGTGAA => 8.2e-05, AGTGAC => 0.00024, 
AGTGAG => 0.00024, AGTGAT => 0.00017, AGTGCA => -3e-05, AGTGCC => 0.00026, AGTGCG => 6.4e-06, AGTGCT => 1.2e-05, 
AGTGGA => 8.9e-05, AGTGGC => 0.00019, AGTGGG => 1.2e-05, AGTGGT => 3.5e-05, AGTGTA => -0.00012, AGTGTC => 2.9e-05, 
AGTGTG => 6.2e-05, AGTGTT => -0.00017, AGTTAA => -0.00032, AGTTAC => -1.9e-05, AGTTAG => -0.00016, AGTTAT => -0.00015, 
AGTTCA => -0.0001, AGTTCC => -1.3e-05, AGTTCG => 6.9e-06, AGTTCT => -9e-05, AGTTGA => -0.00029, AGTTGC => -0.00013, 
AGTTGG => -0.00015, AGTTGT => -0.0002, AGTTTA => -0.00024, AGTTTC => -0.00011, AGTTTG => -0.00014, AGTTTT => -0.00042, 
ATAAAA => -0.00079, ATAAAC => -0.00024, ATAAAG => -0.00019, ATAAAT => -0.00063, ATAACA => -0.00021, ATAACC => -3.9e-05, 
ATAACG => -2.7e-05, ATAACT => -0.00015, ATAAGA => -0.00019, ATAAGC => -8.4e-05, ATAAGG => -0.00011, ATAAGT => -0.00016, 
ATAATA => -0.0005, ATAATC => -0.0001, ATAATG => -0.00018, ATAATT => -0.0004, ATACAA => -0.00027, ATACAC => -0.0002, 
ATACAG => -4.3e-05, ATACAT => -0.00036, ATACCA => -5.6e-05, ATACCC => 2.7e-05, ATACCG => 1e-05, ATACCT => -6.4e-05, 
ATACGA => -3.6e-05, ATACGC => 2.8e-05, ATACGG => -1.6e-05, ATACGT => -5e-05, ATACTA => -0.00011, ATACTC => -3.9e-05, 
ATACTG => -4.1e-05, ATACTT => -0.0002, ATAGAA => -5e-05, ATAGAC => 5.1e-05, ATAGAG => 5.8e-05, ATAGAT => -1.6e-05, 
ATAGCA => -6e-05, ATAGCC => 5.7e-05, ATAGCG => 1e-05, ATAGCT => -5.6e-05, ATAGGA => -5.6e-05, ATAGGC => -1.7e-05, 
ATAGGG => -4.9e-05, ATAGGT => -7.2e-05, ATAGTA => -0.00012, ATAGTC => -3.7e-05, ATAGTG => 1.1e-05, ATAGTT => -0.00016, 
ATATAA => -0.00055, ATATAC => -0.00023, ATATAG => -0.00025, ATATAT => -0.00081, ATATCA => -0.00016, ATATCC => -4.2e-05, 
ATATCG => -1.6e-05, ATATCT => -0.00018, ATATGA => -0.00034, ATATGC => -0.00015, ATATGG => -0.00013, ATATGT => -0.00034, 
ATATTA => -0.00039, ATATTC => -0.00015, ATATTG => -0.00025, ATATTT => -0.00073, ATCAAA => 6.4e-05, ATCAAC => 0.00048, 
ATCAAG => 0.00064, ATCAAT => 0.00016, ATCACA => 0.00011, ATCACC => 0.00041, ATCACG => 0.00012, ATCACT => 0.00012, 
ATCAGA => -0.0001, ATCAGC => 0.00031, ATCAGG => -3.8e-06, ATCAGT => 7.8e-05, ATCATA => -5.2e-05, ATCATC => 0.00051, 
ATCATG => 0.00031, ATCATT => 8.2e-05, ATCCAA => -2.6e-05, ATCCAC => 0.0002, ATCCAG => 0.00054, ATCCAT => -2.1e-05, 
ATCCCA => 2.1e-05, ATCCCC => 0.00012, ATCCCG => 5.4e-05, ATCCCT => 5.8e-05, ATCCGA => 8.5e-05, ATCCGC => 0.00022, 
ATCCGG => 0.00016, ATCCGT => 7.9e-05, ATCCTA => 1.2e-06, ATCCTC => 0.00023, ATCCTG => 0.00064, ATCCTT => -8.5e-07, 
ATCGAA => 0.00012, ATCGAC => 0.00032, ATCGAG => 0.00055, ATCGAT => 0.00025, ATCGCA => 4.6e-05, ATCGCC => 0.00032, 
ATCGCG => 3.1e-05, ATCGCT => 0.00011, ATCGGA => 8.8e-05, ATCGGC => 0.00021, ATCGGG => 8e-05, ATCGGT => 5.6e-05, 
ATCGTA => 2.2e-06, ATCGTC => 0.00017, ATCGTG => 0.00029, ATCGTT => 1.7e-05, ATCTAA => -0.00024, ATCTAC => 0.00033, 
ATCTAG => -0.00013, ATCTAT => 2.7e-05, ATCTCA => -4.7e-05, ATCTCC => 0.00021, ATCTCG => 9.3e-05, ATCTCT => -9.9e-06, 
ATCTGA => -0.0003, ATCTGC => 0.00012, ATCTGG => 3.8e-05, ATCTGT => -8.2e-05, ATCTTA => -0.00013, ATCTTC => 0.00036, 
ATCTTG => 1.1e-05, ATCTTT => -5.8e-05, ATGAAA => -9.3e-05, ATGAAC => 0.00031, ATGAAG => 0.00054, ATGAAT => -1.9e-05, 
ATGACA => -1.2e-05, ATGACC => 0.00027, ATGACG => 8.4e-05, ATGACT => 1e-05, ATGAGA => -7.7e-05, ATGAGC => 0.0002, 
ATGAGG => 4.5e-05, ATGAGT => 3.6e-06, ATGATA => -0.00013, ATGATC => 0.00027, ATGATG => 0.00025, ATGATT => -6.6e-05, 
ATGCAA => -9.1e-05, ATGCAC => 9.1e-05, ATGCAG => 0.00045, ATGCAT => -0.00013, ATGCCA => 3.3e-05, ATGCCC => 0.00021, 
ATGCCG => 0.00012, ATGCCT => 6e-05, ATGCGA => 4.1e-05, ATGCGC => 0.00016, ATGCGG => 0.0001, ATGCGT => 4.3e-05, 
ATGCTA => -3.7e-05, ATGCTC => 0.00012, ATGCTG => 0.00054, ATGCTT => -7.9e-05, ATGGAA => 0.00025, ATGGAC => 0.00054, 
ATGGAG => 0.00091, ATGGAT => 0.00035, ATGGCA => 0.00016, ATGGCC => 0.00052, ATGGCG => 0.00022, ATGGCT => 0.00026, 
ATGGGA => 0.0001, ATGGGC => 0.00035, ATGGGG => 1.7e-05, ATGGGT => 0.00011, ATGGTA => -3.1e-05, ATGGTC => 0.00015, 
ATGGTG => 0.00042, ATGGTT => -1.2e-05, ATGTAA => -0.00039, ATGTAC => 0.00018, ATGTAG => -0.00022, ATGTAT => -0.00023, 
ATGTCA => -5.6e-05, ATGTCC => 0.0002, ATGTCG => 0.00011, ATGTCT => 6.1e-06, ATGTGA => -0.00036, ATGTGC => -4.4e-05, 
ATGTGG => -3.9e-05, ATGTGT => -0.00027, ATGTTA => -0.00022, ATGTTC => 0.00011, ATGTTG => -4.5e-05, ATGTTT => -0.00025, 
ATTAAA => -0.00043, ATTAAC => -6.3e-05, ATTAAG => 5.3e-05, ATTAAT => -0.0003, ATTACA => -0.00018, ATTACC => -2.7e-05, 
ATTACG => -1.4e-05, ATTACT => -0.00013, ATTAGA => -0.00017, ATTAGC => -6.7e-05, ATTAGG => -8.7e-05, ATTAGT => -0.00014, 
ATTATA => -0.00039, ATTATC => -5.7e-06, ATTATG => -0.00011, ATTATT => -0.00048, ATTCAA => -0.00016, ATTCAC => -5.4e-06, 
ATTCAG => 0.00018, ATTCAT => -0.00013, ATTCCA => -2.8e-05, ATTCCC => 3.5e-05, ATTCCG => 4e-05, ATTCCT => 2.5e-05, 
ATTCGA => 1.5e-05, ATTCGC => 8.8e-05, ATTCGG => 3.7e-05, ATTCGT => 2.5e-05, ATTCTA => -0.00013, ATTCTC => -4.3e-05, 
ATTCTG => 0.00013, ATTCTT => -0.00016, ATTGAA => 0.00017, ATTGAC => 0.00029, ATTGAG => 0.00044, ATTGAT => 0.00019, 
ATTGCA => 1e-05, ATTGCC => 0.0003, ATTGCG => 4e-05, ATTGCT => 0.00011, ATTGGA => 0.00014, ATTGGC => 0.0002, 
ATTGGG => 3e-05, ATTGGT => 5.5e-05, ATTGTA => -0.0002, ATTGTC => 9e-05, ATTGTG => 0.00022, ATTGTT => -0.00019, 
ATTTAA => -0.00076, ATTTAC => -0.00015, ATTTAG => -0.00033, ATTTAT => -0.00061, ATTTCA => -0.00032, ATTTCC => -0.00014, 
ATTTCG => -5.2e-05, ATTTCT => -0.00033, ATTTGA => -0.00047, ATTTGC => -0.00021, ATTTGG => -0.00022, ATTTGT => -0.00046, 
ATTTTA => -0.00078, ATTTTC => -0.00036, ATTTTG => -0.00038, ATTTTT => -0.0011, CAAAAA => -0.00047, CAAAAC => -0.00017, 
CAAAAG => -7.5e-06, CAAAAT => -0.00029, CAAACA => -0.00032, CAAACC => -6.4e-05, CAAACG => -3.5e-05, CAAACT => -0.00018, 
CAAAGA => -0.00023, CAAAGC => -0.00012, CAAAGG => -0.00013, CAAAGT => -0.00016, CAAATA => -0.00034, CAAATC => -9.9e-05, 
CAAATG => -0.0001, CAAATT => -0.00027, CAACAA => -0.00013, CAACAC => -8.4e-05, CAACAG => 0.00033, CAACAT => -0.00015, 
CAACCA => -6.4e-05, CAACCC => -2.1e-05, CAACCG => 3e-05, CAACCT => -5.8e-05, CAACGA => -7e-05, CAACGC => 3.8e-06, 
CAACGG => -3.3e-06, CAACGT => -2.6e-05, CAACTA => -6.5e-05, CAACTC => -3.6e-05, CAACTG => 0.00011, CAACTT => -0.00017, 
CAAGAA => 9.5e-05, CAAGAC => 9.7e-05, CAAGAG => 0.00023, CAAGAT => 8.1e-05, CAAGCA => -4.5e-05, CAAGCC => 7.2e-05, 
CAAGCG => 3.2e-06, CAAGCT => 2.2e-05, CAAGGA => 3.8e-05, CAAGGC => 4.8e-05, CAAGGG => -3e-05, CAAGGT => 4.5e-06, 
CAAGTA => -8e-05, CAAGTC => -3.5e-05, CAAGTG => 8.2e-05, CAAGTT => -8.9e-05, CAATAA => -0.00043, CAATAC => -3.3e-05, 
CAATAG => -0.00017, CAATAT => -0.00015, CAATCA => -0.00015, CAATCC => -4.3e-05, CAATCG => 7.7e-06, CAATCT => -8.8e-05, 
CAATGA => -0.00026, CAATGC => -0.00011, CAATGG => -0.00012, CAATGT => -0.00016, CAATTA => -0.00024, CAATTC => -0.00011, 
CAATTG => -8.2e-05, CAATTT => -0.00032, CACAAA => -0.00015, CACAAC => 0.00017, CACAAG => 0.00031, CACAAT => 1.5e-05, 
CACACA => -0.00066, CACACC => 6.3e-05, CACACG => 7.6e-05, CACACT => -6.1e-05, CACAGA => -0.00021, CACAGC => 7.7e-05, 
CACAGG => -9.1e-05, CACAGT => -1.8e-05, CACATA => -0.00013, CACATC => 0.00023, CACATG => 0.00016, CACATT => -0.0001, 
CACCAA => -0.00011, CACCAC => 4e-05, CACCAG => 0.00027, CACCAT => -4.5e-05, CACCCA => -0.00011, CACCCC => -3.5e-05, 
CACCCG => 2.9e-05, CACCCT => -4.5e-05, CACCGA => -1.9e-05, CACCGC => 4.3e-05, CACCGG => 7e-05, CACCGT => -9.7e-06, 
CACCTA => -5.7e-05, CACCTC => -1.3e-05, CACCTG => 0.00032, CACCTT => -9.4e-05, CACGAA => 4.9e-05, CACGAC => 0.00016, 
CACGAG => 0.00028, CACGAT => 0.0001, CACGCA => -4.3e-05, CACGCC => 0.00011, CACGCG => 5.2e-06, CACGCT => -2e-05, 
CACGGA => 5.2e-05, CACGGC => 0.00013, CACGGG => 3.5e-05, CACGGT => -5e-06, CACGTA => -2e-05, CACGTC => 5.9e-05, 
CACGTG => 0.00013, CACGTT => -2.9e-05, CACTAA => -0.0002, CACTAC => 0.00019, CACTAG => -0.00013, CACTAT => 2.6e-05, 
CACTCA => -0.00012, CACTCC => -1.5e-05, CACTCG => 7.8e-05, CACTCT => -7.5e-05, CACTGA => -0.00034, CACTGC => -3.8e-05, 
CACTGG => -7.1e-05, CACTGT => -0.00013, CACTTA => -0.00015, CACTTC => 0.00013, CACTTG => -2.5e-05, CACTTT => -0.00019, 
CAGAAA => 0.00016, CAGAAC => 0.00045, CAGAAG => 0.00086, CAGAAT => 0.00021, CAGACA => 1.6e-05, CAGACC => 0.0003, 
CAGACG => 0.00019, CAGACT => 5.9e-05, CAGAGA => -2.3e-06, CAGAGC => 0.00029, CAGAGG => 8.3e-05, CAGAGT => 0.00012, 
CAGATA => -1.9e-05, CAGATC => 0.00047, CAGATG => 0.00041, CAGATT => 0.00013, CAGCAA => 0.00033, CAGCAC => 0.0003, 
CAGCAG => 0.0017, CAGCAT => 9.3e-05, CAGCCA => 0.00014, CAGCCC => 0.00027, CAGCCG => 0.00025, CAGCCT => 7.4e-05, 
CAGCGA => 0.00011, CAGCGC => 0.00028, CAGCGG => 0.00022, CAGCGT => 8.5e-05, CAGCTA => 4e-05, CAGCTC => 0.00029, 
CAGCTG => 0.0011, CAGCTT => 6.7e-05, CAGGAA => 0.00039, CAGGAC => 0.00056, CAGGAG => 0.0012, CAGGAT => 0.00042, 
CAGGCA => 0.00017, CAGGCC => 0.0006, CAGGCG => 0.00029, CAGGCT => 0.00025, CAGGGA => 7.5e-05, CAGGGC => 0.00037, 
CAGGGG => 1.4e-05, CAGGGT => 0.00014, CAGGTA => -1.1e-05, CAGGTC => 0.00019, CAGGTG => 0.00064, CAGGTT => 5e-05, 
CAGTAA => -0.00025, CAGTAC => 0.00038, CAGTAG => -0.00016, CAGTAT => 0.00015, CAGTCA => -1.1e-05, CAGTCC => 0.00017, 
CAGTCG => 0.00014, CAGTCT => 8.1e-05, CAGTGA => -0.00035, CAGTGC => 9.2e-05, CAGTGG => 3.8e-05, CAGTGT => -3.2e-05, 
CAGTTA => -6.1e-05, CAGTTC => 0.00035, CAGTTG => 0.00012, CAGTTT => 6.7e-05, CATAAA => -0.00024, CATAAC => -6.6e-05, 
CATAAG => 5e-06, CATAAT => -0.00017, CATACA => -0.00021, CATACC => -6.3e-05, CATACG => -1.9e-05, CATACT => -5e-05, 
CATAGA => -0.00012, CATAGC => -5.8e-05, CATAGG => -7.5e-05, CATAGT => -9.3e-05, CATATA => -0.00027, CATATC => -5.7e-05, 
CATATG => -0.00014, CATATT => -0.00027, CATCAA => -0.00013, CATCAC => -2.6e-05, CATCAG => 0.00012, CATCAT => -0.00014, 
CATCCA => -7.4e-05, CATCCC => -1.7e-05, CATCCG => 7.9e-05, CATCCT => -8.3e-05, CATCGA => -1.8e-05, CATCGC => 5.7e-05, 
CATCGG => 7e-06, CATCGT => -9.7e-06, CATCTA => -9.8e-05, CATCTC => -9.8e-05, CATCTG => 4.9e-05, CATCTT => -0.00019, 
CATGAA => 3.3e-05, CATGAC => 8.4e-05, CATGAG => 0.0002, CATGAT => 3.1e-05, CATGCA => -8.4e-05, CATGCC => 0.0001, 
CATGCG => 1e-05, CATGCT => -1.7e-05, CATGGA => 1.4e-05, CATGGC => 8.7e-05, CATGGG => -2.6e-05, CATGGT => -1.4e-05, 
CATGTA => -0.00014, CATGTC => -1.7e-05, CATGTG => 6.6e-05, CATGTT => -0.00016, CATTAA => -0.00035, CATTAC => -3.4e-05, 
CATTAG => -0.00019, CATTAT => -0.00019, CATTCA => -0.00022, CATTCC => -9.7e-05, CATTCG => -2.1e-05, CATTCT => -0.00017, 
CATTGA => -0.00027, CATTGC => -0.00016, CATTGG => -0.00015, CATTGT => -0.00024, CATTTA => -0.00039, CATTTC => -0.00024, 
CATTTG => -0.00024, CATTTT => -0.00066, CCAAAA => -0.00017, CCAAAC => -2.4e-05, CCAAAG => 7.9e-05, CCAAAT => -9.5e-05, 
CCAACA => -2.8e-05, CCAACC => -1.6e-05, CCAACG => 7.6e-06, CCAACT => -3e-05, CCAAGA => -0.00013, CCAAGC => -3.5e-05, 
CCAAGG => -0.00012, CCAAGT => -6.3e-05, CCAATA => -0.00011, CCAATC => -3e-05, CCAATG => 5.6e-05, CCAATT => -8.6e-05, 
CCACAA => -7.8e-05, CCACAC => -0.00012, CCACAG => 0.00017, CCACAT => -9.5e-05, CCACCA => 0.00023, CCACCC => 0.00011, 
CCACCG => 0.00015, CCACCT => 0.00016, CCACGA => 1.4e-05, CCACGC => -9.4e-06, CCACGG => -4.7e-06, CCACGT => -9.1e-07, 
CCACTA => -5.8e-05, CCACTC => -3.9e-05, CCACTG => 0.00012, CCACTT => -0.00011, CCAGAA => 0.0002, CCAGAC => 0.0002, 
CCAGAG => 0.00045, CCAGAT => 0.00022, CCAGCA => -2.2e-05, CCAGCC => 6.9e-05, CCAGCG => 3.3e-05, CCAGCT => 3.1e-05, 
CCAGGA => 0.00011, CCAGGC => 3.9e-05, CCAGGG => -1.4e-05, CCAGGT => 5e-05, CCAGTA => -2.3e-05, CCAGTC => 2.1e-05, 
CCAGTG => 0.00018, CCAGTT => -2.6e-05, CCATAA => -0.00021, CCATAC => 1.6e-05, CCATAG => -0.00015, CCATAT => -5.2e-05, 
CCATCA => -5.3e-05, CCATCC => 1.5e-05, CCATCG => 6.3e-05, CCATCT => -5.2e-05, CCATGA => -0.00022, CCATGC => -0.00011, 
CCATGG => -0.00011, CCATGT => -0.00016, CCATTA => -0.00011, CCATTC => -2.8e-05, CCATTG => -5.4e-05, CCATTT => -0.00021, 
CCCAAA => 0.0002, CCCAAC => 0.00037, CCCAAG => 0.00055, CCCAAT => 0.00021, CCCACA => 0.00011, CCCACC => 0.00022, 
CCCACG => 0.00014, CCCACT => 5.7e-05, CCCAGA => -9.8e-05, CCCAGC => 0.00025, CCCAGG => -0.00011, CCCAGT => 0.00015, 
CCCATA => -2.5e-05, CCCATC => 0.00033, CCCATG => 0.00029, CCCATT => 0.00011, CCCCAA => -8.1e-05, CCCCAC => -7.9e-05, 
CCCCAG => 0.00021, CCCCAT => -5.8e-05, CCCCCA => 5.9e-05, CCCCCC => -5.4e-05, CCCCCG => 6.5e-05, CCCCCT => 4.7e-05, 
CCCCGA => 1.2e-05, CCCCGC => -7.5e-06, CCCCGG => 4.2e-05, CCCCGT => -5.3e-06, CCCCTA => -4.3e-05, CCCCTC => -7.5e-05, 
CCCCTG => 0.00021, CCCCTT => -0.00011, CCCGAA => 0.0001, CCCGAC => 0.00018, CCCGAG => 0.00034, CCCGAT => 0.00014, 
CCCGCA => 1.1e-05, CCCGCC => 0.00011, CCCGCG => 2.5e-05, CCCGCT => 3.9e-05, CCCGGA => 6.9e-05, CCCGGC => 0.00012, 
CCCGGG => -5.2e-06, CCCGGT => 4.2e-05, CCCGTA => 8.2e-06, CCCGTC => 8.5e-05, CCCGTG => 0.00019, CCCGTT => 3.4e-05, 
CCCTAA => -0.00015, CCCTAC => 0.00025, CCCTAG => -0.00012, CCCTAT => 0.00011, CCCTCA => -2.5e-05, CCCTCC => -3.4e-05, 
CCCTCG => 9e-05, CCCTCT => -5.7e-05, CCCTGA => -0.00025, CCCTGC => -4.9e-05, CCCTGG => -7.6e-05, CCCTGT => -0.00011, 
CCCTTA => -6.8e-05, CCCTTC => 0.00019, CCCTTG => 1.5e-05, CCCTTT => -5.5e-05, CCGAAA => -4.4e-05, CCGAAC => 2.5e-05, 
CCGAAG => 6.6e-05, CCGAAT => 3.2e-05, CCGACA => 3.2e-06, CCGACC => 2.7e-05, CCGACG => 3.3e-05, CCGACT => -3.9e-06, 
CCGAGA => -7.1e-05, CCGAGC => -5e-06, CCGAGG => -6e-05, CCGAGT => -8.3e-06, CCGATA => -3.4e-05, CCGATC => 9.3e-06, 
CCGATG => 6.6e-05, CCGATT => -1.6e-05, CCGCAA => 2.2e-05, CCGCAC => 7e-05, CCGCAG => 0.00026, CCGCAT => 2.5e-05, 
CCGCCA => 0.00015, CCGCCC => 0.00014, CCGCCG => 0.0002, CCGCCT => 4.8e-05, CCGCGA => 1.1e-06, CCGCGC => 2.1e-05, 
CCGCGG => -2.1e-05, CCGCGT => -5.6e-06, CCGCTA => 5.5e-06, CCGCTC => 6.1e-05, CCGCTG => 0.00031, CCGCTT => -3.5e-05, 
CCGGAA => 0.00012, CCGGAC => 0.00016, CCGGAG => 0.00035, CCGGAT => 0.00015, CCGGCA => 4.1e-05, CCGGCC => 0.00013, 
CCGGCG => 8.1e-05, CCGGCT => 4.6e-05, CCGGGA => 3.8e-05, CCGGGC => 0.0001, CCGGGG => -4.3e-06, CCGGGT => 3.1e-05, 
CCGGTA => 3.2e-06, CCGGTC => 4.6e-05, CCGGTG => 0.00022, CCGGTT => 2.3e-05, CCGTAA => -6.6e-05, CCGTAC => 0.0001, 
CCGTAG => -4.5e-05, CCGTAT => 6e-05, CCGTCA => 7.2e-06, CCGTCC => 6.7e-05, CCGTCG => 4.6e-05, CCGTCT => 9.6e-06, 
CCGTGA => -9e-05, CCGTGC => -2.4e-05, CCGTGG => -2.8e-05, CCGTGT => -4.4e-05, CCGTTA => -3.1e-05, CCGTTC => 6.4e-05, 
CCGTTG => 1.5e-05, CCGTTT => -2.6e-05, CCTAAA => 9.4e-06, CCTAAC => 1.5e-05, CCTAAG => 2.7e-05, CCTAAT => -2.5e-05, 
CCTACA => -2.7e-05, CCTACC => 3.2e-06, CCTACG => 8.2e-06, CCTACT => -1.6e-06, CCTAGA => -7.4e-05, CCTAGC => -2.6e-06, 
CCTAGG => -6.3e-05, CCTAGT => -1e-05, CCTATA => -7.2e-05, CCTATC => 2.1e-06, CCTATG => -1.1e-06, CCTATT => -2.1e-05, 
CCTCAA => -1.2e-05, CCTCAC => -4.3e-05, CCTCAG => 0.00014, CCTCAT => -6.4e-05, CCTCCA => 0.0002, CCTCCC => -9.8e-05, 
CCTCCG => 6e-05, CCTCCT => 0.0001, CCTCGA => 3.6e-05, CCTCGC => 9.9e-06, CCTCGG => 4.4e-06, CCTCGT => -2e-05, 
CCTCTA => -5.1e-05, CCTCTC => -0.00014, CCTCTG => 8.9e-05, CCTCTT => -0.00013, CCTGAA => 0.00017, CCTGAC => 0.00015, 
CCTGAG => 0.00036, CCTGAT => 0.0002, CCTGCA => 2.8e-05, CCTGCC => 5.3e-05, CCTGCG => 1.8e-05, CCTGCT => 0.00011, 
CCTGGA => 0.00017, CCTGGC => 0.00011, CCTGGG => -4.5e-05, CCTGGT => 0.00011, CCTGTA => -7.6e-05, CCTGTC => 1.9e-08, 
CCTGTG => 0.00022, CCTGTT => -7.6e-06, CCTTAA => -0.00023, CCTTAC => 5.6e-05, CCTTAG => -0.00013, CCTTAT => -4.6e-07, 
CCTTCA => -5.8e-05, CCTTCC => -7.8e-05, CCTTCG => 6.8e-06, CCTTCT => -5.1e-05, CCTTGA => -0.00023, CCTTGC => -9.7e-05, 
CCTTGG => -0.00016, CCTTGT => -0.00015, CCTTTA => -0.00016, CCTTTC => -8.2e-05, CCTTTG => -9.9e-05, CCTTTT => -0.00028, 
CGAAAA => -9e-05, CGAAAC => 5.1e-06, CGAAAG => 7.4e-05, CGAAAT => -5.4e-05, CGAACA => -2.7e-05, CGAACC => 1.9e-05, 
CGAACG => -2.5e-05, CGAACT => -1.9e-05, CGAAGA => -2.4e-05, CGAAGC => 3e-06, CGAAGG => 4.2e-07, CGAAGT => -1e-05, 
CGAATA => -5.2e-05, CGAATC => 1.6e-05, CGAATG => 4.9e-05, CGAATT => -3.3e-05, CGACAA => -1.8e-05, CGACAC => 4.6e-05, 
CGACAG => 0.00011, CGACAT => -1.8e-06, CGACCA => 2.9e-06, CGACCC => 5.7e-05, CGACCG => 1.4e-05, CGACCT => 8.5e-08, 
CGACGA => 2.8e-06, CGACGC => 4.7e-05, CGACGG => 9.5e-06, CGACGT => 1.7e-05, CGACTA => 6.5e-07, CGACTC => 2.6e-05, 
CGACTG => 0.00014, CGACTT => -1.8e-05, CGAGAA => 6.6e-05, CGAGAC => 0.00012, CGAGAG => 0.0002, CGAGAT => 0.00012, 
CGAGCA => 1.5e-05, CGAGCC => 8.6e-05, CGAGCG => -1.1e-05, CGAGCT => 5.3e-05, CGAGGA => 3.2e-05, CGAGGC => 4.1e-05, 
CGAGGG => 1.5e-05, CGAGGT => 2.7e-05, CGAGTA => -5.8e-06, CGAGTC => 1.8e-05, CGAGTG => 0.00011, CGAGTT => -1e-05, 
CGATAA => -9.3e-05, CGATAC => 5.2e-05, CGATAG => -5.7e-05, CGATAT => -1.6e-06, CGATCA => -3.3e-05, CGATCC => 2.7e-05, 
CGATCG => 6.5e-06, CGATCT => 8.9e-06, CGATGA => -0.0001, CGATGC => -2.3e-05, CGATGG => -2.2e-05, CGATGT => -3.5e-05, 
CGATTA => -6e-05, CGATTC => 3.2e-05, CGATTG => 1e-05, CGATTT => -5.4e-05, CGCAAA => 0.00014, CGCAAC => 0.00023, 
CGCAAG => 0.00048, CGCAAT => 0.0001, CGCACA => 1.3e-05, CGCACC => 0.00023, CGCACG => 5.3e-05, CGCACT => 4.8e-05, 
CGCAGA => 6.1e-06, CGCAGC => 0.00019, CGCAGG => 7.7e-05, CGCAGT => 9.6e-05, CGCATA => -7.2e-08, CGCATC => 0.00033, 
CGCATG => 0.00024, CGCATT => 0.00012, CGCCAA => 3.1e-05, CGCCAC => 9e-05, CGCCAG => 0.0003, CGCCAT => 1.5e-05, 
CGCCCA => -3e-05, CGCCCC => -3.9e-05, CGCCCG => -3.1e-05, CGCCCT => -3.4e-05, CGCCGA => -5.6e-06, CGCCGC => 6.6e-05, 
CGCCGG => 4.8e-05, CGCCGT => 2.5e-05, CGCCTA => 2e-05, CGCCTC => 8.7e-05, CGCCTG => 0.00036, CGCCTT => 1.8e-05, 
CGCGAA => 7.2e-05, CGCGAC => 0.00014, CGCGAG => 0.00028, CGCGAT => 0.00011, CGCGCA => -2.6e-05, CGCGCC => 6.8e-05, 
CGCGCG => -3.9e-05, CGCGCT => -2.3e-05, CGCGGA => 8.8e-06, CGCGGC => 8.6e-05, CGCGGG => -2.6e-05, CGCGGT => -2.2e-06, 
CGCGTA => -8.2e-06, CGCGTC => 5.8e-05, CGCGTG => 9.6e-05, CGCGTT => -1.6e-06, CGCTAA => -6.7e-05, CGCTAC => 0.00028, 
CGCTAG => -4.2e-05, CGCTAT => 0.00012, CGCTCA => 1.4e-05, CGCTCC => 0.00019, CGCTCG => 3e-05, CGCTCT => 2.7e-05, 
CGCTGA => -0.00013, CGCTGC => 0.00012, CGCTGG => 0.00013, CGCTGT => 4.1e-05, CGCTTA => -3e-05, CGCTTC => 0.00032, 
CGCTTG => 3e-05, CGCTTT => 6.4e-05, CGGAAA => 7.2e-05, CGGAAC => 9.1e-05, CGGAAG => 0.00019, CGGAAT => 4.4e-05, 
CGGACA => 6.8e-06, CGGACC => 5.6e-05, CGGACG => 9.5e-06, CGGACT => 7.2e-06, CGGAGA => -1.2e-06, CGGAGC => 2.9e-05, 
CGGAGG => 8.4e-06, CGGAGT => 3.6e-07, CGGATA => 6.6e-07, CGGATC => 7.4e-05, CGGATG => 8.1e-05, CGGATT => 2.1e-05, 
CGGCAA => 4.6e-06, CGGCAC => 8.1e-05, CGGCAG => 0.00026, CGGCAT => 1.2e-05, CGGCCA => 7.3e-06, CGGCCC => 0.00011, 
CGGCCG => 7.3e-06, CGGCCT => 2.2e-05, CGGCGA => 1.4e-05, CGGCGC => 6.7e-05, CGGCGG => -4.4e-07, CGGCGT => 1.5e-05, 
CGGCTA => 7.9e-06, CGGCTC => 6.4e-05, CGGCTG => 0.00028, CGGCTT => -1.9e-05, CGGGAA => 0.00015, CGGGAC => 0.00023, 
CGGGAG => 0.0004, CGGGAT => 0.00016, CGGGCA => 4.3e-05, CGGGCC => 0.00017, CGGGCG => -4.8e-06, CGGGCT => 5.2e-05, 
CGGGGA => -1.4e-05, CGGGGC => 7.3e-05, CGGGGG => -1.8e-05, CGGGGT => 1e-05, CGGGTA => 6.3e-06, CGGGTC => 6.4e-05, 
CGGGTG => 0.00013, CGGGTT => -9.5e-07, CGGTAA => -7.1e-05, CGGTAC => 6.5e-05, CGGTAG => -5.7e-05, CGGTAT => 4.6e-06, 
CGGTCA => -1.3e-05, CGGTCC => 2.1e-05, CGGTCG => -1.5e-05, CGGTCT => -9.4e-06, CGGTGA => -8.6e-05, CGGTGC => -4.1e-05, 
CGGTGG => -6.9e-05, CGGTGT => -4.2e-05, CGGTTA => -3.7e-05, CGGTTC => 1.5e-05, CGGTTG => -6e-06, CGGTTT => -1.9e-05, 
CGTAAA => -3.6e-05, CGTAAC => 7.1e-06, CGTAAG => 3.4e-05, CGTAAT => -2.2e-05, CGTACA => -2.7e-05, CGTACC => 1.9e-05, 
CGTACG => -1.5e-06, CGTACT => -8.8e-06, CGTAGA => -1.5e-05, CGTAGC => -9.6e-07, CGTAGG => -1.1e-05, CGTAGT => -1.8e-05, 
CGTATA => -4.9e-05, CGTATC => 3.8e-05, CGTATG => -5.4e-06, CGTATT => -2.7e-05, CGTCAA => 8.5e-06, CGTCAC => 1.5e-05, 
CGTCAG => 0.00014, CGTCAT => -2e-05, CGTCCA => -9.7e-06, CGTCCC => 2.8e-05, CGTCCG => 5.9e-06, CGTCCT => 1.5e-05, 
CGTCGA => 3.2e-06, CGTCGC => 7.5e-05, CGTCGG => -5.7e-06, CGTCGT => 3.1e-05, CGTCTA => -7.3e-06, CGTCTC => 2.6e-05, 
CGTCTG => 0.00017, CGTCTT => -1.8e-05, CGTGAA => 4e-05, CGTGAC => 7.8e-05, CGTGAG => 0.00014, CGTGAT => 4.6e-05, 
CGTGCA => -1e-05, CGTGCC => 0.00011, CGTGCG => -2.6e-05, CGTGCT => 2.8e-05, CGTGGA => 3.8e-05, CGTGGC => 7.6e-05, 
CGTGGG => -1.1e-05, CGTGGT => 2.7e-05, CGTGTA => -2.4e-05, CGTGTC => 1.7e-05, CGTGTG => 9.7e-05, CGTGTT => -5.9e-05, 
CGTTAA => -0.0001, CGTTAC => 3.4e-05, CGTTAG => -5.2e-05, CGTTAT => 1e-05, CGTTCA => -2.7e-05, CGTTCC => 3.9e-05, 
CGTTCG => 6.2e-06, CGTTCT => -1.8e-05, CGTTGA => -8.9e-05, CGTTGC => -2e-05, CGTTGG => -2.9e-05, CGTTGT => -4.8e-05, 
CGTTTA => -8.2e-05, CGTTTC => 1.3e-05, CGTTTG => -1.8e-05, CGTTTT => -0.00013, CTAAAA => -0.00024, CTAAAC => -4.1e-05, 
CTAAAG => 0.00011, CTAAAT => -0.00012, CTAACA => -7.3e-05, CTAACC => -4.6e-07, CTAACG => -2.3e-07, CTAACT => -8.5e-05, 
CTAAGA => -9e-05, CTAAGC => -3.3e-05, CTAAGG => -6.3e-05, CTAAGT => -6.8e-05, CTAATA => -0.00016, CTAATC => -3.2e-05, 
CTAATG => -5.6e-07, CTAATT => -0.00019, CTACAA => -0.0001, CTACAC => -7.8e-05, CTACAG => 0.00012, CTACAT => -0.00014, 
CTACCA => -5.1e-05, CTACCC => 5.5e-06, CTACCG => -3.7e-07, CTACCT => -7.2e-05, CTACGA => -7e-07, CTACGC => 2.4e-05, 
CTACGG => 2e-05, CTACGT => -2.1e-05, CTACTA => -8.3e-05, CTACTC => -2.5e-05, CTACTG => 2.1e-05, CTACTT => -0.00013, 
CTAGAA => 3.9e-05, CTAGAC => 6.2e-05, CTAGAG => 0.00017, CTAGAT => 5e-05, CTAGCA => -4.1e-05, CTAGCC => 3.5e-05, 
CTAGCG => 9.8e-06, CTAGCT => -2.1e-05, CTAGGA => -4.6e-05, CTAGGC => 1e-05, CTAGGG => -3.5e-05, CTAGGT => -4.1e-05, 
CTAGTA => -5.7e-05, CTAGTC => -1.8e-05, CTAGTG => 4.2e-05, CTAGTT => -0.0001, CTATAA => -0.00023, CTATAC => -5.1e-05, 
CTATAG => -0.00015, CTATAT => -0.00017, CTATCA => -8.4e-05, CTATCC => -4.7e-06, CTATCG => 7.3e-06, CTATCT => -0.0001, 
CTATGA => -0.00017, CTATGC => -8.4e-05, CTATGG => -7.6e-05, CTATGT => -0.00014, CTATTA => -0.00017, CTATTC => -6.9e-05, 
CTATTG => -8.6e-05, CTATTT => -0.00029, CTCAAA => 7.1e-05, CTCAAC => 0.00044, CTCAAG => 0.00057, CTCAAT => 0.00018, 
CTCACA => 1.1e-05, CTCACC => 0.00033, CTCACG => 9.1e-05, CTCACT => 2.1e-05, CTCAGA => -0.00015, CTCAGC => 0.00021, 
CTCAGG => -6.1e-05, CTCAGT => 4.8e-05, CTCATA => -5.2e-05, CTCATC => 0.00048, CTCATG => 0.00024, CTCATT => 4.5e-05, 
CTCCAA => -7.3e-05, CTCCAC => -6.1e-05, CTCCAG => 0.00015, CTCCAT => -0.00012, CTCCCA => -0.00015, CTCCCC => -0.00013, 
CTCCCG => -2.9e-05, CTCCCT => -0.00016, CTCCGA => 2.9e-05, CTCCGC => 4.7e-05, CTCCGG => 5.1e-05, CTCCGT => -2.1e-05, 
CTCCTA => -4.8e-05, CTCCTC => -4.4e-05, CTCCTG => 0.00029, CTCCTT => -0.00012, CTCGAA => -1.2e-05, CTCGAC => 9.6e-05, 
CTCGAG => 0.00014, CTCGAT => 6.8e-05, CTCGCA => -3.1e-05, CTCGCC => 7.1e-05, CTCGCG => -1.7e-05, CTCGCT => -3.6e-05, 
CTCGGA => 9e-06, CTCGGC => 6.6e-05, CTCGGG => -2e-06, CTCGGT => 6.4e-07, CTCGTA => -2.2e-05, CTCGTC => 1.7e-05, 
CTCGTG => 8e-05, CTCGTT => -2.9e-05, CTCTAA => -0.00022, CTCTAC => 0.00031, CTCTAG => -0.00016, CTCTAT => 9.3e-05, 
CTCTCA => -8.1e-05, CTCTCC => 6.7e-05, CTCTCG => 4e-05, CTCTCT => -0.00024, CTCTGA => -0.00034, CTCTGC => -5.7e-05, 
CTCTGG => -0.00012, CTCTGT => -0.00021, CTCTTA => -0.00013, CTCTTC => 0.00036, CTCTTG => -6.8e-05, CTCTTT => -9.5e-05, 
CTGAAA => 0.00017, CTGAAC => 0.00039, CTGAAG => 0.001, CTGAAT => 0.00012, CTGACA => 8.3e-05, CTGACC => 0.00042, 
CTGACG => 0.00019, CTGACT => 5.2e-05, CTGAGA => -6.8e-05, CTGAGC => 0.00033, CTGAGG => 7e-05, CTGAGT => 8.2e-05, 
CTGATA => -5.1e-05, CTGATC => 0.00038, CTGATG => 0.00038, CTGATT => 6.9e-05, CTGCAA => 8.8e-05, CTGCAC => 0.00035, 
CTGCAG => 0.0013, CTGCAT => 1.9e-05, CTGCCA => 8.1e-05, CTGCCC => 0.00053, CTGCCG => 0.00021, CTGCCT => 2e-05, 
CTGCGA => 0.00014, CTGCGC => 0.00051, CTGCGG => 0.00033, CTGCGT => 0.00014, CTGCTA => 6e-05, CTGCTC => 0.00052, 
CTGCTG => 0.0015, CTGCTT => -1.3e-05, CTGGAA => 0.00049, CTGGAC => 0.0011, CTGGAG => 0.0019, CTGGAT => 0.00063, 
CTGGCA => 0.00026, CTGGCC => 0.0011, CTGGCG => 0.0003, CTGGCT => 0.0004, CTGGGA => 0.00014, CTGGGC => 0.0007, 
CTGGGG => 9.9e-05, CTGGGT => 0.00014, CTGGTA => 4.7e-05, CTGGTC => 0.00038, CTGGTG => 0.001, CTGGTT => 9.4e-05, 
CTGTAA => -0.00036, CTGTAC => 0.00028, CTGTAG => -0.00022, CTGTAT => 1.5e-05, CTGTCA => -2.9e-05, CTGTCC => 0.00031, 
CTGTCG => 0.00012, CTGTCT => -6.5e-06, CTGTGA => -0.00038, CTGTGC => 8.6e-05, CTGTGG => -2.2e-05, CTGTGT => -0.00017, 
CTGTTA => -0.00012, CTGTTC => 0.00017, CTGTTG => 4.3e-05, CTGTTT => -0.00012, CTTAAA => -0.00024, CTTAAC => -2.9e-05, 
CTTAAG => 1.8e-05, CTTAAT => -0.00011, CTTACA => -9.6e-05, CTTACC => 4.6e-06, CTTACG => -1.2e-05, CTTACT => -6.4e-05, 
CTTAGA => -0.00013, CTTAGC => -3e-05, CTTAGG => -7.9e-05, CTTAGT => -6.3e-05, CTTATA => -0.00012, CTTATC => 2.1e-06, 
CTTATG => -2.6e-05, CTTATT => -0.00015, CTTCAA => -8.9e-05, CTTCAC => -7.1e-05, CTTCAG => 0.00017, CTTCAT => -0.00014, 
CTTCCA => -5.9e-05, CTTCCC => -0.00012, CTTCCG => -1.8e-05, CTTCCT => -0.00013, CTTCGA => 4.3e-05, CTTCGC => 4.5e-05, 
CTTCGG => 4.7e-05, CTTCGT => -1.5e-06, CTTCTA => -9.6e-05, CTTCTC => -0.0001, CTTCTG => 6.7e-05, CTTCTT => -0.00017, 
CTTGAA => -5.1e-05, CTTGAC => 4.5e-05, CTTGAG => 8e-05, CTTGAT => 2.5e-05, CTTGCA => -4.3e-05, CTTGCC => 1.2e-05, 
CTTGCG => -1.5e-05, CTTGCT => -6.2e-05, CTTGGA => -8.8e-06, CTTGGC => 8e-06, CTTGGG => -0.0001, CTTGGT => -4.4e-05, 
CTTGTA => -0.00014, CTTGTC => -4.5e-05, CTTGTG => -6.8e-06, CTTGTT => -0.00018, CTTTAA => -0.00049, CTTTAC => -6e-05, 
CTTTAG => -0.00023, CTTTAT => -0.00022, CTTTCA => -0.0002, CTTTCC => -0.00013, CTTTCG => -3.2e-05, CTTTCT => -0.00026, 
CTTTGA => -0.00037, CTTTGC => -0.00021, CTTTGG => -0.00023, CTTTGT => -0.00031, CTTTTA => -0.00038, CTTTTC => -0.00023, 
CTTTTG => -0.00021, CTTTTT => -0.00061, GAAAAA => -9.7e-05, GAAAAC => 0.00017, GAAAAG => 0.00046, GAAAAT => 1.3e-05, 
GAAACA => -7.6e-05, GAAACC => 9e-05, GAAACG => 5e-05, GAAACT => -6.1e-06, GAAAGA => -0.00011, GAAAGC => 3.5e-05, 
GAAAGG => -3.5e-05, GAAAGT => 1.9e-05, GAAATA => -0.00019, GAAATC => 0.00017, GAAATG => 0.00015, GAAATT => -2.7e-05, 
GAACAA => 5.1e-05, GAACAC => 2.2e-05, GAACAG => 0.00038, GAACAT => -4.4e-05, GAACCA => 8.3e-05, GAACCC => 0.00011, 
GAACCG => 5.9e-05, GAACCT => 8.5e-05, GAACGA => 5.5e-05, GAACGC => 0.00012, GAACGG => 9.2e-05, GAACGT => 3.9e-05, 
GAACTA => 2.5e-05, GAACTC => 9.7e-05, GAACTG => 0.00046, GAACTT => 7.1e-05, GAAGAA => 0.00078, GAAGAC => 0.00054, 
GAAGAG => 0.001, GAAGAT => 0.00065, GAAGCA => 0.00015, GAAGCC => 0.00038, GAAGCG => 7.7e-05, GAAGCT => 0.0003, 
GAAGGA => 0.00017, GAAGGC => 0.00023, GAAGGG => 7e-05, GAAGGT => 0.00012, GAAGTA => 6.3e-05, GAAGTC => 0.00012, 
GAAGTG => 0.00039, GAAGTT => 9.9e-05, GAATAA => -0.00035, GAATAC => 0.00015, GAATAG => -0.00015, GAATAT => -8.3e-06, 
GAATCA => -2.4e-05, GAATCC => 8.2e-05, GAATCG => 4.8e-05, GAATCT => 5.2e-05, GAATGA => -0.00032, GAATGC => 1.9e-05, 
GAATGG => 1.5e-05, GAATGT => 1.4e-05, GAATTA => -0.00011, GAATTC => 8.9e-05, GAATTG => 1.9e-05, GAATTT => -0.00011, 
GACAAA => 0.00031, GACAAC => 0.00053, GACAAG => 0.00085, GACAAT => 0.00028, GACACA => 0.00016, GACACC => 0.00043, 
GACACG => 0.0002, GACACT => 0.00019, GACAGA => 7.4e-06, GACAGC => 0.00054, GACAGG => 9.3e-05, GACAGT => 0.00026, 
GACATA => 8.7e-05, GACATC => 0.00069, GACATG => 0.00058, GACATT => 0.0003, GACCAA => 6.3e-05, GACCAC => 0.00019, 
GACCAG => 0.00053, GACCAT => 5.8e-05, GACCCA => 0.00014, GACCCC => 0.00029, GACCCG => 0.00011, GACCCT => 0.00019, 
GACCGA => 9.8e-05, GACCGC => 0.0002, GACCGG => 0.00016, GACCGT => 8.9e-05, GACCTA => 6.1e-05, GACCTC => 0.00032, 
GACCTG => 0.00084, GACCTT => 0.0001, GACGAA => 0.0002, GACGAC => 0.00047, GACGAG => 0.00076, GACGAT => 0.00033, 
GACGCA => 5.3e-05, GACGCC => 0.00026, GACGCG => 7.6e-05, GACGCT => 9.3e-05, GACGGA => 9.9e-05, GACGGC => 0.00026, 
GACGGG => 0.0001, GACGGT => 8.3e-05, GACGTA => 2.3e-05, GACGTC => 0.00016, GACGTG => 0.00035, GACGTT => 4.6e-05, 
GACTAA => -0.00014, GACTAC => 0.00048, GACTAG => -0.0001, GACTAT => 0.00026, GACTCA => 6.3e-05, GACTCC => 0.00032, 
GACTCG => 0.00018, GACTCT => 0.00018, GACTGA => -0.00026, GACTGC => 0.00021, GACTGG => 0.00021, GACTGT => 7e-05, 
GACTTA => -1.6e-05, GACTTC => 0.00054, GACTTG => 0.00021, GACTTT => 0.00024, GAGAAA => 0.00061, GAGAAC => 0.00079, 
GAGAAG => 0.0017, GAGAAT => 0.0004, GAGACA => 0.00012, GAGACC => 0.0005, GAGACG => 0.00027, GAGACT => 0.00018, 
GAGAGA => -4.5e-05, GAGAGC => 0.00047, GAGAGG => 0.00014, GAGAGT => 0.00023, GAGATA => 4.9e-05, GAGATC => 0.0008, 
GAGATG => 0.00062, GAGATT => 0.00032, GAGCAA => 0.00023, GAGCAC => 0.00038, GAGCAG => 0.0012, GAGCAT => 0.00013, 
GAGCCA => 0.00019, GAGCCC => 0.00041, GAGCCG => 0.00021, GAGCCT => 0.00024, GAGCGA => 0.00018, GAGCGC => 0.00043, 
GAGCGG => 0.00031, GAGCGT => 0.00015, GAGCTA => 0.00013, GAGCTC => 0.00045, GAGCTG => 0.0016, GAGCTT => 0.00021, 
GAGGAA => 0.00095, GAGGAC => 0.0011, GAGGAG => 0.0026, GAGGAT => 0.00094, GAGGCA => 0.00026, GAGGCC => 0.00081, 
GAGGCG => 0.00034, GAGGCT => 0.00042, GAGGGA => 0.00013, GAGGGC => 0.00058, GAGGGG => 7.2e-05, GAGGGT => 0.00023, 
GAGGTA => 6.9e-05, GAGGTC => 0.00032, GAGGTG => 0.00093, GAGGTT => 0.00018, GAGTAA => -0.00018, GAGTAC => 0.00048, 
GAGTAG => -0.00014, GAGTAT => 0.00022, GAGTCA => 8.2e-05, GAGTCC => 0.0003, GAGTCG => 0.00017, GAGTCT => 0.00018, 
GAGTGA => -0.00024, GAGTGC => 0.00023, GAGTGG => 0.00016, GAGTGT => 0.00011, GAGTTA => -1.1e-05, GAGTTC => 0.00053, 
GAGTTG => 0.00022, GAGTTT => 0.0003, GATAAA => 0.00013, GATAAC => 0.00015, GATAAG => 0.00028, GATAAT => 8.3e-05, 
GATACA => -1.5e-05, GATACC => 0.0001, GATACG => 6.2e-05, GATACT => 2.2e-05, GATAGA => -3.8e-05, GATAGC => 8.3e-05, 
GATAGG => -1.2e-05, GATAGT => 2.7e-05, GATATA => -3.6e-05, GATATC => 0.00018, GATATG => 0.00015, GATATT => 6.4e-05, 
GATCAA => 6.3e-05, GATCAC => 9.9e-05, GATCAG => 0.00038, GATCAT => 4.4e-05, GATCCA => 0.00014, GATCCC => 0.00026, 
GATCCG => 0.00011, GATCCT => 0.00016, GATCGA => 7.7e-05, GATCGC => 0.0002, GATCGG => 8.6e-05, GATCGT => 8.4e-05, 
GATCTA => 4.4e-05, GATCTC => 0.00016, GATCTG => 0.00056, GATCTT => 0.00011, GATGAA => 0.00081, GATGAC => 0.00078, 
GATGAG => 0.0012, GATGAT => 0.00074, GATGCA => 0.00026, GATGCC => 0.00064, GATGCG => 0.00016, GATGCT => 0.00037, 
GATGGA => 0.00039, GATGGC => 0.00054, GATGGG => 0.00022, GATGGT => 0.00025, GATGTA => 5.4e-05, GATGTC => 0.00032, 
GATGTG => 0.00074, GATGTT => 0.00018, GATTAA => -0.00026, GATTAC => 0.00014, GATTAG => -0.00014, GATTAT => 2e-05, 
GATTCA => 1.9e-05, GATTCC => 8.3e-05, GATTCG => 8e-05, GATTCT => 3.3e-05, GATTGA => -0.00021, GATTGC => -2.3e-05, 
GATTGG => -1.9e-05, GATTGT => -9.9e-05, GATTTA => -0.00014, GATTTC => 9.5e-05, GATTTG => 7.4e-05, GATTTT => -0.00021, 
GCAAAA => -0.00013, GCAAAC => -6e-05, GCAAAG => 6.8e-05, GCAAAT => -0.00011, GCAACA => -3.7e-05, GCAACC => 9.4e-06, 
GCAACG => -4.2e-06, GCAACT => -1.4e-05, GCAAGA => -0.0001, GCAAGC => -5.7e-05, GCAAGG => -7.1e-05, GCAAGT => -5.7e-05, 
GCAATA => -0.00013, GCAATC => -2.4e-05, GCAATG => 3.5e-05, GCAATT => -0.00011, GCACAA => -5.7e-05, GCACAC => -0.0001, 
GCACAG => 0.00011, GCACAT => -9.9e-05, GCACCA => 5.3e-05, GCACCC => 8.9e-05, GCACCG => 4.5e-05, GCACCT => 4.9e-05, 
GCACGA => 1.4e-05, GCACGC => 9e-06, GCACGG => 2e-05, GCACGT => -1.9e-05, GCACTA => -2.9e-05, GCACTC => -1e-05, 
GCACTG => 0.0002, GCACTT => -8.4e-05, GCAGAA => 0.0002, GCAGAC => 0.00021, GCAGAG => 0.00048, GCAGAT => 0.00022, 
GCAGCA => 0.00014, GCAGCC => 0.00035, GCAGCG => 5.8e-05, GCAGCT => 0.00023, GCAGGA => 8e-05, GCAGGC => 6.8e-05, 
GCAGGG => -2.4e-05, GCAGGT => 6.4e-06, GCAGTA => -2.1e-05, GCAGTC => 1.2e-05, GCAGTG => 0.00014, GCAGTT => -3.3e-05, 
GCATAA => -0.00021, GCATAC => -2.1e-05, GCATAG => -0.00013, GCATAT => -9.7e-05, GCATCA => -2.2e-05, GCATCC => 7.7e-05, 
GCATCG => 2.6e-05, GCATCT => 1.7e-05, GCATGA => -0.0002, GCATGC => -0.00011, GCATGG => -9.6e-05, GCATGT => -0.00014, 
GCATTA => -9.2e-05, GCATTC => -3.6e-05, GCATTG => -4e-05, GCATTT => -0.00022, GCCAAA => 0.00039, GCCAAC => 0.00058, 
GCCAAG => 0.0011, GCCAAT => 0.00038, GCCACA => 0.00025, GCCACC => 0.0006, GCCACG => 0.00022, GCCACT => 0.0002, 
GCCAGA => 2.6e-05, GCCAGC => 0.00048, GCCAGG => 1.8e-05, GCCAGT => 0.00028, GCCATA => 6.4e-05, GCCATC => 0.00074, 
GCCATG => 0.00062, GCCATT => 0.00032, GCCCAA => 7.7e-05, GCCCAC => 0.00016, GCCCAG => 0.00063, GCCCAT => 5.9e-05, 
GCCCCA => 4.5e-05, GCCCCC => 9.7e-05, GCCCCG => 5e-05, GCCCCT => 6.4e-05, GCCCGA => 7.8e-05, GCCCGC => 0.00017, 
GCCCGG => 0.00014, GCCCGT => 7e-05, GCCCTA => 3.4e-05, GCCCTC => 0.00021, GCCCTG => 0.00074, GCCCTT => 4.9e-05, 
GCCGAA => 0.00017, GCCGAC => 0.00024, GCCGAG => 0.00054, GCCGAT => 0.00022, GCCGCA => 8.1e-05, GCCGCC => 0.00045, 
GCCGCG => 5.7e-05, GCCGCT => 0.00014, GCCGGA => 0.00014, GCCGGC => 0.00025, GCCGGG => 1.7e-05, GCCGGT => 9.5e-05, 
GCCGTA => 3.5e-05, GCCGTC => 0.00019, GCCGTG => 0.00034, GCCGTT => 5.3e-05, GCCTAA => -0.00012, GCCTAC => 0.00039, 
GCCTAG => -0.00012, GCCTAT => 0.0002, GCCTCA => 1.8e-05, GCCTCC => 0.00028, GCCTCG => 0.00013, GCCTCT => 6.4e-05, 
GCCTGA => -0.00023, GCCTGC => 0.00013, GCCTGG => 2.6e-05, GCCTGT => -7e-06, GCCTTA => -2.7e-05, GCCTTC => 0.00056, 
GCCTTG => 0.00017, GCCTTT => 0.00015, GCGAAA => -4.4e-05, GCGAAC => 2.3e-05, GCGAAG => 7.2e-05, GCGAAT => 1e-06, 
GCGACA => -9.7e-06, GCGACC => 2.7e-05, GCGACG => 1.9e-05, GCGACT => -1.8e-05, GCGAGA => -6.3e-05, GCGAGC => -6.7e-06, 
GCGAGG => -5.2e-05, GCGAGT => -2e-05, GCGATA => -2.6e-05, GCGATC => 1.6e-05, GCGATG => 6.6e-05, GCGATT => -5.9e-06, 
GCGCAA => -8e-06, GCGCAC => 3.5e-05, GCGCAG => 0.00016, GCGCAT => -1.9e-05, GCGCCA => 1e-05, GCGCCC => 0.00011, 
GCGCCG => 3.3e-05, GCGCCT => 2.4e-06, GCGCGA => -1.1e-05, GCGCGC => -1.1e-05, GCGCGG => -3.4e-05, GCGCGT => -2.3e-05, 
GCGCTA => 4.5e-06, GCGCTC => 6.9e-05, GCGCTG => 0.00034, GCGCTT => -2.5e-05, GCGGAA => 9.1e-05, GCGGAC => 0.00018, 
GCGGAG => 0.00036, GCGGAT => 0.00015, GCGGCA => 9.3e-05, GCGGCC => 0.00033, GCGGCG => 0.00013, GCGGCT => 0.00012, 
GCGGGA => 1.8e-05, GCGGGC => 0.00014, GCGGGG => -1.9e-05, GCGGGT => 4.1e-05, GCGGTA => -5.1e-06, GCGGTC => 6.7e-05, 
GCGGTG => 0.00025, GCGGTT => 1.7e-05, GCGTAA => -7.5e-05, GCGTAC => 6.7e-05, GCGTAG => -5.1e-05, GCGTAT => 1.2e-05, 
GCGTCA => -2.1e-05, GCGTCC => 6.5e-05, GCGTCG => -1.3e-05, GCGTCT => 1.2e-05, GCGTGA => -9.7e-05, GCGTGC => -4e-05, 
GCGTGG => -5.8e-05, GCGTGT => -8.3e-05, GCGTTA => -2.5e-05, GCGTTC => 3e-05, GCGTTG => 1.4e-05, GCGTTT => -2.9e-05, 
GCTAAA => 3.3e-05, GCTAAC => 4e-05, GCTAAG => 0.00011, GCTAAT => -2.5e-05, GCTACA => -1.1e-05, GCTACC => 4.2e-05, 
GCTACG => 1.9e-05, GCTACT => -8.9e-06, GCTAGA => -6.3e-05, GCTAGC => 1.6e-05, GCTAGG => -5.6e-05, GCTAGT => -4.9e-06, 
GCTATA => -4.5e-05, GCTATC => 5.5e-05, GCTATG => 7.6e-05, GCTATT => 6.4e-06, GCTCAA => 3.7e-05, GCTCAC => 1.1e-05, 
GCTCAG => 0.00028, GCTCAT => -6.4e-06, GCTCCA => 0.00013, GCTCCC => 7.3e-05, GCTCCG => 4e-05, GCTCCT => 9.7e-05, 
GCTCGA => 6.5e-05, GCTCGC => 4.9e-05, GCTCGG => 3.6e-05, GCTCGT => 4e-05, GCTCTA => 4e-06, GCTCTC => 1.9e-05, 
GCTCTG => 0.0004, GCTCTT => -5.6e-06, GCTGAA => 0.00022, GCTGAC => 0.00027, GCTGAG => 0.00047, GCTGAT => 0.00024, 
GCTGCA => 0.00017, GCTGCC => 0.00048, GCTGCG => 6.5e-05, GCTGCT => 0.00028, GCTGGA => 0.00019, GCTGGC => 0.00021, 
GCTGGG => -4.1e-05, GCTGGT => 0.0001, GCTGTA => -2e-05, GCTGTC => 0.00011, GCTGTG => 0.00048, GCTGTT => -1.3e-05, 
GCTTAA => -0.00021, GCTTAC => 3.8e-05, GCTTAG => -0.00014, GCTTAT => -2.2e-05, GCTTCA => -5.1e-05, GCTTCC => 1.7e-05, 
GCTTCG => 1.5e-05, GCTTCT => -5.4e-05, GCTTGA => -0.0002, GCTTGC => -7.6e-05, GCTTGG => -0.00012, GCTTGT => -0.00012, 
GCTTTA => -0.00012, GCTTTC => 6.2e-06, GCTTTG => 2.3e-05, GCTTTT => -0.00017, GGAAAA => -0.00011, GGAAAC => 9.7e-05, 
GGAAAG => 0.00018, GGAAAT => -1.9e-05, GGAACA => -4.7e-06, GGAACC => 0.0001, GGAACG => 3.2e-05, GGAACT => 1e-05, 
GGAAGA => -0.00019, GGAAGC => 2.4e-06, GGAAGG => -0.00018, GGAAGT => -1.6e-05, GGAATA => -6.8e-05, GGAATC => 8.9e-05, 
GGAATG => 9.2e-05, GGAATT => -8.4e-06, GGACAA => -2.3e-05, GGACAC => 3.2e-05, GGACAG => 0.00019, GGACAT => -3.2e-05, 
GGACCA => 0.00013, GGACCC => 0.00019, GGACCG => 3.7e-05, GGACCT => 0.00011, GGACGA => 3.6e-05, GGACGC => 7.3e-05, 
GGACGG => 2e-05, GGACGT => 2e-05, GGACTA => -3.6e-05, GGACTC => 4.7e-05, GGACTG => 0.00017, GGACTT => -1.1e-05, 
GGAGAA => 0.00027, GGAGAC => 0.00033, GGAGAG => 0.00045, GGAGAT => 0.00029, GGAGCA => 0.00011, GGAGCC => 0.00026, 
GGAGCG => 3.6e-05, GGAGCT => 0.00015, GGAGGA => 0.00019, GGAGGC => 0.00015, GGAGGG => -0.0001, GGAGGT => 0.0001, 
GGAGTA => 4.7e-06, GGAGTC => 7e-05, GGAGTG => 0.00017, GGAGTT => 2.3e-05, GGATAA => -0.00018, GGATAC => 0.00011, 
GGATAG => -0.00011, GGATAT => -8.1e-06, GGATCA => -3.7e-06, GGATCC => 0.00011, GGATCG => 7.5e-05, GGATCT => 5.5e-05, 
GGATGA => -0.00024, GGATGC => -6.3e-05, GGATGG => -9.6e-05, GGATGT => -0.00012, GGATTA => -0.0001, GGATTC => 0.00012, 
GGATTG => 1.9e-05, GGATTT => -2.2e-05, GGCAAA => 0.00024, GGCAAC => 0.00039, GGCAAG => 0.00066, GGCAAT => 0.00021, 
GGCACA => 0.00011, GGCACC => 0.00041, GGCACG => 9.9e-05, GGCACT => 0.00015, GGCAGA => -0.00011, GGCAGC => 0.00048, 
GGCAGG => -8.3e-05, GGCAGT => 0.00021, GGCATA => 9e-06, GGCATC => 0.00053, GGCATG => 0.00037, GGCATT => 0.00019, 
GGCCAA => -1.4e-05, GGCCAC => 0.00015, GGCCAG => 0.00035, GGCCAT => 9.1e-06, GGCCCA => 3.2e-05, GGCCCC => 0.00013, 
GGCCCG => 1.9e-05, GGCCCT => 7.6e-05, GGCCGA => 5.5e-05, GGCCGC => 0.00017, GGCCGG => 8.3e-05, GGCCGT => 6.7e-05, 
GGCCTA => 1e-05, GGCCTC => 0.00012, GGCCTG => 0.00045, GGCCTT => 2.6e-06, GGCGAA => 0.00012, GGCGAC => 0.00023, 
GGCGAG => 0.00045, GGCGAT => 0.00023, GGCGCA => 3.2e-05, GGCGCC => 0.00021, GGCGCG => 1e-05, GGCGCT => 6.6e-05, 
GGCGGA => 0.00016, GGCGGC => 0.00036, GGCGGG => 6.9e-07, GGCGGT => 0.00017, GGCGTA => 2.4e-05, GGCGTC => 0.00016, 
GGCGTG => 0.00021, GGCGTT => 5.7e-05, GGCTAA => -0.00015, GGCTAC => 0.00044, GGCTAG => -0.00012, GGCTAT => 0.0002, 
GGCTCA => 1.3e-05, GGCTCC => 0.00034, GGCTCG => 7.1e-05, GGCTCT => 0.00011, GGCTGA => -0.00029, GGCTGC => 9e-05, 
GGCTGG => 3.7e-05, GGCTGT => -5.1e-07, GGCTTA => -3.5e-05, GGCTTC => 0.00045, GGCTTG => 6.7e-05, GGCTTT => 9.7e-05, 
GGGAAA => 1.3e-05, GGGAAC => 9.5e-05, GGGAAG => 0.00021, GGGAAT => -2e-06, GGGACA => -3.2e-05, GGGACC => 3.1e-05, 
GGGACG => 1e-05, GGGACT => -3.9e-05, GGGAGA => -0.00016, GGGAGC => -2.1e-05, GGGAGG => -0.00027, GGGAGT => -6.7e-05, 
GGGATA => -5.2e-05, GGGATC => 6.7e-05, GGGATG => 3.5e-05, GGGATT => -4.1e-05, GGGCAA => -6.8e-05, GGGCAC => -2e-05, 
GGGCAG => 6.7e-05, GGGCAT => -6.9e-05, GGGCCA => -6.4e-06, GGGCCC => 9.6e-05, GGGCCG => -1.8e-06, GGGCCT => -3.8e-06, 
GGGCGA => -1.4e-05, GGGCGC => -8.6e-06, GGGCGG => -6.8e-05, GGGCGT => -2.5e-05, GGGCTA => -4.3e-05, GGGCTC => -2.3e-07, 
GGGCTG => 0.00014, GGGCTT => -7.4e-05, GGGGAA => 2.3e-05, GGGGAC => 0.00021, GGGGAG => 0.00025, GGGGAT => 0.00014, 
GGGGCA => 5e-06, GGGGCC => 0.00016, GGGGCG => -3.4e-06, GGGGCT => 2.7e-06, GGGGGA => -8e-05, GGGGGC => 0.00014, 
GGGGGG => -0.00013, GGGGGT => 5.8e-06, GGGGTA => -3.8e-05, GGGGTC => 4.8e-05, GGGGTG => 3.2e-05, GGGGTT => -5.7e-05, 
GGGTAA => -0.00013, GGGTAC => 3e-05, GGGTAG => -0.00011, GGGTAT => -2.8e-05, GGGTCA => -6e-05, GGGTCC => 1e-06, 
GGGTCG => -1.8e-05, GGGTCT => -5e-05, GGGTGA => -0.00018, GGGTGC => -0.0001, GGGTGG => -0.00022, GGGTGT => -0.00014, 
GGGTTA => -7.2e-05, GGGTTC => -5.6e-05, GGGTTG => -6.9e-05, GGGTTT => -0.00014, GGTAAA => -2.7e-05, GGTAAC => 2.4e-05, 
GGTAAG => -1.9e-05, GGTAAT => -3.1e-05, GGTACA => -2.3e-05, GGTACC => 3.1e-05, GGTACG => 2e-06, GGTACT => -7.7e-06, 
GGTAGA => -8.2e-05, GGTAGC => 2.2e-05, GGTAGG => -0.0001, GGTAGT => -1.4e-05, GGTATA => -7.9e-05, GGTATC => 3.6e-05, 
GGTATG => -8.7e-07, GGTATT => -2.5e-05, GGTCAA => 8.9e-06, GGTCAC => 3.5e-05, GGTCAG => 0.00013, GGTCAT => 2.1e-05, 
GGTCCA => 7.4e-05, GGTCCC => 7.5e-05, GGTCCG => 4.1e-05, GGTCCT => 0.0001, GGTCGA => 4.4e-05, GGTCGC => 0.00011, 
GGTCGG => 4.7e-06, GGTCGT => 6.8e-05, GGTCTA => -2.4e-05, GGTCTC => 1.9e-05, GGTCTG => 0.00022, GGTCTT => -1e-05, 
GGTGAA => 0.00017, GGTGAC => 0.00019, GGTGAG => 0.0002, GGTGAT => 0.00019, GGTGCA => 6e-05, GGTGCC => 0.00023, 
GGTGCG => 2.7e-05, GGTGCT => 0.00013, GGTGGA => 0.00021, GGTGGC => 0.00031, GGTGGG => -6.3e-05, GGTGGT => 0.00018, 
GGTGTA => -1.2e-05, GGTGTC => 9.2e-05, GGTGTG => 0.00016, GGTGTT => 2.2e-05, GGTTAA => -0.00019, GGTTAC => 8.2e-05, 
GGTTAG => -0.00014, GGTTAT => 1.7e-05, GGTTCA => -7.1e-05, GGTTCC => 4.1e-05, GGTTCG => 1.2e-05, GGTTCT => -3e-05, 
GGTTGA => -0.00015, GGTTGC => -7.5e-05, GGTTGG => -9.5e-05, GGTTGT => -0.00012, GGTTTA => -0.0001, GGTTTC => 5.4e-06, 
GGTTTG => -4.9e-05, GGTTTT => -0.00023, GTAAAA => -0.00028, GTAAAC => -7.1e-05, GTAAAG => 1.7e-05, GTAAAT => -0.00025, 
GTAACA => -9.2e-05, GTAACC => -9.2e-06, GTAACG => -2.7e-05, GTAACT => -7.5e-05, GTAAGA => -0.00012, GTAAGC => -6.3e-05, 
GTAAGG => -7e-05, GTAAGT => -0.00012, GTAATA => -0.00018, GTAATC => -7.1e-05, GTAATG => -6.4e-05, GTAATT => -0.00021, 
GTACAA => -8.8e-05, GTACAC => -5.2e-05, GTACAG => 4.3e-05, GTACAT => -0.00014, GTACCA => -8.9e-06, GTACCC => 2.9e-05, 
GTACCG => 8.7e-06, GTACCT => -1.7e-05, GTACGA => 1.6e-06, GTACGC => 2e-05, GTACGG => 4.9e-06, GTACGT => -2.5e-05, 
GTACTA => -5.9e-05, GTACTC => -3.2e-05, GTACTG => -2.2e-06, GTACTT => -0.00013, GTAGAA => 3.1e-05, GTAGAC => 6.8e-05, 
GTAGAG => 0.00011, GTAGAT => 5.3e-05, GTAGCA => -1.2e-05, GTAGCC => 6.2e-05, GTAGCG => -3.5e-06, GTAGCT => -3.5e-05, 
GTAGGA => -3.4e-05, GTAGGC => -1.5e-05, GTAGGG => -4.2e-05, GTAGGT => -4.9e-05, GTAGTA => -5.6e-05, GTAGTC => -2.8e-05, 
GTAGTG => 1.6e-05, GTAGTT => -0.0001, GTATAA => -0.00024, GTATAC => -5.2e-05, GTATAG => -0.00012, GTATAT => -0.00025, 
GTATCA => -5.2e-05, GTATCC => 1.6e-05, GTATCG => -1.6e-06, GTATCT => -7.4e-05, GTATGA => -0.00016, GTATGC => -6.6e-05, 
GTATGG => -7.3e-05, GTATGT => -0.00022, GTATTA => -0.00019, GTATTC => -6.2e-05, GTATTG => -0.00011, GTATTT => -0.00036, 
GTCAAA => 0.00012, GTCAAC => 0.00033, GTCAAG => 0.00044, GTCAAT => 0.00013, GTCACA => 5.5e-05, GTCACC => 0.00035, 
GTCACG => 8e-05, GTCACT => 8.5e-05, GTCAGA => -6.6e-05, GTCAGC => 0.00021, GTCAGG => -5.6e-05, GTCAGT => 3.3e-05, 
GTCATA => -3e-05, GTCATC => 0.00047, GTCATG => 0.00022, GTCATT => 0.0001, GTCCAA => -1.8e-05, GTCCAC => 6.7e-05, 
GTCCAG => 0.00028, GTCCAT => -1.8e-05, GTCCCA => -3.8e-05, GTCCCC => -1.5e-06, GTCCCG => 1e-05, GTCCCT => -2.4e-05, 
GTCCGA => 2.3e-05, GTCCGC => 7.6e-05, GTCCGG => 5.1e-05, GTCCGT => 1.4e-05, GTCCTA => -1.9e-05, GTCCTC => 8.4e-05, 
GTCCTG => 0.00033, GTCCTT => -2.8e-05, GTCGAA => 1.4e-06, GTCGAC => 8.3e-05, GTCGAG => 0.00013, GTCGAT => 7e-05, 
GTCGCA => -7.5e-06, GTCGCC => 6.8e-05, GTCGCG => -1.8e-05, GTCGCT => -5.4e-06, GTCGGA => 2e-05, GTCGGC => 7e-05, 
GTCGGG => -1.7e-05, GTCGGT => 2.9e-05, GTCGTA => -3.8e-06, GTCGTC => 6.1e-05, GTCGTG => 8.4e-05, GTCGTT => -7.5e-06, 
GTCTAA => -0.00015, GTCTAC => 0.00023, GTCTAG => -0.00012, GTCTAT => 5.2e-05, GTCTCA => -7e-05, GTCTCC => 0.00012, 
GTCTCG => 2.4e-05, GTCTCT => -7e-05, GTCTGA => -0.00024, GTCTGC => 2.7e-05, GTCTGG => -4.9e-05, GTCTGT => -0.00014, 
GTCTTA => -7.6e-05, GTCTTC => 0.00025, GTCTTG => -2.5e-05, GTCTTT => -4e-05, GTGAAA => 5.5e-05, GTGAAC => 0.0003, 
GTGAAG => 0.00069, GTGAAT => 0.00011, GTGACA => 7.5e-05, GTGACC => 0.00035, GTGACG => 0.00013, GTGACT => 8.1e-05, 
GTGAGA => -8.4e-05, GTGAGC => 0.00013, GTGAGG => -1.1e-05, GTGAGT => -2.6e-05, GTGATA => -2.5e-05, GTGATC => 0.00033, 
GTGATG => 0.0003, GTGATT => 9.2e-05, GTGCAA => 2.5e-05, GTGCAC => 0.00019, GTGCAG => 0.00056, GTGCAT => -5.4e-05, 
GTGCCA => 0.00011, GTGCCC => 0.00041, GTGCCG => 0.00016, GTGCCT => 9.3e-05, GTGCGA => 8.1e-05, GTGCGC => 0.00026, 
GTGCGG => 0.00016, GTGCGT => 4.1e-05, GTGCTA => 2.4e-05, GTGCTC => 0.00027, GTGCTG => 0.00077, GTGCTT => -1.2e-05, 
GTGGAA => 0.00041, GTGGAC => 0.00084, GTGGAG => 0.0013, GTGGAT => 0.00055, GTGGCA => 0.0002, GTGGCC => 0.00089, 
GTGGCG => 0.00019, GTGGCT => 0.00036, GTGGGA => 0.00013, GTGGGC => 0.00054, GTGGGG => 1.8e-05, GTGGGT => 0.00012, 
GTGGTA => 5.2e-05, GTGGTC => 0.0004, GTGGTG => 0.00083, GTGGTT => 0.00013, GTGTAA => -0.00026, GTGTAC => 0.00018, 
GTGTAG => -0.00017, GTGTAT => -4.6e-05, GTGTCA => 2e-06, GTGTCC => 0.00025, GTGTCG => 9.3e-05, GTGTCT => 6.9e-05, 
GTGTGA => -0.00032, GTGTGC => 2.6e-05, GTGTGG => -2e-05, GTGTGT => -0.00057, GTGTTA => -0.00012, GTGTTC => 0.00017, 
GTGTTG => -2.3e-05, GTGTTT => -7.8e-05, GTTAAA => -0.00016, GTTAAC => -1.6e-05, GTTAAG => 2e-05, GTTAAT => -0.00012, 
GTTACA => -7.7e-05, GTTACC => 3.5e-05, GTTACG => -2.3e-06, GTTACT => -4.1e-05, GTTAGA => -8.9e-05, GTTAGC => -1.8e-05, 
GTTAGG => -7e-05, GTTAGT => -9.2e-05, GTTATA => -0.00011, GTTATC => 1.7e-05, GTTATG => -4e-05, GTTATT => -0.00016, 
GTTCAA => -0.00012, GTTCAC => 1e-05, GTTCAG => 0.00015, GTTCAT => -7.5e-05, GTTCCA => 2.3e-05, GTTCCC => 3.2e-05, 
GTTCCG => 2.6e-05, GTTCCT => 3e-05, GTTCGA => -9.3e-06, GTTCGC => 4.7e-05, GTTCGG => 8.3e-06, GTTCGT => 1.9e-05, 
GTTCTA => -7.8e-05, GTTCTC => -3.8e-05, GTTCTG => 0.00011, GTTCTT => -0.00012, GTTGAA => 4.5e-05, GTTGAC => 4.5e-05, 
GTTGAG => 9.7e-05, GTTGAT => 5.6e-05, GTTGCA => -4e-05, GTTGCC => 5.8e-05, GTTGCG => -1.9e-05, GTTGCT => -5.8e-07, 
GTTGGA => 7e-05, GTTGGC => 8.6e-05, GTTGGG => -4.6e-05, GTTGGT => 1.3e-05, GTTGTA => -0.0001, GTTGTC => 1.7e-05, 
GTTGTG => 6.7e-05, GTTGTT => -0.00018, GTTTAA => -0.00041, GTTTAC => -3.2e-05, GTTTAG => -0.00021, GTTTAT => -0.00023, 
GTTTCA => -0.00014, GTTTCC => -8.4e-05, GTTTCG => -2.6e-05, GTTTCT => -0.0002, GTTTGA => -0.00034, GTTTGC => -0.00015, 
GTTTGG => -0.00017, GTTTGT => -0.00034, GTTTTA => -0.0004, GTTTTC => -0.00023, GTTTTG => -0.00027, GTTTTT => -0.00059, 
TAAAAA => -0.0011, TAAAAC => -0.00055, TAAAAG => -0.0005, TAAAAT => -0.00093, TAAACA => -0.00053, TAAACC => -0.00024, 
TAAACG => -0.00014, TAAACT => -0.00037, TAAAGA => -0.00046, TAAAGC => -0.00031, TAAAGG => -0.00031, TAAAGT => -0.00043, 
TAAATA => -0.00081, TAAATC => -0.00036, TAAATG => -0.00056, TAAATT => -0.00066, TAACAA => -0.00036, TAACAC => -0.0002, 
TAACAG => -0.00025, TAACAT => -0.00032, TAACCA => -0.00023, TAACCC => -0.00014, TAACCG => -7.1e-05, TAACCT => -0.00018, 
TAACGA => -9.3e-05, TAACGC => -6.9e-05, TAACGG => -7.7e-05, TAACGT => -9.7e-05, TAACTA => -0.00023, TAACTC => -0.00018, 
TAACTG => -0.00024, TAACTT => -0.00032, TAAGAA => -0.00037, TAAGAC => -0.00017, TAAGAG => -0.00022, TAAGAT => -0.00025, 
TAAGCA => -0.00027, TAAGCC => -0.00016, TAAGCG => -7.5e-05, TAAGCT => -0.00023, TAAGGA => -0.00022, TAAGGC => -0.00016, 
TAAGGG => -0.00015, TAAGGT => -0.00017, TAAGTA => -0.00026, TAAGTC => -0.00015, TAAGTG => -0.00026, TAAGTT => -0.0003, 
TAATAA => -0.0007, TAATAC => -0.00021, TAATAG => -0.00022, TAATAT => -0.00048, TAATCA => -0.00029, TAATCC => -0.00019, 
TAATCG => -8.4e-05, TAATCT => -0.00025, TAATGA => -0.00035, TAATGC => -0.00022, TAATGG => -0.00021, TAATGT => -0.00038, 
TAATTA => -0.00049, TAATTC => -0.00028, TAATTG => -0.00032, TAATTT => -0.00075, TACAAA => -8.4e-05, TACAAC => 0.00033, 
TACAAG => 0.00048, TACAAT => 4.9e-05, TACACA => -6.3e-05, TACACC => 0.00026, TACACG => 0.00012, TACACT => -2.5e-05, 
TACAGA => -3.2e-05, TACAGC => 0.00026, TACAGG => 2.3e-05, TACAGT => -2e-05, TACATA => -0.00025, TACATC => 0.00034, 
TACATG => 0.00023, TACATT => -0.00014, TACCAA => -2.9e-05, TACCAC => 0.00013, TACCAG => 0.00042, TACCAT => -4.7e-06, 
TACCCA => 3.8e-05, TACCCC => 7.6e-05, TACCCG => 5.7e-05, TACCCT => 2.7e-05, TACCGA => 6.7e-05, TACCGC => 0.00021, 
TACCGG => 0.00012, TACCGT => 5e-05, TACCTA => -2.4e-05, TACCTC => 0.00012, TACCTG => 0.00051, TACCTT => -5.2e-05, 
TACGAA => 7.9e-05, TACGAC => 0.00028, TACGAG => 0.00041, TACGAT => 0.00019, TACGCA => 2e-05, TACGCC => 0.00023, 
TACGCG => 5.4e-05, TACGCT => 5.3e-05, TACGGA => 9.9e-05, TACGGC => 0.00025, TACGGG => 7.7e-05, TACGGT => 5.5e-05, 
TACGTA => -2e-05, TACGTC => 0.0001, TACGTG => 0.00021, TACGTT => -1.6e-05, TACTAA => -0.0002, TACTAC => 0.00027, 
TACTAG => -9.2e-05, TACTAT => 3.5e-05, TACTCA => -4.6e-05, TACTCC => 0.00017, TACTCG => 8.9e-05, TACTCT => -2.2e-05, 
TACTGA => -0.00022, TACTGC => 0.00013, TACTGG => 7.8e-05, TACTGT => -9e-05, TACTTA => -0.00017, TACTTC => 0.00028, 
TACTTG => 3.8e-06, TACTTT => -0.00011, TAGAAA => -0.00043, TAGAAC => -0.00016, TAGAAG => -0.00026, TAGAAT => -0.00027, 
TAGACA => -0.0002, TAGACC => -0.00011, TAGACG => -6.8e-05, TAGACT => -0.00016, TAGAGA => -0.00027, TAGAGC => -0.00018, 
TAGAGG => -0.00018, TAGAGT => -0.00018, TAGATA => -0.0002, TAGATC => -0.00013, TAGATG => -0.00023, TAGATT => -0.00025, 
TAGCAA => -0.00023, TAGCAC => -0.00015, TAGCAG => -0.0002, TAGCAT => -0.00021, TAGCCA => -0.00021, TAGCCC => -0.00013, 
TAGCCG => -6.8e-05, TAGCCT => -0.00017, TAGCGA => -6.6e-05, TAGCGC => -5.9e-05, TAGCGG => -5.6e-05, TAGCGT => -5.8e-05, 
TAGCTA => -0.00016, TAGCTC => -0.00017, TAGCTG => -0.00025, TAGCTT => -0.00022, TAGGAA => -0.00025, TAGGAC => -0.00012, 
TAGGAG => -0.00017, TAGGAT => -0.00016, TAGGCA => -0.00016, TAGGCC => -0.0001, TAGGCG => -4.9e-05, TAGGCT => -0.00016, 
TAGGGA => -0.00015, TAGGGC => -0.00011, TAGGGG => -0.00011, TAGGGT => -0.00013, TAGGTA => -0.00014, TAGGTC => -9.7e-05, 
TAGGTG => -0.00015, TAGGTT => -0.00016, TAGTAA => -0.00021, TAGTAC => -0.00012, TAGTAG => -0.00015, TAGTAT => -0.00019, 
TAGTCA => -0.00016, TAGTCC => -0.00013, TAGTCG => -5.9e-05, TAGTCT => -0.00017, TAGTGA => -0.00019, TAGTGC => -0.00014, 
TAGTGG => -0.00015, TAGTGT => -0.00022, TAGTTA => -0.00022, TAGTTC => -0.00018, TAGTTG => -0.0002, TAGTTT => -0.00042, 
TATAAA => -0.00041, TATAAC => -4.2e-05, TATAAG => -1e-05, TATAAT => -0.00026, TATACA => -0.0003, TATACC => -2.6e-05, 
TATACG => -7e-06, TATACT => -0.00011, TATAGA => -0.00018, TATAGC => -5.9e-05, TATAGG => -8.1e-05, TATAGT => -0.00015, 
TATATA => -0.00078, TATATC => -6.6e-05, TATATG => -0.00022, TATATT => -0.00046, TATCAA => -0.00013, TATCAC => 1.6e-06, 
TATCAG => 0.00012, TATCAT => -0.00011, TATCCA => -3.1e-06, TATCCC => 4.5e-05, TATCCG => 6.4e-05, TATCCT => -7.3e-06, 
TATCGA => -3.9e-06, TATCGC => 7.7e-05, TATCGG => 3e-05, TATCGT => -1.5e-06, TATCTA => -0.00014, TATCTC => -4.1e-05, 
TATCTG => 0.00014, TATCTT => -0.00015, TATGAA => 0.00012, TATGAC => 0.00028, TATGAG => 0.00037, TATGAT => 0.00012, 
TATGCA => -2.6e-05, TATGCC => 0.00024, TATGCG => 5.2e-05, TATGCT => 2.8e-05, TATGGA => 0.00013, TATGGC => 0.00019, 
TATGGG => 5.4e-05, TATGGT => 1.1e-05, TATGTA => -0.00031, TATGTC => 4.4e-05, TATGTG => 0.00018, TATGTT => -0.00017, 
TATTAA => -0.00055, TATTAC => -4.7e-05, TATTAG => -0.00022, TATTAT => -0.00042, TATTCA => -0.00023, TATTCC => -6.7e-05, 
TATTCG => -4.3e-05, TATTCT => -0.0002, TATTGA => -0.00029, TATTGC => -0.00013, TATTGG => -0.00012, TATTGT => -0.00035, 
TATTTA => -0.00069, TATTTC => -0.0002, TATTTG => -0.00031, TATTTT => -0.001, TCAAAA => -0.00035, TCAAAC => -0.00011, 
TCAAAG => -7e-05, TCAAAT => -0.00023, TCAACA => -0.00014, TCAACC => -3.9e-05, TCAACG => -1e-05, TCAACT => -0.00013, 
TCAAGA => -0.00018, TCAAGC => -7.9e-05, TCAAGG => -0.00011, TCAAGT => -0.00012, TCAATA => -0.0002, TCAATC => -9.9e-05, 
TCAATG => -9.9e-05, TCAATT => -0.00023, TCACAA => -0.00015, TCACAC => -0.00016, TCACAG => -1.6e-05, TCACAT => -0.0002, 
TCACCA => -5.8e-05, TCACCC => -2.5e-06, TCACCG => 1.2e-05, TCACCT => -7.2e-05, TCACGA => -1.8e-05, TCACGC => -5.6e-06, 
TCACGG => -1.7e-05, TCACGT => -4.4e-05, TCACTA => -9.4e-05, TCACTC => -0.0001, TCACTG => -4.2e-05, TCACTT => -0.0002, 
TCAGAA => 9.8e-06, TCAGAC => 0.00014, TCAGAG => 0.00021, TCAGAT => 8.4e-05, TCAGCA => -7.2e-05, TCAGCC => 4.1e-05, 
TCAGCG => -1.5e-05, TCAGCT => -3.3e-05, TCAGGA => -1.8e-05, TCAGGC => 5.8e-06, TCAGGG => -2.5e-05, TCAGGT => -6.2e-05, 
TCAGTA => -0.0001, TCAGTC => -8.2e-05, TCAGTG => -3.9e-06, TCAGTT => -0.00021, TCATAA => -0.00029, TCATAC => -6.2e-05, 
TCATAG => -0.00017, TCATAT => -0.00022, TCATCA => -0.00012, TCATCC => 1.4e-05, TCATCG => 2.1e-05, TCATCT => -0.00011, 
TCATGA => -0.00024, TCATGC => -0.00013, TCATGG => -0.00012, TCATGT => -0.00023, TCATTA => -0.00024, TCATTC => -0.00014, 
TCATTG => -0.00016, TCATTT => -0.00044, TCCAAA => 8.3e-05, TCCAAC => 0.00028, TCCAAG => 0.00046, TCCAAT => 0.00013, 
TCCACA => 1.6e-05, TCCACC => 0.00026, TCCACG => 0.00016, TCCACT => 2.5e-05, TCCAGA => -0.00011, TCCAGC => 0.00029, 
TCCAGG => -7.6e-05, TCCAGT => 0.00011, TCCATA => -6.5e-05, TCCATC => 0.00029, TCCATG => 0.00025, TCCATT => 3.2e-05, 
TCCCAA => -0.00011, TCCCAC => -1.8e-05, TCCCAG => 0.00016, TCCCAT => -7e-05, TCCCCA => -3.6e-05, TCCCCC => -5e-05, 
TCCCCG => 3.7e-05, TCCCCT => -5.7e-05, TCCCGA => 6.1e-06, TCCCGC => 9.9e-05, TCCCGG => 4.2e-05, TCCCGT => 2e-05, 
TCCCTA => -5.6e-05, TCCCTC => -4.7e-05, TCCCTG => 0.00026, TCCCTT => -0.00013, TCCGAA => 3.2e-05, TCCGAC => 0.00014, 
TCCGAG => 0.00023, TCCGAT => 0.0001, TCCGCA => 1.5e-05, TCCGCC => 0.00014, TCCGCG => 8.8e-06, TCCGCT => 2.5e-05, 
TCCGGA => 5.7e-05, TCCGGC => 0.00013, TCCGGG => 8.4e-06, TCCGGT => 4.7e-05, TCCGTA => -1.2e-05, TCCGTC => 4.7e-05, 
TCCGTG => 0.00017, TCCGTT => 9.9e-07, TCCTAA => -0.0002, TCCTAC => 0.00022, TCCTAG => -0.00014, TCCTAT => 5.2e-05, 
TCCTCA => -4.8e-06, TCCTCC => 0.00013, TCCTCG => 0.00012, TCCTCT => -6.2e-05, TCCTGA => -0.00034, TCCTGC => -6.7e-05, 
TCCTGG => -8.6e-05, TCCTGT => -0.00014, TCCTTA => -8.8e-05, TCCTTC => 0.00013, TCCTTG => 1.3e-05, TCCTTT => -0.00015, 
TCGAAA => -8.7e-05, TCGAAC => 6.4e-06, TCGAAG => 5.6e-05, TCGAAT => -1.9e-05, TCGACA => -2e-05, TCGACC => 2.7e-05, 
TCGACG => 2.7e-05, TCGACT => -2.5e-05, TCGAGA => -5.5e-05, TCGAGC => 1.7e-05, TCGAGG => -2.8e-05, TCGAGT => -3.7e-05, 
TCGATA => -4e-05, TCGATC => -1.4e-06, TCGATG => 4.3e-05, TCGATT => -5.1e-05, TCGCAA => -1.1e-06, TCGCAC => 5e-05, 
TCGCAG => 0.00019, TCGCAT => -5.8e-06, TCGCCA => 6.1e-05, TCGCCC => 0.00015, TCGCCG => 0.0001, TCGCCT => 1.7e-06, 
TCGCGA => 8.2e-06, TCGCGC => 4.7e-05, TCGCGG => 2.4e-06, TCGCGT => -1.4e-05, TCGCTA => -5.6e-07, TCGCTC => 2.4e-05, 
TCGCTG => 0.00029, TCGCTT => -6.2e-05, TCGGAA => 7.7e-05, TCGGAC => 0.00016, TCGGAG => 0.00028, TCGGAT => 0.00012, 
TCGGCA => 4.9e-05, TCGGCC => 0.00016, TCGGCG => 0.00011, TCGGCT => 4.7e-05, TCGGGA => 5.7e-05, TCGGGC => 0.0001, 
TCGGGG => 1.1e-05, TCGGGT => 2.9e-05, TCGGTA => -8.4e-06, TCGGTC => 2.9e-05, TCGGTG => 0.00022, TCGGTT => -2.9e-05, 
TCGTAA => -8.8e-05, TCGTAC => 5.4e-05, TCGTAG => -5.5e-05, TCGTAT => -1.2e-05, TCGTCA => -5.8e-06, TCGTCC => 4.6e-05, 
TCGTCG => 5.7e-05, TCGTCT => -1.3e-05, TCGTGA => -7.8e-05, TCGTGC => -1.8e-05, TCGTGG => -8e-06, TCGTGT => -5.9e-05, 
TCGTTA => -6e-05, TCGTTC => 1e-05, TCGTTG => 1.4e-06, TCGTTT => -0.0001, TCTAAA => -0.00015, TCTAAC => -3e-05, 
TCTAAG => 1.2e-05, TCTAAT => -0.00011, TCTACA => -7.3e-05, TCTACC => -1.8e-05, TCTACG => -1.2e-05, TCTACT => -6.7e-05, 
TCTAGA => -0.00011, TCTAGC => -1.6e-05, TCTAGG => -8.9e-05, TCTAGT => -6.2e-05, TCTATA => -0.00015, TCTATC => -4.8e-05, 
TCTATG => -4.3e-05, TCTATT => -0.00016, TCTCAA => -0.00011, TCTCAC => -7.5e-05, TCTCAG => 9.9e-05, TCTCAT => -0.00014, 
TCTCCA => 7.9e-06, TCTCCC => -7.6e-05, TCTCCG => 5.4e-06, TCTCCT => -6.7e-05, TCTCGA => 1.1e-05, TCTCGC => 1.7e-05, 
TCTCGG => 6.2e-06, TCTCGT => -2.6e-05, TCTCTA => -0.00013, TCTCTC => -0.00027, TCTCTG => 8.1e-05, TCTCTT => -0.00024, 
TCTGAA => 4.4e-05, TCTGAC => 0.00012, TCTGAG => 0.00023, TCTGAT => 9.2e-05, TCTGCA => -2.7e-05, TCTGCC => 5.6e-05, 
TCTGCG => -8.3e-06, TCTGCT => -2.6e-05, TCTGGA => 5.3e-05, TCTGGC => 6.6e-05, TCTGGG => -5.9e-05, TCTGGT => -1.5e-05, 
TCTGTA => -0.00018, TCTGTC => -8.8e-05, TCTGTG => 9.1e-05, TCTGTT => -0.0002, TCTTAA => -0.00037, TCTTAC => -4.4e-05, 
TCTTAG => -0.00019, TCTTAT => -0.00015, TCTTCA => -0.00015, TCTTCC => -9.3e-05, TCTTCG => 2.9e-06, TCTTCT => -0.00016, 
TCTTGA => -0.00029, TCTTGC => -0.00015, TCTTGG => -0.00017, TCTTGT => -0.00025, TCTTTA => -0.00028, TCTTTC => -0.00022, 
TCTTTG => -0.00018, TCTTTT => -0.0005, TGAAAA => -0.00075, TGAAAC => -0.00037, TGAAAG => -0.00042, TGAAAT => -0.00061, 
TGAACA => -0.00039, TGAACC => -0.00021, TGAACG => -0.00012, TGAACT => -0.00036, TGAAGA => -0.00053, TGAAGC => -0.00033, 
TGAAGG => -0.00034, TGAAGT => -0.00036, TGAATA => -0.00039, TGAATC => -0.00025, TGAATG => -0.00043, TGAATT => -0.00046, 
TGACAA => -0.00032, TGACAC => -0.00021, TGACAG => -0.00032, TGACAT => -0.00031, TGACCA => -0.00026, TGACCC => -0.0002, 
TGACCG => -7.8e-05, TGACCT => -0.00027, TGACGA => -0.0001, TGACGC => -9.5e-05, TGACGG => -8.5e-05, TGACGT => -0.00011, 
TGACTA => -0.00019, TGACTC => -0.00023, TGACTG => -0.00033, TGACTT => -0.00035, TGAGAA => -0.00043, TGAGAC => -0.00026, 
TGAGAG => -0.00033, TGAGAT => -0.00031, TGAGCA => -0.00032, TGAGCC => -0.00027, TGAGCG => -0.00012, TGAGCT => -0.00032, 
TGAGGA => -0.00037, TGAGGC => -0.00026, TGAGGG => -0.00025, TGAGGT => -0.00024, TGAGTA => -0.00022, TGAGTC => -0.0002, 
TGAGTG => -0.0003, TGAGTT => -0.00033, TGATAA => -0.00032, TGATAC => -0.00016, TGATAG => -0.00017, TGATAT => -0.00031, 
TGATCA => -0.00025, TGATCC => -0.00019, TGATCG => -6.9e-05, TGATCT => -0.00026, TGATGA => -0.0004, TGATGC => -0.00024, 
TGATGG => -0.00029, TGATGT => -0.00035, TGATTA => -0.00031, TGATTC => -0.00026, TGATTG => -0.00028, TGATTT => -0.00059, 
TGCAAA => -0.00018, TGCAAC => 0.00011, TGCAAG => 0.00019, TGCAAT => -8.6e-05, TGCACA => -0.00016, TGCACC => 8.9e-05, 
TGCACG => 2.3e-05, TGCACT => -0.00012, TGCAGA => -0.00025, TGCAGC => 2e-07, TGCAGG => -0.00016, TGCAGT => -0.00016, 
TGCATA => -0.00018, TGCATC => 0.00012, TGCATG => -2e-05, TGCATT => -0.00022, TGCCAA => -0.00016, TGCCAC => -4.1e-05, 
TGCCAG => 0.00015, TGCCAT => -0.00013, TGCCCA => -0.00013, TGCCCC => -2e-05, TGCCCG => -9.1e-07, TGCCCT => -0.00013, 
TGCCGA => -2.2e-05, TGCCGC => 7.6e-05, TGCCGG => 5.4e-05, TGCCGT => -3.1e-05, TGCCTA => -9.2e-05, TGCCTC => -9.4e-05, 
TGCCTG => 8.4e-05, TGCCTT => -0.00022, TGCGAA => 3.8e-06, TGCGAC => 0.0001, TGCGAG => 0.00019, TGCGAT => 5.7e-05, 
TGCGCA => -6.5e-05, TGCGCC => 7.1e-05, TGCGCG => -3.7e-05, TGCGCT => -3.9e-05, TGCGGA => 1.5e-05, TGCGGC => 9.6e-05, 
TGCGGG => 3.9e-07, TGCGGT => -1.4e-05, TGCGTA => -4.3e-05, TGCGTC => 2.4e-05, TGCGTG => 5.2e-05, TGCGTT => -7.1e-05, 
TGCTAA => -0.00026, TGCTAC => 8.5e-05, TGCTAG => -0.00015, TGCTAT => -5.8e-05, TGCTCA => -0.00015, TGCTCC => 2.6e-05, 
TGCTCG => -1.1e-05, TGCTCT => -0.0002, TGCTGA => -0.00039, TGCTGC => -0.00014, TGCTGG => -0.00022, TGCTGT => -0.00027, 
TGCTTA => -0.00019, TGCTTC => 2e-05, TGCTTG => -0.00013, TGCTTT => -0.00034, TGGAAA => -0.00027, TGGAAC => 2.9e-05, 
TGGAAG => 6.2e-05, TGGAAT => -0.00011, TGGACA => -0.00012, TGGACC => 4.5e-05, TGGACG => -6.6e-06, TGGACT => -0.00013, 
TGGAGA => -0.00029, TGGAGC => -0.0001, TGGAGG => -0.0002, TGGAGT => -0.00016, TGGATA => -0.00016, TGGATC => 3.2e-05, 
TGGATG => -1.7e-05, TGGATT => -0.00016, TGGCAA => -0.00019, TGGCAC => -5.9e-05, TGGCAG => 9e-06, TGGCAT => -0.00014, 
TGGCCA => -0.00022, TGGCCC => -0.0001, TGGCCG => -5.7e-05, TGGCCT => -0.0002, TGGCGA => -3.4e-05, TGGCGC => 2e-05, 
TGGCGG => -3.3e-05, TGGCGT => -3.4e-05, TGGCTA => -9.5e-05, TGGCTC => -6.6e-05, TGGCTG => 9.3e-05, TGGCTT => -0.00022, 
TGGGAA => -0.00012, TGGGAC => 9.8e-05, TGGGAG => 9.7e-05, TGGGAT => 1.6e-05, TGGGCA => -0.00011, TGGGCC => 2.6e-05, 
TGGGCG => -2.8e-05, TGGGCT => -9.6e-05, TGGGGA => -0.00015, TGGGGC => -4.1e-06, TGGGGG => -0.00018, TGGGGT => -0.0001, 
TGGGTA => -9.3e-05, TGGGTC => -3.6e-05, TGGGTG => 2.3e-05, TGGGTT => -0.00017, TGGTAA => -0.00022, TGGTAC => 6.7e-05, 
TGGTAG => -0.00015, TGGTAT => -6.2e-05, TGGTCA => -0.00013, TGGTCC => -3.6e-05, TGGTCG => -3.6e-06, TGGTCT => -0.00013, 
TGGTGA => -0.00028, TGGTGC => -0.0001, TGGTGG => -0.00014, TGGTGT => -0.0002, TGGTTA => -0.00017, TGGTTC => 1.8e-05, 
TGGTTG => -0.00015, TGGTTT => -0.00024, TGTAAA => -0.0004, TGTAAC => -9.5e-05, TGTAAG => -7.8e-05, TGTAAT => -0.00029, 
TGTACA => -0.00025, TGTACC => -6.9e-05, TGTACG => -4.9e-05, TGTACT => -0.00017, TGTAGA => -0.00021, TGTAGC => -0.00013, 
TGTAGG => -0.00014, TGTAGT => -0.00017, TGTATA => -0.00037, TGTATC => -0.0001, TGTATG => -0.00022, TGTATT => -0.00042, 
TGTCAA => -0.00017, TGTCAC => -0.00013, TGTCAG => -4.4e-05, TGTCAT => -0.00022, TGTCCA => -0.00013, TGTCCC => -6.9e-05, 
TGTCCG => -1.1e-05, TGTCCT => -0.00015, TGTCGA => -2.4e-05, TGTCGC => 4.9e-06, TGTCGG => -5.8e-06, TGTCGT => -4.9e-05, 
TGTCTA => -0.00015, TGTCTC => -0.00018, TGTCTG => -0.00012, TGTCTT => -0.00029, TGTGAA => -0.00011, TGTGAC => 9.1e-05, 
TGTGAG => 0.00011, TGTGAT => -5.6e-05, TGTGCA => -0.00019, TGTGCC => 6.6e-05, TGTGCG => -7.5e-05, TGTGCT => -0.00013, 
TGTGGA => -5.9e-05, TGTGGC => 4.6e-05, TGTGGG => -2.3e-05, TGTGGT => -0.00013, TGTGTA => -0.00031, TGTGTC => -0.00011, 
TGTGTG => -0.00054, TGTGTT => -0.00041, TGTTAA => -0.00041, TGTTAC => -0.00011, TGTTAG => -0.00021, TGTTAT => -0.00028, 
TGTTCA => -0.00029, TGTTCC => -0.00013, TGTTCG => -6.1e-05, TGTTCT => -0.0003, TGTTGA => -0.00037, TGTTGC => -0.00023, 
TGTTGG => -0.00027, TGTTGT => -0.00041, TGTTTA => -0.0005, TGTTTC => -0.00029, TGTTTG => -0.00042, TGTTTT => -0.00096, 
TTAAAA => -0.0009, TTAAAC => -0.00028, TTAAAG => -0.00028, TTAAAT => -0.00063, TTAACA => -0.00026, TTAACC => -0.00011, 
TTAACG => -7.4e-05, TTAACT => -0.00023, TTAAGA => -0.00025, TTAAGC => -0.00017, TTAAGG => -0.00016, TTAAGT => -0.00027, 
TTAATA => -0.00047, TTAATC => -0.00019, TTAATG => -0.00026, TTAATT => -0.00057, TTACAA => -0.00026, TTACAC => -0.00015, 
TTACAG => -9.6e-05, TTACAT => -0.0003, TTACCA => -0.0001, TTACCC => -6.3e-05, TTACCG => -2.7e-05, TTACCT => -0.00011, 
TTACGA => -4.8e-05, TTACGC => -2e-05, TTACGG => -1.5e-05, TTACGT => -7.1e-05, TTACTA => -0.00016, TTACTC => -0.00011, 
TTACTG => -0.0001, TTACTT => -0.00027, TTAGAA => -8.4e-05, TTAGAC => -1.4e-05, TTAGAG => -5.7e-06, TTAGAT => -6.3e-05, 
TTAGCA => -0.00011, TTAGCC => -4.3e-05, TTAGCG => -3.7e-05, TTAGCT => -0.00012, TTAGGA => -0.0001, TTAGGC => -6.4e-05, 
TTAGGG => -8.9e-05, TTAGGT => -9.5e-05, TTAGTA => -0.00017, TTAGTC => -9.4e-05, TTAGTG => -0.0001, TTAGTT => -0.00029, 
TTATAA => -0.00049, TTATAC => -0.00018, TTATAG => -0.00025, TTATAT => -0.00051, TTATCA => -0.00021, TTATCC => -9e-05, 
TTATCG => -4.8e-05, TTATCT => -0.00019, TTATGA => -0.00032, TTATGC => -0.00016, TTATGG => -0.00016, TTATGT => -0.00035, 
TTATTA => -0.00057, TTATTC => -0.00025, TTATTG => -0.00032, TTATTT => -0.001, TTCAAA => -0.00019, TTCAAC => 0.00033, 
TTCAAG => 0.00041, TTCAAT => -3.7e-06, TTCACA => -5.5e-05, TTCACC => 0.00034, TTCACG => 8.8e-05, TTCACT => -1.4e-05, 
TTCAGA => -0.00018, TTCAGC => 0.00024, TTCAGG => -6.8e-05, TTCAGT => -5e-05, TTCATA => -0.00019, TTCATC => 0.0004, 
TTCATG => 0.00018, TTCATT => -0.00017, TTCCAA => -0.00012, TTCCAC => 0.00011, TTCCAG => 0.00046, TTCCAT => -0.00011, 
TTCCCA => -0.00013, TTCCCC => 4.6e-05, TTCCCG => 2.3e-05, TTCCCT => -9.3e-05, TTCCGA => 5.8e-05, TTCCGC => 0.00022, 
TTCCGG => 0.00016, TTCCGT => 4.9e-05, TTCCTA => -7.8e-05, TTCCTC => 0.00015, TTCCTG => 0.00055, TTCCTT => -0.00021, 
TTCGAA => -2e-05, TTCGAC => 0.00023, TTCGAG => 0.00032, TTCGAT => 0.00017, TTCGCA => -8.4e-06, TTCGCC => 0.00025, 
TTCGCG => 7.5e-06, TTCGCT => 1.1e-05, TTCGGA => 5.7e-05, TTCGGC => 0.00018, TTCGGG => 6.7e-05, TTCGGT => 1.3e-05, 
TTCGTA => -3.2e-05, TTCGTC => 9.7e-05, TTCGTG => 0.00022, TTCGTT => -6.5e-05, TTCTAA => -0.00034, TTCTAC => 0.0003, 
TTCTAG => -0.00022, TTCTAT => -2.1e-05, TTCTCA => -0.00015, TTCTCC => 0.00016, TTCTCG => 4.1e-05, TTCTCT => -0.00018, 
TTCTGA => -0.0004, TTCTGC => 4e-06, TTCTGG => -3.7e-05, TTCTGT => -0.00029, TTCTTA => -0.00025, TTCTTC => 0.00027, 
TTCTTG => -0.0001, TTCTTT => -0.00036, TTGAAA => -0.00026, TTGAAC => -5e-05, TTGAAG => 7.3e-05, TTGAAT => -0.00019, 
TTGACA => -0.00013, TTGACC => 1.8e-05, TTGACG => -7.4e-06, TTGACT => -0.00012, TTGAGA => -0.0002, TTGAGC => -2.4e-05, 
TTGAGG => -0.0001, TTGAGT => -0.00014, TTGATA => -0.00019, TTGATC => -1e-05, TTGATG => -9.4e-05, TTGATT => -0.00027, 
TTGCAA => -0.00014, TTGCAC => -5.3e-05, TTGCAG => 7.7e-05, TTGCAT => -0.00021, TTGCCA => -5.8e-05, TTGCCC => 4e-05, 
TTGCCG => 2e-05, TTGCCT => -0.00011, TTGCGA => -1.2e-05, TTGCGC => 4.1e-05, TTGCGG => 2e-05, TTGCGT => 3.2e-06, 
TTGCTA => -0.00012, TTGCTC => -6.3e-05, TTGCTG => 9.5e-05, TTGCTT => -0.00025, TTGGAA => 8.7e-05, TTGGAC => 0.00021, 
TTGGAG => 0.00035, TTGGAT => 0.00021, TTGGCA => -2.7e-05, TTGGCC => 0.00022, TTGGCG => 4.7e-05, TTGGCT => 4.4e-05, 
TTGGGA => -6.8e-05, TTGGGC => 0.00012, TTGGGG => -0.00011, TTGGGT => -1.5e-05, TTGGTA => -8.5e-05, TTGGTC => 1.3e-05, 
TTGGTG => 0.00012, TTGGTT => -0.00018, TTGTAA => -0.00047, TTGTAC => -7.3e-05, TTGTAG => -0.00024, TTGTAT => -0.00032, 
TTGTCA => -0.00017, TTGTCC => -3.1e-05, TTGTCG => -4.8e-06, TTGTCT => -0.00016, TTGTGA => -0.00035, TTGTGC => -0.00017, 
TTGTGG => -0.00018, TTGTGT => -0.00039, TTGTTA => -0.0003, TTGTTC => -0.00019, TTGTTG => -0.00029, TTGTTT => -0.00079, 
TTTAAA => -0.00089, TTTAAC => -0.00014, TTTAAG => -0.00012, TTTAAT => -0.00055, TTTACA => -0.00036, TTTACC => -7.3e-05, 
TTTACG => -2.5e-05, TTTACT => -0.00028, TTTAGA => -0.00026, TTTAGC => -0.00011, TTTAGG => -0.00018, TTTAGT => -0.00027, 
TTTATA => -0.00056, TTTATC => -0.00013, TTTATG => -0.00027, TTTATT => -0.0009, TTTCAA => -0.00037, TTTCAC => -0.00018, 
TTTCAG => -6.4e-05, TTTCAT => -0.00036, TTTCCA => -0.00027, TTTCCC => -0.00018, TTTCCG => -3.6e-05, TTTCCT => -0.00029, 
TTTCGA => -4.8e-05, TTTCGC => -1.9e-05, TTTCGG => -3.2e-05, TTTCGT => -7.7e-05, TTTCTA => -0.0003, TTTCTC => -0.00026, 
TTTCTG => -0.00016, TTTCTT => -0.00063, TTTGAA => -1.5e-05, TTTGAC => 0.00026, TTTGAG => 0.00037, TTTGAT => 0.00013, 
TTTGCA => -0.00012, TTTGCC => 0.00022, TTTGCG => 1.5e-05, TTTGCT => -5.2e-05, TTTGGA => 5.6e-05, TTTGGC => 0.00017, 
TTTGGG => -4.5e-05, TTTGGT => -5.2e-05, TTTGTA => -0.00041, TTTGTC => -3.3e-05, TTTGTG => 0.00018, TTTGTT => -0.00062, 
TTTTAA => -0.0012, TTTTAC => -0.00025, TTTTAG => -0.00048, TTTTAT => -0.00085, TTTTCA => -0.00052, TTTTCC => -0.00036, 
TTTTCG => -0.00011, TTTTCT => -0.00068, TTTTGA => -0.00063, TTTTGC => -0.00032, TTTTGG => -0.00037, TTTTGT => -0.00082, 
TTTTTA => -0.0011, TTTTTC => -0.00057, TTTTTG => -0.00065, TTTTTT => -0.0022, );
    
}

1;
