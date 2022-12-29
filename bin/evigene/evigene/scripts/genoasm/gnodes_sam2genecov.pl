#!/usr/bin/perl
# gnodes_sam2genecov.p

=item about gnodes_sam2genecov.pl 

  expansion of gnodes_sam2ucgcov.pl
  ucg: special case of cov depth calc for uniq conserved genes, busco etc;
  C.UCG = Lread * Nread / Width_gene
  with filters to reject genes w/ duplicate reads, and reads with low map qual
  any geneset expansion: keep much of ucgcov, clean up and use for all (cds|mrna) x reads align sam
  .. esp. to measure genes copynum from dna-reads x cds/mrna seq

=item algorithm needs doc

   .. counts read align spans for gene seqs, similar in concept to samtools pileup, but
      adds read duplicate mapping details not in samtools.
   .. simple at basis, complex in some details

=item UPD21DEC:  test gene family read-share matrix
  .. want table of geneid x geneid reads shared, ie multimaps of dupl genes
  .. count how?  grdmat{$topid}{$altid}{nrd}++, {t}{a}{aln}+=$aln ?
  
  .. output separate from genestable
   rows of: geneid glen gnread? @otherids list, w/ shared-read count, cov-align-span/depth? 
  
=cut

use strict;
use Getopt::Long;  

use constant UPD21JUN => 1; #  
use constant UPD21DEC => 1; #  test gene family read-share matrix
use constant UPD22MAR => 1; #  pairmap, same read ids fwd/rev, fix revid

my $debug= $ENV{debug}||1;
my $MINE = $ENV{mine}|| 999; # NM:i:\d+ err, use 2..9, data dependent? dont want?
my $RDLEN=$ENV{rdlen}|| 150; # drop for measure sam col[9] = read
my $MINRDLEN=$ENV{minrdlen}|| 50; #? mainly for variable len w/ tiny read frags in some SRR
my $pMAXDUP= 0.05; #?? need tests

my $ALN_SUB_NMI = 1; #drop: (defined $ENV{nonmi}) ? 0 : 1; # UPD21NOV: TEST for hetero adjust, dont subtract NM:i:nn  from align val
  # ^^ test w/ arath F1 hetz gDNA vs homo gDNA shows no differ effect, both drop MapErr from ~ 1 to 0.1,
  # w/ small increase in genexcopy uniq Cm cover depth vals (eg 1% increase)
  # NOT USEFUL, change to test gene.Merr == NMI + indels as filter for C.UCG  calc
#UPD21NOV: add opts for C.ucg filters: kFullAlignFilter => 16, kUerrFilter => 32, others?
my($optFullAlignFilter,$opt1000filter,$optUerrFilter)=(1, 1, 0); # default fullalignfilter=on, others?
   
#? use  MIN_IDENT as per sam2covtab ? MINALN = MIN_IDENT * RDLEN
my $MIN_IDENT = 0.40; # UPD7e was 0.65; # lo is best? << keep same as gnodes_sam2covtab.pl, or lower for gene x DNA align deficit
my $MIN_DUPIDENT = 0.98; #UPD7e was 0.98; # was .99/1; hi is best? or 1.0; # == ident equal to 1st/top align, lower if desired
   ## gnodes_sam2covtab7e now has  MIN_DUPIDENT=0.999; lower that?, raise this? or what?

my $DOHETZ=0; #UPD19mar, test heterozygosity measure/base/gene, looks not accurate enough
my $HOMIN=0.80; 
my $UCGenePatt='UCG|busco';

use constant { SAMFLAG_nomap => 0x04, SAMFLAG_rev => 0x10, SAMFLAG_read2 => 0x80, SAMFLAG_2ndary => 0x100, SAMFLAG_suppl => 0x800 }; 
my $isPAIRMAP=0; #UPD22MAR

my ($crclassf,$outtab,$inbam,$dosumtab)=("") x 9; 
my @inbam; # opt to allow many -bam inputs, -bam one.bam two.bam three.bam ... or .sam ;  

my $optok= GetOptions( 
  'bam=s', \@inbam, # bam or sam or STDIN ? for samtools view call, add opt many @inbam
  'output|covtab=s',\$outtab, 
  'crclassf|idclassf=s',\$crclassf, # alt table chr => class for savereadids
  'maxerr=s',\$MINE, #?? reuse for NM:i:(err) filter # replace w/ minident ?
  'minident=s',\$MIN_IDENT, #? only used to set MINALN ?
  'maxdup=s', \$pMAXDUP,  #? reuse MIN_DUPIDENT for dup reads/total reads filter |mindupident
  'mindupident=s', \$MIN_DUPIDENT,  # as per sam2covtab; pMAXDUP = 1 - $MIN_DUPIDENT ??
  'minrdlen|minreadlen=i',\$MINRDLEN, # add .. not same as MINALN ??
  'rdlen|readlen=i',\$RDLEN, #  drop
  'UCGeneAnnot=s', \$UCGenePatt, # 
  'filterFullAlign!', \$optFullAlignFilter, # -nofilterFull to turn off
  'filter1000bases!', \$opt1000filter, # -nofilterFull to turn off
  'filterUerr!', \$optUerrFilter, 
  'heterozygosity!', \$DOHETZ, #? drop this
  'HOMIN=s', \$HOMIN, # for DOHETZ
  'summary|sumtab!', \$dosumtab, 
#  'PAIRMAP!', \$isPAIRMAP, #UPD22MAR
  'debug!', \$debug, 
  );

if(@ARGV) { push @inbam, grep(/\.(bam|sam)$/, @ARGV); }
$inbam= shift @inbam if(@inbam and not $inbam);

# my $optinok= (($inbam and -s $inbam) or (defined $MERGE and ($outtab or @ARGV)) );
unless($optok) {
  die "usage: gnodes_sam2genecov.pl -bam uniq_conserved_genes_readmap.bam -output uniq_conserved_genes.covtab
    or: gnodes_sam2genecov.pl < uniq_conserved_genes_readmap.sam > uniq_conserved_genes.covtab
    opts: -idclass genes.idclass : gene id table w/ UCG|busco class flag
    -maxdup=$pMAXDUP portion for uniq class, 
  ";
  # -maxerr=$MINE : max read map err, 
}


my $HOMHET=0.98 - $HOMIN; # DOHETZ options? 
# read map stats: n_mapok == nmapid? sam2cov uses n_mapok for all maps > nmapid
#  nr == n_mapok
#x my $MINALN= $MIN_IDENT * $RDLEN; # but RDLEN can be variable, set per read

my($ngene,  $nmapid, $stw, $scm, $srdlen, $srdlen_nomap, $lrid)= (0) x 9; # $nreadid, $nr
my($n_readid, $n_nomap, $n_mapok, $n_mapbad, $n_dupbad, $n_mapnotucg, $n_rdtooshort, # other: $n_dupbad, $n_partb,  
   $n_intron, $n_insert, $n_delete, $n_softclip, $n_mismatch)= (0) x 19; # globals?
my(%tw, %nrd, %nrdlen, %loc, %locdup, %cmd);
my(%rdshare,%genetab); #UPD21DEC

sub MAIN_stub {}

if($dosumtab){
  ucgcov_sumtab(); 
  exit;
}  

#? allow many -bam inputs, -bam one.bam two.bam three.bam ... or .sam ; does samtools view handle input list? no
my $inh = *STDIN; if($inbam) { open($inh,"samtools view -h $inbam |") or die "FAIL samtools view $inbam"; }
my $outh= *STDOUT; if($outtab) { rename($outtab,"$outtab.old") if(-s $outtab); open($outh, '>', $outtab) or die "writing $outtab"; }

my $topline= "#gnodes_sam2genecov options: minident=$MIN_IDENT, mindupid=$MIN_DUPIDENT, maxdup=$pMAXDUP,  idclass=$crclassf\n"; # maxerr=$MINE,
#NOT this way: $topline =~ s/idclass/hetz.ALN_SUB_NMI=off, idclass/ unless($ALN_SUB_NMI);
#this way: $topline =~ s/idclass/UCG.Merr_filter=on, idclass/ unless($ALN_SUB_NMI);
print $outh $topline; warn $topline if($debug and $outtab);

# assume inbam has all genes, want only UCG/busco subset, need idclassh for that
# * idclassh only used for isUCGene()
my($nidclass,$idclassh,$idclasslist)= ($crclassf) ? read_idclass($crclassf,1) : (0); # cds,te class by id, $idclass->{id} = class

use constant UCGENE_RDFILT => 0; #UPD21Jun04:off, 1 == skip read counts for not-UCGene(id)
my $testUCG= ($nidclass>0 and $UCGenePatt)?1:0;

sub isUCGene { # expand w/ idclass patt
  my($id)= @_;
  if($testUCG) { # allow also for all-are-ucg option
    my $idclass= $idclassh->{$id} || "";
    return($idclass =~ m/$UCGenePatt/i)?1:0; # 
  } else {
    return 1; # no choice
  }  
}

foreach my $itbam ($inbam,@inbam) {
  if($itbam ne $inbam) { close($inh);  open($inh,"samtools view -h $itbam |") or next; }
  warn "#  samtools view $itbam\n" if($debug);

  my($topid,$topidaln, $lrdlen, @savefirst);
  my $fixidPAIRMAP=0;  # can auto-detect PAIRMAP from srid and flags: srid fwd == srid rev
  
  while(<$inh>) {
    if(/^\@/) { 
      if(m/SN:/){ 
        my($td)=m/SN:(\S+)/; my($tw)=m/LN:(\d+)/; 
        #FIXME here? busco gene filter? maybe record all SN/LN, separate UCG later
        if($tw and $td){ $ngene++; $stw += $tw; $tw{$td}=$tw; } # FOXME: ngene++ wrong for @inbam>1, use scalar(%tw)
        #o: if($tw and isUCGene($td)){..}  
      }
      next; } 
    elsif(/^\W/){ next; } 
  
    my @v=split; 
    my($srid,$fl,$td,$tb,$qs,$cig,$rseq)= @v[0,1,2,3,4,5,9];

    #UPD22MAR: use SAMFLAG_read2 NOT SAMFLAG_rev to fix srid from pairmap, same id for fwd/rev reads needs changing *** other gnodes 2
    if($fl & SAMFLAG_read2){ # dont need isPAIRMAP opt
      if($fixidPAIRMAP > 0) { $srid .= "/2"; }
      elsif($fixidPAIRMAP==0) { $fixidPAIRMAP=($srid =~ m,/2,)?-1:1; $srid .= "/2" if($fixidPAIRMAP > 0); }
    }
    # if($isPAIRMAP and ($fl & SAMFLAG_read2)) {
    #   if($fixidPAIRMAP==0) { $fixidPAIRMAP=($srid =~ m,/2,)?-1:1; }
    #   if($fixidPAIRMAP >0) { $srid.="/2"; }
    # }
    
    # UPDmar18: LNtot fix, variable rdlen: add sum_rdlen_notok and n_notok, need sum_rdlen_total= srdlen_ok+srdlen_no & n_ok+n_no for proper total LN
    # srdlen_no,n_no need to include nomap, and notUCGene for flag < 2ndary
    # add opt to skip all rdlen < $MINRDLEN =~ 50..100, ie tiny frag reads in some SRR sets
    # now: $srdlen_tot= $srdlen+$srdlen_nomap; $n_rdtot= $n_mapok + $n_nomap + $n_mapnotucg + $n_mapbad; #  n_rdtot should == n_readid **
    
    my $rdlen=length($rseq);  $rdlen=$lrdlen if($rdlen<3); # empty seq for some, use last rdlen?
    if($rdlen<$MINRDLEN){
      if($fl < SAMFLAG_2ndary){ $n_rdtooshort++; next; }  #BUT check flags < SAMFLAG_2ndary  
      else { } # allow?? is this 2nd part match?
    }
    my $nextid= ($srid ne $lrid);
    if($nextid){ $n_readid++; $topid= $topidaln= 0;}
    $lrid=$srid; # assumes orig map/readid order
    
    #UPD: record all readid mapt to CDS genes, then filt by isUCGene
    #UPD: remove isUCGene() filt in read counts? do only in Cucg calc, table cov for all genes
    if($fl & SAMFLAG_nomap) { $n_nomap++; $srdlen_nomap+=$rdlen; next; }
    elsif( UCGENE_RDFILT and not isUCGene($td) ) { 
      if($nextid and $fl < SAMFLAG_2ndary) { $n_mapnotucg++; $srdlen_nomap+=$rdlen;  } #? $snotucg_rdlen += $rdlen; 
      next;
    }
   
    my($alen)= alignedbases($cig);
    my($nmi)= (m/NM:i:(\d+)/)?$1:0; # $ALN_SUB_NMI and 
    ## ALN_SUB_NMI NOT here: nmi only used for $locdup{$td}{'Uerr'} += $nmi
    my $idaln= $alen - $nmi;
    ##set MINALN for variable rdlen .. dont use fixed MINALN
    my $minaln= $MIN_IDENT * $rdlen; #  not MINALN

    # id pat fix for prefix:id >> SA:Z:([\w\:\-]+), ?? or use td chars
    my($zd,$zb,$zig)= (0,0,0);
    if(m/SA:Z:/){ ($zd,$zb,$zig)=(m/SA:Z:([\w\:\.\-]+),(\d+),.,(\w+)/); }

    use constant DO_FIRSTDUP => 0; # not sure want this addFirst(UorD) change
    # mis-calc: got *lots* of D firstdup when have zero S.Dup for most UCG

    sub addFirst {
      my($UorD, $nextid, $fl,$td,$tb,$cig, $zd,$zb,$zig, $rdlen, $rseq)= @_;

      my($alen,$lenc,$cend,$nlen)= addCigarz($UorD,$fl,$td,$tb,$cig, $zd,$zb,$zig, (($DOHETZ)?$rseq:""));
      $nmapid++ if($nextid); # should be same as $n_mapok, 2nd above, be sure
      $n_mapok++; $srdlen += $rdlen; # if($nextid) ??
      $nrd{$td}++; $nrdlen{$td} += $rdlen;  #<< valuable stats
    }
        
    if($fl >= SAMFLAG_2ndary) { # fixme: ($fl & SAMFLAG_2ndary); 
      #UPD21DEC: problem case ($td eq $topid) .. is this dupl self-align, or 2nd part align? skip or count?
      # if($nextid) { } # problem case: skip 1st align, have 2nd? do what? skip? treat as 1st?
      if($idaln < $minaln or $idaln < $MIN_DUPIDENT*$topidaln) {
        $n_dupbad++;
      } else {
        if($fl & SAMFLAG_2ndary) { # also have SAMFLAG_suppl, other?
        if(@savefirst){ addFirst('D',@savefirst); @savefirst=(); }
        my($alen,$lenc,$cend,$nlen)= addCigarz('D',$fl,$td,$tb,$cig, $zd,$zb,$zig); # skip err 2nd
        if(UPD21DEC and $topid and $td ne $topid) { $rdshare{$topid}{$td}{nrd}++; $rdshare{$topid}{$td}{aln}+= $alen; } 
          # ?? need backlink: $rdshare{$td}{$topid}{xxx} for full parlog counts
          # if($topid eq $td) ?? skip or keep
        }
      }
          
    } else { 
      if(@savefirst){ addFirst('U',@savefirst); @savefirst=(); }
      $topidaln= $idaln; $topid= $td;
      if($idaln < $minaln) {
        $n_mapbad++; $srdlen_nomap+=$rdlen;
        #? if($nmi > 0) { $n_mismatch += $nmi;  $locdup{$td}{'Uerr'} += $nmi; } # count nmi err        
        
      } else {
        # fixme for 1st of dup aligns? assign D after SAMFLAG_2ndary check
        # ** THIS change addFirst(D for U) throws off pMAXDUP calc; not sure wanted
        if(DO_FIRSTDUP) {
        @savefirst=($nextid, $fl,$td,$tb,$cig, $zd,$zb,$zig, $rdlen, (($DOHETZ)?$rseq:""));
        } else { # not DO_FIRSTDUP
        my($alen,$lenc,$cend,$nlen)= addCigarz("U",$fl,$td,$tb,$cig, $zd,$zb,$zig, (($DOHETZ)?$rseq:""));  
        $nmapid++ if($nextid); # should be same as $n_mapok, 2nd above, be sure
        $n_mapok++; $srdlen += $rdlen; # if($nextid) ??
        $nrd{$td}++; $nrdlen{$td} += $rdlen;  #<< valuable stats
        if(UPD21DEC) { $rdshare{$td}{self}{nrd}++; $rdshare{$td}{self}{aln}+= $alen; } #? need self counts
        
        }
        if($nmi>0){ $n_mismatch += $nmi;  $locdup{$td}{'Uerr'} += $nmi; }  #  $alen -= $nmi;
      }
    }

    $lrdlen= $rdlen;
  } # itbam
  if(@savefirst){ addFirst('U',@savefirst); @savefirst=(); }

} # for itbam (inbam,@inbam)

my ($topinfo,$ucmwMedn,$ucmwAve,$ucmwSE)=("",0,0,0);
use constant { kUCGfilter => 1, kIsUniqfilter => 2, k1000filter => 4, k2000filter => 8,
              kFullAlignFilter => 16, kUerrFilter => 32, };

if(UPD21JUN) {
  #>> only if have idclass w/ busco/UCG annots == my $testUCG= ($nidclass>0 and $UCGenePatt)?1:0;

  my $cucg_filt= kIsUniqfilter ; # ?opt for k2000filter
  $cucg_filt += k1000filter if($opt1000filter); #def on
  $cucg_filt += kFullAlignFilter if($optFullAlignFilter); #def on
  $cucg_filt += kUerrFilter if($optUerrFilter);
  my($btopinfo,$bucmwMedn,$bucmwAve,$bucmwSE)= topsummary_UCG("Measured Unique Genes",$cucg_filt); #  >= 1k long
  if($testUCG) {
    $cucg_filt=  kUCGfilter;
    $cucg_filt += k1000filter if($opt1000filter); #def on
    $cucg_filt += kFullAlignFilter if($optFullAlignFilter);
    $cucg_filt += kUerrFilter if($optUerrFilter);
    my($atopinfo,$aucmwMedn,$aucmwAve,$aucmwSE)= topsummary_UCG("UCG Class Genes", $cucg_filt); #  >= 1k long: now label all Filters:
    my $adiff= abs( $aucmwMedn - $aucmwAve);
    my $bdiff= abs( $bucmwMedn - $bucmwAve);
    if($adiff <= $bdiff) {
      $topinfo= $atopinfo . $btopinfo;
      ($ucmwMedn,$ucmwAve,$ucmwSE)= ($aucmwMedn,$aucmwAve,$aucmwSE); 
    } else {
      $topinfo= $btopinfo . $atopinfo;
      ($ucmwMedn,$ucmwAve,$ucmwSE)= ($bucmwMedn,$bucmwAve,$bucmwSE); #<< FIXME need choice from stats
    }
  } else {
    ($topinfo,$ucmwMedn,$ucmwAve,$ucmwSE)= ($btopinfo,$bucmwMedn,$bucmwAve,$bucmwSE);
  }
  print $outh $topinfo; warn $topinfo if($debug and $outtab);
  
} else {
  ($topinfo,$ucmwMedn,$ucmwAve,$ucmwSE)= topsummary_UCG();
  print $outh $topinfo; warn $topinfo if($debug and $outtab);
}

genestable($outh, $ucmwMedn,$ucmwAve,$ucmwSE);

if(UPD21DEC) { 
  my $gfouth= *STDOUT; 
  if($outtab) {
    my $gfout= $outtab; $gfout =~ s/\.\w+$//; $gfout.=".genefamtab"; # replaces .genexcopy ?
    rename($gfout,"$gfout.old") if(-s $gfout); 
    open($outh, '>', $gfout) or die "writing $gfout"; 
  }
  print $gfouth "# genefam_table \n";  
  genefam_table($gfouth, $ucmwMedn,$ucmwAve,$ucmwSE);
}

my $botinfo= join", ", grep /\w/, map{ my $bn=`basename $_ .bam`; chomp($bn); ($bn)?$bn:""; } ($inbam,@inbam);
print $outh "# inbam: ",$botinfo,"\n" if($botinfo);  


#--------------------------------------
use constant C_MEDIAN => 1; # use median not ave, tested both, ave influenced by extremes
use constant USE_CNZMEDIAN => 1; # test, then switch, looks better than cmw .. need for UCG C_MEDIAN also
use constant { sZERO => 9, cZERO => 2,
  spanZERO => 50, # percent of gene span covered, maybe higher: 66%? 75% most appear >90%, poorqual < 75%
  XCOPYisDUP => 1.66,
  cSKEW => 1.85, # C.nz.ave >> C.nz.median
  };

sub _max{ return($_[0] < $_[1])?$_[1]:$_[0]; }
sub _min{ return($_[0] > $_[1])?$_[1]:$_[0]; }
sub _minNot0{ return($_[0] == 0 or $_[0] > $_[1]) ? $_[1] : $_[0]; }

=item genefam_table

try1: dmag7fincds_SRR15012074_b2_upd.genefamtab
.. add CM col to match ,26.9 shared C val
# genefam_table 
Gene_ID	Glen	Nread	tCopy	nPar	Paralogs
Dapma7bEVm000001t1	6030	2035	1.0	1	Dapma7bEVm030496t1,922,26.9
Dapma7bEVm000002t1	6816	1713	1.1	3	Dapma7bEVm027372t1,1071,29.6	Dapma7bEVm017375t1,2,0.3	Dapma7bEVm000002t1,1,0.0
Dapma7bEVm000003t1	5187	1589	1.0	2	Dapma7bEVm027373t1,547,16.0	Dapma7bEVm000003t1,1,0.0
Dapma7bEVm000004t1	4974	1226	1.0	1	Dapma7bEVm027374t1,474,19.5
Dapma7bEVm000005t1	4257	5044	8.8	3	Dapma7bEVm027376t1,4313,141.1	Dapma7bEVm027375t1,2090,91.7	Dapma7bEVm000005t1,12,0.1
    ^^ example real paralog fam: 3 gene reps of ~9 copies total, Nuclear pore complex protein Nup155
Dapma7bEVm000006t1	7260	2519	1.1	1	Dapma7bEVm023601t1,310,19.6
Dapma7bEVm000008t1	6459	1449	0.9	1	Dapma7bEVm027377t1,914,21.0
Dapma7bEVm000009t1	12807	4140	1.1	8	Dapma7bEVm024240t1,1716,27.2	Dapma7bEVm012593t1,29,2.9	Dapma7bEVm014123t1,6,1.3	Dapma7bEVm026420t1,2,0.4	Dapma7bEVm014070t1,1,0.2	Dapma7bEVm000009t1,1,0.0	Dapma7bEVm021698t1,1,0.1	Dapma7bEVm007362t1,1,0.1
   ^^ fixed self and trivial shares: eg: Dapma7bEVm000009t1,1,0.0	Dapma7bEVm021698t1,1,0.1	
=cut

sub genefam_table { # from %rdshare
  my($outh,$ucmwMedn,$ucmwAve,$ucmwSE)= @_;

#  if(UPD21DEC) { $rdshare{$topid}{$td}{nrd}++; $rdshare{$topid}{$td}{aln}+= $alen; } # if($topid eq $td) ?? skip or keep
#  if(UPD21DEC) { $genetab{$td}=[$tw,$rc,$copyNumEst, $CMint, $uok,"$cnz,$cnzmed,$nnz,$cspan"]; } 
  my $Cmap= (C_MEDIAN) ? $ucmwMedn : $ucmwAve; # was: $mcmw : $ucmt; # genestable
  my $pMINRDSHARE= $ENV{minrdshare}|| 0.20; # what? option?
  my $SKIP_NOPAR= $ENV{nonopar} ? 1 : 0; #?? dont need opt, just include all npar == 0
  
  my($nout,$didgfhdr,%topids)=(0,0);
  my @td= sort keys %genetab;
  for my $td (@td) { 
    unless(exists $rdshare{$td}) { 
      # ugh need to scan thru all topid in %rdshare for these? make backlink hash?
      for my $top (@td) {
        if(exists $rdshare{$top} and exists $rdshare{$top}{$td}) {
          $topids{$td}{$top}++;
        }
      }
    }
  }
  
  # all have rdshare{td}{self}?
  for my $td (@td) { 
    my $gtab= $genetab{$td};
    my ($tw,$rc,$copyNumEst, $CMint, $uok, $Cnz)= @$gtab;
    my $minrdshare= $rc * $pMINRDSHARE;
    my ($is2nd,@aid,@parlogs)=(0);
    unless(exists $rdshare{$td}) {  
      my @topa= sort{ $topids{$td}{$b}<=>$topids{$td}{$a} } keys %{$topids{$td}};
      @aid= @topa; # $is2nd=1;
      #^^ this way results in fewer @parlogs, otherwise need to save more in rdshare hash?
      #upd: add info tag for is2nd: count all topa[0] kids?
      if( @topa ) { $is2nd= scalar(  keys %{$rdshare{$topa[0]}} ); } # is this bad?
      $is2nd ||= 1; 
    } else {
      @aid= sort{ $rdshare{$td}{$b}{nrd}<=> $rdshare{$td}{$a}{nrd}}  keys %{$rdshare{$td}};
    }
    
    # my($aself)= grep{ $_ eq 'self' } @aid;#? want this
    #below: @aid= grep{ $_ ne 'self' and $_ ne $td } @aid;
    
    my $selfp=0;
    for my $ad (@aid) {
      #o: next if( $ad eq 'self' or $ad eq $td );
      next if($ad eq $td );
      #upd: maybe add self at end, if have othersl add self to np count, in clean cases np+1 == tCopy
      my $nrdshare=$rdshare{$td}{$ad}{nrd};
      my $alnshare=$rdshare{$td}{$ad}{aln};
      next if($nrdshare < 2 or $nrdshare < $minrdshare); # filter trivial shares
      my $atw= ($ad eq 'self')? $tw{$td} : $tw{$ad}||0; #?$genetab{$ad}->[0];
      #  my($atw,$arc,$acopyNumEst, $aCMint, $auok, $aCnz)= @{$genetab{$ad}};
      my $mtw= _minNot0($atw,$tw);
      my $cmwshare= ($mtw<1) ? 0: $alnshare/$mtw; 
      #NOT HERE?? if(USE_CNZMEDIAN) .....
      my $CMint=  sprintf"%.1f", $cmwshare; # not int($cmw);
      #?? no value? my $copyNumEst= sprintf"%.1f", $cmwshare/$Cmap; #est from cmw or cnz = $sln/$nnz ?
  
      my $parlog="$ad,$nrdshare,$CMint";
      if($ad eq 'self'){ $selfp= $parlog; }
      else { push @parlogs, $parlog; }
    }
    
    my $np= @parlogs; # if $np > NBIG, compress @parlogs ?? drop out low align set?
      #  BUT want to see valid cases of 100 paralogs
    next if($SKIP_NOPAR and $np==0); # skip uniqs..
    unless($didgfhdr++) {
      my @hcols= qw(Gene_ID Glen Nread tCopy CM nPar Paralogs);
      # my @hcols= qw(Gene_ID Glen Nread Rdlen tCopy C.M Uniq Merr C.nz S.Dup,Unq,Nt,Equ,UCG);
      print $outh join("\t",@hcols)."\n";
    } 
    my @cols=($td,$tw,$rc, $copyNumEst, $CMint); #(, $uok,$cnz);  # dont need all genetab cols here
    if($is2nd) { push @parlogs,"minor-of-$is2nd-paralogs"; }
    push @parlogs, $selfp if($selfp);
    print $outh join("\t",@cols,1+$np,@parlogs)."\n"; $nout++; # 1+np here adds self
  }
  return($nout);
}  

sub genestable {
  my($outh,$ucmwMedn,$ucmwAve,$ucmwSE)= @_;
  my $didhdr=0;
  my $Cmap= (C_MEDIAN) ? $ucmwMedn : $ucmwAve; # was: $mcmw : $ucmt;
  
  # main gene data globals: %tw, %nrd %cmd  %loc{$td}[0..tlen]{DUN}
  for my $td (sort keys %tw) { 
    #UPD21DEC: Note cmd{td} sums all align-M spans on gene, but nrd,rclen count reads only for 1st map
    # .. can have dup genes where 2nd near identical has many fewer nrd,rclen than aligns (cmd) indicate
    # .. should nrd,nrdlen change to count 2ndary maps? not yet, needs tests
    my $tw=$tw{$td}; my $cmd=$cmd{$td}||0; 
    my $rc=$nrd{$td}||0; my $rclen=$nrdlen{$td}||0;
    my $cmw= $cmd/$tw; my $crl= $rclen/$tw; # drop int(); was int($rc*$RDLEN/$tw); 
  
    my @lns=(); # median of cov N val for gene span bases, maybe should replace $cmd{td}/tw, cmp to cnz
    my($sld,$slu,$sln,$slerr,$slun,$neq,$nnz, $shet, $shetho)=(0) x 19;
    # ($sld,$slu,$sln,$slerr)= map{ $locdup{$td}{$_} || 0 } qw( D U N Uerr);
    $slerr= $locdup{$td}{'Uerr'}||0;
    #>> change this per gene loc{td}[i] to locdup{td}{(D U N)} ?? need neq, nnz ?
    for(my $i=0; $i<$tw; $i++) { 
      my($ld,$lu,$ln)= map{ $loc{$td}[$i]{$_}||0 } qw(D U N); 
      if($DOHETZ){
        my($ba,$bc,$bg,$bt)= map{ $loc{$td}[$i]{$_}||0 } qw(A C G T); 
        my $bn=$ba+$bc+$bg+$bt; # should == $ln
        my ($b1h,$b2h)=sort{ $b <=> $a} ($ba,$bc,$bg,$bt);
        # UPDhetz: add cds 1,2,3 position stats ($i % 3)
        if($bn > sZERO) {
         my $ishet= ($b1h < $bn*$HOMIN and $b2h >= $bn*$HOMHET)?1:0; # does this hetz stat require b2h be high, > 10%?
         $shet++ if($ishet); $shetho++;
        }
      }
      $sld+=$ld; $slu+=$lu; 
      if($ln>0){ 
        $sln+=$ln; $nnz++; push @lns,$ln; 
        if($lu == $ln){ $slun+=$lu; $neq++; }
      }
    }

    #NOTE: cnz, cnzmed are now best-reliable cov-depth/copynum stats
    my $ceq=($neq<sZERO)?0:sprintf"%.1f", $slun/$neq;
    my $cnz=($nnz<sZERO)?0:sprintf"%.1f", $sln/$nnz;
    my $cnzmed=0; if($nnz>=sZERO){ @lns=sort{$b<=>$a}@lns; $cnzmed= sprintf"%.0f", $lns[int($nnz/2)]; }
    # add spancov = $nnz/$tw == fully covered span .. or calc depthcov over span? for i=0..tw, scspan+= ln - cmw ?
    my $cspan= ($tw<sZERO)?0:sprintf"%.1f", 100*$nnz/$tw;
    
    #o: my $uok= ($sln<sZERO)?"zero":($sld/$sln > $pMAXDUP)?"dupl":"uniq";
    my $merr=($slu<sZERO)?0:sprintf"%.1f",  100*$slerr/$slu;
    my $phet=($shetho<sZERO)?0:sprintf"%.0f,%d/%d",  100*$shet/$shetho,$shet,$shetho;

    my($CMint,$copyNumEst)=(0,0);
    if(USE_CNZMEDIAN) {
      $CMint= $cnzmed; $copyNumEst= sprintf"%.1f", $cnzmed/$Cmap; #est from cmw or cnz = $sln/$nnz ?
    } else {
      $CMint=  sprintf"%.0f", $cmw; # not int($cmw);
      $copyNumEst= sprintf"%.1f", $cmw/$Cmap; #est from cmw or cnz = $sln/$nnz ?
    }
    # my $uok= ($sln<sZERO)?"zero":($sld/$sln > $pMAXDUP or $copyNumEst >= 2)?"dupl":"uniq";

    # UPD21jun22: add 'skew' class, sort of replacing 'dups', skew = C.nz >> C.M (ave >> median)
    # .. this is better indicator of problem genes with subseq of high copy, other of low copy, too messy to call 'gene copy number'
    my $cskewd= ($cnz > cSKEW * $cnzmed)?1:0; # what cSKEW cut off?
    
    # ^^ modify uok class: sld/sln>maxdup == dupls, xcopy >= 1.66 == duplx
    my $uok= ($sln<sZERO or $cspan < spanZERO)?"zero"
      :($cskewd)?"skew"
      :($copyNumEst >= XCOPYisDUP)?"dupx"
      :($sld/$sln > $pMAXDUP)?"dups" #? keep or not
      :"uniq"; 

    my $ucgf="";  # add isUCGene() flag ?? which col? S.D,U,N,E?
    if($testUCG) { $ucgf= (isUCGene($td))?",UCG":""; }
    
    # fixme better output: keep: td,w,rc,crl,cmw, add:uok, maybe: cnz, ceq, maybnot: sld,u,n,un
    # Gene_ID Glen C.LN C.M Uniq C.nz C.eq Nread Rdlen S:dup,unq,nt,equ
    # ** DROP C.LN as useless/noise, replace col?
    # ** replace C.LN = int($crl) w/ xCopy = copyNumEst
    # ** C.nz is good *ave* est, C.M is best median est 
    # ** FIXME: C.M from $cmd{$td} is *not* median score/gene, but ave. Uniq cov depth/gene span
    # .. test cnzmed as replacement for this,
    # Change Uniq col calc: $uok not useful when geneset is all uniq seq, but rdmap cov calls dup count

    unless($didhdr++) {
      # tCopy (total) label corrected was xCopy (excess)
      my @hcols= qw(Gene_ID Glen Nread Rdlen tCopy C.M Uniq Merr C.nz S.Dup,Unq,Nt,Equ,UCG);
      push @hcols,"Hetz" if($DOHETZ);
      print $outh join("\t",@hcols)."\n";
    } 
    my @cols=($td,$tw,$rc,$rclen, $copyNumEst, $CMint, $uok, $merr, "$cnz,$cnzmed,$nnz,$cspan", "$sld,$slu,$sln,$slun$ucgf");  
    push @cols,$phet if($DOHETZ);
    print $outh join("\t",@cols)."\n"; # off: ,"$ceq,$neq"

    if(UPD21DEC) { $genetab{$td}=[$tw,$rc,$copyNumEst, $CMint, $uok,"$cnz,$cnzmed,$nnz,$cspan"]; } 
      # save all? for rdshare table output 
  }  
}


sub topsummary_UCG {
  my($sumlabel, $whichUniqGeneFilter)=@_;
  $whichUniqGeneFilter= kUCGfilter unless($whichUniqGeneFilter);
  $sumlabel ||= "UCG Class Genes";
  my $skipALLsum= (UPD21JUN)?1:0; 
  
  ## fixmed: add median rct from $nrdlen{$td}/$tw{td}
  ## FIXMEd: use only "uniq" set for summary stats .. need to calc after read loop
  ##   my $uok= ($sln<9)?"zero":($sld/$sln > $pMAXDUP)?"dupl":"uniq";

  # UPD21jun22: add 'skew' class, sort of replacing 'dups', skew = C.nz >> C.M (ave >> median)
  # .. this is better indicator of problem genes with subseq of high copy, other of low copy, too messy to call 'gene copy number'
  
  # UPD UCGENE_RDFILT removed isUCGene() filt in read counts? 
  #     add *here* in Cucg calc, table cov for all genes
  # fix @inbam $ngene, and $stw= sum(tw)
  my $ngene= scalar(keys %tw); # tw now includes not-isUCGene
  my $stw=0; for my $td (sort keys %tw){ $stw += $tw{$td}; }

  ##UPD21NOV: retest nonmi=1 == ! ALN_SUB_NMI here: nmi only used for $locdup{$td}{'Uerr'} += $nmi
  ## this retest no good, reduces C.ucg, likely favors *partial* cds aligns.
  ## $topline =~ s/idclass/UCG.Merr_filter=on, idclass/ unless($ALN_SUB_NMI);
  #off: my $ALN_SUB_NMI = (defined $ENV{nonmi}) ? 0 : 1; # UPD21NOV: TEST for hetero adjust, dont subtract NM:i:nn  from align val

  ## re-retest new filter for nonmi=1 : drop *partial* aligns, from cspan of C.nz tuple:
  ## above: "$cnz,$cnzmed,$nnz,$cspan",  my $cspan= ($tw<sZERO)?0:sprintf"%.1f", 100*$nnz/$tw;

  #now opt above: $whichUniqGeneFilter |= kFullAlignFilter if(defined $ENV{nonmi}); # UPD21NOV DEBUG TEST
  my $cspan_min = 90; #? option? dont want too low,
  
  use constant pUERR_MAX =>  0.33; #? test opts? 0.25? 0.50?
  my $uerr_max= 999999;
  if($whichUniqGeneFilter  & kUerrFilter) { # was unless($ALN_SUB_NMI)  # drop this, no good (?)
    my @uerr= sort{ $b <=> $a } map{ $locdup{$_}{'Uerr'}||0 } keys %locdup; #? limit ids to UCG here?
    my $nerr= @uerr;
    $uerr_max= $uerr[ int( pUERR_MAX * $nerr) ]; # drop worst 1/3, 1/4 ?
  }
   
  # totals including dup genes
  my $rct=($stw<sZERO)?0:sprintf "%.1f", $srdlen / $stw; #was $n_mapok * $RDLEN / $stw ;
  my $cmt=($stw<sZERO)?0:sprintf "%.1f", $scm / $stw;  
  my $ardlen= ($n_mapok<sZERO)?0:sprintf "%.0f", $srdlen / $n_mapok; 
  
  #UPD18mar new LN totals: readlen*n_readid
  my $srdlen_tot= $srdlen+$srdlen_nomap; 
  my $nrdlen_tot= $n_mapok + $n_nomap + $n_mapnotucg + $n_mapbad; #  nrdlen_tot should == n_readid **

  my $labelFILT="Filters:";
  $labelFILT .= "UCG," if($whichUniqGeneFilter & kUCGfilter);
  $labelFILT .= "size>=2k," if($whichUniqGeneFilter & k2000filter);
  $labelFILT .= "size>=1k," if($whichUniqGeneFilter & k1000filter);
  $labelFILT .= "covspan>=$cspan_min," if($whichUniqGeneFilter & kFullAlignFilter);
  $labelFILT .= "mErr<$uerr_max," if($whichUniqGeneFilter & kUerrFilter);
  
  my($urdlen,$utw,$ucm,$ucmw,$ssucmw)=(0) x 9; 
  my(@crl,@cmw);
  for my $td (sort keys %tw) { 
    my $tw= $tw{$td} or next; 
    my $gok=1;
    
    # (UPD21JUN) 
    if($whichUniqGeneFilter & kUCGfilter) { $gok=0 unless( isUCGene($td)); }
    elsif($whichUniqGeneFilter & kIsUniqfilter) { } # == ($sln < sZERO or $sld/$sln > $pMAXDUP) below
    if($whichUniqGeneFilter & k2000filter) { $gok= 0 if($tw < 2000); }
    elsif($whichUniqGeneFilter & k1000filter) { $gok= 0 if($tw < 1000); }
    else { $gok=0 if($tw < 200); } #? < 200 lower limit, eg arabid has a few 32 bp cds, other ref prots < 100aa
    
    my $cspan= -1; # always calc, add to sums? No use here, is in genestable()
    if($gok and $whichUniqGeneFilter &  kFullAlignFilter) {
      # ?? is  $locdup{$td}{'N'} same sum N? No, that counts dups, need ln > 0 per base
      my $nnz=0;
      for(my $i=0; $i<$tw; $i++) { 
        my $ln= $loc{$td}[$i]{'N'}||0; 
        $nnz++ if($ln>0);
      }
      $cspan= 100*$nnz/$tw;
      $gok=0 if($cspan < $cspan_min); # kFullAlignFilter
    }
    
    ##UPD21NOV: ALN_SUB_NMI here: nmi only used for $locdup{$td}{'Uerr'} += $nmi
    if($gok and $whichUniqGeneFilter &  kUerrFilter) {
      my $uerr= $locdup{$td}{'Uerr'}||0; # == slerr below
      $gok=0 if($uerr > $uerr_max); # kUerrFilter
    }
    
    next unless($gok);
    
    my $sln= $locdup{$td}{'N'}||0;  
    my $sld= $locdup{$td}{'D'}||0;
    # check td span cover for < spanZERO as bad data, notuniq
    next if($sln < sZERO or $sld/$sln > $pMAXDUP); # $isuniq=0 
    # next unless($isuniq);
    
    my $cmd= $cmd{$td}||0; 
    my $rclen=$nrdlen{$td}||0; # any zero/missing?  
    my $crl= $rclen/$tw;  # dont int() these, lost precision..
    my $cmw= $cmd/$tw; # primary C.m val == read aligns sum/gene from addCigar; 
                       # /tw is problem for partialAlign genes, thus kFullAlignFilter;
                       # alternate fix: use $nnz, ie covered span, but that affects meaning of C.m, can use C.nz tuple vals to correct
    $urdlen += $rclen; $ucm += $cmd; $utw += $tw; 
    $ucmw += $cmw; $ssucmw += $cmw*$cmw; # use ucmw for ave, not ucm/utw
    push @crl,$crl; push @cmw,$cmw;  
  }
  
  @crl= sort{ $b <=> $a } @crl; @cmw= sort{ $b <=> $a } @cmw;
  my $nc= @cmw; my $nch= int($nc/2);
  my $mclr= $crl[$nch]; my $mcmw= $cmw[$nch];
  
  my($ucmwMedn,$ucmwAve,$ucmwSD,$ucmwSE)=($cmw[$nch],0,0,0);
  if($nc > sZERO) {
    $ucmwAve= $ucmw/$nc; 
    $ucmwSD= sqrt(($ssucmw - $ucmwAve*$ucmwAve)/($nc-1)); 
    $ucmwSE= $ucmwSD/sqrt($nc); #? change to report ucmwSD , nc ?
  }
  
  # problem of extremes in  ucm/utw, use instead meanof(@cmw) as better stat
  my $urct=($utw<sZERO)?0:sprintf "%.1f", $urdlen / $utw; #was $n_mapok * $RDLEN / $stw ;
  my $ucmt= $ucmwAve; # old: ($utw<sZERO)?0:sprintf "%.1f", $ucm / $utw;  
  
  my $nmapcdsid= $nmapid+$n_mapnotucg;
  my $pmapucg= ($n_readid<sZERO)?0:sprintf"%.1f", 100*$nmapid/$n_readid; # pmr FIXME should reduce nmapid by filters, ?? here?
  my $pmapcds= ($n_readid<sZERO)?0:sprintf"%.1f", 100*$nmapcdsid/$n_readid;
  
  # UPD add map err stats, maybe mod Gsize est
  # $n_readid, $n_nomap, $n_mapok,  $n_mapnotucg,
  #  base errs:  $n_insert, $n_delete,  $n_mismatch
  #  relative to either $srdlen (nucg reads * rdlen) or $scm = mapt bases
  my $nberr = $n_mismatch + $n_insert + $n_delete; # mid or sid; see also new locdup{td}{Uerr}
  my $pberr = ($scm<sZERO)?0:sprintf "%.2f", 100*$nberr / $scm;  
  
  # FIXME : med: mcmw or ave: ucmt ? see $ucmwAve update, best stat? 
  # with clean data, med =~ ave, w/ unclean data median is more reliable est of true value
  #above: use constant C_MEDIAN => 1; # use median not ave, tested both, ave influenced by extremes
  
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
  my $srdlen_totEST= ($n_mapok>=sZERO) ? $srdlen * $n_readid / $n_mapok : 0; 
  #old: my $dogsalt= (($n_mapok>=sZERO) and ($srdlen_tot == $srdlen * $n_readid / $n_mapok)) ? 0 : 1;  
  # my $MINLENDIF=1; #what ? ratio($srdlen_tot, $srdlen_totEST)>=0.999 ?? << use ratio
  # my $dogsalt= (($n_mapok>=sZERO) and abs($srdlen_tot - $srdlen_totEST) < $MINLENDIF) ? 0 : 1;  
  my $pDIFFLEN= 0.999;
  my $dogsalt= ( $n_mapok>=sZERO and ratio($srdlen_tot, $srdlen_totEST) >= $pDIFFLEN ) ? 0 : 1;
  if($dogsalt){
    my $Genosize_alt= ($Cmap<cZERO) ? 0 :  sprintf"%.1f", $srdlen_tot / $Cmap / 1_000_000; 
    my $ardlen_alt  = ($nrdlen_tot<1)? 0 : sprintf"%.1f", $srdlen_tot / $nrdlen_tot;
    my $GSEeqn_alt  =  sprintf "LN/C= %d / %.1f for Lr: $ardlen_alt, Nr: %d = %d+ %d+ %d+ %d [ok,nomap,noucg,bad]", 
      $srdlen_tot, $Cmap,
      $nrdlen_tot, $n_mapok, $n_nomap, $n_mapnotucg, $n_mapbad; # srdlen_tot sum of all readlens counted, *should* == ave(rdlen) * n_readid
    $addGSEalt=  "# Gsize_alt  = $Genosize_alt Mb of $GSEeqn_alt\n";
  }
  
  #UPD2202: change precision .1f to .2f for the crucial Cucg vals
  my $CmapAveSE= sprintf "%.2f, %.2f +/-%.2f (mdn,ave,sem)",$ucmwMedn,$ucmwAve,$ucmwSE; # replace C.Map/W=$mcmw, $ucmt
  $mclr=int($mclr);
  
  # UPD21JUN: drop C.LN/W=$rct (all), C.LN/W=$mclr (ucg) as uninformative
  #  add a. MIN Glen  filter (MIN_GLEN=1000 default?), b. two UCG est: b1. busco/idclass, b2. 'uniq' measured set
  #  maybe drop 'All Gene Cov Depth' as uninform.

  my $alldepth=($skipALLsum)?"":"# All  Gene Cov Depth n=$ngene, C.Map/W=$cmt ave, for W.genes=$stw, LN.reads=$srdlen\n"; # C.LN/W=$rct, 

  my $topinfo= "# $sumlabel, $nc of $ngene total, $labelFILT --------------\n" . $alldepth . 
  "# Uniq Gene Cov Depth n=$nc, C.Map/W=$CmapAveSE for W.genes=$utw, LN.reads=$urdlen\n" .   # C.LN/W=$mclr,  
  "# Genome_size= $Genosize_mb Mb of $GSEeqn, CDS_size= $CDS_mb Mb,\n" .
  $addGSEalt .
  "#   for Nr.ucgcds=$nmapid ($pmapucg%), Nr.allcds=$nmapcdsid ($pmapcds%), Nr.total=$n_readid, Lrdlen=$ardlen, MapErr=$nberr ($pberr%), Nshort=$n_rdtooshort\n"; 
  #UPD18mar: add new LN.reads=srdlen_tot, others above
  #my $srdlen_tot= $srdlen+$srdlen_nomap; 
  #my $nrdlen_tot= $n_mapok + $n_nomap + $n_mapnotucg + $n_mapbad; #  nrdlen_tot should == n_readid **

  return($topinfo,$ucmwMedn,$ucmwAve,$ucmwSE);
}


sub ratio {
  my($va,$vb)=@_; 
  my($r,$inv)=(0,0);
  if($va > $vb){ ($va,$vb,$inv)=($vb,$va,1); }
  if($vb < 0.01){ $r= 0; } else { $r= $va/$vb;  }
  return (wantarray)? ($r,$inv) : $r;
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

sub alignedbases{ my($cig)=@_; my $alen=0; while($cig =~ m/(\d+)M/g) { $alen+= $1; } return $alen; }
 
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


