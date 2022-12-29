#!/usr/bin/perl
# gnodes_genescov.pl
# from gnodes_covsum.pl for evigene/scripts/genoasm/ 

=item usage

  gnodes_genescov.pl -cdscov xxx_cds.covtab -chrcov xxx_chr.covtab -out xxx_genes.sumtab? \
    -idclass xxx.idclass -anntab xxx.anntab 

  dropse20gnodes/
  $evigene/scripts/genoasm/gnodes_genescov.pl -debug  \
    -title dropse20cdschr -asmid dropse20chrs -metad dropse20chrs.metad \
    -idclass dropse20cdste.idclass -ann dropse20chrs.fa.anntab \
    -cdscov dropse20t1cds_SRR11813283_1_bwa.cdschr7b.covtab -chrcov dropse20chrs_SRR11813283_1_bwa.cdschr7b.covtab

=item fixme do/dont count all genes at same locus ??

  UPD21apr26: readCovtab
   ?? revert to count 1 id/bin locus, more skews gene sums : No, not that either
   .. both ways are wrong: count all ids over-counts some lowqual id-hits, count 1 under-counts some paralogs
   .. likely need better cds.anntab, separate best/2nd cds id aligns, count all best/top per bin-locus, not 2nds

=item fixme  cds.len, stat for dupl/under-cov using cdslen - chr.align bins (in anntab) x chr.cov
  
  .. chrcov act does not measure full cds cov (?), 
  .. maybe want sum(acm) over chr bins vs expected CU (KU) x cdslen/BINSIZE ?
    -- if acm =~ CU and chrbin =~ cdslen/BN, have expected xCopy = 1,
    -- if acm < CU but chrbin > cdslen/BN, have spurious dupls, xCopy < 1, ndup ~ chrbin/cdsbin
    -- if acm =~ CU and chrbin > cdslen/BN, have real dups, xCopy = 1, ndup ~ chrbin/cdsbin
    -- if acm > CU, xCopy > 1, ..
    
=cut    

use strict;
use Getopt::Long;  

my $debug=$ENV{debug}||1;
my @IVAR=(5,6,7); # covtab cols ==  rd_total, rd/mmap, rd_uniq for chrasm x rd map 
my $DOMED= 0; # $ENV{median}||0; * SUCKS UP ALL covtab vals, mem pig
my $DOHIST=1; # replace median from @allvals w/ histo mode(s)/peak, est median, counts for depth vals depth[0..999]?
my $ZEROS= 0; # $ENV{zero}||0;
my $showSUM=1; # $ENV{sum}||0;
my $NTAVE= $ENV{ntave} || 0;
my $MINRD= 0; #? skip all cov rows w/ acovt < MINRD, want to measure zeros? gaps?
my $CBIN=100; # calc from rows ib dist
my $TOTALID="total"; #? "allgenes";
my $KU_OPT=0; # read from metad, or calc from genemeans, unless -CU|kucg=val
my $XHICUT= $ENV{xhicut}|| 1.7; my $XLOCUT= $ENV{xlocut}|| 0.55;# xCopy opts, change def? xhi=1.60 .. 1.66? xlo=0.50..0.54?
my $TRIMAVE= $ENV{trimave}||0; # for ave calc, trim extremes?

use constant kSAMPLE => 1000;  # CBIN size median sample @cbin, should be all same = default $CBIN

use constant UCG_KU_SPECIAL => 1; # special KU cov-depth calc for uniq-conserved-genes: median of median depth
my($nucg,@ucg_med_cov,@ucg_med_ids,@ucg_one_cov,%ucg_id_cov)=(0); # UCG_KU_SPECIAL globals

my $intitle= $ENV{title}||$ENV{name}||""; # was "asource"; 
my $asmid= $intitle; #? intitle for output, asmid for data
my ($sampledata,$anntable,$outmeans,$outsum,$ncpu,$crclassf,$readsetid)=("") x 9;
my ($cdscovtab,$chrcovtab); # input data
my (@lvar);


my $optok= GetOptions( 
  'cdscovtab=s',\$cdscovtab, 'chrcovtab=s',\$chrcovtab,
  'output|outsummary=s', \$outsum,
  'means=s', \$outmeans,
  'title|name=s',\$intitle,
  'asmid=s', \$asmid, # superceeds intitle for finding/using asm info
  'readsetid=s', \$readsetid, # upd21f20
  # 'filter=s',\$filtn, #<< replace w/ @parts loop over filters
  'genomedata|sumdata|sampledata|metadata=s', \$sampledata, # optional input  
  'anntable|annotations=s', \$anntable,
  'crclassf|idclassf=s',\$crclassf, # table of chr/gene id => class ? want for BUSCO, TEfam, : No, rely on anntable having these
  'minread=i', \$MINRD, # 5 default, test, see if it affects over-assembly (Xcopy < 1)
  # 'Rformat!', \$RTABLE, 
  # 'median!', \$DOMED, # -nomedian to save mem  >> change to depth_histogram, count hits of D[0..999],Dk[(1000..10000)/10] ? Dx[>10k]
  # 'histogram|peak|mode!', \$DOHIST, # -nomedian to save mem  >> change to depth_histogram, count hits of D[0..999],Dk[(1000..10000)/10] ? Dx[>10k]
  'CU|kucg=s', \$KU_OPT,
  'XHICUT=s', \$XHICUT, 'XLOCUT=s', \$XLOCUT, # xCopy hi/lo cutoff from 1.0
  'ncpu=i', \$ncpu, # not used here
  'debug!', \$debug, 
  );

die "usage: gnodes_genescov.pl -cdscovtab xxx_cds.covtab -chrcovtab xxx_chr.covtab -anntab xxx.anntab  -name xxxcdschr 
 opts: -metadata xxx.metad -asmid xxx_asm -idclass xxx.idclass -CU=0 -minread=$MINRD -debug 
" unless($optok and ($cdscovtab or $chrcovtab or $outmeans));

if($asmid and not $intitle) { $intitle=$asmid; } # want one of these
elsif($intitle and not $asmid){ $asmid=$intitle; }
unless($readsetid) {
  my $rdset=""; # FIXME leads to messy fname if no SRRnnnn in covtabs
  for my $incov ($cdscovtab,$chrcovtab,@ARGV) {
    if($incov and $incov =~ m/(\w+).(\w+)\.covtab/) { 
      my($asmrd,$vers)=($1,$2); 
      $asmrd =~ s/($asmid|$intitle)//g;
      if($asmrd =~ m/([A-Z]RR\d+)/) { my $na=$1;  $rdset .= $na unless($rdset=~m/$na/); last; }
      elsif($asmrd=~/([A-Za-z0-9]+)/){ my $na=$1;  $rdset .= $na unless($rdset=~m/$na/); last; }
    }
  }
  if($rdset =~ m/\w\w/){ $rdset=substr($rdset,0,19) if(length($rdset>19)); $readsetid=$rdset; }
}

my $title= $intitle;
unless($outsum) { 
  ($outsum=$intitle) =~ s/\W/_/g; # title always ok?
  $outsum.="_".$readsetid if($readsetid and $outsum !~ m/$readsetid/);
  $outsum.="_genesum.txt";
}
unless($outmeans){ ($outmeans=$outsum) =~ s/.genesum.*//; $outmeans.=".genemeans"; }
#add here or putGeneCovsum: outxcopy=outname.genexcopy

warn "#genecov output to sum=$outsum means=$outmeans\n" if $debug;

sub MAINstub {}

my($nsamp,$masmid,%samvals); # masmid == asmid nor not?
($nsamp,$masmid)= readMetad($sampledata);

my($nann,%annvals);
($nann)= ($anntable) ? readAnnots($anntable) : (0); # format from gnodes_annotate.pl, same crID, crLoc, an, ids as covtab

my($nidclass,$idclassh,$idclasslist)= ($crclassf) ? read_idclass($crclassf, 1) : (0); 

my($ntcds,$cdslen,$cdsrdmap)= ($cdscovtab) ? readChrtab($cdscovtab, 1) : 0; 

my($nccds)= ($cdscovtab) ? readCovtab($cdscovtab, 1) : 0;
my($ncchr)= ($chrcovtab) ? readCovtab($chrcovtab, 0) : 0;

# die "Missing input -cdscov $cdscovtab" unless($nccds > 0);
# die "Missing input -chrcov $chrcovtab" unless($ncchr > 0);
if($nccds > 0 or $ncchr > 0) {
  putGeneCovStats($outmeans); 
  %ucg_id_cov=(); # empty hash of all data, use outmeans
}

# step2 usage: read outmeans, write outsum, dont need covtabs
if(-s $outmeans) {
  putGeneCovSum($outsum,$outmeans); 
}

#=======================================================

# my($ntcds,$cdssizes)= ($cdscovtab) ? readChrtab($cdscovtab, 1) : 0; # 
# .cdschr7b.chrtab
# #ChrID	chrlen	nmap	nuniq	nmult	nmread	nomap
# dropse20uc:g117183162t1	489	722	214	508	343	0
# dropse20uc:g117183163t1	3459	2816	1884	932	2088	0

sub readChrtab { # read cds.chrtab if avail for sizes, also chr.chrtab?
  my($covtab, $iscdscov)=@_; #? for both cds,chr covtabs?
  my($nt,%len,%rdmap)=(0);
  (my $chrtab=$covtab) =~ s/covtab/chrtab/;
  open(my $inh,$covtab) or return 0;
  while(<$inh>) {
    next if(/^\W/);
    my @v=split; my($id,$clen,$nrmap,$nuniq,$nmult,$nmread)=@v;
    $len{$id}=$clen; $nt++;
    $rdmap{$id}=[$nrmap,$nuniq,$nmread]; # or all?
  } close($inh);  
  warn "#readChrtab($chrtab) nok=$nt\n" if $debug;
  return($nt,\%len,\%rdmap);
}

sub readCovtab { 
  my($covtab, $iscdscov)=@_; #? for both cds,chr covtabs?
  open(my $inh,$covtab) or die "reading $covtab";
  
  my($nin,$nok,$lcr,$lcb,$lida)= (0) x 9;
  my(@cbin);
  
  while(<$inh>) {
    my @v=split; 
    my $hasivar= (@v >= 8); my $hasann= (@v >= 9);
    if(/^\W/){ 
      if($hasivar and not m/[,=]/ and /Cov/i) { @lvar=@v[@IVAR] unless(@lvar); }  
      next;
    }
    next unless($hasivar); # error? 
    my($cr,$cb,@av)=@v; $nin++;
    
    my($cCDS,$cTE,$cUNK)= splice(@av,0,3); #($ct,$cm,$cu) NOW == CDSrd, TErd, UNKrd
    my($act,$acm,$acu)  = splice(@av,0,3); # == @IVAR, same  rd_total, rd/mmap, rd_uniq for chrasm x rd map 
    my($an,$ids)= (@av>0) ? splice(@av,0,2) : ("",""); # may not exist

    if($nann>0 and my $anidx= $annvals{$cr}{$cb}){
      my($anx,$idx)=split" ",$anidx;
      if($anx){ $an = ($an)? "$an,$anx" : $anx; }
      if($idx){ $ids= ($ids)? "$ids,$idx" : $idx; }
    } 
    # #? use both annvals, idclass? was elsif($nidclass)
    if($nidclass) { #  added from idclassf; use chrid or $ids      
      $ids=$cr unless($ids); #trick for cds.cov; see below ids and busco
      for my $id (split(",",$ids)) {
        if($id and my $anx= $idclassh->{$id} ) { unless($an =~ m/$anx/){ $an = ($an) ? "$an,$anx" : $anx; last; } }
      }
    }

    my $ok= ( $act >= $MINRD )?1:0; # skip? if too low
    next unless($ok);
    $nok++;
    
    # ucg_id_cov => gene_cov table : genecov{id}{cdsm,cdst,chrm,chrt}= @cov_vals (act,acm,acu?)
    # assume 20_000 to 50_000 gene ids, ave cds-length of 2_000/100 bins
    #   %genecov becomes hash-array of 50k x 4 x 20 vals =~ 4 mil vals? 10s of megabytes of mem?
    
    if($iscdscov) {
      # cr == geneid, all together
      my $ida= $cr;
      #? if($ida ne $lida) { putgene(xxx); }  
      #? gene_sum( $cr, $act,$acm,$acu);
      #x push @{$ucg_id_cov{$ida}}, $acm; # keep all? $act,$acm,$acu; rename ucg_id_cov > gene_id_cov
      push @{$ucg_id_cov{$ida}{cdsm}}, $acm; # keep all? $act,$acm,$acu; rename ucg_id_cov > gene_id_cov
      push @{$ucg_id_cov{$ida}{cdst}}, $act; # keep all? $act,$acm,$acu; rename ucg_id_cov > gene_id_cov
      $ucg_id_cov{$ida}{ncds}++;
      $lida= $ida;
      
    } else {
      if($ids and $an =~ /CDS/i) { # require CDS annot? 
        my($ida,@idc); 
        #UPD21apr26: ?? revert to count 1 id/bin locus, more skews gene sums : No, not that either
        # .. both ways are wrong: count all ids over-counts some lowqual id-hits, count 1 under-counts some paralogs
        # .. likely need better cds.anntab, separate best/2nd cds id aligns, count all best/top per bin-locus, not 2nds
        my @ids= split",",$ids; #? add ALL CDS ids ..uniq? ?
        if($nidclass) { @idc= grep{ $idclassh->{$_} =~ m/CDS/i } @ids; }
        else { @idc= grep{ $ucg_id_cov{$_} } @ids; } # this way only? assumes have cdscov, 1st read
        #old2: for $ida (@idc)
        #x ($ida) = shift @idc; # this way gives many more over-asm; 
        #x   .. arath18ccc_SRR10178325 casm<CU	1390 vs casm<CU	57 arath18cab_genesum
        #x ($ida) = sort @idc; # no change from unsort
        # if($ida) # UPD21apr26??
        for $ida (@idc)
        {
          push @{$ucg_id_cov{$ida}{chrm}}, $acm;  
          push @{$ucg_id_cov{$ida}{chrt}}, $act;  
          $ucg_id_cov{$ida}{nchr}++;
          $lida= $ida; # gene ids unordered on chrcov
        }
        
      # } else {        
      #   if($nidclass) { ($ida)= grep{ $idclassh->{$_} =~ m/CDS/i } sort @ids; } # busco|UCG
      #   $ida ||= $ids[0]; # else old way
      #   if($ida) {
      #     push @{$ucg_id_cov{$ida}{chrm}}, $acm;  
      #     push @{$ucg_id_cov{$ida}{chrt}}, $act;  
      #     $ucg_id_cov{$ida}{nchr}++;
      #     $lida= $ida; # gene ids unordered on chrcov
      #   }
      # }        
        
      }
    }

    if($lcr eq $cr and $lcb>0 and $cb > $lcb) { #?
      if(@cbin < kSAMPLE) { 
        my $bspan= $cb - $lcb; push @cbin, $bspan; 
        my $nbin=@cbin; $CBIN= $cbin[int($nbin/2)] if($nbin>100);
        } 
    }
    
    ($lcr,$lcb)= ($cr,$cb);
  } close($inh);

  warn "#readCovtab($covtab,cds=$iscdscov) nok=$nok, nin=$nin\n" if $debug;
  return($nok);
}

=item putGeneCovSum read genecov means table, write summary of classes
  
  sumclass ok, xcopy levels (<=0.5,1,>=2,>=3,..), locov, hicov, tooshort, ..

genemeans:
GeneID	CD_med	CD_ave	CD_nit	CD_sem	CD_rng	Ch_med	Ch_ave	Ch_nit	Ch_sem	Ch_rng	CDt_med	CDt_ave	Cht_med	Cht_ave
dropse20uc:g6903256t1	89	88.42	566	3.75	28,121	92	91.30	621	3.69	62,120	89	88.49	92	91.37
dropse20uc:g6900405t1	98	97.69	458	4.62	52,198	97	97.72	480	4.51	61,183	143	216.65	121	209.09
dropse20uc:g6902868t1	92	91.26	447	4.36	10,134	92	92.72	475	4.29	44,133	92	91.26	92	92.72
...
dropse20uc:g6898477t1	0	0	0	0	0	102	107.00	3	62.02	101,118	0	0	102	107.00
dropse20uc:g117184805t1	0	0	0	0	0	118	115.50	2	0.00	113,118	0	0	118	115.50
allgenes	88	109.06	14327	4.30	22,28894	95	98.03	14326	1.58	36,13623
  ^total now

=cut

sub diffratio {
  my($va,$vb,$mindif)=@_; $mindif||=0.66;
  my($cl,$r)=(0,0);
  return ($cl,$r) if($va<0.01 and $vb<0.01);

  if($vb < 0.01) { $r=($va>=1) ? 9 : 1; }
  else { $r= $va/$vb;  }
  $cl= ($r < 0.57)? -1 : ($r > 1.75) ? 1 : 0;
  # 1/1.75 =~ 0.57; 1/0.55 =~ 1.75
    
  ## old, inv -1  
  # my $inv=0; ($va,$vb,$inv)=($vb,$va,1) if($va>$vb);  
  # $r=$va/$vb;  $cl=($r < $mindif)?1:0; 
  # if($inv and $r > 0.001){ $r= 1/$r; $cl= -$cl; }
  
  #OFF: $r= sprintf"%.4f",$r; #? or not, not printed
  
  return ($cl,$r); 
}

sub putXcopy { 
  my($Oh, $id, @ctn)=@_;
  # @ctn= ($xcds,$txcds,$ncds, $xchr,$txchr,$ncr), only 2 cols here CDS,Chr
 
  sub xfmt{ my($v)=@_; my $d=($v>99)?0:($v>9)?1:2; return sprintf "%.${d}f",$v; }
  # my $XCUT=$ENV{xcut}|| 1.7; my $LCUT= $ENV{lcut}||0.55;
  my($hic,$hiv)=("","");   
  print $Oh $id; 
  while( my @c= splice(@ctn,0,3) ) {
    $hic .= ($c[1]>$XHICUT)?"2":($c[1]<$XLOCUT)?"0":"1"; 
    $hiv .= ($c[0]>$XHICUT)?"x":($c[0]<$XLOCUT)?"z":"n";  
    print $Oh "\t", join(",",xfmt($c[0]),xfmt($c[1]),$c[2]);  
    }
  # add counts of CApatt to outsum, hash it
  my $capatt="c$hic$hiv";
  print $Oh "\t",$capatt,"\n"; 
  return ($capatt);
}

sub putGeneCovSum {
  my($outsum, $outmeans)= @_;

  use constant { kIT => 0, kTOOSHORT => 1, kDIFFBN => 4, kDIFFLEN => 2, 
        kDIFFKU => 8, kDIFFKUn => 64, kDIFFKUC => 16, kDIFFKUCn => 32,
        kTOTKU =>128, kTOTKUn=>256 };
  my @klasses= (kTOOSHORT,kDIFFKU,kDIFFKUn,kDIFFKUC,kDIFFKUCn, kTOTKU,kTOTKUn, kDIFFBN,kDIFFLEN);
  my %klabels= ( 0 => 'item', 1 => 'short',4 => 'dLENbn', 2 => 'dLENcds', 
          8 => 'casm>CU', 64 => 'casm<CU', 16 => 'cds>CU', 32 => 'cds<CU',
          128 => 'ctot>CU', 256 => 'ctot<CU' );
  #? %klabels constant vals are not inserted in hash keys !!
  #? my %klabels= ( kIT => 'item', kTOOSHORT => 'tooshort', kDIFFBN => 'diffbn', kDIFFKU => 'diffku', kDIFFLEN => 'difflen' );
  
  my $MINBN=2; # min bins, ie 2 x 100b
  my $pDIFFBN=0.66; # cdsbins <> chrbins
  my $pDIFFKU=0.66; # cov/KU 

  my $KU= $KU_OPT || $samvals{$asmid}{kucg} || $samvals{$masmid}{kucg}; # should be global -option
  unless($KU) { $KU=  _minnot0($samvals{$asmid}{kulo}, $samvals{$asmid}{kuhi} );   }      # or _max ?
  unless($KU){ 
    open(F,"tail $outmeans |"); while(<F>){ my @v= split; 
    if($v[0] eq $TOTALID){ $KU= _max($v[1], $v[6]);  } } close(F);
  }
  
  # add outxcovtab => name.genemeans => name.genexcopy * or .genexcovtab ?
  (my $outxcopy=$outmeans) =~ s/.genemeans//; $outxcopy.=".genexcopy";
  warn "#putGeneCovSum($outmeans => $outsum, $outxcopy), asmid=$asmid, KU=$KU\n" if $debug;

  open(F,$outmeans) or return;
  rename($outxcopy,"$outxcopy.old") if(-f $outxcopy); 
  open(my $xouth,'>',$outxcopy) or die "writing $outxcopy";  

  # for putXcopy(name.genexcopy)
  my @XCOLS=("GeneID","CDS_xC,tC,nB","Asm_xC,tC,nB","CApatt");
  #above# my $XHICUT=$ENV{xcut}|| 1.7; my $XLOCUT= $ENV{lcut}||0.55;

  print $xouth join("\t",@XCOLS)."\n";

  my($nit);
  my(%class,%sums,%capatt,%capats,@hd,@tcds,@tchr);
  
  while(<F>) {
    next if(/^\W/);
    my @v=split;
    my $id=$v[0];
    if($id =~ /^GeneID/){ @hd=@v; next; }
    elsif($id eq $TOTALID) { 
     @tcds=@v[1..5]; @tchr= @v[6..10]; 
     unless($KU) { $KU= _max($v[1], $v[6]); } # median of medians, chr val most accurate
     next; } # last in?
    
    $nit++;
    my($cdsm,$cdsa,$cdsn,$cdse,$cdsr)= @v[1..5];
    my($crm,$cra,$crn,$cre,$crr)= @v[6..10]; # may be missing?
    my($cdsrt,$crrt)= @v[11,13];
    # for xcopy: id, cdsm/KU,cdsrt/KU,cdsn,  crm/KU, crrt/KU, crn, cpatt: c11nn/c22xn/...
    # FIXME: cpatt==c00zz must ignore zz = xchr,xcds
    
    my $cdsw= ($ntcds) ? $cdslen->{$id} : 0;
    
    my $cl=0; my $tooshort=0;
    if($cdsn < $MINBN and $crn < $MINBN) { $cl |= kTOOSHORT; $tooshort=1; }
    if($tooshort) { } # skip other classes?

    my($clwc,$rwc,$clwa,$rwa)= (0) x 9;
    if($cdsw and not $tooshort) { my $cdswb= $cdsw/$CBIN; 
      ($clwc,$rwc)= diffratio(  $cdsn, $cdswb, $pDIFFBN); 
      ($clwa,$rwa)= diffratio(  $crn, $cdswb, $pDIFFBN); 
      if($clwa or $clwc) { $cl |= kDIFFLEN; }
      push @{$sums{'wchr'}}, $rwa;
      # compare cds-length to readmap len, both cds, chr? measures dups and underasm ?
      }

    # for KU set, count <1 and >1
    my($cldif,$pdifbin)= ($tooshort)?(0,1):diffratio($cdsn, $crn, $pDIFFBN);
    if($cldif) { $cl |= kDIFFBN; } # cov span differs, pdifbin < 1 means crn larger, has dups, pdifbin > 1 means crn missing cds?

    # tCopy vals: add total cov counts, both cdsrt, crrt? or just one?
    # do these before xchr,xcds, if txcds & txchr < LOCUT, skip xchr,xcds test
    my($cltcds,$txcds)=  diffratio($cdsrt, $KU, $pDIFFKU);
    my($cltchr,$txchr)=  diffratio($crrt, $KU, $pDIFFKU);
    if($cltchr<0) { $cl |= kTOTKUn; } # xcov != 1, for chrasm , 
    elsif($cltchr>0) { $cl |= kTOTKU; } # xcov != 1, for chrasm

    ## xCopy vals
    my($noxcopy)= ($tooshort or ($cltchr<0 and $cltcds<0) )?1:0;
    my($clxc,$xcds)= ($noxcopy)?(0,1) : diffratio($cdsm, $KU, $pDIFFKU);
    my($clxa,$xchr)= ($noxcopy)?(0,1) : diffratio($crm, $KU, $pDIFFKU);
    if($clxc<0) { $cl |= kDIFFKUCn; } elsif($clxc>0) { $cl |= kDIFFKUC; } # xcov != 1, for cds
    if($clxa<0) { $cl |= kDIFFKUn; }  elsif($clxa>0) { $cl |= kDIFFKU; } # xcov != 1, for chrasm

    # cds <> chr xcov, ignore this?
    my($clxb,$xcds_chr)=  diffratio($xcds, $xchr, $pDIFFKU); 
    
    # want counts of xcopy per gene: on cds and on chr, count too many, too few, too short
    # sum classes here? ave,medn for xchr,xcds, diff(cds,chr)
    $class{0}++;
    for my $ic (@klasses) { if($cl & $ic) { $class{$ic}++; }  }
    push @{$sums{'difbin'}}, $pdifbin; # was difbin
    push @{$sums{'xcds'}}, $xcds;
    push @{$sums{'xchr'}}, $xchr;
    push @{$sums{'ycds_cr'}}, $xcds_chr;
    push @{$sums{'tchr'}}, $txchr;
    
    #add for xcopy: id, xcds=cdsm/KU,txcds=cdsrt/KU,cdsn,  xchr=crm/KU, txchr=crrt/KU, crn, cpatt: c11nn/c22xn/...
    my($capatt)= putXcopy( $xouth, $id, $xcds, $txcds, $cdsn, $xchr, $txchr, $crn);
    
    $capatt{$capatt}++;
    # for only 4 patt sym, add totals c.2../c.1../c.0.., c...x/c...n/c...z ? 
    # should == casm>CU,casm<CU, ctot>CU,ctot<CU
    # FIXME: c00zz must ignore zz
    my($cap); # capats separate counter from capatt
    ($cap=$capatt) =~ s/c...(.)/c...$1/; $capats{$cap}++; #casm = z/n/x
    ($cap=$capatt) =~ s/c.(.)../c.$1../; $capats{$cap}++; #ctot = 0/1/2
  }
  
  rename($outsum,"$outsum.old") if(-f $outsum); 
  open(my $outh,'>',$outsum) or die "writing $outsum";  
  print $outh "SUMMARY of $outmeans, asmid=$asmid \n";
  print $outh "Params: CU=$KU, xCopy(Cov/CU)<>1.0: hi=$XHICUT, lo=$XLOCUT\n";
 
  #? make 2 cols, add capatt right of Class_counts? max=10 rows
  print $outh "Class_Counts__________\n";
  
  my $it=0;
  for my $ic (kIT,@klasses) {
    my $c= $class{$ic}||0; 
    my $lb= $klabels{$ic}||"na$ic"; # why bad here? hash{kIT} == {'kIT'} not {0}
    print $outh "$lb\t$c"; print $outh ($it % 2 == 0) ? "\t" : "\n"; $it++; 
  } print $outh "\n" if($it % 2 == 1); 
  
  my %clabs=(
    'c.1..' => "\tasm-one-copy",
    'c.2..' => "\tasm-two+copy",
    'c.0..' => "\tasm-miss-copy",
    'c...x' => "\tasm-miss-dup",
    'c...z' => "\tasm-extra-dup",
    'c...n' => "\tasm-norm-dup",
  );
  my @caprows=(); 
  my @cap= sort{ $capatt{$b} <=> $capatt{$a} or $a cmp $b} keys %capatt;
  for my $i (0 .. $#cap) { my $c=$cap[$i]; push @caprows, "$c\t".$capatt{$c}; }
  @cap= sort keys %capats; # put last? first? always only 6?
  for my $i (0 .. $#cap) { my $c=$cap[$i]; my $cl=$clabs{$c}||""; push @caprows, "$c\t".$capats{$c}.$cl; }
  
  my $nct=@caprows;  
  print $outh "\nxCopypat_Counts cpatt:tCDS,tAsm,xCDS,xAsm\n";  
  my($ir,$ic)= ($nct>36)?(9,5) :($nct>27)?(9,4) :($nct>18)?(9,3) : ($nct>9)?(9,2) : ($nct,1);
  for my $j (1 .. $ir) { 
    my @row=();
    for my $k (1 .. $ic) {
      my $kk= $ir * ($k-1) + ($j-1);
      push @row, $caprows[ $kk ] unless($kk>=$nct); 
    }
    print $outh join("\t",@row)."\n";
  }
  
  
  print $outh "\nMedian,Mean,N,StDev,Range________\n";
  for my $sk (sort keys %sums){
    my @v=  @{$sums{$sk}};
    my @stat= meanvar(@v); # $med,$ave,$nit,$se,$range
    print $outh join("\t",$sk,@stat)."\n";
  }
  
  close($outh);
}


sub putGeneCovStats {
  my($outmeans)= @_;
  rename($outmeans,"$outmeans.old") if(-f $outmeans); 
  open(my $outmeanh,'>',$outmeans) or die "writing $outmeans";  
  
  # note: CBIN * n = bases, add Mb vals?
  my($nput,@cdsmd,@chrmd)=(0);
  my @ids= sort{ $ucg_id_cov{$b}{ncds} <=> $ucg_id_cov{$a}{ncds} or $a cmp $b } keys %ucg_id_cov; # by cov size
  my ($didhdr);
  for my $ida (@ids) {
    my $ncds= $ucg_id_cov{$ida}{ncds};
    my $nchr= $ucg_id_cov{$ida}{nchr};
    my @cdsm= ($ncds) ? @{$ucg_id_cov{$ida}{cdsm}} : ();
    my @cdsstat= meanvar(@cdsm); # $med,$ave,$nit,$se,$range
    my @cdsts  = meanvar(($ncds) ? @{$ucg_id_cov{$ida}{cdst}} : ());  
    
    my @chrm= ($nchr) ? @{$ucg_id_cov{$ida}{chrm}} : ();
    my @chrstat= meanvar(@chrm); # $med,$ave,$nit,$se,$range
    my @chrts= meanvar(($nchr) ? @{$ucg_id_cov{$ida}{chrt}} : ());  

    push @cdsmd, $cdsstat[0] if($ncds);
    push @chrmd, $chrstat[0] if($nchr);
    
    unless($didhdr++) {
    my @stath= qw(med ave nit std rng); # std was sem
    my @gstat= map{ "CD_".$_ } @stath;
    my @cstat= map{ "Ch_".$_ } @stath;
    my @gstt = map{ "CDt_".$_ } @stath[0,1];
    my @cstt = map{ "Cht_".$_ } @stath[0,1];
    print $outmeanh join("\t","GeneID",@gstat,@cstat,@gstt,@cstt)."\n" ;
    }
    print $outmeanh join("\t",$ida,@cdsstat,@chrstat,@cdsts[0,1],@chrts[0,1])."\n"; $nput++;
    
    # calc xCopy? need chr,cds sum val?
    # how many xcopy in @cds, how many in @chr, diff?
  }
  
  my @cdstt= meanvar(@cdsmd);
  my @chrtt= meanvar(@chrmd);
  print $outmeanh join("\t",$TOTALID,@cdstt,@chrtt)."\n";
  
  close($outmeanh);
  return($nput);
}

#-----

sub _max{ return($_[0] < $_[1])?$_[1]:$_[0]; }
sub _min{ return($_[0] > $_[1])?$_[1]:$_[0]; }
sub _minnot0{ return ($_[0] == 0) ? $_[1] : ($_[1] == 0) ? $_[0] : _min(@_); }

sub sum_ucg_cov { # is this right for all genes? from UCG_KU_SPECIAL
  my($uid, $ucovs)= @_;
  my $ncov= @$ucovs;
  if($ncov>2) {  
    my @cu= sort{ $b<=>$a } @$ucovs;
    my $cmed= $cu[ int($ncov/2) ];
    push @ucg_med_cov, $cmed; # global median depth for all UCG
    push @ucg_med_ids, $uid; # ? not used, or check for dups?
    $nucg=@ucg_med_cov;
    # if(UCG_KU_MEANS) {  
    #  filt_sum( UCG_KU_MEANS, $cmed,$cmed,$cmed); # ?? $act,$acm,$acu
    # }    
  }
}

sub meanvar { 
  my(@vals)= @_;
  my($med,$nit,$ave,$var,$sd,$se,$s,$range)= (0) x 9;
  $nit= @vals; 
  return($med,$ave,$nit,$se,$range) if($nit<1);
  
  my @sval= sort {$a <=> $b} @vals;
  $med= $sval[ int($nit/2) ]; 
  
  my $nitall=$nit; #??
  if($TRIMAVE and $nit>9) {
    my $lo= $nit * 0.10; my $hi= $nit * 0.90;
    @sval= splice( @sval, $lo, $hi-$lo);
    $nit= @sval;
  }
  for my $i (0..$#sval) { my $v= $sval[$i]; $s += $v; $var += $v*$v; }
  $ave=$s/$nit; if($nit>2){ $var=($var - $ave*$ave)/($nit-1); $sd=sqrt($var); $se=$sd/sqrt($nit); }

  my($rlo,$rhi)= @sval[0,-1]; # or use 1st/3rd quart?
  # $range="$sval[0],$sval[-1]"; # or 1/3 quartile?  $nit*.25, $nit*.75 ?
  map{ my $d=($_ > 99)?0:($_ > 9)?1:2; $_= sprintf"%.${d}f",$_; $_ =~ s/\.0+$//; } ($med,$ave,$sd,$rlo,$rhi);
  
  return($med,$ave,$nit,$sd,"$rlo,$rhi"); # se > $sd ; changed
}

# sub meanvar_orig { 
#   my($outh)= @_;
#   my($med,$nit,$ave,$var,$sd,$se,$s);
#   #x print $outh "$title stats for nt=$nt\n"; 
# 
#   my @hdr=qw(Mean SEM Nitem StDev); 
#   #? new cols: Median Mean SEM Nitem StDev Mode Sum|Var|Range
#   # if($DOHIST) { unshift(@hdr,"Median");  push(@hdr,"Mode"); } 
#   if($DOMED) { unshift(@hdr,"Median"); }
#   push(@hdr,($showSUM)?"Sum":"Var"); # add Range?
#   
#   my $fmt="%.2f\t%.2f\t%d\t%.2f";
#   # if($DOHIST) { $fmt="% 4d\t".$fmt."\t%4d"; } 
#   if($DOMED) { $fmt="% 4d\t".$fmt; }
#   $fmt.=($showSUM)?"\t%d":"\t%.1f"; 
#    
#   #x print $outh join("\t","Item  ",@hdr)."\n"; 
#   for my $i (0..$#sum) { # data cols, drop this? 
#     my $litem= $lvar[$i]||"it$i"; #< @lvar == stat column head
#     $s=$sum[$i]; $nit=$nit[$i]; 
#     ($ave,$var,$sd,$se)=(0) x 4;
#     if($nit>3){ $ave=$s/$nit; $var=($svar[$i] - $ave*$ave)/($nit-1); $sd=sqrt($var); $se=$sd/sqrt($nit); }
#     my @mv=($ave,$se,$nit,$sd); 
#     if($DOMED) { my @ss=sort{$a <=> $b} @{$sval[$i]}; $med= $ss[ int($nit/2) ];  unshift(@mv,$med); }
#     #skip DOHIST: 
#     push @mv,(($showSUM)?$s:$var);
#     
#     #x printf $outh "%-6s\t$fmt\n",$litem,@mv;
#     } 
#   #x print $outh "\n"; #------ ?
# }

sub readAnnots {
  my($anntable)= @_;
  my($nann,$lpt)=(0,0); 
  warn "#read anntable=$anntable\n" if($debug);
  open(F,$anntable); 
  while(<F>){
    next if(/^\W/); my @v= split;
    my($cid,$cib,$an,$ids)= @v; # any more than an,ids ?
    next unless($an or $ids);
    $annvals{$cid}{$cib}="$an\t$ids"; $nann++;
  } close(F);
  return($nann);
}


sub read_idclass { # for case of no anntable, but idclass, eg CDS.covtab
  my($crclassf,$allclasses,$crlenh)=@_; 
  my($nid,$ncl,$ok,$inh)=(0,0,0);
  my %crclass=();   my @crclass=();
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
  if($nid==0 and $CRTPAT and ref($crlenh)){
    for my $cr (sort keys %$crlenh) {
      my($crt)= ($cr =~ m/($CRTPAT)/)? $1 : 'UNK'; 
      $crclass{$cr}= $crt;
    }
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
    warn "#ERR too many idclasses n=$ncl from nid=$nid of $crclassf or sam ids x CRTPAT='$CRTPAT' \n";
    $nid=0; %crclass=(); @crclass=();
  }
  warn "# read nid=$nid, nclass=$ncl from $crclassf\n" if($debug);
  return($nid,\%crclass,\@crclass);
}


sub readMetad {
  my($sampledata)= @_;
  my($nsam,$aid)=(0,0); #? only 1 aid, should have default
  if($sampledata and -f $sampledata) {
    warn "#read sampledata=$sampledata\n" if($debug);
    open(F,$sampledata); # 
    while(<F>){
      next if(/^\W/);
      my($key,$val)= (m/^(\w+)\s*=\s*(.+)$/) ? ($1,$2):(0,0);
      next unless($key);
      
      if($key eq 'pt' or $key eq 'asmid') { $aid=$val; $nsam++;
        map{ $samvals{$aid}{$_}= 0; } qw(flowcyto cytomb asmtotal atotalmb asmname );
    
      } elsif($key eq 'flowcyto') { # /^flowcyto=(.+)$/)  
        # my $flowcyto= $val; # flowcyto=234-391 Mb
        $samvals{$aid}{$key}= $val;
        $samvals{$aid}{cytomb}  = ($val=~m/(\d+)/)?$1:0; # can be range 160-180
    
      } elsif($key eq 'glncformula') { # /^glncformula=(.+)$/ 
        $samvals{$aid}{$key}= $val; # glncformula=nnn Mb
        
      } elsif($key eq 'asmtotal') { # /^asmtotal=(.+)$/)   
        my $asmtotal= $val; 
        $samvals{$aid}{$key} = $asmtotal;
        $samvals{$aid}{atotalmb} = ($asmtotal=~m/(\d+)/)?$1:0; # can be range 160-180
    
      } elsif($key eq 'kulo' or $key eq 'kuhi' or $key eq 'asmname'){ 
        $val=~s/\s.*//; $samvals{$aid}{$key}= $val;  # /^kulo=(\S+)/)  # test hi<lo
        
      } else { # save all valid key = val
        $val=~s/\s.*//; $samvals{$aid}{$key} = $val;
      }
      
    } close(F);
  }
  
  # allow ENV{key} vals ?
  for my $key (qw(asmid asmname asmtotal flowcyto kulo kuhi kucg)) {
    if(my $ev= $ENV{$key}) { 
      if($key eq "asmid"){ $aid=$ev; $nsam++; } 
      $samvals{$aid}{$key}= $ev if($aid); 
      }
  }
  
  return($nsam,$aid); # return ($nsam, \%metavals); #?
}

 
sub openRead {
  my($fn)=@_; my($ok,$inh)=(0); 
  return(0,undef) unless($fn and -f $fn);
  if($fn =~ /\.gz/) { $ok= open($inh,"gunzip -c $fn |"); } else { $ok= open($inh,$fn); }
  # die "cant read $fn" unless ($ok); # warn? leave to caller
  return($ok,$inh);
}

__END__

=item xcopy sum try1n

$evigene/scripts/genoasm/gnodes_genescov.pl -debug -title dropse20cab -asmid dropse20chrs -metad dropse20chrs.metad \
 -idclass dropse20cdste.idclass -ann dropse20chrs.fa.anntab -means dropse20cdschr_SRR11813283.genemeans
 
#genecov output to sum=dropse20cab_genesum.txt means=dropse20cdschr_SRR11813283.genemeans
#read sampledata=dropse20chrs.metad
#read anntable=dropse20chrs.fa.anntab
# read nid=15225, nclass=4 from dropse20cdste.idclass
#putGeneCovSum(dropse20cdschr_SRR11813283.genemeans => dropse20cab_genesum.txt, dropse20cdschr_SRR11813283.genexcopy), asmid=dropse20chrs, KU=95

dropse20cab_genesum.txt
SUMMARY of dropse20cdschr_SRR11813283.genemeans, asmid=dropse20chrs, KU=95
Class_Counts__________
item    14330   short   0
casm>CU 21      casm<CU 33
cds>CU  268     cds<CU  120
ctot>CU 851     ctot<CU 4
dLENbn  1894    dLENcds 0

xCopypat_Counts cpatt:tCDS,tAsm,xCDS,xAsm
c11nn   13187   c22zn   9       c...n   14280   asm-norm-dup
c22nn   539     c11zn   8       c...x   21      asm-miss-dup
c22xn   171     c12zn   4       c...z   29      asm-extra-dup
c12nn   140     c10nz   3       c.0..   4       asm-miss-copy
c21xn   88      c11nz   3       c.1..   13441   asm-one-copy
c21nn   69      c00zz   1       c.2..   885     asm-two+copy
c01zn   65      c01zz   1
c11zz   20      c12nx   1
c22xx   20      c22xz   1

Median,Mean,N,StDev,Range________
difbin  0.84    0.78    14330   0.83    0,9
tchr    1       1.33    14330   4.29    0,294
xcds    0.93    1.15    14330   5.41    0,304
xchr    1       1.03    14330   1.99    0,143
ycds_cr 0.94    1.06    14330   2.36    0,140

dgmm:gnodesf:% head dropse20cdschr_SRR11813283.genexcopy
GeneID  CDS_xC,tC,nB    Asm_xC,tC,nB    CApatt
dropse20uc:g6903256t1   0.94,0.94,566   0.97,0.97,621   c11nn
dropse20uc:g6900405t1   1.03,1.51,458   1.02,1.27,480   c11nn
dropse20uc:g6902868t1   0.97,0.97,447   0.97,0.97,475   c11nn
dropse20uc:g6899878t1   1.11,1.11,273   1.14,1.14,314   c11nn

evigene/scripts/genoasm/gnodes_genescov.pl -debug -title arath18cab -asmid=arath18tair_chr -metad=arath20asm.metad     -idclass arath18tair1cds.idclass -ann arath18tair_chr.anntab -means arath18cdschr_SRR10178325.genemeans
#genecov output to sum=arath18cab_genesum.txt means=arath18cdschr_SRR10178325.genemeans
#read sampledata=arath20asm.metad
#read anntable=arath18tair_chr.anntab
# read nid=27444, nclass=2 from arath18tair1cds.idclass
#putGeneCovSum(arath18cdschr_SRR10178325.genemeans => arath18cab_genesum.txt, arath18cdschr_SRR10178325.genexcopy), asmid=arath18tair_chr, KU=18

arath18cab_genesum.txt
SUMMARY of arath18cdschr_SRR10178325.genemeans, asmid=arath18tair_chr, CU=18
Class_Counts__________
item    27406   short   25
casm>CU 120     casm<CU 57
cds>CU  345     cds<CU  1433
ctot>CU 506     ctot<CU 51
dLENbn  9015    dLENcds 0

xCopypat_Counts cpatt:tCDS,tAsm,xCDS,xAsm
c11nn   25512   c01nn   30      c10zz   4       c...n   27249   asm-norm-dup
c01zn   842     c00nn   28      c22nx   4       c...x   126     asm-miss-dup
c22nn   220     c11zz   13      c10nn   3       c...z   31      asm-extra-dup
c21xn   216     c12nx   13      c11nz   3       c.0..   43      asm-miss-copy
c12nn   113     c12zn   13      c20xz   2       c.1..   26818   asm-one-copy
c11zn   110     c22zn   13      c22nz   2       c.2..   545     asm-two+copy
c22xx   104     c10nz   6       c01zz   1
c21nn   91      c02zn   5       c12zx   1
c22xn   53      c02zx   4       

Median,Mean,N,StDev,Range________
difbin  0.75    0.71    27406   0.82    0,9
tchr    1.06    1.24    27406   9.3     0,790
xcds    0.89    1.05    27406   7.51    0,758
xchr    1       1.10    27406   4.24    0,395
ycds_cr 0.89    0.91    27406   1.65    0,154


arath18cdschr_SRR10178325.genexcopy
GeneID  CDS_xC,tC,nB    Asm_xC,tC,nB    CApatt
AT1G67120t1     0.89,0.89,162   1.06,1.06,204   c11nn
AT3G02260t1     1.00,1.00,153   1.00,1.00,160   c11nn
AT5G23110t1     0.94,0.94,142   1.00,1.00,147   c11nn
AT4G17140t1     0.89,0.89,127   1.00,1.00,156   c11nn
AT1G48090t1     0.89,0.89,126   1.00,1.00,161   c11nn
AT2G17930t1     0.94,0.94,116   1.06,1.06,269   c11nn
AT4G36080t1     1.00,1.00,115   1.06,1.06,262   c11nn
AT1G55860t1     0.89,1.06,111   1.00,1.06,231   c11nn
AT1G70320t1     1.00,1.06,110   1.00,1.06,230   c11nn

=cut

=item eg

  need: cds.covtab, chr.covtab (both?)
        anntab for CDS, ids, other ann?
  maybe: metad, idclass (drop idclass for anntab?)
  minread=0 now default
  
  dropse20gnodes/
  $evigene/scripts/genoasm/gnodes_genescov.pl -debug  \
    -title dropse20cdschr -asmid dropse20chrs -metad dropse20chrs.metad \
    -idclass dropse20cdste.idclass -ann dropse20chrs.fa.anntab \
    -cdscov dropse20t1cds_SRR11813283_1_bwa.cdschr7b.covtab -chrcov dropse20chrs_SRR11813283_1_bwa.cdschr7b.covtab

#genecov output to sum=dropse20cdschr_SRR11813283_genesum.txt means=dropse20cdschr_SRR11813283.genemeans
#read sampledata=dropse20chrs.metad
#read anntable=dropse20chrs.fa.anntab
# read nid=15225, nclass=4 from dropse20cdste.idclass
#readChrtab(dropse20t1cds_SRR11813283_1_bwa.cdschr7b.chrtab) nok=237192
#readCovtab(dropse20t1cds_SRR11813283_1_bwa.cdschr7b.covtab,cds=1) nok=237192, nin=237192
#readCovtab(dropse20chrs_SRR11813283_1_bwa.cdschr7b.covtab,cds=0) nok=1632808, nin=1632808
#putGeneCovSum(dropse20cdschr_SRR11813283_genesum.txt, dropse20cdschr_SRR11813283.genemeans), asmid=dropse20chrs, KU=95

  dromel20gnodes/
  # dromel6rel_ncbigff.anntab	drosmel6gnodncbi.anntab		drosmel6ref_chr.anntab : want ncbigff annots to match pseudog w/ cds dups
  # need a,b reads of SRR11460805, fixme to merge 2+ covtabs
  $evigene/scripts/genoasm/gnodes_genescov.pl -debug \
    -title drosmel6cdschr -asmid drosmel6ref_chr -metad drosmelchr.metad \
    -idclass dromel6relt1cds.idclass -ann drosmel6ref_chr.anntab \
    -cdscov dromel6relt1cds_SRR11460805_1_bwa.cc7d.covtab -chrcov drosmel6ref_chr_SRR11460805_1_bwa.cc7d.covtab
** FIXME: KU=14 too little data, add merged srra,b cds.covtab , have merge.chr.covtab

#genecov output to sum=drosmel6cdschr_SRR11460805_genesum.txt means=drosmel6cdschr_SRR11460805.genemeans
#read sampledata=drosmelchr.metad
#read anntable=drosmel6ref_chr.anntab
# read nid=13954, nclass=2 from dromel6relt1cds.idclass
#readChrtab(dromel6relt1cds_SRR11460805_1_bwa.cc7d.chrtab) nok=219889
#readCovtab(dromel6relt1cds_SRR11460805_1_bwa.cc7d.covtab,cds=1) nok=219889, nin=219889
#readCovtab(drosmel6ref_chr_SRR11460805_1_bwa.cc7d.covtab,cds=0) nok=1437316, nin=1437316
#putGeneCovSum(drosmel6cdschr_SRR11460805_genesum.txt, drosmel6cdschr_SRR11460805.genemeans), asmid=drosmel6ref_chr, KU=14

  aweed20gnodes/
  $evigene/scripts/genoasm/gnodes_genescov.pl -debug \
    -title arath18cdschr -asmid=arath18tair_chr -metad=arath20asm.metad \
    -idclass arath18tair1cds.idclass -ann arath18tair_chr.anntab \
    -cdscov arath18tair1cds_SRR10178325_1_bwa.cc7d.covtab -chrcov arath18tair_chr_SRR10178325_1_bwa.cc7d.covtab 

#genecov output to sum=arath18cdschr_SRR10178325_genesum.txt means=arath18cdschr_SRR10178325.genemeans
#read sampledata=arath20asm.metad
#read anntable=arath18tair_chr.anntab
# read nid=27444, nclass=2 from arath18tair1cds.idclass
#readChrtab(arath18tair1cds_SRR10178325_1_bwa.cc7d.chrtab) nok=337962
#readCovtab(arath18tair1cds_SRR10178325_1_bwa.cc7d.covtab,cds=1) nok=337962, nin=337962
#readCovtab(arath18tair_chr_SRR10178325_1_bwa.cc7d.covtab,cds=0) nok=1196689, nin=1196689
#putGeneCovSum(arath18cdschr_SRR10178325_genesum.txt, arath18cdschr_SRR10178325.genemeans), asmid=arath18tair_chr, KU=18

  apis20gnodes/
  $evigene/scripts/genoasm/gnodes_genescov.pl -debug     -title apis19cdschr -asmid apismel19hav31asm \
    -metad apismel20asm.metad  -idclass amel20cdste.idclass -ann apismel19hav31asm.anntab \
    -cdscov apismel14evg3cds_SRR9108936_1mc_bwa.cc7c.covtab -chrcov apismel19hav31asm_SRR9108936_1mc_bwa.cc7c.covtab

#genecov output to sum=apis19cdschr_SRR9108936_genesum.txt means=apis19cdschr_SRR9108936.genemeans
#read sampledata=apismel20asm.metad
#read anntable=apismel19hav31asm.anntab
# read nid=375870, nclass=5 from amel20cdste.idclass
#readChrtab(apismel14evg3cds_SRR9108936_1mc_bwa.cc7c.chrtab) nok=409891
#readCovtab(apismel14evg3cds_SRR9108936_1mc_bwa.cc7c.covtab,cds=1) nok=409891, nin=409891
#readCovtab(apismel19hav31asm_SRR9108936_1mc_bwa.cc7c.covtab,cds=0) nok=2252366, nin=2252366
#putGeneCovSum(apis19cdschr_SRR9108936_genesum.txt, apis19cdschr_SRR9108936.genemeans), asmid=apismel19hav31asm, KU=24

=cut

=item putGeneCovSum try2...

  .. try sort by total cov: CDt_med or Cht_med, count/class by dups-tcov, single-tcov, low-tcov?
  eg. dropse has ~30-40 diff ku for Ch_med, but 1000+ in Cht_med >= 2 xcopy,
  most w/ Chr_med >= 2x have proper 1x cov on chrasm, expected from genome-total xCopy =~ 1
  .. compare w/ arabid, dromel where xCopy ~> 1 w/ some extra missing dups

drosmel6cdschr_SRR11460805_genesum.txt
SUMMARY of drosmel6cdschr_SRR11460805.genemeans, asmid=drosmel6ref_chr, KU=14
Class_Counts__________
item    13946
short   26
casm>CU 155
casm<CU 117
cds>CU  161
cds<CU  346
ctot>CU 340  *low? check means.tab
ctot<CU 45
dLENbn  639
dLENcds 621

Median,Mean,N,StDev,Range________
difbin  0.92    0.90    13946   0.94    0,9
tchr    1       1.92    13946   9.9     0,142
wchr    1.08    1.47    13902   3.98    0,103
xcds    0.93    0.95    13946   1.05    0,22.9
xchr    1       1.04    13946   1.12    0,8
ycds_cr 0.93    0.93    13946   1       0,10.5


dropse20cdschr_SRR11813283_genesum.txt
SUMMARY of dropse20cdschr_SRR11813283.genemeans, asmid=dropse20chrs, KU=95
Class_Counts__________
item    14330
short   0
casm>CU 21
casm<CU 33
cds>CU  268
cds<CU  120
ctot>CU 851  *? low? check means.tab
ctot<CU 4
dLENbn  1894
dLENcds 1891

Median,Mean,N,StDev,Range________
difbin  0.84    0.78    14330   0.83    0,9
tchr    1       1.33    14330   4.29    0,294
wchr    1.19    2.14    14327   7.78    0,535
xcds    0.93    1.15    14330   5.41    0,304
xchr    1       1.03    14330   1.99    0,143
ycds_cr 0.94    1.06    14330   2.36    0,140


arath18cdschr_SRR10178325_genesum.txt
SUMMARY of arath18cdschr_SRR10178325.genemeans, asmid=arath18tair_chr, KU=18
Class_Counts__________
item    27406
short   25
casm>CU 120
casm<CU 95
cds>CU  345
cds<CU  1496
ctot>CU 506    *??? this seems low, check means tab
ctot<CU 51
dLENbn  9040
dLENcds 8939

Median,Mean,N,StDev,Range________
difbin  0.75    0.71    27406   0.81    0,9
tchr    1.06    1.24    27406   9.3     0,790
wchr    1.33    2.01    27305   7.90    0,751
xcds    0.89    1.05    27406   7.50    0,758
xchr    1       1.09    27406   4.24    0,395
ycds_cr 0.89    0.91    27406   1.65    0,154


apis19cdschr_SRR9108936_genesum.txt
SUMMARY of apis19cdschr_SRR9108936.genemeans, asmid=apismel19hav31asm, KU=24
Class_Counts__________
item    67671
short   8
casm>CU 537     * higher than others, maybe due to TE-cds in evg genes?
casm<CU 2735
cds>CU  1030
cds<CU  14133
ctot>CU 866     * similar to dropse
ctot<CU 2645    * uses apis14evg genes, many lack DNA match, spurious(contam ?)
dLENbn  4240
dLENcds 4081

Median,Mean,N,StDev,Range________
difbin  0.91    0.97    67671   1.32    0,14.8
tchr    1       1.25    67671   16.2    0,2911
wchr    1.10    2.47    67504   25.7    0,1862
xcds    0.79    1.17    67671   18.8    0,4037
xchr    1       1.03    67671   1.45    0,123
ycds_cr 0.82    0.96    67671   5.17    0,673

$evigene/scripts/genoasm/gnodes_genescov.pl -debug     -title daphcar20cdschr -asmid daphcari20chr -metad daphcar20asm.metad -idclass daphsim17evgt1cds.idclass -ann daphcari20chr.anntab   -cdscov daphsim17evgt1cds_SRR10389283_1_bwa.cc7c.covtab -chrcov daphcari20chr_SRR10389283_1_bwa.cc7c.covtab
#genecov output to sum=daphcar20cdschr_SRR10389283_genesum.txt means=daphcar20cdschr_SRR10389283.genemeans
#read sampledata=daphcar20asm.metad
#read anntable=daphcari20chr.anntab
# read nid=38523, nclass=2 from daphsim17evgt1cds.idclass
#readChrtab(daphsim17evgt1cds_SRR10389283_1_bwa.cc7c.chrtab) nok=303332
#readCovtab(daphsim17evgt1cds_SRR10389283_1_bwa.cc7c.covtab,cds=1) nok=303332, nin=303332
#readCovtab(daphcari20chr_SRR10389283_1_bwa.cc7c.covtab,cds=0) nok=1313507, nin=1313507
#putGeneCovSum(daphcar20cdschr_SRR10389283_genesum.txt, daphcar20cdschr_SRR10389283.genemeans), asmid=daphcari20chr, KU=50

 daphcar20cdschr_SRR10389283_genesum.txt
SUMMARY of daphcar20cdschr_SRR10389283.genemeans, asmid=daphcari20chr, KU=50
Class_Counts__________
item    37709
short   19
casm>CU 240     * highish, lower than apis
casm<CU 1272
cds>CU  2856    * high
cds<CU  5277
ctot>CU 1887    * high vs insects,aplant
ctot<CU 535
dLENbn  9087
dLENcds 8543

Median,Mean,N,StDev,Range________
difbin  0.75    0.76    37709   1.12    0,12
tchr    1.02    1.48    37709   9.6     0,1220
wchr    1.33    7.86    37146   71.1    0,9573
xcds    0.84    2.19    37709   105     0,19432
xchr    1       1.18    37709   5.63    0,214
ycds_cr 0.84    1.81    37709   63.2    0,9716


pt=daphpulex_pa42v2; 
$evigene/scripts/genoasm/gnodes_genescov.pl -debug     -title daphplx19cdschr -asmid $pt \
 -metad daphplx20chrs.metad -idclass daphplx17evgt1m_cds.idclass  -ann $pt.anntab \
 -cdscov daphplx17evgt1m_cds_SRR3090572_1_bwa.cc7c.covtab  -chrcov ${pt}_SRR3090572_1_bwa.cc7c.covtab

#genecov output to sum=daphplx18cdschr_SRR3090572_genesum.txt means=daphplx18cdschr_SRR3090572.genemeans
#read sampledata=daphplx20chrs.metad
#read anntable=daphpulex_pa42v2.anntab
# read nid=33037, nclass=2 from daphplx17evgt1m_cds.idclass
#readChrtab(daphplx17evgt1m_cds_SRR3090572_1_bwa.cc7c.chrtab) nok=406382
#readCovtab(daphplx17evgt1m_cds_SRR3090572_1_bwa.cc7c.covtab,cds=1) nok=406382, nin=406382
#readCovtab(daphpulex_pa42v2_SRR3090572_1_bwa.cc7c.covtab,cds=0) nok=1895670, nin=1895670
#putGeneCovSum(daphplx18cdschr_SRR3090572_genesum.txt, daphplx18cdschr_SRR3090572.genemeans), asmid=daphpulex_pa42v2, KU=31

pt=daphplx_gasm16ml; 
env kucg=32 $evigene/scripts/genoasm/gnodes_genescov.pl -debug     -title daphplx16cdschr -asmid $pt \
 -metad daphplx20chrs.metad -idclass daphplx17evgt1m_cds.idclass  -ann $pt.anntab \
 -cdscov daphplx17evgt1m_cds_SRR3090572_1_bwa.cc7c.covtab  -chrcov ${pt}_SRR3090572_1_bwa.cc7c.covtab

* NOTE diff KU, should require same val for all asm here.. 
#genecov output to sum=daphplx16cdschr_SRR3090572_genesum.txt means=daphplx16cdschr_SRR3090572.genemeans
#read sampledata=daphplx20chrs.metad
#read anntable=daphplx_gasm16ml.anntab
# read nid=33037, nclass=2 from daphplx17evgt1m_cds.idclass
#readChrtab(daphplx17evgt1m_cds_SRR3090572_1_bwa.cc7c.chrtab) nok=406382
#readCovtab(daphplx17evgt1m_cds_SRR3090572_1_bwa.cc7c.covtab,cds=1) nok=406382, nin=406382
#readCovtab(daphplx_gasm16ml_SRR3090572_1_bwa.cc7c.covtab,cds=0) nok=1565033, nin=1565033
#putGeneCovSum(daphplx16cdschr_SRR3090572_genesum.txt, daphplx16cdschr_SRR3090572.genemeans), asmid=daphplx_gasm16ml, KU=35
dgmm:dplx20gnodes:% head -33 daphplx16cdschr_SRR3090572_genesum.txt
SUMMARY of daphplx16cdschr_SRR3090572.genemeans, asmid=daphplx_gasm16ml, KU=35

pt=daphplx20maca1b_findeg;
env kucg=32 $evigene/scripts/genoasm/gnodes_genescov.pl -debug     -title daphplx20ma1cdschr -asmid $pt \
 -metad daphplx20chrs.metad -idclass daphplx17evgt1m_cds.idclass  -ann $pt.anntab \
 -cdscov daphplx17evgt1m_cds_SRR3090572_1_bwa.cc7c.covtab  -chrcov ${pt}_SRR3090572_1_bwa.cc7c.covtab

#genecov output to sum=daphplx20ma1cdschr_SRR3090572_genesum.txt means=daphplx20ma1cdschr_SRR3090572.genemeans
#read sampledata=daphplx20chrs.metad
#read anntable=daphplx20maca1b_findeg.anntab
# read nid=33037, nclass=2 from daphplx17evgt1m_cds.idclass
#readChrtab(daphplx17evgt1m_cds_SRR3090572_1_bwa.cc7c.chrtab) nok=406382
#readCovtab(daphplx17evgt1m_cds_SRR3090572_1_bwa.cc7c.covtab,cds=1) nok=406382, nin=406382
#readCovtab(daphplx20maca1b_findeg_SRR3090572_1_bwa.cc7c.covtab,cds=0) nok=2837894, nin=2837894
#putGeneCovSum(daphplx20ma1cdschr_SRR3090572_genesum.txt, daphplx20ma1cdschr_SRR3090572.genemeans), asmid=daphplx20maca1b_findeg, KU=30

---------------

daphplx19cdschr_SRR3090572_genesum.txt
SUMMARY of daphplx18cdschr_SRR3090572.genemeans, asmid=daphpulex_pa42v2, KU=31
Class_Counts__________
item    32981
short   24
casm>CU 554   *? not high
casm<CU 3870
cds>CU  2905
cds<CU  7418
ctot>CU 4899   * high
ctot<CU 819
dLENbn  21091   ?? why so many, aplant next high at 9000
dLENcds 20130   ?? ""

Median,Mean,N,StDev,Range________
difbin  0.37    0.46    32981   0.76    0,27
tchr    1.19    1.93    32981   17.1    0,675   * tchr for daph higher than insect,aplant
wchr    2.60    15.3    31995   49.8    0,714   *? wchr rather high, is this dup gene effect?
xcds    0.90    1.63    32981   27.3    0,2990
xchr    1       1.04    32981   2.39    0,89.1
ycds_cr 0.88    1.23    32981   3.98    0,542

>> redo w/ kucg=32 for all dplx asm

daphplx20ma1cdschr_SRR3090572_genesum.txt
SUMMARY of daphplx20ma1cdschr_SRR3090572.genemeans, asmid=daphplx20maca1b_findeg, KU=32
Class_Counts__________
item    32982
short   25
casm>CU 561    * lower using std KU
casm<CU 3935
cds>CU  2728
cds<CU  8102
ctot>CU 2992
ctot<CU 1587
dLENbn  21945
dLENcds 20983

Median,Mean,N,StDev,Range________
difbin  0.38    0.44    32982   0.65    0,9
tchr    1.09    1.24    32982   1.98    0,84.3
wchr    2.50    22.2    31995   331     0,20069
xcds    0.88    1.58    32982   26.5    0,2896
xchr    0.94    0.96    32982   1.03    0,15.4
ycds_cr 0.89    1.54    32982   20.6    0,2261

daphplx16cdschr_SRR3090572_genesum.txt
SUMMARY of daphplx16cdschr_SRR3090572.genemeans, asmid=daphplx_gasm16ml, KU=32
Class_Counts__________
item    32986
short   29
casm>CU 1544  * higher using std KU
casm<CU 1788
cds>CU  2728
cds<CU  8106
ctot>CU 4407
ctot<CU 617
dLENbn  18430
dLENcds 17468

Median,Mean,N,StDev,Range________
difbin  0.48    0.50    32986   0.70    0,9
tchr    1.19    1.93    32986   20.2    0,944
wchr    2       9.9     31995   29.8    0,434
xcds    0.88    1.58    32986   26.5    0,2896
xchr    1.09    1.29    32986   6.31    0,624
ycds_cr 0.81    0.95    32986   3.26    0,412

>>> WRONG WAY USING SELF-KU depends on chrasm
daphplx16cdschr_SRR3090572_genesum.txt
SUMMARY of daphplx16cdschr_SRR3090572.genemeans, asmid=daphplx_gasm16ml, KU=35
Class_Counts__________
item    32986
short   29
casm>CU 1080   ** 2x of later dplx19ml
casm<CU 2205
cds>CU  2328
cds<CU  8793
ctot>CU 3631   * lower than dplx19ml
ctot<CU 789
dLENbn  18430
dLENcds 17468

Median,Mean,N,StDev,Range________
difbin  0.48    0.50    32986   0.70    0,9
tchr    1.09    1.76    32986   18.5    0,863
wchr    2       9.9     31995   29.8    0,434
xcds    0.80    1.45    32986   24.2    0,2648
xchr    1       1.18    32986   5.77    0,570
ycds_cr 0.81    0.95    32986   3.26    0,412

daphplx20ma1cdschr_SRR3090572_genesum.txt
SUMMARY of daphplx20ma1cdschr_SRR3090572.genemeans, asmid=daphplx20maca1b_findeg, KU=30
Class_Counts__________
item    32982
short   25
casm>CU 906    * higher than 19ml
casm<CU 2868
cds>CU  3111
cds<CU  7419
ctot>CU 3658  
ctot<CU 1154
dLENbn  21945
dLENcds 20983

Median,Mean,N,StDev,Range________
difbin  0.38    0.44    32982   0.65    0,9
tchr    1.17    1.32    32982   2.11    0,90
wchr    2.50    22.2    31995   331     0,20069
xcds    0.93    1.69    32982   28.2    0,3089
xchr    1       1.02    32982   1.10    0,16.4
ycds_cr 0.89    1.54    32982   20.6    0,2260

#-----------------------------------

=item try1c genecovsum 
 .. change SEM to SD 
 ** FIXME: diff ku/kun here are reversed from sense, kun means cov > KU, ku+ means cov < KU
 
dropse20cdschr_genesum.txt
SUMMARY of dropse20cdschr_SRR11813283.genemeans, asmid=dropse20chrs, KU=95
Class_Counts__________
item      14330
short     0
diffku    43
diffkun   21  ** FIXME wrong sense, see above
diffkcds  225
diffkcdsn 320
diffbn    2402
difflen   0 # missing chrtab len

Median,Mean,N,SEM,Range________
difbin  0.8421  0.78    14330   0.01    0.0,2.0000
xcds    0.9263  1.15    14330   0.05    0.0,304
xchr    1.0000  1.03    14330   0.02    0.0,143
ycds_cr 0.9406  1.06    14330   0.02    0.0,140

drosmel6cdschr_SRR11460805_genesum.txt
SUMMARY of drosmel6cdschr_SRR11460805.genemeans, asmid=drosmel6ref_chr, KU=14
Class_Counts__________
item      13946
short     26
diffku    225  # << chr missing gene dups, higher than dpse or arath ? is diffku count bad? 140 by means sort
diffkun   186  ** FIXME wrong sense, see above
diffkcds  884  # << 4x of dropse, ie dups missing from gene set?
diffkcdsn 184
diffbn    819
difflen   775

Median,Mean,N,SEM,Range________
difbin  0.9231  0.90    13946   0.01    0.0,4.
wchr    0.9231  0.90    13902   0.01    0.0,4.
xcds    0.9286  0.95    13946   0.01    0.0,23
xchr    1.0000  1.04    13946   0.01    0.0,8
ycds_cr 0.9286  0.91    13946   0.01    0,10.5

arath18cdschr_SRR10178325_genesum.txt
SUMMARY of arath18cdschr_SRR10178325.genemeans, asmid=arath18tair_chr, KU=18
Class_Counts__________
item      27406
short     25
diffku    130    # chrasm missing gene dups
diffkun   161    # chrasm has too many dups ? or bad gene models?
diffkcds  2223   # genes w/ dups ** FIXME wrong sense, see above
diffkcdsn 459    # below dna cov level means what? mistaken gene model? or excess dupl?
diffbn    11129  # << what do these large nums for cds length mismatches
difflen   11028  # << mean ??

Median,Mean,N,SEM,Range________
difbin  0.7429  0.70    27406   0.00    0.0,4.0
wchr    0.7500  0.70    27305   0.00    0.0,4.0
xcds    0.8889  1.05    27406   0.05    0.0,758
xchr    1.0000  1.09    27406   0.03    0.0,395
ycds_cr 0.8889  0.91    27406   0.01    0.0,154

=cut

=item try1

==> aweed20gnodes/arath18cdschr_SRR10178325.genemeans <==
GeneID	CD_med	CD_ave	CD_nit	CD_sem	CD_rng	Ch_med	Ch_ave	Ch_nit	Ch_sem	Ch_rng	CDt_med	CDt_ave	Cht_med	Cht_ave
AT1G67120t1	16	16.38	162	1.34	3,30	19	18.74	204	1.35	7,33	16	16.38	19	18.74
AT3G02260t1	18	17.69	153	1.48	3,32	18	18.24	160	1.49	6,30	18	17.69	18	18.24
AT5G23110t1	17	17.33	142	1.50	0,28	18	17.95	147	1.52	9,29	17	17.33	18	17.95
AT4G17140t1	16	15.62	127	1.46	3,30	18	17.94	156	1.48	8,34	16	15.62	18	17.94
AT1G48090t1	16	15.67	126	1.47	0,26	18	18.14	161	1.47	4,28	16	15.67	18	18.14
AT2G17930t1	17	17.77	116	1.74	0,42	19	19.27	269	1.21	8,32	17	17.85	19	19.27
AT4G36080t1	18	17.26	115	1.68	0,30	13	15.80	5	7.69	9,26	18	17.32	13	15.80
AT1G55860t1	16	16.95	111	1.68	2,30	18	18.13	231	1.24	10,40	19	18.99	19	19.72
AT1G70320t1	18	17.91	110	1.77	8,34	0	0	0	0	0	19	19.87	0	0

==> dropse20gnodes/dropse20cdschr_SRR11813283.genemeans <==
GeneID	CD_med	CD_ave	CD_nit	CD_sem	CD_rng	Ch_med	Ch_ave	Ch_nit	Ch_sem	Ch_rng	CDt_med	CDt_ave	Cht_med	Cht_ave
dropse20uc:g6903256t1	89	88.42	566	3.75	28,121	92	91.30	621	3.69	62,120	89	88.49	92	91.37
dropse20uc:g6900405t1	98	97.69	458	4.62	52,198	97	97.77	478	4.52	61,183	143	216.65	121	209.60
dropse20uc:g6902868t1	92	91.26	447	4.36	10,134	92	92.72	475	4.29	44,133	92	91.26	92	92.72
dropse20uc:g6899878t1	105	102.61	273	6.27	39,137	108	108.35	314	6.15	78,145	105	102.61	108	108.35
dropse20uc:g6898462t1	91	90.94	267	5.64	0,127	95	95.44	299	5.57	66,133	91	90.94	95	95.44
dropse20uc:g4815526t1	98	95.54	183	7.14	50,137	96	95.70	188	7.04	65,135	100	115.13	101	113.80
dropse20uc:g4818122t1	91	90.64	176	6.91	31,132	93	93.47	202	6.62	67,127	91	90.69	93	93.47
dropse20uc:g6903561t1	90	88.75	173	6.82	48,129	90	89.66	182	6.70	54,125	90	88.75	90	89.66
dropse20uc:g4816080t1	88	87.54	167	6.83	61,111	92	91.02	184	6.75	62,118	88	87.54	92	91.02

wc -l genemeans
   27383 aweed20gnodes/arath18cdschr_SRR10178325.genemeans
   14331 dropse20gnodes/dropse20cdschr_SRR11813283.genemeans

=cut

