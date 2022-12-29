#!/usr/bin/perl
# gnodes2_annotate.pl for evigene/scripts/genoasm/ 

=item usage gnodes2_annotate

  gnodes2_annotate  -chr chr.fasta -cds cds.fasta -te te.fasta 
    opts: -idclasses cdste.idclass : table of class per cds,te ID: CDS,TE,UNK ?RPT (simple repeat) 
            and modifiers: uniq,busco,duplicate, ..
  output: chr.anntab, same chr-locbin rows as sam2covtab for merging
  
=cut

use FindBin;
use lib ("$FindBin::Bin/../"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
our $EVIGENES="$FindBin::Bin/..";  
use strict;
use Getopt::Long;  

use constant UPD21JUN => 1;

use constant { # gnodes_annotate.pl current annot flags : export?
    A_CDS => 'CDS', A_TOP => 'Gtop', A_PART=>'partial', A_LOWQUAL=>'lowqual',
    A_BUSCO => 'busco', A_UCG => 'UCG',
    A_GAP => 'gap',
    A_TE => 'TE', # gff.type=transposon
    A_REPEAT => 'repeat', # gff.type=Simple_repeat or repeat
    A_UNK => 'UNK', # unknown/unclassified type (from repmask only?)
}; 

my $BN= 100; # bin size must match gnodes1_sam2covtab
my $MIN_IDENT= 0.65; # lo is best?
my $MIN_DUPIDENT = 0.98; # was .99/1; hi is best? or 1.0; # == ident equal to 1st/top align, lower if desired
my $evg_chr2agp="$EVIGENES/genoasm/chrfasta2agp.pl";

my $evg_blast2gff="$EVIGENES/genes/blast2evgff.pl"; # TESTs say not right..
my $USEBL2EVGFF= $ENV{usebl2evgff}||0;  # UPD21may12, MABYE use blast2evgff for blast2loctab, drops many valid cds aligns

my $evg_blastscore="$EVIGENES/makeblastscore.pl"; # try this way, maybe good > default?
my $USEBL2CDSLOCS=$ENV{usebl2cdslocs}||1; # UPD21JUN make default, adds copynum.table; -nousebl2cdslocs to turn off

my $debug=$ENV{debug}||1;
my $REUSE=1; # opt out -noreuse
my $ncpu=1;
my($outtab,$chrasm,$teseq,$cdsseq,$crclassf,$buscotsv,$genecovtab)=("") x 19;
my @ingff;

#UPD21apr26: add pMINTOP option .. test high = 0.99 to drop 2ndary aligns
# pMINTOP == MIN_DUPIDENT, 
# MINIDALN not same as MIN_IDENT, need min-ident-span, aln * pident in bases
my $MINIDALN=$ENV{minidal}||49; # opt, shorter for TE than CDS ?
# my $pMINTOP=$ENV{pmintop}||0.50; #? rather low? UPD21apr26: test pmintop=0.99 


# maybe add opt to input busco analysis/full_table_XXXX.tsv for cdseq > idclass table
my $optok= GetOptions( 
  'output=s',\$outtab, 
  'chrasm|assembly|genome=s',\$chrasm,
  'teseq=s',\$teseq,
  'cdsseq=s',\$cdsseq,  
  'angff=s',\@ingff, 
  'crclassf|idclasses=s',\$crclassf,  
  'buscotsv=s',\$buscotsv,  
  'genecovtab=s',\$genecovtab,  
  'ncpu=i', \$ncpu,# 'icpu=i', \$icpu, 
  'minident=s',\$MIN_IDENT, 
  'mindupident=s', \$MIN_DUPIDENT,  
  # 'minalign=s',\$MINALN,  
  'binsize=i', \$BN,
  'REUSE!', \$REUSE, 
  'usebl2cdslocs!', \$USEBL2CDSLOCS,
  'debug!', \$debug, 
  );

# FIXME: allow @gff inputs w/wo cds/te seq
if(@ARGV and my @addf= grep(/\.gff/, @ARGV)){ push @ingff, @addf;}

my $optinok= ($chrasm and ($teseq or $cdsseq or @ingff));

die "usage: gnodes2_annotate  -chr chr.fasta -cds cds.fasta -te te.fasta [ -out chr.anntab ]  
  opts: -idclass cdste.idclass -minident=NNN -minalign=NNN -ncpu 4 ; uses NCBI 'blastn'
"  unless($optok and $optinok);

$ncpu ||= 1;
my $OVTINY= int($BN*0.25); # exclude tiny bin hits for annots
my(%anntype, %busco, %orpar, %cloc, %agpc, ); # orig globals from generdcov4tab
my($nidclass,$idclassh,$idclasslist)=(0, {}, []);

sub MAIN_stub {}

  unless($outtab){ ($outtab=$chrasm) =~ s/\.\w+$//; $outtab.=".anntab"; }
  
  unless($crclassf) {
     my $ncl=0; ($ncl,$crclassf)= make_idclass( $crclassf, $cdsseq, $teseq,$buscotsv); 
  }
  
  ($nidclass,$idclassh,$idclasslist)= ($crclassf) ? read_idclass($crclassf, 1) : (0); # cds,te class by id, $idclass->{id} = class
  
  ## UPD21jun04
  ## add new genecov table from gnodes_sam2genecov.pl
  ## with tCopy, presumed from mRNA/CDS x reads map, as "true" tCopy value
  ## Gene_ID Glen Nread Rdlen tCopy C.M Uniq Merr C.nz S.Dup,Unq,Nt,Equ,UCG
  my($ngenecov,$genecovh)= ($genecovtab) ? read_genescov($genecovtab) : (0); 
  
  my($anngaps,$gapgff)= ann_gaps($chrasm);
  my($anncds,$cdsgff) = ann_cdslocs($chrasm,$cdsseq);
  my($annte,$tegff)   = ann_telocs($chrasm,$teseq);
  unshift @ingff, $gapgff if($anngaps);
  unshift @ingff, $cdsgff if($anncds);
  unshift @ingff, $tegff if($annte);

  # put_anntab FIXME: allow @gff inputs w/wo cds/te seq .. put_anntab(outtab, @allgff); ?
  put_anntab($outtab, @ingff);
  #x put_anntab($outtab, $anngaps,$gapgff, $anncds,$cdsgff, $annte,$tegff);

#------------------------------------------------------------------------  

sub put_anntab {
  my($outtab, @angff )=@_;
  #x my($outtab, $anngaps,$gapgff, $anncds,$cdsgff, $annte,$tegff )=@_;
 
   
  # step1: read each gff into cloc tab, converting to BN bin spans; do all at one go w/ sorted gff?
  # -- see read_annotgff() use for this?
  my $ingff= join" ", grep{ -f $_ } @angff;
  return(0) unless($ingff);
  # my $ingff="";   
  # $ingff .= " $gapgff" if($anngaps);
  # $ingff .= " $cdsgff" if($anncds);
  # $ingff .= " $tegff"  if($annte);
  
  my ($nannot,$no)=(0,0);
  my(%cloc); # global to local here
  
  # UPD21may27: USEBL2CDSLOCS cds.binlocs is in anntab format, merge that replacing cdsgff ..

  sub addIdclassType {
    my($id,$ty,$an)=@_;
    my $idclass= ($nidclass) ? $idclassh->{$id} || "" : "";
    if(not $idclass and $an =~ m/(?:class|type)=(\w+)/) { $idclass=$1; }
    if($idclass) {  
      my $tya="";
      if($idclass =~ /busco/i){ $tya=A_BUSCO; } # "busco"
      elsif($idclass =~ /orlog|ortholog|uniq/i){ $tya="uni"; }
      elsif($idclass =~ /parlog|paralog|dup/i){ $tya="dup"; }
      $ty.=",$tya" unless($ty=~/$tya/);
      unless($tya or $idclass =~ m/$ty/) { $ty.=",$idclass"; } # other annot
      }
    return($ty);
  }
    
  my @binlocs= grep /\.binlocs/, @angff; # UPD21JUN
  if(@binlocs) {
    @angff= grep{ not m/\.binlocs/ } @angff; $ingff= join" ", grep{ -f $_ } @angff;
    my($ok,$inh);
    if(@binlocs>1) {
      $ok= open($inh,"sort -k1,1 -k2,2n -k3,3 @binlocs | ");
    } else {
      ($ok,$inh)= openRead($binlocs[0]);
    }
    if($ok){
    while(<$inh>){
      next if(/^\W/);
      chomp; my($cr,$ib,$ty,$id,@amore)=split"\t",$_;
      #? check binsize, BN == 1+ $ib - $lastib
      
      # FIXME: insert idclass check: busco/UCG other annots of idclass change type
      # id here may be list: ida,idb,.. do all? or first?
      my($id1)= split",",$id;
      ($ty)= addIdclassType($id1,$ty,"");
      
      $anntype{$ty}++; $nannot++;
      my $oa= $cloc{$cr}{$ib}{annot}||""; 
      my $od= $cloc{$cr}{$ib}{ids}||""; 
      unless($oa =~ m/$ty,/) { my $tyo = ($oa) ?  "$oa,$ty" : $ty; $cloc{$cr}{$ib}{annot} = $tyo; }
      unless($od =~ m/$id,/) { my $ids = ($od) ?  "$od,$id" : $id; $cloc{$cr}{$ib}{ids} = $ids; }
    } close($inh);
    }
  }
  # if not $ingff after binlocs ??
  
  open(GFF,"sort -k1,1 -k4,4n -k5,5nr -k3,3 $ingff | ");
  while(<GFF>) {
    next if(/^\W/);
    # my $gff=join("\t",$cr,"evg",$qtype,$cb,$ce,$idal,$co,".","ID=$td;align=$al;ident=$pi;");
    my($cr,$src,$ty,$cb,$ce,$idal,$co,$xx,$an)= split"\t",$_;
    
    # UPD21may12: genes.gff may have mRNA, CDS rows, pick CDS only, Parent=ID
    #o my($id)= $an =~ m/ID=([^;\s]+)/?$1:"";
   
    my($id)= $an =~ m/(?:ID|Parent)=([^;\s]+)/;
    my $gd=0; if($id=~s/_[GC](\d+)$//){ $gd=$1; }  
    
    $ty=~s/transposon/TE/; #? any other class type rewrites?
    $ty=~s/Simple_repeat/repeat/; #? Satellite also
    next if($ty eq 'mRNA'); # skip,  others?
    ## next unless($ty =~ /^(CDS|exon|TE|repeat)/); # |mRNA|CDS .. or what?
    
    ## gff type should have CDS,TE/transposon,gap,repeat, basic class
    ## replace several ID tests w/ idclass modifies gff.type
    ## or/and check $an =~ m/class=(\w+)/
    
    ($ty)= addIdclassType($id,$ty,$an);
    # my $idclass= ($nidclass) ? $idclassh->{$id} || "" : "";
    # if(not $idclass and $an =~ m/(?:class|type)=(\w+)/) { $idclass=$1; }
    # if($idclass) {  
    #   my $tya="";
    #   if($idclass =~ /busco/i){ $tya="busco"; }
    #   elsif($idclass =~ /orlog|ortholog|uniq/i){ $tya="uni"; }
    #   elsif($idclass =~ /parlog|paralog|dup/i){ $tya="dup"; }
    #   $ty.=",$tya" unless($ty=~/$tya/);
    #   unless($tya or $idclass =~ m/$ty/) { $ty.=",$idclass"; } # other annot
    #   }
    
    # old type info .. drop
    my $orpar=$orpar{$id}||""; 
    if($orpar and $ty eq 'CDS') { 
      $ty .= ($orpar =~/parlog/) ? "dup" : ($orpar=~/orlog/)?"uni" : ""; 
      }
      
    $anntype{$ty}++; $nannot++;
    # for my $d (@ids) { if(my $bu= $busco{$d}) { $ty .= ",busco"; last; } } # or busco$bu w/ id?
    my($bb,$be)=map{ int($_/$BN) } ($cb,$ce);  
    for(my $i=$bb; $i<=$be; $i++){ 
      my $ib= $BN*(1+$i); # NOTE 1+i gives *end* of cds span, ie 1..100<, 201..300<
      ## FIXME maybe: require cb,ce to span at least 33%? of ib-BN .. ib, ie exclude tiny bin hits
      next if($cb > $ib - $OVTINY or $ce < ($ib + $OVTINY - $BN));
      # set,append $cloc{$cr}{$ib}{annot,ids}
      my $oa= $cloc{$cr}{$ib}{annot}||""; 
      my $od= $cloc{$cr}{$ib}{ids}||""; 
      unless($oa =~ m/$ty,/) { my $tyo = ($oa) ?  "$oa,$ty" : $ty; $cloc{$cr}{$ib}{annot} = $tyo; }
      unless($od =~ m/$id,/) { my $ids = ($od) ?  "$od,$id" : $id; $cloc{$cr}{$ib}{ids} = $ids; }
    }
  } close(GFF);
  
  # from generdcov4tab.pl
  my @cr= sort keys %cloc; # %$clochash; # want cleanchrtagOUT($cr) for sort, but need IN version for cloc{cr}
  my $ncr=@cr;
  my $clochash= \%cloc;
  
  rename($outtab,"$outtab.old") if(-s $outtab);
  open(OUT,">$outtab");
  # FIXME: @kcols not here
  warn "# write nchr=$ncr\n" if($debug);
  print OUT "#".join("\t","scaf","loc","annot","ids")."\n";  # ,@kcols
  
  for my $cr (@cr){ 
    my @ib=sort{ $a <=> $b } keys %{$cloc{$cr}};  
    for my $ib (@ib) { 
      my @cv=(); # @cv= map{ $cloc{$cr}{$ib}{$_}||0 } @kcols; # dont have @kcols,  
      
      # my $annv= $cloc{$cr}{$ib}{annot};
      my $annv= ($nannot) ? addannot($cr,$ib,$clochash) : "na"; ## addannot uses global %anntype
      my $ids= $cloc{$cr}{$ib}{ids};
      unless($ids) { $ids="noid"; }
      elsif( $ids=~m/,/ ) {
         my %ids=(); map{ $ids{$_}=1 } split",",$ids; 
         $ids=join",",sort keys %ids;
      }
      
      next if($annv eq "na" and $ids eq "noid"); # skip empty rows?
      #?? ($cr)= cleanchrtagOUT($cr); # do before sort?
      print OUT join("\t",$cr,$ib,$annv,$ids)."\n"; $no++; # no ,@cv here, but want table chr,ib to match covtab? or not?
    } 
  } close(OUT);
  return($no);
}


sub ann_gaps {
  my($chrs)= @_;

  # $evg_chr2agp == my $evapp="$EVIGENE/scripts/genoasm/chrfasta2agp.pl";
  (my $cname=$chrs) =~ s/\.fa.*//; $cname=~s/.gz//; 
  my $gapgff= "$cname.gaps.gff";
  my $chragp ="$cname.agp";
  my $SRC="chr2agp";
  if( -s $gapgff and $REUSE) {
    return(1,$gapgff,$chragp);
  }
  
  # opts -contigs $chrctg if want contigs.fasta, -mingap 10 is default, change?
  #      -MINSCLEN 200 -MINCTLEN 50 : defaults
  my $cmd="$evg_chr2agp -mingap 10 -chr $chrs -agp $chragp";
  # want chragp turned into gaps.gff, could pipe thru here, or else keep agp, pull out gaps
  
  my $ok= runcmd($cmd);

  # dmag14bgi:Contig0	169491	174700	7	W	c4_dmag14bgi:Contig0	1	5210	+
  # dmag14bgi:Contig0	174701	177242	8	N	2542	scaffold	yes	paired-ends
  # dmag14bgi:Contig0	177243	184322	9	W	c5_dmag14bgi:Contig0	1	7080	+
  # dmag14bgi:Contig0	184323	185176	10	N	854	scaffold	yes	paired-ends
  my($ngap)= (0);
  if(-s $chragp and open(F,$chragp)) {
    open(OG,">",$gapgff);
    while(<F>) { # see below read_agp
      next unless(/^\w\S+\t\d/);
      my($cid,$cb,$ce,$ci,$FN,$sid,$sb,$se,$so)=split;  
      if($FN =~ m/^[NU]/) { # sid = gap len
         ++$ngap; 
         print OG join("\t",$cid,$SRC,A_GAP,$cb,$ce,$sid,".",".","it=$ci")."\n";
      }
    } close(OG); close(F);
  }
  return($ngap,$gapgff,$chragp);
}

=item subdir bug
genome/*.fa subdir bugs:

#run blastn  -perc_identity 80 -evalue 1e-9 -task dc-megablast -template_type coding -template_length 18 -outfmt 7 -db genome/dmag19skasm -num_threads 24  \
 -out genome/dmag19skasm-genome/dmag7fincds.dcblastn  -query genome/dmag7fincds.fa 
Command line argument error: Argument "out". File is not accessible:  `genome/dmag19skasm-genome/dmag7fincds.dcblastn'
cant read genome/dmag19skasm-genome/dmag7fincds.dcblastn

=cut

use constant NOMAKEDB => 1;

sub ann_cdslocs {
  my($chrs,$cds)= @_;
  my($nloc,$loctab)=(0,0);
  unless($cds){ return(0); }

  # my $MINID_CDS_old= 80; #  # need opt
  # my $blopt_cds_old=" -perc_identity $MINID_CDS -evalue 1e-9 -task dc-megablast -template_type coding -template_length 18 -outfmt 7";
  
  #UPD21AUG16: High copy cds genes in pig suscr4evm are dragging out blastn, other annots.
  #new blopt much reduces blastn runtime for these, also looks okay for other _cds, _ter uses
  # ** blastopts should be user/gene-dataset option(s) 
  # blopt_hixc="-task megablast -perc_identity 90 -evalue 1e-9 -outfmt 7";
  my $MINID_CDS= 90; #was 80; # need opt
  #DAMM -- bug: my $blopt_cds="-task megablast  -perc_identity $MINID_CDS -evalue 1e-9 --outfmt 7";
  my $blopt_cds="-task megablast  -perc_identity $MINID_CDS -evalue 1e-9 -outfmt 7";

  my $isgz= ($cds =~ m/\.gz/);
  (my $cdsn=$cds) =~ s/\.\w+$//;

  my($dbname)= make_chrdb($chrs, NOMAKEDB); # dbname only
  my $blout= basename($dbname)."-".basename($cdsn).".dcblastn";
  my $blcmd= "blastn $blopt_cds -db $dbname -num_threads $ncpu -out $blout ";
  if($isgz) { $blcmd="gunzip -c $cds | $blcmd"; } else { $blcmd.=" -query $cds "; }

  #UPD21jun: USEBL2CDSLOCS change loctab .gff to .binlocs here for 
  ($loctab= $blout) =~ s/\.\w+$//;
  $loctab .= ($USEBL2CDSLOCS) ? ".binlocs" : ".gff"; 
  #o: $loctab.=".gff"; #? .gff.gz ?
  if( -s $loctab and $REUSE) { return(1,$loctab,$blout); }

  unless( -f $blout) { $blout.=".gz" if(-f "$blout.gz"); }
  unless( -s $blout and $REUSE) {
    make_chrdb($chrs); # really make if needed
    my $ok= runcmd($blcmd);
  }
  
  # INSERT parse blout into more usable format
  ($nloc,$loctab)= blast2loctab('CDS',$blout,$loctab,$cds);
  
  return($nloc,$loctab,$blout); # $ok, 
}

=item  3c. cdslocs current method:

  cds=daphmag14t1cds
  chrs=dmag15nwb2fullasm
  -- can this be reused for te.seq? not as good as blastn -task rmblast below 
  -- change dc-megablast to blastn? test 1 shows small increase in aligns, skip?
  
  see latest usage:
     daphmag/gasm20set/dmag20reasm/geneval/run_dcblastcds2locs.sh

  asm=dmag15nwb2fullasm
  orig cds class
  env cds=dmag14t1cds.fa  genome=genome/$asm.fa.gz ncpu=4 maxme=8000  datad=`pwd` ./run_dcblastcds2locs.sh

  $nbin/blastn  -evalue 1e-9 -perc_identity 80 -task dc-megablast -template_type coding -template_length 18 -outfmt 7  \
   -db $chrs -query $cds.fa  -out $chrs-$cds.dcblast1n  -num_threads $ncpu
   
MINIDALN=99; pMINTOP=0.66; 
grep -v '^#' $chrs-$cds.dcblast1n | cut -f1-4,9-10 | \
env minidal=$MINIDALN perl -ne 'BEGIN{ $MINIDALN=$ENV{minidal}||99; } 
($td,$cr,$pi,$al,$cb,$ce)= split; $idal=int($al * $pi/100); 
next if($idal < $MINIDALN); $pi=~s/\..*//; 
$co="f"; if($cb>$ce) { ($cb,$ce)=($ce,$cb); $co="r"; } 
print join("\t",$cr,"CDS",$cb,$ce,$co,$td,$idal,"$al,$pi")."\n";' | \
  sort -k1,1 -k3,3n -k4,4nr | 
env pmintop=$pMINTOP perl -ne 'BEGIN{ $pMINTOP=$ENV{pmintop}||0.50; } 
next if(/^\W/); ($cr,$ty,$cb,$ce,$co,$id,$idal,$alpi)=split; 
if($cr eq $lcr and $cb<$lce and $ce>$lcb) {
   if($lid ne $id and $idal >= $pMINTOP*$topal) { push @x, $id; }
} else { putx() if(@x); $topal=$idal; @x=($cr,$ty,$cb,$ce,$co,$id); }
($lcr,$lcb,$lce,$lid)=($cr,$cb,$ce,$id); END{ putx(); } 
sub putx{ my @xl=splice(@x,0,5); my $ids=join",",@x; print join("\t",@xl,$ids)."\n"; }  ' \
  > $chrs-$cds.cdslocs


=cut

=item repeats -dust no blast

$nbin/blastn $blopt -task blastn -dust no  -db dropse20chrs -query dropse20chrs_repeats.fa | grep ' hits found'
    ^^ use -task blastn -dust no  for most simpel repeats
# 577 hits found
# 98 hits found
# 96 hits found
# 15 hits found
# 4 hits found
# 215 hits found
# 64 hits found
# 2 hits found
# 255 hits found
# 78 hits found

-task rmblastn -dust no
# 589 hits found
# 113 hits found
# 96 hits found
# 4 hits found
# 1 hits found
# 1 hits found
# 9 hits found
# 1 hits found
# 122 hits found
# 15 hits found

versus dust default: -task blastn
# 577 hits found
# 0 hits found
# 96 hits found
# 4 hits found
# 4 hits found
# 0 hits found
# 0 hits found
# 1 hits found
# 0 hits found
# 2 hits found

=item ann_telocs
  v1: expects teseq.fa from species/chrasm, runs blastn to match to chrasm, so-so result
  v2: add 2 opts: 
    a. input tes == teseq.gff already located TE/rept on chrasm, check valid
    b. no tes input, run repeatmasker chrasm, then need opts: -species for rm,
        .. or dont run rm here, but translate rm.out to teseq.gff ?
=cut

sub parse_repmaskout {
  my($inout)= @_;
  my $no=0;
  (my $tegff = $inout) =~ s/\.out.*//; $tegff.=".rmte.gff";
  return(1,$tegff) if(-f $tegff);
  # daphcari20chr.farm.gff << bad name
  my($ok,$inh)= openRead($inout); return(0) unless($ok);
  
  open(G,">$tegff");
  print G "##gff-version 3\n# from RepeatMasker $inout\n";
  while(<$inh>) {
    next unless(/^\s*\d/);
    my @v=split; my($sc,$pv,$pe,$pi,$r,$b,$e,$qleft,$or,$repna,$repcla,$bi,$ei,$lefi,$id)=@v; 
    my $or=($or eq "C") ? "-" : ($or=~/^[+-]/) ? $or : "."; 
    my($tb,$te)=($or eq "-")? ($lefi,$ei) : ($bi,$ei); 
    my($rclass,$rfam)=split "/",$repcla,2; 
    $rclass =~ s/\W/_/g; $rclass=~s/_$//; $rfam="/$rfam" if ($rfam); 
    my $typ=($rclass =~ /Unknown|ARTEFACT/)?"UNK"
          : ($rclass =~ m/DNA|LINE|LTR|RC|SINE/) ? "transposon" : "repeat"; 
    $id= sprintf "rmask%06d",$id;    
    print G join("\t",$r,$rclass,$typ,$b,$e,$sc,$or,".","ID=$id;Name=$rclass$rfam;Target=$repna $tb $te\n");
    $no++;
  }
  close(G); close($inh);
  return($no,$tegff);
}

=item rmask.cat to te.seq

466 27.40 0.48 0.00 NC_037638.1 4243 4450 (27749750) Gypsy-24_DEl-I#LTR/Gypsy 2601 2809 (3339) m_b1s001i1

262 27.59 0.00 0.00 NC_037638.1 20221 20275 (27733925) C MuDR-1_LHu#DNA/MULE-NOF (635) 3923 3869 m_b1s001i6
   (27733925) C MuDR-1_LHu#DNA/MULE-NOF 
              ^ dingbat  extra colm for strand complement
              
 env idp=apimel perl -ne 'if(/^\d+/){ putfa() if($id and $fa);  @v=split;
 $cor=0; if($v[8] eq "C") { ($cor)= splice(@v,8,1); }
 ($cr,$cb,$ce,$te)=@v[4,5,6,8]; ($td,$tc)=split"#",$te;
 $cl=($tc=~/Unknown/)?"UNK":($tc=~/repeat|Sattelite/)?"repeat"
   :($tc=~/LTR|LINE|SINE|DNA|RC/)?"TE":"unk"; $cr=~s/\.\d//;
 $id=join"_",$cr,$cb,$ce;  $fa=""; $infa=1; } 
 elsif($infa){ ($d)=@v=split; if($d=~/^$cr/) { $s=$v[2]; $s=~s/\W//g; $fa.=$s; } } 
 END{ putfa(); } BEGIN{ $IDP=$ENV{idp}||"nada"; $IDP.="_rm"; }
 sub putfa{ $w=length($fa); $fa=~s/(.{100})/$1\n/g;
 print ">$IDP$id type=$cl; len=$w; tecl=$tc; teid=$td;\n",$fa,"\n"; } ' \
  apismel19hav31asm.fa.align > apismel19tefam.fa

=cut

sub repmask_align2teseq {
  my($rmalignf)= @_;
  my $no=0;
  (my $tefasta = $rmalignf) =~ s/\.(align|cat).*//; $tegff.="tefam.fa";
  my $IDP="noname";
  my($ok,$inh)= openRead($rmalignf); return(0) unless($ok);  
  my $outh; open($outh,">$tefasta"); 
  
  sub putfa { 
    my($outh,$id,$fa,$tean)= @_; 
    my $w=length($fa);  $tean=~ s/len=0/len=$w/;
    $fa=~s/(.{100})/$1\n/g; $fa.="\n" unless(substr($fa,-1) eq "\n");
    print $outh ">$id $tean;\n",$fa; return 1;
  }
  
  my($id,$fa,$tean,$atchr,$infa)=("") x 9;
  while(<$inh>){
    if(/^\d+/){ 
      $no += putfa($outh,$id,$fa,$tean) if($id and $fa); 
      my @v=split;
      my $cor=0; if($v[8] eq "C") { ($cor)= splice(@v,8,1); }
      my($cr,$cb,$ce,$te)=@v[4,5,6,8]; 
      my($td,$tc)=split"#",$te;
      my $cl=($tc=~/Unknown/)?"UNK":($tc=~/repeat|Sattelite/)?"repeat"
            :($tc=~/LTR|LINE|SINE|DNA|RC/)?"TE":"unk"; $cr=~s/\.\d//;
      my $id=$IDP . join"_",$cr,$cb,$ce;  $fa=""; $infa=1; $atchr= $cr;
      $tean="type=$cl; len=0; tecl=$tc; teid=$td;";
    } elsif(/^MATRIX/) { 
      $no += putfa($outh,$id,$fa,$tean) if($id and $fa); 
      $infa=0; $id=$fa=$tean="";
    } elsif($infa) { # align lines:  chr pos aaacccggg--aaaccc---ttt pos2
      my @v=split; if($v[0] eq $atchr) { my $s=$v[2]; $s=~s/\W//g; $fa.=$s; } 
    } 
 
  } close($inh); close($outh);
  
  return($no, $tefasta);
}


sub ann_telocs {
  my($chrs,$tes)= @_;
  # $nbin/blastn -task rmblastn  -perc_identity 80 -evalue 1e-9 -outfmt 7 -num_threads 4  \
  # or -task blastn .. either/or not sure

  unless($tes){
    return(0);
  } elsif($tes =~ /\.gff/){
    return(1,$tes,""); # always 
    
  } elsif($tes =~ /\.out/ and $tes =~ /rm/){  # rmout/chrasm.out ..
    my($nte,$tegff)= parse_repmaskout($tes);
    my($ntefa,$tefasta)=(0,"");
    if($nte){
      (my $rmcat=$tes) =~ s/.out/.cat/; 
      $rmcat.=".gz" unless(-f $rmcat); # rmaskr does gzip?
      ($ntefa,$tefasta)= repmask_align2teseq($rmcat) if(-f $rmcat);
    }
    return($nte,$tegff,$tefasta); # FIXME diff 3rd param here
  }
  
  # my $MINID_TE= 80; # opt
  # my $blopt_ter="-task blastn -perc_identity $MINID_TE -evalue 1e-9 -dust no -outfmt 7";

  #UPD21AUG16: High copy cds genes in pig suscr4evm are dragging out blastn, other annots.
  #new blopt much reduces blastn runtime for these, also looks okay for other _cds, _ter uses
  # ** blastopts should be user/gene-dataset option(s) 
  # blopt_hixc="-task megablast -perc_identity 90 -evalue 1e-9 -outfmt 7";
  my $MINID_TE= 90; #was 80; # need opt
  #DAMM -- bug: my $blopt_ter="-task megablast  -perc_identity $MINID_TE -evalue 1e-9 --outfmt 7";
  my $blopt_ter="-task megablast  -perc_identity $MINID_TE -evalue 1e-9 -outfmt 7";

  my $isgz= ($tes =~ m/\.gz/);
  (my $cdsn=$tes) =~ s/\.\w+$//; 
  my ($dbname)= make_chrdb($chrs, NOMAKEDB);
  my $blout= basename($dbname)."-".basename($cdsn).".blastn";
  my $blcmd= "blastn $blopt_ter -db $dbname -num_threads $ncpu -out $blout ";
  if($isgz) { $blcmd="gunzip -c $tes | $blcmd"; } else { $blcmd.=" -query $tes "; }

  my($nloc,$loctab)=(0);
  ($loctab= $blout) =~ s/\.\w+$//; $loctab.=".gff";
  if( -s $loctab and $REUSE) { return(1,$loctab,$blout); }
  
  unless( -s $blout and $REUSE) {
     make_chrdb($chrs); # really make
     my $ok= runcmd($blcmd);
  }
  
  # INSERT parse blout into more usable format
  ($nloc,$loctab)= blast2loctab('TE',$blout,$loctab,$tes);
  
  return($nloc,$loctab,$blout); # $ok, 
}

=item 3t. te.seq method

  see also
  $dmagnew/dmag20reasm/tefind/run_repmask2lib.sh 
  $dmagnew/dmag20reasm/tefind/run_repmod1a.sh 

  test te with run_dcblastcds2locs, 
  env cds=dmag14t1cds.fa  genome=genome/$asm.fa.gz ncpu=4 maxme=8000  datad=`pwd` ./run_dcblastcds2locs.sh
  
  no use blastn -task rmblast instead, better align totals
  $nbin/blastn -task rmblastn  -perc_identity 80 -evalue 1e-9 -outfmt 7 -num_threads 4  \
     -db genome/bldb/$pt  -query dmag20sk4tefam.fa -out $pt-dmag20sk4tefam.i80.rmblastn 
  .. and reduce minidal, 99b too large for TE aligns, reduce CDS also?
  * FIXME: change CDS type to RM type: TER or UNR from  dmag20sk4tefam ids: xxx_(te|un)fnum
  
  -- use CDS_TE.idclass table here also (as for sam2covtab)
  
pt=dmag20sk4maca20ok; 
grep -v '^#' $pt-dmag20sk4tefam.i80.rmblastn  | cut -f1-4,9-10 | env minidal=19 perl -ne \
'BEGIN{ $MINIDALN=$ENV{minidal}||99; } next if(/^\W/);
($td,$cr,$pi,$al,$cb,$ce)= split; $idal=int($al * $pi/100); 
next if($idal < $MINIDALN); $pi=~s/\..*//; 
$tty= ($td =~ m/_(..r)/) ? uc($1) : "UNR"; # TER/UNR
$co="f"; if($cb>$ce) { ($cb,$ce)=($ce,$cb); $co="r"; } 
print join("\t",$cr,$tty,$cb,$ce,$co,$td,$idal,"$al,$pi")."\n";' | \
sort -k1,1 -k3,3n -k4,4nr |  env pmintop=$pMINTOP perl -ne \
'BEGIN{ $pMINTOP=$ENV{pmintop}||0.50; } 
next if(/^\W/); ($cr,$ty,$cb,$ce,$co,$id,$idal,$alpi)=split; 
if($cr eq $lcr and $cb<$lce and $ce>$lcb) {
   if($lid ne $id and $idal >= $pMINTOP*$topal) { push @x, $id; }
} else { putx() if(@x); $topal=$idal; ($al)=split",",$alpi;  @x=($cr,$ty,$cb,$ce,$co,$al,$id); }
($lcr,$lcb,$lce,$lid)=($cr,$cb,$ce,$id); END{ putx(); } 
sub putx{ my @xl=splice(@x,0,5); ($xa)=shift @x; my $ids=join",",@x; print join("\t",@xl,$ids,$xa)."\n"; }  ' \
  > $pt-dmag20sk4tefam.rmtelocs

grep -v '^#' $pt-dmag20sk4tefam.i80.rmblastn  | cut -f1-4,9-10 | head
dmag20skm_unr1f179	dmag20sk4maca_ctg478356	94.393	428	820	394
dmag20skm_unr1f179	dmag20sk4maca_ctg478356	93.519	216	215	1
dmag20skm_unr1f179	dmag20sk4maca_scf238478	94.406	429	555	128
..
dmag20skm_ter1f10	dmag20sk4maca_scf194253	84.356	652	832	190
dmag20skm_ter1f10	dmag20sk4maca_scf194253	86.758	438	1113	676
dmag20skm_ter1f10	dmag20sk4maca_scf194253	85.714	175	1333	1162

-- output this format?
dmag15nwb2fullasm-dmag20sk4tefam.rmtelocs    
dmag24nwb7c_contig00573	UNR	243	280	r	dmag20skm_unr1f189	38
dmag24nwb7c_contig00852	TER	1	388	f	dmag20skm_ter6f2993	388
dmag24nwb7c_contig01062	TER	460	569	r	dmag20skm_ter1f704	110
dmag24nwb7c_contig01062	UNR	825	899	r	dmag20skm_unr1f790	75

=cut

=item UPD21may blast2evgenegff

  USEBL2EVGFF : blast2evgff limits to gene-complete (or most complete) align records
  .. not wanted here, want larger collection of most CDS aligns, including extra exons, pseudos

=cut

sub blast2evgenegff {
  my($btype,$bltab,$loctab)= @_;
  
  #UPD21may12: for gene.cds use blast2evgff instead of simple below 
  # $USEBL2EVGFF= $ENV{usebl2evgff}||0;  # MABYE use, drops many valid cds aligns
  # gunzip -c $pt-daphplx17evgt1m_cds.dcblastn.gz |  \
  # $evigene/scripts/genes/blast2evgff.pl -lowscore 0.80 -alignmax 999 -noaddcds -exont CDS -source evg \
  #  -out $pt-daphplx17bl2evgcdsx.gff

  my $cmd="$evg_blast2gff -lowscore $MIN_DUPIDENT -alignmax 1999 -noaddcds -exontype CDS -source evg  -out $loctab";
  if($bltab =~ m/.gz/){ $cmd="gunzip -c $bltab | $cmd"; } else { $cmd.=" $bltab"; }
  my $ok= runcmd($cmd);
  my $nloc=`grep -c CDS $loctab`; chomp($nloc);
  
  return($nloc, $loctab);
}

=item UPD21may27 blast2bestcdslocs

  USEBL2CDSLOCS: replace simple CDS blast parse w/ makeblastscore3.pl, to exon.gff, 
  tuned to keep p=0.50? low dups, but drop low dups at same locus
  my $USEBL2CDSLOCS=$ENV{usebl2cdslocs}||0;

  pt=daphplx_gasm16ml; ...
  $evigene/scripts/makeblastscore3.pl  -sizes $pt.facount,daphplx17evgt1m_cds.fa.qual  -spans=2 -pMINLOW=0.50 \
    $pt-daphplx17evgt1m_cds.dcblastn.gz -out  $pt-daphplx17evgt1m_cds.btall5x
  
  # btall 2 anntab : change to gff?
 env minpa=33 nolow=1 mindup=0.90 perl -ne \
 'BEGIN{ $MINDUP=$ENV{mindup}||0.90; $NOLO=$ENV{nolow}||0; $MINPA=$ENV{minpa}||33; $BN=100; } 
 ($td,$sc,$bs,$ida,$al,$tw,$sw,$tsp,$ssp)=@v=split; $pal=($tw<1)? 0 : 100*$ida/$tw; next if($pal<$MINPA); 
 ($sor)= ($ssp=~s/:(.)$//)?$1:0; @ssp=split",",$ssp; for $xn (@ssp){ 
  ($xb,$xe)=split"-",$xn; ($xb,$xe)=($xe,$xb) if($xb>$xe); ($ib,$ie)=map{ int($_/$BN) } ($xb,$xe); 
  for($i=$ib; $i<=$ie; $i++){ $ann{$sc}{$i}{$td}=$pal; } $max{$sc}=$ie if($ie>$max{$sc}); 
  }   
  END{ print "#".join("\t",qw(scaf loc annot ids))."\n"; for $sc (sort keys %max) { $imax=$max{$sc}; 
  for ($i=0; $i<=$imax; $i++){ next unless($ann{$sc}{$i}); 
  @ids=sort{ $ann{$sc}{$i}{$b}<=>$ann{$sc}{$i}{$a} or $a cmp $b } keys %{$ann{$sc}{$i}}; $pt=$ann{$sc}{$i}{$ids[0]}; 
  $ids=join",", grep/\w/, map{ my $p=$ann{$sc}{$i}{$_}; my $pf=($p/$pt<$MINDUP)?"-":""; ($pf and $NOLO)? "" : $_.$pf; } @ids; 
  $ib=$BN*($i+1); print join("\t",$sc,$ib,"CDS",$ids)."\n"; } }  } ' \
    $pt-daphplx17evgt1m_cds.btall5x > ${pt}_cds.bt5x33anntab 
   
=cut


sub _max{ return($_[0] < $_[1])?$_[1]:$_[0]; }
sub _min{ return($_[0] > $_[1])?$_[1]:$_[0]; }
sub _minnot0{ return ($_[0] == 0) ? $_[1] : ($_[1] == 0) ? $_[0] : _min(@_); }

sub blast2bestcdslocs {
  my($btype,$bltab,$loctab,$cds)= @_;
  
  # UPD21JUN
  my $MINPA=$ENV{minpa}||33; #? drop to 15 has small effect
  my $MINDUP=$ENV{pmindup}||0.90; # my $NOLO=1; # my $BN=100; 
  my $pMINLOW=$ENV{pminscore} || 0.90; #? same as MINDUP or not
  my $MINFULLCDS= $ENV{pfullcds}||85;
  
  #above: use constant { A_CDS => 'CDS', A_TOP => 'Gtop', A_PART=>'partial', A_LOWQUAL=>'lowqual' }; # annot flags
  
  #? need sizes, for cds at least ** $tw in btall table
  my $sizeopt=""; # my $cdsqual="FIXME.cds.qual";
  if( -f $cds) {
    my $cdsqual=$cds; $cdsqual=~s/\.gz//; $cdsqual.=".qual"; 
    unless( -f $cdsqual) {
      my $cmd="env ismrna=1 off=1 $EVIGENES/prot/aaqual.sh $cds";
      my $ok= runcmd($cmd);
    }
    $sizeopt="-sizes $cdsqual" if(-f $cdsqual);
  }
  
  my $nloc=0;
  $loctab =~ s/.\w+$//; # .gff 
  my $btall=  $loctab . ".btall";
  my $copytab=$loctab . ".genecopyn"; # used .genexcopy, keep that?
  $loctab.=".binlocs"; #<< merge w/ gff to output anntab
  
  # FIXME: ** update opts for makeblastscore.pl -tandupok to get tandem dupl locations **
  # and need to parse col 11 for tandem dup locs on same scaf
  # ?? use $pMINLOW for makeblastscore also ? No, want lowqual/partial aligns here, see below filt
  my $cmd="$evg_blastscore -pMINLOW 0.50 -tandupok -spans=2 $sizeopt -out $btall $bltab";
  
  if( -s $loctab and $REUSE) { 
    $nloc=`grep -c CDS $loctab`; chomp($nloc);
    warn "# reusing $loctab\n# from $cmd\n" if($debug);
    return($nloc, $loctab);
  } elsif( -s $btall and $REUSE) { 
    warn "# reusing $btall\n# from $cmd\n" if($debug);
  } else {
    my $ok= runcmd($cmd);
  }

=item UPD21Jun04 copynum table

  UPD21Jun04: add CopyNum table using this btall parse: both cross-scaf & same-scaf (tandupok)

  copynum++ if(td eq lasttd) and not( (pal < MINPA) and (bits < bmaxcut))
  * Need some test/adjustment on decision of what is a measurable gene copy,
  -- low qual, partial cds x chr aligns may, or not, equate to cds x dna map copy nums;
  -- $bmaxcut= $bs * $pMINLOW; has large effect on tCopy, Nc,s,t counts; pMINLOW=0.90 now .. maybe reduce
    
head daphpulex_pa42v2-cdsaln3g98d.genexcopy
GeneID____________      tCopy   Nc,s,t  Aln     Span    Glen
Daplx7pEVm000023t1      1.0     1,1,0   14018   14215   14268
Daplx7pEVm000024t1      1.0     1,1,0   12367   12336   12317
Daplx7pEVm000055t1      1.0     1,1,0   12711   12708   12708

  add xCopy using genecovtab if avail, xCopy = tCopy / genecov.tCopy
  my($ngenecov,$genecovh)= ($genecovtab) ? read_genescov($genecovtab) : (0); 
  cols: $tcopy,$cmed,$cnz,$tlen,$nread,$isuniq

  ** Add summary stats at end of copynum tab,
  per gCopy levels: 0= gC < 0.66, 1= 0.66<= gC <= 1.55, 2-9= 1.55 < gC <= 9.99, 10-99= 9.99<gC<=99, 100-499, 500+
gCn_level   nGene  xCopy(tCn/gCn)  tC=gC:-,=,+  ave(tCn,gCn)
       1
   2-  9
  10- 99
 100-499
    500+
    
=item eg Arath copynum tables
  
  -- note arath20max is ~10 Mb larger than arath18ta, and likely has some more gene copies
  .. but those stats may be obscured by uncertain filtering of what is a copy
  -- maybe bugs counting tandupl copies: tCopy from Align/Glen, but tandup addition missing?
     .. and tandup copy > 1 missing now from btall table (only 1 tandup tabulated) affects large chrs most      
       
==> arath18tair_chr-arath18tair1cds.genecopyn <==
AT2G05510t1     0       0,0,0   0       0       315     0       0.9
AT2G07623t1     0       0,0,0   0       0       246     0       6.4
AT5G03710t1     0       0,0,0   0       0       246     0       21.6
#---------------------------------------
# Genes Copynum Summary by Copy level, for nGene=27359 found, 3 missed
# tCopy=Chrasm, gCopy=Gene set, xCopy=tCopy/gCopy, x0,xLo,xEq,xHi= xCopy 0,<1,=1,>1
#CnGene nGene   tCopy   xCopy   gCopy   x0,xLo,xEq,xHi
#0      713     1.0     2.3     0.5     0,0,93,620
#1      26212   1.0     1.0     1.0     1,227,25828,156
#2-9    401     1.4     0.5     3.2     1,350,47,3
#10-99  23      1.4     0.1     22.7    1,22,0,0
#99-499 6       1.2     0.0     160.9   0,6,0,0
#500+   7       2.6     0.0     831.5   0,7,0,0
#---------------------------------------

==> arath20max_chr-arath18tair1cds.genecopyn <==
AT5G59616t1     0       0,0,0   0       0       462     0       1.1
AT5G59650t1     0       0,0,0   0       0       2655    0       0.9
AT5G61710t1     0       0,0,0   0       0       468     0       1.8
#---------------------------------------
# Genes Copynum Summary by Copy level, for nGene=27175 found, 176 missed
# tCopy=Chrasm, gCopy=Gene set, xCopy=tCopy/gCopy, x0,xLo,xEq,xHi= xCopy 0,<1,=1,>1
#CnGene nGene   tCopy   xCopy   gCopy   x0,xLo,xEq,xHi
#0      702     1.0     2.3     0.5     0,0,101,601
#1      26212   1.0     1.0     1.0     164,317,25556,175
#2-9    401     1.3     0.5     3.2     11,329,54,7
#10-99  23      1.5     0.1     22.7    1,22,0,0
#99-499 6       1.8     0.0     160.9   0,6,0,0
#500+   7       4.5     0.0     831.5   0,7,0,0
#---------------------------------------
    
=cut
  
  open(F,$btall) or return(0);

  use constant COPYNUM => 1;  # UPD21Jun04
  use constant COPYNUMzeros => 1;  # count chr-align missed genes in ave t/xCopy?
  use constant { C0hi=>0.66, C1hi=>1.55, C2hi=>9.99, C3hi=>99.9, C4hi=>499 };
  my @cn=(); our($ctabh,%cnsum,$ncnsum,%cndid,$ncmiss);

  if(COPYNUM) {
  open( $ctabh,">$copytab"); # save old?
  my @chd=qw(GeneID tCopy Nc,s,t Align Span Glen); my $xlab="";
  if($ngenecov>0) { push @chd, "xCopy","gCopy"; $xlab=" xCopy=tC/gC; gCopy=gene x DNA copy num"; }
  
  print $ctabh "#cols: tCopy = gene-align-spans / gene-length; Nc,s,t = Num copies (scaf,tandem); \n";
  print $ctabh "#cols:  Align, Span= sum align,spans; Glen = gene-length $xlab\n";
  print $ctabh "#".join("\t", @chd)."\n"; 
  }
  
  sub putcn { our($ctabh,%cnsum,$ncnsum,%cndid,$ncmiss);
    my($id,$tlen,$tcopy,$saln,$sspan,$ncp,$nscaf,$ntand)=@_;
    return 0 unless($saln > 0 and $ncp > 0);
    my @ecols=($saln,$sspan,$tlen);
    # if($ngenecov>0 and exists $genecovh->{$id}) 
    if($ngenecov>0) { 
      my $gcopy= (exists $genecovh->{$id}) ? $genecovh->{$id}->[0] : 0; 
      if($gcopy<0.01) {
        push @ecols, 0,0;
      } else {  
        my $xcopy= $tcopy/$gcopy;
        push @ecols, sprintf("%.1f",$xcopy), $gcopy; 
      
        my $gcl= ($gcopy < C0hi)?"0" :(C0hi<=$gcopy and $gcopy<=C1hi)?"1" 
          :(C1hi<$gcopy and $gcopy<=C2hi)?"2-9" :(C2hi<$gcopy and $gcopy<=C3hi)?"10-99"
          :(C3hi<$gcopy and $gcopy<=C4hi)?"100-499":"500+";
        my $ceq= ($xcopy<C0hi)? "xlo" : ($xcopy>C1hi)?"xhi": "xeq"; # add xzero?
        $cnsum{$gcl}{ngene}++; $cnsum{$gcl}{$ceq} ++; $ncnsum++;
        $cnsum{$gcl}{xcopy}+= $xcopy; $cnsum{$gcl}{tcopy}+= $tcopy; $cnsum{$gcl}{gcopy}+= $gcopy;
      }
    }
    $cndid{$id}=$ncp; # or tcopy   
    print $ctabh join("\t",$id,sprintf("%.1f",$tcopy),"$ncp,$nscaf,$ntand",@ecols)."\n";
  }

  sub putcn_missing { our($ctabh,%cnsum,$ncnsum,%cndid,$ncmiss);
    $ncmiss=0;
    if($ngenecov>0) { 
      for my $id (sort keys %{$genecovh}) {
        next if($cndid{$id});
        my($gcopy,$cmed,$cnz,$tlen,$nread,$isuniq)= @{$genecovh->{$id}};
        next if($gcopy < C0hi); # == zero in genecov also
        my @ecols=(0, 0, $tlen, 0, $gcopy); # ($saln,$sspan,$tlen, $xcopy, $gcopy);
        print $ctabh join("\t",$id, 0,"0,0,0",@ecols)."\n"; $ncmiss++;

        my $gcl= ($gcopy < C0hi)?"0" :(C0hi<=$gcopy and $gcopy<=C1hi)?"1" 
          :(C1hi<$gcopy and $gcopy<=C2hi)?"2-9" :(C2hi<$gcopy and $gcopy<=C3hi)?"10-99"
          :(C3hi<$gcopy and $gcopy<=C4hi)?"100-499":"500+";
        my $ceq= "xzero";  # add xzero here
        $cnsum{$gcl}{$ceq} ++; #NO: $ncnsum++;
        if(COPYNUMzeros) { $cnsum{$gcl}{ngene}++; $cnsum{$gcl}{gcopy}+= $gcopy; } #<< skip this sum for ave-copyn of nonzero gene set?
        # my $xcopy=0; my $tcopy=0;
        # zero: $cnsum{$gcl}{xcopy}+= $xcopy; $cnsum{$gcl}{tcopy}+= $tcopy; 
      }
    }
    return $ncmiss;    
  }  
  
  sub putcn_sum { our($ctabh,%cnsum,$ncnsum,%cndid,$ncmiss);
    return "Genes Copynum Summary=0" unless($ncnsum>0);
    my @gcl=sort{ $a<=>$b } keys %cnsum;  
    my $ngt=", for nGene=$ncnsum"; $ngt .= " found, $ncmiss missed" if($ncmiss>0);
    my @ceql= qw(xlo xeq xhi); unshift @ceql, "xzero" if($ncmiss>0);
    my $labxloeqhi= ($ncmiss>0) ? "x0,xLo,xEq,xHi= xCopy 0,<1,=1,>1"
      : "xLo,xEq,xHi= xCopy<1,=1,>1";
    (my $labxleh= $labxloeqhi) =~ s/=.*//;
    
    my $cnbrief="Genes Copynum Summary$ngt, xCopy= ";
    print $ctabh "#---------------------------------------\n";
    print $ctabh "# Genes Copynum Summary by Copy level$ngt\n";
    print $ctabh "# tCopy=Chrasm, gCopy=Gene set, xCopy=tCopy/gCopy, $labxloeqhi\n";
    print $ctabh "#".join("\t",qw(CnGene nGene tCopy xCopy gCopy), $labxleh)."\n";
    for my $gcl (@gcl) {
      my $ng= $cnsum{$gcl}{ngene} or next;  # ng includes xzero for COPYNUMzeros, x/t reduced by zeros
      my($xc,$tc,$gc)= map{ sprintf"%.1f",$cnsum{$gcl}{$_}/$ng } qw( xcopy tcopy gcopy) ;
      my($ceqn)= join",", map{ $cnsum{$gcl}{$_}||0 } @ceql; # qw(xlo xeq xhi); #? add xzero
      $gcl="99-499" if($gcl eq "100-499");
      print $ctabh "#".join("\t",$gcl,$ng,$tc,$xc,$gc,$ceqn)."\n";
      $cnbrief .= $gcl."c:$xc/$ng, " unless($gcl eq "0" or $gcl > 100);
    }
    #daphpulex_pa42v2_SRR13333791_b8_mim cnbrief=  
    #Genes Copynum Summary, for nGene=32408 found, 185 missed, xCopy= 1c:1.2/21151, 2-9c:0.8/3711, 10-99c:0.3/142, 99-499c:0.0/18
    print $ctabh "#---------------------------------------\n";
    return($cnbrief);
  }  
  
  
  my($lasttd,$bmax, $bmaxcut, %ann,%anf,%max);
  while(<F>){ next if(/^\W/);
    #o: my($td,$sc,$bs,$ida,$al,$tw,$sw,$tsp,$ssp)=split; 
    my($td,$sc,$bs,$ida,$al,$tw,$sw,$tsp,$ssp,$tdup,$sdup)=split; 
    
    my $nextgn= ($td ne $lasttd)?1:0;
    if(COPYNUM and $nextgn) {
      putcn(@cn);
      @cn=($td, $tw, (0) x 6); # id, Glen, tCopy, Aln, Span, Ncp, Nscaf, Ntand
    }
    $lasttd= $td;
    my $pal= ($tw<1) ? $ida : 100*$ida/$tw; 
    next if($pal < $MINPA); 
    
    ##? add other lowqual filters here? cuts many partial and some full 2ndary aligns
    ## anntab3e:640041 CDS w/o bmaxcut, anntab3f:557245 w/ bmaxcut 0.90
    ## sort of want most low qual CDS aligns, as evidence of old dups, etc; lower pMINLOW to 0.50 likely wouldnt cut much
    ## .. maybe set aflag to partial?
    if($nextgn){ $bmax=$bs; $bmaxcut= $bs * $pMINLOW; } # dup aligns should be together, for 2ndary locs
    # next if($pal < $MINPA or $bs < $bmaxcut); 

    my $aflag='';
    if($nextgn){ $aflag= A_TOP; } # ? add annot flag for 1st copy == top score
    elsif($bs < $bmaxcut) { $aflag = A_LOWQUAL; } #? A_PART
    
    if(COPYNUM and ($nextgn or $bs >= $bmaxcut)) { #add to @cn
      my($tsb,$tse)= $tsp=~m/^(\d+).(\d+)/;
      my $tspan= ($tse>$tsb)? 1+$tse-$tsb : $al;
      #x $cn[2] += ($tw<1) ? 1 : $tspan/$tw; #tcopy << should be align/tw instead of tspan/tw, tspan is oversized measure
      my $tcopy=  ($tw<1) ? 1 : $al/$tw;
      $cn[2] += $tcopy; # * FIXME tandup needs tcopy addition      
      $cn[3] += $al; 
      $cn[4] += $tspan; 
      $cn[6] ++; $cn[5] ++; # new scanf + all copy
      #?? only one tandup slot?
      if($tdup and $tdup ne "0d") { 
        $cn[7] ++; $cn[5] ++; # tand,all copy
        my($ttsb,$ttse)= $tdup=~m/^(\d+).(\d+)/; 
        my $ttspan= ($ttse>$ttsb)? 1+$ttse-$ttsb : $al;
        $cn[2] +=  ($tw<1) ? $tcopy : _min($ttspan/$tw,$tcopy); # * FIXME tandup needs tcopy addition 
      }
    }
    
    my($sor)= ($ssp=~s/:(.)$//)?$1:0; 
    my @ssp=split",",$ssp; 
    if($sdup and $sdup ne "0d"){ 
      $sdup=~s/:(.)$//; my @sdup = map{ s/$sc://; $_ } split",",$sdup;  push @ssp, @sdup; 
      }
    for my $xn (@ssp){ 
      my($xb,$xe)=split"-",$xn; ($xb,$xe)=($xe,$xb) if($xb>$xe); 
      next unless($xb>0 and $xe>0); # ($xb =~ m/^\d/ and $xe =~ m/^\d/); # no scaf ids here? from @sdup?
      my($ib,$ie)=map{ int($_/$BN) } ($xb,$xe); # $max{$sc}=$ie if($ie>$max{$sc}); 
      for(my $i=$ib; $i<=$ie; $i++){ $ann{$sc}{$i}{$td}=$pal; $anf{$sc}{$i}{$td}=$aflag if($aflag); } 
    }   
  } close(F);
  
  if(COPYNUM) { 
    putcn(@cn); 
    putcn_missing(); 
    my $cnsum= putcn_sum(); 
    close($ctabh); 
    warn "#$cnsum table to $copytab\n" if($debug);
  }
  
  # $loctab=~s/.gff//; $loctab.=".binlocs"; #<< merge w/ gff to output anntab
  open(my $outh,">$loctab");
  print $outh "#".join("\t",qw(scaf loc annot ids))."\n"; 
  for my $sc (sort keys %ann) { 
    # $imax=$max{$sc}; for ($i=0; $i<=$imax; $i++){ next unless($ann{$sc}{$i}); .. }
    my @ibins= sort{$a <=> $b} keys %{$ann{$sc}};
    for my $i (@ibins) {
      my @ids=sort{ $ann{$sc}{$i}{$b}<=>$ann{$sc}{$i}{$a} or $a cmp $b } keys %{$ann{$sc}{$i}}; 
      my $idtop=$ids[0];
      my $ptop= $ann{$sc}{$i}{$idtop}||1; 
      my $ids=join",", grep{ ($ann{$sc}{$i}{$_} / $ptop >= $MINDUP) } @ids; 
      my $ib= $BN*($i+1); 
      my $atype= A_CDS; 
      if(my $aflag= $anf{$sc}{$i}{$idtop}){ $atype.=",$aflag"; }
      $atype .= ','.A_PART if($ptop < $MINFULLCDS);
      print $outh join("\t",$sc,$ib,$atype,$ids)."\n";  $nloc++;
    } 
  } close($outh);
   
  return($nloc, $loctab);
}



sub blast2loctab {
  my($btype,$bltab,$loctab,$qseq)= @_;

  #UPD21apr26: see above, pMINTOP = MIN_DUPIDENT now  
  # my $MINIDALN=$ENV{minidal}||49; # opt, shorter for TE than CDS ?
  # my $pMINTOP=$ENV{pmintop}||0.50; #? rather low? UPD21apr26: test pmintop=0.99 

  unless($loctab) {
    ($loctab= $bltab) =~ s/\.\w+$//; $loctab.=".gff";
  }
  
  warn "#blast2loctab($btype,$bltab,$loctab)\n" if($debug);
  if($USEBL2CDSLOCS and $btype eq 'CDS') { 
    return blast2bestcdslocs($btype,$bltab,$loctab,$qseq); # UPD21may27, MABYE  want this, test
  }
  if($USEBL2EVGFF and $btype eq 'CDS') { 
    return blast2evgenegff($btype,$bltab,$loctab); # UPD21may12, MABYE, not probably, want this, test
  }
  
  my($ok,$inh)= openRead($bltab);
  unless($ok){ warn "cant read $bltab\n"; return 0; }
  my($nloc,$ltd,$lcr,$qid,$qtype,$qtypeb)=(0) x 9;
  my(%locs,%idscore);

=item FIXME TE-reptmodler ID syntax parse

blastn -query repmdlr_output.fa -db chrs ...
.. parse funky "ID#LTR|etc" TE classes

# Query: rnd-1_family-101#LTR/Gypsy ( RepeatScout Family Size = 139, Final Multiple Alignment Size = 100, Localized to 33 out of 35 contigs )
# Database: dropse20chrs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 23 hits found
rnd-1_family-101#LTR/Gypsy	dropse20uc:ChrX	82.519	2763	473	8	1	2756	32641405	32638646	0.0	2938

=cut
  
  while(<$inh>) {
    if(/^\W/){
      # parse Query? for things like teclass=Unknown|LINE|LTR|DNA and type=CDS
      if(/^# Query: (\S+)/){ $qid=$1; $qtype="UNK"; $qtypeb="";
      
        if( m/\stype=(\w+)/ ) { $qtype=$1; }
        elsif(m/class=(\w+)/) { $qtype=$1; }
        #dgg te id syntax: ($td =~ m/_(..r)/) ? uc($1) : "UNR"; # TER/UNR
        elsif($qid =~ m/_(..r)/) { $qtype= uc($1); }
        
        # FIXME btype eq 'TE', convert type/class to simpler set 
        #  transposon, repeat 
        if($btype eq 'TE') { # want this or not? keep detail as class=$classb
          if($qid =~ m/\#(\S+)$/) { $qtype= $1; }# repmodlr ID Query: rnd-1_family-101#LTR/Gypsy 
          if($qtype=~m/DNA|LINE|LTR|RC|SINE/){ $qtypeb=$qtype; $qtype='transposon'; }
          elsif($qtype=~m/repeat|Satellite/i){ $qtypeb=$qtype;  $qtype='repeat'; }# Simple_repeat, ? rRNA, snRNA ?
          elsif($qtype=~/Unknown/i){ $qtypeb=$qtype; $qtype='UNK'; }
          else { $qtypeb=$qtype; $qtype='UNK'; } # what is 'buffer' ??
        }

      }
      next;
    }
    
    my @v=split; my($td,$cr,$pi,$al,$tb,$te,$cb,$ce)= @v[0,1,2,3,6,7,8,9]; # blast -outfmt 6,7
    my $idal=int($al * $pi/100); 
    #UPD21apr26: change al*pi >= MINIDALN to pi >= MIN_IDENT, not so good but what?
    if($idal < $MINIDALN){
      next unless($td eq $ltd and $cr eq $lcr);    
    }
    my $xan=""; $xan .= "bclass=$qtypeb;" if($qtypeb); # extras  
    $pi=~s/\..*//;  my $co="+"; if($cb>$ce) { ($cb,$ce,$co)=($ce,$cb,"-"); } 
    #UPD3mar: add td align spans: "Target=$td $tb $te;" or simpler tspan=$tb-$te; ?
    $xan .= "Target=$td $tb $te;"; # or "span=$tb-$te;"
    my $gff=join("\t",$cr,"evg",$qtype,$cb,$ce,$idal,$co,".","ID=$td;align=$al;ident=$pi;$xan");
    $locs{$cr}{$cb} .= $gff . "\n"; #? add ce to sort largest span? should do contained-in loc filter
    $idscore{$td}= $idal if($idscore{$td} < $idal); #UPD3mar: idscore to filter lower qual aligns
    $nloc++;
    # print join("\t",$cr,"evg",$qtype,$cb,$ce,$idal,$co,".","ID=$td;align=$al;ident=$pi;")."\n"; # to gff format ?

    $ltd=$td; $lcr=$cr;
  } close($inh);
  
  open(OG,">$loctab");
  print OG "##gff-version 3\n";
  my($lcr,$lcb,$lce,$lid,$keepal)=(0) x 9;
  my %idok; # fixme use idok{id} instead of lid, for interspersed repeats
  $nloc=0;
  for my $cr (sort keys %locs) {
    my @cb= sort{$a <=> $b} keys %{$locs{$cr}};
    for my $cb (@cb) {
      my @gff= split"\n",$locs{$cr}{$cb};
      for my $gff (@gff) {
        my($crD,$src,$qtype,$cbD,$ce,$idal,$co,$xx,$an)= split"\t",$gff;
        my($id)= $an =~ m/ID=([^;\s]+)/;
        my $idscore= $idscore{$id}||0;
        if($cr eq $lcr and $cb<$lce and $ce>$lcb) { # simple overlap
          my $ok= ($cb >= $lcb and $ce <= $lce)?0:1; # skip contained-in 
          $ok= 1 if($ok and $idok{$id} or $idal >= $keepal);
          if($ok) { print OG $gff,"\n"; $nloc++; $idok{$id}++; }
          #x if($lid eq $id or $idal >= $keepal) { print OG $gff,"\n"; $nloc++; }
          #o: if($lid ne $id and $idal >= $keepal) { print OG $gff,"\n"; }
        } else { 
          print OG $gff,"\n"; $nloc++; $idok{$id}++;
          $keepal= $MIN_DUPIDENT * $idscore;
          #xb: $keepal= $MIN_DUPIDENT*$idal; # BAD after sort @cb loc, idal here no longer top score
        }
        ($lcr,$lcb,$lce,$lid)=($cr,$cb,$ce,$id); 
      }
    }
  } close(OG);
   
  return($nloc, $loctab);
}

sub make_chrdb {
  my($chrs, $nomake)= @_;
  my $isgz= ($chrs =~ m/\.gz/);
  my $dbname=$chrs; $dbname=~s/.gz//; $dbname=~s/.fa\w*$//;
  
  unless($nomake or -s "$dbname.nsq") {
    my $blcmd="makeblastdb -dbtype nucl -out $dbname -title $dbname";
    if($isgz){ $blcmd= "gunzip -c $chrs | ". $blcmd; } else { $blcmd.=" -in $chrs"; }
    my $ok= runcmd($blcmd);
    return($dbname); # $ok,
  }
  return($dbname);
}

sub runcmd {
  my($cmd)=@_;
  warn "#run $cmd\n" if($debug);
  my $ok= system($cmd);
  return ($ok);
}


#=============================================
#  parts from generdcov4tab.pl --- need updates
#=============================================

  
sub read_annotgff {
  my($ancla,$angff)=@_;
  unless($ancla){ $ancla= ($angff=~/(cdslocs|repmask|gff)/)?$1:"na"; } # use param
  my $nannot=0; # NOT global
  
  my($ok,$inh)= openRead($angff);
  unless($ok){ warn "cant read $ancla:$angff\n"; return 0; }
  
  if($ancla eq "gff"){ # or genegff ?
    # genes.gff, score exons instead of CDS, but rename exon => CDS
    # open(F,$angff); 
    while(<$inh>){ 
    next if(/^\W/); chomp; my @v=split"\t";
    my($cr,$src,$ty,$cb,$ce,$vsc,$co,$vph,$ann)=@v; # YES gff
    
    next if($cr eq 'NOPATH'); # gmap std unmapped
    # NOPATH	100	0	0	0	0	0	0	na	Dapma7bEVm020481t1
    
    ($cr)= cleanchrtagIN($cr);
    $ty =~ s/transposon/TE/;  
    next unless($ty =~ /^(exon|TE)/); # |mRNA|CDS .. or what?
    $ty='CDS' if($ty eq 'exon'); #? change exon to CDS for simplicity?
    
    my($id)= $ann =~ m/(?:ID|Parent)=([^;\s]+)/;
    my $gd=0; if($id=~s/_[GC](\d+)$//){ $gd=$1; } # my($gd)= ($id=~s/_[GC](\d+)$//) ? $1 : 0; 
    # my @ids=split",",$ids; my($id)=$ids[0]; #? sort
    my @ids=($id); my $ids=$id; 
    
    my $orpar=$orpar{$id}||""; 
    if($orpar and $ty eq 'CDS') { 
      $ty .= ($orpar =~/parlog/) ? "dup" : ($orpar=~/orlog/)?"uni" : ""; 
      }
    $anntype{$ty}++; $nannot++;
    for my $d (@ids) { if(my $bu= $busco{$d}) { $ty .= ",".A_BUSCO; last; } } # or busco$bu w/ id?
    my($bb,$be)=map{ int($_/$BN) } ($cb,$ce);  
    for(my $i=$bb; $i<=$be; $i++){ 
      my $ib= $BN*(1+$i); # NOTE 1+i gives *end* of cds span, ie 1..100<, 201..300<
      my $oa= $cloc{$cr}{$ib}{annot}||""; 
      my $od= $cloc{$cr}{$ib}{ids}||""; 
      # FIXME here? off-by-1 diff for CDS(dup,uni) and IDs, due to $ib vs $i calc, $i precedes $ib
      #o: unless($oa =~ m/$ty,/ or ($oa =~ m/CDS/i and $gd > 1)) # need to annotate gd<=1, others?
      unless($oa =~ m/$ty,/) # need to annotate gd<=1, others?
        { 
        if($oa){ $ty = "$oa,$ty"; } $cloc{$cr}{$ib}{annot} = $ty; 
        if($od){ $ids= ($od =~ m/$ids/) ? $od : "$od,$ids"; }
        $cloc{$cr}{$ib}{ids} = $ids;  #?? merge ids? may use both genes.gff, cdslocs
        }
      }
    } 
  
  } elsif($ancla eq "cdslocs"){
    # open(F,$angff); 
    while(<$inh>){ next if(/^\W/); my @v=split;
    my($cr,$ty,$cb,$ce,$co,$ids)=@v; # NOT gff
    ($cr)= cleanchrtagIN($cr);
    $ty=~s/transposon/TE/; # cdslocs type ty may be TE now .. from te gene annots

    # ** annot conflict: TEuni,busco ; n=284 **
    
    my @ids=split",",$ids; my($id)=(@ids)?$ids[0]:""; #? sort
    my($gd)= ($id=~s/_[GC](\d+)$//) ? $1 : 0; 

    # add busco.ids annnot * use at cdslocs w/ idlist

    my $orpar=$orpar{$id}||""; 
    if($orpar and $ty eq "CDS") { 
      $ty .= ($orpar =~/parlog/) ? "dup" : ($orpar=~/orlog/)?"uni" : ""; 
      }
    $anntype{$ty}++; $nannot++;  
    for my $d (@ids) { if(my $bu= $busco{$d}) { $ty .= ",".A_BUSCO; last; } } # or busco$bu w/ id?
    my($bb,$be)=map{ int($_/$BN) } ($cb,$ce);  
    for(my $i=$bb; $i<=$be; $i++){ 
      my $ib= $BN*(1+$i); # NOTE 1+i gives *end* of cds span, ie 1..100<, 201..300<
      ## drop annloc{} for  $cloc{$cr}{$ib}{annot}

      my $oa= $cloc{$cr}{$ib}{annot}||""; 
      my $od= $cloc{$cr}{$ib}{ids}||""; 
      unless($oa =~ m/$ty,/) # need to annotate gd<=1, others?
        { 
        if($oa){ $ty = "$oa,$ty"; } $cloc{$cr}{$ib}{annot} = $ty; 
        if($od){ $ids= ($od =~ m/$ids/) ? $od : "$od,$ids"; } # FIXME many dup ids
        $cloc{$cr}{$ib}{ids} = $ids;  #?? merge ids? may use both genes.gff, cdslocs
        }
        
      }
    }  
      
  } elsif($ancla eq "repmask") {
    # open(F,$angff); 
    while(<$inh>){ next if(/^\W/); 
      my @v=split; my($cr,$ty,$cb,$ce)=@v[0,2,3,4]; next unless($cb and $ce > $cb);
      ($cr)= cleanchrtagIN($cr);
      $ty=~s/transposon/TE/;
      $anntype{$ty}++; $nannot++; # use for addannot() return vector
      my($bb,$be)= map{ int($_/$BN) } ($cb,$ce); #? maybe dont use BIN here?
      #? any te/rep value or just presence?
      for(my $i=$bb; $i<=$be; $i++){ 
        my $ib= $BN*(1+$i); # NOTE 1+i gives *end* of cds span, ie 1..100<, 201..300<
        $ty=",$ty" if($cloc{$cr}{$ib}{annot}); $cloc{$cr}{$ib}{annot} .= $ty;
        } # many ty?
    } 
  }
  close($inh);
  my @ank= sort keys %anntype;
  warn "# read n=$nannot @ank from $angff\n" if($debug);
  return($nannot);
}

sub addannot {
  my($cr,$ib,$clochash)=@_;
  # NOTE: ib here is *end* of bin, ie. 1..100<, 301..400<
  # but annloc{cr}{bb} gives *start* of bin, ie int(99/100) = 0, int(399/100)= 300
  my $ANVAL=2; # any val > 0?
  my %val= map{ $_ => 0 } keys %anntype;
  if( my $anty= $clochash->{$cr}{$ib}{annot} ) {
    for my $ty (grep /\w/, split(",", $anty)) { $val{$ty} += $ANVAL; }
  }

  if(wantarray) { 
    my @vec= map{ $val{$_} } sort keys %val; 
    return @vec; #? only vals? or key=val hash ?
  } else {
    my $anv= join ",", grep { $val{$_} > 0 } sort keys %val; 
    return $anv || "na";
  }
}

=item make_idclass

  grep -h '>' dmag20sk4tefam.fa dmag14t1cds.fa  |  perl -ne \
  '($id)=m/>(\S+)/; $cl="UNK"; if(/type=CDS/){ $cl=(m/Class:Transposon|Name=TE:/)?"TE":"CDS"; } 
  elsif(/teclass=(\w+)/){ $tc=$1; $cl=($tc=~/DNA|RC|LINE|SINE|LTR/)?"TE":"UNK"; } print "$id\t$cl\n"; ' \
   > dmag10cdste.idclass
  # add busco   
  perl -ne '($d,$e)=@v=split; if(/^EOG/){ $bus{$e}=$d; } else { if($bus=$bus{$d}){ s/$/ busco/; } print; }' \
     ../../dmag7finloc_cds.busco.t1ids  dmag20cdste.idclass > dmag20cdste.idclassb
  # or ..
  grep -h '>' dmag20sk4tefam.fa dmag14t1cds.fa  | env busco=../../dmag7finloc_cds.busco.t1ids perl -ne \
  'BEGIN{ if($bf=$ENV{busco}){ open(F,$bf); while(<F>){ ($bd,$td)=split; $bus{$td}=$bd if($td); } close(F); } }
  ($id)=m/>(\S+)/; $cl="UNK"; if(/type=CDS/){ $cl=(m/Class:Transposon|Name=TE:/)?"TE":"CDS"; $cl=" busco" if($bus{$id}); } 
  elsif(/teclass=(\w+)/){ $tc=$1; $cl=($tc=~/DNA|RC|LINE|SINE|LTR/)?"TE":"UNK"; }  print "$id\t$cl\n"; ' \
   > dmag10cdste.idclass

=item add idclass busco table parse?
   .. maybe add opt to input busco analysis/full_table_XXXX.tsv for cdseq > idclass table
   .. look for cdsseq => full_table_*$CDSNAME*.tsv
   
=cut
  
sub readbuscotab { 
  my($butab,$idtab)=@_; 
  my $nbu=0;  my(%buval,%budval,%bd2id);    
  my($ok,$inh)= openRead($butab); return(0) unless($ok);
  while(<$inh>){ 
    # fixme : use Duplicate, pick highest score, for messy trasm w/ alts: option?
    next if(/^\W/);
    my($bd,$bc,$id,$score)=split;
    # FIXME need buval{bd} not id for 1/Dup
    if($bc =~ /Complete/){ $buval{$id}=$score; $budval{$bd}=$id; $bd2id{$bd}= $id; }
    elsif($bc =~ /Duplicate/) { if($score > $budval{$bd}){ 
      if(my $oid= $bd2id{$bd}) { delete $buval{$oid}; }
      $buval{$id}=$score; $budval{$bd}=$score; $bd2id{$bd}= $id; } 
    }
    # if(/^(\w+)\s+Complete\s+(\S+)/){ my($bd,$id)=($1,$2); if($idtab->{$id}){ $idtab->{$id}.=" busco"; $nbu++; } }
  } close($inh);
  for my $id (sort keys %buval) { $idtab->{$id}.=" busco" unless($idtab->{$id}=~/busco/); $nbu++; }
  warn "#readbuscotab: nbu=$nbu from $butab\n" if $debug; # got here twice, where from?
  return ($nbu); 
} 

  ## UPD21jun04
  ## add new genecov table from gnodes_sam2genecov.pl
  ## with tCopy, presumed from mRNA/CDS x reads map, as "true" tCopy value
  ## cols: Gene_ID Glen Nread Rdlen tCopy C.M Uniq Merr C.nz S.Dup,Unq,Nt,Equ,UCG

sub read_genescov {
  my($genecovtab)=@_;
  my($ngc,$ndup,$nuni,%gcovtab)= (0,0,0); 
  my($ok,$inh)= openRead($genecovtab); return(0) unless($ok);
  while(<$inh>){ 
    next if(/^\W/);
    ## cols: Gene_ID Glen Nread Rdlen tCopy C.M Uniq Merr C.nz S.Dup,Unq,Nt,Equ,UCG
    my($id,$tlen,$nread,$rdlen,$tcopy,$cmed,$isuniq,$merr,$cnz,$sdune)=split;
    next unless($tlen>0 and $nread>0);
    $gcovtab{$id}= [$tcopy,$cmed,$cnz,$tlen,$nread,$isuniq]; $ngc++; 
    if($tcopy>=1.75){ $ndup++; } elsif($tcopy>=0.50){ $nuni++; }
  } close($inh);
  warn "#read_genescov: ngene=$ngc ($nuni 1-copy, $ndup 2+copy) from $genecovtab\n" if $debug;  
  return ($ngc,\%gcovtab); 
}
  
sub make_idclass {
  my($outclass, @cdsteseqs)=@_;
  my($no,$ok,$nbu,$inh,$outh,%clt)= (0,0,0);
  unless($outclass) {
    $outclass=$cdsteseqs[0]; $outclass=~s/\.\w+$//; $outclass.=".idclass";
  }

  ## BAD >> overwrote input cds.fa !!!
  #make_idclass: nout=0 () to daphsim17evgt1cds.fa
  # read nid=0, nclass=1 from daphsim17evgt1cds.fa
  # write nchr=441
  
  my $incds=0; my %idtab; my $butab="";
  for my $inf (@cdsteseqs) {  
    $incds=0; my $nin=0;
    
    # teseq may be repmask out here? need out to gff first
    if($inf =~ m/\.out/ and $inf =~ m/rm/) {
      my($nte,$tegff)= parse_repmaskout($inf);
      push @cdsteseqs, $tegff if($nte); 
      
    } elsif($inf =~ m/\.gff/) {
      ## push @cdsteseqs, $tegff; 
      ($ok,$inh)= openRead($inf);
      while(<$inh>) { next if(/^\W/); my @v=split; 
        my($cl,$ida)=@v[2,9]; 
        my($id)= $ida=~m/ID=([^;\s]+)/ or next;
        # check cl is valid class: CDS, TE/transposon, repeat, ??
        $clt{$cl}++;  $nin++;
        if(my $oc= $idtab{$id}){ $cl.=" $oc"; } $idtab{$id}= $cl; 
      } close($inh);
      
    } elsif($inf =~ /\.tsv/) {
      $nbu += readbuscotab($inf,\%idtab) if($inf =~ /full_table/ and not ($nbu>0 and $inf eq $butab));
      #^ maybe re-read of same file from cds test below.. if -busco opt used
      
    } elsif( ($ok,$inh)= openRead($inf) and $ok) {
    while(<$inh>) { 
      # FIXME for te.gff
      if(/^>(\S+)/){ my $id=$1; 
        my $cl="UNK"; 
        # if(/type=CDS/){ $cl=(m/Class:Transposon|Name=TE:/)?"TE":"CDS"; } 
        # elsif(/class=(\w+)/){ my $tc=$1; $cl=($tc=~/DNA|RC|LINE|SINE|LTR/)?"TE":"UNK"; } 
        
        if(/(?:type|class)=(\w+)/i){ my $tc=$1; 
          if($tc eq 'CDS') {
            $cl= (m/Class:Transposon|Name=TE:/)?"TE":"CDS"; $incds++
          } elsif($tc =~ /^(transposon|TE)/i) {
            $cl=(m/DNA|RC|LINE|SINE|LTR/)?"TE":"UNK"; 
          } elsif($tc=~/(repeat|Satellite)/i) { #  Simple_repeat
            $cl='repeat'; # not "RPT"; # <<< this is known as 'repeat' in covsum, but RPT in covtab, change here?
          } elsif($tc=~/^(un)/i){
            $cl="UNK";           
          } 
        }
        $clt{$cl}++;  $nin++;
        if(my $oc= $idtab{$id}){ $cl.=" $oc"; } $idtab{$id}= $cl; 
        # print $outh "$id\t$cl\n"; $no++; #defer output? add busco if find table
        }
      } close($inh);      
    }
    
    if($incds and $nbu == 0) {
      my($ind,$fn)= ($inf=~m,^(.+)/([^\/]+)$,)? ($1,$2):("./",$inf);
      if( opendir(D,$ind) ) { #<< this not working, fullbusco.tsv is inside subdir run_xxx/
        ($butab)= grep /full_table/, map{ chomp; $_; } readdir(D); closedir(D);
        if($butab) { $nbu += readbuscotab($butab,\%idtab); } # if nbu==0 check @butabx ? match fn cds name?
      }
    }
    
  } 
  
  open($outh,">$outclass");
  for my $id (sort keys %idtab) { print $outh "$id\t",$idtab{$id},"\n"; $no++; } close($outh);
  
  my $clt=join", ", map{ "$_:$clt{$_}" } sort keys %clt;
  warn "#make_idclass: nout=$no ($clt) to $outclass\n" if $debug;
  return($no, $outclass);
}

sub read_idclass {
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

# sub OLDread_idclass {
#   my($crclassf,$allclasses)=@_; 
#   my %crclass=(); my $nid=0;  my($ok,$inh);
#   $allclasses||=0;
#   if($crclassf and ($ok,$inh)= openRead($crclassf) and $ok) { 
#     while(<$inh>){ next if(/^\W/); 
#       my($cr,$crclass,@clx)=split;  # may have more columns .. keep all ie CDS,BUSCO 
#       if($allclasses and @clx){ $crclass=join(" ",$crclass,@clx); }
#       $crclass{$cr}= $crclass || 0;  
#     } close($inh); 
#   }
#   my %crcl=(); # make list of values = classes
#   unless(%crclass){ $crcl{'UNK'} = 0; }
#   else { map{ $crcl{$crclass{$_}}++; } keys %crclass; } # or values %crclass > not uniq
#   my @crclass= sort keys %crcl; # sort by count?
#   my $ncl= @crclass;  $nid=scalar(keys %crclass);
#   warn "# read nid=$nid, nclass=$ncl from $crclassf\n" if($debug);
#   return($nid,\%crclass,\@crclass);
# }  

sub read_chragp {
  my($agpf)=@_; my $nagp=0;
  my($ok,$inh)= openRead($agpf);
  unless($ok){ warn "cant read  $agpf\n"; return 0; }
  while(<$inh>){ 
    next unless(/^\w\S+\t\d/);
    my($cid,$cb,$ce,$ci,$FN,$sid,$sb,$se,$so)=split;  
    next if($FN =~ m/^[NU]/ or not($sid=~/\w/ and $se>0));
    # ($sid)= cleanchrtagIN($sid); # avoid this with NOCLEAN==1 w/ read_chragp
    push @{$agpc{$sid}},[$cid,$cb,$ce,$ci,$sid,$sb,$se,$so]; $nagp++;
  } close($inh);
  for my $s (sort keys %agpc) { 
    my @sg=sort{ $a->[5]<=>$b->[5] } @{$agpc{$s}}; $agpc{$s}=\@sg; 
  }
  warn "# read n=$nagp from $agpf\n" if($debug);
  return($nagp);
}

sub basename { my($f,$suf)=@_; $f=`basename $f`; chomp($f); $f=~s/$suf// if($suf); return($f); }

sub openRead {
  my($fn)=@_; my($ok,$inh)=(0); 
  return(0,undef) unless($fn and -f $fn);
  if($fn =~ /\.gz/) { $ok= open($inh,"gunzip -c $fn |"); } else { $ok= open($inh,$fn); }
  # die "cant read $fn" unless ($ok); # warn? leave to caller
  return($ok,$inh);
}

=item try1

$evigene/scripts/genoasm/gnodes_annotate.pl -debug -ncpu 4 \
 -cds dmag14t1cds.fa -te dmag20sk4tefam.fa -idclass dmag20cdste.idclass -chr genome/dmag19skasm.fa
NOTE: all output to genome/dmag19skasm* .. move?

# read nid=0, nclass=3 from dmag20cdste.idclass
#run /Mypathto/evigene/scripts/genoasm/../genoasm/chrfasta2agp.pl -mingap 10 -chr genome/dmag19skasm.fa -agp genome/dmag19skasm.agp
Use of uninitialized value $_ in pattern match (m//) at /Mypathto/evigene/scripts/genoasm/../genoasm/chrfasta2agp.pl line 78, <$inh> line 1543157.
#fa2agp nin=4193, nok=4193, nout=16771, seq=114.6.mb, gaps=8.2.mb, to genome/dmag19skasm.agp, 

#run blastn  -perc_identity 80 -evalue 1e-9 -task dc-megablast -template_type coding -template_length 18 -outfmt 7 -db genome/dmag19skasm.fa -num_threads 4 -out genome/dmag19skasm.fa-dmag14t1cds.dcblastn  -query dmag14t1cds.fa 
cant read genome/dmag19skasm.fa-dmag14t1cds.dcblastn
 ? no blastn on path : forgot export PATH=$PATH:$nbin << NEED error report
 
#run blastn -task blastn -perc_identity 80 -evalue 1e-9 -outfmt 7 -db genome/dmag19skasm.fa -num_threads 4 -out genome/dmag19skasm.fa-dmag20sk4tefam.blastn  -query dmag20sk4tefam.fa 
cant read genome/dmag19skasm.fa-dmag20sk4tefam.blastn

# write nchr=427

=item try2

 export PATH=$PATH:$nbin
BLAST Database error: No alias or index file found for nucleotide database [genome/dmag19skasm.fa] in search path [/Mypathto/biox/bio-grid/daphmag/gasm20set/dmag20reasm/geneval::]

 $evigene/scripts/genoasm/gnodes_annotate.pl -debug -ncpu 4  -cds dmag14t1cds.fa -te dmag20sk4tefam.fa -idclass dmag20cdste.idclass -chr genome/dmag19skasm.fa

# read nid=30730, nclass=3 from dmag20cdste.idclass

#run /Mypathto/evigene/scripts/genoasm/../genoasm/chrfasta2agp.pl -mingap 10 -chr genome/dmag19skasm.fa -agp genome/dmag19skasm.agp
#fa2agp nin=4193, nok=4193, nout=16771, seq=114.6.mb, gaps=8.2.mb, to genome/dmag19skasm.agp, 

#run blastn  -perc_identity 80 -evalue 1e-9 -task dc-megablast -template_type coding -template_length 18 -outfmt 7 -db genome/dmag19skasm.fa -num_threads 4 -out genome/dmag19skasm.fa-dmag14t1cds.dcblastn  -query dmag14t1cds.fa 
BLAST Database error: No alias or index file found for nucleotide database [genome/dmag19skasm.fa] in search path [/Mypathto/biox/bio-grid/daphmag/gasm20set/dmag20reasm/geneval::]

#run blastn -task blastn -perc_identity 80 -evalue 1e-9 -outfmt 7 -db genome/dmag19skasm.fa -num_threads 4 -out genome/dmag19skasm.fa-dmag20sk4tefam.blastn  -query dmag20sk4tefam.fa 
BLAST Database error: No alias or index file found for nucleotide database [genome/dmag19skasm.fa] in search path [/Mypathto/biox/bio-grid/daphmag/gasm20set/dmag20reasm/geneval::]

# write nchr=427

=item try3

$evigene/scripts/genoasm/gnodes_annotate.pl -debug -ncpu 4 -idclass dmag20cdste.idclass \
   -cds dmag14t1cds.fa -te dmag20sk4tefam.fa -chr genome/dmag19skasm.fa

# read nid=30730, nclass=3 from dmag20cdste.idclass
#run /Mypathto/evigene/scripts/genoasm/../genoasm/chrfasta2agp.pl -mingap 10 -chr genome/dmag19skasm.fa -agp genome/dmag19skasm.agp
#fa2agp nin=4193, nok=4193, nout=16771, seq=114.6.mb, gaps=8.2.mb, to genome/dmag19skasm.agp, 

#run makeblastdb -dbtype nucl -out genome/dmag19skasm -title genome/dmag19skasm -in genome/dmag19skasm.fa

Building a new DB, current time: 01/29/2021 15:02:58
New DB name:   /Mypathto/biox/bio-grid/daphmag/gasm20set/dmag20reasm/geneval/genome/dmag19skasm
New DB title:  genome/dmag19skasm
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 4193 sequences in 1.91085 seconds.

#run blastn  -perc_identity 80 -evalue 1e-9 -task dc-megablast -template_type coding -template_length 18 -outfmt 7 -db genome/dmag19skasm -num_threads 4 -out genome/dmag19skasm-dmag14t1cds.dcblastn  -query dmag14t1cds.fa 
...

#run blastn -task blastn -perc_identity 80 -evalue 1e-9 -outfmt 7 -db genome/dmag19skasm -num_threads 4 -out genome/dmag19skasm-dmag20sk4tefam.blastn  -query dmag20sk4tefam.fa 

$evigene/scripts/genoasm/gnodes_annotate.pl -debug -ncpu 4 -cds dmag14t1cds.fa -te dmag20sk4tefam.fa -idclass dmag20cdste.idclass -chr genome/dmag19skasm.fa
# write nchr=3436

head genome/dmag19skasm.anntab
#scaf	loc	annot	ids
dapmag19sk:LG1	100	LTR,TE	dmag20skm_ter1f930
dapmag19sk:LG1	200	LTR,TE	dmag20skm_ter1f930
dapmag19sk:LG1	2100	LTR,TE	dmag20skm_ter6f3990
dapmag19sk:LG1	2200	CDS,LTR,TE	Dapma7bEVm009261t1,Dapma7bEVm025611t1,dmag20skm_ter6f3990
dapmag19sk:LG1	2300	CDS,LTR,TE	Dapma7bEVm009261t1,Dapma7bEVm025611t1,dmag20skm_ter6f3990
dapmag19sk:LG1	2400	CDS,LTR,TE	Dapma7bEVm009261t1,Dapma7bEVm025611t1,dmag20skm_ter6f3990
dapmag19sk:LG1	2500	CDS,LTR,TE	Dapma7bEVm009261t1,Dapma7bEVm025611t1,dmag20skm_ter6f3990
dapmag19sk:LG1	4300	gap	noid
dapmag19sk:LG1	7800	gap	noid

LTR,TE << want LTR or not? .. prob not
LG1	11400	UNK,Unknown << dup class, from where?
LG1	11500	UNK,Unknown

 724 CDS,Simple_repeat,UNK << Simple_repeat should become RPT class

cat  genome/dmag19skasm.anntab | cut -f3| sort | uniq -c | sort -k1,1nr | less
387420 CDS
86961 gap
42234 CDS,UNK,Unknown
35887 UNK,Unknown
27129 CDS,TE
21471 CDS,LTR,TE
9996 CDS,TE,UNK,Unknown
7512 LTR,TE
4738 CDS,LINE,TE
4520 CDS,LTR,TE,UNK,Unknown
4374 CDS,DNA,TE
4268 DNA,TE
3743 CDS,gap
2968 CDS,LINE,TE,UNK,Unknown
1359 UNK,Unknown,gap
1267 CDS,RC,TE
1262 LINE,TE
1098 LTR,TE,UNK,Unknown
1065 CDS,UNK,Unknown,gap
 901 CDS,DNA,TE,UNK,Unknown
 
   45M Jan 29 15:11 genome/dmag19skasm.anntab
  7.5M Jan 29 15:10 genome/dmag19skasm-dmag20sk4tefam.gff
  8.5M Jan 29 15:10 genome/dmag19skasm-dmag20sk4tefam.blastn
   31M Jan 29 15:09 genome/dmag19skasm-dmag14t1cds.gff
   53M Jan 29 15:09 genome/dmag19skasm-dmag14t1cds.dcblastn
  944K Jan 29 15:03 genome/dmag19skasm.nhr
   49K Jan 29 15:03 genome/dmag19skasm.nin
   29M Jan 29 15:03 genome/dmag19skasm.nsq
  2.0M Jan 29 15:02 genome/dmag19skasm.agp
  729K Jan 29 15:02 genome/dmag19skasm.gaps.gff

  357348 genome/dmag19skasm-dmag14t1cds.gff
   83905 genome/dmag19skasm-dmag20sk4tefam.gff
   12578 genome/dmag19skasm.gaps.gff

==> genome/dmag19skasm-dmag14t1cds.gff <==
##gff-version 3
dapmag19sk:LG1	evg	CDS	2134	2362	216	-	.	ID=Dapma7bEVm025611t1;align=231;ident=93;
dapmag19sk:LG1	evg	CDS	2134	2402	258	-	.	ID=Dapma7bEVm009261t1;align=271;ident=95;
dapmag19sk:LG1	evg	CDS	2286	2425	133	-	.	ID=Dapma7bEVm025611t1;align=140;ident=95;
dapmag19sk:LG1	evg	CDS	2286	2425	135	-	.	ID=Dapma7bEVm009261t1;align=140;ident=96;
dapmag19sk:LG1	evg	CDS	11665	11817	129	+	.	ID=Dapma7bEVm008970t1;align=157;ident=82;
dapmag19sk:LG1	evg	CDS	11666	11753	77	+	.	ID=Dapma7bEVm020890t1;align=88;ident=87;
dapmag19sk:LG1	evg	CDS	11718	11825	92	+	.	ID=Dapma7bEVm008970t1;align=109;ident=84;
dapmag19sk:LG1	evg	CDS	12021	12178	146	+	.	ID=Dapma7bEVm008970t1;align=161;ident=91;

==> genome/dmag19skasm-dmag20sk4tefam.gff <==
##gff-version 3
dapmag19sk:LG1	evg	LTR	1	145	128	-	.	ID=dmag20skm_ter1f930;align=148;ident=87;
dapmag19sk:LG1	evg	LTR	2042	2420	368	+	.	ID=dmag20skm_ter6f3990;align=379;ident=97;
dapmag19sk:LG1	evg	LTR	10792	11257	416	-	.	ID=dmag20skm_ter5f2284;align=474;ident=87;
dapmag19sk:LG1	evg	Unknown	11396	11482	74	+	.	ID=dmag20skm_unr1f168;align=88;ident=84;
dapmag19sk:LG1	evg	Unknown	11411	11482	61	-	.	ID=dmag20skm_unr1f201;align=72;ident=86;
dapmag19sk:LG1	evg	Unknown	11575	11768	159	+	.	ID=dmag20skm_unr1f136;align=194;ident=81;
dapmag19sk:LG1	evg	Unknown	11615	11836	197	+	.	ID=dmag20skm_unr1f168;align=223;ident=88;

==> genome/dmag19skasm.gaps.gff <==
dapmag19sk:LG1	chr2agp	gap	4240	4249	10	.	.	it=2
dapmag19sk:LG1	chr2agp	gap	7741	7763	23	.	.	it=4
dapmag19sk:LG1	chr2agp	gap	14309	14324	16	.	.	it=6
dapmag19sk:LG1	chr2agp	gap	17245	17273	29	.	.	it=8
dapmag19sk:LG1	chr2agp	gap	26091	27452	1362	.	.	it=10

=item try4 fixed TE classes, added busco

cut -f3 genome/dmag19skasm.anntab | sort | uniq -c | sort -k1,1nr | less
368872 CDS
86961 gap
59939 CDS,TE
42190 CDS,UNK
36053 UNK
19820 CDS,TE,UNK
18548 CDS,busco   << should be Ku uniq depth set
15289 TE
3633 CDS,gap
2546 TE,UNK
1427 CDS,UNK,repeat
1373 UNK,gap
1187 CDS,TE,UNK,repeat
1068 CDS,UNK,gap
1044 CDS,TE,gap
 697 UNK,repeat
 597 TE,gap
 494 CDS,TE,UNK,gap
 116 CDS,UNK,busco
 110 CDS,busco,gap
  86 TE,UNK,repeat
  72 TE,UNK,gap
  48 CDS,UNK,gap,repeat
  35 CDS,TE,busco
  33 UNK,gap,repeat
  18 CDS,TE,UNK,gap,repeat
   5 CDS,TE,UNK,busco
   2 CDS,TE,busco,gap
   2 TE,UNK,gap,repeat
   
=cut
