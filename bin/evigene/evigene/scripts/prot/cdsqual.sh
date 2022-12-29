#!/bin/bash
# cdsqual.sh == aaqual.sh for seq.cds input, same cols as aa.qual: c2 = aasize, c5 = trsize=cdssize, c6=1-cdsize
# aaqual.sh any*.aa.gz : count protein sizes with faCount (dna), with tr for XXX gap count
# output 3-col table: ID, aatotal, XXgaps
# add on >aalen qual columns; ** Only good for evigene cdna_bestorf.aa files
# >whitefly:vel2k25Loc100074t1 aalen=41,63%,complete; clen=199; strand=-; offs=144-19; 
# ** FIXME2: w/ new aaqual, set col2 == size - gap, not sizetotal, and let user add gaps in
## FIX for ensembl pep.fa with gene: instead of oid= gene:ENSDARG00000030494
## FIX3: add strand to span/offs !!
## FIX4: Selcstop=1 transfer to aaqual string
## FIX4: cdsoff= alias offs=; cxlen=cdslen/trlen alias to clen=; dang also cdsoffs=nnn
## cds: stop_codons = qw(TAA TAG TGA);

export iscds=1; export off=1;
dostat=0; if [ "X$stat" != "X" ]; then dostat=$stat; fi
outd=0; if [ "X$outdir" != "X" ]; then outd=$outdir; fi
if [ "X$span" != "X" ]; then export off=$span; fi; 
if [ "X$outcds" = "X" ]; then export outcds=0; fi;

inaa=$*
osuf=aa; if [ $outcds -gt 0 ]; then osuf=cds; fi

for az in $inaa; do {
  TCAT=cat; 
  nogz=`echo $az | sed 's/.gz//;'`; if [ $az != $nogz ]; then TCAT="gunzip -c"; fi
  nam=`echo $az | sed 's/.gz//; s/\.qual//;'`; 
  if [ $outd != 0 ]; then nam=`basename $nam | sed "s,^,$outd/,;"`; fi
  namout=`echo $nam | sed "s/.cds//; s/$/.$osuf.qual/;"`; 
  
  if [ ! -f $namout ]; then 
  $TCAT $az | perl -ne \
'if(/^>(\S+)/) { puta() if($d); $d=$1; ($al)=m/aalen=([^;\s]+)/; $al||="na"; 
if(/Selcstop=\w/) { $al.=",selcstop" unless($al=~/selc/); } 
($cl)=m/clen=(\d+)/; unless($cl) { ($cdl,$cl)=m/cxlen=(\d+).(\d+)/; } $cl||=0; 
$aas=$aam=$aat=$aag=0; 
unless(($oid)=m/oid=([^;\s]+)/) { ($oid)=m/gene[=:]([^;\s]+)/; } $oid||="noid";
$ofs=0; if($doff){ ($or)=m/strand=(.)/; $or||="."; ($ofs)=m/\b(?:offs|cdsoff\w*)=([\d-]+)/; 
 $clOFF= $cl . ($ofs)?"\t$ofs:$or":"\t0";} 
} elsif(/[\w\*]/){  $lastline=$_;
if(0) { if($iscds) { $aas= (m/(TGA|TAG|TAA)$/i)?1:0; } else { $aas=(m/\*$/)?1:0; } }
if($iscds) { $aat += tr/ACGTacgt/ACGTacgt/; $aag+= tr/Nn/Nn/; if($aam==0){ $aam=(/^ATG/i)?1:-1; } } 
else { $aat += tr/A-WYZa-wyz/A-WYZa-wyz/; $aag += tr/Xx\*/Xx\*/; if($aam==0){ $aam=(/^M/)?1:-1; } } }
END{ puta(); } BEGIN{  $doff=$ENV{off}||0; $doid=$ENV{oid}||0; $iscds=$ENV{iscds}||0; $outcds=$ENV{outcds}||0; }
sub puta { 
 $cdsw=$aat+$aag; $aat=int($cdsw/3); 
 if($iscds) { 
  $aas= ($lastline=~m/(TGA|TAG|TAA)$/i)?1:0; 
  $aat-- if($aas);
 } else { 
  $aas=($lastline=~m/\*$/)?1:0; $aag-- if($aas);
 }
 $cl ||=$cdsw; $ofs ||= "1-$cdsw";
 $clOFF= $cl . (($ofs)?"\t$ofs:$or":"\t0");
 if($al eq "na"){ 
   $part=($aas==1 && $aam==1)?"complete":($aam==1)?"partial3":($aas==1)?"partial5":"partial"; 
   $al=$aat.",99%,$part"; } 
 $clOFF.="\t$oid" if($doid); 
 $aclen=($outcds)?$cdsw:$aat;
 print join("\t",$d,$aclen,$aag,$al,$clOFF)."\n"; }' \
  > $namout

  fi

} done

