#!/usr/bin/env perl
# evg_buscogenesum.pl
# usage: env dotab=1|0 summary=busco.sum.txt evg_buscogenesum.pl buscof/full_table*.tsv
# or:  cat manybusco*/full_table*.tsv | env dotab=1|0 summary=busco.sum.txt evg_buscogenesum.pl

use strict;

## FIXME: stdin not infile, need to check #full_header
# full_table_buplplIAJW.mnpub.tsv 
# # BUSCO version is: 2.0.1 
# # The lineage dataset is: embryophyta_odb9 (Creation date: 2016-02-13, number of species: 30, number of BUSCOs: 1440)
# # To reproduce this run: python /home/ux455375/bio/orthodb/BUSCOv2.py -i plIAJW.mnpub.aa -o buplplIAJW.mnpub -l /home/ux455375/scratchn/chrs/pubdata/orthodb/plants/ -m proteins -c 12 -sp arabidopsis
# #
# # Busco id	Status	Sequence	Score	Length
# EOG09360002	Duplicated	AmeargEVm000042t1	2348.2	1265

my($inf,$nokgene,$nbadgene,$ntrin)=(0) x 9;
my(@hdr,%tab,%did,%didgene,%al,%sc,%bcla,%cla);  
while(<>) { 

  if(/^#/) {
    if(/^# BUSCO/){ 
      putone($inf) if($inf); 
      @hdr=(); $inf=$ARGV; # STDIN $ARGV == '-' always?
      if($inf =~ /^(stdin|\-)/i) { $inf= "busco_full_table.tsv"; }
    }
    push @hdr,$_;
    if(/^# To reproduce this run:/) { 
      my($newinf)= m/ \-o (\S+)/?$1:"busco_full_table.tsv"; #? $ARGV 
      if($inf =~ /^(stdin|\-)/i) { $inf= $newinf; }
    }
  }
  
  next if(/^\W/); my($bd,$cla,$td,$sc,$al)=split; 
  map{ $_ ||= 0; } ($td,$sc,$al); 
  my $gid=$td; $ntrin++ unless($td eq "0");
  if( $gid=~s/t\d+$// ) { $nokgene++; } else { $nbadgene++ unless($td eq "0"); }
  $didgene{$bd}{$gid}++;  
  $tab{$bd}{$gid}{$td}=[$bd,$cla,$td,$sc,$al];
  if($sc > $sc{$bd}){ $sc{$bd}=$sc; $al{$bd}=$al; }
  unless($did{$bd}++){ $cla="Complete" if($cla eq "Duplicated"); $bcla{$bd}=$cla; $cla{$cla}++; }
}

putone($inf) if($inf);

sub putone { 
  #  $okgenes=($nokgene > $nbadgene and $nbadgene < 0.10*$ntrin)?1:0;
  putsum($ENV{summary}); 
  puttab($inf) if($ENV{dotab}); 
  %tab=%did=%didgene=%al=%sc=%bcla=%cla=(); @hdr=();
  $nokgene=$nbadgene=$ntrin=0;
}

sub puttab { 
  # out == input full_table.tsv, dont overwrite
  my($out)=@_; my $oh; 
  if($out){ open($oh,">$out.evg"); } else { return 0; }
  my $okgenes=($nokgene > $nbadgene and $nbadgene < 0.10*$ntrin)?1:0;
  my @bd=sort keys %did; my $nt=@bd;  
  my @hd=qw(BuscoID Status TranscriptID Score Align  GeneID);
  map{ print $oh $_; } grep{ not m/Status/ } @hdr;
  print $oh "#".join("\t",@hd)."\n";
  for my $bd (@bd) { 
    my @gd= sort keys %{$tab{$bd}};  my $ng=@gd;
    for my $gd (@gd) { 
      my @td= keys %{$tab{$bd}{$gd}}; my $ntr=@td;
      @td= sort { $tab{$bd}{$gd}{$b}->[4] <=> $tab{$bd}{$gd}{$a}->[4] or $a cmp $b } @td;
      
      if($okgenes) {
      #* write only top tr/gene, append other trids in row, w/ score?
      my($tda)= shift @td; my @tdb=();
      for my $td (@td) {
        my($bd,$cla,$td,$sc,$al)= @{ $tab{$bd}{$gd}{$td} };
        $td=~s/$gd//; push @tdb,"$td:$al";
      }
      my($bd,$cla,$td,$sc,$al)= @{ $tab{$bd}{$gd}{$tda} };
      my $cc=$cla; $cc="Complete" if($cla eq "Duplicated" and $ng == 1);
      my $tdb= ($gd eq "0")? $gd : (@tdb>0) ? "$gd:${ntr}alts,".join(",",@tdb) : "$gd:${ntr}alt";
      print $oh join("\t",$bd,$cc,$tda,$sc,$al,$tdb)."\n"; 
      } else {
      for my $td (@td) {
      my($bd,$cla,$td,$sc,$al)= @{ $tab{$bd}{$gd}{$td} };
      my $cc=$cla; $cc="Complete" if($cla eq "Duplicated" and $ng == 1);
      print $oh join("\t",$bd,$cc,$td,$sc,$al)."\n"; 
      }
      }
    }
  }
  close($oh);
}

sub putsum { 
  my($out)=@_; my $oh; 
  if($out){ open($oh,">>$out"); } else { $oh=*STDOUT; } # append out
  my $okgenes=($nokgene > $nbadgene and $nbadgene < 0.10*$ntrin)?1:0;
  # map{ print $oh $_; } grep{ not m/Status/ } @hdr;#? not here? 
  map{ print $oh $_; } grep{ m/reproduce/ } @hdr;#? not here? 
  my @cla=qw(Complete Single Duplicated Fragmented Missing);
  my @bd=sort keys %did; my $nt=@bd;  
  my($ina) = $inf =~ m,([\w\.-]+)/[\w\.-]+$,; $ina=~s/run_bu//; 
  my($tsc,$tal)=(0,0); for my $bd (@bd) { $tsc+=$sc{$bd}; $tal+=$al{$bd};  } 
  my($aal, $asc)=map{ sprintf "%.1f", $_/$nt; } ($tal,$tsc); 
  my($ng,$ntr)=(0,0);
  for my $bd (@bd) { 
    my @gd= sort keys %{$didgene{$bd}}; 
    map{ $ng++; $ntr += $didgene{$bd}{$_}; } @gd;
    if($bcla{$bd} eq "Complete") {
      my $cc=(@gd>1)?"Duplicated":"Single"; $cla{$cc}++; 
      }
    }
  my $scla=join",", map{ my $c=$cla{$_}||0; "$_:$c"; } @cla; 
  my $pcla=join",", map{ my $c=$cla{$_}||0; my $l=substr($_,0,1); sprintf("$l:%.1f%%", 100*$c/$nt); } @cla;
  $pcla =~ s/,S/\[S/; $pcla=~s/,F/\],F/;
  my $ngt=($okgenes)?"ngene:$ng,ntr:$ntr":"ntr:$ntr,notgeneids";
  print $oh "$ina:\t","$pcla,n:$nt, $ngt\n";
  print $oh "$ina:\t","$scla,n:$nt, $ngt, align:$aal, score:$asc\n"; 
}