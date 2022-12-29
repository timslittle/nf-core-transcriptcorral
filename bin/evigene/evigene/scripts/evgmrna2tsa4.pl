#!/usr/bin/perl
# evgmrna2tsa4.pl stub to replace   evgmrna2tsa2.pl;   

=item about

  EvidentialGene evgmrna2tsa.pl is OBSOLETE with 2020 March release
  It is replaced by  genes/trclass2pubset.pl and genes/pubset2submit.pl
  
  Usage now is  
    (a) prot/tr2aacds4.pl, which calls subprogram trclass2pubset.pl,
        output to okayset/ with public ids (pubids) and reformated sequences,
        in same fashion of evgmrna2tsa.pl -onlypubset 
    (b) evgpipe_sra2genes4v.pl, calls trclass2pubset.pl and pubset2submit.pl,
        following Step 7 with tr2aacds4, 
        sra2genes Steps 10,11 send outputs to publicset/ (via trclass2pubset) 
        and submitset/ (via pubset2submit)
=cut

use strict;
use constant VERSION => '2020.04.15';  # obsolete

die "EvidentialGene evgmrna2tsa VERSION ",VERSION,"
  is now OBSOLETE.  See replacements via evigene/scripts/evgpipe_sra2genes4v.pl
  or evigene/scripts/prot/tr2aacds4.pl 
  It made mRNA fasta and annotation table for TSA submit via tbl2asn
  using mRNA classes from Evigene tr2aacds, gene product names, vecscreen.

  Usage now is  
    (a) prot/tr2aacds4.pl, which calls subprogram trclass2pubset.pl,
        output to okayset/ with public ids (pubids) and reformated sequences,
        replacing evgmrna2tsa.pl -onlypubset 
        
    (b) evgpipe_sra2genes4v.pl, calls trclass2pubset.pl and pubset2submit.pl,
        following Step 7 with tr2aacds4, 
        sra2genes Steps 10,11 send outputs to publicset/ (via trclass2pubset) 
        and submitset/ (via pubset2submit), replacing evgmrna2tsa.pl -onlysubmit
        
";