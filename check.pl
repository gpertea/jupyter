#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 check.pl interval_start interval_end
/;
my $file='chess2_rep_data.tab';
umask 0002;
getopts('o:') || die($usage."\n");
my ($istart, $iend, $list)=@ARGV;
die($usage."\n") unless $iend;
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my ($excount, $cdcount);
#0     1       2     3     4   5   6     7      8      9    10    11      12    13    14
#id gstatus gtype status type asm chr strand excount exlen cdlen exreps cdreps excov cdcov
# 15=alu_excov, 16=alu_cdcov, 17=line_excov, 18=line_cdscov
open(INF, $file) || die("Error opening $file\n");
while (<INF>) {
 next if m/^id\t/;
 chomp;
 my @t=split(/\t/);
 if (($t[1] eq 'known_fantom' || $t[1] eq 'novel') &&
     ($t[10]>0) ) {
   $excount++ if ($t[13]>=$istart && $t[13]<$iend);
   $cdcount++ if ($t[14]>=$istart && $t[14]<$iend);
   if ($list) {
     print "$_\n" if (lc($list) eq 'alu' && $t[16]>=$istart && $t[16]<$iend);
     print "$_\n" if (lc($list) eq 'line' && $t[18]>=$istart && $t[18]<$iend);
     print "$_\n" if (lc($list) eq 'all' && $t[14]>=$istart && $t[14]<$iend);
   }
 }
}
close(INF);
print STDERR "Count for exon coverage \[$istart, $iend): $excount\n";
print STDERR "Count for CDS coverage  \[$istart, $iend): $cdcount\n";

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

