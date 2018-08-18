#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 check.pl interval_start interval_end
/;
my $file='chess2_rep_moreinfo.tab';
umask 0002;
getopts('o:') || die($usage."\n");
my ($istart, $iend)=@ARGV;
die($usage."\n") unless $iend;
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my $count;
#0     1       2     3     4   5   6     7      8      9    10    11      12    13    14
#id gstatus gtype status type asm chr strand excount exlen cdlen exreps cdreps excov cdcov
open(INF, $file) || die("Error opening $file\n");
while (<INF>) {
 next if m/^id\t/;
 chomp;
 my @t=split(/\t/);
 if ($t[3] eq 'known_fantom' || $t[3] eq 'novel') {
   $count++ if ($t[13]>=$istart && $t[13]<$iend);
 }
}
close(INF);
print "Count for \[$istart, $iend): $count\n";
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

