#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 check.pl [-n] [-p {alu|line|both|all}] [interval_start interval_end]
 Count transcripts in novel genes with a certain repeat content
 Options:
   -n  count non-coding genes instead of coding
   -p  print the entries in the given interval for the 
       given repeat types
   Defaults:
     interval_start=0.1
     interval_end=100.1
/;
my $file='chess2_rep_data.tab';
umask 0002;
getopts('np:o:') || die($usage."\n");
my $noncoding=$Getopt::Std::opt_n;
my $print=$Getopt::Std::opt_p;
if ($print) {
  $print=lc($print);
  die("${usage}Error: invalid -p value '$print'!\n") 
     unless ($print eq 'alu' || $print eq 'line' || $print eq 'both' || $print eq 'all');
}

my ($istart, $iend)=@ARGV;
$istart=0.1 if length($istart)==0;
$iend=100.1 unless $iend;
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
print "id\tgstatus\tgtype\tstatus\ttype\tasm\tchr\tstrand\texcount\texlen\tcdlen\texreps\tcdreps\texcov\tcdcov\talu_excov\talu_cdcov\tline_excov\tline_cdcov\n"
  if $print;
open(INF, $file) || die("Error opening $file\n");
my ($aluix, $lineix)=$noncoding ? (15, 17) : (16, 18);

while (<INF>) {
 next if m/^id\t/;
 chomp;
 my @t=split(/\t/);
 my $pass=0;
 #my $dbg=($t[0] eq 'CHS.696.1');
 if ($noncoding) {
   $pass=($t[1] eq 'known_fantom' || $t[1] eq 'novel') && ($t[2] ne 'protein_coding')
     && ($t[13]>=$istart && $t[13]<$iend);
  #print STDERR "pass=$pass\n" if $dbg;
 } else {
   $pass=($t[1] eq 'known_fantom' || $t[1] eq 'novel') && ($t[2] eq 'protein_coding')
     && ($t[10]>0) && ($t[14]>=$istart && $t[14]<$iend);
 }
 #print STDERR "pass=$pass ; \$t[$aluix]=".$t[$aluix]."(istart=$istart, iend=$iend)\n" if $dbg;
 if ($pass) {
   $excount++ if ($t[13]>=$istart && $t[13]<$iend);
   $cdcount++ if ($t[14]>=$istart && $t[14]<$iend);
   if ($print) {
     if ($print eq 'all') { print "$_\n"; }
     else {
       my $np=1;
       if (($print eq 'alu' || $print eq 'both') && ($t[$aluix]>=$istart) && ($t[$aluix]<$iend)) {
           print "$_\n";
           $np=0;
       }
       print "$_\n" if ($np && ($print eq 'line' || $print eq 'both') && ($t[$lineix]>=$istart) && ($t[$lineix]<$iend));
     }
   }
 }
} #for each INF line
close(INF);
print STDERR "Count for exon coverage \[$istart, $iend): $excount\n";
print STDERR "Count for CDS coverage  \[$istart, $iend): $cdcount\n";

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

