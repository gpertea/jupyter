#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 collect_rep_ovls.pl qexons.bed tabix_QR_output [transcripts.gff]
 
 qexons.bed      : the input that was used as query regions for tabix
 tabix_QR_output : could be '-', is the output of tabix -QR 
 
 The repeat coverage info for each transcript will be printed at stdout
 as tab delimited columns:
 t_id chr strand excount exlen cdslen exovls cdsovls excov% cdscov%
                 Alu-excov% Alu-cdcov% SINE_LINE-excov% SINE_LINE_cdcov%
 
 If transcripts.gff is given, 5 more columns are inserted after t_id:
 
 (gene)STATUS, GENE_TYPE, STATUS, TYPE, ASSEMBLED
 ..as found in the the gff file for the parent gene and the transcript.
/;
umask 0002;
getopts('o:l') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $lonly=$Getopt::Std::opt_l;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
die($usage."\n") if (@ARGV<2 || @ARGV>3);
my $gff;
my %ginfo; #gene_id => [STATUS, GENE_TYPE]
my %tinfo; #t_id=> [ STATUS, TYPE, ASSEMBLED ]
if (@ARGV==3) {
 $gff=pop(@ARGV);
 open(GFF, $gff) || die("Error opening file $gff!\n");
 while(<GFF>) {
  if (m/gene\t/) {
   chomp;
   my ($gid)=(m/\bID=([^;]+)/);
   my ($gs)=(m/\bSTATUS=([^;]+)/);
   my ($gt)=(m/\bGENE_TYPE=([^;]+)/);
   die("Error: GENE_TYPE or STATUS not found for gene $gid\n") if (!$gt || !$gs);
   die("Gene $gid already processed (duplicate?)\n") if exists($ginfo{$gid});
   $ginfo{$gid}=[$gs, $gt];
  }
  elsif (m/\bID=([^;]+)/) {
    my $id=$1;
    chomp;
    my ($s)=(m/\bSTATUS=([^;]+)/);
    my ($t)=(m/\bTYPE=([^;]+)/);
    my ($a)=(m/\bASSEMBLED=([^;]+)/);
    my ($pg)=(m/\bParent=([^;]+)/);
    $t='undetermined' unless $t;
    if ($t eq 'protein_coding') { $t='coding' }
    elsif ($t eq 'non_coding') { $t='noncoding' }
    elsif ($t eq 'protein_coding_uncertain') { $t='coding_maybe' }
    die("Error: transcript $id duplicated in $gff ?\n") if exists($tinfo{$id});
    $tinfo{$id}=[$s, $t, $a, $pg];
  }
 }
 close(GFF);
}

my $qxbed=shift(@ARGV);
die($usage."Error: no BED file $qxbed found!\n") unless -f $qxbed;
open(BED, $qxbed) || die ("Error opening file $qxbed\n");
my @tdata; #currently open transcripts, same order they were encountered
  # [tID, chr, strand, numexons, [exons], [CDS_segs], [repClass_exovls], [repClass_cdovls], [repeats_exovls], [repeats_cdsovls], [Alu_exovls], [Alu_cdovls], [SINE_exovls], [SINE_cdovls]]
  #   0     1     2       3        4[]      5[]            6[]            7[]                  8[]               9[]            10[]           11[]             12[]       13[]
my %th; # tID => $tdata[] entry 
#my $skipreg; #skip repeated region (common exon)
my @toflush; #list of tIDs to flush/write to output
my $debug=0;
#print header
if ($gff) {
print join("\t", 'id','gstatus','gtype', 'status', 'type', 'asm', 'chr', 'strand', 'excount', 'exlen', 
 'cdlen', 'exreps', 'cdreps', 'excov', 'cdcov', 'alu_excov', 'alu_cdcov', 'sine_excov', 'sine_cdcov')."\n";
}
else {
print join("\t", 'id', 'chr', 'strand', 'excount', 'exlen', 
 'cdlen', 'exreps', 'cdreps', 'excov', 'cdcov', 'alu_excov', 'alu_cdcov', 'sine_excov', 'sine_cdcov')."\n";
}
my ($ctd, $ce1, $ce2, $co1, $co2); #current transcript/exon/CDS data
# Example entry in tabix output:
# ##Qry|chrY:57201135-57202145
# chrY	57201330	57201404	SINE|Alu|AluJo	491	-
# chrY	57201487	57201792	SINE|Alu|AluSp	2343	+
while (<>) {
  if (m/^##Qry\|(.+)/) {
    my $qreg=$1;
    writeTData() if (@toflush);
    ($ctd, $ce1, $ce2, $co1, $co2)=getQryData($qreg);
    next;
  }
  # next if ($skipreg);
  #--- process overlap of current exon ($ce1-$ce2) and orf ($co1, $co2) of $ct
  chomp;
  my @r=split(/\t/); # #r[1,2] = repeat coordinates
  $r[1]++; #adjust BED start coordinate
  my ($v1, $v2)=getOvl($ce1, $ce2, $r[1], $r[2]);
  #overlap between exon and repeat region
  if ($v1) {
     my $ovl=$v2-$v1+1;
     addRepOvl($ctd->[6], $r[3], $v1, $v2, $ctd->[8], $ctd->[10], $ctd->[12]) if $ovl>3;
  } else { 
     die("Error: no actual overlap between $ce1-$ce2 and $r[1]-$r[2] ?!\n") 
  }
  if ($co1) {
    ($v1, $v2)=getOvl($co1, $co2, $r[1], $r[2]);
    if ($v1 && $v2-$v1+1>3) {
    #### DEBUG
    #  my $pdebug=$debug;
    #  $debug=0;
       addRepOvl($ctd->[7], $r[3], $v1, $v2, $ctd->[9], $ctd->[11], $ctd->[13]);
    #  $debug=$pdebug;
    }
  }
}
getQryData('');
writeTData() if (@toflush);

close(BED);

exit(0);

#---------- Subroutines -----------

sub getQryData {
  #retrieves BED entries for the current qry interval
  my $qr=$_[0]; # e.g.: chrY:57201135-57202145
  while (<BED>) {
    chomp;
    # chr1	183921	184258	CHS.9.4|183922-184158	4	+	CDS:183922-184158
    my @t=split(/\t/);
    my @ao=split(/\|/,$t[3]); # (t_id, CDS_ovl_coords)
    $t[1]++; #+1 BED start coordinate adjustment
    my $tid=$ao[0];
    my ($o1, $o2);
    if (@ao>1) { #exon has CDS region
        ($o1, $o2)=($ao[1]=~m/^(\d+)\-(\d+)/);
        die("Error parsing tID, orfseg from BED line:(ao[1]=$ao[1])\n$_\n") unless $o2 && $o1 && $o2>=$o1;
    }
    my $rtd=$th{$tid};
    if (!$rtd) {
      #create tdata storage
      $rtd=[ $tid, $t[0], $t[5], $t[4], [], [], [], [], [], [], [], [], [], [] ];
      push(@tdata, $rtd);
      $th{$tid}=$rtd;
    }
    #update transcript data
    push(@{$rtd->[4]}, $t[1].'-'.$t[2]); #store exon coordinates
    push(@{$rtd->[5]}, $o1.'-'.$o2) if $o1; # store CDS coordinates overlapping the exon
    push(@toflush, $tid) if ($rtd->[3] == scalar(@{$rtd->[4]})); # all exons collected?
    
    if ($qr eq $t[0].':'.$t[1].'-'.$t[2]) {
      #found region (exon) match
      return ($rtd, $t[1], $t[2], $o1, $o2); #return transcript data collected so far for this exon
    }
    #skipping exon with no rep overlaps, update exons seen for this tid
  }
  die("Error: EOF reached while looking for $qr in BED file?\n") if $qr;
}

sub writeTData {
  foreach my $tid (@toflush) {
     my $td=$th{$tid};
     die("Error: no tdata found for $tid?!\n") unless $td;
     my $exonlen=getSegLen($$td[4]);
     my $cdslen=getSegLen($$td[5]);
     if ($gff) {
       my @ts=@{$tinfo{$tid}};
       my $gid=pop(@ts);
       my @gd=@{$ginfo{$gid}};
       unshift(@ts, @gd);
       print join("\t", $$td[0], @ts, @{$td}[1..3]);
     } else {
       print join("\t", @{$td}[0..3]); # join(',', @{$$td[4]}), getSegLen($$td[5]), 
     }
     my ($excov, $cdcov)=(showOvl($$td[8], $exonlen, $tid), showOvl($$td[9], $cdslen, $tid));
     my ($aluexcov, $alucdcov)=(showOvl($$td[10], $exonlen, $tid),  showOvl($$td[11], $cdslen, $tid));
     my ($sinexcov, $sincdcov)=(showOvl($$td[12], $exonlen, $tid),  showOvl($$td[13], $cdslen, $tid));
     print "\t".join("\t", $exonlen, $cdslen, 
              showRepOvl($$td[6], $exonlen, $tid), showRepOvl($$td[7]),
              $excov, $cdcov, $aluexcov, $alucdcov, $sinexcov, $sincdcov
              )."\n";
     die("Error: Alu/SINE/LINE exon cov > total repeat exon cov for $tid ?!\n") if ($aluexcov>$excov || $sinexcov>$excov);
     die("Error: Alu/SINE/LINE CDS cov > total repeat CDS cov for $tid ?!\n") if ($alucdcov>$cdcov || $sincdcov>$cdcov);
     delete($th{$tid});
     @tdata= grep { $_->[0] ne $tid } @tdata;
  }
  @toflush=();
}

sub getSegLen {
 my $r=$_[0];
 my $s=0;
 return 0 if @$r<1;
 foreach my $reg (@$r) {
   my ($v1, $v2)=split(/\-/,$reg);
   $s+=int($v2)-int($v1)+1;
 }
 return $s;
}

sub showRepOvl {
 my ($repd, $xlen, $tid)=@_; #$repd = list of [repClass, [[r1,r2], ...]
 return '.' if @$repd == 0;
 my @rlen; #list of repClass:len
 foreach my $rdata (@$repd) {
   my ($rclass, $rlist)=@$rdata;
   my $l;
   foreach my $rr (@$rlist) {
      $l+=$$rr[1]-$$rr[0]+1;
   }
   die ("Error for $tid : class ovl len $l > exon len $xlen !\n")
     if $xlen && $l > $xlen;
   push(@rlen, $rclass.':'.$l);
 }
 return join(',', @rlen);
}

sub showOvl {
 my ($regs, $xlen, $tid)=@_; #$regs = list of [r1,r2]...
 return '0' if ($xlen==0 || @$regs == 0);
 my $rl; #summing up all intervals here
 foreach my $rr (@$regs) {
      $rl+=$$rr[1]-$$rr[0]+1;
 }
 die ("Error for $tid : ovl len $rl > total exon len $xlen !\n")
    if $xlen && $rl > $xlen;
    
 return $lonly ? $rl : sprintf('%.1f',(100.0*$rl/$xlen));
}

sub addRepOvl {
 my ($rset, $rn, $v1, $v2, $arset, $aluset, $sineset)=@_;
 #rset=set of intervals to update for a specific repeat class/type
 #     list of [repname, [ [ovl1-ovl2], ...] ]
 #rn is an entry from tabix overlap report, e.g. "SINE|Alu|AluSp"
 #$v1-$v2 =overlap of exon with this repeat
 #arset = set of intervals overlapped by *any* repeat
 #aluset = set of intervals overlapped by Alu repeats
 #sineset = set of intervals overlapped by LINE/SINE repeats
 my @rps=(split(/\|/, $rn)); # e.g. ('SINE', 'Alu', 'Alusp')
 my $rclass= (@rps>1 && $rps[0] ne $rps[1]) ? join('|',@rps[0,1]) : $rps[0];
 my $added=0;
 #my $r_added=0;
 if ( @$rset>0 ) {
   foreach my $rd (@$rset) {
     #see if any existing should be updated by this v1-v2 overlap
     my ($srclass, $reglist) = @$rd; 
     if ($srclass eq $rclass) {
       #get union interval ($r1, $r2)U($s1,$s2)
       ovlRegList($v1, $v2, $reglist);
       $rd=[$srclass, $reglist];
       $added=1;
       last;
     }
   }
 }
 unless ($added) { #new overlap class
   #if ($debug) { print STDERR "adding new overlap $rclass: $v1-$v2\n"; }
   push(@$rset, [$rclass, [[$v1, $v2]]] );
 }
 ovlRegList($v1, $v2, $arset) if ($arset); #update list of any repeat regions
 ovlRegList($v1, $v2, $aluset) if ($aluset && $rclass=~m/\|Alu/); #update list of Alu repeat regions
 ovlRegList($v1, $v2, $sineset) if ($sineset && $rclass=~m/^[LS]INE/);
}

sub ovlRegList {
 my ($v1, $v2, $reglist)=@_;
 if (@$reglist==0) { #new repeat overlap list
   push(@$reglist, [$v1, $v2]);
   return;
 }
 my @tomerge; #indexes of @$reglist entries to merge
 my $i=0;
 my $ins=-1;
 #if ($debug) {
 #  print STDERR "add/merge overlap $v1-$v2 to existing list: ", 
 #    join(',', (map { $$_[0].'-'.$$_[1] } @$reglist)),"\n";
 #}
 foreach my $reg (@$reglist) {
   my ($r1,$r2)=@$reg;
   #if ($debug) {
   #  print STDERR "     loopDBG: testing $v1-$v2 vs $r1-$r2 (last if $r1>$v2)\n";
   # }
   $ins=$i if $r1>$v1;
   last if $r1>$v2;
   if ($r1 <= $v2 && $v1 <= $r2) {
     #overlap with this region
     ##-DEBUG
     ##my ($p1, $p2)=($r1,$r2);
     push(@tomerge, $i);
     $r1 = $v1 if $v1<$r1;
     $r2 = $v2 if $v2>$r2;
     $reg=[$r1, $r2];
     ##if ($debug) {
     ##  print STDERR " ..merging $p1-$p2 with $v1-$v2 = $r1-$r2\n";
     ##}
   }
   $i++;
 }
 if (@tomerge>0) {
   #join all these overlapping intervals
   if (@tomerge>1) {
     my $bi=shift(@tomerge);
     foreach my $di (@tomerge) {
       $$reglist[$bi]->[0]=$$reglist[$di]->[0] if $$reglist[$di]->[0]<$$reglist[$bi]->[0];
       $$reglist[$bi]->[1]=$$reglist[$di]->[1] if $$reglist[$di]->[1]>$$reglist[$bi]->[1];
     }
     splice(@$reglist, $tomerge[0], scalar(@tomerge));
   }
 }
 else { #no overlaps, add this region in the right place
   if ($ins>=0) {
      splice(@$reglist, $ins, 0, [$v1, $v2]);
   }
   else {
      push(@$reglist, [$v1, $v2]);
   }
 }
 #if ($debug) { 
 #  print STDERR "  => resulting list: ", 
 #    join(',', (map { $$_[0].'-'.$$_[1] } @$reglist)),"\n";
 #}
}

# --------------------- 
sub getOvl { 
  #returns coordinates of an overlap region between e1-e2 and o1-o2
  my ($e1, $e2, $o1, $o2)=@_;
  #my ($e1, $e2)=($r->[0], $r->[1]);
  if ($o1 <= $e2 && $e1 <= $o2) {
     return ( ($o1<$e1 ? $e1 : $o1), ($o2>$e2 ? $e2 : $o2) );
  }
  else { return (0,0); }
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

