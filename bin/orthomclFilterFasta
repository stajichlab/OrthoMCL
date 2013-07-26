#!/usr/bin/perl

use strict;

&usage() unless scalar(@ARGV) == 3 || scalar(@ARGV) == 5;

my $inputDir = $ARGV[0];
my $minLength = $ARGV[1];
my $maxStopPercent = $ARGV[2];
my $goodProteinsFile = $ARGV[3]? $ARGV[3] : 'goodProteins.fasta' ;
my $poorProteinsFile = $ARGV[4]? $ARGV[4] : 'poorProteins.fasta';

opendir(DIR, $inputDir) || die "Can't open input directory '$inputDir'\n";
my @files = readdir(DIR);
closedir(DIR);

die "Input directory $inputDir does not contain any files" unless scalar(@files);

my $rejectRates = [];
open(GOOD, ">$goodProteinsFile");
open(BAD, ">$poorProteinsFile");
foreach my $file (@files) {
  next unless $file =~ /\.fasta$/;
  $file =~ /(\w+)\.fasta$/ || die "File '$file' does not have a name in xxxx.fasta format\n";;
  my $abbrev = $1;
  open(F, "$inputDir/$file") || die "Can't open input file '$file'\n";
  print STDERR "processing file $file\n";
  my $seqCount;
  my $rejectSeqCount;
  my $currentSeq;
  my $currentLen;
  my $currentStopCnt;

  # process lines of one file
  while (<F>) {
    chomp;
    # handle prev seq
    if (/\>/) {
      (/\>([^|]+)|/ && $1 eq $abbrev) ||
	die "The ID on def line '$_' is missing the prefix '$abbrev|' '$1'\n" unless $1 eq $abbrev;
      if ($currentSeq) {
	die "Error: zero length sequence in file $file.  Look near line '$_'\n" if $currentLen == 0;
	$seqCount++;
	$rejectSeqCount += &handleSeq($currentSeq, $currentLen, $currentStopCnt);
	$currentSeq = "";
	$currentLen = 0;
	$currentStopCnt = 0;
      }
    } else {
      $currentLen += length($_);
      $currentStopCnt += tr/[^A-Za-z]//; # this removes the stop codon from $_
    }
    $currentSeq .= "$_\n";
  }
  $rejectSeqCount += &handleSeq($currentSeq, $currentLen, $currentStopCnt);
  $seqCount++;

  # add file stats to reject count if it qualifies
  if ($rejectSeqCount) {
    my $pct = $rejectSeqCount/$seqCount * 100;
    if ($pct > 10) {
      push(@$rejectRates, [$file, $pct]);
    }
  }
  close(F);
}

if (scalar(@$rejectRates)) {
  print "\nProteomes with > 10% poor proteins:\n";
  my @sortedRR = sort {$b->[1] <=> $a->[1]} @$rejectRates;
  foreach my $reject (@sortedRR) {
    my $intPct = int($reject->[1]);
    print "  $reject->[0]\t$intPct%\n";
  }
}

sub handleSeq {
  my ($seq, $len, $stopCnt) = @_;
  my $isBad = 0;
  my $stopPercent = (($len - $stopCnt)/$len)* 100;

  if ($len < $minLength || $stopPercent > $maxStopPercent) {
    print BAD $seq;
    $isBad = 1;
  } else {
    print GOOD $seq;
  }
  return $isBad;
}

sub usage {
  print STDERR "
Create goodProteins.fasta containing all good proteins and rejectProteins.fasta containing all rejects.  Input is a directory containing a set of compliant input .fasta files (as produced by orthomclAdjustFasta).

Usage:
  orthomclFilterFasta input_dir min_length max_percent_stops [good_proteins_file poor_proteins_file]

where:
  input_dir:               a directory containing a set of .fasta files
  min_length:              minimum allowed length of proteins.  (suggested: 10)
  max_percent_stop:        maximum percent stop codons.  (suggested 20)
  good_proteins_file:      optional.  By default goodProteins.fasta in the current dir.
  poor_proteins_file:      optional.  By default poorProteins.fasta in the current dir.

EXAMPLE: orthomclSoftware/bin/orthomclFilterFasta my_orthomcl_dir/compliantFasta 10 20

";
  exit(1);
}
