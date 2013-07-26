#!/usr/bin/perl

use strict;

&usage() unless scalar(@ARGV) == 2;

my $fastafile = $ARGV[0];
my $taxa_str = $ARGV[1];

if ($fastafile =~ /\.gz$/) {
  open(IN, "zcat $fastafile|") || die "Can't open fasta file '$fastafile'\n";
} else {
  open(IN, $fastafile) || die "Can't open fasta file '$fastafile'\n";
}
my @t = split(/,/, $taxa_str);
scalar(@t) > 1 || die "Invalid taxa specification '$taxa_str'\n";

my $taxa;
map {$taxa->{$_} = 1} @t;

my $currentSeq;
my $currentTaxon;
# process lines of one file
while (<IN>) {
  chomp;
  # handle prev seq
  if (/\>([^\|]+)/) {
    print $currentSeq if ($currentSeq && $taxa->{$currentTaxon});
    $currentTaxon = $1;
    $currentSeq = "";
  }
  $currentSeq .= "$_\n" if $taxa->{$currentTaxon};
}
print $currentSeq if $taxa->{$currentTaxon};

sub usage {
print STDERR "
Reduce a fasta file by taxon.  Input is a fasta file and a set of taxa.
Output is a fasta file that contains only those taxa

Usage:
  orthomclReduceFasta fasta_file taxa

where:
  fasta_file: a standard orthomcl compatible fasta file. (.gz file is ok)
  taxa:        a comma delimited list of taxon abbreviations

EXAMPLE: orthomclSoftware/bin/orthomclReduceFasta proteins.fasta hsa,pfa,txo

";
exit(1);
}
