#!/usr/bin/perl

use strict;

&usage() unless scalar(@ARGV) == 3;

my $taxoncode = $ARGV[0];
my $inputfile = $ARGV[1];
my $idField = $ARGV[2];

open(IN, $inputfile) || die "Can't open input file '$inputfile'\n";
open(OUT, ">$taxoncode.fasta") || die "Can't open output file '$taxoncode.fasta'\n";

my %ids;
while(<IN>) {
  if (/\>/) {
    s/^\>\s*//;
    s/\s+/ /g;
    s/\s*\|\s*/\|/g;
    my @a = split(/[\s\|]/);
    my $id = $a[$idField-1];
    die "Fasta file '$inputfile' contains a duplicate id: $id\n" if $ids{$id};
    $ids{$id} = 1;
    print OUT ">$taxoncode|$id\n";
  } else {
    print OUT $_;
  }
}



sub usage {
print STDERR "
Create an OrthoMCL compliant .fasta file, by adjusting definition lines.

Usage:
  orthomclAdjustFasta taxon_code fasta_file id_field

where:
  taxon_code:  a three or four letter unique abbreviation for the taxon
  fasta_file:  the input fasta file
  id_field:    a number indicating what field in the definition line contains
               the protein ID.  Fields are separated by either ' ' or '|'. Any
               spaces immediately following the '>' are ignored.  The first
               field is 1. For example, in the following definition line, the
               ID (AP_000668.1) is in field 4:  >gi|89106888|ref|AP_000668.1|

Input file requirements:
  (1) .fasta format
  (2) a unique id is provided for each sequence, and is in the field specified
      by id_field

Output file format:
  (1) .fasta format
  (2) definition line is of the form:
         >taxoncode|unique_protein_id

The output file is named taxoncode.fasta

Note: if your input files do not meet the requirements, you can do some simple perl or awk processing of them to create the required input files to this program, or the required output files.  This program is provided as a convenience, but OrthoMCL users are expected to have the scripting skills to provide OrthoMCL compliant .fasta files.

EXAMPLE: orthomclSoftware/bin/orthomclAdjustFasta hsa Homo_sapiens.NCBI36.53.pep.all.fa 1

";
exit(1);
}
