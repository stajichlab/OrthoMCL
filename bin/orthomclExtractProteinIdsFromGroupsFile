#!/usr/bin/perl
use strict;

usage() unless scalar(@ARGV) == 1;

open(F, $ARGV[0]) || die "Could not open file '$ARGV[0]'\n";

while(<F>) {
  chomp;

  my @a = split;
  shift @a;
  print join("\n", @a) . "\n";
}

sub usage {

die "
Extract protein IDs from an orthomcl groups file.

Usage orthomclExtractProteinIdsFromGroupsFile groups_file

Output to standard out.
";

}
