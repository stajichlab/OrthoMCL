#!/usr/bin/perl
use strict;

usage() unless scalar(@ARGV) == 1;

open(F, $ARGV[0]) || die "Could not open file '$ARGV[0]'\n";

while(<F>) {
  chomp;

  my @a = split;
  shift @a;
  foreach my $member1 (@a) {
    foreach my $member2 (@a) {
      next unless $member1 lt $member2;
      print "$member1 $member2\n";
    }
  }
}

sub usage {

die "
Extract protein ID pairss from an orthomcl groups file.

Usage orthomclExtractProteinIdsFromGroupsFile groups_file

Output to standard out.
";

}
