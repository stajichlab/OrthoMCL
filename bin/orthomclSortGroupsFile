#!/usr/bin/perl

use strict;

&usage() unless scalar(@ARGV) == 1;

my $groupsFile = $ARGV[0];

open(F, $groupsFile) || die "Can't open groups file '$groupsFile'\n";

my %groups;
while(<F>) {
  chomp;
  my @a = split(/\s/);
  my $groupId = shift(@a);
  my @b = sort(@a);
  my $sortedMembers = join(" ", @b);
  $groups{$sortedMembers} = $groupId;
}

foreach my $grp (sort(keys(%groups))) {
  print "$groups{$grp} $grp\n"
}

sub usage {
    die "

Produce a groups file in which the members of each group are sorted.  Useful
for comparing groups files.

Output to standard out.

NOTE: output files is not sorted, only the members of each group are.  Pass to sort for final sorting.

Usage:  orthomclSortGroupsFile groups.txt | sort > groups_sorted.txt

";
}

