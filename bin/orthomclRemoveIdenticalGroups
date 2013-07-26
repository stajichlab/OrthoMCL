#!/usr/bin/perl

use strict;

&usage() unless scalar(@ARGV) == 2;

my $groupsFile = $ARGV[0];
my $refGroupsFile = $ARGV[1];

open(R, $refGroupsFile) || die "Can't open ref groups file '$refGroupsFile'\n";

my %refs;
while(<R>) {
  chomp;
  my @a = split(/\s/);
  my $groupId = shift(@a);
  my @b = sort(@a);
  my $sortedMembers = join(" ", @b);
  $refs{$sortedMembers} = 1;
}
close(R);

open(G, $groupsFile) || die "Can't open groups file '$groupsFile'\n";
while(<G>) {
  my $grp = $_;
  chomp;
  my @a = split(/\s/);
  shift(@a);
  my @b = sort(@a);
  my $sortedMembers = join(" ", @b);
  print "$grp" unless $refs{$sortedMembers};
}

sub usage {
    die "

Produce a groups file that has had removed from it those groups that are 
identical in a second, reference, groups file.   

The group IDs are ignored by the filter.   Only the members of groups are
considered.

The input files must be sorted (see orthomclSortGroupsFile).
Usage:  orthomclRemoveIdenticalGroups groups.txt ref_groups.txt

Output to standard out.

";
}

