#!/usr/bin/perl

use strict;

&usage() unless scalar(@ARGV) == 2;

my $groupsfile = $ARGV[0];
my $taxa_str = $ARGV[1];

if ($groupsfile =~ /\.gz$/) {
  open(IN, "zcat $groupsfile|") || die "Can't open groups file '$groupsfile'\n";
} else {
  open(IN, $groupsfile) || die "Can't open groups file '$groupsfile'\n";
}
my @t = split(/,/, $taxa_str);
scalar(@t) > 1 || die "Invalid taxa specification '$taxa_str'\n";

my $taxa;
map {$taxa->{$_} = 1} @t;


while(<IN>) {
    chomp;
    my @group = split(/\s/);
    my $group_id = shift(@group);  
    my $filteredGroup;
    foreach my $id (@group) {
	my ($taxon, $dontcare) = split(/\|/, $id);
	push(@$filteredGroup, $id) if $taxa->{$taxon};
    }
    if ($filteredGroup) {
	my $fg = join(" ", @$filteredGroup);
	print "$group_id $fg\n";
    }
}



sub usage {
print STDERR "
Reduce a groups file by taxon.  Input is a groups file and a set of taxa.
Output is a groups file that contains only those taxa

Usage:
  orthomclReduceGroups groups_file taxa

where:
  groups_file: a standard orthomcl groups file. (.gz file is ok)
  taxa:        a comma delimited list of taxon abbreviations

EXAMPLE: orthomclSoftware/bin/orthomclReduceGroups groups.txt hsa,pfa,txo

";
exit(1);
}
