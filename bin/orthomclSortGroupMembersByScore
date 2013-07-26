#!/usr/bin/perl

use strict;

&usage() unless scalar(@ARGV) == 2;

my $groupsFile = shift(@ARGV);
my $pairsDir = shift(@ARGV);

my $orthologPairsFile = "$pairsDir/orthologs.txt";
my $paralogPairsFile = "$pairsDir/inparalogs.txt";
my $coorthologPairsFile = "$pairsDir/coorthologs.txt";

die "groups file '$groupsFile' does not exist" unless -e $groupsFile;
die "ortholog pairs file '$orthologPairsFile' does not exist" unless -e $orthologPairsFile;
die "paralog pairs file '$paralogPairsFile' does not exist" unless -e $paralogPairsFile;
die "coortholog pairs file '$coorthologPairsFile' does not exist" unless -e $coorthologPairsFile;

# pass through groups, making hashes both ways
my %protein2groupId;
my %groupId2proteins;
print STDERR "Scanning groups file\n";
open(INGROUPS, $groupsFile)  || die "Can't open groups file '$groupsFile' for reading";
while (<INGROUPS>) {
  chomp;
  my @group = split(/\s+/, $_);
  my $groupId = shift @group;
  $groupId2proteins{$groupId} = \@group;
  foreach my $protein (@group) {
      $protein2groupId{$protein} = $groupId;
  }
}
close(INGROUPS);

# iterate through score pairs (from ortholog, paralog and coortholog files)
# for each pair whose proteins are in the same group, add the score to each proteins' sum
print STDERR "Scanning pairs files\n";
my %protein2scoreSum;  # protein --> its score sum
addScoresToProteins(\%protein2scoreSum, $orthologPairsFile, \%protein2groupId);
addScoresToProteins(\%protein2scoreSum, $paralogPairsFile, \%protein2groupId);
addScoresToProteins(\%protein2scoreSum, $coorthologPairsFile, \%protein2groupId);

# iterate through groups.  for each, sort by score
print STDERR "Sorting groups\n";
foreach my $groupId (keys(%groupId2proteins)) {
    my @sortedGroup = sort { $protein2scoreSum{$b} <=> $protein2scoreSum{$a} } @{$groupId2proteins{$groupId}};
    print join(" ", $groupId, @sortedGroup) . "\n";
}

#######################################################################################

# iterate through score pairs (from pairs file)
# for each pair whose proteins are in the same group, add the score to each proteins' sum
sub addScoresToProteins {
  my ($protein2scoreSum, $pairsFile, $protein2group) = @_;

  open(PAIRS, $pairsFile) || die "can't open pairs file '$pairsFile'";
  while(<PAIRS>) {
    my @a = split(/\s+/);
    my $protein1 = $a[0];
    my $protein2 = $a[1];
    my $score = $a[2];
    if ($protein2group->{$protein1} eq $protein2group->{$protein2}) {
      $protein2scoreSum->{$protein1} += $score;
      $protein2scoreSum->{$protein2} += $score;
    }
  }
  close(PAIRS);
}

sub usage {
  die "
Sort the groups in an OrthoMCL groups file so that each group has its members sorted by score, with the highest first.  (The order of the groups in the file is unchanged.) 

The score for each protein is made by summing its pairwise scores with the other members of its group.  The pairwise scores a provided in the orthologs.txt, paralogs.txt and coorthologs.txt files provided in the input direcgtory.

Usage: orthomclSortGroupMembersByScore groups_file pairs_dir

Input:
 - groups_file: standard orthomcl groups file
 - pairs_dir: a directory containing these files:
    orthologs.txt
    paralogs.txt
    coorthologs.txt

Output:
 - standard out:  sorted groups file

";

}
