#!/usr/bin/perl

use strict;

usage() unless (scalar(@ARGV) == 2 || scalar(@ARGV) == 3);
my $proteinsFile = $ARGV[0];
my $groupsFile = $ARGV[1];
my $idPrefix = $ARGV[2];

my $startIdNumber;

open(F, $groupsFile) || die "Can't open groups file '$groupsFile'\n";
my $groups;
while (<F>) {
    my @a = split(/\s+/);

    if ($idPrefix) {
      if ($a[0] =~ /${idPrefix}([0-9]+):/ ) {
	$startIdNumber = $1 if $1 > $startIdNumber;
      } else {
	die "group ID \"" . $a[0] . "\" does not begin with the supplied prefix \"$idPrefix\"";
      }
    }

    shift(@a);
    map {$groups->{$_} = 1} @a;
}
close(F);

$startIdNumber++ if $startIdNumber;

open(F, $proteinsFile) || die "Can't open proteins file '$proteinsFile'\n";
while(<F>) {
    next unless /\>\s*(\S+)/;
    next if $groups->{$1};
    my $id;
    if ($startIdNumber) {
	$id = "$idPrefix$startIdNumber: ";
	$startIdNumber++;
    }
    print "$id$1\n";
}

sub usage {
 print STDERR "
Find proteins that are not in the groups.txt file.  Output those protein IDs to standard out.

usage: orthomclSingletons good_proteins groups_file [id_prefix]

If optional id_prefix is provided, start each line with a generated group
id that begins with the start_id_number, in the format of a standard orthomcl groups file, like this:
  OG_124445: some_protein_id

EXAMPLE: orthomclSoftware/bin/orthomclSingletons my_orthomcl_dir/goodProteins.fasta my_orthomcl_dir/groups.txt > my_orthomcl_dir/singletons.txt

";
 exit(1);
}
