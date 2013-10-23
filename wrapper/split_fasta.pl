#!/usr/bin/perl 

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use warnings;
use strict "vars";
use strict "refs";

use FindBin;
use Bio::SeqIO; 
use lib "$FindBin::Bin/../lib";
use vars qw($LOG $CMD_ARGS);


my $usage = "
Usage:

        split_fasta [count] <fasta_input>

        This script splits a multi-fasta file into the number of files specified by count.

";
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
if (@ARGV != 2){
    print $usage;
    exit;
}

#variables that are persistent outside of try block
my $count = shift @ARGV;
die "ERROR: the count nust be set to a value greater than 1.\n"if ($count <= 1);
my $infile = shift @ARGV;
die "ERROR: The file \'$infile\' does not exist.\n" if(! -e $infile);

my $seqIO_obj = Bio::SeqIO->new(-file=>$infile, -format=>'fasta');
my $total_seq_count = `grep -c '>' $infile`;
chomp $total_seq_count;
my $average_seqs_per_file =  (  $total_seq_count / $count );
my $file_count_with_more = $count;
#my $file_count_with_less = 0;
my $seqs_per_file = $average_seqs_per_file;
if ($average_seqs_per_file > int $average_seqs_per_file){
  $seqs_per_file = (int $average_seqs_per_file) + 1; 
  $file_count_with_more =  $count * ($average_seqs_per_file - int $average_seqs_per_file);
  #$file_count_with_less = $count - $file_count_with_more;
}
print "$infile: newFileCount=$count; totalSeqCount=$total_seq_count; seqsPerFile=$average_seqs_per_file\n";

#for each new file
for (my $i = 0 ; $i < $count ; $i++){
  #print $seqs_per_file num of seqs
  my $outfile = $infile;
  $outfile =~ s/\.fasta$//;
  $outfile .= "_part_$i.fasta";
  open OUT, ">$outfile";
  if ($i >= $file_count_with_more ){
    $seqs_per_file = int $average_seqs_per_file;
  }
  for (my $j=0 ; $j < $seqs_per_file ; $j++){
    my $seq_obj = $seqIO_obj->next_seq();
    last if !defined $seq_obj;
    my $id = $seq_obj->id;
    my $seq = $seq_obj->seq;
    my $desc = $seq_obj->desc;
    #$seq =~ s/(.{1,80})/$1\n/g;
    print OUT ">$id $desc\n$seq\n";
  }
  close OUT;
}
sub round_up {
  my $num = shift;
  my $int = int $num;
  if ($num > $int ){
    $num = ++$int;
  }
  return $num;
}
