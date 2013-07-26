#!/usr/bin/perl

use DBI;
use FindBin;
use lib "$FindBin::Bin/../lib/perl";
use OrthoMCLEngine::Main::Base;
use strict;


usage() unless (@ARGV >= 1);
my $configFile = $ARGV[0];
my $sqlLog = $ARGV[1];
my $suffix = $ARGV[2];

my $base = OrthoMCLEngine::Main::Base->new($configFile);
my $dbh = $base->getDbh();

if ($sqlLog) {
  open (LOGFILE, ">$sqlLog");
}


runSql("drop table " . $base->getConfig("similarSequencesTable") . $suffix);
runSql("drop table " . $base->getConfig("inParalogTable") . $suffix);
runSql("drop table " . $base->getConfig("orthologTable") . $suffix);
runSql("drop table " . $base->getConfig("coOrthologTable") . $suffix);
runSql("drop view " . $base->getConfig("interTaxonMatchView") . $suffix);

##############################################################


sub runSql {
 my $sql = $_[0];
 if ($sqlLog) {
     logSql($sql);
    }
  my $stmt = $dbh->prepare($sql) or die DBI::errstr;
      $stmt->execute() or die DBI::errstr;
}


sub logSql {
  my $sql = $_[0];
  print LOGFILE "\n$sql";
}

sub usage {
 print STDERR "
Drop OrthoMCL schema

usage: orthomclDropSchema config_file sql_log_file

where:
  config_file : see below
  sql_log_file : optional log of sql executed

EXAMPLE: orthomclDropSchema my_orthomcl_dir/orthomcl.config my_orthomcl_dir/drop_schema.log

NOTE: the database login in the config file must have the required database privileges on the tables specified in the config file.

Sample Config File:

dbVendor=oracle  (or mysql)
dbConnectString=dbi:Oracle:orthomcl
dbLogin=my_db_login
dbPassword=my_db_password
oracleIndexTablespace=
similarSequencesTable=SimilarSequences
orthologTable=Ortholog
inParalogTable=InParalog
coOrthologTable=CoOrtholog
interTaxonMatchView=InterTaxonMatch

";
 exit(1);
}

