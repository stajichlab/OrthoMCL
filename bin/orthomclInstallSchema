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

my $dbVendor;
my $intType = ($dbVendor eq 'oracle') ? 'NUMBER' : 'INT';

my $oracleNoLogging;
my $oracleIndexTablespace;

if ($base->getConfig("dbVendor") eq 'oracle') {
  $oracleNoLogging = " NOLOGGING ";
  $oracleIndexTablespace = $base->getConfig("oracleIndexTablespace");
  $oracleIndexTablespace = " TABLESPACE " . $oracleIndexTablespace
    if $oracleIndexTablespace;
}

createSimilarSequencesTable();
createInParalogTable();
createOrthologTable();
createCoOrthologTable();
createInterTaxonMatchView();

##############################################################

sub createSimilarSequencesTable {

  my $sst = $base->getConfig("similarSequencesTable");

  my $sql = "
CREATE TABLE $sst$suffix (
 QUERY_ID                 VARCHAR(60),
 SUBJECT_ID               VARCHAR(60),
 QUERY_TAXON_ID           VARCHAR(40),
 SUBJECT_TAXON_ID         VARCHAR(40),
 EVALUE_MANT              FLOAT,
 EVALUE_EXP               $intType,
 PERCENT_IDENTITY         FLOAT,
 PERCENT_MATCH            FLOAT
) $oracleNoLogging
";
  runSql($sql);

  $sql = "
CREATE INDEX ss_qtaxexp_ix$suffix
ON $sst$suffix(query_id, subject_taxon_id,
evalue_exp, evalue_mant,
query_taxon_id, subject_id) $oracleNoLogging $oracleIndexTablespace
";
  runSql($sql);

  $sql = "
CREATE INDEX ss_seqs_ix$suffix
ON $sst$suffix(query_id, subject_id,
evalue_exp, evalue_mant, percent_match) $oracleNoLogging $oracleIndexTablespace
";
  runSql($sql);
}


sub createInParalogTable {

 my $ipt = $base->getConfig("inParalogTable");
 my $sql = "
CREATE TABLE $ipt$suffix (
 SEQUENCE_ID_A           VARCHAR(60),
 SEQUENCE_ID_B           VARCHAR(60),
 TAXON_ID                VARCHAR(40),
 UNNORMALIZED_SCORE      FLOAT,
 NORMALIZED_SCORE        FLOAT
)
";
  runSql($sql);

$sql = "
CREATE INDEX inparalog_seqa_ix$suffix
ON $ipt$suffix(sequence_id_a) $oracleNoLogging $oracleIndexTablespace
";
  runSql($sql);


$sql = "
CREATE INDEX inparalog_seqb_ix$suffix
ON $ipt$suffix(sequence_id_b) $oracleNoLogging $oracleIndexTablespace
";
  runSql($sql);

}


sub createOrthologTable  {

 my $olt = $base->getConfig("orthologTable");
 my $sql = " 
CREATE TABLE $olt$suffix (
 SEQUENCE_ID_A           VARCHAR(60),
 SEQUENCE_ID_B           VARCHAR(60),
 TAXON_ID_A              VARCHAR(40),
 TAXON_ID_B              VARCHAR(40),
 UNNORMALIZED_SCORE      FLOAT,
 NORMALIZED_SCORE        FLOAT
)
";
  runSql($sql);

$sql = "
CREATE INDEX ortholog_seq_a_ix$suffix
ON $olt$suffix(sequence_id_a) $oracleNoLogging $oracleIndexTablespace
";
  runSql($sql);

$sql = "
CREATE INDEX ortholog_seq_b_ix$suffix
ON $olt$suffix (sequence_id_b) $oracleNoLogging $oracleIndexTablespace
";
  runSql($sql);
}


sub createCoOrthologTable  {

 my $cot = $base->getConfig("coOrthologTable");
 my $sql = "
CREATE TABLE $cot$suffix (
 SEQUENCE_ID_A           VARCHAR(60),
 SEQUENCE_ID_B           VARCHAR(60),
 TAXON_ID_A              VARCHAR(40),
 TAXON_ID_B              VARCHAR(40),
 UNNORMALIZED_SCORE      FLOAT,
 NORMALIZED_SCORE        FLOAT
)
";
  runSql($sql);
$sql = "
CREATE INDEX coortholog_seq_a_ix$suffix
ON $cot$suffix(sequence_id_a) $oracleNoLogging $oracleIndexTablespace
";
  runSql($sql);

$sql = "
CREATE INDEX coortholog_seq_b_ix$suffix
ON $cot$suffix (sequence_id_b) $oracleNoLogging $oracleIndexTablespace
";
  runSql($sql);
}

sub createInterTaxonMatchView {

 my $sst = $base->getConfig("similarSequencesTable");
 my $itv = $base->getConfig("interTaxonMatchView");
 my $sql = "
CREATE OR REPLACE VIEW $itv$suffix
	AS SELECT ss.query_id, ss.subject_id, ss.subject_taxon_id,
	ss.evalue_mant, ss.evalue_exp
	FROM $sst$suffix ss
	WHERE ss.subject_taxon_id != ss.query_taxon_id
";

## added by sofia
 if( $base->getConfig("dbVendor") eq 'sqlite'){
   $sql = "
CREATE VIEW $itv$suffix
	AS SELECT ss.query_id, ss.subject_id, ss.subject_taxon_id,
	ss.evalue_mant, ss.evalue_exp
	FROM $sst$suffix ss
	WHERE ss.subject_taxon_id != ss.query_taxon_id
";
  } 
## end add by sofia
 runSql($sql);
}

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
Create OrthoMCL schema in an Oracle or Mysql database.

usage: orthomclInstallSchema config_file sql_log_file table_suffix

where:
  config_file : see below
  sql_log_file : optional log of sql executed
  table_suffix : optional string to append to database object names

EXAMPLE: orthomclSoftware/bin/orthomclInstallSchema my_orthomcl_dir/orthomcl.config my_orthomcl_dir/install_schema.log

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

