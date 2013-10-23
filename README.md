Authors:
Sofia Robb sofia.robb-at-ucr.edu

Overview:
Modifying OrthoMCL to run with SQLite in addition to the current functionality with MySQL and Oracle.

Progress: (Overview of step from http://orthomcl.org/common/downloads/software/v2.0/UserGuide.txt)
===================================================
=========== Overview of steps =====================
===================================================

This is an overview of the thirteen steps to run orthomcl.  Details are in the next sections. 


WORKS with SQLITE3 (4) run orthomclInstallSchema to install the required schema into the database.

WORKS with SQLITE3 (5) run orthomclAdjustFasta (or your own simple script) to generate protein fasta files in the required format.

WORKS with SQLITE3 (6) run orthomclFilterFasta to filter away poor quality proteins, and optionally remove alternative proteins. Creates a single large goodProteins.fasta file (and a poorProteins.fasta file)

WORKS with SQLITE3 (7) run all-v-all NCBI BLAST on goodProteins.fasta (output format is tab delimited text).

WORKS with SQLITE3 (8) run orthomclBlastParser on the NCBI BLAST tab output to create a file of similarities in the required format

WORKS with SQLITE3 (9) run orthomclLoadBlast to load the output of orthomclBlastParser into the database.

WORKS with SQLITE3 (10) run the orthomclPairs program to compute pairwise relationships. 

WORKS with SQLITE3 (11) run the orthomclDumpPairsFiles program to dump the pairs/ directory from the database

WORKS with SQLITE3 (12) run the mcl program on the mcl_input.txt file created in Step 11.

WORKS with SQLITE3 (13) run orthomclMclToGroups to convert mcl output to groups.txt
