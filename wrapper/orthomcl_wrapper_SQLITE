#!/bin/bash -l

source /usr/local/Modules/3.2.9/init/bash

module load ncbi-blast

SCRIPTS="$( cd "$(dirname "$0")" ; pwd -P )"
ORTHOMCL_BIN=$ORTHO ##set environmental variable in .bashrc: export ORTHO=bin_dir_for_orthoMCL: source .bashrc
OPTS=$(getopt -o p:u:s:o: -l "protein_dir:,uniqid:,split:,ortho_db:" -n orthoWrapper.sh -- "$@")
if [ $? -ne 0 ]
then
    exit 1
fi

eval set -- "$OPTS"

while true ; do
    case "$1" in
        -p|--protein_dir) 
          shift; 
          if [ -n "$1" ]; then
            echo "setting PROTEIN_FASTA_DIR to be $1";
            PROTEIN_FASTA_DIR=`realpath $1`;
            shift;
          fi
          ;;
        -u|--uniqid) 
          shift; 
          if [ -n "$1" ]; then
            echo "setting UNIQID to be $1";
            UNIQID=$1;
            shift;
          fi
          ;;
        -o|--ortho_db) 
          shift;      
          if [ -n "$1" ]; then
            echo "setting ORTHO_DB to be $1";
            ORTHO_DB=$1;
            shift;
          fi
          ;;
        -s|--split) 
          shift;      
          if [ -n "$1" ]; then
            echo "setting SPLIT to be $1";
            SPLIT=$1;
            shift;
          fi
          ;;
         --)
           shift;
           break;
           ;;

    esac
done

if [ ! $PROTEIN_FASTA_DIR  ] ; then
  echo "Need to provide protein fasta dir"
  echo "$0 -p [protein dir] -u [fasta header uniq id col] -o [ortho db name]"
  echo "   
protein_dir  :     a directory of FASTA files. One FASTA for each species.
                   Each FSATA file must have a name in the form 'xxxx.fasta' 
                   where xxxx is a three or four letter unique taxon code.  
                   For example: hsa.fasta or eco.fasta
                   [no default]
db_name      :     input a name for the new database
                   [default=orthoMCL.DB.`date +%Y.%m.%d`]
uniq_id      :     a number indicating what field in the definition line contains
                   the protein ID.  Fields are separated by either ' ' or '|'. Any
                   spaces immediately following the '>' are ignored.  The first
                   field is 1. For example, in the following definition line, the
                   ID (AP_000668.1) is in field 4:  >gi|89106888|ref|AP_000668.1|
                   [default=2]
split        :     How many blast jobs do you want to run? 
		   [default=500]"
  exit 1

fi
if [ ! $UNIQID  ] ; then
  ## guess
  UNIQID=2
fi
if [ ! $ORTHO_DB  ] ; then
  ORTHO_DB="orthoMCL.DB.`date +%Y.%m.%d`"; 
fi

if [ ! $SPLIT  ] ; then 
  SPLIT=500
fi

echo "PROTEIN_FASTA_DIR is $PROTEIN_FASTA_DIR"
echo "UNIQID id $UNIQID"
echo "ORTHO_DB is $ORTHO_DB"
echo "SPLIT is $SPLIT"


DIR=`pwd`
cd $DIR

if [ ! -e orthomcl.config ] ; then
  echo "dbVendor=sqlite 
dbConnectString=DBI:SQLite:database=${ORTHO_DB}.sqlite
dbLogin=none
dbPassword=none
similarSequencesTable=SimilarSequences
orthologTable=Ortholog
inParalogTable=InParalog
coOrthologTable=CoOrtholog
interTaxonMatchView=InterTaxonMatch
percentMatchCutoff=50
evalueExponentCutoff=-5
oracleIndexTblSpc=NONE" > orthomcl.config
fi

# drop db to get a clean slate
if [ -e "$ORTHO_DB.sqlite" ] ; then
  echo "Deleting $ORTHO_DB.sqlite"
  rm $ORTHO_DB.sqlite
fi


$ORTHOMCL_BIN/orthomclInstallSchema orthomcl.config



#make a copy of the original FASTA into dir called proteomes
if [ ! -d proteomes ]; then
 mkdir proteomes
 cd proteomes
 for file in $PROTEIN_FASTA_DIR/*.fasta
 do
  echo $file
  base=`basename $file .fasta`
  cp $file $base.pep 
 done
 cd ..
fi

# let orthoMCL clean the FASTAs
if [ ! -d cleanseq ] ; then
  mkdir cleanseq
  cd cleanseq
  for file in $DIR/proteomes/*.pep
   do
   base=`basename $file .pep`
   ## which field has a uniq ID? => 2
   $ORTHOMCL_BIN/orthomclAdjustFasta $base $file $UNIQID 
  done
  cd $DIR
  $ORTHOMCL_BIN/orthomclFilterFasta cleanseq 10 10
fi

## format goodProteins BLAST DB
if [ ! -e $DIR/goodProteins.fasta.pin ] ; then
  makeblastdb -in $DIR/goodProteins.fasta -title "OrthoMCLPeps" -dbtype prot
fi

## split up goodProteins into smaller FASTA
MAX=`expr $SPLIT - 1`
if [ ! -d $DIR/good_proteins_split ] ; then
  mkdir $DIR/good_proteins_split
fi
if [ ! -e $DIR/good_proteins_split/goodProteins_part_${MAX}.fasta ] ; then
  cd $DIR/good_proteins_split
  cp -s $DIR/goodProteins.fasta .
  $SCRIPTS/split_fasta.pl $SPLIT goodProteins.fasta
  rm -f $DIR/good_proteins_split/goodProteins.fasta
fi 
cd $DIR 

## get ready for blast
if [ ! -d $DIR/blastp_out ] ; then
  mkdir $DIR/blastp_out
fi
BLAST_Q=0
## consider adding check for qsub:  if [ `which qsub` != '' ] ; then
if [ ! -d $DIR/blastp_out ] ; then
  ## make array BLAST job
  BLAST_Q=1
  echo "#PBS -t 0-${MAX}
#PBS -q js
##PBS -l nodes=1:ppn=4
module load ncbi-blast
cd $DIR
blastp -query good_proteins_split/goodProteins_part_\${PBS_ARRAYID}.fasta -db goodProteins.fasta -num_threads 4 -outfmt 6 -out $DIR/blastp_out/goodProteins_part_\${PBS_ARRAYID}_vs_goodProteins.BLASTP.tab -evalue 1e-3" > $DIR/run_split_blast_array.sh


  #run array BLAST job
  BLASTJOB=`qsub run_split_blast_array.sh`
  echo "JOBS SUBMITTED AND HOPEFULLY RUNNING ON YOUR QUEUE:"
  echo "$BLASTJOB = qsub $DIR/run_split_blast_array.sh"
fi


echo "
#module load orthoMCL
#module load ncbi-blast

cd $DIR 

# This will reformat the BLAST data into something to load into the database
if [ ! -f goodProteins.BLASTP.bpo ]; then
 cat blastp_out/*tab > goodProteins.BLASTP.tab
 $ORTHOMCL_BIN/orthomclBlastParser goodProteins.BLASTP.tab cleanseq > goodProteins.BLASTP.bpo
fi

# Load the data into the DB
$ORTHOMCL_BIN/orthomclLoadBlast orthomcl.config goodProteins.BLASTP.bpo

# now run the ortholog/paralog initial finding
rm -rf pairs pairs.log 
$ORTHOMCL_BIN/orthomclPairs orthomcl.config  pairs.log cleanup=no

# Dump out the ortholog groups and the mclInput file
$ORTHOMCL_BIN/orthomclDumpPairsFiles orthomcl.config

# Run mcl for clustering
mcl mclInput  --abc -I 1.5 -o mclOutput.I15.out

# convert the MCL clusters into OrthoMCL groups
$ORTHOMCL_BIN/orthomclMclToGroups OG 1 < mclOutput.I15.out > mclGroups.I15.table
" > $DIR/finish.orthomcl.sh

if [ $BLAST_Q -eq 1 ] ; then
  qsub -W depend=afterokarray:$BLASTJOB $DIR/finish.orthomcl.sh
  echo "qsub -W depend=afterokarray:$BLASTJOB $DIR/finish.orthomcl.sh" 
else
  sh $DIR/finish.orthomcl.sh
fi
