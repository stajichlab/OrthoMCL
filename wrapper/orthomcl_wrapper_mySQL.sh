#!/bin/bash -l

source /usr/local/Modules/3.2.9/init/bash

module load orthoMCL
module load ncbi-blast

#"database=orthomcl_andrii;host=puffball.fungalgenomes.org"
SCRIPTS=./
ORTHOMCL_BIN=$ORTHO ##set environmental variable in .bashrc: export ORTHO=bin_dir_for_orthoMCL: source .bashrc

OPTS=$(getopt -o p:u:s:o:d:h:x:n: -l "protein_dir:,uniqid:,split:,db_name:,db_host:,db_psswrd:,db_name" -n $0 -- "$@")
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
        -s|--split) 
          shift;      
          if [ -n "$1" ]; then
            echo "setting SPLIT to be $1";
            SPLIT=$1;
            shift;
          fi
          ;;
        -d|--db_name)
          shift;
          if [ -n "$1" ]; then
            echo "setting DB_NAME to be $1";
            DB_NAME=$1;
            shift;
          fi
          ;;
         -h|--db_host)
          shift;
          if [ -n "$1" ]; then
            echo "setting DB_HOST to be $1";
            DB_HOST=$1;
            shift;
          fi
          ;;
          -x|--db_psswrd)
          shift;
          if [ -n "$1" ]; then
            echo "setting DB_PSSWRD to be $1";
            DB_PSSWRD=$1;
            shift;
          fi
          ;;
          -n|--db_user)
          shift;
          if [ -n "$1" ]; then
            echo "setting DB_USER to be $1";
            DB_USER=$1;
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
  echo "$0 -p [protein dir] -u [fasta header uniq id col] -d [ortho db name] -h [db_host] -x [db_psswrd] -n [db_user]"
  echo "
protein_dir  :     a directory of FASTA files. One FASTA for each species.
                   Each FSATA file must have a name in the form 'xxxx.fasta'
                   where xxxx is a three or four letter unique taxon code.
                   For example: hsa.fasta or eco.fasta
                   [no default:required]
mySQL_db_name :    mySQL database name
                   For example: orthomcl_andrii
                   [no default:required]
mySQL_username:    mySQL user name
                   For example: Orhtomcl 
                   [no default:required]
mySQL_Psswrd :     mySQL password
                   For example: P@SSW0RD!
                   [no default:required]
mySQL_host   :     mySQL host
                   For example: puffball.fungalgenomes.org 
                   [no default]
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

if [ ! $SPLIT  ] ; then 
  SPLIT=500
fi

echo "PROTEIN_FASTA_DIR is $PROTEIN_FASTA_DIR"
echo "UNIQID id $UNIQID"
echo "SPLIT is $SPLIT"


DIR=`pwd`
cd $DIR
if [ !$PSSWRD ] ; then
  DB_PSSWRD="0rthMC1"
fi
if [ !$DB_USER ] ; then
  DB_USER="orthomcl"
fi

if [ $HOST  ] ; then
  CONNECT_HOST=";host=$HOST"
  ADMIN_HOST="-h $HOST"
else
  CONNECT_HOST=""
  ADMIN_HOST=""
fi
if [ ! -e orthomcl.config ] ; then
  echo "dbVendor=mysql 
dbConnectString=DBI:mysql:database=$DB_NAME$CONNECT_HOST
dbLogin=orthomcl
dbPassword=$PASS
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
mysqladmin drop -f $ADMIN_HOST -u $DB_USER -p$PASS $DB_NAME
mysqladmin create  $ADMIN_HOST -u $DB_USER -p$PASS $DB_NAME


orthomclInstallSchema orthomcl.config



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
fi
cd cleanseq
for file in $DIR/proteomes/*.pep
do
 base=`basename $file .pep`
 ## which field has a uniq ID? => 2
 orthomclAdjustFasta $base $file $UNIQID 
done
cd $DIR

orthomclFilterFasta cleanseq 10 10


## format goodProteins BLAST DB
if [ ! -e $DIR/goodProteins.fasta.pin ] ; then
  makeblastdb -in $DIR/goodProteins.fasta -title "OrthoMCLPeps" -dbtype prot
fi

## split up goodProteins into smaller FASTA
MAX=`expr $SPLIT - 1`
if [ ! -d $DIR/good_proteins_split ] ; then
  mkdir $DIR/good_proteins_split
fi

if [ ! -e $DIR/goodProteins_part_${MAX}.fasta ] ; then
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
 
if [ ! -e goodProteins.BLASTP.bpo ] ; then 
  ## make array BLAST job
  BLAST_Q=1
  echo "#PBS -t 0-${MAX}
  #PBS -l nodes=1:ppn=4
  module load ncbi-blast
  cd $DIR
  blastp -query good_proteins_split/goodProteins_part_\${PBS_ARRAYID}.fasta -db goodProteins.fasta -num_threads 4 -outfmt 6 -out $DIR/blastp_out/goodProteins_part_\${PBS_ARRAYID}_vs_goodProteins.BLASTP.tab -evalue 1e-3" > $DIR/run_split_blast_array.sh


  #run array BLAST job
  BLASTJOB=`qsub run_split_blast_array.sh`
  echo "JOBS SUBMITTED AND HOPEFULLY RUNNING ON YOUR QUEUE:"
  echo "$BLASTJOB = qsub $DIR/run_split_blast_array.sh"

  echo "
  module load orthoMCL
  module load ncbi-blast

  cd $DIR 

  # This will reformat the BLAST data into something to load into the database
  if [ ! -f goodProteins.BLASTP.bpo ]; then
   cat blastp_out/*tab > goodProteins.BLASTP.tab
   orthomclBlastParser goodProteins.BLASTP.tab cleanseq > goodProteins.BLASTP.bpo
  fi

fi


# Load the data into the DB
orthomclLoadBlast orthomcl.config goodProteins.BLASTP.bpo

# now run the ortholog/paralog initial finding
rm -rf pairs pairs.log 
orthomclPairs orthomcl.config  pairs.log cleanup=no

# Dump out the ortholog groups and the mclInput file
orthomclDumpPairsFiles orthomcl.config

# Run mcl for clustering
mcl mclInput  --abc -I 1.5 -o mclOutput.I15.out

# convert the MCL clusters into OrthoMCL groups
orthomclMclToGroups OG 1 < mclOutput.I15.out > mclGroups.I15.table
" > $DIR/finish.orthomcl.sh

if [ $BLAST_Q -eq 0 ]; then
  qsub -W depend=afterokarray:$BLASTJOB $DIR/finish.orthomcl.sh
  echo "qsub -W depend=afterokarray:$BLASTJOB $DIR/finish.orthomcl.sh"
else
  sh $DIR/finish.orthomcl.sh
fi
