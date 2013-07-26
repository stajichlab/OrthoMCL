USE orthomcl;

LOAD DATA
LOCAL INFILE "~/orthoPackage/data/orthomcl_all_3.lst"
REPLACE INTO TABLE orthomcl.SimilarSequences
FIELDS TERMINATED BY '\t'
