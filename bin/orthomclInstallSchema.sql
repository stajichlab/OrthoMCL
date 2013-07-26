USE orthomcl;

--Remove schema if it already exists

DROP TABLE IF EXISTS orthomcl.SimilarSequences;
DROP TABLE IF EXISTS orthomcl.InParalog;
DROP TABLE IF EXISTS orthomcl.Ortholog;
DROP TABLE IF EXISTS orthomcl.CoOrtholog;
DROP TABLE IF EXISTS orthomcl.BestInterTaxonScore;
DROP TABLE IF EXISTS orthomcl.BestQueryTaxonScore;
DROP TABLE IF EXISTS orthomcl.InterTaxonMatch;

--Create schema

CREATE TABLE orthomcl.SimilarSequences (
 QUERY_ID                 VARCHAR(15),
 SUBJECT_ID               VARCHAR(15),
 QUERY_TAXON_ID           VARCHAR(15),
 SUBJECT_TAXON_ID         VARCHAR(15),
 EVALUE_MANT              BIGINT(20),
 EVALUE_EXP               BIGINT(20),
 PERCENT_IDENTITY         FLOAT,
 PERCENT_MATCH            FLOAT  
);



CREATE INDEX ss_qtaxexp_ix ON orthomcl.SimilarSequences(query_id, subject_taxon_id, evalue_exp, evalue_mant, query_taxon_id, subject_id);
CREATE INDEX ss_seqs_ix ON orthomcl.SimilarSequences(query_id, subject_id, evalue_exp, evalue_mant);


-----------------------------------------------------------

CREATE TABLE orthomcl.InParalog (
 SEQUENCE_ID_A           VARCHAR(15),
 SEQUENCE_ID_B           VARCHAR(15),
 TAXON_ID                VARCHAR(15),
 UNNORMALIZED_SCORE      DOUBLE,
 NORMALIZED_SCORE        DOUBLE    
);


------------------------------------------------------------

CREATE TABLE orthomcl.Ortholog (
 SEQUENCE_ID_A           VARCHAR(15),
 SEQUENCE_ID_B           VARCHAR(15),
 TAXON_ID_A              VARCHAR(15),
 TAXON_ID_B              VARCHAR(15),
 UNNORMALIZED_SCORE      DOUBLE,
 NORMALIZED_SCORE        DOUBLE    
);

CREATE INDEX orthomcl.ortholog_seq_a_ix on orthomcl.ortholog(sequence_id_a);
CREATE INDEX orthomcl.ortholog_seq_b_ix on orthomcl.ortholog(sequence_id_b);


-------------------------------------------------------------
 
CREATE TABLE orthomcl.CoOrtholog (
 SEQUENCE_ID_A           VARCHAR(15),
 SEQUENCE_ID_B           VARCHAR(15),
 TAXON_ID_A              VARCHAR(15),
 TAXON_ID_B              VARCHAR(15),
 UNNORMALIZED_SCORE      DOUBLE,
 NORMALIZED_SCORE        DOUBLE    
);




CREATE OR REPLACE VIEW orthomcl.InterTaxonMatch 
	AS SELECT ss.query_id, ss.subject_id, ss.subject_taxon_id, 
	ss.evalue_mant, ss.evalue_exp 
	FROM orthomcl.SimilarSequences ss 
	WHERE ss.subject_taxon_id != ss.query_taxon_id;




-- exit;
