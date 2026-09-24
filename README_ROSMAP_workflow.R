## Reprossing RNASeq data for AMP-AD1.0
# J Lundin
# July 14 2026

## WORKFLOW


DE_STEP01_synapse_metadata.R
  Loads and cleans metadata
  output: final md_all file (syn76260726) (n=2807)

DE_STEP02_synapse_counts_tpms_from_synapse.R
  Loads and filters counts and tpms by batch
  output: 4 individual files of filtered batches (b1-b4) (output to local server, temp files)
  
  Moved to DE_STEP02 file: DE_STRP03a_QC_preprocessing_prep_sep_by_tissue.R
  Loads and combines filtered files from step 2 (counts filtered for both counts and tpms thresholds were output)
  output: combined_full_all4.csv (syn74910836)
  
  Check for additional count and tpm filter by tissue and sex
  output: list of gene_ids to remove based on tissue and sex specific filtering (syn75279942)

Technical variables
  technical_stats_multiqc_star.R (syn76227881) ROSMAP_multiqc_star_technical_stats.csv
  technical_stats_fastqc.R (syn76283403) [run on AWS] (ROSMAP_fq_stats.rds contains: basic_stats.txt, phred_per_base.txt, base_content.txt) 
                                                 
ROSMAP_QC_CQN_USEME3_b1_b4.Rmd
  Run checks on RIN thresholds, inferred sex, fast-qc, multiqc, calcuate technical covariate for RNA metrics (plus technical outliers check) -- filter final md and count files
  output: ROSMAP_md_counts_DLPFC_CN_PCC.rds ("metadata", "counts") (syn76226100)
  output: ROSMAP_md_counts_DLPFC_CN_PCC_TC_FC.rds ("metadata", "counts") (syn76425058)
  
  Ran CQN normalization on QC file
  output: ROSMAP_md_counts_cqn_DLPFC_CN_PCC.rds ("metadata", "counts", "dge_cqn") (syn76226123)  
  output: ROSMAP_md_counts_cqn_DLPFC_CN_PCC_FC_TC.rds ("metadata", "counts", "dge_cqn") (syn76425182)  

  Filtered additional samples flagged as outliers for RNA alignment and QC metrics QC
  output: ROSMAP_md_counts_cqn_DLPFC_CN_PCC_FINAL.rds ("metadata", "counts", "dge_cqn") (synxx) ## USED FOR Differential Expression models

                                                 
#here for reference but not longer used 
DE_QC_b1tob4_24July2026.Rmd  
   Outlier detection using PCA (crude)
   Calculate SVs - no longer using bc of PCA of RNA metrics
   Variance partitioning visualization - figure saved
   Outlier detection using PCA (residualized with incorporation of RNA metrics and SVs) - use AWS version instead (on residualized data)
   #output: ROSMAP_md_sv_final_for_DE_models.csv [final metadata file, includes final metadata with RNA metrics and SVs] ] ("md_sv") (synxxx)
