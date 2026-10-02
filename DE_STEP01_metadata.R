
## Differential expression analyses on RNAseq data from ROSMAP
# STEP 01: metadata from synapse
# J Lundin
# 24 July 2026

## load libraries
pacman::p_load(tidyverse, limma, edgeR, biomaRt, DESeq2, vsn, sva, pamr)
pacman::p_load(synapser,dplyr,purrr,readr,lubridate,stringr,tibble,ggplot2)
synLogin()


## pulling metadata from synapse ----
#assay metadata
dlpfcCovObj2 <- synapser::synGet('syn21088596')
metadata_assay <- read.csv(dlpfcCovObj2$path,stringsAsFactors = F)

#biospecimen metadata
dlpfcCovObj2 <- synapser::synGet('syn21323366')
metadata_biosp <- read.csv(dlpfcCovObj2$path,stringsAsFactors = F)
table(metadata_biosp$assay)
metadata_biosp <- metadata_biosp[metadata_biosp$assay == "rnaSeq"& metadata_biosp$tissue != 'blood',]

#clinical metadata
dlpfcCovObj2old <- synapser::synGet('syn3191087') #ROSMAP_clinical.csv
metadata_clinical_old <- read.csv(dlpfcCovObj2old$path,stringsAsFactors = F)
dlpfcCovObj2 <- synapser::synGet('syn73713768') #clinical_harmonized from Jaclyns metadata harmonization project 02/2026 #ROSMAP_clinical_harmonized.csv
metadata_clinical <- read.csv(dlpfcCovObj2$path,stringsAsFactors = F)
table(metadata_clinical$projid)

metadata_temp <- merge(metadata_assay, metadata_biosp, by="specimenID")
metadata_temp2 <- merge(metadata_temp, metadata_clinical, by = "individualID")
metadata_temp2b <- metadata_temp2 %>% filter(organ == "brain")
table(metadata_temp2b$tissue)

comb <- metadata_temp2b


# Batch IDs
# Collapse weird "0, 6, 7" into a single number batch
comb[ comb$libraryBatch %in% "0, 6, 7", ]$sequencingBatch <- 9
comb[ comb$libraryBatch %in% "0, 6, 7", ]$libraryBatch <- 9
table(comb$sequencingBatch)
# Remove letters from data cuts one and two
comb$sequencingBatch <- gsub('NYGC', '', comb$sequencingBatch)
comb$sequencingBatch <- gsub('RISK_', '',comb$sequencingBatch)

colnames(comb)[colnames(comb)=='notes'] <- 'data_contribution'
comb$data_contribution <- gsub('data contribution batch ', '', comb$data_contribution)

comb$final_batch <- paste0(comb$data_contribution, '_', comb$sequencingBatch)
table(comb$final_batch)
table(comb$dataContributionGroup)
table(comb$sex)
table(comb$apoeGenotype, comb$apoe4Status)

# diagnosis
comb$diagnosis <- 'OTHER'
comb[ (comb$Braak == "Stage IV" | comb$Braak == "Stage V" |comb$Braak == "Stage VI" ) &
        (comb$amyCerad == "Frequent/Definite/C3" | comb$amyCerad == "Moderate/Probable/C2" ) &
        comb$cogdx == 4  , ]$diagnosis <- 'AD'
comb[ (comb$Braak == "Stage III" |comb$Braak == "Stage II"|comb$Braak == "Stage I"|comb$Braak == "None") &
  (comb$amyCerad == "Sparse/Possible/C1" | comb$amyCerad == "None/No AD/C0" ) &
        comb$cogdx == 1  , ]$diagnosis <- 'CT'
table(comb$diagnosis, useNA="always")
table(comb$diagnosis, comb$sex, useNA="always")


## New AD definition for the AMP_AD1.0 re-analysis using the DivCo definitions. 
# From "Diagnostic harmonization" subsection (2.4) of Reddy et al 
# AD diagnosis: Braak stage ≥ IV and 
#   CERAD measure equal to moderate/probable AD or frequent/definite AD. 
# Control diagnosis was assigned to individuals with Braak stage ≤ III and 
#   CERAD measure equal to none/no AD or sparse/possible AD. 
# Any donors who did not fall under these criteria were assigned as “other.”

# diag2
comb$diag2 <- 'OTHER2'
comb[ (comb$Braak == "Stage IV" | comb$Braak == "Stage V" |comb$Braak == "Stage VI" ) &
        (comb$amyCerad == "Frequent/Definite/C3" | comb$amyCerad == "Moderate/Probable/C2" ), ]$diag2 <- 'AD2'
comb[ (comb$Braak == "Stage III" |comb$Braak == "Stage II"|comb$Braak == "Stage I"|comb$Braak == "None") &
        (comb$amyCerad == "Sparse/Possible/C1" | comb$amyCerad == "None/No AD/C0" ), ]$diag2 <- 'CT2'
table(comb$diag2, useNA="always")
table(comb$diag2, comb$sex, useNA="always")

vars_b123 <- c("individualID", "specimenID", 'projid', 'cohort', "tissue",
             'diagnosis', 'diag2', 'apoeGenotype', 'apoe4Status', 
             'PMI', 'Braak', 'amyCerad', 'cogdx', 'dcfdx_lv', 'sex', 
             'educ', 'race', 'isHispanic', 'age_at_visit_max', 'age_first_ad_dx',
             'ageDeath', 'cts_mmse30_first_ad_dx', 'cts_mmse30_lv', 'RIN', 
             'libraryBatch', 'sequencingBatch', 'dataContributionGroup', 
             'libraryPrep', 'libraryPreparationMethod','final_batch')

comb <- comb %>% dplyr::select(all_of(vars_b123))

# Sequencing Statistics - UPDATED 7-10-2026
metrics <- read.csv(synapser::synGet('syn76227881')$path, stringsAsFactors = F)

metrics2 <- metrics %>% dplyr::select(specimenID, "picard_UNPAIRED_READS_EXAMINED"  ,          
                               "picard_READ_PAIRS_EXAMINED" ,                  
                               "picard_UNMAPPED_READS"  ,                         
                               "picard_PERCENT_DUPLICATION"   ,                   
                               "picard_READS_UNMAPPED"     ,                      
                               "rsem_alignable_percent" ,                         
                               "rsem_uniquely_aligned_percent"  ,   
                               "samtools_reads_mapped", 
                               "samtools_reads_duplicated",
                               "samtools_error_rate",                             
                               "samtools_average_length"  ,                
                               "samtools_average_quality"  ,                       
                               "samtools_insert_size_average" ,       
                               "samtools_percentage_of_properly_paired_reads_...",
                               "samtools_reads_mapped_percent"  ,                  
                               "samtools_reads_mapped_and_paired_percent"  ,      
                               "samtools_reads_properly_paired_percent"   ,       
                               "samtools_reads_duplicated_percent"     ,                 
                               "samtools_percent_mapped_X",                        
                               "samtools_percent_mapped_Y"      ,                 
                               "cutadapt_percent_trimmed_R2" ,                     
                               "cutadapt_percent_trimmed_R1",                     
                               "cutadapt_mean_percent_trimmed" ,          
                               "Average.input.read.length"  ,                     
                               "Uniquely.mapped.reads.."   ,                      
                               "Number.of.splices..GT.AG" ,                       
                               "Number.of.splices..GC.AG" ,                        
                               "Number.of.splices..AT.AC" ,                       
                               "Mismatch.rate.per.base..." ,                   
                               "Deletion.rate.per.base"   ,                                      
                               "Insertion.rate.per.base"  ,                        
                               "X..of.reads.mapped.to.multiple.loci"  ,           
                               "X..of.reads.mapped.to.too.many.loci"    ,         
                               "X..of.reads.unmapped..too.many.mismatches"   ,    
                               "X..of.reads.unmapped..too.short" ,                
                               "X..of.reads.unmapped..other" ,                    
                               "X..of.chimeric.reads" ) 
comb2 <- merge(comb, metrics2, by.x="specimenID", all.x=T)
comb2 <- comb2[comb2$final_batch != "NA_NA",]
md <- comb2


## batch 4 ----
# Biospecimen 
dlpfcCovObj2 <- synapser::synGet('syn24185681')
metadata_b4_biosp <- read.csv(dlpfcCovObj2$path,stringsAsFactors = F)
metadata_b4_biosp <- metadata_b4_biosp[metadata_b4_biosp$assay == "rnaSeq",]

# assay
dlpfcCovObj3 <- synapser::synGet('syn24185682')
metadata_b4_assay <- read.csv(dlpfcCovObj3$path,stringsAsFactors = F)

table(metadata_b4_biosp$organ)
table(metadata_b4_biosp$tissue)
table(metadata_b4_biosp$assay, metadata_b4_biosp$tissue)

md_b4 <- merge(metadata_b4_biosp, metadata_b4_assay, by="specimenID" )
md_b4 <- merge(metadata_clinical, md_b4, by="individualID" )


# diagnosis
md_b4$diagnosis <- 'OTHER'
md_b4[ (md_b4$Braak == "Stage IV" | md_b4$Braak == "Stage V" |md_b4$Braak == "Stage VI" ) &
        (md_b4$amyCerad == "Frequent/Definite/C3" | md_b4$amyCerad == "Moderate/Probable/C2" ) &
        md_b4$cogdx == 4  , ]$diagnosis <- 'AD'
md_b4[ (md_b4$Braak == "Stage III" |md_b4$Braak == "Stage II"|md_b4$Braak == "Stage I"|md_b4$Braak == "None") &
        (md_b4$amyCerad == "Sparse/Possible/C1" | md_b4$amyCerad == "None/No AD/C0" ) &
        md_b4$cogdx == 1  , ]$diagnosis <- 'CT'
table(md_b4$diagnosis, useNA="always")
table(md_b4$diagnosis, md_b4$sex, useNA="always")

## New AD definition for the AMP_AD1.0 re-analysis using the DivCo definitions. 
# From "Diagnostic harmonization" subsection (2.4) of Reddy et al 
# AD diagnosis: Braak stage ≥ IV and 
#   CERAD measure equal to moderate/probable AD or frequent/definite AD. 
# Control diagnosis was assigned to individuals with Braak stage ≤ III and 
#   CERAD measure equal to none/no AD or sparse/possible AD. 
# Any donors who did not fall under these criteria were assigned as “other.”

# diag2
md_b4$diag2 <- 'OTHER2'
md_b4[ (md_b4$Braak == "Stage IV" | md_b4$Braak == "Stage V" |md_b4$Braak == "Stage VI" ) &
        (md_b4$amyCerad == "Frequent/Definite/C3" | md_b4$amyCerad == "Moderate/Probable/C2" ), ]$diag2 <- 'AD2'
md_b4[ (md_b4$Braak == "Stage III" |md_b4$Braak == "Stage II"|md_b4$Braak == "Stage I"|md_b4$Braak == "None") &
        (md_b4$amyCerad == "Sparse/Possible/C1" | md_b4$amyCerad == "None/No AD/C0" ), ]$diag2 <- 'CT2'
table(md_b4$diag2, useNA="always")
table(md_b4$diag2, md_b4$sex, useNA="always")

md_b4$final_batch <- paste0('4_', md_b4$sequencingBatch)
#table(md_b4$final_batch)

vars_b4 <- c("individualID", "specimenID", 'projid', 'cohort', "tissue",
             'diagnosis', 'diag2', 'apoeGenotype', 'apoe4Status', 
             'PMI', 'Braak', 'amyCerad', 'cogdx', 'dcfdx_lv', 'sex', 
             'educ', 'race', 'isHispanic', 'age_at_visit_max', 'age_first_ad_dx',
             'ageDeath', 'cts_mmse30_first_ad_dx', 'cts_mmse30_lv', 'RIN', 
             'libraryBatch', 'sequencingBatch', 'dataContributionGroup', 
             'libraryPrep', 'libraryPreparationMethod','final_batch')

md_b4_temp <- md_b4 %>% dplyr::select(all_of(vars_b4))

md_b4_metrics <- merge(md_b4_temp , metrics2, by.x="specimenID", all.x=T)

md_b4 <- md_b4_metrics


## BRINGING DATA TOGETHER BY BRAIN REGION ----
## metadata files for batch 1 through batch 4 ----

md_all <- rbind(md, md_b4)

md_all <- md_all %>% separate_wider_delim(final_batch, delim = "_", names = c("procbatch", "seqbatch"), cols_remove = FALSE) 
md_all <- md_all %>% mutate(across(c(procbatch, seqbatch), as.numeric)) 
md_all <- md_all %>% arrange(procbatch, seqbatch)
table(md_all$final_batch)
summary(md$RIN)
summary(md_b4$RIN)

md_all$age_cat[md_all$ageDeath == "90+"] <- "90+"
md_all$age_cat[as.numeric(md_all$ageDeath[md_all$ageDeath != "90+"])<85] <- "<85"
md_all$age_cat[as.numeric(md_all$ageDeath[md_all$ageDeath != "90+"])>=85] <- "ge85lt90"
md_all$age_cat[md_all$ageDeath == "90+"] <- "90+"
table(md_all$age_cat, useNA="always")


md_all <- md_all %>% mutate(tissue2 = recode(tissue, 
                                             "dorsolateral prefrontal cortex" = "DLPFC",
                                             "frontal cortex"  = "FC",
                                             "Head of caudate nucleus"  = "CN",
                                             "posterior cingulate cortex"   = "PCC",
                                             "temporal cortex"   = "TC"))  

## send raw metadata to synapse ----
temp_path_md_all <- tempfile(fileext = ".csv")
write_csv(md_all, temp_path_md_all)
syn_md_all <- File(path = temp_path_md_all, name = "md_all.csv", parent = "syn75192265")
synStore(syn_md_all, forceVersion = FALSE)


