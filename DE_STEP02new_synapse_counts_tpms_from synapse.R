## Differential expression analyses on RNAseq data from ROSMAP
# STEP 02: count and tpm data from synapse
# J Lundin
# 24 July 2026

## Load libraries
pacman::p_load(tidyverse, limma, edgeR, biomaRt, DESeq2, vsn, sva, pamr)
pacman::p_load(synapser,dplyr,purrr,readr,lubridate,stringr,tibble,ggplot2)
pacman::p_load(httr, purrr)
synLogin()

# get functions from git
urls <- c("https://raw.githubusercontent.com/jessica-lundin-sage/pers_Lundin/main/RNASeq_DE/AMP_AD1/functions/functions_filter_low_count_genes.R")
walk(urls, function(u) {
  resp <- GET(u, add_headers(Authorization = paste("token", Sys.getenv("GITHUB_PAT"))))
  stop_for_status(resp)
  tmp <- tempfile(fileext = ".R")
  writeLines(content(resp, "text", encoding = "UTF-8"), tmp)
  source(tmp)
})

# metadata
md_all <- read.csv(synapser::synGet('syn76260726')$path, header=TRUE)


## filtering on low counts ----
filtered_genes_b1 <- filter_gene_expression(
  tpm_file   = read.csv(synapser::synGet("syn69369139")$path, sep="\t", header=TRUE, check.names = FALSE),
  reads_file  = read.csv(synapser::synGet("syn69369134")$path, sep="\t", header=TRUE, check.names = FALSE),
  metadata  = read.csv(synapser::synGet('syn76260726')$path, header=TRUE),
  synid_outfile = c("syn75192265"),
  study_var = c("ROSMAP_b1")
)

filtered_genes_b2 <- filter_gene_expression(
  tpm_file   = read.csv(synapser::synGet("syn69681788")$path, sep="\t", header=TRUE, check.names = FALSE),
  reads_file  = read.csv(synapser::synGet("syn69681787")$path, sep="\t", header=TRUE, check.names = FALSE),
  metadata  = read.csv(synapser::synGet('syn76260726')$path, header=TRUE),
  synid_outfile = c("syn75192265"),
  study_var = c("ROSMAP_b2")
)

filtered_genes_b3 <- filter_gene_expression(
  tpm_file   = read.csv(synapser::synGet("syn69369207")$path, sep="\t", header=TRUE, check.names = FALSE),
  reads_file  = read.csv(synapser::synGet("syn69369206")$path, sep="\t", header=TRUE, check.names = FALSE),
  metadata  = read.csv(synapser::synGet('syn76260726')$path, header=TRUE),
  synid_outfile = c("syn75192265"),
  study_var = c("ROSMAP_b3")
)

filtered_genes_b4 <- filter_gene_expression(
  tpm_file    = read.csv(synapser::synGet("syn69369251")$path, sep="\t", header=TRUE, check.names = FALSE),
  reads_file  = read.csv(synapser::synGet("syn69369249")$path, sep="\t", header=TRUE, check.names = FALSE),
  metadata  = read.csv(synapser::synGet('syn76260726')$path, header=TRUE),
  synid_outfile = c("syn75192265"),
  study_var = c("ROSMAP_b4")
)


## previously DE_STEP03a_QC_preprocessing_prep_sep_by_tissue.R

## filtered gene count data  ----
# create gene reads dataset for COMBINED based on filtered genes across 4 batches
rosmap_counts_b1_filtered <- read.table(synapser::synGet('syn77437138')$path, sep="\t", check.names = FALSE, header=T)
rosmap_counts_b2_filtered <- read.table(synapser::synGet('syn77454077')$path, sep="\t", check.names = FALSE, header=T)
rosmap_counts_b3_filtered <- read.table(synapser::synGet('syn77454156')$path, sep="\t", check.names = FALSE, header=T)
rosmap_counts_b4_filtered <- read.table(synapser::synGet('syn77454242')$path, sep="\t", check.names = FALSE, header=T)

# Find IDs present in all three
common_ids <- intersect(intersect(rosmap_counts_b1_filtered$gene_id, rosmap_counts_b2_filtered$gene_id), rosmap_counts_b3_filtered$gene_id)
common_ids_all <- intersect(common_ids, rosmap_counts_b4_filtered$gene_id)

# Subset the data frames
batch1_clean <- rosmap_counts_b1_filtered[rosmap_counts_b1_filtered$gene_id %in% common_ids_all, ]
batch2_clean <- rosmap_counts_b2_filtered[rosmap_counts_b2_filtered$gene_id %in% common_ids_all, ]
batch3_clean <- rosmap_counts_b3_filtered[rosmap_counts_b3_filtered$gene_id %in% common_ids_all, ]
batch4_clean <- rosmap_counts_b4_filtered[rosmap_counts_b4_filtered$gene_id %in% common_ids_all, ]

batch1_clean <- batch1_clean[,-2]
batch2_clean <- batch2_clean[,-2]
batch3_clean <- batch3_clean[,-2]
batch4_clean <- batch4_clean[,-2]


# 1. Put datasets into a list
list_of_batches <- list(batch1_clean, batch2_clean, batch3_clean, batch4_clean)

# 2. Join them all together by "gene_id"
combined_full_all4 <- list_of_batches %>% 
  purrr::reduce(full_join, by = "gene_id")

# Total NAs in the entire dataset
sum(is.na(combined_full_all4 ))

## send count data filtered for low gene/tpm counts to synapse ----
temp_path_combined_full_all4 <- tempfile(fileext = ".csv")
write_csv(combined_full_all4, temp_path_combined_full_all4)
syn_full_all4 <- File(path = temp_path_combined_full_all4, name = "combined_full_all4.csv", parent = "syn74831376")
synStore(syn_full_all4, forceVersion = FALSE)


## Create a filtered gene tpm data from STEP02 (file above filtered for both counts and tpms, but only saved counts) ----

## filtered gene tpm data  ----
# create gene tpms dataset for COMBINED based on filtered genes across 4 batches
rosmap_tpm_b1_filtered <- read.table(synapser::synGet('syn77437145')$path, sep="\t", check.names = FALSE, header=T)
rosmap_tpm_b2_filtered <- read.table(synapser::synGet('syn77454103')$path, sep="\t", check.names = FALSE, header=T)
rosmap_tpm_b3_filtered <- read.table(synapser::synGet('syn77454165')$path, sep="\t", check.names = FALSE, header=T)
rosmap_tpm_b4_filtered <- read.table(synapser::synGet('syn77454243')$path, sep="\t", check.names = FALSE, header=T)


# create gene reads dataset for COMBINED based on filtered genes across 4 batches
# Find IDs present in all three
common_ids_tpm <- intersect(intersect(rosmap_tpm_b1_filtered$gene_id, rosmap_tpm_b2_filtered$gene_id), rosmap_tpm_b3_filtered$gene_id)
common_ids_tpm_all <- intersect(common_ids_tpm, rosmap_tpm_b4_filtered$gene_id)

# Subset the data frames
batch1_tpm_clean <- rosmap_tpm_b1_filtered[rosmap_tpm_b1_filtered$gene_id %in% common_ids_tpm_all, ]
batch2_tpm_clean <- rosmap_tpm_b2_filtered[rosmap_tpm_b2_filtered$gene_id %in% common_ids_tpm_all, ]
batch3_tpm_clean <- rosmap_tpm_b3_filtered[rosmap_tpm_b3_filtered$gene_id %in% common_ids_tpm_all, ]
batch4_tpm_clean <- rosmap_tpm_b4_filtered[rosmap_tpm_b4_filtered$gene_id %in% common_ids_tpm_all, ]

batch1_tpm_clean <- batch1_tpm_clean[,-2]
batch2_tpm_clean <- batch2_tpm_clean[,-2]
batch3_tpm_clean <- batch3_tpm_clean[,-2]
batch4_tpm_clean <- batch4_tpm_clean[,-2]


# 1. Put datasets into a list
list_of_batches_tpm <- list(batch1_tpm_clean, batch2_tpm_clean, batch3_tpm_clean, batch4_tpm_clean)

# 2. Join them all together by "gene_id"
combined_full_tpm_all4 <- list_of_batches_tpm %>% 
  purrr::reduce(full_join, by = "gene_id")

combined_full_tpm_all4 <- combined_full_tpm_all4[combined_full_tpm_all4$gene_id %in% combined_full_all4$gene_id,] # filter from above already applied to counts, now matching gene_ids for tpms


# Total NAs in the entire dataset
sum(is.na(combined_full_all4 ))

## send count data filtered for low gene/tpm counts to synapse ----
temp_path_combined_full_tpm_all4 <- tempfile(fileext = ".csv")
write_csv(combined_full_tpm_all4, temp_path_combined_full_tpm_all4)
syn_full_tpm_all4 <- File(path = temp_path_combined_full_tpm_all4, name = "combined_full_tpm_all4.csv", parent = "syn75192265")
synStore(syn_full_tpm_all4, forceVersion = FALSE)



## process tpms and counts filtering again but stratified by tissue and sex ----
library(tidyverse)

sex_var <- c("female","male")
tissue2_var <- c("DLPFC","PCC","CN","FC","TC")
mat1<-NULL
mat2<-NULL
for (i in 1:length(tissue2_var)){
  for (j in 1:length(sex_var)){
    combined_full_all4 <- read.csv(synapser::synGet('syn74910836')$path, stringsAsFactors = F, check.names = FALSE)
    combined_full_tpm_all4 <- read.csv(synapser::synGet('syn75276608')$path, stringsAsFactors = F, check.names = FALSE)
    md_all <- read.csv(synapser::synGet('syn76260726')$path, stringsAsFactors = F)
    
    
    read_thresh = 6
    prop_samples = 0.2
    tpm_thresh = 0.1
    
    # 1. Load Data
    gene_reads <- combined_full_all4
    gene_tpm   <- combined_full_tpm_all4 
    pheno      <- md_all %>% filter(sex == sex_var[j] & tissue2 == tissue2_var[i])
    
    # 2. Match IDs
    # Get list of sample IDs present in phenotype file
    RNA_link <- pheno$specimenID
    
    # Filter columns: Keep 'Name', 'Description', and any columns matching the phenotype IDs
    keep_cols_tpm   <- intersect(names(gene_tpm), RNA_link)
    keep_cols_tpm   <- c("gene_id", keep_cols_tpm)
    
    keep_cols_reads <- intersect(names(gene_reads), RNA_link)
    keep_cols_reads <- c("gene_id", keep_cols_reads)
    
    gene_tpm   <- gene_tpm[, keep_cols_tpm]
    gene_reads <- gene_reads[, keep_cols_reads, drop = FALSE]
    
    
    # 4. Filter by Reads Threshold
    rownames(gene_reads) <- gene_reads$gene_id
    read_matrix <- as.matrix(gene_reads[, 2:ncol(gene_reads)])
    passed_reads <- rowMeans(read_matrix >= read_thresh) >= prop_samples
    cat(paste("Table for sex and tissue specific test of counts >=6 in >=20% of samples for ", sex_var[j]," ",tissue2_var[i]))
    cat(" ")
    table(passed_reads)
    cat(table(passed_reads))
    remove_for_strat <- names(passed_reads)[!passed_reads]
    
    # 3. Filter by TPM Threshold
    # Logic: rowMeans on a logical matrix is much faster than sapply/transpose
    rownames(gene_tpm) <- gene_tpm$gene_id
    tpm_matrix <- as.matrix(gene_tpm[, 2:ncol(gene_tpm)])
    passed_tpm <- rowMeans(tpm_matrix >= tpm_thresh) >= prop_samples
    cat(paste("Table for sex and tissue specific test of tpm >=0.1 in >=20% of samples for ", sex_var[j]," ",tissue2_var[i]))
    cat(" ")
    table(passed_tpm)
    cat(table(passed_tpm))
    remove_for_strat_tpm <- names(passed_tpm)[!passed_tpm]
    
    final_genes_remove <- unique(names(passed_tpm)[!passed_tpm], 
                                 names(passed_reads)[!passed_reads])
    
    if (length(final_genes_remove) != 0){
      mat2 <- data.frame(gene_id = final_genes_remove, stringsAsFactors = FALSE)
      mat2$sex <- sex_var[j]
      mat2$tissue2 <- tissue2_var[i]
      
      mat1 <- rbind(mat1, mat2)
    }
    
    # write to synapse
    file_path <- paste("Remove gene_ids from strat analysis reads and tpm_", sex_var[j],"_",tissue2_var[i],".txt",sep="")
    write.table(final_genes_remove, file = file_path, row.names = FALSE)
    synapser::synStore(synapser::File(path = file_path, parent = "syn75192265"))
    
  }
}

table(mat1$sex, mat1$tissue2)

mat1 <- as.data.frame(mat1)

mat1$gene_id_clean <- sub("\\..*", "", mat1$gene_id)

# 1. Install and load
BiocManager::install("biomaRt")
library(biomaRt)

# 2. Connect to the Ensembl human database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 3. Define your list of target Gene IDs
my_ensembl_ids <- unique(mat1$gene_id_clean)

# 4. Fetch the chromosome, start, and end positions
gene_positions <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"),
  filters = "ensembl_gene_id",
  values = my_ensembl_ids,
  mart = mart
)

mat1 <- merge(mat1, gene_positions, by.x="gene_id_clean", by.y="ensembl_gene_id", all=T)


# write to synapse
file_path <- paste("Remove gene_ids from stratified analysis reads and tpm.csv",sep="")
write.table(mat1, file = file_path, row.names = FALSE, sep="\t")
synapser::synStore(synapser::File(path = file_path, parent = "syn75192265"))


