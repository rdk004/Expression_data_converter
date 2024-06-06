# ---------------------------------------------------------------------------- #
# SECTION 1
# Let us first load all the libraries required 

library(readr)
library(data.table)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer")
BiocManager::install("GenomicRanges")
library(rtracklayer)
library(GenomicRanges)
library(dplyr)

# ---------------------------------------------------------------------------- #
# SECTION 2
# Import the .tsv file

tsv_df <- read_tsv('your file path')
# ---------------------------------------------------------------------------- #
# SECTION 3
# Apply the antilog transformation to count data (Assuming Ensembl_ID is the identifier column)

anti_log_transformed_df <- tsv_df %>% 
  mutate(across(where(is.numeric), ~ 2^. - 1))

# Create a dataframe of Ensembl Gene IDs only

# Let us create a vector of Ensembl Gene IDs from the original dataframe

ensembl_ids_cancername_genes<- tsv_df["Ensembl_ID"]
# ---------------------------------------------------------------------------- #
# SECTION 4
# Load the GTF file and extract the gene lengths of correspoding Ensembl IDs

gtf <- import('path of gtf file in system')

# Extract the gene annotations
gene_annotations <- gtf[gtf$type == "gene"]

# Calculate gene lengths
gene_lengths <- width(gene_annotations)

# Create a data frame with gene IDs and their lengths
gene_df <- data.frame(
  gene_id = gene_annotations$gene_id,
  length = gene_lengths
)

colnames(gene_df)[1] <- "Ensembl_ID"

# Match Ensembl IDs to gene lengths
gene_lengths_ensembl_id <- gene_df[match(ensembl_ids_cancername_genes$Ensembl_ID, gene_df$Ensembl_ID), ]

# We want gene lengths in kb, while the GTF files have it in b. So we divide by 1000

gene_lengths_ensembl_id$length <- gene_lengths_ensembl_id$length / 1000 # Now we have the lengths of every gene
# ---------------------------------------------------------------------------- #
# SECTION 5
# TPM transformations

# Formula description: # TPM = x/sum(x) * 1e6 (where x = c/l; c is the count value (antilog transformed) and l is the length of the gene in kb; sum(x) is all x values summed for a sample)

# Step 1: Let us apply the first transformation - normalize each element with corresponding gene length

# Assuming your data frames are named anti_log_transformed_tsv_df and gene_lengths_ensembl_id
# Extract numeric columns (excluding gene IDs)
count_matrix <- anti_log_transformed_tsv_df[, -1]

# Initialize an empty data frame for normalized counts
normalized_counts_df <- data.frame(matrix(NA, nrow = nrow(count_matrix), ncol = ncol(count_matrix)))
colnames(normalized_counts_df) <- colnames(count_matrix)

# Normalize each sample (column) by gene lengths
for (col in colnames(count_matrix)) {
  normalized_counts_df[, col] <- count_matrix[, col] / gene_lengths_ensembl_id$length
}

# Round each value to 2 decimal places
gene_len_norm_twodec <- round(normalized_counts_df, digits = 2)
class(gene_len_norm_twodec)
class(ensembl_id_luad_genes)

# Lets add the ensembl ids as rownames back to the dataframe using cbind()

gene_len_norm_df <- cbind(ensembl_id_luad_genes, gene_len_norm_twodec)
# ---------------------------------------------------------------------------- #

# Step 2: Let us apply the second transformation - divide each every element by column_sum/10^6

# Assuming your dataframe is named 'my_df'
numeric_cols <- sapply(gene_len_norm_df, is.numeric)
col_sum_vector <- 1000000/(colSums(gene_len_norm_df[, numeric_cols], na.rm = TRUE))
class(col_sum_vector)

# Remove the first column (Ensembl gene IDs)
numeric_df <- gene_len_norm_df[, -1]

# Construct a diagonal matrix from the vector
diag_matrix <- diag(col_sum_vector)

# Multiply the data frame by the diagonal matrix
result_matrix <- as.matrix(numeric_df) %*% diag_matrix

# Convert the result back to a data frame
result_df <- as.data.frame(result_matrix)

tpm_df <- round(result_df, digits = 2)

# Final TPM data

final_tpm_df <- cbind(ensembl_id_luad_genes, tpm_df)

# Access column names from df1
column_names <- colnames(gene_len_norm_df)

colnames(final_tpm_df) <- column_names

View(final_tpm_df) # This is your converted TPM dataframe
# ---------------------------------------------------------------------------- #
# SECTION 6
# Let us export this TPM dataframe:
write.csv(final_tpm_df, file = '/home/nsclab/Rishabh_Kulkarni/Lab_project_2024/Datasets/TPM files/LUAD/luad_tcga_tpm.csv', row.names = FALSE)
# ---------------------------------------------------------------------------- #
