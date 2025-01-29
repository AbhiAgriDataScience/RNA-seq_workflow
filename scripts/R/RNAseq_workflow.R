# RNA seq analysis
# This script is a second phase of the RNA seq analysis after featureCounts analysis

#### ======================= Step 0: Set up directories================================================
# Get the current working directory
getwd()
# "C:/Users/abhis/OneDrive/Documents/RNAseq_analysis"

# Set the working directory to the same location of where the other initial half analysis was carried out
# Could not do that since the files are in linux system so copied the folder to a current directory

# Set a new working directory
setwd("~/RNAseq_analysis/RNAseq/RNA_seq_workflow")

# Create a folder for all R analysis
new_folder_path = "~/RNAseq_analysis/RNAseq/RNA_seq_workflow/R analysis"
dir.create(new_folder_path)

####===============  Step 1: Load required libraries===================================================
BiocManager::install('edgeR', ask = FALSE, update = TRUE)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(dplyr)

BiocManager::install('pheatmap', ask = FALSE, update = TRUE)
library(pheatmap)

BiocManager::install('clusterProfiler', ask = FALSE, update = TRUE)
library(clusterProfiler)

####=============== Step 2: Load data ==========================================
counts = read.table('~/RNAseq_analysis/RNAseq/RNA_seq_workflow/output/4_final_counts/gene_counts.txt',
                    header = TRUE, # use first row as column names
                    comment.char = '#', # skip lines starting with '#'
                    sep = '\t', # tab-separated file
                    check.names = FALSE # prevent R from altering column names
                    )
head(counts)

# Extract relevant columns
gene_lengths = counts$Length # Gene length for normalization
head(gene_lengths)
rownames(counts) = counts$Geneid # Use Gene IDs as row names
count_data = counts[, 7:ncol(counts)] # Extract row counts

# Rename column names of count_data
colnames(count_data) = c("barley_control 1",
                         "barley_control 2",
                         "barley_control 3",
                         "barley_control 4",
                         "barley_Pst 1",
                         "barley_Pst 2",
                         "barley_Pst 3",
                         "barley_Pst 4"
)
head(count_data)

####===============Step 3: Data Quality Assessment =======================================
# 3.1 Basic statistics
library_sizes = colSums(count_data)
summary(library_sizes)

# Boxplot of raw counts
boxplot(log2(count_data + 1), las = 2, main = 'Raw Counts Distribution', col = 'skyblue')

# PCA on raw counts
pca = prcomp(t(log2(count_data + 1)))

# Create a PCA plot
plot(pca$x[,1], pca$x[,2], xlab = 'PC1', ylab = 'PC2', main = 'PCA of Raw Counts', pch = 19, col = 'dodgerblue')

# Add sample names as labels
text(pca$x[,1], pca$x[,2],
     labels = colnames(count_data), # use column names as labels
     pos = 4, # poistion the labels to the right of points
     cex= 0.7, # adjust lable size
     col = 'red')

####===============Step 4: Normalization ==================================================
# 1. TMM normalization (for Differential Expression)

# Create DGEList object (edgeR package)
dge = DGEList(counts = count_data)

# Normalize using TMM
dge = calcNormFactors(dge)

# Access normalized counts
tmm_counts = cpm(dge, normalized.lib.sizes = TRUE) # counts per million

head(tmm_counts)

# 2. Fragments Per Kilobase of transcript per Million mapped reads (FPKM): Normalizes library size and gene length
# Convert counts to FPKM
counts_to_fpkm = function(counts, gene_length, library_sizes) {
  fpkm = t(t(counts)/ library_sizes) * 1e6 # Normalize by library size
  fpkm = fpkm/(gene_length/1e3) # Normnalize by gene length (kb)
  return (fpkm)
}

#  Calculate FPKM
fpkm_counts = counts_to_fpkm(count_data, gene_lengths, library_sizes)

# 3. Transcripts per million (TPM) Calculation (we won't use this): A more modern metric similar to FPKM but easier to interpret and compare across samples.
counts_to_tpm = function(counts, gene_lengths) {
  rpk = counts/(gene_lengths/1e3) # reads per kilobase
  scaling_factors = colSums(rpk)
  tpm = t(t(rpk) / scaling_factors) * 1e6
  return(tpm)
}

# Calculate TPM
tpm_counts = counts_to_tpm(count_data, gene_lengths)

####===============Step 5: Differential Expression Analysis ===============================
# Experimental design
group = factor(c('Control', 'Control', 'Control', 'Control', 'Pst', 'Pst','Pst','Pst'))

# Add experimetnal design to DGE object
dge$samples$group = group


# Estimate dispersions
dge = estimateDisp(dge, design = model.matrix(~ group))

# Perform differential expression analysis
fit = glmFit(dge, design = model.matrix(~ group))
lrt = glmLRT(fit) # genewise likelihood ratio test

# Extract results
results = topTags(lrt, n = Inf)$table
write.table(results, "R analysis/differential_expression_results.txt", sep = "\t", quote = FALSE, row.names = TRUE)


####===============Step 6: Data Visualization ===================================================
# MA plot
plotMD(lrt, main = "MA Plot", xlim = c(-5, 5))

# Volcano plot
results$Significant <- results$FDR < 0.05
ggplot(results, aes(x = logFC, y = -log10(FDR), color = Significant)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log Fold Change", y = "-log10(FDR)")

# Heatmap of significant genes
library(pheatmap)
sig_genes = rownames(results[results$FDR < 0.05, ])

pheatmap(log2(tmm_counts[sig_genes, ] + 1), cluster_rows = TRUE, cluster_cols = TRUE, main = "Significant Genes Heatmap")

list(sig_genes)

#### ===============Step 7: Functional Enrichment Analysis =====================================================================

# This isnt needed either
# We don't require these steps since rownames have already been assigned with gene_id. Therefore, we can directly proceed with Functional Enrichment Analysis
# Install AnnotationHub for Plant Data
#library(AnnotationHub)

# Search for barley annotations
#ah = AnnotationHub()
#query(ah, c("Hordeum", "Morex"))

# nothing in Annotation hub

# Therefore, we will proceed with our own custom Annotation (GFF3)

#library(clusterProfiler)

#if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")
#library(org.Hv.eg.db) # didnt work

###didnt work. Will start again from codes below

# Step 1: Extract relevant information from GFF3 file

# The GFF3 file contains detailed gene annotations, including GO terms. You can parse this file in R to extract the required information.

# Install required packages if not already installed
#if (!requireNamespace("rtracklayer", quietly = TRUE)) 
 # install.packages("BiocManager"); BiocManager::install("rtracklayer")

# Load the required library
#library(rtracklayer)

# # Load GFF3 file
# gff_file <- "annotation/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.60.chr.gff3"
# gff_data <- import.gff3(gff_file)
# 
# 
# head(gff_data)
# 
# #Extract Gene IDs
# # Filter rows for genes
# genes <- gff_data[gff_data$type == "gene", ]
# head(genes)
# 
# # Create a data frame of Gene IDs and Names
# gene_info <- data.frame(
#   GeneID = genes$gene_id,
#   Description = genes$description
# )
# 
# # Preview extracted data
# head(gene_info)
# 
# # Handling missing values (GO.term.accession)
# gene_info = gene_info[!is.na(gene_info$GeneID), ]
# 
# # Rename GeneID to GID in gene_info
# colnames(gene_info)[colnames(gene_info) == 'GeneID'] = "GID"
# 
# # Step 2: Map Gene IDs to GO terms
# # Load GO mapping data
# go_mapping <- read.table("annotation/mart_export.txt", header = TRUE, sep = ",", fill = TRUE, comment.char = "")  # Replace with actual file
# head(go_mapping)
# 
# # Rename Gene.stable.ID to GeneID in go_mapping
# colnames(go_mapping)[colnames(go_mapping) == 'Gene.stable.ID'] = "GID"
# 
# # Merge with gene_info
# gene_go <- merge(gene_info, go_mapping, by = "GeneID", all.x = TRUE)
# 
# # Preview the merged data
# head(gene_go)
# 
# # Handling missing values (GO.term.accession)
# gene_go = gene_go_terms[!is.na(gene_go_terms$GO.term.accession), ]
# 
# 
# #####  Create a custom org.Hv.eg.db package
# # Format GO data for OrgDb creation
# gene_go_formatted <- data.frame(
#   GID = gene_go$GeneID,
#   GO = gene_go$GO.term.accession,
#   EVIDENCE = "IEA"  # Example evidence code: Inferred from Electronic Annotation
# )
# 
# # Remove duplicates from gene_go_formatted
# gene_go_formatted <- gene_go_formatted[!duplicated(gene_go_formatted), ]
# 
# # Verify there are no duplicates
# any(duplicated(gene_go_formatted))  # Should return FALSE
# 
# # Replace NULL or NA values in the Description column
# gene_info$Description[is.na(gene_info$Description)] <- "No description available"
# 
# # Step 5: Create Custom OrgDb Package
# # BiocManager::install('AnnotationForge', update = TRUE)
# library(AnnotationForge)
# 
# makeOrgPackage(
#   gene_info = gene_info,
#   go = gene_go_formatted,
#   version = "0.1",
#   maintainer = "Abhishek Shrestha <abhishek.shrestha39@gmail.com>",
#   author = "Abhishek Shrestha <abhishek.shrestha39@gmail.com>",
#   outputDir = ".",
#   tax_id = "4513",  # Taxonomy ID for Hordeum vulgare
#   genus = "Hordeum",
#   species = "vulgare",
#   goTable = "go"
# )


######  Perform GO Enrichment Analysis


#### GO Enrichment analysis part II
# Step 1: Extract relevant information from GFF3 file

# The GFF3 file contains detailed gene annotations, including GO terms. You can parse this file in R to extract the required information.

# Install required packages if not already installed
#if (!requireNamespace("rtracklayer", quietly = TRUE)) 
 # install.packages("BiocManager"); BiocManager::install("rtracklayer")

# Load the required library
library(rtracklayer)

# Load GFF3 file
#gff_file <- "annotation/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.60.chr.gff3"
#gff_data <- import.gff3(gff_file)


#head(gff_data)

#Extract Gene IDs
# Filter rows for genes
#genes <- gff_data[gff_data$type == "gene", ]
#head(genes)

# Create a data frame of Gene IDs and Names
#gene_info <- data.frame(
 # GeneID = genes$gene_id,
  #Description = genes$description)

# Preview extracted data
#head(gene_info)

# Handling missing values (GO.term.accession)
#gene_info = gene_info[!is.na(gene_info$GeneID), ]

# Rename GeneID to GID in gene_info
#colnames(gene_info)[colnames(gene_info) == 'GeneID'] = "GeneID"

# # Step 2: Map Gene IDs to GO terms
# # Load GO mapping data
# gene_go_terms <- read.table("annotation/mart_export.txt", header = TRUE, sep = ",", fill = TRUE, comment.char = "")  # Replace with actual file
# head(gene_go_terms)
# 
# # Rename Gene.stable.ID to GeneID in go_mapping
# colnames(gene_go_terms)[colnames(gene_go_terms) == 'Gene.stable.ID'] = "GeneID"
# 
# # Create a named list of GO terms for each gene
# go_list = split(gene_go_terms$GO.term.accession, gene_go_terms$GeneID)
# 
# # Input: Gene of interest
# gene_list = sig_genes
# 
# # Perform enrichment analysis using clusterProfiler
# enrichment_result <- enricher(
#   gene = gene_list,
#   TERM2GENE = gene_go_terms[, c("GO.term.accession", "GeneID")]
# )
# 
# # View enrichment results
# head(enrichment_result)
# 
# # Visualize top GO terms
# library(ggplot2)
# dotplot(enrichment_result, showCategory = 10)  # Adjust showCategory as needed


######------- Case 2: Remove barley_Pst1 and barley_control 2 samples

# # Load the counts data
# counts = read.table('~/RNAseq_analysis/RNAseq/RNA_seq_workflow/output/4_final_counts/gene_counts.txt',
#                     header = TRUE,
#                     comment.char = '#',
#                     sep = '\t',
#                     check.names = FALSE)
# 
# # Extract relevant columns
# gene_lengths = counts$Length
# rownames(counts) = counts$Geneid
# count_data = counts[, 7:ncol(counts)]
# 
# # Rename column names of count_data
# colnames(count_data) = c("barley_control 1",
#                          "barley_control 2",
#                          "barley_control 3",
#                          "barley_control 4",
#                          "barley_Pst 1",
#                          "barley_Pst 2",
#                          "barley_Pst 3",
#                          "barley_Pst 4")
# 
# # Drop `barley_control 2` and `barley_Pst 1`
# count_data = count_data[, -c(2, 5)]
# colnames(count_data)
# 
# ####===============Step 3: Data Quality Assessment =======================================
# # 3.1 Basic statistics
# library_sizes = colSums(count_data)
# summary(library_sizes)
# 
# # Boxplot of raw counts
# boxplot(log2(count_data + 1), las = 2, main = 'Raw Counts Distribution', col = 'skyblue')
# 
# # PCA on raw counts
# pca = prcomp(t(log2(count_data + 1)))
# 
# # Create a PCA plot
# plot(pca$x[, 1], pca$x[, 2], xlab = 'PC1', ylab = 'PC2', main = 'PCA of Raw Counts', pch = 19, col = 'dodgerblue')
# text(pca$x[, 1], pca$x[, 2], labels = colnames(count_data), pos = 4, cex = 0.7, col = 'red')
# 
# ####===============Step 4: Normalization ==================================================
# # 1. TMM normalization
# dge = DGEList(counts = count_data)
# dge = calcNormFactors(dge)
# tmm_counts = cpm(dge, normalized.lib.sizes = TRUE)
# 
# # 2. Calculate FPKM
# counts_to_fpkm = function(counts, gene_length, library_sizes) {
#   fpkm = t(t(counts) / library_sizes) * 1e6
#   fpkm = fpkm / (gene_length / 1e3)
#   return(fpkm)
# }
# fpkm_counts = counts_to_fpkm(count_data, gene_lengths, library_sizes)
# 
# ####===============Step 5: Differential Expression Analysis ===============================
# group = factor(c('Control', 'Control', 'Control', 'Pst', 'Pst', 'Pst'))
# dge$samples$group = group
# dge = estimateDisp(dge, design = model.matrix(~ group))
# fit = glmFit(dge, design = model.matrix(~ group))
# lrt = glmLRT(fit)
# 
# results = topTags(lrt, n = Inf)$table
# write.table(results, "R analysis/differential_expression_results_dropped.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# 
# ####===============Step 6: Data Visualization ===================================================
# # MA plot
# plotMD(lrt, main = "MA Plot", xlim = c(-5, 5))
# 
# # Volcano plot
# results$Significant <- results$FDR < 0.05
# ggplot(results, aes(x = logFC, y = -log10(FDR), color = Significant)) +
#   geom_point() +
#   theme_minimal() +
#   labs(title = "Volcano Plot", x = "Log Fold Change", y = "-log10(FDR)")
# 
# # Heatmap of significant genes
# library(pheatmap)
# sig_genes = rownames(results[results$FDR < 0.05, ])
# pheatmap(log2(tmm_counts[sig_genes, ] + 1), cluster_rows = TRUE, cluster_cols = TRUE, main = "Significant Genes Heatmap")
# 
# ####===============Step 7: GO Enrichment Analysis =========================================
# # Load GO mapping data
# gene_go_terms <- read.table("annotation/mart_export.txt", header = TRUE, sep = ",", fill = TRUE, comment.char = "")
# colnames(gene_go_terms)[colnames(gene_go_terms) == 'Gene.stable.ID'] = "GeneID"
# 
# # Create a named list of GO terms for each gene
# go_list = split(gene_go_terms$GO.term.accession, gene_go_terms$GeneID)
# 
# # Input: Gene of interest
# gene_list = sig_genes
# 
# # Perform enrichment analysis using clusterProfiler
# enrichment_result <- enricher(
#   gene = gene_list,
#   TERM2GENE = gene_go_terms[, c("GO.term.accession", "GeneID")]
# )
# 
# 
# # View enrichment results
# head(enrichment_result)
# 
# # Visualize top GO terms
# dotplot(enrichment_result, showCategory = 4)
# 
# # Debugging steps
# overlap_genes <- intersect(gene_list, gene_go_terms$GeneID)
# length(overlap_genes)  # Number of overlapping genes
# if (length(overlap_genes) == 0) {
#   stop("No overlapping genes found between gene_list and TERM2GENE. Check for mismatches.")
# }
# 
# head(gene_go_terms)
# summary(gene_go_terms$GO.term.accession)
# 
# if (length(overlap_genes) > 0) {
#   enrichment_result <- enricher(
#     gene = overlap_genes,
#     TERM2GENE = gene_go_terms[, c("GO.term.accession", "GeneID")]
#   )
# }
# 
# 
# # Gene Set enrichment Analysis
# # Prepare a ranked gene list (e.g., based on log fold change)
# ranked_genes <- results$logFC
# names(ranked_genes) <- rownames(results)
# 
# # Remove NA values (if any)
# ranked_genes <- na.omit(ranked_genes)
# 
# # Sort the list in decreasing order
# ranked_genes <- sort(ranked_genes, decreasing = TRUE)
# 
# # Check the structure of the ranked list
# head(ranked_genes)
# 
# # Perform GSEA
# gsea_result <- GSEA(
#   geneList = ranked_genes,
#   TERM2GENE = gene_go_terms[, c("GO.term.accession", "GeneID")],
#   pvalueCutoff = 0.05
# )
# 
# # View results
# head(gsea_result)
# dotplot(gsea_result, showCategory = 10)
# 
# # Load GO.db
# if (!requireNamespace("GO.db", quietly = TRUE)) {
#   BiocManager::install("GO.db")
# }
# library(GO.db)
# 
# # Map GO IDs to ontology categories
# go_categories <- AnnotationDbi::select(
#   GO.db,
#   keys = gsea_result$ID,
#   columns = "ONTOLOGY",
#   keytype = "GOID"
# )
# 
# # Merge ontology categories with GSEA results
# gsea_result_annotated <- merge(
#   gsea_result,
#   go_categories,
#   by.x = "ID",
#   by.y = "GOID",
#   all.x = TRUE
# )
# 
# # Filter for Molecular Functions (MF)
# mf_terms <- gsea_result_annotated[gsea_result_annotated$ONTOLOGY == "MF", ]
# 
# # View top Molecular Functions
# head(mf_terms)
# 
# 
# library(ggplot2)
# 
# GeneRatio = count / setSize
# 
# 
# # Add GeneRatio to mf_terms
# mf_terms$GeneRatio <- mf_terms$Count / mf_terms$setSize
# 
# # Create a custom dot plot with ggplot2
# ggplot(mf_terms, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
#   geom_point() +
#   scale_color_gradient(low = "blue", high = "red") +
#   theme_minimal() +
#   labs(
#     title = "Top Enriched Molecular Functions",
#     x = "Gene Ratio",
#     y = "GO Term Description",
#     size = "Gene Count",
#     color = "Adjusted P-value"
#   )


###############Case 3 #######-------------------------------------------------------------
# Step 1. Import gene counts into R studio
# 1.1 Import required R-libraries

library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(KEGG.db)
library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)

# Step 1.2: Import featureCounts output
countdata = read.table('~/RNAseq_analysis/RNAseq/RNA_seq_workflow/output/4_final_counts/gene_counts.txt',
                    header = TRUE, # use first row as column names
                    skip = 1, # skip first row
                    row.names = 1, # make row names the gene identifiers
                    #sep = '\t', # tab-separated file
                    #check.names = FALSE # prevent R from altering column names
)
head(countdata)
# Remove char/length columns
count_data = countdata[, 6:ncol(countdata)] # Extract row counts
View(count_data)

# Rename column names of count_data
colnames(count_data) = c("SRR8437484",# control 1
                         "SRR8437485",# control 2
                         "SRR8437482",# control 3
                         "SRR8437483",# control 4
                         "SRR8437480",# Pst 48h 1
                         "SRR8437481",# Pst 48h 2
                         "SRR8437478",# Pst 48h 3
                         "SRR8437479" # Pst 48h 4
)

# Drop `barley_control 2` and `barley_Pst 1`
count_data = count_data[, -c(2, 5)]
# make sure IDs are correct
head(count_data)

# Create metadata
metadata <- data.frame(
  SampleID = c("SRR8437484", "SRR8437482", "SRR8437483",
               "SRR8437481", "SRR8437478", "SRR8437479"),
  Condition = c("Control", "Control", "Control",
                "Pst", "Pst", "Pst"),
  TimePoint = rep("48h", 6) # All collected at 48h
)

# Write to a file
write.table(metadata, "R analysis/metadata_filtered.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Step 1.3: Import metadata text file. The SampleID's must be the first column

metadata = read.delim('R analysis/metadata_filtered.txt', row.names = 1, header = TRUE)

# Add sampleID's to the mapping file
metadata$sampleid <- row.names(metadata)

# Reorder sampleID's to match featureCounts column order. 
metadata <- metadata[match(colnames(count_data), metadata$sampleid), ]
head(metadata)

# Ensure Condition is a factor
metadata$Condition = as.factor(metadata$Condition)

# Step 1.4: Make DESeq2 object from counts and metadata
ddsMat <- DESeqDataSetFromMatrix(countData = count_data,
                                 colData = metadata,
                                 design = ~ Condition)

# # Find differential expressed genes
ddsMat <- DESeq(ddsMat)
# Step 1.5: Get basic statistics about the number of significant genes

# Get results from testing with FDR adjust pvalues
results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05)

# Generate summary of testing. 
summary(results)

mcols(results, use.names = T)

# Subset for only significant genes (q < 0.05)
results_sig <- subset(results, padj < 0.05) # 117 rows 6 columns


# Step 2: Annotate gene symbols
# Using biomaRt
library(biomaRt)

# 1. Connect to Ensembl Plants (contains barley genome annotations)
mart <- useMart(
  biomart = "plants_mart", 
  dataset = "hvulgare_eg_gene", 
  host = "https://plants.ensembl.org"
)

# 2. Fetch GO terms and KEGG pathways 
# 2.1 Extract gene IDs from results_sig1
gene_ids = rownames(results_sig)
head(gene_ids)
length(gene_ids) # 117

# 2.2 Query BioMart for GO terms and KEGG pathways

# Fetch GO terms and KEGG annotations
annotations <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "description",
    "entrezgene_id",
    "go_id",
    "namespace_1003"  # To specify BP, MF, or CC
  ),
  filters = "ensembl_gene_id",
  values = gene_ids,  # List of your significant gene IDs
  mart = mart
)

nrow(results) # 35106
nrow(annotations) # 263

# Match rows of annotations to results
#annotations <- annotations[match(rownames(results), annotations$ensembl_gene_id), ]


# Convert Row Names to a Column in results_sig
results_sig$GeneID <- rownames(results_sig)

# Create a data frame from annotations
annotations_df <- data.frame(
  GeneID = annotations$ensembl_gene_id,
  symbol = annotations$external_gene_name,
  description = annotations$description,
  entrez = annotations$entrezgene_id,
  go_id = annotations$go_id,
  namespace = annotations$namespace_1003,
  stringsAsFactors = FALSE
)

# Check the structure of annotations_df
str(annotations_df)

# Convert results_sig to dataframe
results_sig_df <- as.data.frame(results_sig)

# Perform the merge
results_sig1 <- merge(
  results_sig_df,          # Converted results_sig data frame
  annotations_df,          # Annotation data frame
  by = "GeneID",           # Merge by GeneID
  all.x = TRUE             # Keep all rows from results_sig_df
)


# Check the result
head(results_sig1)
# Save the annotated results to a file
write.table(results_sig1, "R analysis/annotated_significant_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 3. Perform GO/KEGG enrichment (Start from here tomorrow)
library(clusterProfiler)

# 3.1 Prepare TERM2GENE for Enrichment
# Create TERM2GENE data frame
term2gene = annotations[, c('go_id', 'ensembl_gene_id')]
# 263 rows

# head(go_enrichment) # Got 0 rows 

# Resolving the issue
# Verify term2gene formatting
# Check for any empty or malformed GO terms
term2gene_clean <- term2gene[grep("^GO:\\d+$", term2gene$go_id), ]
print(nrow(term2gene_clean))  # Count rows with valid GO IDs

# Validate gene_ids vector
# Check if gene_ids matches term2gene
valid_genes <- intersect(gene_ids, term2gene_clean$ensembl_gene_id)
print(length(valid_genes))  # Count matching genes

term2gene_clean <- term2gene_clean[term2gene_clean$ensembl_gene_id %in% gene_ids, ]
print(nrow(term2gene_clean))  # Rows remaining after filtering

library(clusterProfiler)

go_enrichment <- enricher(
  gene = valid_genes,  # Use filtered genes
  TERM2GENE = term2gene,
  pvalueCutoff = 0.05
)

# Check if enrichment results exist and are non-empty
if (!is.null(go_enrichment) && nrow(as.data.frame(go_enrichment)) > 0) {
  
  # Convert enrichment results to a data frame
  enrichment_results <- as.data.frame(go_enrichment)
  
  # Select relevant columns (e.g., Description, GeneRatio, p-value, adjusted p-value, gene list)
  selected_columns <- enrichment_results[, c("ID", "Description", "GeneRatio", "pvalue", "p.adjust", "geneID")]
  
  # Write to a text file
  write.table(selected_columns, file = "R analysis/GO_enrichment_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
  print("Enrichment results saved to 'GO_enrichment_results.txt'")
} else {
  print("No significant enrichment results to save.")
}

##### Step 4: Write all important results to .txt files
# Save normalized gene counts to a .txt file
write.table(x = as.data.frame(counts(ddsMat), normalized = T), 
            file = 'R analysis/normalized_counts.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Save significant normalized gene counts to a .txt file
write.table(x = as.data.frame(counts(ddsMat[row.names(results_sig)], normalized = T)), 
            file = 'R analysis/normalized_counts_significant.txt', 
            sep = '\t', 
            quote = F, 
            col.names = NA)

# Save DESeq2 results (all genes) to a .txt file
write.table(x = as.data.frame(results), 
            file = "R analysis/results_gene_annotated.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Save significant annotated results table to a .txt file
write.table(x = as.data.frame(results_sig), 
            file = "R analysis/results_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Save GO Enrichment results to a .txt file
write.table(x = as.data.frame(go_enrichment),
            file = "R_analysis/GO_enrichment_resuls.txt",
            sep = "\t",
            quote = F,
            row.names = F)

# Step 5: Plotting gene expression data

# Step 5.1: PCA plot
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Plot PCA by column variable
plotPCA(ddsMat_rlog, intgroup = "Condition", ntop = 50) +
  theme_bw() + # remove default ggplot2 theme
  geom_point(size = 3) + # Increase point size
  scale_y_continuous(limits = c(-10, 10)) + # change limits to fix figure dimensions
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "Top 50 most variable genes") 

# Step 5.2: Heatmap
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)
head(ddsMat_rlog)

# Gather 30 significant genes and make matrix
mat <- assay(ddsMat_rlog[row.names(results_sig)])[1:30, ]

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(ddsMat_rlog)$Condition)
)

# Specify colors you want to annotate the columns by.
ann_colors = list(Group = c(LoGlu = "lightblue", HiGlu = "darkorange"))

# Make Heatmap with pheatmap function.
library(RColorBrewer)
pheatmap(mat = mat, 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = annotation_col, # Add multiple annotations to the samples
         annotation_colors = ann_colors,# Change the default colors of the annotations
         fontsize = 6, # Make fonts smaller
         cellwidth = 50, # Make the cells wider
         cellheight = 10,
         show_colnames = F)

# step 5.3. Volcano Plot
# Gather Log-fold change and FDR-corrected pvalues from DESeq2 results
## - Change pvalues to -log10 (1.3 = 0.05)
data <- data.frame(gene = row.names(results),
                   pval = -log10(results$padj), 
                   lfc = results$log2FoldChange)

# Remove any rows that have NA as an entry
data <- na.omit(data)

# Color the points which are up or down
## If fold-change > 0 and pvalue > 1.3 (Increased significant)
## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
data <- mutate(data, color = case_when(data$lfc > 0 & data$pval > 1.3 ~ "Increased",
                                       data$lfc < 0 & data$pval > 1.3 ~ "Decreased",
                                       data$pval < 1.3 ~ "nonsignificant"))

# Make a basic ggplot2 object with x-y values
vol <- ggplot(data, aes(x = lfc, y = pval, color = color))

# Add ggplot2 layers
vol +   
  ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
  geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "#008B00", Decreased = "#CD4F39", nonsignificant = "darkgray")) +
  theme_bw(base_size = 14) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression(log[2]("LoGlu" / "HiGlu"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
  scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values

# Step 5.4: MA plot (Visualize distribution of RNA-seq data)
plot(results$baseMean, results$log2FoldChange, 
     log = "", ylim = c(-5, 5),
     xlab = "Mean of Normalized Counts",
     ylab = "Log2 Fold Change",
     main = "MA Plot without Log Scaling")

# Step 5.5: Dispersions plot
plotDispEsts(ddsMat)

# Step 5.6: Single gene plot
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Get gene with highest expression
top_gene <- rownames(results)[which.min(results$log2FoldChange)]

# Plot single gene
plotCounts(dds = ddsMat, 
           gene = top_gene, 
           intgroup = "Condition", 
           normalized = T, 
           transform = T,
           col = as.numeric(ddsMat$Condition),
           main = paste("Expression of", top_gene))


# Get list of installed packages
installed_packages = installed.packages()[, "Package"]

# Save package list to a file
writeLines(installed_packages, "requirements_r.txt")

