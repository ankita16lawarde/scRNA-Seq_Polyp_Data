# differential expression analysis between healthy and polyps samples

## endometrium samples scRNA-seq analysis
setwd("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/scRNA-Seq_polyp_data/analysis")


options(future.globals.maxSize = 4000 * 1024^2)


library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

## marker identification
seurat_integrated <- readRDS(file = "clusters_polypdata_integrated_from_0.2.rds")
seurat_integrated


# Assign identity of clusters : 
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

DimPlot(seurat_integrated, 
        label = TRUE)  #+ NoLegend()

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

## merged data
# Cluster name	  CLUSTER NUMBERS 
# Epithelial	    0,5,9
# Stromal	        1,7,10
# Perivascular	  2,19
# Endothelial	    3
# Immune	        4,6,8,15,16,17,18,20
# Macrophage	    11
# Ciliated	      14
# B cells	        13
# Unknown 	        12


## cell type for control and healthy data analysis
#0,6,9,13 Stromal
#1,5,12 immune
#2 Perivascular
#4,14 Epithelial
#7 endothelial
#10 macrophage
#11 ciliated
#8 Others
#3 Excluded

## celltypes for subcell analysis

#Stromal 0
#Perivascular/SMC 1,8
#Immune  2, 6, 14
#Endothelial 3
#Epithelial 4
#Excluded/unknown 5, 7, 10, 13
#Macrophage 9
#B cells 11
#Ciliated 12

# new.cluster.ids <- c("Stromal",
#                      "PerivascularSMC",
#                      "Immune",
#                      "Endothelial",
#                      "Epithelial",
#                      "Excluded",
#                      "Immune",
#                      "Excluded",
#                      "PerivascularSMC",
#                      "Macrophage",
#                      "Excluded",
#                      "Bcells",
#                      "Ciliated",
#                      "Excluded",
#                      "Immune")

new.cluster.ids <- c("Epithelial",
                     "Stromal",
                     "Perivascular",
                     "Endothelial",
                     "Immune",
                     "Epithelial",
                     "Immune",
                     "Stromal",
                     "Immune",
                     "Epithelial",
                     "Unknown1",
                     "Macrophage",
                     "Unknown2",
                     "Bcells",
                     "Ciliated",
                     "Immune",
                     "Immune",
                     "Immune",
                     "Immune",
                     "Perivascular",
                     "Immune")



levels(seurat_integrated)

names(new.cluster.ids) <- levels(seurat_integrated)

seurat_integrated <- RenameIdents(seurat_integrated, new.cluster.ids)
levels(seurat_integrated)

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE, 
        label.size = 4)

dev.off()

DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE, 
        label.size = 4,
        split.by = "sampleType")

dev.off()

#########

seurat_integrated@meta.data
seurat_integrated@active.ident

# Rename "Healthy" to "Endometrium"
levels(seurat_integrated$sampleType)[levels(seurat_integrated$sampleType) == "Healthy"] <- "adEN"
levels(seurat_integrated$sampleType)[levels(seurat_integrated$sampleType) == "Polyps"] <- "EPs"

levels(seurat_integrated$sampleType)
#metadata <- seurat_integrated@meta.data

#metadata$cell_types <- seurat_integrated@active.ident
#head(metadata)

#seurat_integrated@meta.data <- metadata
#seurat_integrated@meta.data


## pseudobuld DEG analysis using DEseq2
#https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html

# Load libraries
#BiocManager::install("apeglm")

library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)


# Extract raw counts and metadata to create SingleCellExperiment object
seurat_integrated
seurat_integrated@assays$RNA@counts
counts <- seurat_integrated@assays$RNA@counts 
#head(counts)

scrna <- JoinLayers(seurat_integrated)
scrna
counts <- LayerData(scrna, assay = "RNA", layer = "counts")
counts
scrna@meta.data

rm(scrna)

metadata <- seurat_integrated@meta.data
seurat_integrated@active.ident
#scrna@active.ident
factor(seurat_integrated@active.ident)

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat_integrated@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)
sce


# Explore the raw counts for the dataset
## Check the assays present
assays(sce)

## Check the counts matrix
dim(counts(sce))
counts(sce)[1:6, 1:6]

#Explore the cellular metadata for the dataset

dim(colData(sce))

head(colData(sce))


##Preparing the single-cell dataset for pseudobulk analysis


# Extract unique names of clusters (= levels of cluster_id factor variable)
levels(colData(sce)$cluster_id)
cluster_names <- levels(colData(sce)$cluster_id)
cluster_names

# Total number of clusters
length(cluster_names)


# Extract unique names of samples (= levels of sample_id factor variable)
colData(sce)$sampleType
#colData(sce)$sampleType <- gsub("_", "", colData(sce)$sampleType)
colData(sce)$sampleType <- factor(colData(sce)$sampleType)
colData(sce)$sampleType


#colData(sce)$DataSet <- gsub("_", "", colData(sce)$DataSet) 
colData(sce)$DataSet <- factor(colData(sce)$DataSet)
colData(sce)$DataSet

levels(colData(sce)$DataSet)
levels(colData(sce)$sampleType)


sample_names <- levels(colData(sce)$DataSet)
sample_names

# Total number of samples
length(sample_names)
colData(sce)

# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups <- colData(sce)[, c("cluster_id", "DataSet")]
head(groups)
groups$cluster_id
groups$DataSet
colnames(groups)

#rm(aggr_counts)

# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                groupings = groups, fun = "sum") 

# Explore output matrix
class(aggr_counts)
dim(aggr_counts)
aggr_counts[1:6, 1:6]

# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts <- t(aggr_counts)
aggr_counts[1:6, 1:6]
dim(aggr_counts)


# Understanding tstrsplit()

## Exploring structure of function output (list)
tstrsplit(colnames(aggr_counts), "_") %>% str()

## Comparing the first 10 elements of our input and output strings
head(colnames(aggr_counts), n = 10)
head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 10)
tstrsplit(colnames(aggr_counts), "_")[[1]]


# Using which() to look up tstrsplit() output
#b_cell_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == "Bcells")
#b_cell_idx

#colnames(aggr_counts)[b_cell_idx]
#aggr_counts[1:10, b_cell_idx]


# As a reminder, we stored our cell types in a vector called cluster_names
cluster_names
cluster_names[20]

# Loop over all cell types to extract corresponding counts, and store information in a list
1:length(cluster_names)
## Initiate empty list
counts_ls <- list()

for (i in 1:length(cluster_names)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
  
}

# Explore the different components of the list
str(counts_ls)
counts_ls

counts_ls[[20]]

# Remove the 20th element (index 19) from the list
counts_ls <- counts_ls[-20]

# Verify the removal
str(counts_ls)


# Reminder: explore structure of metadata
head(colData(sce))

# Extract sample-level variables
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  select(DataSet, sampleType)

dim(metadata)
head(metadata)
metadata

# Exclude duplicated rows
metadata <- metadata[!duplicated(metadata), ]

dim(metadata)
head(metadata)
metadata

# Rename rows
rownames(metadata) <- metadata$DataSet
head(metadata)


# Number of cells per sample and cluster
t <- table(colData(sce)$DataSet,
           colData(sce)$cluster_id)
t[1:6, 1:6]
rownames(t)
colnames(t)
head(t)


# Creating metadata list

## Initiate empty list
metadata_ls <- list()
length(counts_ls)
colnames(counts_ls[[1]])
colnames(t)
head(metadata)
names(metadata) <- c("sample_id", "sampleType")
metadata

counts_ls


for (i in 1:length(counts_ls)) {
  #print(i)
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
  
}

# Explore the different components of the list
str(metadata_ls)
metadata_ls

#Creating a DESeq2 object

# Select cell type of interest
cluster_names



###########################################################################################
##########################################################################################
# Create directories to save results if they don't already exist:
if (!dir.exists("DESeq2")) { dir.create("DESeq2") }
if (!dir.exists("DESeq2/merged_unknown2")) { dir.create("DESeq2/merged_unknown2") }
setwd("DESeq2/merged_unknown2")


# Function to run DESeq2 Wald Test and get results for any cluster:
## clustx is the name of the cluster (cell type) on which to run the function
## A is the sample group to compare (e.g. stimulated condition)
## B is the sample group to compare against (base/control level)
## padj_cutoff defines the ajusted p-value cutoff for significance (set to 0.05 by default)

## This function assumes the counts matrices and metadata for all clusters have been prepared
## and arranged in matching named lists (as illustrated in tutorial above)
## This function assumes the contrast (e.g. stim vs. control) is stored in a variable named "group_id"
metadata

get_dds_resultsAvsB <- function(clustx, A, B, padj_cutoff = 0.05) {
  
  print(clustx) # useful for debugging
  
  # Extract counts matrix and metadata for cluster x
  idx <- which(names(counts_ls) == clustx)
  cluster_counts <- counts_ls[[idx]]
  cluster_metadata <- metadata_ls[[idx]]
  
  cluster_metadata
  # Print error message if sample names do not match
  if ( all(colnames(cluster_counts) != rownames(cluster_metadata)) ) {
    print("ERROR: sample names in counts matrix columns and metadata rows do not match!")
  }
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ sampleType)
  
  # Transform counts for data visualization
  rld <- rlog(dds, blind = TRUE)
  save(rld, file = paste0(idx,"_rlog_transformed_data_HC_scRNA-Seq_data.RData"))
  
  
  # Generate QC plots
  
  ## Plot and save PCA plot
  DESeq2::plotPCA(rld, intgroup = "sampleType")
  if (!dir.exists("results")) { dir.create("results") }
  ggsave(paste0("results/", clustx, "_specific_PCAplot.png"))
  
  ## Extract rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  ## Plot and save heatmap
  png(paste0("results/", clustx, "_specific_heatmap.png"),
      height = 6, width = 7.5, units = "in", res = 300)
  pheatmap(rld_cor, annotation = cluster_metadata[, c("sampleType"), drop = FALSE])
  dev.off()
  
  
  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  
  ## Plot dispersion estimates
  png(paste0("results/", clustx, "_dispersion_plot.png"),
      height = 5, width = 6, units = "in", res = 300)
  plotDispEsts(dds)
  dev.off()
  
  ## Output and shrink results of Wald test for contrast A vs B
  contrast <- paste(c("sampleType", "Polyps", "vs", "Healthy"), collapse = "_")
  print(resultsNames(dds))
  
  res <- results(dds, name = contrast, alpha = 0.05)
  res <- lfcShrink(dds, coef = contrast, res = res)
  
  ## Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  
  write.csv(res_tbl,
            paste0("results/", clustx, "_", contrast, "_all_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  ## Subset the significant results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  write.csv(sig_res,
            paste0("results/", clustx, "_", contrast, "_signif_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  
  # Generate results visualization plots
  
  ## Extract normalized counts from dds object
  normalized_counts <- counts(dds, normalized = TRUE)
  save(normalized_counts, file = paste0(idx,"_normalized_data_HC_scRNA-Seq_data.RData"))
  
  ## Extract top 20 DEG from resLFC (make sure to order by padj)
  top20_sig_genes <- sig_res %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene) %>%
    head(n = 20)
  
  ## Extract matching normalized count values from matrix
  top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
  
  ## Convert wide matrix to long data frame for ggplot2
  top20_sig_df <- data.frame(top20_sig_counts)
  top20_sig_df$gene <- rownames(top20_sig_counts)
  
  top20_sig_df <- melt(setDT(top20_sig_df), 
                       id.vars = c("gene"),
                       variable.name = "cluster_sample_id") %>% 
    data.frame()
  
  #head(top20_sig_df)
  ## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
  top20_sig_df$cluster_sample_id <- gsub("\\X", "", top20_sig_df$cluster_sample_id)
  #top20_sig_df$cluster_sample_id <- gsub("\\  ", "", top20_sig_df$cluster_sample_id)
  
  ## Join counts data frame with metadata
  top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds)),
                             by = "cluster_sample_id")
  
  ## Generate plot
  ggplot(top20_sig_df, aes(y = value, x = sampleType, col = sampleType)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6) +
    geom_jitter(height = 0, width = 0.15, size = 2, alpha = 0.6) +  # ⬅️ Increased point size
    scale_y_continuous(trans = 'log10') +
    ylab("log10 of normalized expression level") +
    #ggtitle("Top 20 Significant DE Genes") +
    theme_bw(base_size = 14) +  # ⬅️ Base font size for theme
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 13, face = "bold"),  # facet labels
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
      panel.border = element_blank(),
      panel.grid.major = element_blank()
    ) +
    scale_x_discrete(labels = c("Polyps" = "EPs", "Healthy" = "adEN")) +
    scale_color_manual(
      name = "Sample Type",
      values = c("Polyps" = "#E41A1C", "Healthy" = "#377EB8"),
      labels = c("Polyps" = "EPs", "Healthy" = "adEN")
    ) +
    facet_wrap(~ gene, scales = "free_y")
  
  
  ggsave(paste0("results/", clustx, "_", contrast, "_top20_DE_genes.png"))
  
}

resultsNames(dds)
#get_dds_resultsAvsB("Bcell", A = "Peritoniallesion", B = "Endometrium")

# Run the script on all clusters comparing stimulated condition relative to control condition
#map("Bcells", get_dds_resultsAvsB, A = "Healthy", B = "Control")

map("Perivascular", get_dds_resultsAvsB, A = "Polyps", B = "Healthy")
map("Unknown1", get_dds_resultsAvsB, A = "Polyps", B = "Healthy")
