# MONOCLE analysis
install.packages("devtools")
devtools::install_github("cole-trapnell-lab/monocle-release@develop")

remotes::install_github("satijalab/seurat-wrappers")

setwd("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/scRNA-Seq_polyp_data/analysis")

options(future.globals.maxSize = 4000 * 1024^2)

#BiocManager::install("TSCAN")

library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(scater)
library(monocle3)
library(patchwork)
library(SeuratWrappers)
library(ggrastr)
library(terra)

## load the saved seurat object
## marker identification
seurat_integrated <- readRDS(file = "clusters_polypdata_integrated_from_0.2.rds")
seurat_integrated@meta.data

# Assign identity of clusters : 
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

DimPlot(seurat_integrated, 
        label = TRUE)  #+ NoLegend()

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

## merged data
# Cluster name	  CLUSTER NUMBERS 
# Epithelial	    0,5,9
# Stromal	        1,7
# Perivascular	  2,19
# Endothelial	    3
# Immune	        4,6,8,15,16,17,18,20
# Macrophage	    11
# Ciliated	      14
# B cells	        13
# Unknown 	      10,  12


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

seurat_integrated@meta.data$celltypes <- seurat_integrated@active.ident

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE, 
        label.size = 4)
dev.off()

## split the seurat data into two seurat object, according to the condition
table(seurat_integrated@meta.data$sampleType)

### polyp samples trajectory
# Create a Seurat object for condition "B"
seurat_polyps <- subset(seurat_integrated, subset = sampleType == "Polyps")

# Remove Immune, Bcells, and Macrophage clusters
seurat_polyps_filtered <- subset(seurat_polyps, idents = c("Immune", "Bcells", "Macrophage"), invert = TRUE)

rm(seurat_polyps)

DimPlot(seurat_polyps,
        reduction = "umap",
        label = TRUE, 
        label.size = 4)
dev.off()


DimPlot(seurat_polyps,
        reduction = "pca",
        label = TRUE, 
        label.size = 4)


#seurat_polyps_filtered@active.ident

#seurat_healthy@meta.data$celltypes <- seurat_healthy@active.ident
seurat_polyps_filtered@meta.data$celltypes <- seurat_polyps_filtered@active.ident

# Convert Seurat object to a CellDataSet
cds_polyps <- as.cell_data_set(seurat_polyps_filtered)

cds_polyps

head(colData(cds_polyps))
fData(cds_polyps)
rownames(fData(cds_polyps))[1:10]

fData(cds_polyps)$gene_short_name <- rownames(fData(cds_polyps))
head(fData(cds_polyps))

head(counts(cds_polyps))

#1. Assign partitions
length(cds_polyps@colData@rownames)
cds_polyps@colData@rownames
recreate.partitions <- c(rep(1, length(cds_polyps@colData@rownames)))
names(recreate.partitions) <- cds_polyps@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds_polyps@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

#Assign cluster information
list.cluster <- seurat_polyps_filtered@active.ident
cds_polyps@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

cds_polyps
head(colData(cds_polyps))

#Assign UMAP coordinates
cds_polyps@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seurat_polyps_filtered@reductions$umap@cell.embeddings

cluster.before.traj <-plot_cells(cds_polyps, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                 group_label_size = 5) + theme(legend.position = "right")

cluster.before.traj


cds_polyps1 <- learn_graph(cds_polyps, use_partition = TRUE)

monocle_graph <- "m3_graph_polyps_unknown_1_2_noimmune.jpeg"

# Open a JPEG device with 600 DPI resolution
jpeg(monocle_graph, width = 12, height = 10, units = "in", res = 300)


plot_cells(cds_polyps1, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)
dev.off()


plot_cells(cds_polyps1, color_cells_by = "partition", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)

principal_points <- "m3_princ_pt_polyps_unknown_1_2.jpeg"

# Open a JPEG device with 600 DPI resolution
jpeg(principal_points, width = 12, height = 10, units = "in", res = 300)

plot_cells(cds_polyps1,
           color_cells_by = "celltypes",
           label_cell_groups=FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_principal_points = TRUE,       # set this to TRUE
           graph_label_size=3)

dev.off()

get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                             (which.max(table(closest_vertex[cell_ids,]))))]       
  root_pr_nodes
}

# call the function to automatically find the node in the principal graph where our DN cells reside
stromal_node_id = get_correct_root_state(cds_polyps1, 
                                         cell_phenotype = 'celltypes', "Stromal")

stromal_node_id      # check the node found ## "Y_85"



cds_polyps1 <- order_cells(cds_polyps1, reduction_method = "UMAP", 
                            root_pr_nodes = stromal_node_id)


pseudo_polyp <- "m3_pseudotime_polyps_unknown_1_2_noimmune.jpeg"

# Open a JPEG device with 600 DPI resolution
jpeg(pseudo_polyp, width = 12, height = 10, units = "in", res = 300)


plot_cells(cds_polyps1, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)


dev.off()

# access pseudotime calculated for each cell and store it alongside cell metadata
pseudotime <- pseudotime(cds_polyps1) 	
cds_polyps1@colData$pseudotime <- pseudotime 


# make the subset CDS
test_genes = c("HES4", "TTLL10", "TNFRSF18", "TNFRF4")
cds_subset2 <- cds_polyps1[rowData(cds_polyps1)$gene_short_name %in% test_genes,]

# produce violin plots
plot_genes_violin(cds_subset2, group_cells_by="celltypes", ncol=2)

plot_genes_in_pseudotime(cds_subset2)

p1 <- plot_genes_in_pseudotime(cds_subset2,color_cells_by="celltypes") 
min(p1$data$expectation)
p2 <-plot_genes_in_pseudotime(cds_subset2,color_cells_by="celltypes",min_expr=min(p1$data$expectation)) 
p2


a <- cds_subset2
cds_exprs <- SingleCellExperiment::counts(a)
cds_exprs <- Matrix::t(counts(a))/size_factors(a)
p2$data$expression<-as.numeric(cds_exprs)
p2


head(pseudotime(cds_polyps1), 10)

cds_polyps1$monocle3_pseudotime <- pseudotime(cds_polyps1)
data.pseudo <- as.data.frame(colData(cds_polyps1))
data.pseudo

ggplot(data.pseudo, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot()

ggplot(data.pseudo, aes(monocle3_pseudotime, celltypes, fill = celltypes)) + geom_boxplot()

pseudo_polyp_order <- "m3_pseudotime_polyps_cellorder_U1_2.jpeg"
# Open a JPEG device with 600 DPI resolution
jpeg(pseudo_polyp_order, width = 12, height = 10, units = "in", res = 300)

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(celltypes, monocle3_pseudotime), fill = celltypes)) + geom_boxplot()
dev.off()

save(cds_polyps1, file = "cds_polyps_M3_single_start_Y85_U1_2.RData")

deg_Y95 <- graph_test(cds_polyps1, neighbor_graph = "principal_graph")
deg_Y95
#save(deg_Y95, file = "M3_stromal_start_polyps_U1_2.RData")
deg_Y95 %>% arrange(q_value) %>% filter(status == "OK") %>% head()
write.csv(deg_Y95, file = "DEG_polyps_U1_2_noommune.csv")

save(cds_polyps1, file = "cds_polyps_Y85_start_trajectory_u12.RData")

FeaturePlot(seurat_polyps, features = c("HES4", "TTLL10", "TNFRSF18", "TNFRF4"))


####################################################################################

# Cluster cells if needed (Monocle3 can cluster cells based on UMAP or PCA)
#set.seed(45)
#cds_polyps <- cluster_cells(cds_polyps, reduction_method = "UMAP", clsuter_method = "louvain")

# Transfer Seurat cluster identities to Monocle3
#cds_polyps$cell_clusters <- seurat_polyps$celltypes


get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                             (which.max(table(closest_vertex[cell_ids,]))))]       
  root_pr_nodes
}
# call the function to automatically find the node in the principal graph where our DN cells reside
DN_node_id = get_correct_root_state(cds_polyps, cell_phenotype = 'cluster', "Stromal")

DN_node_id      # check the node found

# specifying root cells: `root_pr_nodes` argument - check the principal points
plot_cells(cds_healthy,
           color_cells_by = "celltypes",
           label_cell_groups=FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_principal_points = TRUE,       # set this to TRUE
           graph_label_size=3)



# Learn the trajectory graph
cds_polyps <- learn_graph(cds_polyps, use_partition = FALSE)

plot_cells(cds_polyps, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)

cds_polyps <- order_cells(cds_polyps, reduction_method = "UMAP", root_cells = colnames(cds_polyps[, clusters(cds_polyps) == "Stromal"]))
plot_cells(cds_polyps, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)

cds_polyps2 <- order_cells(cds_polyps, reduction_method = "UMAP")

plot_cells(cds_polyps2, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)


head(pseudotime(cds_polyps), 10)
cds_polyps$monocle3_pseudotime <- pseudotime(cds_polyps)
data.pseudo_polyps <- as.data.frame(colData(cds_polyps))
data.pseudo_polyps

ggplot(data.pseudo_polyps, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()
ggplot(data.pseudo_polyps, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot()

ggplot(data.pseudo_polyps, aes(monocle3_pseudotime, celltypes, fill = celltypes)) + geom_boxplot()
ggplot(data.pseudo_polyps, aes(monocle3_pseudotime, reorder(celltypes, monocle3_pseudotime), fill = celltypes)) + geom_boxplot()

cds_polyps2$monocle3_pseudotime <- pseudotime(cds_polyps2)
data.pseudo2 <- as.data.frame(colData(cds_polyps2))
data.pseudo2

ggplot(data.pseudo2, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()
ggplot(data.pseudo2, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot()

ggplot(data.pseudo2, aes(monocle3_pseudotime, celltypes, fill = celltypes)) + geom_boxplot()
ggplot(data.pseudo2, aes(monocle3_pseudotime, reorder(celltypes, monocle3_pseudotime), fill = celltypes)) + geom_boxplot()

save(cds_polyps, file = "cds_polyps_M3.RData")

deg <- graph_test(cds_polyps, neighbor_graph = "principal_graph")
deg
save(deg, file = "deg_M3_stromal_start_polyps.RData")
deg %>% arrange(q_value) %>% filter(status == "OK") %>% head()

FeaturePlot(seurat_polyps, features = c("HES4", "TTLL10", "TNFRSF18", "TNFRF4"))

#Add pseudotime values into the seuratobject
seurat_polyps$pseudotime <- pseudotime(cds_polyps)
FeaturePlot(seurat_polyps, features = "pseudotime")

# make the subset CDS
test_genes=c("HES4", "TTLL10", "TNFRSF18", "TNFRF4")
cds_subset1 <- cds_polyps[rowData(cds_polyps)$gene_short_name %in% test_genes,]

# produce violin plots
plot_genes_violin(cds_subset1, group_cells_by="celltypes", ncol=2)

plot_genes_in_pseudotime(cds_subset1, color_cells_by="celltypes", min_expr=0.5)

# get gene IDs for which q-value < 0.05
pr_deg_ids <- row.names(subset(deg, q_value < 0.05))
pr_deg_ids

# group genes into modules
gene_module_df <- find_gene_modules(cds_polyps[pr_deg_ids,], resolution=1e-2) 	

####################################################
# Plot the graph to visualize cell clustering and trajectory
plot_cells(cds_polyps, color_cells_by = "cluster", label_principal_points = TRUE)  # Use appropriate metadata for coloring

# Plot UMAP in Monocle3
plot_cells(cds_polyps, color_cells_by = "ident", reduction_method = "UMAP")

# Plot PCA in Monocle3
plot_cells(cds_polyps, color_cells_by = "ident", reduction_method = "PCA")


# Order cells by pseudotime
cds_polyps <- order_cells(cds_polyps)

# Plot the cells along the trajectory, colored by pseudotime
plot_cells(cds_polyps, color_cells_by = "pseudotime", label_cell_groups = TRUE, label_leaves = TRUE)

# Perform graph-based differential expression analysis along pseudotime
de_genes_polyps <- graph_test(cds_polyps, cores = 4)

de_genes_polyps

# Check the results
head(de_genes)

save(de_genes_polyps, file = "de_genes_polyps_data.RData")


de_genes_polyps_ordered <- de_genes_polyps[order(de_genes_polyps$q_value), ]

head(de_genes_polyps_ordered)

rownames(de_genes_polyps_ordered)[1:10]

rowData(cds_polyps)
rowData(cds_polyps)$gene_short_name <- row.names(rowData(cds_polyps))

# Plot genes with significant differential expression along pseudotime
plot_genes_by_group(cds_polyps, markers = rownames(de_genes_polyps_ordered)[1:50], 
                    group_cells_by = "cell_clusters")

colData(cds)
gene_fits_polyps <- fit_models(cds_polyps, model_formula_str = "~cell_clusters")
gene_fits

fit_coefs <- coefficient_table(gene_fits)
fit_coefs$term

sampleType_terms <- fit_coefs %>% filter(term == "sampleTypePolyps")
sampleType_terms

sampleType_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)


cds_subset <- cds[row.names(subset(rowData(cds),
                                   gene_short_name %in% c("SAMD11", "HES4", "AGRN", "TTLL10"))),]

plot_genes_violin(cds_subset, group_cells_by="cell_clusters", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

plot_genes_in_pseudotime(cds_subset)

# Plot cells colored by gene expression
plot_cells(cds, color_cells_by = "cluster", genes = c("SAMD11", "HES4", "AGRN", "TTLL10"))

evaluate_fits(gene_fits)


cds_polyps_test_res <- graph_test(cds_polyps, neighbor_graph="principal_graph", cores=4)
cds_polyps_test_res

save(cds_polyps_test_res, file = "cds_polyps_test_res.RData")
load(file = "cds_polyps_test_res.RData")

cds_polyps_test_res$gene_short_name <- rownames(cds_polyps_test_res)
writexl::write_xlsx(cds_polyps_test_res, path = "Diff_genes_along_trajectory_polyps.xlsx")

#save(cds, file = "cds_polys_healthy.RData")

load("cds_polys_healthy.RData")
load("cds_per_test_res.RData")

de_genes_ordered_cds <- cds_polyps_test_res[order(cds_polyps_test_res$q_value), ]
de_genes_ordered_cds

plot_cells(cds_polyps, genes=c("FAAP20", "PRKCZ"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)


pr_deg_ids <- row.names(subset(cds_polyps_test_res, q_value < 0.05))
pr_deg_ids

#gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))
colData(cds)
traject_genes <- c("FAAP20", "PRKCZ")
traject_lineage_cds <- cds_polyps[rowData(cds_polyps)$gene_short_name %in% traject_genes]
traject_lineage_cds <- order_cells(traject_lineage_cds)
traject_lineage_cds

colData(traject_lineage_cds)

#colData(cds)$color_cells_by <- col_date(cds)$cell_cluster

plot_genes_in_pseudotime(traject_lineage_cds,
                         color_cells_by ="cell_clusters",
                         min_expr=0.5)


cds_subset <- choose_cells(cds)
cds_subset
subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
subset_pr_test_res
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
pr_deg_ids

traject_genes <- c("HES4",  "ISG15", "AGRN" , "RNF223")
traject_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% traject_genes]
traject_lineage_cds <- order_cells(traject_lineage_cds)
traject_lineage_cds

plot_genes_in_pseudotime(traject_lineage_cds,
                         color_cells_by ="cell_clusters",
                         min_expr=0.5)



