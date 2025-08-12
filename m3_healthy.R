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
#seurat_integrated@meta.data

# Assign identity of clusters : 
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

table(seurat_integrated@meta.data$DataSet)


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

table(seurat_integrated@meta.data$sampleType, seurat_integrated@meta.data$celltypes)

saveRDS(seurat_integrated, file = "seurat_integrated_with_celltype.rds")

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE, 
        label.size = 4)
dev.off()


### polyp samples trajectory
# Create a Seurat object for condition "B"
seurat_healthy <- subset(seurat_integrated, subset = sampleType == "Healthy")

DimPlot(seurat_healthy,
        reduction = "umap",
        label = TRUE, 
        label.size = 4)
dev.off()

#seurat_healthy@meta.data$celltypes <- seurat_healthy@active.ident
seurat_healthy@meta.data$celltypes
levels(seurat_healthy)
seurat_healthy@meta.data$celltypes <- seurat_healthy@active.ident


############################################################
## remove the immune, b cells and machogpage cell clusters 
############################################################

# Remove Immune, Bcells, and Macrophage clusters
seurat_healthy_filtered <- subset(seurat_healthy, idents = c("Immune", "Bcells", "Macrophage"), invert = TRUE)

# Check remaining clusters
levels(seurat_healthy_filtered)

seurat_healthy_filtered@meta.data$celltypes


# Convert Seurat object to a CellDataSet
#cds_healthy <- as.cell_data_set(seurat_healthy)

cds_healthy <- as.cell_data_set(seurat_healthy_filtered)

cds_healthy

head(colData(cds_healthy))
fData(cds_healthy)
rownames(fData(cds_healthy))[1:10]

fData(cds_healthy)$gene_short_name <- rownames(fData(cds_healthy))
head(fData(cds_healthy))

head(counts(cds_healthy))

#1. Assign partitions
length(cds_healthy@colData@rownames)
cds_healthy@colData@rownames
recreate.partitions <- c(rep(1, length(cds_healthy@colData@rownames)))
names(recreate.partitions) <- cds_healthy@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds_healthy@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

#Assign cluster information
list.cluster <- seurat_healthy@active.ident
cds_healthy@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
#cds_healthy@clusters@listData[["PCA"]][["clusters"]] <- list.cluster

# seurat_healthy@reductions$umap@cell.embeddings
# seurat_healthy@reductions$pca@cell.embeddings
#Assign UMAP coordinates
cds_healthy@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seurat_healthy_filtered@reductions$umap@cell.embeddings
cds_healthy

#cds_healthy@int_colData@listData[["reducedDims"]]@listData[["PCA"]] <- seurat_healthy_filtered@reductions$pca@cell.embeddings

cds_healthy@clusters
colData(cds_healthy)$celltypes

cluster.before.traj <-plot_cells(cds_healthy, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                 group_label_size = 5) + theme(legend.position = "right")

cluster.before.traj

# Cluster cells if needed (Monocle3 can cluster cells based on UMAP or PCA)
#set.seed(45)
cds_healthy <- cluster_cells(cds_healthy)

# Transfer Seurat cluster identities to Monocle3
#cds_healthy$cell_clusters <- seurat_healthy$celltypes

# Learn the trajectory graph
cds_healthy <- learn_graph(cds_healthy, use_partition = F)

plot_cells(cds_healthy, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5, label_principal_points = F)

plot_cells(cds_healthy, color_cells_by = "celltypes", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5, label_principal_points = F)

monocle_graph_healthy <- "m3_graph_healthy_no_immune.jpeg"

jpeg(monocle_graph_healthy, width = 12, height = 10, units = "in", res = 600)


plot_cells(cds_healthy, color_cells_by = "celltypes", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5, label_principal_points = F)

dev.off()




colnames(cds_healthy)
unique(clusters(cds_healthy))

get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                             (which.max(table(closest_vertex[cell_ids,]))))]       
  root_pr_nodes
}

DN_node_id = get_correct_root_state(cds_healthy, cell_phenotype = 'celltypes', "Stromal")

DN_node_id

get_correct_root_state(cds_healthy, cell_phenotype = 'celltypes', "Ciliated")

get_correct_root_state(cds_healthy, cell_phenotype = 'celltypes', "Endothelial")


cds_healthy <- order_cells(cds_healthy, reduction_method = "UMAP",  
                           root_pr_nodes = DN_node_id)

#cds_healthy <- order_cells(cds_healthy, reduction_method = "PCA", root_cells = colnames(cds_healthy[, clusters(cds_healthy) == "Stromal"]))
dev.off()

plot_psuedo <- "m3_pseudotime_healthy_noimmune.jpeg"
jpeg(plot_psuedo, width = 12, height = 10, units = "in", res = 600)


plot_cells(cds_healthy, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)

dev.off()


head(cds_healthy@int_colData$reducedDims$UMAP)
cds_healthy@colData@rownames
cds_healthy@reducedDims$UMAP[ , "partition"]
# Extract branch information
branch_info <- cds_healthy@reducedDims$UMAP[ , "partition"]  # The partition corresponds to the branch ID
cds_healthy$branch_id <- branch_info


#cds_healthy2 <- order_cells(cds_healthy, reduction_method = "UMAP")

#plot_cells(cds_healthy2, color_cells_by = "pseudotime", label_groups_by_cluster = T,
#           label_branch_points = T, label_roots = F, label_leaves = F)


head(pseudotime(cds_healthy), 10)

cds_healthy$monocle3_pseudotime <- pseudotime(cds_healthy)
data.pseudo <- as.data.frame(colData(cds_healthy))
data.pseudo

ggplot(data.pseudo, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot()

ggplot(data.pseudo, aes(monocle3_pseudotime, celltypes, fill = celltypes)) + geom_boxplot()
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(celltypes, monocle3_pseudotime), fill = celltypes)) + geom_boxplot()

cds_healthy2$monocle3_pseudotime <- pseudotime(cds_healthy2)
data.pseudo2 <- as.data.frame(colData(cds_healthy2))
data.pseudo2

ggplot(data.pseudo2, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()
ggplot(data.pseudo2, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot()

ggplot(data.pseudo2, aes(monocle3_pseudotime, celltypes, fill = celltypes)) + geom_boxplot()
ggplot(data.pseudo2, aes(monocle3_pseudotime, reorder(celltypes, monocle3_pseudotime), fill = celltypes)) + geom_boxplot()

save(cds_healthy, file = "cds_healthy_M3_no_immune.RData")

deg <- graph_test(cds_healthy, neighbor_graph = "principal_graph")
deg

writexl::write_xlsx(deg, path = "DEG_along_trajectory_Healthy_noimmune.xlsx")

save(deg, file = "deg_M3_stromal_start_healthy_no_immune.RData")
deg %>% arrange(q_value) %>% filter(status == "OK") %>% head()

FeaturePlot(seurat_healthy_filtered, features = c("HES4", "TTLL10", "TNFRSF18", "TNFRF4"))

#Add pseudotime values into the seuratobject
seurat_healthy$pseudotime <- pseudotime(cds_healthy)
FeaturePlot(seurat_healthy, features = "pseudotime")

# make the subset CDS
test_genes=c("HES4", "TTLL10", "TNFRSF18", "TNFRF4")
cds_subset1 <- cds_healthy[rowData(cds_healthy)$gene_short_name %in% test_genes,]

# produce violin plots
plot_genes_violin(cds_subset1, group_cells_by="celltypes", ncol=2)

plot_genes_in_pseudotime(cds_subset1, color_cells_by="monocle3_pseudotime")

p1<-plot_genes_in_pseudotime(cds_subset1,color_cells_by="celltypes") 
min(p1$data$expectation)
p2<-plot_genes_in_pseudotime(cds_subset1,color_cells_by="celltypes",min_expr=min(p1$data$expectation)) 
p2


a <- cds_subset1
cds_exprs <- SingleCellExperiment::counts(a)
cds_exprs <- Matrix::t(counts(a))/size_factors(a)
p2$data$expression<-as.numeric(cds_exprs)
p2

########cds_subset1
deg_knn <- graph_test(cds_healthy, neighbor_graph = "knn")
deg_knn
save(deg_knn, file = "deg_knn_M3_stromal_start_healthy.RData")
deg_knn %>% arrange(q_value) %>% filter(status == "OK") %>% head()

FeaturePlot(seurat_healthy, features = c("HES4", "TTLL10", "TNFRSF18", "TNFRF4"))

#Add pseudotime values into the seuratobject
pseudotime(cds_healthy)
seurat_healthy$pseudotime <- pseudotime(cds_healthy)
FeaturePlot(seurat_healthy, features = "pseudotime")

# make the subset CDS
test_genes=c("HES4", "TTLL10", "TNFRSF18", "TNFRF4")
cds_subset1 <- cds_healthy[rowData(cds_healthy)$gene_short_name %in% test_genes,]

# produce violin plots
plot_genes_violin(cds_subset1, group_cells_by="celltypes", ncol=2)

plot_genes_in_pseudotime(cds_subset1, color_cells_by="celltypes", min_expr=0.5)


######################
# a helper function to identify the root principal points
get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                             (which.max(table(closest_vertex[cell_ids,]))))]       
  root_pr_nodes
}
# call the function to automatically find the node in the principal graph where our DN cells reside
DN_node_id = get_correct_root_state(cds_healthy, cell_phenotype = 'celltypes', "Stromal")

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


##### with default parameters

cds_healthy1 <- learn_graph(cds_healthy, use_partition = TRUE)

monocle_graph <- "m3_graph_healthy_unknown_1_2.jpeg"

# Open a JPEG device with 600 DPI resolution
jpeg(monocle_graph, width = 12, height = 10, units = "in", res = 300)


plot_cells(cds_healthy1, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)

dev.off()

plot_cells(cds_healthy1, color_cells_by = "partition", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)

principal_points <- "m3_princ_pt_healthy_unknown_1_2.jpeg"

# Open a JPEG device with 600 DPI resolution
jpeg(principal_points, width = 12, height = 10, units = "in", res = 300)

plot_cells(cds_healthy1,
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
stromal_node_id = get_correct_root_state(cds_healthy1, 
                                         cell_phenotype = 'celltypes', "Stromal")

stromal_node_id      # check the node found ## "Y_190"



cds_healthy1 <- order_cells(cds_healthy1, reduction_method = "UMAP", 
                            root_pr_nodes = stromal_node_id)

pseudo_healthy <- "m3_pseudotime_healthy_unknown_1_2.jpeg"

# Open a JPEG device with 600 DPI resolution
jpeg(pseudo_healthy, width = 12, height = 10, units = "in", res = 300)


plot_cells(cds_healthy1, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)

dev.off()

# access pseudotime calculated for each cell and store it alongside cell metadata
pseudotime <- pseudotime(cds_healthy1) 	
cds_healthy1@colData$pseudotime <- pseudotime 


# make the subset CDS
test_genes = c("HES4", "TTLL10", "TNFRSF18", "TNFRF4")
cds_subset2 <- cds_healthy1[rowData(cds_healthy1)$gene_short_name %in% test_genes,]

# produce violin plots
plot_genes_violin(cds_subset2, group_cells_by="celltypes", ncol=2)

plot_genes_in_pseudotime(cds_subset2)

p1 <- plot_genes_in_pseudotime(cds_subset2,color_cells_by="celltypes") 
min(p1$data$expectation)
p2<-plot_genes_in_pseudotime(cds_subset2,color_cells_by="celltypes",min_expr=min(p1$data$expectation)) 
p2


a <- cds_subset2
cds_exprs <- SingleCellExperiment::counts(a)
cds_exprs <- Matrix::t(counts(a))/size_factors(a)
p2$data$expression<-as.numeric(cds_exprs)
p2


head(pseudotime(cds_healthy1), 10)

cds_healthy1$monocle3_pseudotime <- pseudotime(cds_healthy1)
data.pseudo <- as.data.frame(colData(cds_healthy1))
data.pseudo

ggplot(data.pseudo, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot()

ggplot(data.pseudo, aes(monocle3_pseudotime, celltypes, fill = celltypes)) + geom_boxplot()

pseudo_healthy_order <- "m3_pseudotime_healthy_cellorder_U1_2.jpeg"
# Open a JPEG device with 600 DPI resolution
jpeg(pseudo_healthy_order, width = 12, height = 10, units = "in", res = 300)

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(celltypes, monocle3_pseudotime), fill = celltypes)) + geom_boxplot()
dev.off()

save(cds_healthy1, file = "cds_healthy_M3_single_start_Y190.RData")

deg_Y212 <- graph_test(cds_healthy1, neighbor_graph = "principal_graph")
deg_Y212
save(deg_Y212, file = "deg_Y190_M3_stromal_start_healthy_U12.RData")
deg_Y212 %>% arrange(q_value) %>% filter(status == "OK") %>% head()
write.csv(deg_Y212, file = "deg_Y190_stromal_healthy_u12.csv")

save(cds_healthy1, file = "cds_healthy_Y190_start_trajectory_u12.RData")

FeaturePlot(seurat_healthy, features = c("HES4", "TTLL10", "TNFRSF18", "TNFRF4"))

