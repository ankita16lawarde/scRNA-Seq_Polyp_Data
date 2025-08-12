## trajectory analysis and differential expression in trajectory

## polyp data

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
library(TSCAN)


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
# Unknown1 	      10
#Unknown2         12


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

## custum bubble plots
features <- c("MTCL1","KIAA1324", "SMAD9","CPM", "CTNNA2",
              "ECM1", "CPXM1", "IGF1", "TOP2A",
              "MUSTN1", "MCAM", "MYLK",
              "NR5A2", "SLCO2A1", "EMCN", "CD34",
              "CPVL", "CSF1R", "CD68",
              "ARMC3", "DNAH7",
              "IL7R", "GZMK", "IL32", "CD3D",
              "BANK1", "CD37", "BLNK")

dev.off()
# Open a PNG device with the specified filename, dimensions, and resolution
png(filename = "DotPlot_Seurat.png", 
    width = 10,    # Width in inches
    height = 8,    # Height in inches
    units = "in",  # Unit of measurement
    res = 600)     # Resolution of the image

# Generate the DotPlot with modified y-axis label
DotPlot(seurat_integrated, 
        features = features, 
        cols = c("lightgrey", "darkblue")) +
  RotatedAxis() + 
  ylab("Cell types") +      # Replace with your desired label
  xlab("Marker genes")               # Optional: Customize x-axis label too
  

# Close the device to write the plot to file
dev.off()


##### trajectory analysis for Polyps samples

# Create a Seurat object for condition "B"
seurat_polyps <- subset(seurat_integrated, subset = sampleType == "Healthy")

DimPlot(seurat_polyps,
        reduction = "umap",
        label = TRUE, 
        label.size = 4)
dev.off()

#seurat_healthy@meta.data$celltypes <- seurat_healthy@active.ident
#seurat_polyps@meta.data$celltypes <- seurat_polyps@active.ident

# Set up metadata as desired for aggregation and DE analysis
seurat_polyps@meta.data$celltypes <- factor(seurat_polyps@active.ident)

# Remove Immune, Bcells, and Macrophage clusters
seurat_polyps_filtered <- subset(seurat_polyps, idents = c("Immune", "Bcells", "Macrophage"), invert = TRUE)

rm(seurat_polyps)
rm(seurat_integrated)

# Create single cell experiment object
sce_polyps <- as.SingleCellExperiment(seurat_polyps_filtered)
sce_polyps

seurat_polyps[['umap']]@cell.embeddings
reducedDim(sce_polyps, "PCA", withDimnames = TRUE) <- seurat_polyps_filtered[['pca']]@cell.embeddings
reducedDim(sce_polyps, "UMAP", withDimnames = TRUE) <- seurat_polyps_filtered[['umap']]@cell.embeddings

colData(sce_polyps)$celltypes

## basic steps
by.cluster <- aggregateAcrossCells(sce_polyps, ids=colData(sce_polyps)$celltypes)
by.cluster
centroids <- reducedDim(by.cluster, "PCA")
centroids

# Set clusters=NULL as we have already aggregated above.
library(TSCAN)
mst <- createClusterMST(centroids, clusters=NULL)
mst

line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="UMAP")
line.data

##############################################################################
# Load libraries
library(scater)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

# -----------------------
# Step 1: Cell type color map (for UMAP points)
# -----------------------
# celltype_colors <- c(
#   "Epithelial"    = "#F781BF",
#   "Stromal"       = "#FF7F00",
#   "Perivascular"  = "#BDB76B",
#   "Endothelial"   = "#4DAF4A",
#   "Immune"        = "#00CED1",
#   "Macrophage"    = "#1E90FF",
#   "Unknown1"      = "#00FFFF",
#   "Unknown2"      = "#9370DB",
#   "Bcells"        = "#DA70D6",
#   "Ciliated"      = "#CD5C5C"
# )

celltype_colors <- c(
  "Epithelial"    = "#F781BF",
  "Stromal"       = "#FF7F00",
  "Perivascular"  = "#BDB76B",
  "Endothelial"   = "#4DAF4A",
  "Unknown1"      = "#00FFFF",
  "Unknown2"      = "#9370DB",
  "Ciliated"      = "#CD5C5C"
)

# -----------------------
# Step 2: Base UMAP plot
# -----------------------
umap_plot <- plotUMAP(sce_polyps, colour_by = "celltypes") +
  scale_color_manual(values = celltype_colors, name = "Celltypes") +
  guides(color = guide_legend(override.aes = list(size = 3)))

umap_plot
# -----------------------
# Step 3: Process MST edge data
# -----------------------
line.data_clean <- line.data %>%
  separate(edge, into = c("from", "to"), sep = "--")

edge_summary <- line.data_clean %>%
  group_by(from, to) %>%
  summarize(
    from_UMAP1 = first(UMAP_1),
    from_UMAP2 = first(UMAP_2),
    to_UMAP1 = last(UMAP_1),
    to_UMAP2 = last(UMAP_2),
    mid_x = mean(c(first(UMAP_1), last(UMAP_1))),
    mid_y = mean(c(first(UMAP_2), last(UMAP_2))),
    edge_label = paste(from, "→", to),
    .groups = "drop"
  )

# -----------------------
# Step 4: Terminal nodes
# -----------------------
terminal_nodes <- c("Perivascular", "Endothelial", "Ciliated")

terminal_edges <- edge_summary %>%
  filter(to %in% terminal_nodes | from %in% terminal_nodes) %>%
  mutate(
    terminal_is_to = to %in% terminal_nodes,
    final_from = ifelse(terminal_is_to, from, to),
    final_to = ifelse(terminal_is_to, to, from),
    final_from_x = ifelse(terminal_is_to, from_UMAP1, to_UMAP1),
    final_from_y = ifelse(terminal_is_to, from_UMAP2, to_UMAP2),
    final_to_x = ifelse(terminal_is_to, to_UMAP1, from_UMAP1),
    final_to_y = ifelse(terminal_is_to, to_UMAP2, from_UMAP2)
  ) %>%
  distinct(final_to, .keep_all = TRUE)

# Assign path labels manually
terminal_edges <- terminal_edges %>%
  mutate(path_label = case_when(
    final_to == "Perivascular" ~ "Path 1",
    final_to == "Endothelial"  ~ "Path 2",
    final_to == "Ciliated"     ~ "Path 3",
    TRUE ~ NA_character_
  ))

# -----------------------
# Step 5: Segment MST into segments
# -----------------------
line.data_segments <- line.data %>%
  separate(edge, into = c("from", "to"), sep = "--", remove = FALSE) %>%
  group_by(edge) %>%
  mutate(
    x_start = lag(UMAP_1),
    y_start = lag(UMAP_2),
    x_end = UMAP_1,
    y_end = UMAP_2
  ) %>%
  filter(!is.na(x_start)) %>%
  ungroup()

# -----------------------
# Step 6: Annotate segments by path
# -----------------------
line.data_segments <- line.data_segments %>%
  mutate(path_id = case_when(
    edge %in% c("Unknown1--Stromal", "Unknown1--Perivascular") ~ "Path 1",
    edge %in% c("Unknown1--Stromal", "Unknown2--Unknown1", "Unknown2--Endothelial") ~ "Path 2",
    edge %in% c("Unknown1--Stromal", "Unknown2--Epithelial", "Ciliated--Epithelial") ~ "Path 3",
    TRUE ~ NA_character_
  ))

view(line.data_segments)

terminal_edges <- terminal_edges %>%
  mutate(path_id = path_label)

# -----------------------
# Step 7: Define path colors
# -----------------------

path_colors <- c(
  "Path 1" = "#d7191c",  # red
  "Path 2" = "#2E8B57",  # blue  # orange
  "Path 3" = "#DDA0DD"   # green
)

# -----------------------
# Step 8: Final plot construction
# ----------------------- 
# Add UMAP cell points as background
# Prepare UMAP coordinates
umap_df <- as.data.frame(reducedDim(sce_polyps, "UMAP"))
umap_df$celltypes <- colData(sce_polyps)$celltypes
umap_df
# Plot from scratch using ggplot
umap_plot <- ggplot() +
  
  # 1. UMAP cell points by celltype
  geom_point(data = umap_df, aes(x = UMAP_1, y = UMAP_2, color = celltypes), size = 1.5, alpha = 0.7) +
  scale_color_manual(name = "Celltypes", values = celltype_colors) +
  guides(
    color = guide_legend(title = "Celltypes", override.aes = list(size = 3)),
    fill = guide_legend(title = "Trajectory Paths", override.aes = list(size = 3))
  )+
  
  # 2. Trajectory segments colored by path
  geom_segment(data = line.data_segments,
               aes(x = x_start, y = y_start,
                   xend = x_end, yend = y_end, color = path_id),
               size = 1.2, lineend = "round", show.legend = TRUE) +
  
  # # 3. Add black arrowhead outline
  # geom_segment(data = terminal_edges,
  #              aes(x = final_from_x, y = final_from_y,
  #                  xend = final_to_x, yend = final_to_y),
  #              arrow = arrow(length = unit(0.3, "inches"), type = "closed"),
  #              color = "black", size = 2.2, lineend = "round") +
  # 
  # 4. Add colored arrowhead fill (overlaid)
  geom_segment(data = terminal_edges,
               aes(x = final_from_x, y = final_from_y,
                   xend = final_to_x, yend = final_to_y,
                   color = path_id),
               arrow = arrow(length = unit(0.25, "inches"), type = "closed"),
               size = 1, lineend = "round") +
  
  # 5. Add path labels at endpoints
  geom_text_repel(data = terminal_edges,
                  aes(x = final_to_x, y = final_to_y, label = path_label),
                  size = 4, fontface = "bold", color = "black",
                  box.padding = 0.4, segment.color = NA) +
  
  # 6. Apply both legends cleanly
  scale_color_manual(
    name = "Celltypes",
    values = c(celltype_colors, path_colors),
    breaks = c(names(celltype_colors), names(path_colors)),
    labels = c(names(celltype_colors), names(path_colors)),
    guide = guide_legend(override.aes = list(size = 3))
  ) +
  
  # 7. Theme tweaks
  theme_minimal() +
  labs(title = "adEN") + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+

  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(color = guide_legend(override.aes = list(size = 3)))


umap_plot
dev.off()

ggsave(filename = "Healthy_TSCAN_with_path_info_noimmune.png", 
       plot = umap_plot, 
       width = 10,   # Increased width
       height = 8,   # Increased height
       dpi = 600,
       bg = "white")


umap_plot

dev.off()

#################################################################################



mst_healthy <- "mst_healthy_unknown_1_2.jpeg"
# Open a JPEG device with 600 DPI resolution
jpeg(mst_healthy, width = 7, height = 7, units = "in", res = 600)


plotUMAP(sce_polyps, colour_by="celltypes") + 
  geom_line(data=line.data, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))

dev.off()

colLabels(sce_polyps) <- colData(sce_polyps)$celltypes
colLabels(sce_polyps)

save(sce_polyps, file = "sce_healthy_u1_2_noimmune.RData")
load("sce_healthy_u1_2.RData")

map.tscan <- mapCellsToEdges(sce_polyps, mst=mst, use.dimred="PCA")
map.tscan
tscan.pseudo <- orderCells(map.tscan, mst)
tscan.pseudo
#save(tscan.pseudo, file = "tscan_pseudo_polyps.RData")
save(tscan.pseudo, file = "tscan_pseudo_healthy_U1_2.RData")
load("tscan_pseudo_healthy_U1_2.RData")

colnames(tscan.pseudo)
tscan.pseudo@metadata
pathnames(tscan.pseudo)
colData(tscan.pseudo[,2])


common.pseudo <- averagePseudotime(tscan.pseudo) 
common.pseudo

pseudo_healthy <- "pseudotime_healthy_unknown_1_2_no_labels.jpeg"
# Open a JPEG device with 600 DPI resolution
jpeg(pseudo_healthy, width = 12, height = 10, units = "in", res = 300)

plotUMAP(sce_polyps, colour_by=I(common.pseudo)) + 
         #text_by="celltypes", text_colour="red") +
  geom_line(data=line.data, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))

dev.off()

pseudo.all <- quickPseudotime(sce_polyps, use.dimred="UMAP")
head(pseudo.all$ordering)
pseudo.all$ordering

# pseudo.og <- quickPseudotime(sce_polyps, use.dimred="UMAP", outgroup=TRUE)
# set.seed(10101)
# plot(pseudo.og$mst)
# 
# 
# pseudo.mnn <- quickPseudotime(sce_polyps, use.dimred="UMAP", with.mnn=TRUE)
# #save(pseudo.mnn, file = "pseudo_mnn_polyps.RData")
# save(pseudo.mnn, file = "pseudo_mnn_healthy.RData")
# 
# 
# mnn.pseudo <- averagePseudotime(pseudo.mnn$ordering)
# plotUMAP(sce_polyps, colour_by=I(mnn.pseudo), text_by="celltypes", text_colour="red") +
#   geom_line(data=pseudo.mnn$connected$UMAP, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))
# 
# 
# library(slingshot)
# sce.sling <- slingshot(sce_polyps, reducedDim='UMAP', clusterLabels = colData(sce_polyps)$celltypes)
# head(sce.sling$slingPseudotime_1)
# sce.sling$slingPseudotime_2
# 
# #save(sce.sling, file = "sce_sling_polyps.RData")
# save(sce.sling, file = "sce_sling_healthy.RData")
# 
# embedded <- embedCurves(sce.sling, "UMAP")
# embedded <- slingCurves(embedded)[[1]] # only 1 path.
# embedded <- data.frame(embedded$s[embedded$ord,])
# 
# plotUMAP(sce.sling, colour_by="slingPseudotime_1") +
#   geom_path(data=embedded, aes(x=UMAP_1, y=UMAP_2), linewidth =1.2)
# 
# sce.sling2 <- slingshot(sce_polyps, cluster=colLabels(sce_polyps), reducedDim='UMAP')
# pseudo.paths <- slingPseudotime(sce.sling2)
# head(pseudo.paths)
# 
# #save(sce.sling2, file = "sce_sling2_polyps.RData")
# save(sce.sling2, file = "sce_sling2_healthy.RData")
# 
# # Taking the rowMeans just gives us a single pseudo-time for all cells. Cells
# # in segments that are shared across paths have similar pseudo-time values in 
# # all paths anyway, so taking the rowMeans is not particularly controversial.
# 
# shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)
# 
# # Need to loop over the paths and add each one separately.
# gg <- plotUMAP(sce.sling2, colour_by=I(shared.pseudo))
# embedded <- embedCurves(sce.sling2, "UMAP")
# embedded <- slingCurves(embedded)
# 
# for (path in embedded) {
#   embedded <- data.frame(path$s[path$ord,])
#   gg <- gg + geom_path(data=embedded, aes(x=UMAP_1, y=UMAP_2), linewidth = 1.2)
# }
# 
# gg
# 
# 
# curve.assignments <- slingBranchID(sce.sling2)
# table(curve.assignments)
# 
# 
# sce.sling3 <- slingshot(sce_polyps, cluster=colLabels(sce_polyps), 
#                         reducedDim='UMAP', approx_points=100)
# 
# pseudo.paths3 <- slingPseudotime(sce.sling3)
# head(pseudo.paths3)
# 
# sce.sling4 <- slingshot(sce_polyps, cluster=colLabels(sce_polyps), 
#                         reducedDim='UMAP', approx_points=100, omega=TRUE)
# 
# pseudo.paths4 <- slingPseudotime(sce.sling4)
# head(pseudo.paths4)
# head(tscan.pseudo)
# 
# shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)
# 
# gg <- plotUMAP(sce.sling4, colour_by=I(shared.pseudo))
# embedded <- embedCurves(sce.sling4, "UMAP")
# embedded <- slingCurves(embedded)
# for (path in embedded) {
#   embedded <- data.frame(path$s[path$ord,])
#   gg <- gg + geom_path(data=embedded, aes(x=UMAP_1, y=UMAP_2), linewidth =1.2)
# }
# gg


##############################  Characterizing trajectories #####################################################
                               #Changes along a trajectory
#################################################################################################################

tscan.pseudo
pathnames(tscan.pseudo)
tscan.pseudo[,2]

# pseudo_perivascular <- testPseudotime(sce_polyps, pseudotime=tscan.pseudo[,2])[[1]]
# pseudo_ciliated <- testPseudotime(sce_polyps, pseudotime=tscan.pseudo[,1])[[1]]
# pseudo_endothelial <- testPseudotime(sce_polyps, pseudotime=tscan.pseudo[,3])[[1]]
# pseudo_macrophage <- testPseudotime(sce_polyps, pseudotime=tscan.pseudo[,4])[[1]]

pseudo_perivascular <- testPseudotime(sce_polyps, pseudotime=tscan.pseudo[,1])[[1]]
pseudo_endothelial <- testPseudotime(sce_polyps, pseudotime=tscan.pseudo[,2])[[1]]
#pseudo_macrophage <- testPseudotime(sce_polyps, pseudotime=tscan.pseudo[,3])[[1]]
pseudo_ciliated <- testPseudotime(sce_polyps, pseudotime=tscan.pseudo[,3])[[1]]

pseudo_perivascular[order(pseudo_perivascular$p.value),]
table(pseudo_perivascular$p.value < 0.05)

pseudo_ciliated[order(pseudo_ciliated$p.value),]
table(pseudo_ciliated$p.value < 0.05)

pseudo_endothelial[order(pseudo_endothelial$p.value),]
table(pseudo_endothelial$p.value < 0.05)


pseudo_macrophage[order(pseudo_macrophage$p.value),]
table(pseudo_macrophage$p.value < 0.05)


pathStat(tscan.pseudo)[,1]
pathStat(tscan.pseudo)[,2]
pathStat(tscan.pseudo)[,3]
#pathStat(tscan.pseudo)[,4]

#load("sce_healthy_u1_2.RData")
load("sce_healthy_u1_2_noimmune.RData")
load("tscan_pseudo_healthy_U1_2.RData")


sce_polyps$TSCAN.first <- pathStat(tscan.pseudo)[,1]
sce_polyps$TSCAN.second <- pathStat(tscan.pseudo)[,2]
sce_polyps$TSCAN.third <- pathStat(tscan.pseudo)[,3]
#sce_polyps$TSCAN.fourth <- pathStat(tscan.pseudo)[,4]

colData(sce_polyps)
str(sce_polyps)
str(tscan.pseudo)
# Testing against the first path again.
pseudo2 <- testPseudotime(sce_polyps, pseudotime=sce_polyps$TSCAN.first)
pseudo2
#save(pseudo2, file = "pseudo_ciliated_polyps.RData")
save(pseudo2, file = "pseudo_perivascular_healthy_U1_2.RData")
load("pseudo_perivascular_healthy_U1_2.RData")

pseudo3 <- testPseudotime(sce_polyps, pseudotime=sce_polyps$TSCAN.second)
pseudo3
#save(pseudo3, file = "pseudo_perivascular_polyps.RData")
save(pseudo3, file = "pseudo_endothelial_healthy_U1_2.RData")
load("pseudo_endothelial_healthy_U1_2.RData")

# pseudo4 <- testPseudotime(sce_polyps, pseudotime=sce_polyps$TSCAN.third)
# pseudo4
# #save(pseudo4, file = "pseudo_endothelial_polyps.RData")
# save(pseudo4, file = "pseudo_macrophage_healthy_U1_2.RData")
# load("pseudo_macrophage_healthy_U1_2.RData")


pseudo5 <- testPseudotime(sce_polyps, pseudotime=sce_polyps$TSCAN.third)
pseudo5
#save(pseudo5, file = "pseudo_macrophage_polyps.RData")
save(pseudo5, file = "pseudo_ciliated_healthy_U1_2.RData")
load("pseudo_ciliated_healthy_U1_2.RData")


pseudo2
pseudo_ciliated

# write.csv(pseudo_ciliated, file = "pseudo_ciliated_polyps.csv")
# write.csv(pseudo_endothelial, file = "pseudo_endothelial_polyps.csv")
# write.csv(pseudo_macrophage, file = "pseudo_macrophage_polyps.csv")
# write.csv(pseudo_perivascular, file = "pseudo_perivascular_polyps.csv")

write.csv(pseudo_ciliated, file = "pseudo_ciliated_healthy_U1_2.csv")
write.csv(pseudo_endothelial, file = "pseudo_endothelial_healthy_U1_2.csv")
#write.csv(pseudo_macrophage, file = "pseudo_macrophage_healthy_U1_2.csv")
write.csv(pseudo_perivascular, file = "pseudo_perivascular_healthy_U1_2.csv")


###########################################################################################################################
                                             ## Heat map 
############################################################################################################################

genelist <- readxl::read_excel(path = "Trajectory marker list_27.2.25.xlsx", 
                               sheet = "EMT-MET genes_mesenchymal")


genelist
#colnames(genelist) <- "genes"
genelist <- genelist$`EMT-MET genes`
genelist

#genelist <- gsub("PDGFRβ","PDGFRB", genelist)

#pseudo5, pseudo4
table(is.na(match(rownames(pseudo3), genelist)))


## custom genelist,
genelist <- c("MECOM", "EYA2", "MAGI1", "MYH11",
              "ADIRF", "TAGLN", "MYL9", "TPM2")

#PATJ, EYA2, KIAA1324, PAX8
genelist2 <- c("MECOM", "RHEX", "CTNNA2", "LDB2",
               "VWF", "A2M")

#factor(rownames(genelist_sorted), levels =levels(factor(rownames(genelist_sorted))))
genelist_sorted <- pseudo3[!is.na(match(rownames(pseudo3), genelist2)),]
genelist_sorted
on.path <- !is.na(sce_polyps$TSCAN.second)
#length(on.path)

valid_genes <- intersect(rownames(sce_polyps), rownames(genelist_sorted))
#valid_genes
length(valid_genes)

plot <- "heatmap_path2_endothelial_healthy.jpeg"
# Open a JPEG device with 600 DPI resolution
jpeg(plot, width = 7, height = 7, 
     units = "in", res = 300, pointsize = 15)


plotHeatmap(sce_polyps[,on.path], order_columns_by="TSCAN.second", 
            colour_columns_by="celltypes", features= valid_genes,
            center=TRUE)

dev.off()

###############################################################################
sorted <- pseudo5[order(pseudo5$p.value),]
sorted


up.left <- sorted[sorted$logFC < 0,]
head(up.left, 10)


best <- head(rownames(up.left), 10)
best

plotExpression(sce_polyps, features=best,
               x="TSCAN.fourth", colour_by="celltypes")

dev.off()

up.right <- sorted[sorted$logFC > 0,]
head(up.right, 10)

best <- head(rownames(up.right), 10)
best

plotExpression(sce_polyps, features=best,
               x="TSCAN.fourth", colour_by="celltypes")

on.first.path <- !is.na(sce_polyps$TSCAN.third)
on.first.path
dev.off()
plotHeatmap(sce_polyps[,on.first.path], order_columns_by="TSCAN.third", 
            colour_columns_by="celltypes", features=head(rownames(up.right), 50),
            center=TRUE)

plotHeatmap(sce_polyps[,on.first.path], order_columns_by="TSCAN.third", 
            colour_columns_by="celltypes", features=head(rownames(up.left), 50),
            center=TRUE)


