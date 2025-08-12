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
# Unknown 	       10, 12


# new.cluster.ids <- c("0 Epithelial",
#                      "1 Stromal",
#                      "2 Perivascular",
#                      "3 Endothelial",
#                      "4 Immune",
#                      "5 Epithelial",
#                      "6 Immune",
#                      "7 Stromal",
#                      "8 Immune",
#                      "9 Epithelial",
#                      "10 Unknown1",
#                      "11 Macrophage",
#                      "12 Unknown2",
#                      "13 Bcells",
#                      "14 Ciliated",
#                      "15 Immune",
#                      "16 Immune",
#                      "17 Immune",
#                      "18 Immune",
#                      "19 Perivascular",
#                      "20 Immune")

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

seurat_integrated@meta.data$celltypes <- factor(seurat_integrated@active.ident)
seurat_integrated@meta.data$celltypes
seurat_integrated$integrated_snn_res.0.4

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

# Assuming 'seurat_obj' is your Seurat object
seurat_integrated$sampleType
seurat_integrated$sampleType <- as.factor(seurat_integrated$sampleType)  # Ensure 'sampletype' is a factor

# Rename "Healthy" to "Endometrium"
levels(seurat_integrated$sampleType)[levels(seurat_integrated$sampleType) == "Healthy"] <- "adEN"
levels(seurat_integrated$sampleType)[levels(seurat_integrated$sampleType) == "Polyps"] <- "EPs"

levels(seurat_integrated$sampleType)


# Now, create the UMAP plot again
p1 <- DimPlot(seurat_integrated,
              reduction = "umap",
              label = FALSE, 
              label.size = 4, 
              pt.size = 0.8, 
              split.by = "sampleType") 

ggsave(filename = "withoutlabels_umap_clusterlabels.png", 
       plot = p1, 
       width = 12,   # Increased width
       height = 8,   # Increased height
       dpi = 600)

dev.off()




#######################################################

cluster_labels <- c(
  "0" = "0: Epithelial", "1" = "1: Stromal", "2" = "2: Perivascular",
  "3" = "3: Endothelial", "4" = "4: Immune", "5" = "5: Epithelial",
  "6" = "6: Immune", "7" = "7: Stromal", "8" = "8: Immune",
  "9" = "9: Epithelial", "10" = "10: Unknown1", "11" = "11: Macrophage",
  "12" = "12: Unknown2", "13" = "13: B cells", "14" = "14: Ciliated",
  "15" = "15: Immune", "16" = "16: Immune", "17" = "17: Immune",
  "18" = "18: Immune", "19" = "19: Perivascular", "20" = "20: Immune"
)

Idents(seurat_integrated) <- "integrated_snn_res.0.4"  # or whatever the original cluster identities are

# Get unique cluster levels
cluster_levels <- levels(Idents(seurat_integrated))

# Assign labels for scale_color_manual
new_labels <- cluster_labels[cluster_levels]

# Assign custom colors if desired (optional)
custom_colors <- scales::hue_pal()(length(cluster_levels))
names(custom_colors) <- cluster_levels

dev.off()

p <- DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4,
        pt.size = 0.8, #+
        split.by = "sampleType") +
  scale_color_manual(
    values = custom_colors,
    labels = new_labels
  ) +
  labs(color = "Cluster: Cell Type") +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        panel.grid = element_blank(),
        panel.background = element_blank())

# Save the plot with larger size and 600 dpi  umap_integrated_custom_legeds.png, umap_integrated_largepoints.png
ggsave(filename = "sepratedby_sampletype.png", 
       plot = p, 
       width = 12,   # Increased width
       height = 8,   # Increased height
       dpi = 600)

###################################################



# Generate the UMAP plot with increased point size
umap_plot <- DimPlot(seurat_integrated,
                     reduction = "umap",
                     label = TRUE, 
                     label.size = 4, 
                     pt.size = 0.8, 
                     split.by = "sampleType",
                     group.by = "cluster_label") + 
  labs(color = "Cluster (Cell Type)") + 
  theme_minimal()

umap_plot

# Save the UMAP plot with 600 dpi resolution
ggsave(filename = "umap_plot_new.png", 
       plot = umap_plot, 
       width = 8, 
       height = 6, 
       dpi = 600)

library(ggplot2)

# Generate the UMAP plot
umap_plot <- DimPlot(seurat_integrated,
                     reduction = "umap",
                     label = TRUE, 
                     label.size = 4, 
                     pt.size = 0.8) +
                     #split.by = "sampleType") +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  xlim(-15, 15) +  # Adjust according to your data spread
  ylim(-15, 15)

# Save the plot with larger size and 600 dpi
ggsave(filename = "umap_integrated.png", 
       plot = umap_plot, 
       width = 10,   # Increased width
       height = 8,   # Increased height
       dpi = 600)


library(grid)
library(ggplot2)

# Function to extract legend
extract_legend <- function(ggplot_obj) {
  ggplot_grob <- ggplotGrob(ggplot_obj)
  legend <- ggplot_grob$grobs[which(sapply(ggplot_grob$grobs, function(x) x$name) == "guide-box")]
  return(legend)
}

# Extract the legend
legend_plot <- extract_legend(umap_plot)

# Save the legend separately as PNG
png("legend_plot_corner_axis.png", width = 400, height = 800, res = 600)
grid.draw(legend_plot)
dev.off()
# Save the legend separately as PNG with 600 dpi
png("legend_plot_corner_axis.png", width = 400, height = 800, res = 600)
grid.newpage()
grid.draw(legend_plot)
dev.off()



##### trajectory analysis for Polyps samples

# Create a Seurat object for condition "B"
seurat_polyps <- subset(seurat_integrated, subset = sampleType == "EPs")

DimPlot(seurat_polyps,
        reduction = "umap",
        label = TRUE, 
        label.size = 4)
dev.off()


# Remove Immune, Bcells, and Macrophage clusters
seurat_polyps_filtered <- subset(seurat_polyps, idents = c("Immune", "Bcells", "Macrophage"), invert = TRUE)

rm(seurat_polyps)

#seurat_healthy@meta.data$celltypes <- seurat_healthy@active.ident
#seurat_polyps@meta.data$celltypes <- seurat_polyps@active.ident

seurat_polyps_filtered@meta.data$celltypes
# Set up metadata as desired for aggregation and DE analysis
#seurat_polyps@meta.data$celltypes <- factor(seurat_polyps@active.ident)


# Create single cell experiment object
sce_polyps <- as.SingleCellExperiment(seurat_polyps_filtered)
sce_polyps$celltypes

seurat_polyps_filtered[['umap']]@cell.embeddings
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

#################################################################################

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
terminal_nodes <- c("Perivascular","Ciliated", "Endothelial")

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
    final_to == "Ciliated" ~ "Path 2",
    final_to == "Endothelial"  ~ "Path 3",
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

line.data_segments

# -----------------------
# Step 6: Annotate segments by path
# -----------------------
line.data_segments <- line.data_segments %>%
  mutate(path_id = case_when(
    edge %in% c("Unknown1--Stromal", "Unknown1--Perivascular", "Unknown1--Epithelial") ~ "Path 1",
    edge %in% c("Unknown1--Stromal", "Unknown1--Epithelial", "Ciliated--Epithelial") ~ "Path 2",
    edge %in% c("Unknown1--Stromal", "Unknown2--Epithelial", "Unknown2--Endothelial") ~ "Path 3",
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
  "Path 2" = "#DDA0DD",  # blue
  "Path 3" = "#2E8B57"  # orange
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
  
  # 3. Add black arrowhead outline
  # geom_segment(data = terminal_edges,
  #              aes(x = final_from_x, y = final_from_y,
  #                  xend = final_to_x, yend = final_to_y),
  #              arrow = arrow(length = unit(0.3, "inches"), type = "closed"),
  #              color = "black", size = 0.5, lineend = "round") +

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
  labs(title = "EPs") + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(color = guide_legend(override.aes = list(size = 3)))

umap_plot
dev.off()

ggsave(filename = "Polyps_TSCAN_with_path_info_noimmune.png", 
       plot = umap_plot, 
       width = 10,   # Increased width
       height = 8,   # Increased height
       dpi = 600,
       bg = "white")


#############################################################################



# Plot
# Construct the filename for saving the heatmap
mst_polyps <- "MST_polyps_unknown_1_2_noimmune.jpeg"
# Open a JPEG device with 600 DPI resolution
jpeg(mst_polyps, width = 7, height = 7, units = "in", res = 600)


plotUMAP(sce_polyps, colour_by="celltypes") + 
  geom_line(data=line.data, 
            mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))

dev.off()


####################################################
library(dplyr)
library(igraph)


# Step 1: Extract edges from the MST
mst_edges <- as_data_frame(mst, what = "edges")

# Rename columns for clarity
colnames(mst_edges) <- c("from", "to", "weight", "gain")

# Display the edges (Optional)
head(mst_edges)

rm(mst_edges_expanded)

# Step 1: Join UMAP coordinates for 'from'
mst_edges_expanded <- mst_edges %>%
  left_join(umap_df, by = c("from" = "celltypes")) %>%
  mutate(UMAP_1.x = UMAP_1, UMAP_2.x = UMAP_2) %>%
  select(-UMAP_1, -UMAP_2)  # Remove the original columns

# Step 2: Join UMAP coordinates for 'to'
mst_edges_expanded <- mst_edges_expanded %>%
  left_join(umap_df, by = c("to" = "celltypes")) %>%
  mutate(UMAP_1.y = UMAP_1, UMAP_2.y = UMAP_2) %>%
  select(-UMAP_1, -UMAP_2)  # Remove the original columns

  
# Display the result
head(mst_edges_expanded)

library(ggplot2)
library(scater)

# Plot UMAP with cell types
p <- plotUMAP(sce_polyps, colour_by = "celltypes") +
  
  # Plot MST lines with arrows
  geom_segment(data = mst_edges_expanded, 
               aes(x = UMAP_1.x, y = UMAP_2.x, 
                   xend = UMAP_1.y, yend = UMAP_2.y,
                   group = interaction(from, to)), 
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", 
               size = 0.6) +
  
  # Optional: Add labels for start and end points
  geom_text(data = mst_edges_expanded, 
            aes(x = UMAP_1.x, y = UMAP_2.x, label = from),
            color = "blue", size = 3, vjust = -1) +
  geom_text(data = mst_edges_expanded, 
            aes(x = UMAP_1.y, y = UMAP_2.y, label = to),
            color = "red", size = 3, vjust = -1) +
  
  theme_minimal() +
  theme(legend.position = "right")

ggsave("UMAP_MST_Trajectory.png", plot = p, width = 10, height = 8, dpi = 300)


#####################################################

colLabels(sce_polyps) <- colData(sce_polyps)$celltypes
colLabels(sce_polyps)

save(sce_polyps, file = "sce_polyps_u1_2_noimmune.RData")
load("sce_polyps_u1_2.RData")

map.tscan <- mapCellsToEdges(sce_polyps, mst=mst, use.dimred="PCA")
map.tscan
tscan.pseudo <- orderCells(map.tscan, mst)
tscan.pseudo
save(tscan.pseudo, file = "tscan_pseudo_polyps_unknown1_2_noimmune.RData")
load("tscan_pseudo_polyps_unknown1_2.RData")
tscan.pseudo

common.pseudo <- averagePseudotime(tscan.pseudo) 
common.pseudo

pseudo_polyps <- "pseudotime_polyps_unknown_1_2_nolabels.jpeg"
# Open a JPEG device with 600 DPI resolution
jpeg(pseudo_polyps, width = 12, height = 10, units = "in", res = 600)


plotUMAP(sce_polyps, colour_by=I(common.pseudo)) + 
         #text_by="celltypes", text_colour="red") +
  geom_line(data=line.data, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))

dev.off()


pseudo.all <- quickPseudotime(sce_polyps, use.dimred="UMAP")

#############################################################################################################
                                  #Characterizing trajectories
                                  #Changes along a trajectory
##############################################################################################################

tscan.pseudo
pathnames(tscan.pseudo)
tscan.pseudo[,2]

pseudo_perivascular <- testPseudotime(sce_polyps, pseudotime=tscan.pseudo[,1])[[1]]
pseudo_ciliated <- testPseudotime(sce_polyps, pseudotime=tscan.pseudo[,2])[[1]]
pseudo_endothelial <- testPseudotime(sce_polyps, pseudotime=tscan.pseudo[,3])[[1]]
#pseudo_macrophage <- testPseudotime(sce_polyps, pseudotime=tscan.pseudo[,4])[[1]]

pseudo_perivascular[order(pseudo_perivascular$p.value),]
table(pseudo_perivascular$p.value < 0.05)

pseudo_ciliated[order(pseudo_ciliated$p.value),]
table(pseudo_ciliated$p.value < 0.05)

pseudo_endothelial[order(pseudo_endothelial$p.value),]
table(pseudo_endothelial$p.value < 0.05)


# pseudo_macrophage[order(pseudo_macrophage$p.value),]
# table(pseudo_macrophage$p.value < 0.05)

pathStat(tscan.pseudo)[,1]
pathStat(tscan.pseudo)[,2]
pathStat(tscan.pseudo)[,3]
#pathStat(tscan.pseudo)[,4]



load("sce_polyps_u1_2_noimmune.RData")
#load("tscan_pseudo_polyps_unknown1_2.RData")
load("tscan_pseudo_polyps_unknown1_2_noimmune.RData")



sce_polyps$TSCAN.first <- pathStat(tscan.pseudo)[,1]
sce_polyps$TSCAN.second <- pathStat(tscan.pseudo)[,2]
sce_polyps$TSCAN.third <- pathStat(tscan.pseudo)[,3]
#sce_polyps$TSCAN.fourth <- pathStat(tscan.pseudo)[,4]


# Testing against the first path again.
pseudo2 <- testPseudotime(sce_polyps, pseudotime=sce_polyps$TSCAN.first)
pseudo2
save(pseudo2, file = "pseudo_perivascular_polyps_U1_2.RData")
load("pseudo_perivascular_polyps_U1_2.RData")

pseudo3 <- testPseudotime(sce_polyps, pseudotime=sce_polyps$TSCAN.second)
pseudo3
save(pseudo3, file = "pseudo_ciliated_polyps_U1_2.RData")
load("pseudo_ciliated_polyps_U1_2.RData")

pseudo4 <- testPseudotime(sce_polyps, pseudotime=sce_polyps$TSCAN.third)
pseudo4
save(pseudo4, file = "pseudo_endothelial_polyps_U1_2.RData")
load("pseudo_endothelial_polyps_U1_2.RData")

# pseudo5 <- testPseudotime(sce_polyps, pseudotime=sce_polyps$TSCAN.fourth)
# pseudo5
# save(pseudo5, file = "pseudo_macrophage_polyps_U1_2.RData")
# load("pseudo_macrophage_polyps_U1_2.RData")

pseudo2
pseudo_ciliated

write.csv(pseudo_ciliated, file = "pseudo_ciliated_polyps_U1_2.csv")
write.csv(pseudo_endothelial, file = "pseudo_endothelial_polyps_U1_2.csv")
#write.csv(pseudo_macrophage, file = "pseudo_macrophage_polyps_U1_2.csv")
write.csv(pseudo_perivascular, file = "pseudo_perivascular_polyps_U1_2.csv")

#pseudo5

#######################################################################################################################
                                                      #Heat map
########################################################################################################################


# read the genelist file to make plots,
genelist <- readxl::read_excel(path = "Trajectory marker list_27.2.25.xlsx", 
                               sheet = "intermediate EMT")


genelist
genelist <- genelist$`MET-EMT intermediate state`
genelist

#genelist <- gsub("PDGFRβ","PDGFRB", genelist)

#pseudo5, pseudo4
table(is.na(match(rownames(pseudo2), genelist)))


# genelist_order <-  sort(rownames(genelist_sorted))
# genelist_order
# 
# genelist_sorted <- genelist_sorted[genelist_order,]
# genelist_sorted


## custom genelist
genelist <- c("MECOM", "EYA2", "MAGI1", "MYH11",
              "ADIRF", "TAGLN", "MYL9", "TPM2")

#PATJ, EYA2, KIAA1324, PAX8
genelist2 <- c("MECOM", "RHEX", "CTNNA2", "LDB2",
               "VWF", "A2M")

pseudo2
#factor(rownames(genelist_sorted), levels =levels(factor(rownames(genelist_sorted))))
genelist_sorted <- pseudo4[!is.na(match(rownames(pseudo4), genelist2)),]
genelist_sorted
on.path <- !is.na(sce_polyps$TSCAN.third)
on.path
#length(on.path)

valid_genes <- intersect(rownames(sce_polyps), rownames(genelist_sorted))
#valid_genes
length(valid_genes)

valid_genes

plot <- "heatmap_path3_endothelial_polyps.jpeg"
# Open a JPEG device with 600 DPI resolution
jpeg(plot, width = 7, height = 7, 
     units = "in", res = 300, pointsize = 15)


plotHeatmap(sce_polyps[,on.path], order_columns_by="TSCAN.third", 
            colour_columns_by="celltypes", features= valid_genes,
            center=TRUE)

dev.off()


####################################################################################
sorted <- pseudo4[order(pseudo4$p.value),]
sorted


up.left <- sorted[sorted$logFC < 0,]
head(up.left, 10)
head(rownames(up.right), 50)

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
plotHeatmap(sce_polyps[,on.first.path], order_columns_by="TSCAN.third", 
            colour_columns_by="celltypes", features=head(rownames(up.right), 50),
            center=TRUE)

plotHeatmap(sce_polyps[,on.first.path], order_columns_by="TSCAN.third", 
            colour_columns_by="celltypes", features=head(rownames(up.left), 50),
            center=TRUE)


