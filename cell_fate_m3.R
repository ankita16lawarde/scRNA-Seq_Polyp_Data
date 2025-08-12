
##############################################################################
seurat_healthy_filtered
seurat_polyps_filtered

cds_polyps <- as.cell_data_set(seurat_polyps_filtered)

#seurat_healthy_filtered@reductions$pca@cell.embeddings

# Extract PCA embeddings from Seurat
#pca_coords <- seurat_polyps_filtered@reductions$pca@cell.embeddings

pca_coords <- seurat_polyps_filtered@reductions$pca@cell.embeddings

# Ensure cell order matches
pca_coords <- pca_coords[colnames(cds_polyps), ]

# Assign PCA to Monocle3 reducedDims slot
reducedDims(cds_polyps)$PCA <- pca_coords



cds_polyps <- reduce_dimension(cds_polyps, reduction_method = "UMAP", 
                               preprocess_method = "PCA")



cds_polyps <- cluster_cells(cds_polyps, reduction_method = "UMAP")

cds_polyps <- learn_graph(cds_polyps, use_partition = F)
cds_polyps$celltypes

root_cell <- colnames(cds_polyps)[cds_polyps$celltypes == "Stromal"][1]
root_cell
cds_polyps <- order_cells(cds_polyps, root_cells = root_cell)


plot_cells(cds_polyps,
           color_cells_by = "celltypes",         # Color by cell type
           label_groups_by_cluster = FALSE,      # No cluster labels
           label_leaves = FALSE,                 # No terminal node labels
           label_branch_points = TRUE,           # Show numeric branch labels
           graph_label_size = 6,                 # Label font size
           show_trajectory_graph = TRUE          # Show black graph lines
) +
  xlab("Component 1") +
  ylab("Component 2")

#cell_graph <- "m3_graph_PCA_polyps.jpeg"

cell_graph <- "m3_graph_PCA_healthy.jpeg"

# Open a JPEG device with 600 DPI resolution
jpeg(cell_graph, width = 9, height = 7, units = "in", res = 600)


plot_cells(cds_polyps,
           color_cells_by = "celltypes",
           label_groups_by_cluster = F,
           label_leaves = F,
           label_branch_points = F,
           graph_label_size = 6,
           group_label_size = 6,
           show_trajectory_graph = TRUE
) +
  scale_color_manual(values = c(
    "Stromal" = "#E41A1C",
    "Epithelial" = "#4DAF4A",
    "Endothelial" = "#377EB8",
    "Ciliated" = "#984EA3",
    "Perivascular" = "#FF7F00",
    "Unknown1" = "#999999",
    "Unknown2" = "#A65628"
  )) +
  xlab("Component 1") +
  ylab("Component 2")

dev.off()


#cell_pseudo <- "m3_psudotime_polyps.jpeg"

cell_pseudo <- "m3_psudotime_healthy.jpeg"

# Open a JPEG device with 600 DPI resolution
jpeg(cell_pseudo, width = 9, height = 7, units = "in", res = 600)


plot_cells(cds_polyps,
           color_cells_by = "pseudotime",          # Now coloring by pseudotime
           label_groups_by_cluster = FALSE,
           label_leaves = F,
           label_branch_points = F,
           graph_label_size = 6,
           show_trajectory_graph = TRUE,
           label_roots = T,
          # label_principal_points = T
) +
  xlab("Component 1") +
  ylab("Component 2") +
  scale_color_viridis_c(option = "C")   # Optional: Use a perceptually uniform colormap

dev.off()

# Step 1: Get principal graph from Monocle3
graph_obj <- principal_graph(cds_polyps)$UMAP

# Step 2: Get graph cluster memberships using igraph
graph_clusters <- igraph::clusters(graph_obj)

# Step 3: Map cells to closest graph node
closest_vertex <- cds_polyps@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
head(closest_vertex)
# Just extract the first (and only) column as a numeric vector
closest_vertex_vector <- closest_vertex[, 1]
closest_vertex_vector

# Ensure vertex names match (like "Y_1", "Y_2", ...)
vertex_names <- names(membership_vector)

# Match closest vertex IDs to vertex names
closest_vertex_names <- vertex_names[closest_vertex_vector]

# Assign membership (cell fate state)
cds_polyps$cell_fate_state <- membership_vector[closest_vertex_names]

unique(cds_polyps$cell_fate_state)

# Optional: label states as State 1, State 2, etc.
cds_polyps$cell_fate_state <- factor(cds_polyps$cell_fate_state,
                                     labels = paste0("State ", sort(unique(cds_polyps$cell_fate_state))))


plot_cells(cds_polyps,
           color_cells_by = "cell_fate_state",
           label_groups_by_cluster = FALSE,
           label_branch_points = TRUE,
           label_leaves = TRUE,
           graph_label_size = 6,
           show_trajectory_graph = TRUE
) +
  xlab("Component 1") +
  ylab("Component 2") +
  scale_color_manual(values = c("State 1" = "#1f78b4",
                                "State 2" = "#33a02c",
                                "State 3" = "#e31a1c"))  # Add more colors if needed

# Load igraph (if not already loaded)
library(igraph)

# 1. Get principal graph from Monocle3
g <- principal_graph(cds_polyps)$UMAP

# 2. Identify leaf nodes as nodes with degree == 1
leaf_node_indices <- which(degree(g) == 1)

# 3. Get the names of leaf nodes
leaf_node_names <- names(leaf_node_indices)

# Preview
leaf_node_names

# 1. Extract closest vertex for each cell
closest_vertex <- cds_polyps@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex[, 1]

# 2. Match closest vertex to leaf node names
cell_terminal_state <- ifelse(paste0("Y_", closest_vertex) %in% leaf_node_names,
                              paste0("State_", match(paste0("Y_", closest_vertex), leaf_node_names)),
                              NA)

# 3. Assign to metadata
cds_polyps$terminal_state <- factor(cell_terminal_state)
cds_polyps$terminal_state

plot_cells(cds_polyps,
           color_cells_by = "terminal_state",
           label_groups_by_cluster = FALSE,
           label_branch_points = TRUE,
           label_leaves = TRUE,
           graph_label_size = 2,
           show_trajectory_graph = TRUE
) +
  scale_color_manual(values = c(
    "State_1" = "#1f78b4",  # blue
    "State_2" = "#33a02c",  # green
    "State_3" = "#e31a1c",  # red
    "State_4" = "#ff7f00",  # orange
    "State_5" = "#6a3d9a",  # purple
    "State_6" = "#b15928",  # brown
    "State_7" = "#a6cee3",  # light blue
    "State_8" = "#b2df8a"   # light green
  )) +
  xlab("Component 1") +
  ylab("Component 2")


# Manually assign states to specific cell types
cds_polyps$cell_fate_state <- NA

cds_polyps$cell_fate_state[cds_polyps$celltypes == "Perivascular"] <- "State_1"
cds_polyps$cell_fate_state[cds_polyps$celltypes == "Ciliated"]     <- "State_2"
cds_polyps$cell_fate_state[cds_polyps$celltypes == "Endothelial"]  <- "State_3"

# Convert to factor
cds_polyps$cell_fate_state <- factor(cds_polyps$cell_fate_state,
                                     levels = c("State_1", "State_2", "State_3"))


plot_cells(cds_polyps,
           color_cells_by = "cell_fate_state",
           label_groups_by_cluster = FALSE,
           label_branch_points = TRUE,
           label_leaves = TRUE,
           graph_label_size = 6,
           show_trajectory_graph = TRUE
) +
  xlab("Component 1") +
  ylab("Component 2") +
  scale_color_manual(values = c(
    "State_1" = "#1f78b4",  # Perivascular
    "State_2" = "#33a02c",  # Ciliated
    "State_3" = "#e31a1c"   # Endothelial
  ))

# Extract pseudotime values
pseudotime_values <- pseudotime(cds_polyps)

# Step 1: Extract pseudotime values
pt <- pseudotime(cds_polyps)

# Step 2: Create quantile cutoffs for early, mid, late pseudotime
quantile_cutoffs <- quantile(pt, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
quantile_cutoffs

# Show quantile values
quantile_cutoffs

# State 1 (early)
cat("State_1 (early):", quantile_cutoffs[1], "to", quantile_cutoffs[2], "\n")

# State 2 (mid)
cat("State_2 (mid):", quantile_cutoffs[2], "to", quantile_cutoffs[3], "\n")

# State 3 (late)
cat("State_3 (late):", quantile_cutoffs[3], "to", quantile_cutoffs[4], "\n")

## for healthy
# State_1 (early): 0 to 13.08904 
# > # State 2 (mid)
#   > cat("State_2 (mid):", quantile_cutoffs[2], "to", quantile_cutoffs[3], "\n")
# State_2 (mid): 13.08904 to 23.80901 
# > # State 3 (late)
#   > cat("State_3 (late):", quantile_cutoffs[3], "to", quantile_cutoffs[4], "\n")
# State_3 (late): 23.80901 to 34.58028 

#for polyps
# State_1 (early): 0 to 10.80244 
# > 
#   > # State 2 (mid)
#   > cat("State_2 (mid):", quantile_cutoffs[2], "to", quantile_cutoffs[3], "\n")
# State_2 (mid): 10.80244 to 25.81128 
# > 
#   > # State 3 (late)
#   > cat("State_3 (late):", quantile_cutoffs[3], "to", quantile_cutoffs[4], "\n")
# State_3 (late): 25.81128 to 51.12953 

# Step 3: Initialize vector to store fate state
cell_fate_state <- rep(NA, length(pt))

# Step 4: Assign states based on pseudotime intervals
cell_fate_state[pt >= quantile_cutoffs[1] & pt <= quantile_cutoffs[2]] <- "State_1"  # Early
cell_fate_state[pt > quantile_cutoffs[2] & pt <= quantile_cutoffs[3]]  <- "State_2"  # Mid
cell_fate_state[pt > quantile_cutoffs[3] & pt <= quantile_cutoffs[4]]  <- "State_3"  # Late

cell_fate_state

# Step 5: Add to colData and convert to factor
cds_polyps$cell_fate_state <- factor(cell_fate_state, levels = c("State_1", "State_2", "State_3"))

table(cds_polyps$cell_fate_state, cds_polyps$celltypes)

plot_cells(cds_polyps,
           color_cells_by = "cell_fate_state",
           group_cells_by = "celltypes",
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_leaves = F,
           show_trajectory_graph = TRUE,
           label_cell_groups = F
) +
  scale_color_manual(values = c(
    "State_1" = "#1f78b4",  # early
    "State_2" = "#33a02c",  # middle
    "State_3" = "#e31a1c"   # late
  )) +
  xlab("Component 1") +
  ylab("Component 2")


library(dplyr)
library(ggplot2)
library(ggrepel)
# 1. Extract UMAP coordinates and metadata
umap_coords <- as.data.frame(reducedDims(cds_polyps)$UMAP)
umap_coords$celltypes <- cds_polyps$celltypes
umap_coords$cell_fate_state <- cds_polyps$cell_fate_state

# 2. Compute cell type label positions (e.g., median position for each type)
label_coords <- umap_coords %>%
  group_by(celltypes) %>%
  summarize(UMAP_1 = median(V1), UMAP_2 = median(V2), .groups = "drop")


#cell_fates <- "m3_cell_fate_by_psudotime_polyps.jpeg"

cell_fates <- "m3_cell_fate_by_psudotime_healthy.jpeg"

# Open a JPEG device with 600 DPI resolution
jpeg(cell_fates, width = 9, height = 7, units = "in", res = 600)


# 3. Plot with adjusted label size and point size
plot_cells(cds_polyps,
           color_cells_by = "cell_fate_state",
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           show_trajectory_graph = TRUE,
           group_cells_by = "celltypes",
           cell_size = 0.5  # 🔼 Increase point size here
) +
  scale_color_manual(
    name = "Cell Fate",
    values = c(
      "State_1" = "#a6cee3", # early 
      "State_2" = "#b2df8a", # mid
      "State_3" = "#fcae91" # late
    ),
    labels = c(
      "State_1" = "1: Early",
      "State_2" = "2: Mid",
      "State_3" = "3: Late"
    )
  ) + 
  geom_text_repel(
    data = label_coords,
    aes(x = UMAP_1, y = UMAP_2, label = celltypes),
    size = 5,
    fontface = "plain",
    box.padding = 0.4,
    point.padding = 0.3,
    segment.size = 0.3
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.position = "right") +
  xlab("Component 1") +
  ylab("Component 2")


dev.off()


# Create a data frame with pseudotime and cell type information
pt_df <- data.frame(
  Pseudotime = pseudotime(cds_polyps),
  CellType = cds_polyps$celltypes
)

# Optional: Remove NA pseudotime values (e.g., unassigned cells)
pt_df <- pt_df[!is.na(pt_df$Pseudotime), ]

# Reorder CellType factor based on median pseudotime
pt_df$CellType <- factor(pt_df$CellType,
                         levels = names(sort(tapply(pt_df$Pseudotime, pt_df$CellType, median))))


library(ggplot2)


#cell_order <- "m3_cell_order_by_psudotime_polyps.jpeg"

cell_order <- "m3_cell_order_by_psudotime_healthy.jpeg"

# Open a JPEG device with 600 DPI resolution
jpeg(cell_order, width = 9, height = 7, units = "in", res = 600)

ggplot(pt_df, aes(x = CellType, y = Pseudotime, fill = CellType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = "Cell Type",
       y = "Pseudotime") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_brewer(palette = "Set2")  # Or use scale_fill_manual for custom colors

dev.off()

