## cell cluster ordering plot
setwd("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/scRNA-Seq_polyp_data/analysis")
# Load required libraries
library(TSCAN)
library(ggplot2)
library(dplyr)
library(tidyr)

#load("sce_healthy_u1_2.RData")
load("sce_healthy_u1_2_noimmune.RData")
load("tscan_pseudo_healthy_U1_2.RData")

#load("sce_polyps_u1_2.RData")
load("sce_polyps_u1_2_noimmune.RData")
#load("tscan_pseudo_polyps_unknown1_2.RData")
load("tscan_pseudo_polyps_unknown1_2_noimmune.RData")


tscan.pseudo


sce_polyps$TSCAN.first <- pathStat(tscan.pseudo)[,1]
sce_polyps$TSCAN.second <- pathStat(tscan.pseudo)[,2]
sce_polyps$TSCAN.third <- pathStat(tscan.pseudo)[,3]
#sce_polyps$TSCAN.fourth <- pathStat(tscan.pseudo)[,4]


# Extract pseudotime and cell type information
pseudotime_df <- data.frame(
  cell_id = colnames(sce_polyps),
  celltype = sce_polyps$celltypes,  # Ensure this exists
  TSCAN.first = sce_polyps$TSCAN.first,
  TSCAN.second = sce_polyps$TSCAN.second,
  TSCAN.third = sce_polyps$TSCAN.third#,
  #TSCAN.fourth = sce_polyps$TSCAN.fourth
)


pseudotime_df

# Reshape data: Convert wide format to long format
pseudotime_long <- pseudotime_df %>%
  pivot_longer(cols = starts_with("TSCAN"), names_to = "path", values_to = "pseudotime") %>%
  filter(!is.na(pseudotime))  # Remove NA values (cells not assigned to a path)




# Rename paths to match biological interpretation ## healthy samples
pseudotime_long$path <- recode(pseudotime_long$path,
                               "TSCAN.first"  = "Perivascular",
                               "TSCAN.second" = "Endothelial",
                               #"TSCAN.third"  = "Macrophage",
                               "TSCAN.third" = "Ciliated"
)

# Rename paths to match biological interpretation ## polyps samples
pseudotime_long$path <- recode(pseudotime_long$path,
                               "TSCAN.first"  = "Perivascular",
                               "TSCAN.second" = "Ciliated",
                               "TSCAN.third"  = "Endothelial"#,
                               #"TSCAN.fourth" = "Macrophage"
)

ggplot(pseudotime_long, aes(x = reorder(celltype, pseudotime), y = pseudotime, fill = celltype)) +
  geom_boxplot(outlier.shape = NA) +  # Removes outliers for clarity
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~ path, scales = "free_y") +  # Separate plots for each path
  theme_minimal() +
  labs(x = "Cell Types", y = "TSCAN Pseudotime") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Generate the vertically stacked boxplot
ggplot(pseudotime_long, aes(x = reorder(celltype, pseudotime), y = pseudotime, fill = celltype)) +
  geom_boxplot(outlier.shape = NA) +  # Removes outliers for clarity
  scale_fill_brewer(palette = "Set3") +
  facet_grid(path ~ ., scales = "free_y") +  # Stacks plots vertically
  theme_minimal() +
  labs(x = "Cell Types", y = "TSCAN Pseudotime") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0)  # Ensures path labels are horizontal
  )

# Generate the Boxplot with Pseudotime on X-Axis
ggplot(pseudotime_long, aes(x = pseudotime, y = reorder(celltype, pseudotime), fill = celltype)) +
  geom_boxplot(outlier.shape = NA) +  # Removes extreme outliers
  scale_fill_brewer(palette = "Set3") +
  facet_grid(path ~ ., scales = "free_y") +  # Stacks plots vertically for each path
  theme_minimal() +
  labs(x = "TSCAN Pseudotime", y = "Cell Types") +
  theme(
    strip.text.y = element_text(angle = 0),  # Ensures facet labels are horizontal
    legend.position = "none"  # Removes redundant legend
  )

dev.off()

# Remove Perivascular cells from the Ciliated and Endothelial paths
pseudotime_filtered <- pseudotime_long %>%
  filter(!(celltype == "Perivascular" & path %in% c("Ciliated", "Endothelial")))

# Set factor levels in reverse order for bottom-to-top appearance
pseudotime_filtered$path <- factor(pseudotime_filtered$path,
                                   levels = c("Endothelial", "Ciliated", "Perivascular"))



# # Set factor levels in reverse order to control bottom-to-top appearance
# pseudotime_long$path <- factor(pseudotime_long$path,
#                                levels = c("Endothelial", "Ciliated", "Perivascular"))

# Create the updated plot
orderplot_filtered <- ggplot(pseudotime_filtered, aes(x = pseudotime, y = reorder(celltype, pseudotime), fill = celltype)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_brewer(palette = "Set3") +
  facet_grid(path ~ ., scales = "free_y") +
  theme_minimal(base_size = 16) +
  labs(x = "TSCAN Pseudotime", y = "Cell Types") +
  theme(
    strip.text.y = element_text(angle = 0, size = 16),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  )

ggsave(filename = "pseudotime_Cellordering_healthy_noimmune.png", 
       plot = orderplot_filtered, 
       width = 12,   # Increased width
       height = 8,   # Increased height
       dpi = 600,
       bg = "white")


orderplot_filtered

dev.off()

# Create the plot
orderplot <- ggplot(pseudotime_long, aes(x = pseudotime_filtered, y = reorder(celltype, pseudotime), fill = celltype)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_brewer(palette = "Set3") +
  facet_grid(path ~ ., scales = "free_y") +  # Facets ordered bottom-up
  theme_minimal(base_size = 16) +
  labs(x = "TSCAN Pseudotime", y = "Cell Types") +
  theme(
    strip.text.y = element_text(angle = 0, size = 16),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  )

orderplot

ggsave(filename = "pseudotime_Cellordering_polyps_noimmune.png", 
       plot = orderplot, 
       width = 12,   # Increased width
       height = 8,   # Increased height
       dpi = 600,
       bg = "white")

# Define custom colors for the paths
path_colors <- c(
  "Perivascular" = "#d7191c",  # Red
  "Endothelial"  = "#2E8B57",  # Blue
  #"Macrophage"   = "#2b83ba",  # Green
  "Ciliated"     = "#DDA0DD"   # Purple
)

# path_colors <- c(
#   "Path 1" = "#d7191c",  # red
#   "Path 2" = "#2E8B57",  # blue
#   "Path 3" = "#2b83ba",  # orange
#   "Path 4" = "#DDA0DD"   # green
# )

# Generate the single plot
ggplot(pseudotime_long, aes(x = pseudotime, y = reorder(celltype, pseudotime), fill = path)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Slight transparency for clarity
  scale_fill_manual(values = path_colors) +  # Custom colors for each path
  theme_minimal() +
  labs(x = "TSCAN Pseudotime", y = "Cell Types", fill = "Path") +
  theme(
    legend.position = "right",  # Keep legend for path colors
    axis.text.y = element_text(size = 10),  # Adjust text size if needed
    axis.text.x = element_text(size = 10)
  )
