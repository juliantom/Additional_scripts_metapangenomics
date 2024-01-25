# Install and load necessary packages
#install.packages(c("ggplot2", "dplyr", "tidyr", "scales"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(scales)
library(ape)
library(ggtree)

# ADJUST the following lines as needed
############################################

# Variables for input files (matrix and newick files are expected to be in the same path)
path_to_files="/PATH/TO/INPUT_FILES" # Path to files
ani_file <- "FILENAME_ANI_MATRIX.txt" # Should be a TSV file with header (row 1) and rownames (column 1)
newick_file <- "FILENAME_ANI_MATRIX.newick" # Should be a newick tree

# Variables for ouput plots
path_to_output="/PATH/TO/OUTPUT" # Path
plot_file <- "PLOT_NAME.png" # Plot name with no dendogram
plot_file_dendrogram <- "PLOT_NAME_DENDROGRAM.png" # Plot name with dendogram

##########################################

# The following code contains the color palette that will be used for the heatmap values.
# It is always suggested to use not more than 10 color gradients.
# display.brewer.pal(n = 10, name = 'RdBu')
# brewer.pal(n = 10, name = 'RdBu')
# "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"

# Read table as matrix
full_path_ani_matrix=full_path <- paste(path_to_files, ani_file, sep = "/")
ani_matrix <- as.matrix(read.csv(full_path_ani_matrix, sep = "\t", header = TRUE, row.names = 1))

# Read the dendrogram from the Newick file
full_path_ani_tree=full_path <- paste(path_to_files, newick_file, sep = "/")
dendrogram <- read.tree(full_path_ani_tree)

# Melt matrix so that headers and rownames become "variables" (Var1 and Var2) and percent of identity "values"
melted_cormat <- melt(ani_matrix)

# Transform values to percentages
melted_cormat <- melted_cormat %>% mutate(value = value * 100)

# Extract labels from dendrogram
dendrogram_labels <- dendrogram$tip.label

# Identify common labels between dendrogram and melted_cormat
common_labels <- intersect(dendrogram_labels, unique(c(melted_cormat$Var1, melted_cormat$Var2)))

# Reorder the melted_cormat based on common labels to ensure the plot will be a specular image in diagonal
melted_cormat <- melted_cormat %>%
  filter(Var1 %in% common_labels, Var2 %in% common_labels) %>%
  mutate(Var1 = factor(Var1, levels = rev(dendrogram_labels)),
         Var2 = factor(Var2, levels = rev(dendrogram_labels))) %>%
  arrange(Var1, Var2)

# Specify 10 breaks for the color gradient for value 90%-100%
custom_breaks <- seq(90, 100, length.out = 3)

# Create a tree object from the dendrogram
tree <- as.phylo(dendrogram)

# Create the heatmap plot with ggplot
p <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    scale_fill_gradientn(colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"), 
                         breaks = custom_breaks, 
                         na.value = "#f0f0f0", 
                         limits = c(90, 100)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 2, hjust = 0, color = "black"), 
          axis.text.y = element_text(angle = 0, vjust = 0.5, size = 2, hjust = 0, color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          axis.title = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
    labs( fill = "ANI (%)") +
    coord_fixed()

# Plot heatmap
p

# Save the heatmap as a PNG file with 300 dpi
full_path_plot=full_path <- paste(path_to_output, plot_file, sep = "/")
ggsave(full_path_plot, plot = p, dpi = 300)

# Create a tree object from the dendrogram to be added to the plot
tree <- as.phylo(dendrogram)

# Add dendrogram to the heatmap plot
p2 <- ggtree(tree) + geom_tree() + p

# Plot heatmap with tree
p2

# Save the heatmap + dendrogram as a PNG file with 300 dpi
full_path_plot_dendrogram=full_path <- paste(path_to_output, plot_file_dendrogram, sep = "/")
ggsave(full_path_plot_dendrogram, plot = p2, dpi = 300)
