# -----------------------------------------------------------------------------------------
# Mapping Glycoprotein Structure Reveals Defining Events in the Evolution of the Flaviviridae
# -----------------------------------------------------------------------------------------
# Authors: Jonathon C.O. Mifsud, Spyros Lytras, Michael Oliver, Kamilla Toon, Vincenzo A. Costa,
#          Edward C. Holmes, Joe Grove
# -----------------------------------------------------------------------------------------
#
# Script Purpose:
# This script is designed take the ns5b trees and run the treespace analysis 
# including the generation of Supplementary Figure 1.
#
# The script also generates a *BONUS* plot not included in the manuscript but is 
# a fun way to explore the varitation between all the tree varitations
#
# -----------------------------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(treespace)
library(adegenet)
library(adegraphics)
library(rgl)
library(ape)
library(phytools)
library(factoextra)
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(cowplot)
library(readxl)

# Create an empty list to store plots
plot_list <- list()

setwd("/Users/jonmifsud/Documents/GitHub/phd/p6_flavi_stucture/final/for_upload/trees/ns5b/all_variations/")

################################################################################
#                   Analysis for phylo set with ClustalO                       #                                   
################################################################################

outgroup <- c("TOMB_Cucumber_leaf_spot_virus", "TOMB_Red_clover_necrotic_mosaic_virus", "TOMB_Turnip_crinkle_virus")

# Initialize empty list to store the phylogenies
phylogenies <- list()

# Set the directory path where the phylogeny files are located
dir_path <- "./"

# Get the file names of all the phylogeny files in the directory
file_names <- list.files(dir_path, pattern = "\\.treefile$", full.names = TRUE)

# Loop through the file names and read the phylogenies
for (file_name in file_names) {
  # Read the phylogeny from file
  phy <- read.tree(file_name)
  phy_rooted <- midpoint.root(phy)
  node <- findMRCA(phy_rooted, tips = outgroup)
  phy_rooted_outgroup <- root(phy_rooted, node = node, resolve.root = TRUE)
  
  # Add the phylogeny to the list
  phylogenies[[file_name]] <- phy_rooted_outgroup
}

class(phylogenies)<-"multiPhylo"

for (i in 1:length(phylogenies)) {
  # Access the i-th phylogeny
  phylo <- phylogenies[[i]]
  
  # Check if all node.labels are populated (non-empty)
  if (all(nchar(phylo$node.label) > 0)) {
    print(paste("Phylogeny", i, "has all node labels populated."))
  } else {
    print(paste("Phylogeny", i, "does not have all node labels populated."))
  }
}

# Phylogeny name cleaning
new_names <- str_remove(names(phylogenies),
                         ".//flaviviridae_ns5_")
new_names <- str_remove(new_names,
                        "_20230630.treefile")
new_names <- str_remove(new_names,
                        "trimmed_")
new_names <- str_replace_all(new_names,
                        "MUSCLE", "MU")
new_names <- str_replace_all(new_names,
                             "MAFFT", "MA")
new_names <- str_replace_all(new_names,
                             "ClustalOmega", "CO")
new_names <- str_replace_all(new_names,
                             "gappyout", "go")
new_names <- str_replace_all(new_names,
                             "_trimmed", "")
new_names <- str_replace_all(new_names,
                             ".fasta", "")
new_names <- str_replace_all(new_names,
                             "cons175", "cons17.5")
new_names <- str_replace_all(new_names,
                             "cons75", "cons7.5")
new_names <- str_replace_all(new_names,
                             "cons125", "cons12.5")
# Rename the list elements with the updated names
names(phylogenies) <- new_names

# here we use treespace to compare the trees using treeVec a Kendall Colijn metric vector (for rooted trees)
# in this case we retain 2 principle components 
res <- treespace(phylogenies, nf=2, processors = 8)

# here we plot out the trees in to examine their distribution interactively.
plotGrovesD3(res$pco, treeNames=names(phylogenies))

# make a df of the phylogeny construction information e.g., aligner, trimal consenus, gapthresold and aa models
row_names <- tibble(full = row.names(res$pco$l1)) %>% 
  mutate(cons = word(full, 1, sep = "_"),
         gt = word(full, 2, sep = "_"),
         model =  word(full, 3, sep = "_"),
         aliger = word(full, 4, sep = "_"))

# turn phylo information into factors to use in plots.
row_names$full <- factor(row_names$full, levels = unique(row_names$full))
row_names$cons <- factor(row_names$cons, levels = unique(row_names$cons))
row_names$gt <- factor(row_names$gt, levels = unique(row_names$gt))
row_names$model <- factor(row_names$model, levels = unique(row_names$model))
row_names$aliger <- factor(row_names$aliger, levels = unique(row_names$aliger))

# Here we plot out the same thing as above but we are colouring by aligner
# We highlight that Mafft and Mucle generally cluster together while Clustalo
# has outliers
plotGrovesD3(res$pco, groups=row_names$aliger)

################################################################################
#                   Analysis for phylo set without ClustalO                    #                                   
################################################################################
# Here we rerun the analysis without clustalo
setwd("without_clustal/")
dir_path2 <- "./"
outgroup <- c("TOMB_Cucumber_leaf_spot_virus", "TOMB_Red_clover_necrotic_mosaic_virus", "TOMB_Turnip_crinkle_virus") 
file_names <- list.files(dir_path2, pattern = "\\.treefile$", full.names = TRUE)
phylogenies_without_clustal <- list()

# Loop through the file names and read the phylogenies
for (file_name in file_names) {
  
  # Read the phylogeny from file
  phy <- read.tree(file_name)
  phy_rooted <- midpoint.root(phy)
  node <- findMRCA(phy_rooted, tips = outgroup)
  
  # root on tombus
  phy_rooted_outgroup <- root(phy_rooted, node = node, resolve.root = TRUE)
  
  # Add the phylogeny to the list
  phylogenies_without_clustal[[file_name]] <- phy_rooted_outgroup
}

class(phylogenies_without_clustal)<-"multiPhylo"

# Clean phylo names
new_names_without_clustal <- str_remove(names(phylogenies_without_clustal),
                        "./.flaviviridae_ns5_")
new_names_without_clustal <- str_remove(new_names_without_clustal,
                        "_20230630.treefile")
new_names_without_clustal <- str_remove(new_names_without_clustal,
                        "trimmed_")
new_names_without_clustal <- str_replace_all(new_names_without_clustal,
                             "MUSCLE", "MU")
new_names_without_clustal <- str_replace_all(new_names_without_clustal,
                             "MAFFT", "MA")
new_names_without_clustal <- str_replace_all(new_names_without_clustal,
                             "ClustalOmega", "CO")
new_names_without_clustal <- str_replace_all(new_names_without_clustal,
                             "gappyout", "go")
new_names_without_clustal <- str_replace_all(new_names_without_clustal,
                             "_trimmed", "")
new_names_without_clustal <- str_replace_all(new_names_without_clustal,
                             ".fasta", "")
new_names_without_clustal <- str_replace_all(new_names_without_clustal,
                             "cons175", "cons17.5")
new_names_without_clustal <- str_replace_all(new_names_without_clustal,
                             "cons75", "cons7.5")
new_names_without_clustal <- str_replace_all(new_names_without_clustal,
                             "cons125", "cons12.5")
# Rename the list elements with the updated names
names(phylogenies_without_clustal) <- new_names_without_clustal

# use treespace
res_without_clustal <- treespace(phylogenies_without_clustal, nf=2, processors = 8)
row_names_without_clustal <- tibble(full = row.names(res_without_clustal$pco$l1)) %>% 
  mutate(cons = word(full, 1, sep = "_"),
         gt = word(full, 2, sep = "_"),
         model =  word(full, 3, sep = "_"),
         aliger = word(full, 4, sep = "_"))

# assign parameters to factors
row_names_without_clustal$full <- factor(row_names_without_clustal$full, levels = unique(row_names_without_clustal$full))
row_names_without_clustal$cons <- factor(row_names_without_clustal$cons, levels = unique(row_names_without_clustal$cons))
row_names_without_clustal$gt <- factor(row_names_without_clustal$gt, levels = unique(row_names_without_clustal$gt))
row_names_without_clustal$model <- factor(row_names_without_clustal$model, levels = unique(row_names_without_clustal$model))
row_names_without_clustal$aliger <- factor(row_names_without_clustal$aliger, levels = unique(row_names_without_clustal$aliger))
row_names_without_clustal <- row_names_without_clustal %>% 
  mutate(opacity = case_when(cons == "cons5" ~ 0.3,
                             cons == "cons7.5" ~ 0.7,
                             TRUE ~ 1))

# Here we use the hierarchical clustering on principal components output to identity 4 clusters
# The three clusters look good but the 4th is different. 
phylo_groves_without_clustal <- findGroves(res_without_clustal, nclust = 4)

# Lets overlay some of the parameters (e.g, aligner) as shapes
# aligner
plotGrovesD3(phylo_groves_without_clustal, tooltip_text=names(phylogenies_without_clustal), legend_width=100, col_lab="Cluster", symbol_var = row_names_without_clustal$aliger)
# aa sub model
plotGrovesD3(phylo_groves_without_clustal, tooltip_text=names(phylogenies_without_clustal), legend_width=100, col_lab="Cluster", symbol_var = row_names_without_clustal$model)
# trimal gap threshold
plotGrovesD3(phylo_groves_without_clustal, tooltip_text=names(phylogenies_without_clustal), legend_width=100, col_lab="Cluster", symbol_var = row_names_without_clustal$gt)

# trimal cons value
# there are more conservation values than there are shapes in R so I added in an opacity column so that these can be told apart
# then in post R I've replaced these with different shapes
plotGrovesD3(phylo_groves_without_clustal, tooltip_text=names(phylogenies_without_clustal), legend_width=100, col_lab="Cluster", symbol_var = row_names_without_clustal$cons,
             opacity_var = row_names_without_clustal$opacity)

# Create median trees for each of the clusters
combined_without_clustal <- interaction(phylo_groves_without_clustal$groups, drop = TRUE)
med_trees <- medTree(phylogenies_without_clustal, combined_without_clustal)
med_trees_extracted <- lapply(med_trees, function(e) ladderize(e$trees[[1]]))

# We can do the same but create a median tree for the interaction between aligner and cluster
# combined_without_clustal <- interaction(phylo_groves_without_clustal$groups, row_names_without_clustal$aliger, drop = TRUE)
# med_trees <- medTree(phylogenies_without_clustal, combined_without_clustal)
# med_trees_extracted <- lapply(med_trees, function(e) ladderize(e$trees[[i]]))

# plot trees
par(mfrow=c(2,2))
for(i in 1:length(med_trees_extracted)) plot.phylo(med_trees_extracted[[i]], main=paste(names(med_trees_extracted)[i]),cex=1.5, show.tip.label = FALSE)

# write the trees
# for(i in 1:length(med_trees_extracted)) write.tree(med_trees_extracted[[i]], paste0("med_tree_without_clustal_", i, ".newick"))

#saveRDS(med_trees_extracted, "median_trees.RDS")
#saveRDS(phylo_groves_without_clustal, "phylo_groves_without_clustal.RDS")
#saveRDS(res_without_clustal, "res_without_clustal.RDS")
#saveRDS(phylo_groves, "phylo_groves.RDS")
#saveRDS(res, "res.RDS")


# Extra plots ==================================================================

library(ggtree)
library(phytools)

 #If you would like to visualise all of the tree varitations and how they differ
 #in the major clades this is the script section for you
 #There wasnt any room for this figure in the paper but it can be found on
 #Zenodo /trees/ns5b/all_varitations/facet_phylo_plots.pdf
setwd("/Users/jonmifsud/Documents/GitHub/phd/p6_flavi_stucture/final/for_upload/trees/ns5b/all_variations/")
tips <- read.tree("../flaviviridae_ns5_trimmed_cons5_gt0.9_LG+F+R10_MUSCLE_20230630.treefile")
tips <- as_tibble(tips$tip.label) %>% mutate(clade = word(value, 1, sep = "_"))

# Get the file names of all the phylogeny files in the directory
# let change the working dir to the all_varitations (inc clustal) folder

file_names <- list.files(dir_path, pattern = "\\.treefile$", full.names = TRUE)

# make a df of the phylo information e.g., aligner, trimal consenus
row_names <- tibble(full = row.names(res$pco$l1)) %>% 
  mutate(cons = word(full, 1, sep = "_"),
         gt = word(full, 2, sep = "_"),
         model =  word(full, 3, sep = "_"),
         aliger = word(full, 4, sep = "_"))

row_names$orig_name <- file_names


orthoflavi <- tips %>% filter(str_detect(value, "FJ..")) %>% filter(!str_detect(clade, "FJUN"))
jingmen <- (orthoflavi %>% filter(str_detect(clade, "FJJI")))$value
ortho <- (orthoflavi %>% filter(!str_detect(clade, "FJJI")))$value

pestilgf <- tips %>% filter(str_detect(clade, "PL.."))%>% filter(!str_detect(clade, "PLUN"))
pesti <- (pestilgf %>% filter(str_detect(clade, "PLPV")))$value
lgf <- (pestilgf %>% filter(str_detect(clade, "PLLG")))$value

hepacipegi <- tips %>% filter(str_detect(clade, "HP.."))%>% filter(!str_detect(clade, "HPUN"))
hepaci <- (hepacipegi %>% filter(str_detect(clade, "HPHV")))$value
pegi <- (hepacipegi %>% filter(str_detect(clade, "HPPV")))$value

groupings <- list(jingmen=jingmen,
                  ortho=ortho,
                  pesti=pesti,
                  lgf=lgf,
                  hepaci=hepaci,
                  pegi=pegi)

# Define colors for groups
group_colors <- c(jingmen = "#9CD6F2",
                  ortho = "#0F783C",
                  pesti = "#43AA99",
                  lgf = "#DCCC77",
                  hepaci = "#892155",
                  pegi = "#CD6777",
                  Other = "grey")

# Loop over each tree and create a plot with smaller size and add title
for (i in 1:nrow(row_names)) {
  tree_file <- row_names$full[i]  # Adjust if it's just an identifier
  phy <- read.tree(file = row_names$orig_name[i])  # Adjust the file path and extension
  phy_rooted <- midpoint.root(phy)
  node <- findMRCA(phy_rooted, tips = outgroup)
  phy_rooted_outgroup <- root(phy_rooted, node = node, resolve.root = TRUE)
  tree_grouped <- groupOTU(phy_rooted_outgroup, groupings)
  
  # Convert group to factor and use group_colors for coloring
  p <- ggtree(tree_grouped, size = 0.5, aes(color = group)) + 
    scale_color_manual(values = group_colors) +
    labs(title = tree_file) +  # Add tree name as title
    theme_void() +
    theme(plot.title = element_text(size = 8, hjust = 0.5),
          plot.margin = margin(10, 10, 10, 10),
          legend.position = "none")  # Hide the legend
  plot_list[[i]] <- list(plot = p, aliger = row_names$aliger[i])
}

# Sort plot_list based on aliger
plot_list <- plot_list[order(sapply(plot_list, function(x) x$aliger))]

# Define the number of plots per page
plots_per_page <- 15

# Function to create a page of plots
create_page <- function(plot_list, start_index, end_index) {
  subset_plots <- lapply(plot_list[start_index:end_index], function(x) x$plot)
  return(plot_grid(plotlist = subset_plots, ncol = 5, align = 'v'))  # Adjust ncol as needed
}

# Open a PDF device for multiple pages
pdf("facet_phylo_plots.pdf", width = 8.3, height = 11.7)

# Create and save each page
num_pages <- ceiling(length(plot_list) / plots_per_page)
for (page in 1:num_pages) {
  start_index <- (page - 1) * plots_per_page + 1
  end_index <- min(start_index + plots_per_page - 1, length(plot_list))
  print(create_page(plot_list, start_index, end_index))
}

# Close the PDF device
dev.off()


# Code Graveyard ===============================================================

# Figure plot B - full treespace analysis but for all trees including clustalO.
# this was ran prior to the decision to remove clustalo alignments

# phylo_groves <- findGroves(res, nclust = 3)
# plotGrovesD3(phylo_groves, tooltip_text=names(phylogenies), legend_width=100, col_lab="Cluster", symbol_var = row_names$aliger)
# plotGrovesD3(phylo_groves, tooltip_text=names(phylogenies), legend_width=100, col_lab="Cluster", symbol_var = row_names$model)
# plotGrovesD3(phylo_groves, tooltip_text=names(phylogenies), legend_width=100, col_lab="Cluster", symbol_var = row_names$gt)
# plotGrovesD3(phylo_groves, tooltip_text=names(phylogenies), legend_width=100, col_lab="Cluster", symbol_var = row_names$cons)

# Combined group and aligner phylogenies
# combined <- interaction(row_names$aliger, phylo_groves$groups, drop = TRUE)
# med_trees <- medTree(phylogenies, combined)
# med_trees_extracted <- lapply(med_trees, function(e) ladderize(e$trees[[1]]))

# plot trees
# par(mfrow=c(2,7))
# for(i in 1:length(med_trees_extracted)) plot.phylo(med_trees_extracted[[i]], main=paste(names(med_trees_extracted)[i]),cex=1.5, show.tip.label = FALSE)


# Combined group and aligner phylogenies
#med_trees_group <- medTree(phylogenies, phylo_groves$groups)
#med_trees_group_extracted <- lapply(med_trees_group, function(e) ladderize(e$trees[[1]]))

# plot trees
#par(mfrow=c(2,3))
#for(i in 1:length(med_trees_group_extracted)) plot.phylo(med_trees_group_extracted[[i]], main=paste(names(med_trees_group_extracted)[i]),cex=1.5, show.tip.label = FALSE)
