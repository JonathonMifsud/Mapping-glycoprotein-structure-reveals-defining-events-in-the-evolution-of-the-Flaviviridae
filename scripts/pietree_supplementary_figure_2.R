# -----------------------------------------------------------------------------------------
# Supplementary Tree Visualization for "Mapping Glycoprotein Structure Reveals Defining
# Events in the Evolution of the Flaviviridae"
# -----------------------------------------------------------------------------------------
# Authors: Jonathon C.O. Mifsud, Spyros Lytras, Michael Oliver, Kamilla Toon, Vincenzo A. Costa,
#          Edward C. Holmes, Joe Grove
# -----------------------------------------------------------------------------------------
#
# Purpose:
# This script produces the phylogenetics tree for in Supplementary Figure 2
# Much of this script deals with adding pie charts to nodes based on bootstrap and SH-aLRT support values.

# -----------------------------------------------------------------------------------------

# Load necessary libraries
library(ape)
library(phytools)
library(ggtree)
library(tidyverse)
library(ggimage) # For adding images to ggplot objects


setwd("~/Documents/GitHub/phd/p6_flavi_stucture/final/") # Set working directory

# Load the raw Tree 18 in
tree <- ape::read.tree("../flaviviridae_folding_draft_070823/trees/ns5b/flaviviridae_ns5_trimmed_cons5_gt0.9_LG+F+R10_MUSCLE_20230630.treefile")

# Remove two of three tombus outgroup tips tips from the tree for clarity
tombus_tips_to_remove <- c("TOMB_Cucumber_leaf_spot_virus", "TOMB_Turnip_crinkle_virus")
tree_single_tombus <- drop.tip(tree, tombus_tips_to_remove)

# Reroot the tree for better visualization
tombus_tip <- "TOMB_Red_clover_necrotic_mosaic_virus"
tombus_node <- filter(ggtree(tree_single_tombus)$data, label == tombus_tip)$node

# ladderize the tree
tombus_rooted_tree <- ladderize(phytools::reroot(tree_single_tombus, node = tombus_node, position = 1), right = FALSE)

# Initial tree visualization setup
p1 <- ggtree(tombus_rooted_tree) +
  ggtree::geom_tiplab(size = 6, family = "Helvetica") + # add Tip labels
  ggtree::theme_tree() + # Basic tree theme
  theme(text = element_text(family = "Helvetica")) +
  geom_treescale(family = "Helvetica") +
  geom_rootedge(rootedge = 3)

# Prepare data for pie charts based on SHaLRT and UFboot support values
p1_data <- p1$data
p1_data <- p1_data[!p1_data$isTip, ]
p1_data <- separate(data = p1_data, col = label, sep = "/", into = c("SHaLRT", "UFboot"))
p1_data$SHaLRT <- as.numeric(p1_data$SHaLRT)
p1_data$UFboot <- as.numeric(p1_data$UFboot)
p1_data <- p1_data %>%
  mutate(pie = case_when(
    (SHaLRT >= 80 & UFboot >= 95) ~ "B", # full circle
    (SHaLRT < 80 & UFboot >= 95) ~ "R", # half circle coloured black on the right half
    (SHaLRT >= 80 & UFboot < 95) ~ "L", # half circle coloured black on the left half
    (SHaLRT < 80 & UFboot < 95) ~ "N", # White circle
    TRUE ~ ""
  ))

# Create pie chart data, filtering out nodes with low support
# The reason why there is a z in zBlackRight is due to alphabetical ordering effecting plotting order
# specifically we want BlackRight to end up on the right side of the pie graph
# There is surely a better way to do this manually but I am lazy
pie_data <- p1_data %>%
  select(node, pie) %>%
  mutate(
    BlackLeft = ifelse(pie == "B", 1, ifelse(pie == "L", 0.5, ifelse(pie == "R", 0, 0))),
    WhiteRight = ifelse(pie == "B", 0, ifelse(pie == "L", 0.5, ifelse(pie == "R", 0, 1))),
    WhiteLeft = ifelse(pie == "B", 0, ifelse(pie == "L", 0, ifelse(pie == "R", 0.5, 0))),
    zBlackRight = ifelse(pie == "B", 0, ifelse(pie == "L", 0, ifelse(pie == "R", 0.5, 0)))
  ) %>%
  select(node, BlackLeft, zBlackRight, WhiteLeft, WhiteRight, pie) %>%
  column_to_rownames("node")

# Filter out nodes marked as 'N' to avoid adding circles for low support
# Instead we will just show these as blank i.e., no pie
pie_data_filtered <- pie_data %>%
  filter(pie != "N", pie != "")

# Assign colors to each pie chart component
my_palette <- c(
  BlackLeft = "black",
  zBlackRight = "black",
  WhiteLeft = "white",
  WhiteRight = "white"
)

# Add pie charts to the tree visualization
pies <- nodepie(pie_data_filtered, cols = 2:5, outline.color = "black", outline.size = 0.5)
pies <- lapply(pies, function(g) g + scale_fill_manual(values = my_palette))
p2 <- p1 + geom_inset(pies, width = .0125, height = .0125)

# Save the tree visualization with standard labels
ggsave(plot = p2, filename = "flaviviridae_ns5_std_labels.svg", height = 140, width = 50, limitsize = FALSE, dpi = 1000)


# splitting up into clades -----------------------------------------------------
library(readxl)

# Now we have a full family tree lets do the same thing and split the tree up into the 3 subclades
# so that we can better visualize it across multiple pages.

# this is Supplementary Table 1
sequence_summary_table <- read_xlsx("../ns5b_sequence_summary_table.xlsx") %>%
  mutate(tree_tip_with_accession = paste0(tree_tip_label, " ", "(", nucl_sequence_accession, ")"))

orthoflavi <- (sequence_summary_table %>% filter(str_detect(clade, "FJ..")))$tree_tip_with_accession
pestilgf <- (sequence_summary_table %>% filter(str_detect(clade, "PL..")))$tree_tip_with_accession
hepacipegi <- (sequence_summary_table %>% filter(str_detect(clade, "HP..")))$tree_tip_with_accession

orthoflavi_mrca <- findMRCA(tombus_rooted_tree_with_accession, orthoflavi)
pestilgf_mrca <- findMRCA(tombus_rooted_tree_with_accession, pestilgf)
hepacipegi_mrca <- findMRCA(tombus_rooted_tree_with_accession, hepacipegi)

orthoflavi_tree <- extract.clade(tombus_rooted_tree_with_accession, node = orthoflavi_mrca)
orthoflavi_p1_with_accession <- ggtree(orthoflavi_tree) +
  ggtree::geom_tiplab(size = 3, family = "Helvetica") +
  # ggtree::geom_tippoint(size = tipPointSize, ggtree::aes(color = group)) +
  # ggtree::scale_color_manual(values = colors) +
  ggtree::theme_tree() +
  theme(text = element_text(family = "Helvetica")) +
  geom_treescale(family = "Helvetica") +
  geom_rootedge(rootedge = 2)

orthoflavi_data_with_accession <- orthoflavi_p1_with_accession$data
orthoflavi_data_with_accession <- orthoflavi_data_with_accession[!orthoflavi_data_with_accession$isTip, ]
orthoflavi_data_with_accession <- orthoflavi_data_with_accession[!orthoflavi_data_with_accession$label == "Root", ]
orthoflavi_data_with_accession <- separate(data = orthoflavi_data_with_accession, col = label, sep = "/", into = c("SHaLRT", "UFboot"))
orthoflavi_data_with_accession$SHaLRT <- as.numeric(orthoflavi_data_with_accession$SHaLRT)
orthoflavi_data_with_accession$UFboot <- as.numeric(orthoflavi_data_with_accession$UFboot)
orthoflavi_data_with_accession <- orthoflavi_data_with_accession %>%
  mutate(pie = case_when(
    (SHaLRT >= 80 & UFboot >= 95) ~ "B",
    (SHaLRT < 80 & UFboot >= 95) ~ "R",
    (SHaLRT >= 80 & UFboot < 95) ~ "L",
    (SHaLRT < 80 & UFboot < 95) ~ "N",
    TRUE ~ ""
  ))

orthoflavi_data_with_accession <- orthoflavi_data_with_accession[!orthoflavi_data_with_accession$isTip, ] %>%
  select(node, pie) %>%
  mutate(
    BlackLeft = ifelse(pie == "B", 1, ifelse(pie == "L", 0.5, ifelse(pie == "R", 0, 0))),
    WhiteRight = ifelse(pie == "B", 0, ifelse(pie == "L", 0.5, ifelse(pie == "R", 0, 1))),
    WhiteLeft = ifelse(pie == "B", 0, ifelse(pie == "L", 0, ifelse(pie == "R", 0.5, 0))),
    zBlackRight = ifelse(pie == "B", 0, ifelse(pie == "L", 0, ifelse(pie == "R", 0.5, 0))),
    node2 = node
  ) %>%
  select("node", "BlackLeft", "zBlackRight", "WhiteLeft", "WhiteRight", "pie", "node2") %>%
  column_to_rownames("node2")

# As this is a massive tree I don't really want to include circles for those with low support i.e., N
# So lets subset the dataset
orthoflavi_data_with_accession_filtered <- orthoflavi_data_with_accession %>%
  filter(!pie == "N") %>%
  filter(!pie == "")

# assign colours for each va;ue
my_palette <- c(
  BlackLeft = "black",
  zBlackRight = "black",
  WhiteLeft = "white",
  WhiteRight = "white"
)

orthoflavi_pies_with_accession <- nodepie(orthoflavi_data_with_accession_filtered, cols = 2:5, outline.color = "black", outline.size = 0.4)
orthoflavi_pies_with_accession <- lapply(orthoflavi_pies_with_accession, function(g) g + scale_fill_manual(values = my_palette))
orthoflavi_p2_with_accession <- orthoflavi_p1_with_accession + geom_inset(orthoflavi_pies_with_accession, width = .023, height = .023, hjust = 0.005, vjust = 0.11)
# dim are a1 width /2
ggsave(plot = orthoflavi_p2_with_accession, filename = "orthoflavi_ns5_with_accessions.svg", height = 23.4, width = 16.5, dpi = 1000)

# PESTI ------------------------------------------------------------------------------
pestilgf_tree <- extract.clade(tombus_rooted_tree_with_accession, node = pestilgf_mrca)
pestilgf_p1_with_accession <- ggtree(pestilgf_tree) +
  ggtree::geom_tiplab(size = 4, family = "Helvetica") +
  # ggtree::geom_tippoint(size = tipPointSize, ggtree::aes(color = group)) +
  # ggtree::scale_color_manual(values = colors) +
  ggtree::theme_tree() +
  theme(text = element_text(family = "Helvetica")) +
  geom_treescale(family = "Helvetica") +
  geom_rootedge(rootedge = 3)

pestilgf_data_with_accession <- pestilgf_p1_with_accession$data
pestilgf_data_with_accession <- pestilgf_data_with_accession[!pestilgf_data_with_accession$isTip, ]
pestilgf_data_with_accession <- pestilgf_data_with_accession[!pestilgf_data_with_accession$label == "Root", ]
pestilgf_data_with_accession <- separate(data = pestilgf_data_with_accession, col = label, sep = "/", into = c("SHaLRT", "UFboot"))
pestilgf_data_with_accession$SHaLRT <- as.numeric(pestilgf_data_with_accession$SHaLRT)
pestilgf_data_with_accession$UFboot <- as.numeric(pestilgf_data_with_accession$UFboot)
pestilgf_data_with_accession <- pestilgf_data_with_accession %>%
  mutate(pie = case_when(
    (SHaLRT >= 80 & UFboot >= 95) ~ "B",
    (SHaLRT < 80 & UFboot >= 95) ~ "R",
    (SHaLRT >= 80 & UFboot < 95) ~ "L",
    (SHaLRT < 80 & UFboot < 95) ~ "N",
    TRUE ~ ""
  ))

# The z in zBlackRight is just used to effect the plotting order so that BlackRight
# ends up on the right side of the pie graph
# there is a way to do this manually but I am lazy
pestilgf_data_with_accession <- pestilgf_data_with_accession[!pestilgf_data_with_accession$isTip, ] %>%
  select(node, pie) %>%
  mutate(
    BlackLeft = ifelse(pie == "B", 1, ifelse(pie == "L", 0.5, ifelse(pie == "R", 0, 0))),
    WhiteRight = ifelse(pie == "B", 0, ifelse(pie == "L", 0.5, ifelse(pie == "R", 0, 1))),
    WhiteLeft = ifelse(pie == "B", 0, ifelse(pie == "L", 0, ifelse(pie == "R", 0.5, 0))),
    zBlackRight = ifelse(pie == "B", 0, ifelse(pie == "L", 0, ifelse(pie == "R", 0.5, 0))),
    node2 = node
  ) %>%
  select("node", "BlackLeft", "zBlackRight", "WhiteLeft", "WhiteRight", "pie", "node2") %>%
  column_to_rownames("node2")

# As this is a massive tree I don't really want to include circles for those with low support i.e., N
# So lets subset the dataset
pestilgf_data_with_accession_filtered <- pestilgf_data_with_accession %>%
  filter(!pie == "N") %>%
  filter(!pie == "")

# assign colours for each va;ue
my_palette <- c(
  BlackLeft = "black",
  zBlackRight = "black",
  WhiteLeft = "white",
  WhiteRight = "white"
)

pestilgf_pies_with_accession <- nodepie(pestilgf_data_with_accession_filtered, cols = 2:5, outline.color = "black", outline.size = 0.4)
pestilgf_pies_with_accession <- lapply(pestilgf_pies_with_accession, function(g) g + scale_fill_manual(values = my_palette))
pestilgf_p2_with_accession <- pestilgf_p1_with_accession + geom_inset(pestilgf_pies_with_accession, width = .025, height = .025, hjust = 0.0055, vjust = 0.11)
# dim are a1 width /2
ggsave(plot = pestilgf_p2_with_accession, filename = "pestilgf_ns5_with_accessions.svg", height = 23.4, width = 16.5, dpi = 1000)


# HEPACI -----------------------------------------------------------------------

hepacipegi_tree <- extract.clade(tombus_rooted_tree_with_accession, node = hepacipegi_mrca)
hepacipegi_p1_with_accession <- ggtree(hepacipegi_tree) +
  ggtree::geom_tiplab(size = 3.5, family = "Helvetica") +
  # ggtree::geom_tippoint(size = tipPointSize, ggtree::aes(color = group)) +
  # ggtree::scale_color_manual(values = colors) +
  ggtree::theme_tree() +
  theme(text = element_text(family = "Helvetica")) +
  geom_treescale(family = "Helvetica") +
  geom_rootedge(rootedge = 3)

hepacipegi_data_with_accession <- hepacipegi_p1_with_accession$data
hepacipegi_data_with_accession <- hepacipegi_data_with_accession[!hepacipegi_data_with_accession$isTip, ]
hepacipegi_data_with_accession <- hepacipegi_data_with_accession[!hepacipegi_data_with_accession$label == "Root", ]
hepacipegi_data_with_accession <- separate(data = hepacipegi_data_with_accession, col = label, sep = "/", into = c("SHaLRT", "UFboot"))
hepacipegi_data_with_accession$SHaLRT <- as.numeric(hepacipegi_data_with_accession$SHaLRT)
hepacipegi_data_with_accession$UFboot <- as.numeric(hepacipegi_data_with_accession$UFboot)
hepacipegi_data_with_accession <- hepacipegi_data_with_accession %>%
  mutate(pie = case_when(
    (SHaLRT >= 80 & UFboot >= 95) ~ "B",
    (SHaLRT < 80 & UFboot >= 95) ~ "R",
    (SHaLRT >= 80 & UFboot < 95) ~ "L",
    (SHaLRT < 80 & UFboot < 95) ~ "N",
    TRUE ~ ""
  ))

# The z in zBlackRight is just used to effect the plotting order so that BlackRight
# ends up on the right side of the pie graph
# there is a way to do this manually but I am lazy
hepacipegi_data_with_accession <- hepacipegi_data_with_accession[!hepacipegi_data_with_accession$isTip, ] %>%
  select(node, pie) %>%
  mutate(
    BlackLeft = ifelse(pie == "B", 1, ifelse(pie == "L", 0.5, ifelse(pie == "R", 0, 0))),
    WhiteRight = ifelse(pie == "B", 0, ifelse(pie == "L", 0.5, ifelse(pie == "R", 0, 1))),
    WhiteLeft = ifelse(pie == "B", 0, ifelse(pie == "L", 0, ifelse(pie == "R", 0.5, 0))),
    zBlackRight = ifelse(pie == "B", 0, ifelse(pie == "L", 0, ifelse(pie == "R", 0.5, 0))),
    node2 = node
  ) %>%
  select("node", "BlackLeft", "zBlackRight", "WhiteLeft", "WhiteRight", "pie", "node2") %>%
  column_to_rownames("node2")

# As this is a massive tree I don't really want to include circles for those with low support i.e., N
# So lets subset the dataset
hepacipegi_data_with_accession_filtered <- hepacipegi_data_with_accession %>%
  filter(!pie == "N") %>%
  filter(!pie == "")

# assign colours for each va;ue
my_palette <- c(
  BlackLeft = "black",
  zBlackRight = "black",
  WhiteLeft = "white",
  WhiteRight = "white"
)

hepacipegi_pies_with_accession <- nodepie(hepacipegi_data_with_accession_filtered, cols = 2:5, outline.color = "black", outline.size = 0.4)
hepacipegi_pies_with_accession <- lapply(hepacipegi_pies_with_accession, function(g) g + scale_fill_manual(values = my_palette))
hepacipegi_p2_with_accession <- hepacipegi_p1_with_accession + geom_inset(hepacipegi_pies_with_accession, width = .027, height = .027, hjust = 0.0055, vjust = 0.11)
# dim are a1 width /2
ggsave(plot = hepacipegi_p2_with_accession, filename = "hepacipegi_ns5_with_accessions.svg", height = 23.4, width = 16.5, dpi = 1000)
