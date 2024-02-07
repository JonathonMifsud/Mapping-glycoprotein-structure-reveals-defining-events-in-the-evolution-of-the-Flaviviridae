# -----------------------------------------------------------------------------------------
# Mapping Glycoprotein Structure Reveals Defining Events in the Evolution of the Flaviviridae
# -----------------------------------------------------------------------------------------
# Authors: Jonathon C.O. Mifsud, Spyros Lytras, Michael Oliver, Kamilla Toon, Vincenzo A. Costa,
#          Edward C. Holmes, Joe Grove
# -----------------------------------------------------------------------------------------
#
# Script Purpose:
# This script is designed to obtain host information for Flaviviridae sequences.
# Used in Figure 2C
# It utilizes functions from a custom R function collection (in development) named Phyparse.
#
# Note:
# This script is not functional without access to the Phyparse package, which contains
# custom functions used throughout. If you are interested in the early code of Phyparse,
# please get in contact
# -----------------------------------------------------------------------------------------

library(tidyverse)
library(readxl)

# Load custom Phyparse package
devtools::load_all("~/Documents/GitHub/phyParse/")

# Set working directory to project folder
setwd("~/Documents/GitHub/phd/p6_flavi_stucture/")

# Reading sequence summary table and preprocessing
df <- read_xlsx("flaviviridae_folding_draft_070823/tables/sequence_summary_table.xlsx")

# Extract distinct tree tip labels, append accession numbers, and select relevant columns
tree_tips <- df %>%
  distinct(tree_tip_label, .keep_all = TRUE) %>%
  mutate(tree_tips_with_accessions = paste0(nucl_sequence_accession, "_", tree_tip_label)) %>%
  select(tree_tips_with_accessions, nucl_sequence_accession, tree_tip_label)

# Validate nucleotide accessions and join to phylogeny tips
checked_accessions <- phyparse::check_accessions(tree_tips$nucl_sequence_accession,
  tip_labels = tree_tips,
  join_results_to_tips = TRUE,
  entrez_api_key = "your_entrez_api_key",
  accession_col = "nucl_sequence_accession",
  tree_tip_label_col = "tree_tips_with_accessions"
)

# Pull NCBI summary information for each sequence
summary_info <- phyparse::pull_ncbi_summary_information(checked_accessions, entrez_api_key = "your_entrez_api_key")

# Check and standardize host names using NCBI information
host_names <- phyparse::check_host_names(summary_info, entrez_api_key = "your_entrez_api_key")

# Rename columns for clarity and consistency
checked_accessions_renamed <- checked_accessions %>%
  rename(orig_label = tree_tip_labels)

# Extract taxonomic lineage for each host and join with tip information
host_names_with_lineage <- phyparse::pull_lineage(
  host_table = host_names,
  host_id_column = "id",
  host_column = "host",
  organism_column = "organism",
  join_to_tips = TRUE,
  virus_accession_column = "accession",
  checked_accessions = checked_accessions_renamed
)


# Define a vector of novel hostnames to check
hostnames <- c(
  "Rhipicephalus_sanguineus", "Insects", "Bemisia_tabaci", "Gecarcinus_lateralis",
  "Locusta_migratoria", "Semibalanus_balanoides", "Artemia_sinica",
  "Phenacoccus_solenopsis", "Plodia_interpunctella", "Brachymeria_lasus",
  "Parasteatoda_tepidariorum"
)

# We do the same thing but for those sequences that are novel / aren't on GenBank
# Check and retrieve lineage for novel hosts
novel_hosts <- phyparse::check_host_names_simple(hostnames, entrez_api_key = "your_entrez_api_key")
novel_host_lineage <- phyparse::pull_lineage(novel_hosts)
write.csv(novel_host_lineage, "novel_host_lineage.csv")
