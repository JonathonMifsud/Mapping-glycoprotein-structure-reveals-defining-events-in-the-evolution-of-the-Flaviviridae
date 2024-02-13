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
#
# Parts of this script are heavily based upon Liam Brierley wonderful blog post examining
# Coronaviruses and their hosts https://ropensci.org/blog/2020/11/10/coronaviruses-and-hosts/
# https://doi.org/10.59350/qenh9-cyj40
# -----------------------------------------------------------------------------------------

library(tidyverse)
library(rentrez)
library(taxize)
library(readxl)

# Reading sequence summary table and preprocessing
df <- read_xlsx("Supplementary_table_1.xlsx", skip = 1)

# Extract distinct tree tip labels, append accession numbers, and select relevant columns
tree_tips <- df %>%
  distinct(tree_tip_label, .keep_all = TRUE) %>%
  mutate(tree_tips_with_accessions = paste0(nucl_sequence_accession, "_", tree_tip_label)) %>%
  select(tree_tips_with_accessions, nucl_sequence_accession, tree_tip_label) %>% 
  filter(!str_detect(nucl_sequence_accession, "hou_et_al|novel")) # remove tips that have no GenBank data associated

accessions <- tree_tips$nucl_sequence_accession

# lets store results here
nucl_seq_result <- vector("list", length(accessions))

# Split the accessions into batches of 50 (NCBI API limit?)
query_index <- split(seq(1, length(accessions)), ceiling(seq_along(seq(1, length(accessions))) / 50))

for (i in 1:length(query_index)) {
  suppressWarnings(
    nucl_seq_result[[i]] <- rentrez::entrez_summary(
      db = "nucleotide",
      id = accessions[unlist(query_index[i])],
      api_key = entrez_api_key
    )
  )
  Sys.sleep(1)
}

# Combine results if there are more than one batch which there is
if (length(accessions) == 1) {
  nucl_seq_result <- nucl_seq_result
} else {
  nucl_seq_result <- nucl_seq_result %>%
    purrr::flatten() %>%
    unname() %>%
    structure(class = c("list", "esummary_list"))
}

summary_info <- nucl_seq_result

# define two functions to parse esummary results 
# for the purposes of extracting host information
extract_subfeatures <- function(esummary_records) {
  
  # Check that the input is a list of esummary records
  if (!is.list(esummary_records) || length(esummary_records) == 0) {
    stop("ERROR: Input must be a non-empty list of esummary records.")
  }
  
  # Extract the accession and subtype information from the esummary records
  all_accession <- rentrez::extract_from_esummary(esummary_records, elements = c("accessionversion"))
  all_subname_list <- rentrez::extract_from_esummary(esummary_records, elements = c("subname", "subtype"), simplify = FALSE)
  
  # Convert the subtype information to a data frame and set column names
  all_subname_list <- all_subname_list %>%
    lapply(function(x) {
      subname_row <- x["subname"] %>%
        stringr::str_split("\\|") %>%
        unlist() %>%
        matrix(nrow = 1, byrow = FALSE) %>%
        as.data.frame()
      subtype_names <- x["subtype"] %>%
        stringr::str_split("\\|") %>%
        unlist()
      magrittr::set_colnames(subname_row, subtype_names)
    })
  
  # Combine the accession and subtype information into a single data frame
  all_metadata_df <-  data.frame(accession = all_accession,
                                 suppressMessages(dplyr::bind_rows(all_subname_list)))
  
  return(all_metadata_df)
}
extract_organism <- function(esummary_records, collapse_sequence_versions = TRUE) {
  # Check if input is a non-empty list of esummary records
  if (!is.list(esummary_records) || length(esummary_records) == 0) {
    stop("ERROR: Input must be a non-empty list of esummary records. Please use 'check_accessions()' to obtain the correct input.")
  }
  
  # Extract organism information from esummary records
  all_accession <- as.data.frame(t(rentrez::extract_from_esummary(esummary_records, elements = c("accessionversion", "organism")))) %>%
    dplyr::mutate(ncbi_accession = stringr::str_remove(accessionversion, stringr::regex("\\.[[:digit:]]"))) %>%
    tidyr::unnest(cols = c(accessionversion, organism)) %>%
    dplyr::rename(ncbi_accession_version = accessionversion) %>%
    dplyr::select(ncbi_accession, ncbi_accession_version, organism)
  
  # Collapse sequence versions if requested
  if (collapse_sequence_versions) {
    all_accession <- all_accession %>%
      dplyr::distinct(ncbi_accession, organism, .keep_all = TRUE)
  }
  
  # Check if number of unique rows matches the length of esummary records
  if (nrow(all_accession) != length(esummary_records)) {
    message("WARNING: The number of unique rows in input does not match the length of output. This is likely due to collapse_sequence_versions and is safe to proceed.")
  }
  
  # Return data frame with organism information
  return(all_accession)
}

# run these functions
subfeatures_df <- extract_subfeatures(summary_info)
organism_df <- extract_organism(summary_info)

# it would be good to keep organism information so lets join organism and subfeatures to create input df
input_df <- subfeatures_df[, c("host", "accession")] %>%
  dplyr::left_join(organism_df, by = c("accession" = "ncbi_accession_version")) %>%
  dplyr::select(-ncbi_accession)

# now lets check which of these host names are valid
seq_result <- vector("list", length(input_df[, "host"]))
query_index <- split(seq(1,length(input_df[, "host"])), ceiling(seq_along(seq(1,length(input_df[, "host"])))/50))
for (i in 1:length(query_index)) {
  seq_result[[i]]<-taxize::get_uid_(
    sci_com = stringr::str_trim(unique(input_df$host[unlist(query_index[i])])),
    ask = T,
    key = entrez_api_key
  )
  Sys.sleep(0.5)
}

uid_raw <- seq_result %>%
    purrr::flatten()

# Lets see which are none valid
message(
    cat(paste0('"', paste(names(uid_raw[lapply(uid_raw, length) == 0]), collapse="\", \""), '"')),
    " was not found in NCBI")

# store these in a df
info_df <-
  input_df %>% dplyr::filter(input_df[, "host"] %in% names(uid_raw[lapply(uid_raw, length) == 0])) %>%
  dplyr::left_join(organism_df, by = c("accession" = "ncbi_accession_version", "organism")) %>%
  dplyr::select(accession, host, organism)
     
# I then take this table that contains the tips without a valid host, for those 
# with a host but potentially a synonym or incorrectly entered I correct this using NCBI
# taxonomy. For those with no host information I search through associated papers 
# in an attempt to identify the host. Lastly I take the valid hosts from subfeatures_df
# to create a table with the columns "host" = virus host, "host_id" = ncbi_taxid for virus host, 
# accession = virus_accession, tree_tip_label = virus organism
# and manually add in hosts for the novel sequences. 

# This information is present in Supplementary table 3 in its final form
# we will load it in here. 

# The last step is to get the lineage information for each host

supt3 <- read_xlsx("Supplementary_table_3.xlsx", skip = 1) %>% 
  select(host,id,accession,tree_tip_label)

input_df <- supt3[, c("id", "host")]
input_df <- setNames(input_df, c("id", "host"))

if (!is.null("accession") && "accession" %in% colnames(checked_accessions)) {
  input_df <- cbind(input_df, supt3[, "accession"])
  colnames(input_df)[3] <- "accession"
}

output <- input_df %>%
  dplyr::bind_cols(
    input_df$id %>%
      taxize::classification(db = "ncbi", apikey = entrez_api_key) %>%
      lapply(function(x) {
        x %>%
          data.frame() %>%
          dplyr::distinct() %>%
          dplyr::filter(
            rank %in% c("class", "order", "family", "genus", "species", "kingdom", "phylum")
          ) %>%
          reshape2::melt(id = "rank") %>%
          tidyr::unite(col = rank_var, c("rank", "variable")) %>%
          tidyr::spread(key = rank_var, value = value)
      }) %>%
      dplyr::bind_rows() %>%
      dplyr::select(
        kingdom_name,
        kingdom_id,
        phylum_name,
        phylum_id,
        class_name,
        class_id,
        order_name,
        order_id,
        family_name,
        family_id,
        genus_name,
        genus_id,
        species_name,
        species_id
      )
  )
