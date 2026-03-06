# =====================================================================
# Tiara vs FCS-GX comparison and visualization script
# =====================================================================

#install package
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)

# ------------------------------
# Step 0: Define file paths and output directory
# ------------------------------
fcs_path <- "/home/niall/FSP/contamination_detection/tiara_vs_fcs/debug/EGP017_25_Com_1_FCS_GX.tsv" 
tiara_path <- "/home/niall/FSP/contamination_detection/tiara_vs_fcs/debug/tiara_Com_1.txt"
output_dir <- "/home/niall/FSP/contamination_detection/tiara_vs_fcs/debug"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# ------------------------------
# Step 1: Extract sample ID from input filename
# ------------------------------
# This assumes the sample ID is the first underscore-separated block that looks like EGP###_##_###
# (You can adjust this regex if needed.)
sample_id <- str_extract(
  basename(fcs_path),
  "EGP[0-9]+_[0-9]+_B[0-9]+_[0-9]+"
)

if (is.na(sample_id)) {
  sample_id <- tools::file_path_sans_ext(basename(fcs_path))
}

cat(paste("Processing sample:", sample_id, "\n"))


# ------------------------------
# Step 2: Read files
# ------------------------------
fcs <- read.delim(fcs_path, stringsAsFactors = FALSE, check.names = TRUE)
tiara <- read.delim(tiara_path, stringsAsFactors = FALSE, check.names = TRUE)


# ------------------------------
# Step 3: Process FCS-GX results
# ------------------------------
df_fcs <- fcs %>%
  select(
    seq_id = seq.id,        # contig ID
    seq_len = seq.len,      # sequence length
    div1 = div.1,           # domain-level classification (e.g. fung:basidiomycetes)
    species_fcs = tax.name.1 # species-level assignment (e.g. Boletus reticuloceps)
  ) %>%
  mutate(
    seq_len = as.numeric(seq_len),
    div1 = ifelse(
      is.na(div1) | div1 == "" | div1 == "unassigned" | div1 == "unknown",
      "Unknown",
      div1
    ),
    species_fcs = ifelse(
      is.na(species_fcs) | species_fcs == "",
      "Unknown",
      species_fcs
    ),
    domain_fcs = case_when(
      grepl("^prok", div1, ignore.case = TRUE) | grepl("^bact", div1, ignore.case = TRUE) ~ "bacteria",
      grepl("^arch", div1, ignore.case = TRUE) ~ "archaea",
      grepl("^fung", div1, ignore.case = TRUE) |
        grepl("^plnt", div1, ignore.case = TRUE) |
        grepl("^anml", div1, ignore.case = TRUE) ~ "eukarya",
      div1 == "Unknown" ~ "Unknown",
      TRUE ~ "Unknown"
    )
  )

# ------------------------------
# Step 4: Process Tiara results 
# ------------------------------
df_tiara <- tiara %>%
  mutate(
    domain_tiara = ifelse(
      class_fst_stage == "" | class_fst_stage == "unknown",
      "Unknown",
      class_fst_stage
    )
  ) %>%
  select(seq_id = sequence_id, domain_tiara)




# ------------------------------
# Step 5: Merge and compare
# ------------------------------
#due to an error, when joining, all contigs present in FCS-GX but absent in tiara were dropped from the final table.
#This is now ammended 



library(dplyr)
library(tidyr)

df_compare <- df_fcs %>%
  left_join(df_tiara %>% select(seq_id, domain_tiara), by = "seq_id") %>%
  replace_na(list(domain_tiara = "Unknown")) %>%
  mutate(match = ifelse(domain_tiara == domain_fcs, "match", "mismatch"))

# ------------------------------
# Save the merged comparison table
# ------------------------------
out_file <- file.path(output_dir, paste0(sample_id, "tiara_vs_fcs_compare.tsv"))
write_tsv(df_compare, out_file)

# ------------------------------
# Step 6: Summary by sequence length
# ------------------------------
summary_bp <- df_compare %>%
  group_by(match) %>%
  summarise(total_bp = sum(seq_len, na.rm = TRUE)) %>%
  mutate(percent_bp = round(100 * total_bp / sum(total_bp), 2))

cat("Summary of matches by total sequence length (organelles excluded):\n")



# ------------------------------
# Step 7: Breakdown by domain pair
# ------------------------------
domain_bp <- df_compare %>%
  mutate(pair = paste(domain_tiara, "→", domain_fcs)) %>%
  group_by(pair) %>%
  summarise(total_bp = sum(seq_len, na.rm = TRUE)) %>%
  mutate(percent_bp = round(100 * total_bp / sum(total_bp), 2))

cat("\nBreakdown by domain pair (as % of total bp):\n")

# ------------------------------
# create  blobtags
# ------------------------------

df_compare <- df_compare %>%
  mutate(
    # Standardize unknowns first
    species_fcs = ifelse(is.na(species_fcs) | species_fcs == "Unknown", "unknown", species_fcs),
    
    
    # Replace spaces with underscores
    species_fcs = gsub(" ", "_", species_fcs),
    
    # Create blob_tag depending on match status
    blob_tag = ifelse(
      match == "match",
      species_fcs,  # if match → use species name
      paste0(domain_tiara, "_", domain_fcs)  # if mismatch → combine domain_species
    ),
    
    # Override blob_tag for contigs <1kb
    blob_tag = ifelse(seq_len < 1000, "Unknown", blob_tag)
  )


# ------------------------------
# Step 9: Export comparison CSV
# ------------------------------
#now all contigs >1kbp are labelled as "unknown"
#this is as tiara will not test contigs below 1kbp and FCS-GX is not accurate below 1kbp.

blob_taxonomy <- df_compare %>%
  select(seq_id, blob_tag) %>%
  rename(taxonomy = blob_tag)





# Write it to the output directory
write_tsv(blob_taxonomy, file.path(output_dir, paste0(sample_id, "blobtools_taxonomy.tsv")))
