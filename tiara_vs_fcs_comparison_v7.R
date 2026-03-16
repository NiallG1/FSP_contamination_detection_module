# =====================================================================
# Tiara vs FCS-GX comparison and visualization script version 7
# =====================================================================

#install package
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)

# ----------------------------------------------
# Step 0: Define file paths and output directory
# ----------------------------------------------
fcs_path <- "" 
tiara_path <- ""
output_dir <- ""

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# ---------------------------------------------
# Step 1: Extract sample ID from input filename
# ---------------------------------------------
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
    div1 = div,           # domain-level classification (e.g. fung:basidiomycetes)
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


# --------------------------------------------
# Step 3B: Process FCS-GX results for chimeras
# --------------------------------------------
library(dplyr)
library(stringr)

df_fcs_collapsed <- df_fcs %>%
  mutate(contig_base = str_remove(seq_id, "~.*")) %>%
  group_by(contig_base) %>%
  summarise(
    seq_len = sum(seq_len),

    n_domains = n_distinct(domain_fcs[domain_fcs != "Unknown"]),
    n_species = n_distinct(species_fcs[species_fcs != "Unknown"]),

    species_list = paste(
      unique(species_fcs[species_fcs != "Unknown"]),
      collapse = "; "
    ),

    species_fcs = case_when(
      n_domains > 1 ~ "possible cross-domain chimera",
      n_species == 1 ~ first(species_fcs[species_fcs != "Unknown"]),
      n_species > 1 ~ "possible chimera",
      TRUE ~ "Unknown"
    ),

    chimera_species = case_when(
      species_fcs %in% c("possible chimera", "possible cross-domain chimera") ~ species_list,
      TRUE ~ NA_character_
    ),

    domain_fcs = case_when(
      n_domains == 1 ~ first(domain_fcs[domain_fcs != "Unknown"]),
      n_domains > 1 ~ "mixed",
      TRUE ~ "Unknown"
    ),

    div1 = first(div1),
    .groups = "drop"
  ) %>%
  rename(seq_id = contig_base)
df_fcs_collapsed <- df_fcs_collapsed %>%
  arrange(as.numeric(stringr::str_extract(seq_id, "(?<=NODE_)\\d+")))

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


library(dplyr)
library(tidyr)

# ------------------------------
# Step 5: Merge and compare new
# ------------------------------

df_compare <- df_fcs_collapsed %>%
  left_join(
    df_tiara %>% select(seq_id, domain_tiara),
    by = "seq_id"
  ) %>%
  replace_na(list(domain_tiara = "Unknown")) %>%
  mutate(
    match = case_when(
      species_fcs %in% c("possible chimera", "possible cross-domain chimera") ~ "chimera",
      domain_tiara == domain_fcs ~ "match",
      TRUE ~ "mismatch"
    )
  )

# ----------------------------------------
# step 6: Save the merged comparison table
# ----------------------------------------
out_file <- file.path(output_dir, paste0(sample_id, "tiara_vs_fcs_compare.tsv"))
write_tsv(df_compare, out_file)



# ------------------------------
# step 7: create  blobtags
# ------------------------------

df_compare <- df_compare %>%
  mutate(
    # Standardize unknowns first
    species_fcs = ifelse(is.na(species_fcs) | species_fcs == "Unknown", "unknown", species_fcs),

    # Replace spaces with underscores
    species_fcs = gsub(" ", "_", species_fcs),

    # Create blob_tag depending on match status
    blob_tag = case_when(
      seq_len < 1000 ~ "Unknown",                # override small contigs
      match == "match" ~ species_fcs,            # agreement → species
      match == "chimera" ~ "possible_chimera",   # new chimera case
      TRUE ~ paste0(domain_tiara, "_", domain_fcs)  # mismatch
    )
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
