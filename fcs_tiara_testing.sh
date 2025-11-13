# ============================================================
# Tiara–FCS-GX Comparison and BlobTools Assignment Generator
# ============================================================

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)

# ------------------------------
# Step 1: Define input paths
# ------------------------------
# Example: adjust these as needed
fcs_path   <- "path/to/fcs_output.txt"
tiara_path <- "path/to/tiara_output.tsv"
sample_id  <- "Sample_X"
output_dir <- "results"

dir.create(output_dir, showWarnings = FALSE)

# ------------------------------
# Step 2: Read FCS-GX results
# ------------------------------
# The FCS-GX format is space/tab separated with fields like:
# seq_id, seq_len, ..., species_name, taxid, lineage
# We'll keep seq_id, seq_len, species_name, lineage

fcs <- read.table(fcs_path, sep = "\t", quote = "", fill = TRUE, comment.char = "",
                  stringsAsFactors = FALSE)

# Extract essential columns using regex-safe indexing
# Typical useful columns:
# V1 = seq_id, V2 = seq_len, V9 = species_name, V10 = taxid, V11 = lineage
df_fcs <- fcs %>%
  transmute(
    seq_id = V1,
    seq_len = as.numeric(V2),
    species_fcs = ifelse(is.na(V9) | V9 == "", "Unknown", V9),
    lineage = ifelse(is.na(V11) | V11 == "", "Unknown", V11),
    domain_fcs = case_when(
      grepl("prok", lineage, ignore.case = TRUE) | grepl("bact", lineage, ignore.case = TRUE) ~ "bacteria",
      grepl("arch", lineage, ignore.case = TRUE) ~ "archaea",
      grepl("fung", lineage, ignore.case = TRUE) |
        grepl("plnt", lineage, ignore.case = TRUE) |
        grepl("anml", lineage, ignore.case = TRUE) ~ "eukarya",
      TRUE ~ "Unknown"
    )
  )

# ------------------------------
# Step 3: Read and process Tiara results
# ------------------------------
tiara <- read.delim(tiara_path, stringsAsFactors = FALSE, check.names = TRUE)

df_tiara <- tiara %>%
  filter(!grepl("organelle|mitochondria|plastid", class_fst_stage, ignore.case = TRUE)) %>%
  mutate(
    domain_tiara = ifelse(class_fst_stage == "" | class_fst_stage == "unknown", "Unknown", class_fst_stage)
  ) %>%
  select(seq_id = sequence_id, domain_tiara)

# ------------------------------
# Step 4: Merge and compare assignments
# ------------------------------
df_compare <- df_tiara %>%
  inner_join(df_fcs, by = "seq_id") %>%
  mutate(
    match = ifelse(domain_tiara == domain_fcs, "match", "mismatch")
  )

# ------------------------------
# Step 5: Generate final assignment table (for BlobTools)
# ------------------------------
df_final <- df_compare %>%
  mutate(
    tax_name = ifelse(
      match == "match",
      species_fcs,                                    # use FCS species if domains match
      paste(domain_tiara, domain_fcs, sep = "/")      # combine domains if mismatch
    ),
    tax_rank = ifelse(match == "match", "species", "domain_pair"),
    source = ifelse(match == "match", "fcs-gx", "combined")
  ) %>%
  select(seq_id, tax_name, tax_rank, source)

final_path <- file.path(output_dir, paste0(sample_id, "_BlobTools_Assignments.tsv"))
write.table(df_final, final_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste("\n✅ Final BlobTools assignment table written to:", final_path, "\n"))

# ------------------------------
# Step 6: Summaries by sequence length
# ------------------------------
summary_bp <- df_compare %>%
  group_by(match) %>%
  summarise(total_bp = sum(seq_len, na.rm = TRUE)) %>%
  mutate(percent_bp = round(100 * total_bp / sum(total_bp), 2))

cat("\nSummary of matches by total sequence length (organelles excluded):\n")
print(summary_bp)

# ------------------------------
# Step 7: Breakdown by domain pair
# ------------------------------
domain_bp <- df_compare %>%
  mutate(pair = paste(domain_tiara, "→", domain_fcs)) %>%
  group_by(pair) %>%
  summarise(total_bp = sum(seq_len, na.rm = TRUE)) %>%
  mutate(percent_bp = round(100 * total_bp / sum(total_bp), 2))

cat("\nBreakdown by domain pair (as % of total bp):\n")
print(domain_bp)

# ------------------------------
# Step 8: Plots
# ------------------------------
bar_plot <- ggplot(domain_bp, aes(x = reorder(pair, percent_bp), y = percent_bp)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = paste("Tiara vs FCS-GX:", sample_id, "\n% Total Sequence Length per Domain ≥1 kb"),
    x = "Tiara → FCS-GX",
    y = "% of total sequence length"
  ) +
  theme_minimal(base_size = 13)

bar_path <- file.path(output_dir, paste0(sample_id, "_Tiara_FCSGX_domain_barplot.png"))
ggsave(bar_path, bar_plot, width = 7, height = 5, dpi = 300)
cat(paste("Bar plot saved to:", bar_path, "\n"))

# Stacked bar plot: match vs mismatch
stacked_plot <- ggplot(summary_bp, aes(x = "Sample", y = percent_bp, fill = match)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  scale_fill_manual(values = c("match" = "#2ca02c", "mismatch" = "#d62728")) +
  labs(
    title = paste("Tiara vs FCS-GX Classification Agreement (≥1 kb)\n", sample_id),
    x = "",
    y = "% of total sequence length",
    fill = "Comparison"
  ) +
  geom_text(aes(label = paste0(percent_bp, "%")),
            position = position_stack(vjust = 0.5),
            color = "white", size = 5) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

stacked_path <- file.path(output_dir, paste0(sample_id, "_FCSGX_Tiara_MatchMismatch_StackedBar.png"))
ggsave(stacked_path, stacked_plot, width = 5, height = 4, dpi = 300)
cat(paste("Stacked bar plot saved to:", stacked_path, "\n"))

# ============================================================
# End of script
# ============================================================
