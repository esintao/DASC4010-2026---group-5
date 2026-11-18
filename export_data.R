if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("curatedMetagenomicData", quietly = TRUE)) BiocManager::install("curatedMetagenomicData")
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) BiocManager::install("SummarizedExperiment")

library(curatedMetagenomicData)
library(SummarizedExperiment)

# Load YachidaS_2019
se_list <- curatedMetagenomicData("YachidaS_2019.relative_abundance", dryrun = FALSE)
se <- se_list[[1]]

# 1. METADATA TABLE
metadata_df <- as.data.frame(colData(se))
metadata_df <- cbind(sample_id = rownames(metadata_df), metadata_df)

# 2. OTU TABLE (Samples x Features)
raw_assay <- assay(se, "relative_abundance")
tax_strings <- rownames(raw_assay) # Original taxonomic strings

# Give features clean IDs (e.g., feat_1, feat_2...) to keep table lightweight
feature_ids <- paste0("feat_", seq_len(nrow(raw_assay)))
rownames(raw_assay) <- feature_ids

otu_df <- as.data.frame(t(raw_assay))
otu_df <- cbind(sample_id = rownames(otu_df), otu_df)

# 3. TAXONOMY TABLE (Features x Ranks)
parse_taxonomy <- function(tax_str) {
  parts <- strsplit(tax_str, "\\|")[[1]]
  ranks <- c(Kingdom = NA, Phylum = NA, Class = NA, Order = NA, Family = NA, Genus = NA, Species = NA)
  for (p in parts) {
    if (grepl("^k__", p)) ranks["Kingdom"] <- sub("^k__", "", p)
    if (grepl("^p__", p)) ranks["Phylum"] <- sub("^p__", "", p)
    if (grepl("^c__", p)) ranks["Class"] <- sub("^c__", "", p)
    if (grepl("^o__", p)) ranks["Order"] <- sub("^o__", "", p)
    if (grepl("^f__", p)) ranks["Family"] <- sub("^f__", "", p)
    if (grepl("^g__", p)) ranks["Genus"] <- sub("^g__", "", p)
    if (grepl("^s__", p)) ranks["Species"] <- sub("^s__", "", p)
  }
  return(as.list(ranks))
}

tax_list <- lapply(tax_strings, parse_taxonomy)
tax_df <- do.call(rbind, lapply(tax_list, as.data.frame))
rownames(tax_df) <- feature_ids
tax_df <- cbind(feature_id = rownames(tax_df), tax_df)

# Export all three as tab-delimited text files
write.table(metadata_df, "metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(otu_df, "otu_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tax_df, "tax_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cat("Successfully exported metadata.txt, otu_table.txt, and tax_table.txt!\n")