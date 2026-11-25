library(curatedMetagenomicData)
library(mia)
library(phyloseq)

tse_list <- curatedMetagenomicData("YachidaS_2019.relative_abundance", dryrun = FALSE, counts = FALSE)
tse <- tse_list[[1]]

ps_yachida <- convertToPhyloseq(tse, assay.type = "relative_abundance")

# 2. Verify the phyloseq object
ps_yachida


# 1. Export OTU / Abundance Table (taxa as rows, samples as columns)
otu_df <- as.data.frame(otu_table(ps_yachida))
write.table(otu_df, file = "otu_table.txt", sep = "\t", quote = FALSE, col.names = NA)

# 2. Export Sample Metadata Table
# Wrapping in data.frame() converts phyloseq's S4 sample_data class into a standard data frame
meta_df <- data.frame(sample_data(ps_yachida), check.names = FALSE)
write.table(meta_df, file = "sample_data.txt", sep = "\t", quote = FALSE, col.names = NA)

# 3. Export Taxonomy Table
tax_df <- as.data.frame(tax_table(ps_yachida))
write.table(tax_df, file = "tax_table.txt", sep = "\t", quote = FALSE, col.names = NA)