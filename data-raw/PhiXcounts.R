## code to prepare `PhiXcounts` dataset

# Count k_mers of PhiX fasta
phiX_fasta <- fs::path_package("extdata", "phiX.fasta.gz", package = "kmer2tns")

count_matrix <- Biostrings::readDNAStringSet(phiX_fasta) %>%
  Biostrings::oligonucleotideFrequency(width = 8)

# Convert to data.table and only include count > 1
PhiXcounts <- data.table(kmer = colnames(count_matrix), count = as.numeric(count_matrix[1, ]))[count > 1]

# Add PhiXcounts to data
usethis::use_data(PhiXcounts, overwrite = TRUE)

utils::write.csv(PhiXcounts, "inst/extdata/PhiXcounts.csv")
