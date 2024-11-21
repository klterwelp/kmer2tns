## code to prepare `ecoliCounts` dataset

# Read Ecoli fasta
ecoli_fasta <- fs::path_package("extdata", "ecoliJE86-ST05.fasta.gz", package = "kmer2tns")
# Count k_mers of Ecoli fasta
ecoliCounts <- count_kmers_positions(ecoli_fasta)

# Add ecoliCounts to data
usethis::use_data(ecoliCounts, overwrite = TRUE, compress = "xz")

