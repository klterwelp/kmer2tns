## code to prepare `testCounts` dataset

# Read test_sequence fasta
test_fasta <- fs::path_package("extdata", "test_sequence.fasta", package = "kmer2tns")
# Count k_mers of Ecoli fasta
testCounts <- count_kmers_positions(test_fasta)

# Add ecoliCounts to data
usethis::use_data(testCounts, overwrite = TRUE, compress = "xz")

