# Check that fs is installed to run unit tests
if (requireNamespace("fs", quietly = TRUE)) {
  library(fs)
} else {
  stop("Package 'fs' is required for these tests.")
}

# Test for file not found
test_that("count_kmers handles file not found error", {
  expect_error(count_kmers(seq_path = "non_existent_file.fasta", k = 7), "File not found")
})

# Test for file not Fasta
test_that("count_kmers handles file not fasta error", {
  csv <- fs::path_package("extdata", "phiXcounts.csv", package = "kmer2tns")
  expect_error(count_kmers(seq_path = csv, k = 7))
})

# Test for invalid k value
dummy_fasta <- fs::path_package("extdata", "phiX.fasta.gz", package = "kmer2tns")
test_that("count_kmers handles invalid k value", {
  expect_error(count_kmers(seq_path = dummy_fasta , k = 0), "Selected k is too small")
  expect_error(count_kmers(seq_path = dummy_fasta, k = 501), "Selected k is too large")
})

# Test for default k value
test_that("count_kmers uses default k value if not provided", {
  result <- count_kmers(seq_path = dummy_fasta)
  default_k <- round(log(5386, base = 4)) + 1  # 5386 = length of PhiX sequence
  expect_equal(nchar(result$kmer[1]), default_k)
})

# Test for correct k-mer counting
test_that("count_kmers correctly counts kmers", {
    result <- count_kmers(dummy_fasta, 8)
    expect_equal(result, PhiXcounts)
  })


# Test that output is a data.table
test_that("count_kmers correctly counts kmers", {
  result <- count_kmers(dummy_fasta, 8)
  expect_s3_class(result, "data.table")
})

