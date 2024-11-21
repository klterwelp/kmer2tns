# tests for cluster_kmers
# Helper function to create test data
create_test_kmer_table <- function(kmers, counts, positions) {
  data.table(kmer = kmers, count = counts, positions = positions)
}

test_that("identify_top_kmer identifies the correct top k-mer", {
  # Test in case of singal max k-mer
  dt1 <- create_test_kmer_table(c("AAAA", "CCCC", "GGGG"),
    c(10, 5, 3),
    list(c(1:9), c(10:14), c(800:802)))
  expect_equal(identify_top_kmer(dt1)$kmer, "AAAA")
  expect_equal(identify_top_kmer(dt1)$count, 10)

  # Test multiple kmers with same max count (should pick alphabetically first)
  dt2 <- create_test_kmer_table(
    c("CCCC", "AAAA", "GGGG"),
    c(10, 10, 3),
    list(c(1:9), c(10:19), c(800:802))
  )
  expect_equal(identify_top_kmer(dt2)$kmer, "AAAA")
  expect_equal(identify_top_kmer(dt2)$count, 10)

})


test_that("new_identify_core_kmers identifies related sequences correctly", {
  # Test exact matches and hamming distance cases
  dt <- create_test_kmer_table(
    c("AAAA", "AAAT", "AACA", "GGGG", "TTTC"),
    c(10, 8, 7, 2, 1),
    list(c(1:9), c(10:17), c(40:46), c(80:81), c(800))
  )

  # Test with max_diff = 1
  result1 <- new_identify_core_kmers(dt, 1)
  expect_true("AAAA" %in% result1$kmer)
  expect_true("AAAT" %in% result1$kmer)
  expect_false("GGGG" %in% result1$kmer)

  # Test with max_diff = 2
  result2 <- new_identify_core_kmers(dt, 2)
  expect_true("AACA" %in% result2$kmer)
  expect_false("GGGG" %in% result2$kmer)

  # Test reverse complement inclusion
  dt_rev <- create_test_kmer_table(
    c("AAAA", "TTTT", "GGGG"),
    c(10, 8, 2),
    list(c(1:9), c(10:17), c(800:801))
  )
  result_rev <- new_identify_core_kmers(dt_rev, 0)
  expect_true("TTTT" %in% result_rev$kmer)
})

test_that("new_expand_core_kmers expands clusters correctly", {
  # Create test data
  main_table <- create_test_kmer_table(
    c("AAAA", "AAAT", "AACC", "ATAA", "GGGG", "GGGC", "GGCC"),
    c(250, 8, 7, 6, 220, 6, 5),# Core kmers above secondary_cutoff (200)
    list(c(1:249), c(1:7), c(1:6), c(1:5), c(200:419), c(1000:1005), c(2000:2004))
  )

  core_list <- list(
    create_test_kmer_table(c("AAAA"), c(250), list(c(1:249))),  # Above secondary_cutoff
    create_test_kmer_table(c("GGGG"), c(220), list(c(1:219)))   # Above secondary_cutoff
  )

  cutoff_list <- list(
    lower_cutoff = 2,
    core_cutoff = 10,
    primary_cutoff = 20,
    secondary_cutoff = 200,
    tertiary_cutoff = 2000
  )

  result <- new_expand_core_kmers(main_table, core_list, cutoff_list)

  # Should now include kmers with hamming distance 2
  expect_true("AACC" %in% result[[1]]$kmer)  # Distance 2 from AAAA
  expect_true("GGCC" %in% result[[2]]$kmer)  # Distance 2 from GGGG

  # Should still include distance 1 kmers
  expect_true("AAAT" %in% result[[1]]$kmer)
  expect_true("GGGC" %in% result[[2]]$kmer)
})

test_that("cluster_kmers handles edge cases correctly", {
  # Test invalid inputs
  expect_error(cluster_kmers("10", create_test_kmer_table(c("AAAA"), c(10), list(c(10)))))
  expect_error(cluster_kmers(4, create_test_kmer_table(c("AAAA"), c(10), list(c(10)))))
  expect_error(cluster_kmers(10, create_test_kmer_table(c("AAAA"), c(8), list(c(1)))))

  # Test minimal valid input
  min_table <- create_test_kmer_table(
    c("AAAA", "AAAT", "GGGG"),
    c(10, 8, 6),
    list(c(1:9), c(10:17), c(800:805))
  )
  result <- cluster_kmers(5, min_table)
  expect_s3_class(result, "data.table")

  # Test cutoff scaling
  big_table <- create_test_kmer_table(
    c("AAAA", "AAAT", "AACC", "GGGG"),
    c(1000, 800, 700, 600),
    list(c(1:999), c(1000:1799), c(2000:2599), c(9000:9599))
  )
  result_big <- cluster_kmers(100, big_table)
  expect_s3_class(result_big, "data.table")

  # Verify each cluster contains related sequences
  expect_true(result_big[, {
    expect_true(all(count >= 10))
    if (.N > 1) {
      top_kmer <- identify_top_kmer(.SD)$kmer
      other_kmers <- kmer[kmer != top_kmer]
      ham_dists <- stringdist::stringdist(top_kmer, other_kmers, method = "hamming")
      expect_true(all(ham_dists <= 3))
    }
  }, by = p_cloud]$V1)

})
