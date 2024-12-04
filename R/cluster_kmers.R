# Function to set top_kmer to the kmer with max_count and lowest in alphabet

#' Identify top kmer in a kmer_count table
#'
#' @param table kmer_count table
#'
#' @return kmer with highest count and its count
#' @export
#'
#' @examples
#' identify_top_kmer(ecoliCounts)
#' @importFrom magrittr %>%
#' @import data.table
identify_top_kmer <- function(table) {
  # set top_kmer to kmer with the max_count
  max_count <- table[["count"]] %>% max()
  top_kmers <- table[table[["count"]] == max_count]
  # requires with to force "count" colname to be explicitly used (long debug)
  if (length(top_kmers$kmer) > 1) { # if there's more than one kmer with max count
    # set top_kmer to alphabetically smallest
    top_kmer_name <- min(top_kmers$kmer) # pick kmer earliest in alphabet
    top_kmer <- top_kmers[top_kmers[["kmer"]] == top_kmer_name]
  } else {
    # set top_kmer to top_kmers
    top_kmer <- top_kmers
    top_kmer_name <- top_kmer$kmer
  }
  return(top_kmer)
}


#' Identify core kmers from a kmer_count table
#'
#' @param table kmer_count table
#' @param max_diff max allowed differences between kmers
#'
#' @return kmer_count table of just the most repetitive kmer and its closely related sequences
#' @export
#'
#' @examples
#' new_identify_core_kmers(ecoliCounts, 3)
new_identify_core_kmers <- function(table, max_diff) {
  # identify top_kmer and its count
  top_kmer <- identify_top_kmer(table)
  top_kmer_name <- top_kmer$kmer
  # add top_kmer to core_table
  core_table <- top_kmer
  # remove top_kmer from kmer_counts
  table <-table[table[["kmer"]] != top_kmer$kmer]
  # add other kmers to core_table if they have fewer diffs than max_diff
  table_ham <- table %>%
    dplyr::mutate(hamming_f = stringdist::stringdist(top_kmer_name, table$kmer, method = "hamming"))
  # add kmers that are less than max_diff into core_table
  filt_table <- table_ham[table_ham[["hamming_f"]] <= max_diff]
  filt_table <- filt_table %>% dplyr::select(-hamming_f)
  core_table <- data.table::rbindlist(list(core_table, filt_table))
  # compare reverse_complement of top_kmer w/ sequences that did not match forward
  rev_top_kmer <- Biostrings::reverseComplement(Biostrings::DNAString(top_kmer_name)) %>% as.character()
  filt_table <- table_ham[table_ham[["hamming_f"]] > max_diff] # only sequences > max_diff
  filt_table <- filt_table %>%
    dplyr::mutate(hamming_r = stringdist::stringdist(rev_top_kmer, filt_table$kmer, method = "hamming"))
  filt_table <- filt_table[filt_table[["hamming_r"]] <= max_diff] # find any sequences less than max_diff for reverse_complement
  filt_table <- filt_table %>% dplyr::select(-hamming_f, -hamming_r) # remove extra cols
  core_table <- data.table::rbindlist(list(core_table, filt_table)) # combine tables
  return(core_table)
}

#' Expand core kmers
#'
#' @param table kmer_counts used to generate core_table list
#' @param core_list list of core_tables generated
#' @param cutoff_list cutoffs generated automatically from core_cutoff
#'
#' @return list of expanded p_clouds
#' @export
#'
new_expand_core_kmers <- function(table, core_list, cutoff_list) {
  # filter table by core_cutoff (only need to expand for non-core kmers)
  table <- table[table[["count"]] < cutoff_list$core_cutoff]
  # start empty expanded_p_cloud list and p_cloud_name
  expanded_p_clouds <- list()
  p_cloud_name <- 0
  for (core_table in core_list) {
    # add 1 to p_cloud_name
    p_cloud_name <- p_cloud_name + 1
    # use top_kmer_count to identify max_diff allowed for expansion
    top_kmer_count <- identify_top_kmer(core_table)$count
    if (top_kmer_count > cutoff_list$tertiary_cutoff) {
      max_diff <- 3
    } else if (top_kmer_count > cutoff_list$secondary_cutoff) {
      max_diff <- 2
    } else if (top_kmer_count > cutoff_list$primary_cutoff) {
      max_diff <- 1
    } else {
      max_diff <- 0
    }
    # add core_table to expanded_core_table
    expanded_core_table <- core_table
    # only add expanded kmers if top_kmer_count is greater than primary cutoff
    if (max_diff > 0) {
      # add all core kmers to expanded core_table
      # remove any kmers that are already in the expanded_p_clouds
      if (length(expanded_p_clouds) > 0) {
        p_cloud_kmers <- unlist(lapply(expanded_p_clouds, function(dt) dt$kmer)) # find all kmers
        table <- table[!table[["kmer"]] %in% p_cloud_kmers] # remove any kmers already in table
      }
      # compare to kmer_count table
      hamming_matrix <- stringdist::stringdistmatrix(table$kmer, core_table$kmer, useNames = TRUE, method = "hamming")
      # convert to data.table for faster querying
      hamming_dt <- as.data.frame(hamming_matrix, row.names = rownames(hamming_matrix)) %>%
        data.table::as.data.table(keep.rownames=TRUE)
      # filter out any values that are less than max_diff and find the kmer names (colnames)
      expanded_kmer_names <-  hamming_dt[rowSums(hamming_dt <= max_diff) > 0]$rn
      # if any of kmers have a hamming distance less than or equal to max_diff, add to expanded kmer name
      expanded_kmers <- table[table[["kmer"]] %in% expanded_kmer_names]
      log_message("Adding: ", length(expanded_kmer_names), " to P-cloud: ", p_cloud_name)
      # add expanded_kmers to expanded_core_table
      expanded_core_table <- data.table::rbindlist(list(expanded_kmers, expanded_core_table))
    }

    expanded_p_clouds[[p_cloud_name]] <- expanded_core_table
  }
  return(expanded_p_clouds)
}


#' Cluster kmers into similar repetitive sequences
#'
#' @param core_cutoff any integer value greater than 5
#' @param kmer_counts counts table from count_kmers()
#'
#' @return table of kmer clusters (P-cloud object)
#' @export
#'
#' @examples
#' cluster_kmers(10, ecoliCounts)
cluster_kmers <- function(core_cutoff, kmer_counts) {
  # Check if user given variables are in the correct format
  log_message("Checking if core_cutoff is numeric...")

  if (!is.numeric(core_cutoff)) {
    stop("core_cutoff must be a numeric value.")
  }

  if (core_cutoff < 5) {
    stop("core_cutoff must be greater than 5")
  }

  log_message("Checking if kmer_counts is in the kmer_count table format...")

  if (!validate_kmer_count_table(kmer_counts)) {
    stop("kmer_count must be in kmer_count data.table format")
  }

  log_message("Checking if core_cutoff value makes sense...")

  if (core_cutoff >= max(kmer_counts[["count"]])) { # greater than maximum count in given table
    stop("core_cutoff must be less than max count in kmer_counts")
  }

  log_message(max(kmer_counts[["count"]]))

  # SET variables based on core_cutoff:
  log_message("Setting variables based on core_cutoff...")

  # setting lower_cutoff
  if (core_cutoff <= 10) {
    lower_cutoff <- 2 # set lower_cutoff to 2
  } else {
    lower_cutoff <- round(core_cutoff/10) # set lower_cutoff to 1/10 core_cutoff
  }

  # setting primary-tertiary cutoffs
  primary_cutoff <- 2*core_cutoff
  secondary_cutoff <- 20*core_cutoff
  tertiary_cutoff <- 200*core_cutoff

  cutoff_list <- list("lower_cutoff" = lower_cutoff,
                      "core_cutoff" = core_cutoff,
                      "primary_cutoff" = primary_cutoff,
                      "secondary_cutoff" = secondary_cutoff,
                      "tertiary_cutoff" = tertiary_cutoff)

  log_message("Cutoff variables assigned:", "\n",
              "lower_cutoff: ", lower_cutoff, "\n",
              "primary_cutoff: ", primary_cutoff, "\n",
              "secondary_cutoff: ", secondary_cutoff, "\n",
              "tertiary_cutoff: ", tertiary_cutoff, "\n")


  log_message("Identifying core kmers...")
  p_clouds_list <- list()
  p_cloud_name <- 0
  log_message("Filtering kmer_counts by core_cutoff...")
  core_pool <- kmer_counts[kmer_counts[["count"]] > core_cutoff]
  log_message("Filtered kmer_counts by core_cutoff.")

  while (nrow(core_pool) > 0) { #while there's still kmers in the core_pool
    log_message("There are ", nrow(core_pool), "rows")
    core_table <- new_identify_core_kmers(core_pool, 3) # identify core_kmers
    max_count <- max(core_table[["count"]])
    p_cloud_name <- p_cloud_name + 1
    log_message("P_cloud: ", p_cloud_name)
    p_clouds_list[[p_cloud_name]] <- core_table
    # add core_table to list with top_core_kmer as its name
    core_pool <- core_pool[!core_pool[["kmer"]] %in% core_table[["kmer"]]]
  }
  log_message("Finished identifying core kmers.")
  log_message("Expanding p-clouds...")

  expanded_p_clouds_list <- new_expand_core_kmers(table = kmer_counts,
                        core_list = p_clouds_list,
                        cutoff_list = cutoff_list)

  log_message("Done expanding p-clouds.")
  log_message("Converting p-cloud list to data.table")
  expanded_p_clouds_tbl <- rbindlist(expanded_p_clouds_list, idcol = "p_cloud")

  return(expanded_p_clouds_tbl)

}








