# Purpose: Function to cluster k-mer counts into P-clouds of repetitive k-mers
#' @importFrom magrittr %>%

# Establish variables based on user_defined core_cutoff
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

  if (core_cutoff >= max(kmer_counts$count)) { # greater than maximum count in given table
    stop("core_cutoff must be less than max count in kmer_counts")
  }

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

  log_message("Cutoff variables assigned:", "\n",
              "lower_cutoff: ", lower_cutoff, "\n",
              "primary_cutoff: ", primary_cutoff, "\n",
              "secondary_cutoff: ", secondary_cutoff, "\n",
              "tertiary_cutoff: ", tertiary_cutoff, "\n")

  # Function to compare sequences with reverse complement in mind
  compare_seqs_rev <- function(seq_1, seq_2, max_diff) {
    # convert character strings into vector of characters
    string2vec <- function(seq) {
      unlist(strsplit(seq, split = ""))
    }
    vec_1 <- string2vec(seq_1)
    vec_2 <- string2vec(seq_2)
    # check that character strings are the same length
    if (length(vec_1) != length(vec_2)) {
      stop("kmers are not the same length")
    }
    # compare_seqs to max_diff function
    compare_seqs <- function(vec_1, vec_2, max_diff) {
      diff <- 0
      for ( i in 1:length(vec_1)) { # for each character in vec_1
        if (vec_1[i] != vec_2[i]) { # if character isn't the same
          diff <- diff + 1 # add one to diff
          if (diff > max_diff) {
            break
          }
        }
      }
      return(diff) # return differences
    }

    # compare sequences including reverse complement of seq_2
    fw_diff <- compare_seqs(vec_1, vec_2, max_diff) # set forward difference
    if (fw_diff > max_diff) { # if fw_diff > max_diff, check reverse complement differences
      rev_seq_2 <- Biostrings::reverseComplement(Biostrings::DNAString(seq_2)) %>%
        as.character()
      rev_vec_2 <- string2vec(rev_seq_2) # turn to vector of characters
      rc_diff <- compare_seqs(vec_1, rev_vec_2, max_diff) # compare reverse complement
      diff <- rc_diff # set differences to reverse complement differences
    } else {
      diff <- fw_diff # set differences to forward complement differences
    }
    return(diff)

  }

  # Function to identify core_kmers
  identify_core_kmers <- function(kmer_counts, max_diff) {
    # create empty core_table
    core_table <- data.table::data.table(kmer = character(), count = numeric())
    # set top_kmer to kmer with the max_count
    max_count <- max(kmer_counts$count)
    top_kmers <- kmer_counts[count == max_count]
    # return table of top_kmers if there's multiple
    if (length(top_kmers$kmer) > 1) { # if there's more than one kmer with max count
      # set top_kmer to alphabetically smallest
      top_kmer_name <- min(top_kmers$kmer)
      top_kmer <- top_kmers[kmer == top_kmer_name]
    } else {
      # set top_kmer to top_kmers
      top_kmer <- top_kmers
      top_kmer_name <- top_kmer$kmer
    }
    # add top_kmer to core_table
    core_table <- rbind(core_table, top_kmer)
    # remove top_kmer from kmer_counts
    kmer_counts <-kmer_counts[kmer != top_kmer$kmer]
    # add other kmers to core_table if they have fewer diffs than max_diff
    for (kmer_name in kmer_counts$kmer) { # for each kmer in kmer_counts
      diff <- compare_seqs_rev(top_kmer$kmer, kmer_name, max_diff) # compare top_kmer vs kmer
      if(diff <= max_diff) { # if differences less than max_diff
        kmer_table <- kmer_counts[kmer == kmer_name] # acquire kmer counts info
        core_table <- rbind(core_table, kmer_table) # add to core table
      }
    }
    return(core_table)
  }

  # Identify core_kmers and return a p_clouds_core list
  # Create empty p_clouds_list
  log_message("Identifying core kmers...")
  p_clouds_list <- list()
  log_message("Filtering kmer_counts by core_cutoff...")
  core_pool <- kmer_counts[count > core_cutoff]

  while (nrow(core_pool) > 0) { #while there's still kmers in the core_pool
    core_table <- identify_core_kmers(kmer_counts = core_pool, max_diff = 3) # identify core_kmers
    max_count <- max(core_table$count)
    top_core_kmers <- core_table[count == max_count]$kmer # find name of top_core_kmer
    top_core_kmer_name <- min(top_core_kmers)
    p_clouds_list[[top_core_kmer_name]] <- core_table
    # add core_table to list with top_core_kmer as its name
    core_pool <- core_pool[!kmer %in% core_table$kmer]
  }
}


# Function for Comparing Sequences
# Initalize counter of differences
# for each character in the first sequence
# if letter at position i isn't equal in sequence one and two
# add one to differences
# if differences is greater than max_diff
# exit loop
# return differences

# Function for P cloud Core Identification
# top_core_kmer = k-mer with highest count in core_pool
# top_core_kmer = Sort k-mers by count, take top sequence
# Add top_core_kmer named data.table
## data.table will contain top_core_kmer and count
# For each kmer in core_pool
  # Compare sequence of top_core_kmer with kmer
    # if differences greater than max_diff
    # Reverse complement kmer
    # Compare sequence of top_core_kmer and reverse complement kmer
      # if differences greater than max_diff
        # do nothing
      # else if < max_diff, add kmer to table
  # else add kmer to top_core_kmer table in P-cloud list

# Initiate empty P-cloud list object
# core_pool = Filter k-mers > core_cutoff
# while core_pool is NOT empty do;
  # Apply P-cloud Core ID (with max_diff = 3) function to core_pool
    # Return first core data.table
  # Add first core data.table to P-cloud list object
  # Filter the filtered core_pool by kmers inside first core data.table
# Repeat until core_pool is empty
  # Returns a list of core data.tables in P-cloud list

# outer_layer_pool = Filter k-mers < core_cutoff and > lower cutoff
# for each data.table in P-cloud core list
  # if top core count > tertiary cutoff
    # set max_diff to 3
    # else if top core count > secondary cutoff
      # set max_diff to 2
        # else if top core count > primary cutoff
          # set max_diff to 1
            # else if top core count < primary cutoff
              # max_diff is 0 and move onto next data.table
    # for each k-mer in data.table
      # Apply comparing sequence function w/ set max_diff to list of outer_layer_pool
        # For k-mers < max_diff add to data.table and remove from outer_layer_pool
        # else do nothing
# Return








