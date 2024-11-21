
#' Count kmers with position information
#'
#' @param seq_path path to fasta file
#' @param k length of k-mers
#'
#' @return data.table with three columns: kmers, count, and positions
#' @export
#'
#' @examples
#' seq_path <- fs::path_package("extdata", "phiX.fasta.gz", package = "kmer2tns")
#' count_kmers_positions(seq_path = seq_path)
count_kmers_positions <- function(seq_path, k = NA) {
  # Check file exists + file format is correct
  log_message("Checking fasta file...")

  # Check that file exists
  ## if file doesn't exist, quit function
  if (!file.exists(seq_path)) {
    stop("File not found: ", seq_path)
  }

  # Check that file is in fasta format
  ## if not in fasta format, quit function
  tryCatch(
    {
      sequence <- as.character(Biostrings::readDNAStringSet(seq_path))
    },
    warning = function(w) {
      message("Warning: ", conditionMessage(w))
    },
    error = function(e) {
      message(
        "Error: ", conditionMessage(e), "\n",
        "Check that file is in fasta format."
      )
      stop(call. = FALSE)
    }
  )
  ## readDNAStringSet will only read fasta/fastq files


  # Check that file contains only one sequence, concatenate if not
  tryCatch(
    {
      num_seqs <- length(sequence)
    },
    warning = function(w) {
      message("Warning: ", conditionMessage(w))
    },
    error = function(e) {
      message(
        "Error: ", conditionMessage(e), "\n",
        "Failed to calculate number of sequences in fasta file"
      )
      stop(call. = FALSE)
    }
  )

  if (num_seqs > 1) {
    warning(
      "There were: ", num_seqs, " sequences in the fasta file.", "\n",
      "Will concatenate sequences and count k-mers."
    )
    ## concatenate sequences
    tryCatch(
      {
        sequence <- paste(sequence, collapse = "")
      },
      warning = function(w) {
        message("Warning: ", conditionMessage(w))
      },
      error = function(e) {
        message(
          "Error: ", conditionMessage(e), "\n",
          "Issue concatenating sequences"
        )
        stop(call. = FALSE)
      }
    )
  }

  # Check that the length of the given sequence makes sense
  tryCatch(
    {
      seq_length <- nchar(sequence)
      if (seq_length <= 0) {
        message("Sequence length must be greater than 0.")
        stop(call. = FALSE)
      }
      if (is.na(seq_length)) {
        message("Sequence length exceeds maximum limit for R.")
        stop(call. = FALSE)
      }
      if (seq_length < 100) {
        message("Sequence length must be at least 100.") # may change MIN later
        stop(call. = FALSE)
      }
    },
    error = function(e) {
      message("Error calculating sequence length: ", conditionMessage(e))
      stop(call. = FALSE)
    }
  )
  # nchar() may already throw an error already for these checks

  log_message("Finished checking fasta file.")

  log_message("Checking or assigning value of k...")
  # Check that k is a number or assign k if left blank
  tryCatch(
    {
      # convert k to nearest numeric integer
      k <- round(as.numeric(k))
      # Calculate default k
      default_k <- round(log(seq_length, base = 4)) + 1
      # check that k is not NA
      # if k is NA, use default k
      # default k based on: https://academic.oup.com/bioinformatics/article/21/5/582/219976
      if (is.na(k)) {
        k <- default_k
        message(
          "Invalid input for k. Defaulting to k = ", k, "\n"
        )
      }
    },
    error = function(e) { # if as.numeric or round has an error, use default k
      message("Error: ", conditionMessage(e))
      k <- default_k
      message("Error converting input for k to numeric. Defaulting to k = ", k)
    }
  )


  #Check that k is not too large or too small
  if (k < default_k) {
    stop("Selected k is too small. Pick a value that is at least: ", default_k)
  } else if (k > 500) {
    stop("Selected k is too large. Pick a value that is less than: ", 500)
    # around 4^511 will be larger than the double max for R
    # even k=31 will be over R's integer limit for 4^k
  }

  log_message("Finished checking value of k...")
  log_message("k is now: ", k)

  # Determine hash table size
  log_message("Determining hash table size...")

  tryCatch({
    hash_size1 <- 4^k
    # this is the theoretical space of all k-mers of size k
    hash_size2 <- seq_length - k + 1
    # this is the number of k-mers that will be counted from a sequence
    # if all k-mers were unique, this would be exact # of k-mers counted

    # assign hash size to the smallest of the two options
    hash_size <- min(hash_size1, hash_size2)

    # check that the hash size is not larger than integer limit
    if (hash_size > .Machine$integer.max) {
      message("Error: hash size is larger than integer max.", "\n",
              "To fix this: use a smaller genome or value of k")
      stop(call. = FALSE)
    }
  }, error = function(e) { # if calculating hash size has an error
    message("Error calculating hash size: ", conditionMessage(e))
    stop(call. = FALSE)
  })

  log_message("Hash size will be: ", hash_size)

  # Generate hash table
  log_message("Creating hash table...")
  # Hash table to store counts greater > min_length
  hash_kmers <- new.env(hash = TRUE, parent = emptyenv(), size = hash_size)
  # Hash table for initial count
  pre_count <- new.env(hash = TRUE, parent = emptyenv(), size = hash_size)
  log_message("Hash table generation complete.")

  # Count k-mers loop [Direct-count method]
  log_message("Counting k-mers...")

  # Function to count k-mers with positions
  count_kmers_with_positions <- function(seq_length, k, sequence, hash_kmers, pre_count) {
    tryCatch({
      for (i in 1:(seq_length - k + 1)) { # for each kmer along DNA sequence
        kmer <- substr(sequence, i, i + k - 1) # identify kmer as substring between i and i + k - 1
        if (exists(kmer, envir = hash_kmers)) { # check if kmer exists in main hash table
          # If k-mer exists, update count and append position
          kmer_data <- hash_kmers[[kmer]]
          kmer_data$count <- kmer_data$count + 1
          kmer_data$positions <- c(kmer_data$positions, i)
          hash_kmers[[kmer]] <- kmer_data
        } else {
          if (exists(kmer, envir = pre_count)) { # check if kmer exists in pre count hash table
            # If k-mer exists, update count and append position to hash_kmers
            kmer_data <- pre_count[[kmer]]
            kmer_data$count <- kmer_data$count + 1
            kmer_data$positions <- c(kmer_data$positions, i)
            hash_kmers[[kmer]] <- kmer_data
          } else {
            # If k-mer doesn't exist, initialize with count = 1 and positions = [position] in pre_count
            pre_count[[kmer]] <- list(count = 1, positions = i)
          }
        }
      }

    }, error = function(e) { # if counting k-mers has an error
      message("Error counting k-mers: ", conditionMessage(e))
      stop(call. = FALSE)
    })
  }
  # Count k-mers and store their positions
  count_kmers_with_positions(seq_length, k, sequence, hash_kmers, pre_count)

  log_message("Finished counting k-mers.")


  log_message("Converting hash table to data.table object...")
  # Extract kmers and length of hash table
  kmers <- ls(hash_kmers)
  num_kmers <- length(kmers)
  # Preallocate vector size for kmer data
  kmer <- character(length = num_kmers)
  count <- integer(length = num_kmers)
  positions <- vector(mode = "list", length = num_kmers)
  # Loop over kmers and extract values for position and count
  for (i in seq_along(kmers)) {
    kmer_data <- hash_kmers[[kmers[i]]] # only extracts once
    kmer[i] <- kmers[i]
    count[i] <- kmer_data$count
    positions[[i]] <- kmer_data$positions
  }

  hash_kmers_df <- data.table(kmer = kmer, count = count, positions = positions)

  log_message("Converted hash table to data.table object.")

  return(hash_kmers_df)
}


