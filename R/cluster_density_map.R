#' Demarcate repetitive regions from kmer clusters
#'
#' @param cluster_table output from cluster_kmers. Table with 4 cols: p_cloud, kmer, count, positions
#' @param threshold value from 1-10 that determines what's allowed in the smoothing algorithm.
#'
#' @return kmer density table with 10 columns: start_position, p_cloud, kmer, count, p_cloud_kmers, rep_contig_group, rep_contig_length, max_p_cloud_kmers, rep_contig_start, rep_contig_end
#' @export
#'
#' @examples
#' testClusters <- cluster_kmers(core_cutoff = 5, kmer_counts = testCounts)
#' identify_cluster_density <- identify_cluster_density(testClusters, 8)
identify_cluster_density <- function(cluster_table, threshold = 8){
  # Find length of kmers in given table
  k <- cluster_table[,unique(nchar(kmer))]

  # Set position table to data.table where start_position is unlist of positions
  position_table <- cluster_table[, .(start_position = unlist(positions),
         p_cloud = p_cloud,
         kmer = kmer,
         count = count), by = .I][,I:=NULL]

  # Sort by position
  position_table <- position_table[order(start_position)]
  # Calculate kmer density window
  density_cluster_table <- position_table[, p_cloud_kmers := runner::runner(
    x = kmer,
    k = 10,
    idx = start_position,
    f = function(x) length(unique(x)))]

  # Calculate number of continuous kmer stretches in p_clouds
  density_cluster_table <- density_cluster_table[, rep_contig_group :=
                                                   cumsum(c(1, diff(start_position) != 1))]

  # Calculate length of continuous kmer stretches in p_clouds
  density_cluster_table <- density_cluster_table[, rep_contig_length := .N, by = rep_contig_group]

  # Calculate max_p_cloud_kmers per rep_contig_group
  density_cluster_table <- density_cluster_table[, max_p_cloud_kmers := max(p_cloud_kmers), by = rep_contig_group]

  # Filter table based on threshold. Removes contigs < threshold length AND aren't in a high density region
  density_cluster_table <- density_cluster_table[(rep_contig_length >= threshold) | (max_p_cloud_kmers >= threshold)]

  # Add start and end of each contig group after filtering
  density_cluster_table <- density_cluster_table[, `:=` (rep_contig_start = min(start_position),
                                                       rep_contig_end = max(start_position) + k - 1),
                                                       by = rep_contig_group]
  # Identify gaps and combine these small groups to larger proximal contig groups
  # Rows that are continuous less than threshold
  small_clusters <- density_cluster_table[rep_contig_length < threshold, unique(rep_contig_group)]
  for(cluster in small_clusters) {
    # Proximal groups to small clusters
    merging_clusters <- c(cluster + 1, cluster - 1)
    # Not all proximal clusters will be large contig groups, may have already been filtered out
    existing_clusters <- density_cluster_table[rep_contig_group %in% merging_clusters, unique(rep_contig_group)]
    # assign small cluster to the nearest existing_cluster
    if (length(existing_clusters) > 1) {
      cluster_start <- density_cluster_table[rep_contig_group == cluster, unique(rep_contig_start)]
      cluster_end <- density_cluster_table[rep_contig_group == cluster, unique(rep_contig_end)]
      above_cluster_start <- density_cluster_table[ rep_contig_group == merging_clusters[1], unique(rep_contig_start)]
      below_cluster_end <- density_cluster_table[rep_contig_group == merging_clusters[2], unique(rep_contig_end)]
      above_diff <- above_cluster_start - cluster_end
      below_diff <- cluster_start - below_cluster_end
      if (above_diff < below_diff) {
        # replace with above_cluster
        density_cluster_table[rep_contig_group == cluster, rep_contig_group := merging_clusters[1]]
      } else {
        # replace with below_cluster
        density_cluster_table[rep_contig_group == cluster, rep_contig_group := merging_clusters[2]]
      }
    } else if (length(existing_clusters) > 0) {
      # assign small cluster to the only proximal existing_cluster
      density_cluster_table[rep_contig_group == cluster, rep_contig_group := existing_clusters]
    } else if (length(existing_clusters) == 0) {
      # remove small cluster from table (no proximal big groups to join)
      density_cluster_table <- density_cluster_table[rep_contig_group != cluster]
    }
  }

  # with the updated groups, find the new start and end
  density_cluster_table <- density_cluster_table[, `:=` (rep_contig_start = min(start_position),
                                                         rep_contig_end = max(start_position) + k - 1),
                                                 by = rep_contig_group]

  return(density_cluster_table)

}

export_clusters_fasta <- function(density_cluster_table, seq_path, out_path) {
  # Import fasta sequence as biostrings object
  genome_sequence <- Biostrings::readDNAStringSet(seq_path)
  # Select only columns of interest
  rep_contigs <- density_cluster_table[,.(rep_contig_group, rep_contig_start, rep_contig_end)] %>%
    unique() # get rid of duplicates

  fasta_sequences <- Biostrings::DNAStringSet()

  # Extract sequences for each contig group and add to the DNAStringSet
  for (i in 1:nrow(rep_contigs)) {
    contig <- rep_contigs[i]
    seq <- Biostrings::subseq(genome_sequence, start = contig$rep_contig_start, end = contig$rep_contig_end)
    names(seq) <- paste0("rep_seq_", contig$rep_contig_group)
    fasta_sequences <- append(fasta_sequences, seq)
  }

  # Export sequences to out_path
  Biostrings::writeXStringSet(fasta_sequences, filepath = out_path)
  # Return sequences as biostrings object
  return(fasta_sequences)

}
