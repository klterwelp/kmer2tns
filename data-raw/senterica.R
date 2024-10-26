## code to prepare `senterica` fasta file internal tests
# download fasta from NCBI
# list of NCBI IDs for Salmonella Enterica SA20021456
ids <- c("CP030219.1", "CP030220.1", "CP030221.1", "CP030222.1")
# function to fetch sequences
fetch_sequences <- function(id) {
  rentrez::entrez_fetch(db = "nucleotide", id = id, rettype = "fasta")
}
# fetch all sequences
sequences <- lapply(ids, fetch_sequences)
# combine all sequences into one string
combined_sequences <- paste(unlist(sequences), collapse = "\n")
# add to extdata
writeLines(combined_sequences, "inst/extdata/sentericaSA20021456.fasta")
# compress to save space
system("gzip inst/extdata/sentericaSA20021456.fasta")
