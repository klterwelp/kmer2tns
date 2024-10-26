## code to prepare phiX fasta file internal tests
# download fasta from NCBI
phiX_fasta <- rentrez::entrez_fetch(db = "nucleotide", id = "NC_001422.1", rettype = "fasta")
writeLines(phiX_fasta, "inst/extdata/phiX.fasta")
# compress to save space
system("gzip inst/extdata/phiX.fasta")


