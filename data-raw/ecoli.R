## code to prepare ecoli fasta file internal tests
# download fasta from NCBI
ecoli_fasta <- rentrez::entrez_fetch(db = "nucleotide", id = "AP022815.1", rettype = "fasta")
writeLines(ecoli_fasta, "inst/extdata/ecoliJE86-ST05.fasta")
# compress to save space
system("gzip inst/extdata/ecoliJE86-ST05.fasta")
