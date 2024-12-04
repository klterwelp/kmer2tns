## code to prepare test sequence with known repetitive elements as a fasta

# generate random sequences with set seed
set.seed(10)
random_sequence <- stringi::stri_rand_strings(1, 1*10^6, '[ACGT]')
set.seed(5)
random_repeat <- stringi::stri_rand_strings(1, 100, '[ACGT]')


# positions in sequence to be replaced with random_repeat
replace_positions <- c(1, 200, 500, 1000, 2000, 3000, 5000, 8000, 10000, 11000,
                       30000, 80000, 50000, 100000, 900000)


# replace the sequence at specified positions
for (pos in replace_positions) {
  end_pos <- pos + 99
  random_sequence <- stringi::stri_sub_replace(random_sequence, from = pos,
                                               to = end_pos,
                                               value = random_repeat)
}

# export random_sequence into fasta
random_dna <- Biostrings::DNAStringSet(x = random_sequence)
Biostrings::writeXStringSet(random_dna,
                            filepath = "inst/extdata/test_sequence.fasta")
