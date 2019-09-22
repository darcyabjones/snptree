#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(DECIPHER)))

# Round iranges down to nearest multiple of 3
offset_start <- function(x) {
  start(x) - ((start(x) - 1) %% 3)
}

# Round iranges up to nearest multiple of 3
offset_end <- function(x) {
  end(x) - (end(x) %% -3)
}


args = commandArgs(trailingOnly=TRUE)

infile <- args[1]
outfile <- args[2]


dna <- readDNAStringSet(infile)

# De-replicate the input sequences.
u_dna <- unique(dna)
index <- match(dna, u_dna)

# Do codon alignment of the dereplicated sequences
U_DNA <- AlignTranslation(u_dna)

# Re-replicate the sequences.
DNA <- U_DNA[index]
names(DNA) <- names(dna)

# Mask alignments where entropy < 1 or gap fraction < 0.2
MDNA <- MaskAlignment(DNA, showPlot = FALSE)

# Fix the mask to make sure codon phase is preserved.
out_mask <- as(colmask(MDNA), "IRanges")
start(out_mask) <- offset_start(out_mask)  # Bumps start to start of codon
end(out_mask) <- offset_end(out_mask)  # Bumps end to end of codon.
colmask(MDNA) <- asNormalIRanges(out_mask, force=TRUE)

MASKED <- as(MDNA, "DNAStringSet")

# Write to file.
writeXStringSet(MASKED, outfile)
