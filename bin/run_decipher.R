#!/usr/bin/env Rscript

VERSION = "0.0.1"

suppressPackageStartupMessages(library("optparse"))
suppressWarnings(suppressPackageStartupMessages(library("Biostrings")))
suppressWarnings(suppressPackageStartupMessages(library("DECIPHER")))

option_list <- list(
    make_option(
        c("-i", "--infile"),
        type="character",
        action="store",
        help="The input fasta to align (required)."
    ),
    make_option(
        c("-o", "--outfile"),
        type="character",
        action="store",
        help="The output file to write to (required)."
    ),
    make_option(
        c("-g", "--gencode"),
        type="character",
        action="store",
        default="1",
        help="The genetic code to use to translate the codons. Default 1 (standard)"
    ),
    make_option(
        c("-m", "--maxgap"),
        type="numeric",
        action="store",
        default=0.2,
        help="The maximum proportion of gaps allowed in a column. Default 0.2."
    ),
    make_option(
        c("-e", "--minentropy"),
        type="numeric",
        action="store",
        default=1,
        help="The minimum entropy score for a column in the alignment. Default 1."
    ),
    make_option(
        "--version",
        type="logical",
        action="store_true",
        default=FALSE,
        help="Print version and exit.",
    )
)

parser <- OptionParser(
    usage = "%prog --infile in.fasta --outfile out.fasta",
    option_list = option_list
)

args <- parse_args(parser)

log_stderr <- function(...) {
  cat(sprintf(...), sep='', file=stderr())
}

quit_with_err <- function(...) {
  log_stderr(...)
  quit(save = "no", status = 1, runLast = FALSE)
}

validate_file <- function(path) {
  if (is.null(path)) {
    quit_with_err("Please provide required file")
  }
}

# Round iranges down to nearest multiple of 3
offset_start <- function(x) {
  start(x) - ((start(x) - 1) %% 3)
}

# Round iranges up to nearest multiple of 3
offset_end <- function(x) {
  end(x) - (end(x) %% -3)
}


main <- function(args) {
  if (args$version) {
    cat(VERSION, file=stdout())
    quit(save = "no", status = 0, runLast = FALSE)
  }

  validate_file(args$infile)
  validate_file(args$outfile)

  # Select a genetic code to use for translation
  GC <- Biostrings::getGeneticCode(args$gencode)

  dna <- readDNAStringSet(args$infile)

  # De-replicate the input sequences.
  u_dna <- unique(dna)
  index <- match(dna, u_dna)

  # Do codon alignment of the dereplicated sequences
  sink(stderr(), type = "output")
  U_DNA <- AlignTranslation(u_dna, geneticCode = GC)

  # Re-replicate the sequences.
  DNA <- U_DNA[index]
  names(DNA) <- names(dna)

  # Mask alignments where entropy < 1 or gap fraction < 0.2
  MDNA <- MaskAlignment(
    DNA,
    showPlot = FALSE,
    maxFractionGaps = args$maxgap,
    threshold = args$minentropy
  )

  # Fix the mask to make sure codon phase is preserved.
  out_mask <- as(colmask(MDNA), "IRanges")
  start(out_mask) <- offset_start(out_mask)  # Bumps start to start of codon
  end(out_mask) <- offset_end(out_mask)  # Bumps end to end of codon.
  colmask(MDNA) <- asNormalIRanges(out_mask, force=TRUE)

  MASKED <- as(MDNA, "DNAStringSet")

  # Write to file.
  writeXStringSet(MASKED, args$outfile)
}

main(args)
