#!/usr/bin/env python3

import sys
import argparse

from Bio import AlignIO

def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        From a codon aligned sequence return a raxml style partition file.
        """
    )

    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="Input fasta file.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output file path.",
    )

    parser.add_argument(
        "-n", "--name",
        type=str,
        default="gene",
        help="The basename to prefix each partition name with."
    )

    return parser.parse_args(args)


def main():

    try:
        args = cli(sys.argv[0], sys.argv[1:])
        inmsa = AlignIO.read(args.infile, format="fasta")

        length = inmsa.get_alignment_length()
        print(f"DNA, {args.name}_1 = 1-{length}\\3", file=args.outfile) 
        print(f"DNA, {args.name}_2 = 2-{length}\\3", file=args.outfile) 
        print(f"DNA, {args.name}_3 = 3-{length}\\3", file=args.outfile) 
    except Exception as e:
        print("Error: infile={args.infile} outfile={args.outfile}", file=sys.stderr)
        print(e)
        sys.exit(255)


if __name__ == "__main__":
    main()
