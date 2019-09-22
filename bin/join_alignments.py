#!/usr/bin/env python3

import os
import sys
import argparse
from collections import defaultdict, namedtuple

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


Partition = namedtuple("Partition", ["og", "start", "end"])


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        Creates a fasta and grouping file for iqtree given a bed-like file.
        """
    )

    parser.add_argument(
        "infiles",
        type=str,
        nargs="+",
        help="Input fasta files.",
    )

    parser.add_argument(
        "-o", "--outfasta",
        required=True,
        type=argparse.FileType('w'),
        help="Output file path.",
    )

    parser.add_argument(
        "-p", "--outpartition",
        required=True,
        type=argparse.FileType('w'),
        help="Output file path.",
    )

    parser.add_argument(
        "-t", "--trim",
        action="store_true",
        default=False,
        help="Trim alignments to a multiple of three for codon alignments.",
    )

    parser.add_argument(
        "-e", "--exclude",
        action="store_true",
        default=False,
        help=("Don't include any alignments that don't have "
              "lengths in multiples of 3."),
    )

    parser.add_argument(
        "--type",
        type=str,
        default="DNA",
        help="The molecule type to include in the partition file."
    )

    return parser.parse_args(args)


def decide_if_exclude(align):
    return align.get_alignment_length() % 3 != 0


def decide_trim_amount(align):
    return align.get_alignment_length() % 3


def parse_alignments(alignment_files, trim, exclude):
    isolates = set()
    alignments = list()

    for alignment_file in alignment_files:
        alignment = AlignIO.read(alignment_file, format="fasta")

        if exclude and decide_if_exclude(alignment):
            continue

        if trim:
            trim_amount = decide_trim_amount(alignment)
        else:
            trim_amount = 0

        og = os.path.splitext(alignment_file)[0]
        row = {}
        for seq in alignment:
            isolates.add(seq.id)

            if trim_amount == 0:
                row[seq.id] = str(seq.seq)
            else:
                last_pos = alignment.get_alignment_length() - trim_amount
                row[seq.id] = str(seq.seq)[:last_pos]

            alignments.append((og, row))

    return isolates, alignments


def make_blank_seq(length):
    return "-" * length


def join_alignments(isolates, alignments):
    partitions = list()
    joined_alignments = defaultdict(list)

    i = 0
    for og, alignment in alignments:
        length = None

        for isolate, seq in alignment.items():
            joined_alignments[isolate].append(seq)
            length = len(seq)

        missing_isolates = isolates.difference(alignment.keys())
        for isolate in missing_isolates:
            seq = make_blank_seq(length)
            joined_alignments[isolate].append(seq)

        start = i
        i += length
        end = i
        partition = Partition(og, start, end)
        partitions.append(partition)

    joined_seqs = []
    for isolate, seqs in joined_alignments.items():
        seq = Seq("".join(seqs))
        sr = SeqRecord(id=isolate, name=isolate, description=isolate, seq=seq)
        joined_seqs.append(sr)

    return MultipleSeqAlignment(joined_seqs), partitions


def write_partition(partitions, handle, kind="DNA"):
    strings = [
        f"{kind}, {p.og} = {p.start + 1}-{p.end}"
        for p
        in partitions
    ]
    handle.write("\n".join(strings))
    return


def main():
    args = cli(sys.argv[0], sys.argv[1:])
    isolates, alignments = parse_alignments(
        args.infiles,
        args.trim,
        args.exclude
    )

    msa, partitions = join_alignments(isolates, alignments)
    AlignIO.write(msa, args.outfasta, format="fasta")
    write_partition(partitions, args.outpartition, kind=args.type)
    return


if __name__ == "__main__":
    main()
