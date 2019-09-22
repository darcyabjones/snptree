#!/usr/bin/env python3

import sys
import argparse
from statistics import mean

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        From a codon aligned sequence, removes sequence after the stop codons.
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
        "-g", "--gencode",
        type=str,
        default=1,
        help="The genetic code to use to find stop codons."
    )

    parser.add_argument(
        "-m", "--max-missing",
        dest="max_missing",
        default=1,
        type=float,
        help=(
            "Filter out codon columns where the proportion of gaps is "
            "this number or higher. 0 will filter out columns with any "
            "gaps, 1 will filter out columns that are all gaps (default), "
            "> 1 will do no filtering."
        )
    )

    return parser.parse_args(args)


class GapInCodonError(Exception):

    def __init__(self, msg):
        self.msg = msg
        return

    def __str__(self):
        return self.msg


class NoValidPhase(Exception):

    def __init__(self, msg):
        self.msg = msg
        return

    def __str__(self):
        return f"Could not find a valid phase for the codons. {self.msg}"


def iter_codons(seq):
    for i in range(0, len(seq), 3):
        j = i + 3
        codon = seq[i:j]
        if len(codon) != 3:
            break

        yield codon
    return


def truncate_after_stop(seq, gencode=1, gap="-"):
    seen_first = False
    seen_stop = False
    internal_gap = None
    out_codons = []
    for codon in iter_codons(seq):
        # This avoids translation error
        if codon == "---":
            pass
        elif "-" in codon:
            if internal_gap is not None:
                raise GapInCodonError(str(internal_gap))

            if seen_first:
                internal_gap = codon
            else:
                seen_first = True

        # We've already seen a gap, so that one wasn't the last.
        elif internal_gap is not None:
            raise GapInCodonError(str(internal_gap))

        elif codon.translate(table=gencode) == "*":
            seen_stop = True

        if seen_stop or "-" in codon:
            out_codons.append(gap * 3)
        else:
            seen_first = True
            out_codons.append(str(codon))

    return out_codons


def filter_codon_columns(msa_codons, max_missing=0.1):
    passed_columns = []
    for column in zip(*msa_codons):
        missing = mean(c == "---" for c in column)
        if missing <= max_missing and missing < 1:
            passed_columns.append(column)

    return zip(*passed_columns)


def filter_msa(msa, gencode=1, gap="-", max_missing=0.1, offset=0):
    msa_codons = [
        truncate_after_stop(s.seq[offset:], gencode=gencode, gap=gap)
        for s
        in msa
    ]

    filtered_codons = filter_codon_columns(
        msa_codons,
        max_missing=max_missing
    )

    filtered_seqs = []
    for orig_seq, codons in zip(msa, filtered_codons):
        new_seq = Seq("".join(codons))
        sr = SeqRecord(
            id=orig_seq.id,
            name=orig_seq.name,
            description=orig_seq.description,
            seq=new_seq
        )
        filtered_seqs.append(sr)

    return MultipleSeqAlignment(filtered_seqs)


def score_msa(codons):
    score = 0
    for column in zip(*codons):
        missing = mean("-" not in c for c in column)
        score += missing

    return score


def select_best_offset(msa, gencode=1, gap="-"):
    scores = []
    errors = []

    for offset in [0, 1, 2]:
        try:
            codons = [
                truncate_after_stop(seq.seq[offset:], gencode=gencode, gap=gap)
                for seq
                in msa
            ]

        except GapInCodonError as e:
            error = f"Invalid phase: {offset}. Codon: {str(e)}."
            errors.append(error)
            continue

        scores.append((offset, score_msa(codons)))

    if len(scores) == 0:
        joined = ' '.join(errors)
        raise NoValidPhase(joined)

    return max(scores, key=lambda x: x[1])[0]


def main():
    args = cli(sys.argv[0], sys.argv[1:])
    inmsa = AlignIO.read(args.infile, format="fasta")

    assert 0 <= args.max_missing, "Max missing must be >= 0."
    try:
        offset = select_best_offset(inmsa, gencode=args.gencode)

        if offset != 0:
            print(
                f"Using offset {offset} in file {args.infile.name}.",
                file=sys.stderr
            )

        outmsa = filter_msa(
            inmsa,
            offset=offset,
            gencode=args.gencode,
            max_missing=args.max_missing
        )

        if inmsa.get_alignment_length() > (outmsa.get_alignment_length() + 9):
            print(
                f"Trimmed more than 3 codons in {args.infile.name}",
                file=sys.stderr
            )

    except GapInCodonError as e:
        print(
            f"Encountered internal gap in {args.infile.name}, Codon: {str(e)}",
            file=sys.stderr
        )
        # 255 because it's the only exit code xargs respects.
        sys.exit(255)
    except NoValidPhase as e:
        print(
            f"Could not find valid phase in {args.infile.name}.",
            f"{str(e)}",
            file=sys.stderr
        )
        # 255 because it's the only exit code xargs respects.
        sys.exit(255)

    AlignIO.write(outmsa, args.outfile, "fasta")
    return


if __name__ == "__main__":
    main()
