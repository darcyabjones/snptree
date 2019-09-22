#!/usr/bin/env python3

import os
import sys
import argparse

import pandas as pd
from Bio import SeqIO


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        Selects orthologous groups to align based on busco results.
        """
    )

    parser.add_argument(
        "indir",
        type=str,
        help="The directory containing busco results.",
    )

    parser.add_argument(
        "-o", "--outdir",
        default="og_seqs",
        type=str,
        help="Where to store the sequences.",
    )

    parser.add_argument(
        "-m", "--max-missing",
        default=0.1,
        type=float,
        help="The maximum number of isolates that the OG can be missing in.",
    )

    parser.add_argument(
        "-q", "--quiet",
        default=False,
        action="store_true",
        help="Don't print the selected OGs.",
    )

    return parser.parse_args(args)


def read_busco_results(isolates, path):
    """ """

    dfs = []
    for isolate in isolates:
        ft_path = os.path.join(path, isolate, f"full_table_{isolate}.tsv")
        df = pd.read_csv(
            ft_path,
            sep="\t",
            comment="#",
            names=["id", "status", "contig", "start", "end", "score", "length"]
        )

        df["isolate"] = isolate
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)


def select_buscos(isolates, df, max_missing, quiet=False):
    """ """

    passed = []
    for id_, subdf in df.groupby("id"):
        if (subdf["status"] == "Complete").all():
            passed.append(id_)
            continue

        duplicated_isolates = (
            subdf
            .loc[subdf["status"] == "Duplicated", "isolate"]
            .unique()
        )

        if len(duplicated_isolates) > 0:
            if not quiet:
                print(
                    "OG",
                    id_,
                    "was duplicated in isolates:",
                    ", ".join(duplicated_isolates),
                    file=sys.stdout
                )
            continue

        nmissing = (subdf["status"].isin(["Missing", "Fragmented"])).sum()
        if nmissing / len(isolates) < max_missing:
            passed.append(id_)

    if not quiet:
        print("Selected OGs:")
        print("\n".join(passed))

    return passed


def seqs_not_all_same(seqs):
    any_false = False
    for seq in seqs[1:]:
        if seq.seq != seqs[0].seq:
            any_false = True
            break

    return any_false


def get_selected_seqs(ogs, isolates, path):
    """ """

    for og in ogs:
        seqs = []

        for isolate in isolates:
            seq_path = os.path.join(
                path,
                isolate,
                "single_copy_busco_sequences",
                f"{og}.fna"
            )

            if not os.path.isfile(seq_path):
                continue

            seq = None
            for s in SeqIO.parse(seq_path, format="fasta"):
                if seq is None:
                    seq = s
                elif len(s) > len(seq):
                    seq = s

            assert seq is not None
            seq.id = isolate
            seq.name = isolate
            seq.description = og

            seqs.append(seq)

        if seqs_not_all_same(seqs):
            yield (og, seqs)
    return


def main():

    args = cli(sys.argv[0], sys.argv[1:])

    isolates = [
        d
        for d
        in os.listdir(args.indir)
        if os.path.isdir(os.path.join(args.indir, d))
    ]

    reports = read_busco_results(isolates, args.indir)
    ogs = select_buscos(
        isolates,
        reports,
        args.max_missing,
        args.quiet
    )

    os.makedirs(args.outdir, exist_ok=True)
    for og, seq in get_selected_seqs(ogs, isolates, args.indir):
        out_path = os.path.join(args.outdir, f"{og}.fasta")
        SeqIO.write(seq, out_path, format="fasta")

    return


if __name__ == "__main__":
    main()
