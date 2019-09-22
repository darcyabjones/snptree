#!/usr/bin/env python3

import sys
import argparse
from collections import Counter
from collections import defaultdict


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        Creates a fasta and grouping file for iqtree given a bed-like file.
        """
    )

    parser.add_argument(
        "infile",
        default=sys.stdin,
        type=argparse.FileType('r'),
        help="Input fasta files. Use '-' for stdin.",
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

    return parser.parse_args(args)


HEADER = ["seqid", "start", "end"]
GROUP_HEADER = ["gseqid", "gstart", "gend", "gene"]


class Record(object):

    def __init__(self, seqid, start, end, group, samples):
        self.seqid = str(seqid)
        self.start = int(start)
        self.end = int(end)
        self.samples = list(samples)
        self.group = str(group)
        return

    def is_useful(self):
        return any((s != self.samples[0]) for s in self.samples)

    def impute(self):
        counter = Counter(s for s in self.samples if s != ".")
        base, freq = counter.most_common(1)[0]

        self.replace_missing(base, missing=".")
        return

    def replace_missing(self, replacement, missing="."):
        for i, sample in enumerate(self.samples):
            if sample == missing:
                self.samples[i] = replacement

        return


class Reader(object):

    def __init__(self, samples, group_col="gene"):
        self.samples = samples
        self.group_col = group_col
        return

    @classmethod
    def from_header(cls, header, group_col="gene"):
        hcols = header.strip("# ").split("\t")
        return cls(hcols[len(HEADER):], group_col=group_col)

    def read_line(self, line):
        sline = line.strip().split("\t")
        kwargs = {}

        kwargs.update(dict(zip(HEADER, sline)))
        group_info = dict(zip(GROUP_HEADER, sline[-len(GROUP_HEADER):]))
        samples = sline[len(HEADER): -len(GROUP_HEADER)]
        return Record(
            group=group_info[self.group_col],
            samples=samples,
            **kwargs
        )


def loop_infile(reader, infile):
    for line in infile:
        rec = reader.read_line(line)
        if rec.is_useful():
            rec.replace_missing("-", missing=".")
            rec.replace_missing("-", missing="*")
            yield rec
    return


def group_records(records):
    groups = defaultdict(list)
    for record in records:
        groups[record.group].append(record)
    return groups


def groups_to_fasta(groups, sample_order):
    fasta = defaultdict(list)
    spans = dict()

    # Format is 1 based inclusive ends
    i = 1
    for group, records in groups.items():
        records.sort(key=lambda x: x.start)
        start = i
        for record in records:
            for sample_name, sample_base in zip(sample_order, record.samples):
                fasta[sample_name].append(sample_base)
            i += 1

        end = i - 1
        spans[group] = (start, end)

    return format_fasta(fasta), format_partition(spans)


def format_fasta(fasta_dict):
    for key, value in fasta_dict.items():
        yield ">{}\n{}".format(key, "".join(value))
    return


def format_partition(span_dict):
    for group, (start, end) in span_dict.items():
        yield "DNA, {} = {}-{}".format(group, start, end)
    return


def main():
    args = cli(sys.argv[0], sys.argv[1:])

    header_line = next(args.infile)
    reader = Reader.from_header(header_line)
    records = loop_infile(reader, args.infile)
    grouped = group_records(records)

    fasta, partitions = groups_to_fasta(grouped, reader.samples)

    args.outfasta.write("\n".join(fasta))
    args.outpartition.write("\n".join(partitions))
    return


if __name__ == "__main__":
    main()
