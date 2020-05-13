#!/usr/bin/env python

from argparse import ArgumentParser
from Bio.Seq import reverse_complement


def revcomp(seq):
    return reverse_complement(seq)


def main(args):
    r = revcomp(args.seq)
    print(r, end="")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("seq", type=str,
                        help="Input fasta string")
    args = parser.parse_args()
    main(args)
