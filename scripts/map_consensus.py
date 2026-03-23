#!/usr/bin/env python3

import argparse

from magnumopus import map_reads_to_ref

def parse_args():

    p = argparse.ArgumentParser(
        description="Extract 16S sequences read mapping"
        )
    p.add_argument(
        "-1", "--read1",
        required=True,
        help="Path to read1"
        )
    p.add_argument(
        "-2", "--read2",
        required=True,
        help="Path to read2"
        )
    p.add_argument(
        "-r", "--ref-seqs",
        required=True,
        help="Path to the reference sequences file for use in read mapping"
        )
    p.add_argument(
        "-s", "--seq",
        required=False,
        help="Name of sequence for which to print the consensus"
        )
    

    args = p.parse_args()
    return args


def main():
    args = parse_args()
    
    sam = map_reads_to_ref(
        ref=args.ref_seqs,
        r1=args.read1,
        r2=args.read2
    )
    if args.seq:
        print(sam.consensus(args.seq, fasta=True), end="")
    else:
        print(sam.best_consensus(fasta=True), end="")


if __name__ == "__main__":
    main()
