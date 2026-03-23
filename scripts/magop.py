#!/usr/bin/env python3

import argparse
import os
import re
import sys
from tempfile import TemporaryDirectory
from collections import OrderedDict

import numpy as np
from biotite.sequence.phylo import neighbor_joining

from magnumopus.mapping import map_reads_to_ref
from magnumopus.nw import needleman_wunsch
from magnumopus.ispcr import step_one, step_two, step_three

# Scoring for alignment
MATCH, MISMATCH, GAP = 1, -1, -1


def ispcr(primers, template, max_size):
    hits = step_one(primers, template)
    pairs = step_two(hits, max_size)
    return step_three(pairs, template)


def parse_fasta(text):
    header, seq = None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                yield header, "".join(seq)
            header = line[1:].strip()
            seq = []
        else:
            seq.append(line)
    if header:
        yield header, "".join(seq)


def group_reads(files):
    pairs = {}
    for f in files:
        m = re.match(r"(.+)_([12])\.(fastq|fq)(\.gz)?$", os.path.basename(f))
        if not m:
            raise ValueError("Bad filename: {}".format(f))
        sid, num = m.group(1), m.group(2)
        pairs.setdefault(sid, {})[num] = f
    return {k: (v["1"], v["2"]) for k, v in pairs.items()}


def clean_label(raw):
    h = raw.split()[0]
    h = h.split(";", 1)[0]
    h = re.sub(r"_[0-9]+-[0-9]+$", "", h)
    h = h.split(":", 1)[0]
    return h


def revcomp(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def distance(a, b):
    (a1, b1), _ = needleman_wunsch(a, b, MATCH, MISMATCH, GAP)
    mismatches_fwd = sum(x != y for x, y in zip(a1, b1))

    b_rc = revcomp(b)
    (a2, b2), _ = needleman_wunsch(a, b_rc, MATCH, MISMATCH, GAP)
    mismatches_rc = sum(x != y for x, y in zip(a2, b2))

    d = float(min(mismatches_fwd, mismatches_rc))
    if d < 0:
        d = 0.0
    return d


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--assemblies", nargs="*")
    parser.add_argument("-p", "--primers", required=True)
    parser.add_argument("-r", "--reads", nargs="*")
    parser.add_argument("-s", "--ref_seqs")
    parser.add_argument("-m", "--max_size", type=int, default=400)
    args = parser.parse_args()

    if args.reads and not args.ref_seqs:
        sys.exit("ERROR: -s/--ref_seqs is required when using -r/--reads")

    label_to_seq = OrderedDict()

    # Assemblies
    if args.assemblies:
        for a in args.assemblies:
            fa = ispcr(args.primers, a, args.max_size)
            for h, s in parse_fasta(fa):
                label = clean_label(h)
                if label not in label_to_seq:
                    label_to_seq[label] = s

    # Reads (ERR/SRR)
    if args.reads:
        pairs = group_reads(args.reads)
        for sid, (r1, r2) in pairs.items():
            sam = map_reads_to_ref(ref=args.ref_seqs, r1=r1, r2=r2)
            cons = sam.best_consensus(fasta=False)

            if not cons or not cons.strip():
                print(f"WARNING: no consensus for {sid}, skipping", file=sys.stderr)
                continue

            with TemporaryDirectory() as t:
                fp = os.path.join(t, "consensus.fasta")
                with open(fp, "w") as f:
                    f.write(f">{sid}\n{cons}\n")

                try:
                    fa = ispcr(args.primers, fp, args.max_size)
                except RuntimeError as e:
                    print(f"WARNING: isPCR failed for {sid}: {e}", file=sys.stderr)
                    continue

                for _, s in parse_fasta(fa):
                    label = sid
                    if label not in label_to_seq:
                        label_to_seq[label] = s
                    break

    # Reference sequences
    if args.ref_seqs:
        fa = ispcr(args.primers, args.ref_seqs, args.max_size)
        for h, s in parse_fasta(fa):
            label = clean_label(h)
            if label not in label_to_seq:
                label_to_seq[label] = s

    if not label_to_seq:
        sys.exit("ERROR: No amplicons found")

    labels = list(label_to_seq.keys())
    seqs = list(label_to_seq.values())

    labels, seqs = zip(*sorted(zip(labels, seqs), key=lambda x: x[0]))
    labels = list(labels)
    seqs = list(seqs)

    n = len(seqs)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = distance(seqs[i], seqs[j])
            if abs(d) < 1e-8:
                d = 0.0
            D[i][j] = D[j][i] = d

    tree = neighbor_joining(D)
    print(tree.to_newick(labels=labels, round_distance=2))


if __name__ == "__main__":
    main()