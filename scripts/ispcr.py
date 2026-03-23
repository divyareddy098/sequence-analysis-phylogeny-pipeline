#!/usr/bin/env python3

import subprocess
import tempfile
from typing import List, Tuple


# BLAST: find primer hits 

def step_one(primer_file: str, assembly_file: str) -> List[List[str]]:
    """
    Find primer matches using BLAST.
    """
    cmd = [
        "blastn",
        "-task", "blastn-short",
        "-word_size", "6",
        "-penalty", "-2",
        "-query", primer_file,
        "-subject", assembly_file,
        "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qlen",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"blastn failed (code {result.returncode})\n"
            f"cmd: {' '.join(cmd)}\n"
            f"stderr:\n{result.stderr}"
        )

    hits: List[List[str]] = []

    for line in result.stdout.splitlines():
        cols = line.strip().split("\t")
        if len(cols) < 13:
            continue

        pident = float(cols[2])
        length = int(cols[3])
        qlen = int(cols[12])

        # full-length match of the primer, >= 80% identity
        if length == qlen and pident >= 80.0:
            hits.append(cols)

    return hits


#  Pair forward & reverse primers

def step_two(hits: List[List[str]], max_size: int) -> List[Tuple[str, int, int]]:
    """
    Pair forward and reverse primers.
    """
    per_contig = {}  

    for h in hits:
        qid = h[0]      
        contig = h[1]
        sstart = int(h[8])
        send = int(h[9])

        strand = "+" if sstart <= send else "-"
        s5 = min(sstart, send)
        s3 = max(sstart, send)

        per_contig.setdefault(contig, []).append(
            {"qid": qid, "strand": strand, "s5": s5, "s3": s3}
        )

    paired_regions: List[Tuple[str, int, int]] = []

    for contig, hits_list in per_contig.items():
        fwd_hits = [h for h in hits_list if h["strand"] == "+"]
        rev_hits = [h for h in hits_list if h["strand"] == "-"]

        for f in fwd_hits:
            f_3p = f["s3"]  
            for r in rev_hits:
                r_3p = r["s5"]  
                if r_3p <= f_3p:
                    continue
                length = r_3p - f_3p
                if 0 < length <= max_size:
                    # 0-based, half-open interval
                    start0 = f_3p
                    end0 = r_3p
                    paired_regions.append((contig, start0, end0))

    return paired_regions


#  Extract amplicons with seqtk 

def step_three(hit_pairs: List[Tuple[str, int, int]],
               assembly_file: str) -> str:
    """
    Extract amplicons from the assembly.
    """
    if not hit_pairs:
        return ""

    bed_lines = []
    for contig, start0, end0 in hit_pairs:
        if start0 < end0:
            bed_lines.append(f"{contig}\t{start0}\t{end0}")

    if not bed_lines:
        return ""

    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as bed_file:
        bed_file.write("\n".join(bed_lines))
        bed_path = bed_file.name

    cmd = ["seqtk", "subseq", assembly_file, bed_path]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"seqtk failed (code {result.returncode})\n"
            f"cmd: {' '.join(cmd)}\n"
            f"stderr:\n{result.stderr}"
        )

    return result.stdout