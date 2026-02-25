# Copyright (C) 2026 Berkem Toprak Elmacı
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License.

import sys
import argparse
from pathlib import Path

from fasta_analyzer import FastaAnalyzer

# HELPERS

_COMP       = str.maketrans("ACGT", "TGCA")
STOP_CODONS = {"TAA", "TAG", "TGA"}


def reverse_complement(seq):
    return seq.translate(_COMP)[::-1]


def parse_fasta(filepath):
    records = []
    with open(filepath, "r") as fh:
        header, seq_parts = None, []
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header, seq_parts = line[1:], []
            else:
                seq_parts.append(line)
        if header is not None:
            records.append((header, "".join(seq_parts).upper()))
    return records

# ORF

def find_orfs(seq, min_len=100):
    """
    Search all 6 reading frames for ORFs.
    Returns list of dicts sorted by length (descending).
    Coordinates are 1-based on the forward strand.
    """
    orfs    = []
    seq_len = len(seq)
    rc_seq  = reverse_complement(seq)

    for strand, working_seq in [("+", seq), ("-", rc_seq)]:
        for frame in range(3):
            i = frame
            while i + 3 <= len(working_seq):
                if working_seq[i:i+3] == "ATG":
                    for j in range(i, len(working_seq) - 2, 3):
                        if working_seq[j:j+3] in STOP_CODONS:
                            orf_len = j + 3 - i
                            if orf_len >= min_len:
                                if strand == "+":
                                    start, end = i + 1, j + 3
                                else:
                                    end   = seq_len - i
                                    start = seq_len - (j + 3) + 1
                                orfs.append({
                                    "strand": strand,
                                    "frame":  frame + 1,
                                    "start":  start,
                                    "end":    end,
                                    "length": orf_len,
                                })
                            i = j + 3
                            break
                    else:
                        break
                i += 3

    return sorted(orfs, key=lambda x: -x["length"])

# DISPLAY

def run(records, min_len=100, all_frames=False):
    print("=" * 60)
    print(f"  ORF FINDER  (min length: {min_len} bp, 6 frames)")
    print("=" * 60)

    for header, seq in records:
        print(f"\n>{header}  [{len(seq)} bp]")
        orfs = find_orfs(seq, min_len=min_len)

        if not orfs:
            print("  No ORFs found.\n")
            continue

        longest = orfs[0]
        print(f"  Total ORFs    : {len(orfs)}")
        print(f"  Longest ORF   : {longest['length']} bp  "
              f"(strand {longest['strand']}, frame {longest['frame']}, "
              f"pos {longest['start']}–{longest['end']})")
        print()
        print(f"  {'#':<4} {'Strand':<7} {'Frame':<6} {'Start':<10} {'End':<8} {'Length'}")
        print("  " + "-" * 48)

        display = orfs if all_frames else orfs[:20]
        for idx, orf in enumerate(display, 1):
            print(f"  {idx:<4} {orf['strand']:<7} {orf['frame']:<6} "
                  f"{orf['start']:<10} {orf['end']:<8} {orf['length']} bp")

        if not all_frames and len(orfs) > 20:
            print(f"  ... ({len(orfs) - 20} more ORFs — use --all-frames to show all)")

# CLI

def build_parser():
    parser = argparse.ArgumentParser(
        prog="orf_analyzer.py",
        description="ORF Finder — 6 reading frames (3 forward + 3 reverse complement)")
    parser.add_argument("file",        nargs="?",  help="Input .fasta file")
    parser.add_argument("--min",       type=int,   default=100, metavar="BP",
                        help="Minimum ORF length in bp (default: 100)")
    parser.add_argument("--all-frames", action="store_true",
                        help="Show all ORFs (default: top 20 per sequence)")
    return parser


if __name__ == "__main__":
    parser = build_parser()
    args   = parser.parse_args()

    if args.file:
        p = Path(args.file)
        if not p.exists():
            sys.exit(f"ERROR: File not found: {args.file}")
        filepath = p
    else:
        filelist = list(Path().glob("*.fasta"))
        if len(filelist) == 1:
            filepath = filelist[0]
        elif len(filelist) == 0:
            sys.exit("ERROR: No .fasta file found.")
        else:
            sys.exit("ERROR: Multiple .fasta files found. Please specify one.")

    records = parse_fasta(filepath)
    if not records:
        sys.exit("ERROR: No sequences found in file.")

    print(f"\nFile      : {filepath}")
    print(f"Sequences : {len(records)}\n")
    run(records, min_len=args.min, all_frames=args.all_frames)
