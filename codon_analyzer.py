# Copyright (C) 2026 Berkem Toprak Elmacı
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License.

import sys
import argparse
from pathlib import Path
from collections import Counter

# HELPERS

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

# CODON TABLE

CODON_TABLE = {
    "TTT":"Phe","TTC":"Phe","TTA":"Leu","TTG":"Leu",
    "CTT":"Leu","CTC":"Leu","CTA":"Leu","CTG":"Leu",
    "ATT":"Ile","ATC":"Ile","ATA":"Ile","ATG":"Met",
    "GTT":"Val","GTC":"Val","GTA":"Val","GTG":"Val",
    "TCT":"Ser","TCC":"Ser","TCA":"Ser","TCG":"Ser",
    "CCT":"Pro","CCC":"Pro","CCA":"Pro","CCG":"Pro",
    "ACT":"Thr","ACC":"Thr","ACA":"Thr","ACG":"Thr",
    "GCT":"Ala","GCC":"Ala","GCA":"Ala","GCG":"Ala",
    "TAT":"Tyr","TAC":"Tyr","TAA":"Stop","TAG":"Stop",
    "CAT":"His","CAC":"His","CAA":"Gln","CAG":"Gln",
    "AAT":"Asn","AAC":"Asn","AAA":"Lys","AAG":"Lys",
    "GAT":"Asp","GAC":"Asp","GAA":"Glu","GAG":"Glu",
    "TGT":"Cys","TGC":"Cys","TGA":"Stop","TGG":"Trp",
    "CGT":"Arg","CGC":"Arg","CGA":"Arg","CGG":"Arg",
    "AGT":"Ser","AGC":"Ser","AGA":"Arg","AGG":"Arg",
    "GGT":"Gly","GGC":"Gly","GGA":"Gly","GGG":"Gly",
}

# CODON LOGIC

def codon_usage(seq):
    """Count in-frame codons starting at position 0."""
    counts = Counter()
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3 and all(b in "ACGT" for b in codon):
            counts[codon] += 1
    return counts

# DISPLAY


def run(records):
    print("=" * 60)
    print("  CODON USAGE TABLE")
    print("=" * 60)

    for header, seq in records:
        print(f"\n>{header}  [{len(seq)} bp]")
        counts = codon_usage(seq)
        total  = sum(counts.values())

        if total == 0:
            print("  No valid codons found.\n")
            continue

        aa_groups = {}
        for codon, count in sorted(counts.items()):
            aa = CODON_TABLE.get(codon, "???")
            aa_groups.setdefault(aa, []).append((codon, count))

        print(f"  {'Codon':<6} {'AA':<6} {'Count':<8} {'Freq/1000'}")
        print("  " + "-" * 35)
        for aa in sorted(aa_groups):
            for codon, count in sorted(aa_groups[aa], key=lambda x: -x[1]):
                print(f"  {codon:<6} {aa:<6} {count:<8} {count / total * 1000:.2f}")

        print(f"\n  Total codons analyzed: {total}")

# CLI

def build_parser():
    parser = argparse.ArgumentParser(
        prog="codon_analyzer.py",
        description="Codon Usage Table — standard genetic code, frame +1")
    parser.add_argument("file", nargs="?", help="Input .fasta file")
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
    run(records)
