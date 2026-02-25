# Copyright (C) 2026 Berkem Toprak Elmacı
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License.

import sys
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

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

# GC WINDOW

def sliding_gc(seq, window, step):
    """Return (positions, gc_percentages) for a sliding window scan."""
    positions, gc_values = [], []
    for i in range(0, len(seq) - window + 1, step):
        chunk = seq[i:i+window]
        gc = (chunk.count("G") + chunk.count("C")) / window * 100
        positions.append(i + window // 2)
        gc_values.append(gc)
    return positions, gc_values

# DISPLAY + PLOT

def run(records, window=100, step=10, out_path="gc_window.png"):
    print("=" * 60)
    print(f"  SLIDING WINDOW GC  (window={window} bp, step={step} bp)")
    print("=" * 60)

    fig, axes = plt.subplots(len(records), 1,
                             figsize=(12, 3.5 * len(records)),
                             squeeze=False)

    for idx, (header, seq) in enumerate(records):
        ax = axes[idx][0]

        if len(seq) < window:
            print(f"\n  SKIPPED >{header}: sequence shorter than window "
                  f"({len(seq)} bp < {window} bp)")
            ax.set_visible(False)
            continue

        positions, gc_vals = sliding_gc(seq, window, step)
        overall = (seq.count("G") + seq.count("C")) / len(seq) * 100

        ax.fill_between(positions, gc_vals, overall,
                        where=[v >= overall for v in gc_vals],
                        alpha=0.4, color="#e74c3c", label="Above average")
        ax.fill_between(positions, gc_vals, overall,
                        where=[v < overall for v in gc_vals],
                        alpha=0.4, color="#3498db", label="Below average")
        ax.plot(positions, gc_vals, color="#2c3e50", linewidth=0.8)
        ax.axhline(overall, color="black", linewidth=1, linestyle="--",
                   label=f"Overall GC ({overall:.1f}%)")

        ax.set_title(f">{header[:60]}", fontsize=9, loc="left")
        ax.set_xlabel("Position (bp)", fontsize=8)
        ax.set_ylabel("GC content (%)", fontsize=8)
        ax.set_ylim(0, 100)
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(axis="y", alpha=0.3)

        print(f"\n  >{header}")
        print(f"    Length       : {len(seq):,} bp")
        print(f"    Overall GC   : {overall:.2f}%")
        print(f"    GC range     : {min(gc_vals):.2f}% – {max(gc_vals):.2f}%")

    plt.suptitle("Sliding Window GC Analysis", fontsize=11, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\n  Plot saved → {Path(out_path).resolve()}")

# CLI

def build_parser():
    parser = argparse.ArgumentParser(
        prog="gc_analyzer.py",
        description="Sliding Window GC Analysis — generates gc_window.png")
    parser.add_argument("file",      nargs="?", help="Input .fasta file")
    parser.add_argument("--window",  type=int,  default=100, metavar="BP",
                        help="Window size in bp (default: 100)")
    parser.add_argument("--step",    type=int,  default=10,  metavar="BP",
                        help="Step size in bp (default: 10)")
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
    run(records, window=args.window, step=args.step)
