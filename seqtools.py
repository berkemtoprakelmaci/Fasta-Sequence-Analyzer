# Copyright (C) 2026 Berkem Toprak Elmacı
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License.

import sys
import argparse
from pathlib import Path

try:
    from fasta_analyzer import FastaAnalyzer
    import orf_analyzer
    import gc_analyzer
    import codon_analyzer
except ImportError as e:
    sys.exit(f"ERROR: Could not import required module — {e}\n"
             f"Make sure all .py files are in the same directory.")


# FILE RESOLVER

def resolve_input(positional_file):
    if positional_file:
        p = Path(positional_file)
        if not p.exists():
            sys.exit(f"ERROR: File not found: {positional_file}")
        return p
    filelist = list(Path().glob("*.fasta"))
    if len(filelist) == 1:
        return filelist[0]
    elif len(filelist) == 0:
        sys.exit("ERROR: No .fasta file found. Please specify the file as an argument.")
    else:
        names = ", ".join(f.name for f in filelist)
        sys.exit(f"ERROR: Multiple .fasta files found ({names}).\n"
                 f"Usage: python seqtools.py <command> file.fasta")

# SUBCOMMAND HANDLERS

def cmd_stats(filepath, args):
    FastaAnalyzer(filepath)

def cmd_orf(filepath, args):
    records = orf_analyzer.parse_fasta(filepath)
    orf_analyzer.run(records, min_len=args.min, all_frames=args.all_frames)

def cmd_gc_window(filepath, args):
    records = gc_analyzer.parse_fasta(filepath)
    gc_analyzer.run(records, window=args.window, step=args.step)

def cmd_codon(filepath, args):
    records = codon_analyzer.parse_fasta(filepath)
    codon_analyzer.run(records)

# CLI

def build_parser():
    parser = argparse.ArgumentParser(
        prog="seqtools.py",
        description="FASTA Sequence Analysis Toolkit",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python seqtools.py stats      genome.fasta
  python seqtools.py orf        genome.fasta --min 300
  python seqtools.py orf        genome.fasta --all-frames
  python seqtools.py gc_window  genome.fasta --window 200 --step 50
  python seqtools.py codon      genome.fasta
        """)

    sub = parser.add_subparsers(dest="command")

    p_stats = sub.add_parser("stats", help="Basic sequence statistics (via fasta_analyzer.py)")
    p_stats.add_argument("file", nargs="?")

    p_orf = sub.add_parser("orf", help="ORF finder — 6 reading frames (via orf_analyzer.py)")
    p_orf.add_argument("file",         nargs="?")
    p_orf.add_argument("--min",        type=int, default=100, metavar="BP",
                       help="Minimum ORF length in bp (default: 100)")
    p_orf.add_argument("--all-frames", action="store_true",
                       help="Show all ORFs (default: top 20 per sequence)")

    p_gc = sub.add_parser("gc_window", help="Sliding window GC analysis + plot (via gc_analyzer.py)")
    p_gc.add_argument("file",     nargs="?")
    p_gc.add_argument("--window", type=int, default=100, metavar="BP",
                      help="Window size in bp (default: 100)")
    p_gc.add_argument("--step",   type=int, default=10,  metavar="BP",
                      help="Step size in bp (default: 10)")

    p_codon = sub.add_parser("codon", help="Codon usage table (via codon_analyzer.py)")
    p_codon.add_argument("file", nargs="?")

    return parser


if __name__ == "__main__":
    parser = build_parser()

    args_raw = sys.argv[1:]
    if args_raw and args_raw[0].endswith(".fasta"):
        args_raw = ["stats"] + args_raw

    args = parser.parse_args(args_raw)

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    filepath = resolve_input(args.file)

    print(f"\nFile      : {filepath}")
    print(f"Sequences : {sum(1 for line in open(filepath) if line.startswith('>'))}\n")

    dispatch = {
        "stats":     cmd_stats,
        "orf":       cmd_orf,
        "gc_window": cmd_gc_window,
        "codon":     cmd_codon,
    }
    dispatch[args.command](filepath, args)
