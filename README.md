# FASTA Sequence Analysis Toolkit

A modular command-line toolkit for FASTA sequence analysis.  
Each analysis module works standalone or through the unified `seqtools.py` entry point.

---

## File Structure

```
.
├── seqtools.py         # Main CLI — delegates to all modules
├── fasta_analyzer.py   # Basic sequence statistics (standalone or imported)
├── orf_analyzer.py     # ORF detection — 6 reading frames
├── gc_analyzer.py      # Sliding window GC analysis + plot
├── codon_analyzer.py   # Codon usage frequency table
├── environment.yml     # Conda environment
└── README.md
```

---

## Installation

```bash
conda env create -f environment.yml
conda activate fasta-toolkit
```

---

## Usage

### Via seqtools.py (recommended)

All modules are accessible through a single entry point:

```bash
python seqtools.py stats      genome.fasta
python seqtools.py orf        genome.fasta [--min 100] [--all-frames]
python seqtools.py gc_window  genome.fasta [--window 100] [--step 10]
python seqtools.py codon      genome.fasta
```

If only one `.fasta` file exists in the directory, the filename can be omitted.

---

### fasta_analyzer.py — Basic Statistics

Computes per-sequence nucleotide composition and frequency metrics.

```bash
python fasta_analyzer.py
python fasta_analyzer.py genome.fasta
```

**Output:**
- Total sequence length
- GC / AT content
- Nucleotide frequencies (A, G, C, T)
- Dinucleotide frequencies (all 16 combinations)
- GC skew, AT skew
- CpG Observed/Expected ratio

---

### orf_analyzer.py — ORF Finder

Scans all 6 reading frames (3 forward + 3 reverse complement) for open reading frames.

```bash
python orf_analyzer.py genome.fasta
python orf_analyzer.py genome.fasta --min 300
python orf_analyzer.py genome.fasta --all-frames
```

| Flag | Default | Description |
|------|---------|-------------|
| `--min` | 100 bp | Minimum ORF length filter |
| `--all-frames` | off | Show all ORFs (default: top 20 per sequence) |

**Output:** Strand, frame, start/end coordinates (1-based, forward strand), length; longest ORF summary per sequence.

---

### gc_analyzer.py — Sliding Window GC

Computes regional GC variation and saves a matplotlib figure as `gc_window.png`.

```bash
python gc_analyzer.py genome.fasta
python gc_analyzer.py genome.fasta --window 200 --step 50
```

| Flag | Default | Description |
|------|---------|-------------|
| `--window` | 100 bp | Window size |
| `--step` | 10 bp | Step size |

Red regions exceed the overall GC average; blue regions fall below.

---

### codon_analyzer.py — Codon Usage Table

Counts in-frame codons (frame +1 from position 0) and reports frequency per 1000 codons, grouped by amino acid.

```bash
python codon_analyzer.py genome.fasta
```

---

## Running Modules Standalone

Every module can also be run independently without `seqtools.py`:

```bash
python orf_analyzer.py    genome.fasta --min 300
python gc_analyzer.py     genome.fasta --window 150 --step 25
python codon_analyzer.py  genome.fasta
```

---

## Requirements

| Package | Version |
|---------|---------|
| Python | ≥ 3.9 |
| matplotlib | ≥ 3.7 |

All other dependencies (`sys`, `pathlib`, `collections`, `argparse`) are part of the Python standard library.

---

## Example Output

**stats**
```
ANALYSIS REPORT

>NM_007294 BRCA1 mRNA
Total sequence length: 5711 bp
GC content: 43.62%
AT content: 56.38%

Nucleotide Frequency
A content: 28.84% , 1647 bp
G content: 14.78% ,  844 bp
C content: 14.78% ,  844 bp
T content: 41.60% , 2376 bp

GC skew: -0.092461
AT skew: +0.141803
CpG Observed/Expected ratio: 0.441295
```

**orf**
```
ORF FINDER  (min length: 300 bp, 6 frames)

>NM_007294  [5711 bp]
  Total ORFs    : 12
  Longest ORF   : 5592 bp  (strand +, frame 1, pos 1–5592)

  #    Strand  Frame  Start      End      Length
  ------------------------------------------------
  1    +       1      1          5592     5592 bp
  2    -       2      183        1289     1107 bp
  ...
```

**codon**
```
CODON USAGE TABLE

>NM_007294  [5711 bp]
  Codon  AA     Count    Freq/1000
  -----------------------------------
  GCA    Ala    58       30.42
  GCC    Ala    71       37.24
  ...
  Total codons analyzed: 1904
```
