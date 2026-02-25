# Copyright (C) 2026 Berkem Toprak Elmacı
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License.

import sys
from pathlib import Path
from collections import Counter


class FastaAnalyzer:

    def __init__(self, filename):
        filelist=list(Path().glob("*.fasta"))

        with open(filename, "r") as f:
            content=f.read()

        sequence_number=content.count(">")

        sequence_dict={
            f"sequence{i}":None
            for i in range(1,sequence_number+1)
            }
        report=["ANALYSIS REPORT"]

        entries = content.split(">")[1:] 
        for i, entry in enumerate(entries, 1):
            lines = entry.splitlines()
            sequence_id = ">" + lines[0]
            sequence_dict[f"sequence{i}"] = "".join(lines[1:]).upper()
            report.append(sequence_id)

        for i in range(1, sequence_number+1):
            seq = sequence_dict[f"sequence{i}"]

            mono = Counter(seq)
            A, G, C, T = mono["A"], mono["G"], mono["C"], mono["T"]
            
            di = Counter(seq[j:j+2] for j in range(len(seq)-1))
            GC, CG = di["GC"], di["CG"]
            CA, CT = di["CA"], di["CT"]
            GA, GT = di["GA"], di["GT"]
            AG, AC = di["AG"], di["AC"]
            AT, TA = di["AT"], di["TA"]
            TC, TG = di["TC"], di["TG"]
            TT, AA = di["TT"], di["AA"]
            GG, CC = di["GG"], di["CC"]

            N = A+G+C+T
            GC_content = (G+C)/N if N != 0 else 0
            AT_content = (A+T)/N if N != 0 else 0
            A_ratio = A/N if N != 0 else 0
            G_ratio = G/N if N != 0 else 0
            C_ratio = C/N if N != 0 else 0
            T_ratio = T/N if N != 0 else 0
            GC_skew = (G-C)/(G+C) if (G+C) != 0 else 0
            AT_skew = (A-T)/(A+T) if (A+T) != 0 else 0

            CpG = (CG*N)/(C*G) if (C!=0 and G!=0) else 0
            
            report[i] = report[i]+f"""
Total sequence length:{N} bp
GC content: {GC_content * 100:.2f}%
AT content: {AT_content * 100:.2f}%

Nucleotide Frequency
A content: {A_ratio * 100:.2f}% , {A} bp
G content: {G_ratio * 100:.2f}% , {G} bp
C content: {C_ratio * 100:.2f}% , {C} bp
T content: {T_ratio * 100:.2f}% , {T} bp

Dinucleotide Frequency
CG: {CG/(N-1):.6f}
GC: {GC/(N-1):.6f}
CA: {CA/(N-1):.6f}
CT: {CT/(N-1):.6f}
GA: {GA/(N-1):.6f}
GT: {GT/(N-1):.6f}
AG: {AG/(N-1):.6f}
AC: {AC/(N-1):.6f}
AT: {AT/(N-1):.6f}
TA: {TA/(N-1):.6f}
TC: {TC/(N-1):.6f}
TG: {TG/(N-1):.6f}
TT: {TT/(N-1):.6f}
AA: {AA/(N-1):.6f}
GG: {GG/(N-1):.6f}
CC: {CC/(N-1):.6f}

GC skew: {GC_skew:.6f}
AT skew: {AT_skew:.6f}

CpG Observed/Expected ratio: {CpG:.6f}

"""

        for i in range(len(report)):
            print(report[i]+"\n")


if __name__ == "__main__":
    filelist = list(Path().glob("*.fasta"))

    if len(sys.argv)==1 and len(filelist)==1:
        filename=filelist[0]
    elif len(sys.argv)==2 and sys.argv[1].endswith(".fasta"):
        filename=sys.argv[1]
    else:
        print("""ERROR: There are multiple .fasta files in the directory, the program cannot run automatically.
Please try `python fasta_analyzer.py sequencefile.fasta` in bash""")
        sys.exit()

    FastaAnalyzer(filename)
