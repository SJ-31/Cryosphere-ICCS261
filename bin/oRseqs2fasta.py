#!/usr/bin/env python

import sys
import pandas as pd

filepath = sys.argv[1]
print(filepath)
name = filepath[: filepath.find('_')]
seqs = pd.read_csv(f'{sys.argv[1]}', sep='\t')
fasta: list = []
for index, row in seqs.iterrows():
    fasta.append(f'>{name} Count:{row["Count"]} Percent:{row["Percentage"]}\n')
    fasta.append(f'{row["#Sequence"]}\n')
with open(f'{name}_oRseqs.fasta', 'w') as f:
    f.write(''.join(fasta))
