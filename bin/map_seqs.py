#!/usr/bin/env python

import argparse
import subprocess as sp

import numpy as np
import pandas as pd
from Bio import SearchIO, SeqIO

parser = argparse.ArgumentParser(description='For every sequence in a fasta'
                                 'file, map it to its top hit from the'
                                 'specified database')
parser.add_argument('--database', nargs=1, required=True,
                    help='The name of the blast database (include path)'
                    )
parser.add_argument('--query', nargs=1, required=True,
                    help='A file in fasta format containing the query sequences'
                    )
parser.add_argument('--output', nargs=1, required=True,
                    help='Name of the output file containing the sequences'
                    'and their top hits'
                    )
args = parser.parse_args()

db = args.database[0]
to_parse = args.query[0]
output = args.output[0]
parse_type = "fasta"
blastdb: str = "data/blastdb/16sCombined"
params: str = (
    "-outfmt 5 " +
    "-max_target_seqs 10 " +
    "-perc_identity 90 "
)

seq_frame = pd.DataFrame()
sequences: np.array = np.array([])
ids: np.array = np.array([])
top_hits: np.array = np.array([])
# Search the fasta file
for query in SeqIO.parse(to_parse, parse_type):
    temp = open('temp.fasta', 'w+')
    temp.write(f'>{query.id}\n{query.seq}')
    t = temp.read()
    temp.close
    blastsearch: str = (f"blastn {params} -query temp.fasta -db {blastdb}"
                        "> temp.txt")
    sp.run(blastsearch, shell=True)
    parse = SearchIO.read('temp.xml', 'blast-xml')
    if not len(parse):
        continue
    sequences = np.append(sequences, str(query.seq))
    top_hits = np.append(top_hits, parse[0].id)
    ids = np.append(ids, query.id)

# Construct the dataframe
seq_frame['Query id'] = ids
seq_frame['Query sequence'] = sequences
seq_frame['Top hit'] = top_hits
seq_frame.to_csv(output)
