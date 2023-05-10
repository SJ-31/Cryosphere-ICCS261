
import subprocess as sp

import argparse
import numpy as np
import pandas as pd
from Bio import SearchIO, SeqIO

parser = argparse.ArgumentParser(description='Change maker control file'
                                 'params with'
                                 'command line arguments')
parser.add_argument('params', nargs='+', metavar='F_P',
                    help='F is the file the parameter is located in, '
                    'P the parameter (no spaces)\n'
                    'Ex: exe_RepeatMasker=/usr/bin/RepeatMasker/RepeatMasker '
                    'will change the location of the RepeatMasker executable'
                    'in the "maker_exe.ctl" file'
                    )
args = parser.parse_args()

# blastdb: str = "data/blastdb/16sCombined"
# params: str = (
#     "-outfmt 5 " +
#     "-max_target_seqs 10 " +
#     "-perc_identity 90 "
# )
# to_parse = "testblast.fasta"
# parse_type = "fasta"
#
# seq_frame = pd.DataFrame()
# sequences: np.array = np.array([])
# ids: np.array = np.array([])
# top_hits: np.array = np.array([])
# # Search the fasta file
# for query in SeqIO.parse(to_parse, parse_type):
#     temp = open('temp.fasta', 'w+')
#     temp.write(f'>{query.id}\n{query.seq}')
#     t = temp.read()
#     temp.close
#     blastsearch: str = (f"blastn {params} -query temp.fasta -db {blastdb}"
#                         "> temp.txt")
#     sp.run(blastsearch, shell=True)
#     parse = SearchIO.read('temp.txt', 'blast-xml')
#     if not len(parse):
#         continue
#     sequences = np.append(sequences, str(query.seq))
#     top_hits = np.append(top_hits, parse[0].id)
#     ids = np.append(ids, query.id)
#
# # Construct the dataframe
# seq_frame['Query id'] = ids
# seq_frame['Query sequence'] = sequences
# seq_frame['Top hit'] = top_hits
# seq_frame.to_csv('output.csv')
