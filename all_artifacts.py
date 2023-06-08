from qiime2 import Artifact
from qiime2 import Metadata
import pandas as pd
metadata = pd.read_csv('./ds_metadata.tsv', sep='\t')

# Principle component analyses
pcoaBC = Artifact.load('results/8-Analysis/Merged/Merged-PCOA_B_braycurtis.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaJ = Artifact.load('results/8-Analysis/Merged/Merged-PCOA_B_jaccard.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaUU = Artifact.load('results/8-Analysis/Merged/Merged-PCOA_BPhylo_unweighted_unifrac.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaNU = Artifact.load('results/8-Analysis/Merged/Merged-PCOA_BPhylo_weighted_normalized_unifrac.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaWU = Artifact.load('results/8-Analysis/Merged/Merged-PCOA_BPhylo_weighted_unifrac.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')