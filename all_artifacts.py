from qiime2 import Artifact
from qiime2 import Metadata
import pandas as pd
metadata = pd.read_csv('./ds_metadata.tsv', sep='\t')

# Principle component analyses
pcoaBC_2D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-2D_B_braycurtis.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaJ_2D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-2D_B_jaccard.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaUU_2D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-2D_BPhylo_unweighted_unifrac.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaNU_2D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-2D_BPhylo_weighted_normalized_unifrac.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaWU_2D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-2D_BPhylo_weighted_unifrac.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaBC_3D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-3D_B_braycurtis.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaJ_3D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-3D_B_jaccard.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaUU_3D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-3D_BPhylo_unweighted_unifrac.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaNU_3D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-3D_BPhylo_weighted_normalized_unifrac.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaWU_3D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-3D_BPhylo_weighted_unifrac.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')