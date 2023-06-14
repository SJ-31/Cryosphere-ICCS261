from qiime2 import Metadata, Artifact
from skbio.tree import TreeNode
import pandas as pd
metadata = pd.read_csv('./ds_metadata.tsv', sep='\t')
id_key: dict = {
    "BhM": "Bihor mountains", "BrH": "Barrow mountain high",
    "BrL": "Barrow mountain low", "CaS": "Catriona snow",
    "CeY": "Central Yakutia", "GrI": "Greenland ice",
    "EaI": "East Iceland glaciers", "CrC": "Cryoconite",
    "HaS": "Hailstone", "NzS": "New Zealand soil",
    "SvG": "Sverdrup glacier", "StR": "Storglaciaren",
    "ViS": "Villum station"
}

# Principle component analyses
pcoaBC_2D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-2D_B_braycurtis.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaJ_2D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-2D_B_jaccard.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaUU_2D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-2D_BPhylo_unweighted_unifrac_FastTree.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaNU_2D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-2D_BPhylo_weighted_normalized_unifrac_FastTree.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaWU_2D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-2D_BPhylo_weighted_unifrac_FastTree.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaBC_3D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-3D_B_braycurtis.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaJ_3D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-3D_B_jaccard.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaUU_3D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-3D_BPhylo_unweighted_unifrac_FastTree.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaNU_3D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-3D_BPhylo_weighted_normalized_unifrac_FastTree.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')
pcoaWU_3D = Artifact.load('results/8-Analysis/Merged/Merged-PCOA-3D_BPhylo_weighted_unifrac_FastTree.qza').view(Metadata).to_dataframe().merge(metadata, left_index=True, right_on='sample-id')

# Taxonomy
#   BLAST
blast = {}
for loc in id_key:
    blast_result = f'./results/3-Classified/{loc}-BLAST_ALL.qza'
    blast_top = f'./results/3-Classified/{loc}-BLAST_TopHits.qza'
    blast[loc] = [Artifact.load(r).view(Metadata).to_dataframe() for r in
                (blast_result, blast_top)]

# Phylogenetic trees, imported as Newick strings
def tree_dict(tree_type: str, ids: dict):
    t_dict = {}
    for loc in ids:
        tree = f'./results/6-RootedTrees/{loc}-{tree_type}_RootedTree.qza'
        t_dict[loc] = Artifact.load(tree).view(TreeNode).__str__()
    return t_dict

fasttree = tree_dict("FastTree", id_key)
iqtree = tree_dict("IQTREE", id_key)