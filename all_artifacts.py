from qiime2 import Metadata, Artifact
from skbio.tree import TreeNode
from skbio.stats.distance import DistanceMatrix
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

def to_distance(artifact_path: str):
    return Artifact.load(artifact_path).view(DistanceMatrix)

def to_df(artifact_path: str):
    return Artifact.load(artifact_path).view(Metadata).to_dataframe()

# Principle component analyses
beta_metrics: dict = { "bc": "B_braycurtis", "ja": "B_jaccard",
                    "uu": "BPhylo_unweighted_unifrac",
                    "wn": "BPhylo_weighted_normalized_unifrac",
                    "wu": "BPhylo_weighted_unifrac"}
alpha_metrics: dict = {"sh": "shannon", "pi": "pielou_e",
                    "si": "simpson", "se": "simpson_e"}
alpha = {}
for loc in id_key:
    path = "results/7-Diversity"
    alpha[loc] = {}
    for name, metric in alpha_metrics.items():
        alpha[loc][name] = Artifact.load(f"{path}/{loc}/{loc}-A_{metric}.qza").view(Metadata).to_dataframe()
beta_dm = {}
pcoa2D = {}
pcoa3D = {}
for name, metric in beta_metrics.items():
        dm_path = "results/7-Diversity/Merged/Merged-"
        pcoa_path = "results/8-Analysis/Merged/Merged-PCOA"
        if not name in {"bc", "ja"}:
            ext = "_FastTree"
        else:
            ext = ""
        beta_dm[name] = to_distance(f"{dm_path}{metric}{ext}.qza")
        pcoa2D[name] = to_df(f"{pcoa_path}-2D_{metric}{ext}.qza").merge(metadata, left_index=True, right_on="sample-id")
        pcoa3D[name] = to_df(f"{pcoa_path}-3D_{metric}{ext}.qza").merge(metadata, left_index=True, right_on="sample-id")

# Taxonomy
#   BLAST
blast = {}
for loc in id_key:
    blast_result = f'./results/3-Classified/{loc}-BLAST_All.qza'
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