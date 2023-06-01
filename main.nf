// raw_ch = Channel.fromPath(params.raw_dir)

/*
 * Module imports
 */

include { VSEARCH_CLUSTER_DE-NOVO } from './modules/vsearch_cluster_de-novo.nf'
include { CLASSIFY-CONSENSUS-BLAST } from './modules/classify-consensus-blast.nf'
include { CLASSIFY-CONSENSUS-VSEARCH } from './modules/classify-consensus-vsearch.nf'
include { MAFFT } from './modules/mafft.nf'
include { RAXML-RAPID-BOOTSTRAP } from './modules/raxml_boostrap.nf'
include { IQTREE-ULTRAFAST-BOOTSTRAP } from './modules/iqtree_ultrafast_bootstrap.nf'
include { FASTTREE } from './modules/fasttree.nf'
include { MIDPOINT-ROOT } from './modules/midpoint-root.nf'

/*
 * Main workflow
 */

Channel.fromPath()

workflow {
    VSEARCH_CLUSTER_DE-NOVO(table_ch, seq_ch,
        params.outdirOTU, "0.99").set { otu_ch }
    //
    CLASSIFY-CONSENSUS-BLAST(otu_ch, params.blast_args, params.refSeqs,
        params.refIDs, params.outdirClassified).set { blast_ch }
    CLASSIFY-CONSENSUS-VSEARCH(otu_ch, params.vsearch_args, params.refSeqs,
        params.refIDs, params.outdirClassified).set { vsearch_ch }

    // Phylogeny
    MAFFT(otu_ch).set { aligned_ch }
    FASTTREE(aligned_ch, params.outdirTrees).set { fasttree_ch }
    RAXML-RAPID-BOOTSTRAP(aligned_ch, params.outdirTrees,
        '1000', 'GTRGAMMA').set { raxml_ch }
    IQTREE-ULTRAFAST-BOOTSTRAP(aligned_ch, params.outdirTrees,
        '1000').set { iqtree_ch }
    MIDPOINT-ROOT(iqtree_ch.mix(raxml_ch).mix(fasttree_ch),
        params.outdirRooted).set { rooted_ch }
}


