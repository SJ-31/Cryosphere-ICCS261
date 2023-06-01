// raw_ch = Channel.fromPath(params.raw_dir)

/*
 * Module imports
 */

include { VSEARCH_CLUSTER_DENOVO } from './modules/vsearch_cluster_de-novo.nf'
include { CLASSIFY_CONSENSUS_BLAST } from './modules/classify-consensus-blast.nf'
include { CLASSIFY_CONSENSUS_VSEARCH } from './modules/classify-consensus-vsearch.nf'
include { MAFFT } from './modules/mafft.nf'
include { RAXML_RAPID_BOOTSTRAP } from './modules/raxml_bootstrap.nf'
include { IQTREE_ULTRAFAST_BOOTSTRAP } from './modules/iqtree_ultrafast_bootstrap.nf'
include { FASTTREE } from './modules/fasttree.nf'
include { MIDPOINTROOT } from './modules/midpoint-root.nf'

/*
 * Main workflow
 */

Channel.fromPath(params.denoisedTable)
    .map{ it -> ["$params.name", it ] }.set{ denoisedTable_ch }
Channel.fromPath(params.denoisedSeqs)
    .set { denoisedSeqs_ch }


workflow {
    VSEARCH_CLUSTER_DENOVO(denoisedTable_ch, denoisedSeqs_ch,
    params.outdirOTU, "0.99")
        .set { otu_ch }
    //
    CLASSIFY_CONSENSUS_BLAST(otu_ch.seqs, params.blast_args,
    params.refSeqs, params.refIDs, params.outdirClassified)
        .set { blast_ch }
    CLASSIFY_CONSENSUS_VSEARCH(otu_ch.seqs, params.vsearch_args, params.refSeqs, params.refIDs, params.outdirClassified)
        .set { vsearch_ch }

    // Phylogeny
    MAFFT(otu_ch.seqs, params.outdirAligned)
        .set { aligned_ch }
    FASTTREE(aligned_ch, params.outdirTrees)
        .set { fasttree_ch }
    RAXML_RAPID_BOOTSTRAP(aligned_ch, params.outdirTrees,
    '1000', 'GTRGAMMA')
        .set { raxml_ch }
    IQTREE_ULTRAFAST_BOOTSTRAP(aligned_ch, params.outdirTrees,
    '1000')
        .set { iqtree_ch }
    MIDPOINTROOT(iqtree_ch.mix(raxml_ch).mix(fasttree_ch),
    params.outdirRooted)
        .set { rooted_ch }
}


