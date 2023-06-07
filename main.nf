/*
 * Module imports
 */
include { DADA2 } from './modules/dada2'
include { CUTADAPT } from './modules/cutadapt'
include { QSCORE } from './modules/qscore'
include { MERGE } from './modules/merge'
include { BETADIVERSITY } from './modules/betadiv.nf'
include { ALPHADIVERSITY } from './modules/alphadiv.nf'
include { VSEARCH_CLUSTER_DENOVO } from './modules/vsearch_cluster_de-novo.nf'
include { PCOA } from './modules/pcoa'
include { CLASSIFY_CONSENSUS_BLAST } from './modules/classify-consensus-blast.nf'
include { CLASSIFY_CONSENSUS_VSEARCH } from './modules/classify-consensus-vsearch.nf'
include { MAFFT } from './modules/mafft.nf'
include { RAXML_RAPID_BOOTSTRAP; RAXML } from './modules/raxml.nf'
include { IQTREE_ULTRAFAST_BOOTSTRAP; IQTREE } from './modules/iqtree.nf'
include { FASTTREE } from './modules/fasttree.nf'
include { MIDPOINTROOT } from './modules/midpoint-root.nf'

/*
 * Main workflow
 */

def separate_merge = branchCriteria {
    merged: it[0] =~ /Merged/
    extra: it[0] =~ /$params.ignore/
    the_rest: !(it[0] =~ /Merged/)
}

iqtree_ch = Channel.empty()
raxml_ch = Channel.empty()

Channel.fromPath('samples.tsv')
    .splitCsv(header: true, sep: "\t" )
    .map { it -> [ it.Name, it.Type, it.Trimto, it.Path ] }
    .set { samples_ch }

workflow {
    // Clean and merge fastq
    CUTADAPT(samples_ch, params.outdirClean)
        .set { trimmed_ch }
    QSCORE(trimmed_ch, params.outdirClean)
        .set { quality_ch }
    DADA2(quality_ch.seqs, params.outdirClean)
        .set { dd_ch }
    dd_ch.table.flatten().filter{ !(it =~ /$params.ignore/) }.filter( ~/.*qza/ ).collect().map { it -> ['Merged', it ]}
        .set { mfreqs }
    dd_ch.seqs.flatten().filter{ !(it =~ /$params.ignore/) }.filter( ~/.*qza/ ).collect().map { it -> ['Merged', it ]}
        .set { mseqs }
    MERGE(mfreqs, mseqs, params.outdirClean)
        .set { merged_ch }
    VSEARCH_CLUSTER_DENOVO(dd_ch.table.mix(merged_ch.table), dd_ch.seqs.mix(merged_ch.seqs),params.outdirOTU, "0.99")
        .set { all_otu }
    all_otu.seqs.branch(separate_merge)
        .set { otu_seqs }

    // Classify taxonomy
    CLASSIFY_CONSENSUS_BLAST(otu_seqs.the_rest, params.blast_args,
    params.refSeqs, params.refIDs, params.outdirClassified)
        .set { blast_ch }
    // CLASSIFY_CONSENSUS_VSEARCH(otu_seqs.other, params.vsearch_args, params.refSeqs, params.refIDs, params.outdirClassified)
    //     .set { vsearch_ch } // Not done yet

    // Construct phylogeny
    MAFFT(otu_seqs.merged.mix(otu_seqs.extra).mix(otu_seqs.the_rest),
    params.outdirAligned).branch(separate_merge)
        .set { aligned_ch }
    FASTTREE(aligned_ch.merged.mix(aligned_ch.the_rest).mix(aligned_ch.extra),
    params.outdirTrees)
        .set { fasttree_ch }
    // RAXML_RAPID_BOOTSTRAP(aligned_ch, params.outdirTrees,
    // '1000', 'GTRGAMMA')
    RAXML(aligned_ch.the_rest, params.outdirTrees, 'GTRGAMMA')
        .tap { raxml_ch }
    // IQTREE_ULTRAFAST_BOOTSTRAP(aligned_ch, params.outdirTrees,
    // '100')
    //     .set { iqtree_ch }
    IQTREE(aligned_ch.the_rest, params.outdirTrees)
        .tap { iqtree_ch }
    MIDPOINTROOT(iqtree_ch.mix(raxml_ch).mix(fasttree_ch),
    params.outdirRooted).branch {
        fasttree: it[1] =~ /FastTree/
        iqtree: it[1] =~ /IQTREE/
        raxml: it[1] =~ /RAxML/
        all: true
    }.set { rooted_ch }
    all_otu.freqs.join(rooted_ch.fasttree)
        .branch(separate_merge)
        .set { freqs_trees }
        // todo: You could change the tree to base the phylogeny on...

    // Compute diversity
    BETADIVERSITY(freqs_trees.merged.mix(freqs_trees.extra), params.beta, params.outdirDiversity)
        .transpose()
        .set { distance_matrices }
    ALPHADIVERSITY(freqs_trees.the_rest.mix(freqs_trees.extra),
    params.alpha, params.outdirDiversity)

    // Analyze diversity
    PCOA(distance_matrices, params.pcoa_dimensions,
    params.outdirAnalysis)


}


