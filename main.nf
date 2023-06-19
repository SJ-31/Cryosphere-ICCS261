/*
 * Module imports
 */
include { DADA2 } from './modules/dada2'
include { EXPORT } from './modules/export'
include { CUTADAPT } from './modules/cutadapt'
include { QSCORE } from './modules/qscore'
include { MERGE } from './modules/merge'
include { PICRUST2 } from './modules/picrust2'
include { CLASSIFY_SKLEARN } from './modules/classify_sklearn'
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
fasttree_ch = Channel.empty()

def get_channel(path) {
    return Channel
            .fromPath(path)
            .map{ it -> [ it.baseName.replaceAll(/-.*/, ''), it ] }
}

Channel.fromPath('samples.tsv')
    .splitCsv(header: true, sep: "\t" )
    .map { it -> [ it.Name, it.Type, it.Trimto, it.Path ] }
    .set { samples_ch }

workflow clean_cluster {
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
        .set { otu_ch }
    EXPORT(otu_ch.freqs, otu_ch.seqs, params.outdirOTUExport)

}

workflow phylogeny {
    get_channel("$params.outdirOTU/${params.target}*Seqs*")
        .tap { all_seq_ch }
        .branch(separate_merge)
        .set { seq_ch }
    if ( params.blast ) {
        CLASSIFY_CONSENSUS_BLAST(seq_ch.the_rest.mix(seq_ch.extra),
        params.blast_args, params.refSeqs, params.refIDs,
        params.outdirClassified)
    }
    if ( params.vsearch ) {
        CLASSIFY_CONSENSUS_VSEARCH(seq_ch.the_rest.mix(seq_ch.extra),
        params.vsearch_args, params.refSeqs, params.refIDs,
        params.outdirClassified)
    }
    if ( params.sklearn ){
        CLASSIFY_SKLEARN(seq_ch.the_rest.mix(seq_ch.extra),
        params.classifier,
        params.outdirClassified
        )
    }
    // Construct phylogeny
        MAFFT(all_seq_ch, params.outdirAligned)
            .set { aligned_ch }
        if ( params.fasttree ) {
            FASTTREE(aligned_ch, params.outdirTrees)
                .tap { fasttree_ch }
        }
        if ( params.raxml ) {
            RAXML(aligned_ch, params.outdirTrees, 'GTRGAMMA')
                .tap { raxml_ch }
        }
        if ( params.iqtree ) {
            IQTREE(aligned_ch, params.outdirTrees)
                .tap { iqtree_ch }
        }
        MIDPOINTROOT(iqtree_ch.mix(raxml_ch).mix(fasttree_ch),
        params.outdirRooted)
}

workflow diversity {
    tree_ch = Channel.empty()
    get_channel("$params.outdirRooted/*")
        .map { it -> [ it[0], (it =~ /.*-(.*)_.*/)[0][1], it[1] ]}
        // Extract the tree builder automatically
        .tap { all_ch }.branch(separate_merge)
        .set { tree_ch }
    // Compute diversity
    BETADIVERSITY(tree_ch.merged.mix(tree_ch.extra), params.beta, params.outdirOTU, params.outdirDiversity)
        .transpose()
        .set { distance_matrices }
    ALPHADIVERSITY(tree_ch.the_rest.mix(tree_ch.extra),
    params.alpha, params.outdirOTU, params.outdirDiversity)

    // Analyze diversity
    PCOA(distance_matrices, params.pcoa_dimensions,
    params.outdirAnalysis)
}

workflow function_annotation {
    get_channel("$params.outdirOTUExport/*biom")
    .join(get_channel("$params.outdirOTUExport/*fasta"))
        .set { freqs_seqs_ch }
    PICRUST2(freqs_seqs_ch, params.outdirFunctions)
}
