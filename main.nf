// raw_ch = Channel.fromPath(params.raw_dir)

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

Channel.fromPath('samples.tsv')
    .splitCsv(header: true, sep: "\t" )
    .map { it -> [ it.Name, it.Type, it.Trimto, it.Path ] }
    .set { samples_ch }

workflow {
    // Clean fastq
    CUTADAPT(samples_ch, params.outdirClean)
        .set { trimmed_ch }
    QSCORE(trimmed_ch, params.outdirClean)
        .set { quality_ch }
    DADA2(quality_ch.seqs, params.outdirClean)
        .set { dd_ch }
    dd_ch.table.flatten().filter{ !(it =~ /$params.ignore/) }.filter( ~/.*qza/ )
        .collect().map { it -> ['Merged', it ]}.set { mfreqs }
    dd_ch.seqs.flatten().filter{ !(it =~ /$params.ignore/) }.filter( ~/.*qza/ )
        .collect().map { it -> ['Merged', it ]}.set { mseqs }
    MERGE(mfreqs, mseqs, params.outdirClean)
        .set { merged_ch }
    VSEARCH_CLUSTER_DENOVO(dd_ch.table.mix(merged_ch.table), dd_ch.seqs.mix(merged_ch.seqs),
    params.outdirOTU, "0.99")
        .set { otu_ch }
    otu_ch.view()
    // otu_ch.branch {

    // }

    // Classify taxonomy
    // CLASSIFY_CONSENSUS_BLAST(otu_ch.seqs, params.blast_args,
    // params.refSeqs, params.refIDs, params.outdirClassified)
    //     .set { blast_ch }
    // CLASSIFY_CONSENSUS_VSEARCH(otu_ch.seqs, params.vsearch_args, params.refSeqs, params.refIDs, params.outdirClassified)
    //     .set { vsearch_ch } // Not done yet
    // Merge afterwards or else you die

    // BETADIVERSITY()
    // ALPHADIVERSITY()
    // BETAPHYLO()
    // ALPHAPHYLO()

    // // Construct phylogeny
    // MAFFT(otu_ch.seqs, params.outdirAligned)
    //     .set { aligned_ch }
    // FASTTREE(aligned_ch, params.outdirTrees)
    //     .set { fasttree_ch }
    // RAXML_RAPID_BOOTSTRAP(aligned_ch, params.outdirTrees,
    // '1000', 'GTRGAMMA')
    //     .set { raxml_ch }
    // IQTREE_ULTRAFAST_BOOTSTRAP(aligned_ch, params.outdirTrees,
    // '1000')
    //     .set { iqtree_ch }
    // MIDPOINTROOT(iqtree_ch.mix(raxml_ch).mix(fasttree_ch),
    // params.outdirRooted)
    //     .set { rooted_ch }
}


