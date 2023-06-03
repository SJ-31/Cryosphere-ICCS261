process CLASSIFY_CONSENSUS_VSEARCH {
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'
    publishDir "$outdir", mode: "copy"

    input:
    tuple val(name), path(otus)
    val(args)
    val(refSeqs)
    val(refIDs)
    val(outdir)
    //
    output:
    tuple val(name), val('Vsearch'), path("*All*"), emit: all
    // FeatureData[Taxonomy]
    tuple val(name), val('Vsearch'), path("*Top*"), emit: top
    // FeatureData[BLAST6]
    script:
    """
    qiime feature-classifier classify-consensus-vsearch \
        ${args} \
        --i-query $otus \
        --i-reference-reads $refSeqs \
        --i-reference-taxonomy $refIDs \
        --o-classification ${name}-Vsearch_All.qza  \
        --o-search-results ${name}-Vsearch_TopHits.qza
    """
    //
}