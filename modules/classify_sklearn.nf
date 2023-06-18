process CLASSIFY_SKLEARN {
    publishDir "$outdir", mode: "copy"
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'

    input:
    tuple val(name), path(otus)
    val(classifier)
    val(outdir)
    //
    output:
    tuple val(name), path("*Sklearn.qza"), emit: freqs // FeatureData[Taxonomy]
    //
    script:
    """
    qiime feature-classifier classify-sklearn \
        --i-reads $otus \
        --i-classifier $classifier \
        --o-classification ${name}-Sklearn.qza
    """
    //
}