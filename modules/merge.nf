process MERGE {
    publishDir "$outdir", mode: "copy"
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'

    input:
    tuple val(name), path(tables)
    tuple val(name), path(seqs)
    val(outdir)
    //
    output:
    path("mergedSeqs.qza"), emit: seqs
    // FeatureTable[Sequences]
    tuple val(name), path("mergedFreqs.qza"), emit: table
    // FeatureTable[Frequency]
    script:
    """
    qiime feature-table merge-seqs \
        --i-data $seqs \
        --o-merged-data mergedSeqs.qza
    qiime feature-table merge \
        --i-tables $tables \
        --o-merged-table mergedFreqs.qza
    """
    //
}