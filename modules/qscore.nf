process QSCORE {
    publishDir "$outdir", pattern: "*-QSeqs.qza", mode: "copy"
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'

    input:
    tuple val(name), val(type), val(trunc), path(artifact)
    val(outdir)
    //
    output:
    tuple val(name), val('single'), val(trunc), path("*-QSeqs.qza"), emit: seqs
    // SampleData[SequencesWithQuality]
    tuple val(name), path("*-QStats.qza"), emit: stats
    // QualityFilterStats
    script:
    """
    qiime quality-filter q-score  \
        --i-demux $artifact  \
        --p-min-quality 4 \
        --p-quality-window 3 \
        --o-filtered-sequences ${name}-QSeqs.qza \
        --o-filter-stats ${name}-QStats.qza
    """
    //
}