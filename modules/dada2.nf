process DADA2_PAIRED {
    publishDir "$outdir", mode: 'symlink', pattern: "*-denoised.qza"
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'

    input:
    tuple val(name), path(artifact)
    val(mode)
    val(outdir)
    //
    output:
    tuple val(name), path("*Table*"), emit: table
    // FeatureTable[Frequency]
    path("*Stats*"), emit: stats
    // SampleData[DADA2Stats]
    path("*Seqs*"), emit: seqs
    // FeatureData[Sequence]
    script:
    if  (mode == 'paired' )
        """
        qiime dada2 denoise-paired \
        --i-demultiplexed-seqs  \
        --p-trim-left \
        --p-trim -right \
        --p-trunc-len-r \
        --p-trunc-len-f \
        --o-table ${name}-denoisedSeqs.qza \
        --o-repesentative-sequences ${name}-denoisedTable.qza \
        --o-denoising-stats ${name}-denosiedStats.qza
        """
    else if ( mode == 'pyro')
    //
}