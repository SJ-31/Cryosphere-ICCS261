process DADA2 {
    tag "Trimming $name"
    publishDir "$outdir", mode: "copy", pattern: "*.qza"
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'

    input:
    tuple val(name), val(type), val(trunc), path(artifact)
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
    if  ( type == "paired" )
        """
        qiime dada2 denoise-paired \
        --i-demultiplexed-seqs $artifact \
        --p-trunc-len-r $trunc \
        --p-trunc-len-f $trunc \
        --o-table ${name}-denoisedTable.qza \
        --o-representative-sequences ${name}-denoisedSeqs.qza \
        --o-denoising-stats ${name}-denoised.qza
        """
    else if ( type == "single")
        """
        qiime dada2 denoise-single \
        --i-demultiplexed-seqs $artifact \
        --p-trunc-len $trunc \
        --o-table ${name}-denoisedTable.qza \
        --o-representative-sequences ${name}-denoisedSeqs.qza \
        --o-denoising-stats ${name}-denoisedStats.qza
        """
    //
}