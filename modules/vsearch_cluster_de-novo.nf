process VSEARCH_CLUSTER_DENOVO {
    tag "Clustering $name, Table: $table, Seqs: $seqs"
    publishDir "$outdir", pattern: "*otu*", mode: "copy"
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'
    memory '2 GB'

    input:
    tuple val(name), path(table)
    path(seqs)
    val(outdir)
    val(percent_identity)
    //
    output:
    tuple val(name), path("*Freqs.qza"), emit: freqs // FeatureData[Frequency]
    tuple val(name), path("*Seqs.qza"), emit: seqs // FeatureData[Sequence]
    //
    script:
    """
    qiime vsearch cluster-features-de-novo \
        --i-sequences $seqs  \
        --i-table $table \
        --p-perc-identity $percent_identity \
        --o-clustered-table $name-otuFreqs.qza \
        --o-clustered-sequences $name-otuSeqs.qza
    """
    //
}