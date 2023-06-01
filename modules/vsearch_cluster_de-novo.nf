process VSEARCH_CLUSTER_DE-NOVO {
    tag "Clustering $name, Table: $table, Seqs: $seqs"
    publishDir "$outdir", pattern: "*otu*"
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'

    input:
    tuple val(name), path(table)
    path(seqs)
    val(outdir)
    val(percent_identity)
    //
    output:
    tuple val(name), path("*Freqs.qza"), emit: freqs
    path("*Seqs.qza"), emit: seqs
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