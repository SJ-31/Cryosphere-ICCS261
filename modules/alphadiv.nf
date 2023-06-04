process ALPHADIVERSITY {
    publishDir "$outdir/$name", mode: 'copy'
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'

    input:
    tuple val(name), path(table)
    tuple val(name), path(phylogeny)
    map(metrics)
    val(outdir)
    //
    output:
    tuple val(name), path("-*A_*.qza")
    // DistanceMatrix
    script:
    for ( metric in metrics["non-phylogenetic"] ){
    """
    qiime diversity alpha \
        --i-table $table \
        --p-metric $metric \
        --o-distance-matrix ${name}-A_${metric}.qza
    """
    }
    for ( metric in metrics["phylogenetic"] ){
    """
    qiime diversity alpha-phylogenetic \
        --i-table $table \
        --i-phylogeny $phylogeny \
        --o-distance-matrix ${name}-APhylo_${metric}.qza
    """
    }
    //
}