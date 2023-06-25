process ALPHADIVERSITY {
    publishDir "$outdir/$name", mode: 'copy'
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'

    input:
    tuple val(name), val(builder), path(phylogeny)
    val(metrics)
    val(otus)
    val(outdir)
    //
    output:
    tuple val(name), path("*-A*.qza")
    // DistanceMatrix
    shell:
    table = "${otus}/${name}-otuFreqs.qza"
    '''
    for metric in !{metrics["non-phylogenetic"].join(" ")}
        do
        qiime diversity alpha \
            --i-table !{table} \
            --p-metric $metric \
            --o-alpha-diversity !{name}-A_${metric}.qza
        done
    for metric in !{metrics["phylogenetic"].join(" ")}
        do
        qiime diversity alpha-phylogenetic \
            --i-table !{table} \
            --p-metric $metric \
            --i-phylogeny !{phylogeny} \
            --o-alpha-diversity !{name}-APhylo_${metric}.qza
        done
    '''
    //
}