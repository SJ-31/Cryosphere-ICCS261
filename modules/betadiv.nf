process BETADIVERSITY {
    publishDir "$outdir/$name", mode: 'copy'
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'

    input:
    tuple (val(name), val(builder), path(phylogeny))
    val(metrics)
    val(otus)
    val(outdir)
    //
    output:
    tuple val(name), path("*-B*.qza")
    // DistanceMatrix
    shell:
    table = "${otus}/${name}-otuFreqs.qza"
    '''
    for metric in !{metrics["non-phylogenetic"].join(" ")}
        do
        qiime diversity beta \
            --i-table !{table} \
            --p-metric $metric \
            --o-distance-matrix \
            !{name}-B_${metric}.qza
        done
    for metric in !{metrics["phylogenetic"].join(" ")}
        do
        qiime diversity beta-phylogenetic \
            --i-table !{table} \
            --p-metric $metric \
            --i-phylogeny !{phylogeny} \
            --o-distance-matrix \
            !{name}-BPhylo_${metric}_!{builder}.qza
        done
    '''
    // You can't parallelize unifrac or else it dies
}