process PCOA  {
    publishDir "$outdir/$name", mode: 'copy'
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'

    input:
    tuple val(name), path(matrix)
    each dimensions
    val(outdir)
    //
    output:
    tuple (val(name),
    path("${name}-PCOA-${dimensions}D_${metric}.qza"))
    //
    script:
    metric = matrix.baseName.replaceAll(/.*-/, '')
    """
    qiime diversity pcoa \
        --i-distance-matrix $matrix \
        --p-number-of-dimensions $dimensions \
        --o-pcoa ${name}-PCOA-${dimensions}D_${metric}.qza
    """
    //
}