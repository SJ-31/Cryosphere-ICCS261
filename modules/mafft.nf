process MAFFT {
    publishDir "$outdir"
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'

    input:
    tuple val(name), path(sequences)
    val(outdir)
    //
    output:
    tuple val(name), path("*Aligned*")
    //
    script:
    """
    qiime alignment mafft \
        --i-sequences $sequences \
        --o-alignment ${name}-Aligned.qza
    """
    //
}
