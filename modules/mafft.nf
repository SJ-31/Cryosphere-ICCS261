process MAFFT {
    publishDir "$outdir"
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'

    input:
    tuple val(name), val(builder), path(unrooted_tree)
    val(outdir)
    //
    output:
    tuple val(name), path("*Rooted*")
    //
    script:
    """
    qiime phylogeny midpoint-root \
        -i-tree $unrooted_tree \
        -o-rooted-tree ${name}-${builder}_RootedTree.qza
    """
    //
}