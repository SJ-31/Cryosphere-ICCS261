process FASTTREE {
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'
    publishDir "$outdir", mode: 'symlink'

    input:
    tuple val(name) path(aligned_sequence)
    val(outdir)
    //
    output:
    tuple val(name), val('FastTree'), path("*Unrooted*")
    //
    script:
    """
    qiime phylogeny fasttree \
        --i-alignment $aligned_sequence \
        --o-tree ${name}-FastTree_UnrootedTree.qza
    """
    //
}