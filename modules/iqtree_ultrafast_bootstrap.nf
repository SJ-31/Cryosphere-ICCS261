process IQTREE_ULTRAFAST_BOOTSTRAP {
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'
    publishDir "$outdir", mode: "copy"

    input:
    tuple val(name), path(aligned_sequence)
    val(outdir)
    val(replicates)
    //
    output:
    tuple val(name), val('IQTREE'), path("*UnrootedTree.qza")
    // Phylogeny[Unrooted]
    script:
    """
    qiime phylogeny iqtree-ultrafast-bootstrap \
        --i-alignment $aligned_sequence \
        --p-bootstrap-replicates $replicates \
        --o-tree ${name}-IQTREE_UnrootedTree.qza
    """
    //
}