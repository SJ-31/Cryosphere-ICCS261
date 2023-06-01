process RAXML-RAPID-BOOTSTRAP {
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'
    publishDir "$outdir", mode: 'symlink'

    input:
    tuple val(name) path(aligned_sequence)
    val(outdir)
    val(replicates)
    val(model)
    //
    output:
    tuple val(name), val('RAxML'), path("*Unrooted*)
    //
    script:
    """
    qiime phylogeny raxml-rapid-bootstrap \
    --i-alignment $aligned_sequence \
    --p-bootstrap-replicates $replicates \
    --p-substitution-model $model \
    --o-tree ${name}-RAxML_UnrootedTree.qza
    """
    // Default model is GTRGAMMA
    // DEfault bootstrap replicate number is 100
    //
}