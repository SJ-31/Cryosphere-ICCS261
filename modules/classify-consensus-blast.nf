process CLASSIFY-CONSENSUS-BLAST {
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'
    publishDir "$outdir", mode: 'symlink'

    input:
    tuple val(name) path(otus)
    val(args)
    val(refSeqs)
    val(refIDs)
    val(outdir)
    //
    output:
    tuple val(name), val('BLAST'), path("*All*"), emit: all
    tuple val(name), val('BLAST'), path("*Top*"), emit: top
    //
    script:
    """
    qiime feature-classifier classify-consensus-blast \
        ${args} \
        --i-query $otus \
        --i-reference-reads $refSeqs \
        --i-reference-taxonomy $refIDs \
        --o-classification ${name}-BLAST_All.qza  \
        --o-search-results ${name}-BLAST_TopHits.qza
    """
    //
}