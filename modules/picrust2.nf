process PICRUST2 {
    publishDir "$outdir", mode: 'copy'
    conda '/home/sc31/Bio_SDD/miniconda3/envs/picrust2'

    input:
    tuple val(name), path(biom), path(fasta)
    val(outdir)
    //
    output:
    tuple val(name), path("*PICRUST2")
    //
    script:
    """
    picrust2_pipeline.py \
        -s $fasta \
        -i $biom \
        -p 1 \
        --remove_intermediate \
        -o ${name}_PICRUST2
    gunzip ${name}_PICRUST2/*/*
    """
    //
}