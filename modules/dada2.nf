process DADA2 {
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'

    input:
    path()
    //
    output:
    path()
    //
    script:
    """
    qiime dada2 denoise-paired \
    --i-demultiplexed-seqs <>.qza \
    --p-trim-left \
    --p-trim -right \
    --p-trunc-len-r \
    --p-trunc-len-f \
    --o-repesentative-sequences  \
    --o-table \
    --o-denoising-stats
    """
    //
}