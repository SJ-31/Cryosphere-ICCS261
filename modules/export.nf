process EXPORT {
    publishDir "$outdir", mode: 'copy'
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'

    input:
    tuple val(name), path(freqs)
    tuple val(name), path(seqs)
    val(outdir)
    //
    output:
    tuple val(name), path("*biom")
    tuple val(name), path("*fasta")
    //
    script:
    """
    qiime tools export \
    --input-path $freqs \
    --output-path temp_OTU
    cp temp_OTU/feature-table.biom ./${name}-otuFreqs.biom
    qiime tools export \
    --input-path $seqs \
    --output-path temp_seqs
    cp temp_seqs/dna-sequences.fasta ./${name}-otuSeqs.fasta
    """
    //
}