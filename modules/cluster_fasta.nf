process CLUSTER_BLAST {
    publishDir "$outdir", mode: 'copy'

    input:
    path(files)
    val(outdir)
    //
    output:
    path("*.txt")
    //
    script:
    """
    cat $files > all.fasta
    cd-hit -i all.fasta -o clustered.fasta
    blastn -query clustered.fasta -db $projectDir/data/blastdb/16SAll > oRseqs_hits.txt
    """
    //
}