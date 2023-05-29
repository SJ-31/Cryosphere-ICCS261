process FASTQC {
    publishDir "$outdir", mode: 'copy', pattern: "*.html"

    input:
    tuple val(name), path(reads)
    val(outdir)
    //
    output:
    path("*_fastqc"), emit: dirs
    path("*.html"), emit: html
    //
    script:
    """
    fastqc $reads \
    --extract
    """
    //
}

process EXTRACT_FASTQC {
    publishDir "$outdir", mode: 'copy', pattern: "*fasta"

    input:
    path(fastqc)
    val(outdir)
    //
    output:
    path("*Adapter_Content*"), emit: aContent
    path("*fasta"), emit: oSeqs
    //
    script:
    def name = fastqc.baseName.replaceAll(/_.*/, '')
    """
    parse_fastqc.py $fastqc
    oRseqs2fasta.py ${name}_Overrepresented_sequences.txt
    """
    //
}