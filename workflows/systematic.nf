Channel.fromPath(params.sys_data)
    .map { it -> [ it.baseName, it ] }
    .set { reads_ch }

/*
 * Module imports
 */

include { FASTQC; EXTRACT_FASTQC } from '../modules/fastqc'
include { CLUSTER_BLAST } from '../modules/cluster_fasta'

workflow sys {
    EXTRACT_FASTQC(FASTQC(reads_ch, params.sys_results).dirs).set { fastqc_ch }
    CLUSTER_BLAST(fastqc_ch.oSeqs.collect(), params.sys_results)
}