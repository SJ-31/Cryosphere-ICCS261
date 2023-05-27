// raw_ch = Channel.fromPath(params.raw_dir)

/*
 * Main workflow
 */

include { sys } from './workflows/systematic'

workflow {
    if ( params.run_sys )
        sys()
}
