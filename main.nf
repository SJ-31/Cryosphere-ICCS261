// raw_ch = Channel.fromPath(params.raw_dir)

/*
 * Main workflow
 */

include { sys } from './workflows/systematic'
include { IQTREE-ULTRAFAST-BOOTSTRAP } from './modules/iqtree_ultra'
include { sys } from './workflows/systematic'
include { sys } from './workflows/systematic'

workflow {
    if ( params.run_sys )
        sys()

}
