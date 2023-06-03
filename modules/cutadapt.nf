process CUTADAPT {
    conda '/home/sc31/Bio_SDD/miniconda3/envs/qiime2-2023.2'
    publishDir "$outdir", mode: "copy"

    input:
    tuple val(name), val(type), val(trunc), val(artifact_path)
    val(outdir)
    //
    output:
    tuple val(name), val(type), val(trunc), path("*-Trimmed.qza")
    // Either SampleData[SequencesWithQuality]
    // or SampleData[PairedEndSequencesWithQuality]
    script:
    if ( type == "paired" ) {
    """
    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences $artifact \
        --o-trimmed-sequences ${name}-Trimmed.qza
    """
    }
    else if (type == "single") {
    """
    qiime cutadapt trim-single \
        --i-demultiplexed-sequences $artifact
        --o-trimmed-sequences ${name}-Trimmed.qza
    """
    }
    //
}