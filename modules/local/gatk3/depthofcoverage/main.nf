process GATK3_DEPTHOFCOVERAGE {
    tag "${meta.id}_${meta.ref}${meta.type}_${meta.ref2}${meta.type2}_${meta.consensus}"
    label 'process_medium'

    // not confirmed
    container 'broadinstitute/gatk3:3.8-1'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta2), path(dict)
    tuple val(meta3), path(bed)

    output:
    path "*.sample_summary", emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals = bed ? "-L ${bed}" : ""
    """
    java -jar \$gatk_dir/GenomeAnalysisTK.jar \\
        -T DepthOfCoverage \\
        -R $fasta \\
        -I $bam \\
        -o ${prefix}.coverage \\
        ${args} \\
        ${intervals}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prfiex}.coverage.sample_summary
    """
}

