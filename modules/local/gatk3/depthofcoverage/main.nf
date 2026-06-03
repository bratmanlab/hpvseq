process GATK3_DEPTHOFCOVERAGE {
    tag "${meta.id}_DEPTHOFCOVERAGE"
    label 'process_higher'

    // not confirmed
    container 'broadinstitute/gatk3:3.8-1'

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai), path(dict), path(bed)
//    tuple val(meta2), path(fasta), path(fai)
//    tuple val(meta2), path(dict)
//    tuple val(meta3), path(bed)

    output:
    path "*.sample_summary", emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals = ( bed && bed.size() > 0) ? "-L ${bed}" : ""
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

