process INTEGRATION_SEARCHPV {
    tag "$meta.id"
    label 'process_high'
    publishDir "results", mode: 'copy', saveAs: { filename ->
        if (filename == "HPVfusionPointContig.txt") {
            return "${meta.id}_${filename}" // Renames the folder
        } else if (filename.startsWith("{meta.id}")) {
            return filename                      // Keeps original name for anything else
        }
    }

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta_hg), path(fai_hg)
    tuple val(meta3), path(fasta_virus), path(fai_virus)

    output:
    tuple val(meta), path("${meta.id}__${meta3.virus}"), emit: searchpv_dir
    path("${meta.id}__${meta3.virus}/call_fusion_virus/HPVfusionPointContig.txt"), emit: HPVfusionPointContig

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-mem ${task.memory.giga} -nt ${task.cpus}'
    """
    searcHPV -fastq1 ${reads[0]} \
        -fastq2 ${reads[1]} \
	-humRef ${fasta_hg} \
	-virRef ${fasta_virus} \
	-output ${meta.id}__${meta3.virus} \
        ${args}
    """

    stub:
    """
    mkdir -p ${meta.id}__${meta3.virus}/call_fusion_virus
    touch ${meta.id}__${meta3.virus}/call_fusion_virus/HPVfusionPointContig.txt
    """
}

