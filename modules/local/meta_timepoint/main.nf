process META_TIMEPOINT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reads), path(index)
    path(library)

    output:
    tuple val(meta), path(reads), path("timepoint.txt"), emit: reads

    script:
    def pool = ""
    if (params.key_pool) { 
        pool = meta.id.tokenize(params.pool_delim).find { it.toLowerCase().startsWith('pool') } 
    }
 
    """
    label=\$(awk '{print \$1"_"\$2}' $index)
    if [[ "${pool}" != "" ]]; then
        label=\$label"_"$pool
    fi
    awk -F'\t' -v label=\${label} '{if(\$1==label) print \$2}' ${library} > timepoint.txt
    """

    stub:
    """
    touch timepoint.txt 
    """
}

