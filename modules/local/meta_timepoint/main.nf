process META_TIMEPOINT {
    tag "${meta.id}_METATIMEPOINT"
    label 'process_single'

    input:
    tuple val(meta), path(reads), path(index)
    path(library)

    output:
    tuple val(meta), path(reads), path("timepoint.txt"), emit: reads

    script:
    def pool = ""
    if (params.key_pool) { 
        pool = meta.id.tokenize(params.pool_delim).find { it.toLowerCase().startsWith('pool') }
        if (pool.toLowerCase() == "pool") {
            def ext = meta.id.split("${pool}${params.pool_delim}")[1].split(params.pool_delim)[0]
            pool = "${pool}_${ext}"
        }
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

