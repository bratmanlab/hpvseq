process SUMMARIZE_INTEGRATION {
    tag "summarize_integration"

    input:
    path integration_files 

    output:
    path("*integration.summ.txt"), emit: summ_integration

    script:
    integration_file1       = integration_files[0]
    summ_integration_file   = "integration.summ.txt"
    """
    ## integration
    head -1 ${integration_file1} | awk -F'\t' '{OFS="\t"; print "Sample", \$0}' > ${summ_integration_file} 
    for file in ${integration_files};do
        samp=\${file/_HPVfusionPointContig*}
        awk -F'\t' -v samp=\$samp 'NR>1{OFS="\t"; print samp, \$0}' \$file >> ${summ_integration_file}
        done
    """

    stub:
    """
    touch integration.summ.txt
    """
}

