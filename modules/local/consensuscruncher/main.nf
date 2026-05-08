process CONSENSUSCRUNCHER {
    tag "${meta.id}_${meta.ref}${meta.type}_${meta.ref2}${meta.type2}_${meta.consensus}"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(cytoband)

    output:
    tuple val(meta), path("${bam.baseName}"), emit: consensus_dir
    path("${bam.baseName}/*.read_families.txt"), optional: true, emit: optional_read_families
    tuple val(meta), path("${bam.baseName}/dcs_sc/*.dcs.sc.sorted.bam"), emit: bam_dcssc
    tuple val(meta), path("${bam.baseName}/dcs_sc/*.dcs.sc.sorted.bam.bai"), emit: bai_dcssc
    tuple val(meta), path("${bam.baseName}/dcs_sc/*.all.unique.dcs.sorted.bam"), emit: bam_alluniquedcs
    tuple val(meta), path("${bam.baseName}/dcs_sc/*.all.unique.dcs.sorted.bam.bai"), emit: bai_alluniquedcs
    tuple val(meta), path("${bam.baseName}/dcs/*.dcs.sorted.bam"), emit: bam_dcs
    tuple val(meta), path("${bam.baseName}/dcs/*.dcs.sorted.bam.bai"), emit: bai_dcs
    tuple val(meta), path("${bam.baseName}/sscs/*.sscs.sorted.bam"), emit: bam_sscs
    tuple val(meta), path("${bam.baseName}/sscs/*.sscs.sorted.bam.bai"), emit: bai_sscs
    tuple val(meta), path("${bam.baseName}/sscs_sc/*.sscs.sc.sorted.bam"), emit: bam_sscssc
    tuple val(meta), path("${bam.baseName}/sscs_sc/*.sscs.sc.sorted.bam.bai"), emit: bai_sscssc

    script:
    def samp = bam.baseName
    def config = "config_${meta.id}_${meta.ref}${meta.type}_${meta.ref2}${meta.type2}.ini"
    """
    echo "[fastq2bam]" > $config
    echo "fastq1 = skip" >> $config
    echo "fastq2 = skip" >> $config
    echo "output = skip" >> $config
    echo "name = ${samp}" >> $config
    echo "bwa = skip" >> $config
    echo "ref = skip" >> $config
    echo "samtools = ${params.samtools_path}" >> $config
    echo "bpattern = skip" >> $config
    echo "[consensus]" >> $config
    echo "bam = ${bam}" >> $config
    echo "c_output = ." >> $config
    echo "bedfile = ${cytoband}" >> $config
    python3 ${params.consensuscruncher_dir}/ConsensusCruncher.py -c $config consensus
    """

    stub:
    def samp = bam.baseName
    def config = "config_${meta.id}_${meta.ref}${meta.type}_${meta.ref2}${meta.type2}.ini"
    """
    mkdir -p $samp/dcs_sc
    mkdir -p $samp/dcs
    mkdir -p $samp/sscs
    mkdir -p $samp/sscs_sc
    touch $config
    touch $samp/${samp}.read_families.txt 
    touch $samp/dcs_sc/${samp}.dcs.sc.sorted.bam 
    touch $samp/dcs_sc/${samp}.dcs.sc.sorted.bam.bai 
    touch $samp/dcs_sc/${samp}.all.unique.dcs.sorted.bam 
    touch $samp/dcs_sc/${samp}.all.unique.dcs.sorted.bam.bai 
    touch $samp/dcs/${samp}.dcs.sorted.bam 
    touch $samp/dcs/${samp}.dcs.sorted.bam.bai 
    touch $samp/sscs/${samp}.sscs.sorted.bam 
    touch $samp/sscs/${samp}.sscs.sorted.bam.bai 
    touch $samp/sscs_sc/${samp}.sscs.sc.sorted.bam 
    touch $samp/sscs_sc/${samp}.sscs.sc.sorted.bam.bai 
    """
}

