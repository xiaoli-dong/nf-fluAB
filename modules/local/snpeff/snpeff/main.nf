process SNPEFF_SNPEFF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1' :
        'biocontainers/snpeff:A5.2--hdfd78af_1' }"

    input:
    tuple val(meta), path(vcf)
    val   db
    //tuple val(meta2), path(cache)
    path(snpeff_config)
    path(dataDir)

    output:
    tuple val(meta), path("*.snpeff.vcf"),   emit: vcf
    tuple val(meta), path("*.stats.csv"),       emit: report
    tuple val(meta), path("*.html"),      emit: summary_html
    tuple val(meta), path("*.genes.txt"), emit: genes_txt
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 6144
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    //def cache_command = cache ? "-dataDir \${PWD}/${cache}" : ""
    //-XX:-UsePerfData Disabling perfdata feature. suppresses the creation of the hsperfdata_userid directories. 
    """
    snpEff -XX:-UsePerfData \\
        $db \\
        -csvStats ${prefix}.stats.csv \\
        -config $snpeff_config \\
        -dataDir $dataDir \\
        -i vcf -o vcf \\
        $vcf > ${prefix}.snpeff.vcf
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ann.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """

}
