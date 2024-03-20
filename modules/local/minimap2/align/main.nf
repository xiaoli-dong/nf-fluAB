process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "bioconda::minimap2=2.26"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2%3A2.26--he4a0461_2' :
        'quay.io/biocontainers/minimap2%3A2.26--he4a0461_2' }"
    input:
    tuple val(meta), path(reads)
    path reference
    val sam_format
    
    output:
    tuple val(meta), path("*.paf"), optional: true, emit: paf
    tuple val(meta), path("*.sam"), optional: true, emit: sam
    path "versions.yml"           , emit: versions

    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //def sam_output = sam_format ? "-a | samtools sort | samtools view -@ ${task.cpus} -b -h -o ${prefix}.bam" : "-o ${prefix}.paf"
    def sam_output = sam_format ? "-a -o ${prefix}.sam" : "-o ${prefix}.paf"

    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        ${reference ?: reads} \\
        $reads \\
        $sam_output


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
