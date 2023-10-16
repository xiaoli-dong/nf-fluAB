process MERGE_SEGMENTS {
    tag "$meta.id"
    label 'process_me'

    conda (params.enable_conda ? "bioconda::seqkit=2.3.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit%3A2.3.1--h9ee0642_0':
        'quay.io/biocontainers/seqkit%3A2.3.1--h9ee0642_0' }"


    input:

    tuple val(meta), path(fa)

    output:
    tuple val(meta), path("*_working_reference.fa"), emit: working_reference
    //path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat ${fa} > ${prefix}_working_reference.fa

    """
}
