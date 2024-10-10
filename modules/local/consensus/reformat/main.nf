process CONSENSUS_REFORMAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9':
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.consensus.fa'), emit: fasta

    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    consensus_reformat.py \\
        $args \\
        --fasta ${fasta} \\
        --prefix ${prefix} \\
        > ${prefix}.consensus.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
