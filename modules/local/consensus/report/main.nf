process CONSENSUS_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(stats)
    tuple val(meta), path(cov)
    tuple val(meta), path(typing)
    tuple val(meta), val(nextclade_csv) 

    output:
    tuple val(meta), path("*.consensus_summary.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    consensus_summary.py \\
        --c-stats-file ${stats} \\
        --c-coverage-file ${cov} \\
        --c-typing-file ${typing} \\
        --c-nextclade-files ${nextclade_csv} \\
        > ${prefix}.consensus_summary.csv

   
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
