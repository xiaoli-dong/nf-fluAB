process CONSENSUS_REPORT{
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
    tuple val(meta), val(nextclade_dbnames)
    tuple val(meta), path(ref_screen)
    tuple val(meta), val(dbver)
    tuple val(meta), val(pipelinever)

    output:
    tuple val(meta), path("${prefix}.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def stats_file = stats ? "--c-stats-file ${stats}" : ""
    def cov_file = cov ? " --c-coverage-file ${cov}" : ""
    def typing_file = typing ? "--c-typing-file ${typing}" : ""
    def nextclade_csv_files = nextclade_csv ? "--c-nextclade-files ${nextclade_csv}" : ""
    def nextclade_db_files = nextclade_dbnames ? "--c-nextclade-dbnames ${nextclade_dbnames}" : ""
    def mashscreen_file = ref_screen ? "--c-mashscreen-file ${ref_screen}" : ""

    """
    consensus_summary.py \\
        ${mashscreen_file} \\
        ${stats_file} \\
        ${cov_file} \\
        ${typing_file} \\
        ${nextclade_csv_files} \\
        ${nextclade_db_files} \\
        --db-ver ${dbver} \\
        --pipeline-ver ${pipelinever} \\
        > ${prefix}.csv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
