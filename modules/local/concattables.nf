process CONCATTABLES {
    tag "$meta.id"
    label 'process_low'

    //conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas%3A2.2.1' :
        'quay.io/biocontainers/pandas%3A2.2.1' }"

    input:
    tuple val(meta), path(input_tsv)

    output:
    tuple val(meta), path("${prefix}.${output_ext}"), emit: output_tsv
    path "versions.yml"

    script: // This script is bundled with the pipeline, in nf-core/influenza/bin/
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    output_ext = task.ext.output_ext ?: 'tsv'

    """
    concat_tables.py $args -o ${prefix}.${output_ext} ${input_tsv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
