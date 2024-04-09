process CONSENSUS_REHEADER {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    // e.g: >KX351456 Human|2|PB2|H3N2|USA|A/Rochester/0091/2013|na|na|na|na
    tuple val(meta), path(tsv)

    output:
    //barcode52_minimap2_clair3_segment_7-ref_accession=AB509035
    tuple val(meta), path('*.filtered.txt'), emit: txt

    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    consensus_reheader.py \\
        $args \\
        --input ${tsv} \\
        --prefix ${prefix} \\
        --output ${prefix}.filtered.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
