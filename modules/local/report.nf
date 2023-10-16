process REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(raw_stats), path(qc_stats), path(consensus_stats), path(blastn_outfmt6), path(mapping_summary), path(nextclade_csv) 
    

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.segments_oneline.csv"), emit: oneline_csv
    tuple val(meta), path("*.segments.csv"), emit: multiline_csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //print(nextclade_csv.getClass())
    //print(nextclade_csv.toString())
    tsvList = nextclade_csv.join(',')
    //print(tsvList)
    
    """
    get_report.py \\
        -s ${meta.id} \\
        -r ${raw_stats} \\
        -q ${qc_stats} \\
        -f ${consensus_stats} \\
        -b ${blastn_outfmt6} \\
        -m ${mapping_summary} \\
        -n ${tsvList} \\
        --prefix ${prefix}

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
