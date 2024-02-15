process MAPPING_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(screen_file) //mash screen ouput format
    tuple val(meta), path(cov_file) //samtool coverage produced txt file
        
    output:
    tuple val(meta), path("*.mapping_summary.csv"), emit: mapping_summary
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    mappingSummary.py \\
        -s ${meta.id} \\
        -m ${screen_file} \\
        -c ${cov_file} \\
        -p ${prefix} \\
        -t ${meta.id}.mapping_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

process REHEADER {
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
    tuple val(meta), path('*.filtered.txt'), emit: txt

    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    reheader.py \\
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

process SEGMENT2TYPEDATA {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    
    tuple val(meta), path(fasta)
    tuple val(meta), path(tsv)

    output:
    
    tuple val(meta), path("*.seg2typedata.tsv"),  optional:true, emit: out_tsv
    tuple val(meta), path("*.segcontig.fa"), optional:true, emit: fasta 
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def nextclade_db_dir = params.nextclade_dataset_base
    def gzipped = fasta.toString().endsWith('.gz')
    def outfile = gzipped ? file(fasta.baseName).baseName : fasta.baseName
    def command = gzipped ? 'zcat' : 'cat'

    """
    $command $fasta > ${outfile}.fixed.fa

    prepare_input_nextclad.py \\
        -t ${tsv} \\
        -f ${outfile}.fixed.fa \\
        -d ${nextclade_db_dir} \\
        -i ${meta.id} \\
        > ${prefix}.seg2typedata.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

process FILTERMASH {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(screen_file)
    
    
    output:
    tuple val(meta), path("*.filtered.screen"), emit: screen
    //tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    filter_mash.py -i ${screen_file} -p ${prefix}.filtered 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}


