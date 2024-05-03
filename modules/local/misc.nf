
process SEGMENT2DATASET {
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
    
    tuple val(meta), path("*.seg2dataset.tsv"),  optional:true, emit: out_tsv
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

    seg2dataset.py \\
        -t ${tsv} \\
        -f ${outfile}.fixed.fa \\
        -d ${nextclade_db_dir} \\
        -i ${meta.id} \\
        > ${prefix}.seg2dataset.tsv

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
    val min_shared_hashes
    
    
    output:
    tuple val(meta), path("*keep_best.screen"), emit: screen
    //tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    filter_mash.py -i ${screen_file} -p ${prefix}-keep_best -m ${min_shared_hashes}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}


process SAMCLIP {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::perl=5.26.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2':
        'biocontainers/perl:5.26.2' }"

    input:
    tuple val(meta), path(sam)
    tuple val(meta), path(fasta), path(fai)
    
    
    output:
    tuple val(meta), path("*.samclip.sam"), emit: sam
    //tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
   
    samclip --ref ${fasta} ${args} ${sam} > ${prefix}.samclip.sam 2>${prefix}.error.log 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl -v 2>&1) | sed -n \'2 p\' | sed 's/^.*?(v//g; s/^).*//g;' ))
    END_VERSIONS
    """
}

process RENAMECHROM {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(bed)
    tuple val(meta), path(header)
    
    
    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    rename_bed_chrom.py -b ${bed} -i ${header} -o ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
process PLOT_DEPTH {
     tag "$meta.id"
    label 'process_single'


    conda "conda-forge::r-base=4.0.3 conda-forge::r-optparse=1.6.6 conda-forge::r-ggplot2=3.3.3 "
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0fb10cdfa7:05678da05b8e5a7a5130e90a9f9a6c585b965afa-0' :
        'quay.io/biocontainers/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0fb10cdfa7:05678da05b8e5a7a5130e90a9f9a6c585b965afa-0' }"

    input:
    tuple val(meta), path(genomecov)

    output:
    tuple val(meta), path('*.pdf'), emit: pdf
    //tuple val(meta), path('*.jpeg'), emit: png
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plot_depth.R $genomecov $prefix 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
