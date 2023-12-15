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
        -t ${prefix}.mapping_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}



process FORMAT_CONSENSUS {
    tag "$sample"
    conda 'bioconda::shiptv=0.4.0'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/shiptv:0.4.0--pyh5e36f6f_0'
    } else {
    container 'quay.io/biocontainers/shiptv:0.4.0--pyh5e36f6f_0'
    }

    input:
    // e.g: >KX351456 Human|2|PB2|H3N2|USA|A/Rochester/0091/2013|na|na|na|na
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.consensus_with_name.fasta'), emit: consensus_with_name
    tuple val(meta), path('*.consensus.fasta'), emit: consensus_fasta

    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    formatFasta.py \\
    --sample-name ${prefix} \\
    --output1-fasta ${prefix}.consensus_with_name.fasta \\
    --output2-fasta ${prefix}.consensus.fasta \\
    $fasta

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

process PROCESSGVCF {
    tag "$meta.id"
    label 'process_medium'

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     "https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:1744f68fe955578c63054b55309e05b41c37a80d-0" :
    //     "quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:1744f68fe955578c63054b55309e05b41c37a80d-0" }"
    conda (params.enable_conda ? "bioconda::pysam=0.19.0 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-57736af1eb98c01010848572c9fec9fff6ffaafd:402e865b8f6af2f3e58c6fc8d57127ff0144b2c7-0' :
        'quay.io/biocontainers/mulled-v2-57736af1eb98c01010848572c9fec9fff6ffaafd:402e865b8f6af2f3e58c6fc8d57127ff0144b2c7-0' }"

    input:

    tuple val(meta), path(vcf)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.mask.txt"), emit: mask
    tuple val(meta), path("*.variants.vcf"), emit: variants
    tuple val(meta), path("*.consensus.vcf"), emit: consensus
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    process_gvcf.py \\
        ${args} \\
        -m ${prefix}.mask.txt \\
        -v ${prefix}.variants.vcf \\
        -c ${prefix}.consensus.vcf \\
        ${vcf}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        processgvcf: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}

process SPLITVCF {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::tabix=1.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix%3A1.11--hdfd78af_0':
        'quay.io/biocontainers/tabix%3A1.11--hdfd78af_0' }"

    input:

    tuple val(meta), path(vcf)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.vcf"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    zcat $vcf | awk $args '\$0 ~ /^#/ || \$0 ~ vartag' > ${prefix}.vcf
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


