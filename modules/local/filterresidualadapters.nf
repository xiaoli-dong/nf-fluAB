process FILTERRESIDUALADAPTERS {
    tag "$meta.id"
    label 'process_medium'

    // //conda (params.enable_conda ? "conda-forge::python=3.8.3 bioconda::pysam=0.16.01" : null)
    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     "https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:1744f68fe955578c63054b55309e05b41c37a80d-0" :
    //     "quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:1744f68fe955578c63054b55309e05b41c37a80d-0" }"

    conda (params.enable_conda ? "bioconda::pysam=0.19.0 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-57736af1eb98c01010848572c9fec9fff6ffaafd:402e865b8f6af2f3e58c6fc8d57127ff0144b2c7-0' :
        'quay.io/biocontainers/mulled-v2-57736af1eb98c01010848572c9fec9fff6ffaafd:402e865b8f6af2f3e58c6fc8d57127ff0144b2c7-0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_posttrim_filter.fq.gz"), emit: reads
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    [ ! -f  ${prefix}_1.fastp.fq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastp.fq.gz
    [ ! -f  ${prefix}_2.fastp.fq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastp.fq.gz

    filter_residual_adapters.py \\
        --input_R1 ${prefix}_1.fastp.fq.gz \\
        --input_R2 ${prefix}_2.fastp.fq.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_residual_adapters.py: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
