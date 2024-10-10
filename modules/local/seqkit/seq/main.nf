process SEQKIT_SEQ {
    tag "$meta.id"
    label 'process_me'

    conda (params.enable_conda ? "bioconda::seqkit=2.3.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit%3A2.3.1--h9ee0642_0':
        'quay.io/biocontainers/seqkit%3A2.3.1--h9ee0642_0' }"


    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    seqkit seq --name ${fasta} -o ${prefix}.txt 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit | sed '3!d; s/Version: //' )
    END_VERSIONS
    """
}
