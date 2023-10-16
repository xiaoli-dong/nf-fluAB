process SEQKIT_GREP {
    tag "$meta.id"
    label 'process_me'

    conda (params.enable_conda ? "bioconda::seqkit=2.3.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit%3A2.3.1--h9ee0642_0':
        'quay.io/biocontainers/seqkit%3A2.3.1--h9ee0642_0' }"


    input:

    tuple val(meta), path(screen)
    path refdb

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*_working_reference.fa"), emit: fasta
    tuple val(meta), path("*_segment_*.fa"), emit: segments
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //cut -f3 ${screen} | seqkit grep -f - ${refdb} > ${prefix}_working_reference.fa
    """

    while read segment accession
    do
	    seqkit grep --pattern \$accession ${refdb} > ${prefix}_working_segment_\$segment.fa

    done < <(cut -f1,3 ${screen})

    cat ${prefix}_working_segment_*.fa > ${prefix}_working_reference.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit | sed '3!d; s/Version: //' )
    END_VERSIONS
    """
}
