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
