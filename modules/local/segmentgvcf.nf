process SEGMENTGVCF {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::freebayes=1.3.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.6--hb089aa1_0':
        'quay.io/biocontainers/freebayes:1.3.6--hb089aa1_0' }"

    input:

    tuple val(meta), path(fasta)
    tuple val(meta), path(fai)
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    tuple val(meta), val(tlist)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.gvcf"), emit: gvcf
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def cmd = ''
    for( int i=0; i< tlist.size(); i++ ) {

        cmd += "freebayes $args --fasta-reference ${prefix}_working_segment_${tlist[i].segid}.fa ${prefix}_working_segment_${tlist[i].segid}.bam > ${prefix}_working_segment_${tlist[i].segid}.gvcf;\n"
    }
    cmd

    """
    $cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segmentgvcf: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}
