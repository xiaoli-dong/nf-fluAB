process SPLITBAM2SEGMENTS {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:

    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    tuple val(meta), val(tlist)
    //tuple val(meta), path(screen)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*_working_segment_*.bam"), emit: bam
    tuple val(meta), path("*_working_segment_*.bai"), emit: bai
    tuple val(meta), path("*.cov"), emit: cov
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def cmd = ''
    for( int i=0; i< tlist.size(); i++ ) {
        cmd += "samtools view -h ${bam} ${tlist[i].acc} -b -o ${prefix}_working_segment_${tlist[i].segid}.bam;\n"
        cmd += "samtools index ${prefix}_working_segment_${tlist[i].segid}.bam;\n"
        cmd += " samtools depth -d 0 -a ${prefix}_working_segment_${tlist[i].segid}.bam > ${prefix}_${tlist[i].segid}.cov; \n"
    }
    cmd

    """
    $cmd
    cat ${prefix}_*.cov > ${prefix}.cov

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        splitbam2segments: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
