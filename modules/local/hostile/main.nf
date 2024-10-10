process HOSTILE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::hostile=0.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hostile%3A0.4.0--pyhdfd78af_0':
        'biocontainers/hostile%3A0.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    val(aligner)
    path(ref) //directory, index has same prefix as the directory name

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "--fastq1 ${reads[0]}" : "--fastq1 ${reads[0]} --fastq2 ${reads[1]}"
    def index = aligner =~ /bowtie2/ ? "${ref}/${ref}" : "${ref}"

   // if (aligner =~ /bowtie2/) {
        """
        hostile clean \\
            --aligner ${aligner} \\
            --out-dir ./ \\
            --index ${index} \\
            --threads $task.cpus \\
            ${input} \\
            $args \\
            >& ${prefix}.hostile_${aligner}.log
        
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hostile: \$(hostile --version)
        END_VERSIONS
        """
    /* } else if (aligner =~ /minimap2/) {
        
        """
        hostile clean \\
            --aligner ${aligner} \\
            --out-dir ./ \\
            --index ${index} \\
            --threads $task.cpus \\
            ${input} \\
            $args \\
            >& ${prefix}.hostile_${aligner}.log
        
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hostile: \$(hostile --version)
        END_VERSIONS
        """
    } */
}
