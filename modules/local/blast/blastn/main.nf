process BLAST_BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::blast=2.14.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast%3A2.14.1--pl5321h6f7f691_0':
        'quay.io/biocontainers/blast%3A2.14.1--pl5321h6f7f691_0' }"

    input:
    tuple val(meta), path(fasta)
    path  db //fasta reference file

    output:
    tuple val(meta), path('*.typing_output.tsv'), emit: tsv
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    

    """
    
    if [[ $fasta == *.gz ]] 
    then
        gzip -dc $fasta | blastn -num_threads $task.cpus -subject ${db}  $args -out ${prefix}.blastn.tsv
        awk '\$NF > ${args2} {if (!seen[\$1]++) print \$0}' ${prefix}.blastn.tsv > ${prefix}.typing_output.tsv
    else
        blastn -num_threads $task.cpus -subject ${db} -query $fasta $args -out ${prefix}.blastn.tsv
        awk '\$NF > ${args2} {if (!seen[\$1]++) print \$0}' ${prefix}.blastn.tsv > ${prefix}.typing_output.tsv
    fi
    
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """

}
