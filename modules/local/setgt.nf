process SETGT {
    tag "$meta.id"
    label 'process_single'

   
    conda "bioconda::vcfpy=0.13.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcfpy%3A0.13.8--pyhdfd78af_0':
        'biocontainers/vcfpy:0.13.8--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    
    output:
    tuple val(meta), path("*.setGT.vcf"), emit: vcf
    
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
   
    """
    vcf_setGT.py \\
        ${args} \\
        --input ${vcf} \\
        --tbi ${tbi} \\
        --output-vcf ${prefix}.setGT.vcf \\
        

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfpy: 0.13.8
    END_VERSIONS
    """
}
