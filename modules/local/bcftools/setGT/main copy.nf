process BCFTOOLS_SETGT {
    tag "$meta.id"
    label 'process_medium'

    conda 'modules/nf-core/bcftools/consensus/environment.yml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path('*.vcf.gz'), emit: vcf
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    //bcftools +setGT S10/variants/freebayes/S10.norm.vcf -- -t q -i 'strlen(ref) !=strlen(alt) && FORMAT/AD[:1]/FORMAT/DP >= 0.5' -n 'c:1/1' | bcftools +setGT  -- -t q -i 'strlen(ref) !=strlen(alt) && FORMAT/AD[:1]/FORMAT/DP < 0.5' -n 'c:0/0'  | bcftools +setGT -- -t q -i 'GT="1" && FORMAT/AD[:1]/FORMAT/DP < 0.75' -n 'c:0/1' | bcftools +setGT -o  mytest.setGT.vcf.gz -- -t q -i 'GT="1" && FORMAT/AD[:1]/FORMAT/DP >=0.75'  -n 'c:1/1' 

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''

    def prefix = task.ext.prefix ?: "${meta.id}"

    /*
    args: if vaf >= params.upper_ambiguity_freq, set GT="1/1", carry alt
    args2: if vaf < params.upper_ambiguity_freq, set GT="0/1", carray both ref+alt
    args3: for indel, if the vaf >= 0.5, set GT="1/1", 
    args 4, for indel, if the vaf < 0.5, set GT="0/0", carry reference
    */
    """
    #bcftools +setGT ${vcf} ${args} | bcftools +setGT ${args2} | bcftools +setGT ${args3} | bcftools +setGT -o ${prefix}.setgt.vcf.gz ${args4}
    bcftools +setGT ${vcf} ${args} | bcftools +setGT ${args2} | bcftools +setGT ${args3}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
