process CLAIR3 {
  tag "$meta.id"
  label 'process_medium'

  conda 'bioconda::clair3==1.0.8'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/clair3%3A1.0.8--py39hf5e1c6e_0'
  } else {
    container 'quay.io/biocontainers/clair3%3A1.0.8--py39hf5e1c6e_0'
  }

  input:
  tuple val(meta), path(bam), path(bai) //indexed bam
  tuple val(meta), path(fasta), path(fai) //indexed ref
  path model_path 

  output:
  tuple val(meta), path(vcf), path(tbi),  optional:true, emit: vcf_tbi
  tuple val(meta), path(fasta), emit: fasta
  tuple val(meta), path(fasta), path(fai), emit: fasta_fai
  path (clair3_dir), emit: output_dir
  path (clair3_log), emit: log
  path "versions.yml" , emit: versions

  when:
    task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"

  vcf         = "${prefix}.clair3.gvcf.gz"
  tbi         = "${prefix}.clair3.gvcf.gz.tbi"  
  clair3_dir   = "${prefix}.clair3"
  clair3_log   = "${clair3_dir}/run_clair3.log"
  
  """
  run_clair3.sh \\
      --gvcf \\
      --bam_fn=${bam} \\
      --ref_fn=$fasta \\
      --model_path=${model_path} \\
      --threads=${task.cpus} \\
      --output=${clair3_dir} \\
      $args 
  
  if [ -f ${clair3_dir}/merge_output.gvcf.gz ]; then
        ln -s ${clair3_dir}/merge_output.gvcf.gz ${vcf}
        ln -s ${clair3_dir}/merge_output.gvcf.gz.tbi ${tbi}
  fi
 
  
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      clair3: \$(head -n1 ${clair3_dir}/run_clair3.log | sed 's/^.*CLAIR3 VERSION: v//; s/ .*\$//')
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}