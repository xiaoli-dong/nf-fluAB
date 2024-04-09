include {
    CLAIR3
} from '../../modules/local/clair3/main'

include {
    TABIX_TABIX as TABIX_TABIX_VCF;
} from '../../modules/nf-core/tabix/tabix'

workflow VARIANTS_NANOPORE {   

    take:
        bam_bai //[meta, bam, bai]
        fasta_fai // [meta, fasta, fai]
        //model_path
    main:
        ch_versions = Channel.empty()   
        vcf = Channel.empty()

        if(params.nanopore_variant_caller == 'clair3'){
            bam_bai.join(fasta_fai).multiMap {
                it ->
                    bam_bai: [it[0], it[1], it[2]]
                    fasta_fai: [it[0], it[3], it[4]]
            }
            .set{
                ch_input
            }

            CLAIR3(
                ch_input.bam_bai, 
                ch_input.fasta_fai, 
                Channel.value(file(params.clair3_variant_model))
            )
            ch_versions = ch_versions.mix(CLAIR3.out.versions)

            TABIX_TABIX_VCF(CLAIR3.out.vcf)
            vcf_tbi = CLAIR3.out.vcf.join(TABIX_TABIX_VCF.out.tbi)
           
           
        }
       
    emit:
        vcf_tbi  //[meta, vcf.gz]
        versions = ch_versions
        
}
