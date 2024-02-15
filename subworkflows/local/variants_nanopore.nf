include {
    CLAIR3
} from '../../modules/local/clair3/main'

workflow variants_clair3 {   

    take:
        bam_bai //[meta, bam, bai]
        fasta_fai // [meta, fasta, fai]
        model_path
    main:
        ch_versions = Channel.empty()

        CLAIR3(
            bam_bai, 
            fasta_fai, 
            model_path
        )
        ch_versions = ch_versions.mix(CLAIR3.out.versions)

        CLAIR3.out.vcf.join(fasta_fai).multiMap{
            it ->
                vcf: [it[0], it[1], it[2]]
                fasta: [it[0], it[3]]
        }.set{
            ch_input
        }

        
       
    emit:
        vcf = CLAIR3.out.vcf //[meta, vcf.gz]
        fasta_fai
        versions = ch_versions
        
}
