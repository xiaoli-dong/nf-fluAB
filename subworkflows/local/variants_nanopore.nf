include {
    CLAIR3
} from '../../modules/local/clair3/main'

workflow VARIANTS_NANOPORE {

    take:
        bam_bai //[meta, bam, bai]
        fasta_fai // [meta, fasta, fai]
        //model_path
    main:
        ch_versions = Channel.empty()
        vcf = Channel.empty()
        /*
            Clair3 was seeing no read at 2265. Clair3 filters the alignments with the following four flags:
            read unmapped (0x4)
            mate unmapped (0x8)
            not primary alignment (0x100)
            supplementary alignment (0x800)
            Using samtools, the alignments can be filtered with samtools view -F 2316

            Clair3 doesn't call variants in the first 16bp and last 16bp of a sequence because of
            1) algorithmic limit, and
            2) usually degenerated coverage and alignment performance in the head and tail of a sequence
            that makes variant calling unreliable.
            A solution to do the variants calling at the end is to add some 'N' to the tail of your reference genome
            so the end gets out of the tail 16bp limit.
        */
        vcf_tbi = Channel.empty()
        if(params.variant_caller == 'clair3'){
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
            vcf_tbi = CLAIR3.out.vcf_tbi
            ch_versions = ch_versions.mix(CLAIR3.out.versions)



        }

    emit:
        vcf_tbi //[meta, vcf.gz]
        versions = ch_versions

}
