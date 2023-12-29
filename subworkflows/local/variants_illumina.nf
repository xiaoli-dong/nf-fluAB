include {
    BCFTOOLS_MPILEUP
} from '../../modules/nf-core/bcftools/mpileup/main.nf'

workflow variants_bcftools {   

    take:
        bam_bai //[meta, bam]
        fasta_fai //[meta, fasta]
    main:
        ch_versions = Channel.empty()
        bam_bai.join(fasta_fai).multiMap{
            it ->
                bam_interval: [it[0], it[1], []]
                fasta: [it[0], it[3]]
        }.set{
            ch_input
        }
        
        BCFTOOLS_MPILEUP(ch_input.bam_interval, ch_input.fasta, false)
        vcf_tbi = BCFTOOLS_MPILEUP.out.vcf.join(BCFTOOLS_MPILEUP.out.tbi)
        ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)
       
    emit:
        vcf_tbi //[meta, vcf.gz, vcf.gz.tbi]
        fasta_fai //[meta, fasta, fai]
        versions = ch_versions
        
}

include {
    FREEBAYES
} from '../../modules/nf-core/freebayes/main.nf'

include {
    TABIX_TABIX
} from '../../modules/nf-core/tabix/tabix'

include {
    BCFTOOLS_SORT
} from '../../modules/nf-core/bcftools/sort'

workflow variants_freebayes {   

    take:
        bam_bai //[meta, bam, bai]
        fasta_fai // [meta, fasta, fai]
    main:
        ch_versions = Channel.empty()
        bam_bai.join(fasta_fai)
            .multiMap{
                it ->
                    bam_bai: [it[0], it[1], it[2], [], [], []]
                    fasta: it[3]
                    fai: [it[4]]
            }
            .set{
                ch_input
            }

        FREEBAYES (ch_input.bam_bai, ch_input.fasta, ch_input.fai, [], [], [])
        ch_versions.mix(FREEBAYES.out.versions)
       
       
        //without sorting it gives error: [E::hts_idx_push] Unsorted positions on seque
        BCFTOOLS_SORT(FREEBAYES.out.vcf) //emit vcf.gz
        ch_versions.mix(BCFTOOLS_SORT.out.versions)

        TABIX_TABIX(BCFTOOLS_SORT.out.vcf)  //emit [meta, tbi]
        ch_versions.mix(TABIX_TABIX.out.versions)

        vcf_tbi = BCFTOOLS_SORT.out.vcf.join(TABIX_TABIX.out.tbi)

      emit:
        vcf_tbi  //[meta, vcf.gz, tbi]
        fasta_fai
        versions = ch_versions
        
}
