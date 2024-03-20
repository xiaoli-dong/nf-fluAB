include {
    SAMTOOLS_VIEW;
    SAMTOOLS_VIEW as SAMTOOLS_VIEW_2SAM;
} from '../../modules/nf-core/samtools/view/main'

include {
    SAMCLIP;  
} from '../../modules/local/misc'

include {
    SAMTOOLS_INDEX;
    SAMTOOLS_INDEX as  SAMTOOLS_INDEX_CLEAN;
} from '../../modules/nf-core/samtools/index/main'

include {
    SAMTOOLS_SORT as SAMTOOLS_SORT_CLEAN;
} from '../../modules/nf-core/samtools/sort/main'

include {
    SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_CLEAN_RMDUP
} from '../../modules/nf-core/samtools/coverage/main'

include {
    BAM_MARKDUPLICATES_PICARD;
} from '../nf-core/bam_markduplicates_picard'

workflow PREPROCESS_BAM {   
    take:
        bam_bai
        fasta_fai //[meta, fasta, fai] reference used to produce bam_bai

    main:
        ch_versions = Channel.empty()
        //referene: https://wiki.bits.vib.be/index.php/Call_variants_with_samtools_1.0

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            filter out secondary and supplementary 
            samtools view ext.args = " -h -F0x900 --output-fmt bam -b"
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        
        
        SAMTOOLS_VIEW(
            bam_bai, 
            fasta_fai.map{
                meta, fasta, fai ->[meta, fasta]
                }, 
            []
        )
        SAMTOOLS_INDEX(
             SAMTOOLS_VIEW.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())
        
        ch_bam_bai = SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_INDEX.out.bai).view()
        
        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            filter sam for soft and hard clipped alignment 
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        

        if(!params.skip_rm_soft_hard_clip){

            ch_bam_bai.join(fasta_fai).multiMap{
            it->
                bam_bai: [it[0], it[1], it[2]]
                fasta: [it[0], it[3]]
            }.set{
                ch_input
            }
           

            //convert bam to SAM
            SAMTOOLS_VIEW_2SAM(
                ch_input.bam_bai, 
                ch_input.fasta,
                []
            )

            SAMTOOLS_VIEW_2SAM.out.sam.join(fasta_fai).multiMap{
                it ->
                    sam: [it[0], it[1]]
                    fasta_fai: [it[0], it[2], it[3]]
            }.set{ 
                ch_input 
            }

            SAMCLIP(ch_input.sam, ch_input.fasta_fai)
            ch_versions = ch_versions.mix(SAMCLIP.out.versions.first())

            SAMTOOLS_SORT_CLEAN (
                SAMCLIP.out.sam
            )
            ch_versions = ch_versions.mix(SAMTOOLS_SORT_CLEAN.out.versions.first())

            SAMTOOLS_INDEX_CLEAN(
                SAMTOOLS_SORT_CLEAN.out.bam
            )

            ch_bam_bai = SAMTOOLS_SORT_CLEAN.out.bam.join(SAMTOOLS_INDEX_CLEAN.out.bai)
        }

         /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            remove duplicates
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        ch_bam_bai.join(fasta_fai).multiMap{
            it ->
                bam: [it[0], it[1]]
                fasta: [it[0], it[3]]
                fai: [it[0], it[4]]
        }.set{ ch_input }
    
        //all the duplicates will be removed from the output
        BAM_MARKDUPLICATES_PICARD(ch_input.bam, ch_input.fasta, ch_input.fai)

        ch_bam_bai = BAM_MARKDUPLICATES_PICARD.out.bam.join(BAM_MARKDUPLICATES_PICARD.out.bai)
        
        SAMTOOLS_COVERAGE_CLEAN_RMDUP(ch_bam_bai)
        coverage = SAMTOOLS_COVERAGE_CLEAN_RMDUP.out.coverage
        
        
    emit:
        bam_bai = ch_bam_bai // bam file with secondary, supplementary, soft/hard clipped reads, duplicate reads removed
        coverage = SAMTOOLS_COVERAGE_CLEAN_RMDUP.out.coverage
        versions = ch_versions
        
}