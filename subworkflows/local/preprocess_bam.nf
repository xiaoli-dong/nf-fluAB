include {
    SAMTOOLS_VIEW as SAMTOOLS_VIEW_RM_SEC_SUP;
    SAMTOOLS_VIEW as SAMTOOLS_VIEW_BAM2SAM;
} from '../../modules/nf-core/samtools/view/main'

include {
    RENAMECHROM as RENAMECHROM_PERBASE;
    RENAMECHROM as RENAMECHROM_PERBASE_RM_SEC_SUP;
    RENAMECHROM as RENAMECHROM_PERBASE_SAMCLIP;
    RENAMECHROM as RENAMECHROM_PERBASE_RMDUP;
    SAMCLIP; 
    PLOT_DEPTH as PLOT_DEPTH_INPUT;
    PLOT_DEPTH as PLOT_DEPTH_RM_SEC_SUP;
    PLOT_DEPTH as PLOT_DEPTH_SAMCLIP;
    PLOT_DEPTH as PLOT_DEPTH_RMDUP;

} from '../../modules/local/misc'

include {
    SAMTOOLS_INDEX as SAMTOOLS_INDEX_RM_SEC_SUP;
    SAMTOOLS_INDEX as SAMTOOLS_INDEX_SAMCLIP;
} from '../../modules/nf-core/samtools/index/main'

include {
    SAMTOOLS_SORT as SAMTOOLS_SORT_SAMCLIP;
} from '../../modules/nf-core/samtools/sort/main'

include {
    SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_RM_SEC_SUP;
    SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_SAMCLIP;
    SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_RMDUP;
} from '../../modules/nf-core/samtools/coverage/main'

include {
    BAM_MARKDUPLICATES_PICARD as PICARD_MARKDUPLICATES;
} from '../nf-core/bam_markduplicates_picard'

include {
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_PERBASE;
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_PERBASE_RM_SEC_SUP;
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_PERBASE_SAMCLIP;
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_PERBASE_RMDUP;
} from '../../modules/nf-core/bedtools/genomecov/main.nf'


workflow PREPROCESS_BAM {   
        take:
            bam_bai
            fasta_fai //[meta, fasta, fai] reference used to produce bam_bai
            ref_header //[meta, txt] reference header

        main:
            ch_versions = Channel.empty()
            //ref_header.view()
            //get perbase depth file and plot it
            ch_input = bam_bai.map{ it -> [it[0], it[1], 1] } //[meta, bam, scale]
            BEDTOOLS_GENOMECOV_PERBASE(ch_input, [], "bed")

            BEDTOOLS_GENOMECOV_PERBASE.out.genomecov
                .join(ref_header)
                .multiMap{ 
                    it -> 
                        bed: [it[0], it[1]] //[meta, bed]
                        header: [it[0], it[2]] //[meta, header] 
                }.set{ ch_input }
            //ch_input.header.view()
            RENAMECHROM_PERBASE(ch_input.bed, ch_input.header)

            PLOT_DEPTH_INPUT(RENAMECHROM_PERBASE.out.tsv)
            //PLOT_DEPTH_INPUT(BEDTOOLS_GENOMECOV_PERBASE.out.genomecov)
            ch_versions = ch_versions.mix(PLOT_DEPTH_INPUT.out.versions.first())


        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            filter out secondary and supplementary 
            samtools view ext.args = " -h -F0x900 --output-fmt bam -b"
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        SAMTOOLS_VIEW_RM_SEC_SUP(
            bam_bai, 
            fasta_fai.map{
                meta, fasta, fai ->[meta, fasta]
                }, 
            []
        )
        SAMTOOLS_INDEX_RM_SEC_SUP(
             SAMTOOLS_VIEW_RM_SEC_SUP.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW_RM_SEC_SUP.out.versions.first())
        ch_bam_bai = SAMTOOLS_VIEW_RM_SEC_SUP.out.bam.join(SAMTOOLS_INDEX_RM_SEC_SUP.out.bai)//.view()
        
        SAMTOOLS_COVERAGE_RM_SEC_SUP(ch_bam_bai)
        coverage_rm_sec_sup = SAMTOOLS_COVERAGE_RM_SEC_SUP.out.coverage

        ch_input = ch_bam_bai.map{ it -> [it[0], it[1], 1] } //[meta, bam, scale]
        BEDTOOLS_GENOMECOV_PERBASE_RM_SEC_SUP(ch_input, [], "bed")

        BEDTOOLS_GENOMECOV_PERBASE_RM_SEC_SUP.out.genomecov
                .join(ref_header)
                .multiMap{ 
                    it -> 
                        bed: [it[0], it[1]] //[meta, bed]
                        header: [it[0], it[2]] //[meta, header] 
                }.set{ ch_input }
            //ch_input.header.view()
        RENAMECHROM_PERBASE_RM_SEC_SUP(ch_input.bed, ch_input.header)

        PLOT_DEPTH_RM_SEC_SUP(RENAMECHROM_PERBASE_RM_SEC_SUP.out.tsv)
        
        
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
            SAMTOOLS_VIEW_BAM2SAM(
                ch_input.bam_bai, 
                ch_input.fasta,
                []
            )

            SAMTOOLS_VIEW_BAM2SAM.out.sam.join(fasta_fai).multiMap{
                it ->
                    sam: [it[0], it[1]]
                    fasta_fai: [it[0], it[2], it[3]]
            }.set{ 
                ch_input 
            }

            SAMCLIP(ch_input.sam, ch_input.fasta_fai)
            ch_versions = ch_versions.mix(SAMCLIP.out.versions.first())

             SAMTOOLS_SORT_SAMCLIP (
                SAMCLIP.out.sam
            )
            ch_versions = ch_versions.mix(SAMTOOLS_SORT_SAMCLIP.out.versions.first())

             SAMTOOLS_INDEX_SAMCLIP(
                 SAMTOOLS_SORT_SAMCLIP.out.bam
            )

            ch_bam_bai =  SAMTOOLS_SORT_SAMCLIP.out.bam.join(SAMTOOLS_INDEX_SAMCLIP.out.bai)
            SAMTOOLS_COVERAGE_SAMCLIP(ch_bam_bai)
            coverage_SAMCLIP = SAMTOOLS_COVERAGE_SAMCLIP.out.coverage

            ch_input = ch_bam_bai.map{ it -> [it[0], it[1], 1] } //[meta, bam, scale]
            BEDTOOLS_GENOMECOV_PERBASE_SAMCLIP(ch_input, [], "bed")
            //PLOT_DEPTH_SAMCLIP(BEDTOOLS_GENOMECOV_PERBASE_SAMCLIP.out.genomecov)

            BEDTOOLS_GENOMECOV_PERBASE_SAMCLIP.out.genomecov
                .join(ref_header)
                .multiMap{ 
                    it -> 
                        bed: [it[0], it[1]] //[meta, bed]
                        header: [it[0], it[2]] //[meta, header] 
                }.set{ ch_input }
            //ch_input.header.view()
            RENAMECHROM_PERBASE_SAMCLIP(ch_input.bed, ch_input.header)

            PLOT_DEPTH_SAMCLIP(RENAMECHROM_PERBASE_SAMCLIP.out.tsv)
        
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
        PICARD_MARKDUPLICATES(ch_input.bam, ch_input.fasta, ch_input.fai)
        ch_bam_bai = PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai)
        
        SAMTOOLS_COVERAGE_RMDUP(ch_bam_bai)
        coverage_rmdup = SAMTOOLS_COVERAGE_RMDUP.out.coverage

        ch_input = ch_bam_bai.map{ it -> [it[0], it[1], 1] } //[meta, bam, scale]
        BEDTOOLS_GENOMECOV_PERBASE_RMDUP(ch_input, [], "bed")

        //PLOT_DEPTH_RMDUP(BEDTOOLS_GENOMECOV_PERBASE_RMDUP.out.genomecov)

        BEDTOOLS_GENOMECOV_PERBASE_RMDUP.out.genomecov
                .join(ref_header)
                .multiMap{ 
                    it -> 
                        bed: [it[0], it[1]] //[meta, bed]
                        header: [it[0], it[2]] //[meta, header] 
                }.set{ ch_input }
            //ch_input.header.view()
        RENAMECHROM_PERBASE_RMDUP(ch_input.bed, ch_input.header)

        PLOT_DEPTH_RMDUP(RENAMECHROM_PERBASE_RMDUP.out.tsv)
        
        
    emit:
        bam_bai = ch_bam_bai // bam file with secondary, supplementary, soft/hard clipped reads, duplicate reads removed
        coverage_rm_sec_sup
        //coverage_SAMCLIP
        coverage_rmdup
        versions = ch_versions
        
}