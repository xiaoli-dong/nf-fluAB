include {
    SAMTOOLS_VIEW as SAMTOOLS_FILTER;
    SAMTOOLS_VIEW as SAMTOOLS_VIEW_BAM2SAM;
} from '../../modules/local/samtools/view/main'

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
    SAMTOOLS_SORT as SAMTOOLS_SORT_SAMCLIP;
} from '../../modules/local/samtools/sort/main'

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

//improvement

workflow PREPROCESS_BAM {   
        take:
            bam_bai //coordinate ordered bam
            fasta_fai //[meta, fasta, fai] reference used to produce bam_bai
            ref_header //[meta, txt] reference header

        main:
            ch_versions = Channel.empty()
            /*
            https://www.htslib.org/workflow/wgs-call.html
            */
               
        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        exclude alignments:
            read unmapped (0x4)
            mate unmapped (0x8)
            not primary alignment (0x100)
            supplementary alignment (0x800) 
            samtools view ext.args = " -h -F 2316 -b"
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        SAMTOOLS_FILTER(
            bam_bai, 
            fasta_fai.map{
                meta, fasta, fai ->[meta, fasta]
                }, 
            []
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FILTER.out.versions.first())
        SAMTOOLS_FILTER.out.bam.view()
        SAMTOOLS_FILTER.out.bai.view()
        ch_bam_bai = SAMTOOLS_FILTER.out.bam.join(SAMTOOLS_FILTER.out.bai).view()

         /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           remove  PCR and optical duplicates
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
        coverage = SAMTOOLS_COVERAGE_RMDUP.out.coverage
        versions = ch_versions
        
}