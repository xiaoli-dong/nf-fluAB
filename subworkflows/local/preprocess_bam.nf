
include {
    RENAMECHROM;
    PLOT_DEPTH;
} from '../../modules/local/misc'

include {
    SAMTOOLS_COVERAGE;
} from '../../modules/nf-core/samtools/coverage/main'

include {
    SAMTOOLS_VIEW;
} from '../../modules/nf-core/samtools/view/main'

include {
    BAM_MARKDUPLICATES_PICARD;
} from '../nf-core/bam_markduplicates_picard'

include {
    BEDTOOLS_GENOMECOV;
} from '../../modules/nf-core/bedtools/genomecov/main.nf'


/*
references:
    https://www.htslib.org/workflow/wgs-call.html
    https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00791-w
*/

workflow PREPROCESS_BAM {
    take:
        bam_bai //coordinate ordered bam
        fasta //[meta, fasta] reference used to produce bam_bai
        fai
        ref_header //[meta, txt] reference header

    main:
        ch_versions = Channel.empty()



        /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        remove  PCR and optical duplicates
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    SAMTOOLS_VIEW(bam_bai, fasta, [])
    ch_bam_bai = bam_bai
    SAMTOOLS_VIEW.out.bam
        .join(fasta)
        .join(fai)
        .multiMap{
            it ->
                bam_bai: [it[0], it[1]]
                fasta: [it[0], it[2]]
                fai: [it[0], it[3]]
        }
        .set{ch_input}


    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
    if(params.platform == 'illumina'){
        //all the duplicates will be removed from the output
        BAM_MARKDUPLICATES_PICARD(ch_input.bam_bai, ch_input.fasta, ch_input.fai)
        ch_bam_bai = BAM_MARKDUPLICATES_PICARD.out.bam.join( BAM_MARKDUPLICATES_PICARD.out.bai)
    }
    SAMTOOLS_COVERAGE(ch_bam_bai)
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions)

    ch_input = ch_bam_bai.map{ it -> [it[0], it[1], 1] } //[meta, bam, scale]
    BEDTOOLS_GENOMECOV(ch_input, [], "bed")
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)
    //PLOT_DEPTH(BEDTOOLS_GENOMECOV.out.genomecov)

    BEDTOOLS_GENOMECOV.out.genomecov
            .join(ref_header)
            .multiMap{
                it ->
                    bed: [it[0], it[1]] //[meta, bed]
                    header: [it[0], it[2]] //[meta, header]
            }.set{ ch_input }
        //ch_input.header.view()
    RENAMECHROM(ch_input.bed, ch_input.header)
    ch_versions = ch_versions.mix(RENAMECHROM.out.versions)
    PLOT_DEPTH(RENAMECHROM.out.tsv)
    ch_versions = ch_versions.mix(PLOT_DEPTH.out.versions)

emit:
    bam_bai = ch_bam_bai // bam file with secondary, supplementary, soft/hard clipped reads, duplicate reads removed
    coverage = SAMTOOLS_COVERAGE.out.coverage
    versions = ch_versions

}
