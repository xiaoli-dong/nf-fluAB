include {
    BEDTOOLS_MASKFASTA
} from '../../modules/nf-core/bedtools/maskfasta/main.nf'

include {
    BEDTOOLS_GENOMECOV
} from '../../modules/nf-core/bedtools/genomecov/main.nf'

include {
    SAMTOOLS_FAIDX as  SAMTOOLS_FAIDX_MASK;
} from '../../modules/nf-core/samtools/faidx/main'

//mask fasta sequences for those low depth region < mindepth cutoff
workflow MASK_FASTA {   

    take:
        bam_bai //[meta, bam, bai]
        fasta // [meta, fasta]
    main:
        ch_versions = Channel.empty()

        ch_input = bam_bai.map{
            it -> [it[0], it[1], 1] //[meta, bam, scale]
        }
        //mask low depth regions of the reference
        BEDTOOLS_GENOMECOV(ch_input, [], "bed")
        ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)

        BEDTOOLS_GENOMECOV.out.genomecov.join(fasta).multiMap{
            it ->
                bed: [it[0], it[1]]
                fasta: it[2]
        }.set{
            ch_input
        }

        BEDTOOLS_MASKFASTA(ch_input.bed, ch_input.fasta)
        ch_versions.mix(BEDTOOLS_MASKFASTA.out.versions)

        SAMTOOLS_FAIDX_MASK(BEDTOOLS_MASKFASTA.out.fasta)
        ch_versions.mix(SAMTOOLS_FAIDX_MASK.out.versions)

        fasta_fai = BEDTOOLS_MASKFASTA.out.fasta.join(SAMTOOLS_FAIDX_MASK.out.fai)


    emit:
        fasta_fai
        versions = ch_versions
        
}

//produce bed file for those low depth region < mindepth cutoff
workflow LOW_DEPTH_BED_FROM_BAM {   

    take:
        bam_bai

    main:
        ch_versions = Channel.empty()
        ch_input = bam_bai.map{
            it -> [it[0], it[1], 1] //[meta, bam, scale]
        }
        //bed file of the low depth regions < params.mindepth
        BEDTOOLS_GENOMECOV(ch_input, [], "bed")
        ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)

    emit:
        bed = BEDTOOLS_GENOMECOV.out.genomecov
        versions = ch_versions
        
}
