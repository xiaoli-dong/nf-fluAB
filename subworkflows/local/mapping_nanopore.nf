include { 
    MINIMAP2_ALIGN
} from '../../modules/local/minimap2/align/main'

include {
    SAMTOOLS_INDEX
} from '../../modules/nf-core/samtools/index/main'

include {  
    SAMTOOLS_FIXMATE       
} from '../../modules/local/samtools/fixmate/main'

include {
    SAMTOOLS_SORT
} from '../../modules/local/samtools/sort/main'

include {   
    SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_MAPPING;
} from '../../modules/nf-core/samtools/coverage/main'

workflow MAPPING_NANOPORE {   

    take:
        reads
        fasta //reference [meta, fasta]

    main:
        ch_versions = Channel.empty()
        
        sam = Channel.empty()
        
        if ( params.mapping_tool == 'minimap2' ){
            sam_format = true
           
            reads.join(fasta).multiMap{
                it ->
                    reads: [it[0], it[1]]
                    fasta: it[2]
            }.set{ ch_input }

            MINIMAP2_ALIGN(ch_input.reads, ch_input.fasta, sam_format)
            ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
            sam = MINIMAP2_ALIGN.out.sam
            
        }

        /*
        Because minimap2 can sometimes leave unusual FLAG information on SAM records, 
        it is helpful when working with many tools to first clean up read pairing 
        information and flags:
        */
        SAMTOOLS_FIXMATE(sam)
        ch_versions = ch_versions.mix(SAMTOOLS_FIXMATE.out.versions)
       
        SAMTOOLS_SORT (SAMTOOLS_FIXMATE.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
        
        // SAMTOOLS_INDEX (SAMTOOLS_SORT.out.bam)
        // ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
        // bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)
        
        SAMTOOLS_COVERAGE_MAPPING(SAMTOOLS_SORT.out.bam_bai)
        ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE_MAPPING.out.versions)
        
        
    emit:
        bam_bai = SAMTOOLS_SORT.out.bam_bai
        coverage = SAMTOOLS_COVERAGE_MAPPING.out.coverage
        versions = ch_versions
}