include { 
    MINIMAP2_ALIGN
} from '../../modules/local/minimap2/align/main'

include {
    SAMTOOLS_INDEX
} from '../../modules/nf-core/samtools/index/main'

include {
    SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_MAPPING;  
} from '../../modules/nf-core/samtools/coverage/main'

include {
    SAMTOOLS_SORT
} from '../../modules/nf-core/samtools/sort/main'

include {
    FILTERMASH;
} from '../../modules/local/misc'

workflow MAPPING_NANOPORE {   

    take:
        reads
        fasta //reference [meta, fasta]

    main:
        ch_versions = Channel.empty()
        
        sam = Channel.empty()
        
        if ( params.nanopore_reads_mapping_tool == 'minimap2' ){
            sam_format = true
           
            reads.join(fasta).multiMap{
                it ->
                    reads: [it[0], it[1]]
                    fasta: it[2]
            }.set{ ch_input }

            MINIMAP2_ALIGN(ch_input.reads, ch_input.fasta, sam_format)
            ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())
            sam = MINIMAP2_ALIGN.out.sam
            
        }
     
       
        SAMTOOLS_SORT (sam)
        SAMTOOLS_INDEX (SAMTOOLS_SORT.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
        
        bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)
        
        SAMTOOLS_COVERAGE_MAPPING(bam_bai)
        ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE_MAPPING.out.versions.first())
        
    emit:
        bam_bai 
        coverage = SAMTOOLS_COVERAGE_MAPPING.out.coverage
        versions = ch_versions
}