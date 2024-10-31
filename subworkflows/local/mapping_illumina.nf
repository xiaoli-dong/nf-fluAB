include {   
    BWAMEM2_INDEX       
} from '../../modules/local/bwamem2/index/main'

include {   
    BWAMEM2_MEM         
} from '../../modules/local/bwamem2/mem/main'

include {   
    MINIMAP2_ALIGN      
} from '../../modules/local/minimap2/align/main'

include {   
    SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_MAPPING;
} from '../../modules/nf-core/samtools/coverage/main'

include {  
    SAMTOOLS_FIXMATE       
} from '../../modules/local/samtools/fixmate/main'

include {  
    SAMTOOLS_SORT       
} from '../../modules/local/samtools/sort/main'

workflow MAPPING_ILLUMINA {   

    take:
        reads //short reads
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
            
            //MINIMAP2_ALIGN(reads, fasta, sam_format)
            ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())
            sam = MINIMAP2_ALIGN.out.sam
        }
        else if(params.mapping_tool == 'bwa'){
           
            BWAMEM2_INDEX(fasta)
            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())

            reads.join(BWAMEM2_INDEX.out.index).multiMap{
                it ->
                    reads: [it[0], it[1]]
                    index: [it[0], it[2]]
            }.set{ ch_input }

            BWAMEM2_MEM ( ch_input.reads, ch_input.index) 
            //BWAMEM2_MEM ( reads, BWAMEM2_INDEX.out.index)
            ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
            sam = BWAMEM2_MEM.out.sam
        }
        /*
        Because BWA can sometimes leave unusual FLAG information on SAM records, 
        it is helpful when working with many tools to first clean up read pairing 
        information and flags:
        */
        SAMTOOLS_FIXMATE(sam)
        ch_versions = ch_versions.mix(SAMTOOLS_FIXMATE.out.versions.first())

        //sort from name order to coordinate order
        SAMTOOLS_SORT (SAMTOOLS_FIXMATE.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())
        
        /* SAMTOOLS_COVERAGE_MAPPING(SAMTOOLS_SORT.out.bam_bai)
        ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE_MAPPING.out.versions.first())
 */
        
    emit:
        bam_bai = SAMTOOLS_SORT.out.bam_bai
        //coverage = SAMTOOLS_COVERAGE_MAPPING.out.coverage
        versions = ch_versions
        
}