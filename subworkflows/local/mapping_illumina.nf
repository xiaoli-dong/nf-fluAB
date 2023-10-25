
include {SAMTOOLS_FAIDX} from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_SORT} from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX} from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_DEPTH} from '../../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_COVERAGE} from '../../modules/nf-core/samtools/coverage/main'
include { BWAMEM2_INDEX} from '../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM} from '../../modules/nf-core/bwamem2/mem/main'


workflow MAPPING_ILLUMINA {   

    take:
        illumina_reads
        fasta
    main:
        ch_versions = Channel.empty()

        SAMTOOLS_FAIDX(fasta)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())
        
        BWAMEM2_INDEX(fasta)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())

        illumina_reads.join(BWAMEM2_INDEX.out.index).multiMap{
            it ->
            reads: [it[0], it[1]]
            index: [it[0], it[2]]
        }.set{
            input_bwa_ch
        }

        //only keep primary alignment
        //bwa mem | samtools view -h -F260
        BWAMEM2_MEM (
            input_bwa_ch.reads,
            input_bwa_ch.index,
            false
        )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

        SAMTOOLS_SORT (
            BWAMEM2_MEM.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

        //SAMTOOLS_SORT.out.bam.view()

        SAMTOOLS_INDEX (
            SAMTOOLS_SORT.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())
        //SAMTOOLS_INDEX.out.bai.view()
        
        SAMTOOLS_DEPTH(SAMTOOLS_SORT.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

        input_ch = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)
        //input_ch.view()
        SAMTOOLS_COVERAGE(input_ch)
        
        
    emit:
        bam = SAMTOOLS_SORT.out.bam
        bai = SAMTOOLS_INDEX.out.bai
        fasta_fai = SAMTOOLS_FAIDX.out.fai
        fasta
        coverage = SAMTOOLS_COVERAGE.out.coverage
        versions = ch_versions
        
}