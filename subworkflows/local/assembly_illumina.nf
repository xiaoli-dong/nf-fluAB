
include {fetchRefs_illumina} from './fetchRefs_illumina'
include {process_bam} from './process_bam'
include { make_consensus } from './make_consensus'
include {MAPPING_REPORT} from '../../modules/local/mapping_report'
include { BWAMEM2_INDEX} from '../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM} from '../../modules/nf-core/bwamem2/mem/main'
include { MINIMAP2_ALIGN } from '../../modules/local/minimap2/align/main'
include {SAMTOOLS_FAIDX} from '../../modules/nf-core/samtools/faidx/main'
include { FREEBAYES } from '../../modules/nf-core/freebayes/main.nf'

workflow ASSEMBLY_ILLUMINA {   

    take:
        illumina_reads
        ch_flu_db_msh
        ch_flu_db_fasta

    main:
        ch_versions = Channel.empty()
        //illumina_reads.view()
        fetchRefs_illumina(illumina_reads, ch_flu_db_msh, ch_flu_db_fasta)
        fasta = fetchRefs_illumina.out.fasta
        screen = fetchRefs_illumina.out.screen
        ch_versions.mix(fetchRefs_illumina.out.versions)

        SAMTOOLS_FAIDX(fasta)
        fai = SAMTOOLS_FAIDX.out.fai

        
        bam = null
        

        if ( params.illumina_reads_mapping_tool == 'minimap2' ){
            bam_format = true
            cigar_paf_format = false
            cigar_bam = false
           
            illumina_reads.join(fasta).multiMap{
                it ->
                    reads: [it[0], it[1]]
                    fasta: it[2]
            }
            .set{
                ch_input
            }

            MINIMAP2_ALIGN(ch_input.reads, ch_input.fasta, bam_format, cigar_paf_format, cigar_bam)
            ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())
            bam = MINIMAP2_ALIGN.out.bam
        }
        else if(params.illumina_reads_mapping_tool == 'bwa'){
           
            BWAMEM2_INDEX(fasta)
            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())

            illumina_reads.join(BWAMEM2_INDEX.out.index).multiMap{
                it ->
                reads: [it[0], it[1]]
                index: [it[0], it[2]]
            }.set{
                ch_input
            }

            //only keep primary alignment
            //bwa mem | samtools view -h -F260
            BWAMEM2_MEM (
                ch_input.reads,
                ch_input.index,
                false
            )
            ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
            bam = BWAMEM2_MEM.out.bai
        }
        //MAPPING_ILLUMINA(illumina_reads, fasta)
        process_bam(bam)
        coverage = process_bam.out.coverage
        ch_versions.mix(process_bam.out.versions)

        screen.join(coverage)
            .multiMap {
                it ->
                    screen: [it[0], it[1]]
                    coverage: [it[0], it[2]]
            }.set{
                ch_input
            }
        MAPPING_REPORT(ch_input.screen, ch_input.coverage)
        MAPPING_REPORT.out.csv.view()

        //variant calling
        //todo: need join before calling to make sure fasta fail are all for the correspoing samples
        process_bam.out.bam
            .join(process_bam.out.bai)
            .join(fasta)
            .join(fai)
            .multiMap {
                it ->
                    mapping: [it[0], it[1], it[2], [], [], []]
                    fasta: it[3]
                    fai: it[4]
                    
            }.set{
                ch_input_freebayes
            }

        samples     = []
        populations = []
        cnv         = []
        FREEBAYES (
            ch_input_freebayes.mapping, 
            ch_input_freebayes.fasta, 
            ch_input_freebayes.fai, 
            samples, 
            populations, 
            cnv
        )

        //for making consensus
        FREEBAYES.out.vcf.join(fasta).multiMap{
                it ->
                    vcf: [it[0], it[1]]
                    fasta: [it[0], it[2]]
        }.set{
                ch_input
        }

        //callConsensus_illumina(ch_input.vcf, ch_input.fasta)
        make_consensus(ch_input.vcf, ch_input.fasta)
        ch_versions.mix(make_consensus.out.versions)
        
    emit:
        consensus = make_consensus.out.fasta
        consensus_stats = make_consensus.out.stats
        mapping_summary = MAPPING_REPORT.out.csv 
        versions = ch_versions
        
}