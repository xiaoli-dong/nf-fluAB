
include {fetchRefs_illumina} from './fetchRefs_illumina'
include {MAPPING_ILLUMINA} from './mapping_illumina'
include {callConsensus_illumina} from './callConsensus_illumina'
include {MAPPING_REPORT} from '../../modules/local/mapping_report'

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

        MAPPING_ILLUMINA(illumina_reads, fasta)
        coverage = MAPPING_ILLUMINA.out.coverage
        ch_versions.mix(MAPPING_ILLUMINA.out.versions)

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

        callConsensus_illumina(MAPPING_ILLUMINA.out.bam, MAPPING_ILLUMINA.out.bai, MAPPING_ILLUMINA.out.fasta, MAPPING_ILLUMINA.out.fasta_fai)
        ch_versions.mix(callConsensus_illumina.out.versions)
        
    emit:
        consensus = callConsensus_illumina.out.fasta
        consensus_stats = callConsensus_illumina.out.stats
        mapping_summary = MAPPING_REPORT.out.csv 
        versions = ch_versions
        
}