
include {
    BBMAP_BBDUK;
} from      '../../modules/nf-core/bbmap/bbduk/main'

include {
    FASTP
} from '../../modules/nf-core/fastp/main'

include {
    HOSTILE
} from '../../modules/local/hostile/main'

include {
    SEQKIT_STATS as INPUT_STATS;
    SEQKIT_STATS as TRIMMED_STATS;
    SEQKIT_STATS as HOSTILE_STATS;
} from '../../modules/local/seqkit/stats/main'

include {
    CSVTK_CONCAT as CONCAT_INPUT_STATS;
    CSVTK_CONCAT as CONCAT_TRIMMED_STATS;
    CSVTK_CONCAT as CONCAT_HOSTILE_STATS;
} from '../../modules/nf-core/csvtk/concat/main'

workflow QC_ILLUMINA {   

    take:
        reads
        adapter_fasta
        hostile_ref
    main:
        ch_versions = Channel.empty()
        trimmed_reads = Channel.empty()
        qc_reads = Channel.empty()

        //concat stats format
        in_format = "tsv"
        out_format = "tsv"

        INPUT_STATS(reads)
        ch_versions = ch_versions.mix(INPUT_STATS.out.versions)
        CONCAT_INPUT_STATS(INPUT_STATS.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "illumina_reads.input_seqstats"], files)}, in_format, out_format )
        
        //default
        if ( params.illumina_reads_qc_tool == 'bbduk' ){
            BBMAP_BBDUK(reads, adapter_fasta)
            ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions)

            //get rid of zero size contig file and avoid the downstream crash
            BBMAP_BBDUK.out.reads
                .filter {meta, reads -> reads[0].size() > 0 && reads[0].countFastq() > 0}
                .set { trimmed_reads }
        }
        else if ( params.illumina_reads_qc_tool == 'fastp'){
            save_trimmed_fail = false
            save_merged       = false
            FASTP ( reads, adapter_fasta, save_trimmed_fail, save_merged )
            ch_versions = ch_versions.mix(FASTP.out.versions)
           
            //get rid of zero size contig file and avoid the downstream crash
            FASTP.out.reads
                .filter {meta, reads -> reads[0].size() > 0 && reads[0].countFastq() > 0}
                .set { trimmed_reads }
        }

        TRIMMED_STATS(trimmed_reads)
        ch_versions = ch_versions.mix(TRIMMED_STATS.out.versions)

        CONCAT_TRIMMED_STATS(TRIMMED_STATS.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "illumina_reads.trimmed_seqstats"], files)}, in_format, out_format )
        ch_versions = ch_versions.mix(CONCAT_TRIMMED_STATS.out.versions)

        HOSTILE(trimmed_reads, "bowtie2", hostile_ref)
        HOSTILE.out.reads
            .filter {meta, reads -> reads[0].size() > 0 && reads[0].countFastq() > 0}
            .set { qc_reads }
        ch_versions = ch_versions.mix(HOSTILE.out.versions)

        HOSTILE_STATS(qc_reads)
        ch_versions = ch_versions.mix(HOSTILE_STATS.out.versions)

        CONCAT_HOSTILE_STATS(HOSTILE_STATS.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "illumina_reads.dehost_seqstats"], files)}, in_format, out_format )
        ch_versions = ch_versions.mix(CONCAT_HOSTILE_STATS.out.versions)
        
    emit:
        input_stats = INPUT_STATS.out.stats
        qc_reads = HOSTILE.out.reads
        qc_stats = HOSTILE_STATS.out.stats
        versions = ch_versions
}
