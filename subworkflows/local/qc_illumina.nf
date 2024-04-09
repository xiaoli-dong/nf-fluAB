
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
        in_format = "tsv"
        out_format = "tsv"
        trimmed_reads = reads
        qc_reads = reads

        // FASTQC_INPUT(reads)
        // ch_versions = ch_versions.mix(FASTQC_INPUT.out.versions.first())
        INPUT_STATS(reads)
        ch_versions = ch_versions.mix(INPUT_STATS.out.versions.first())
        CONCAT_INPUT_STATS(INPUT_STATS.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "illumina_reads.input_seqstats"], files)}, in_format, out_format )
        
        //default
        if ( params.illumina_reads_qc_tool == 'bbduk' ){
            BBMAP_BBDUK(reads, adapter_fasta)
            ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
            //get rid of zero size contig file and avoid the downstream crash
            /* BBMAP_BBDUK.out.reads
                .filter { meta, reads -> reads.countFastq() > 0 }
                .set { qc_reads } */
            trimmed_reads = BBMAP_BBDUK.out.reads
            //FASTQC_QC(trimmed_reads)
           TRIMMED_STATS(trimmed_reads)
        }
        else if ( params.illumina_reads_qc_tool == 'fastp'){
            save_trimmed_fail = false
            save_merged       = false
            FASTP ( reads, adapter_fasta, save_trimmed_fail, save_merged )
            ch_versions = ch_versions.mix(FASTP.out.versions.first())
            //get rid of zero size contig file and avoid the downstream crash
            FASTP.out.reads
                .filter {meta, reads -> reads[0].size() > 0 && reads[0].countFastq() > 0}
                //.filter { meta, reads -> reads[0].countFastq() > 0 }
                .set { trimmed_reads }
            trimmed_reads = FASTP.out.reads
            //FASTQC_QC(trimmed_reads)
           TRIMMED_STATS(trimmed_reads)
        }
        CONCAT_TRIMMED_STATS(TRIMMED_STATS.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "illumina_reads.trimmed_seqstats"], files)}, in_format, out_format )

        HOSTILE(trimmed_reads, "bowtie2", hostile_ref)
        //FASTQC_QC_HOSTILE(HOSTILE.out.reads)
        HOSTILE_STATS(HOSTILE.out.reads)
        CONCAT_HOSTILE_STATS(HOSTILE_STATS.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "illumina_reads.dehost_seqstats"], files)}, in_format, out_format )
        
        
    emit:
        qc_reads = HOSTILE.out.reads
        input_stats = INPUT_STATS.out.stats
        qc_stats = HOSTILE_STATS.out.stats
        versions = ch_versions
        
}
