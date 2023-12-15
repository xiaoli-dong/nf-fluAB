
include {
    FASTQC as FASTQC_INPUT;
    FASTQC as FASTQC_QC;
} from     '../../modules/nf-core/fastqc/main'

include {BBMAP_BBDUK} from      '../../modules/nf-core/bbmap/bbduk/main'
include {FASTP} from '../../modules/nf-core/fastp/main'
include {HOSTILE_BOWTIE2} from '../../modules/local/hostile/bowtie2/main'
include {
    SEQKIT_STATS as SEQKIT_STATS_INPUT;
    SEQKIT_STATS as SEQKIT_STATS_TRIMMED;
    SEQKIT_STATS as SEQKIT_STATS_HOSTIE;
} from '../../modules/local/seqkit/stats/main'

include {
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_INPUT;
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_QC;
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_HOSTIE;
} from '../../modules/nf-core/csvtk/concat/main'

workflow QC_ILLUMINA {   

    take:
        reads
        adapter_fasta
        hostile_human_ref_bowtie2
    main:
        ch_versions = Channel.empty()
        in_format = "tsv"
        out_format = "tsv"
        trimmed_reads = reads
        qc_reads = reads

        FASTQC_INPUT(reads)
        ch_versions = ch_versions.mix(FASTQC_INPUT.out.versions.first())
        SEQKIT_STATS_INPUT(reads)
        ch_versions = ch_versions.mix(SEQKIT_STATS_INPUT.out.versions.first())
        CSVTK_CONCAT_STATS_INPUT(SEQKIT_STATS_INPUT.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "illumina_reads.input_seqstats"], files)}, in_format, out_format )
        
        //default
        if ( params.illumina_reads_qc_tool == 'bbduk' ){
            BBMAP_BBDUK(reads, adapter_fasta)
            ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
            //get rid of zero size contig file and avoid the downstream crash
            /* BBMAP_BBDUK.out.reads
                .filter { meta, reads -> reads.countFastq() > 0 }
                .set { qc_reads } */
            trimmed_reads = BBMAP_BBDUK.out.reads
            SEQKIT_STATS_TRIMMED(trimmed_reads)
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
            SEQKIT_STATS_TRIMMED(trimmed_reads)
        }
        CSVTK_CONCAT_STATS_QC(SEQKIT_STATS_TRIMMED.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "illumina_reads.trimmed_seqstats"], files)}, in_format, out_format )

        HOSTILE_BOWTIE2(trimmed_reads, hostile_human_ref_bowtie2)
        FASTQC_QC(HOSTILE_BOWTIE2.out.reads)
        SEQKIT_STATS_HOSTIE(HOSTILE_BOWTIE2.out.reads)
        CSVTK_CONCAT_STATS_HOSTIE(SEQKIT_STATS_HOSTIE.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "illumina_reads.dehost_seqstats"], files)}, in_format, out_format )
        
        
    emit:
        qc_reads = HOSTILE_BOWTIE2.out.reads
        input_stats = SEQKIT_STATS_INPUT.out.stats
        qc_stats = SEQKIT_STATS_HOSTIE.out.stats
        versions = ch_versions
        
}
