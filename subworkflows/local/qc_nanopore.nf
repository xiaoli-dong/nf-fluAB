//nanopore
include {
    PORECHOP_PORECHOP
}  from '../../modules/nf-core/porechop/porechop/main.nf'

include {
    CHOPPER
} from '../../modules/local/chopper/main.nf'

include {
    HOSTILE_CLEAN
} from '../../modules/local/hostile/clean'

include {
    SEQKIT_STATS as INPUT_STATS
    SEQKIT_STATS as PORECHOP_STATS;
    SEQKIT_STATS as CHOPPER_STATS;
    SEQKIT_STATS as HOSTILE_STATS;
} from '../../modules/local/seqkit/stats/main'

include {
    CSVTK_CONCAT as CONCAT_INPUT_STATS;
    CSVTK_CONCAT as CONCAT_PORECHOP_STATS;
    CSVTK_CONCAT as CONCAT_CHOPPER_STATS;
    CSVTK_CONCAT as CONCAT_HOSTILE_STATS;
} from '../../modules/nf-core/csvtk/concat/main'
//fastplong
//--do_read_splitting, --detect_mid_strand_adapter, --detect_mid_strand_barcodes, --require_barcodes_both_end
workflow QC_NANOPORE {

    take:
        reads
        //hostile_ref
    main:

        ch_versions = Channel.empty()
        INPUT_STATS(reads)
        ch_versions = ch_versions.mix(INPUT_STATS.out.versions)

        // QC
        PORECHOP_PORECHOP(reads)
        ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions)

        PORECHOP_PORECHOP.out.reads
                .filter {meta, reads -> reads.size() > 0 && reads.countFastq() > 0}
                .set { trimmed_reads }
        //trimmed_reads = PORECHOP_PORECHOP.out.reads
        PORECHOP_STATS(trimmed_reads)
        ch_versions = ch_versions.mix(PORECHOP_STATS.out.versions)

        CHOPPER(trimmed_reads)
        ch_versions = ch_versions.mix(CHOPPER.out.versions)


        CHOPPER.out.fastq
            .filter {meta, fastq -> fastq.size() > 0 && fastq.countFastq() > 0}
            //.filter { meta, reads -> reads[0].countFastq() > 0 }
            .set { filtered_reads }
        CHOPPER_STATS(filtered_reads)

        //HOSTILE(filtered_reads, "minimap2", hostile_ref)
        HOSTILE_CLEAN(filtered_reads, [params.hostile_ref_name_nanopore, params.hostile_ref_dir])
        ch_versions = ch_versions.mix(HOSTILE_CLEAN.out.versions)
        HOSTILE_STATS(HOSTILE_CLEAN.out.fastq)

        qc_reads = HOSTILE_CLEAN.out.fastq //gzip compressed


        in_format = "tsv"
        out_format = "tsv"
        CONCAT_INPUT_STATS(INPUT_STATS.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads.input_seqstats"], files)}, in_format, out_format )
        CONCAT_PORECHOP_STATS(PORECHOP_STATS.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads.porechop_seqstats"], files)}, in_format, out_format )
        CONCAT_CHOPPER_STATS(CHOPPER_STATS.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads.chopper_seqstats"], files)}, in_format, out_format )
        CONCAT_HOSTILE_STATS(HOSTILE_STATS.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads.dehost_seqstats"], files)}, in_format, out_format )

    emit:

        input_stats = INPUT_STATS.out.stats
        qc_stats = HOSTILE_STATS.out.stats
        qc_reads
        versions = ch_versions

}
