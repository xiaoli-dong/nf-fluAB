//nanopore
include {
    PORECHOP_PORECHOP
}  from '../../modules/nf-core/porechop/porechop/main.nf'  

include {
    CHOPPER
} from '../../modules/local/chopper/main.nf'

include {
    HOSTILE
} from '../../modules/local/hostile/main'

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

workflow QC_NANOPORE {

    take:
        reads
        hostile_ref
    main:

        ch_versions = Channel.empty()
        INPUT_STATS(reads)
        ch_versions = ch_versions.mix(INPUT_STATS.out.versions.first())
        
        // QC
        PORECHOP_PORECHOP(reads)
        ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions.first())
        PORECHOP_STATS(PORECHOP_PORECHOP.out.reads)
        ch_versions = ch_versions.mix(PORECHOP_STATS.out.versions.first())
        
        CHOPPER(PORECHOP_PORECHOP.out.reads)
        ch_versions = ch_versions.mix(CHOPPER.out.versions.first())
        CHOPPER_STATS(CHOPPER.out.fastq)

        HOSTILE(CHOPPER.out.fastq, "minimap2", hostile_ref)
        ch_versions = ch_versions.mix(HOSTILE.out.versions.first())
        HOSTILE_STATS(HOSTILE.out.reads)

        qc_reads = HOSTILE.out.reads //gzip compressed
               
        
        in_format = "tsv"
        out_format = "tsv"
        CONCAT_INPUT_STATS(INPUT_STATS.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads.input_seqstats"], files)}, in_format, out_format )
        CONCAT_PORECHOP_STATS(PORECHOP_STATS.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads.porechop_seqstats"], files)}, in_format, out_format ) 
        CONCAT_CHOPPER_STATS(CHOPPER_STATS.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads.chopper_seqstats"], files)}, in_format, out_format ) 
        CONCAT_HOSTILE_STATS(HOSTILE_STATS.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads.dehost_seqstats"], files)}, in_format, out_format ) 

    emit:
        qc_reads
        input_stats = INPUT_STATS.out.stats
        qc_stats = HOSTILE_STATS.out.stats
        versions = ch_versions

}
