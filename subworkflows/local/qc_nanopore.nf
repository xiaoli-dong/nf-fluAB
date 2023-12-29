//nanopore
include {
    NANOPLOT as NANOPLOT_INPUT;
    NANOPLOT as NANOPLOT_QC;
} from '../../modules/nf-core/nanoplot/main' 

include {PORECHOP_PORECHOP}  from '../../modules/nf-core/porechop/porechop/main.nf'  
include {CHOPPER} from '../../modules/local/chopper/main.nf'
include {HOSTILE} from '../../modules/local/hostile/main'
include {
    SEQKIT_STATS as SEQKIT_STATS_INPUT
    SEQKIT_STATS as SEQKIT_STATS_PORECHOP;
    SEQKIT_STATS as SEQKIT_STATS_CHOPPER;
    SEQKIT_STATS as SEQKIT_STATS_HOSTILE;
} from '../../modules/local/seqkit/stats/main'

include {
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_INPUT;
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_PORECHOP;
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_CHOPPER;
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_HOSTILE;
} from '../../modules/nf-core/csvtk/concat/main'

workflow QC_NANOPORE {

    take:
        reads
        hostile_human_ref
    main:

        ch_versions = Channel.empty()
        //reads.view()
        NANOPLOT_INPUT(reads)
        ch_versions = ch_versions.mix(NANOPLOT_INPUT.out.versions.first())
        //reads.view()
        SEQKIT_STATS_INPUT(reads)
        ch_versions = ch_versions.mix(SEQKIT_STATS_INPUT.out.versions.first())
        
        // QC
        PORECHOP_PORECHOP(reads)
        SEQKIT_STATS_PORECHOP(PORECHOP_PORECHOP.out.reads)
        ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions.first())
        
        CHOPPER(PORECHOP_PORECHOP.out.reads)
        ch_versions = ch_versions.mix(CHOPPER.out.versions.first())
        SEQKIT_STATS_CHOPPER(CHOPPER.out.fastq)

        HOSTILE(CHOPPER.out.fastq, "minimap2", hostile_human_ref)
        ch_versions = ch_versions.mix(HOSTILE.out.versions.first())
        SEQKIT_STATS_HOSTILE(HOSTILE.out.reads)

        qc_reads = HOSTILE.out.reads //gzip compressed
               
        NANOPLOT_QC(qc_reads)
        
        in_format = "tsv"
        out_format = "tsv"
        CSVTK_CONCAT_STATS_INPUT(SEQKIT_STATS_INPUT.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads.input_seqstats"], files)}, in_format, out_format )
        CSVTK_CONCAT_STATS_PORECHOP(SEQKIT_STATS_PORECHOP.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads.porechop_seqstats"], files)}, in_format, out_format ) 
        CSVTK_CONCAT_STATS_CHOPPER(SEQKIT_STATS_CHOPPER.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads.chopper_seqstats"], files)}, in_format, out_format ) 
        CSVTK_CONCAT_STATS_HOSTILE(SEQKIT_STATS_HOSTILE.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads.dehost_seqstats"], files)}, in_format, out_format ) 

    emit:
        qc_reads
        input_stats = SEQKIT_STATS_INPUT.out.stats
        qc_stats = SEQKIT_STATS_HOSTILE.out.stats
        versions = ch_versions

}
