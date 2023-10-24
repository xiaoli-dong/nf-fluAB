/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowIllumina.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
def checkPathParamList = [ params.input]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/* ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
 */
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include {RUN_ILLUMINA_QC} from '../subworkflows/local/qc_illumina'

include {FETCH_REFERENCES} from '../subworkflows/local/fetch_references'
include {MAPPING} from '../subworkflows/local/mapping'
include {callConsensus} from '../subworkflows/local/callConsensus'
include {RUN_TYPING} from '../subworkflows/local/run_typing'
include {callClade} from '../subworkflows/local/callClade'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include {
    CSVTK_CONCAT as CONCAT_STATS_SHORT_RAW;
    CSVTK_CONCAT as CONCAT_STATS_SHORT_QC;
    CSVTK_CONCAT as CONCAT_STATS_SUMMARY;
    
} from '../modules/nf-core/csvtk/concat/main'

include {MAPPING_REPORT} from '../modules/local/mapping_report'
include {REPORT} from '../modules/local/report'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []


workflow ILLUMINA {

    ch_versions = Channel.empty()
    Channel
        .value(file( "${params.flu_db_msh}" ))
        .set { ch_flu_db_msh}
    Channel
        .value(file( "${params.flu_db_fasta}" ))
        .set { ch_flu_db_fasta}

    Channel
    .value(file( "${params.exclude}" ))
    .set { ch_exclude_list}

    Channel
        .value(file( "${params.typingdb}" ))
        .set { ch_typingdb}

    
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    short_reads = INPUT_CHECK.out.reads
    in_format = "tsv"
    out_format = "tsv"

     if(!params.skip_short_reads_qc){
        RUN_ILLUMINA_QC(short_reads)
        ch_versions = ch_versions.mix(RUN_ILLUMINA_QC.out.versions)
        //RUN_ILLUMINA_QC.out.qc_reads.view()
        
        //get rid of zero size contig file and avoid the downstream crash
        RUN_ILLUMINA_QC.out.qc_reads
            .filter {meta, reads -> reads[0].size() > 0 && reads[0].countFastq() > 0}
            .set { short_reads }

        
        CONCAT_STATS_SHORT_RAW(RUN_ILLUMINA_QC.out.raw_stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "short_reads_raw_seqstats"], files)}, in_format, out_format )
        CONCAT_STATS_SHORT_QC(RUN_ILLUMINA_QC.out.qc_stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "short_reads_qc_seqstats"], files)}, in_format, out_format )
    }

    FETCH_REFERENCES(short_reads, ch_flu_db_msh, ch_flu_db_fasta, ch_exclude_list)
    ch_versions.mix(FETCH_REFERENCES.out.versions)
    ref_fasta = FETCH_REFERENCES.out.fasta

    MAPPING(short_reads, ref_fasta)
    ch_versions.mix(MAPPING.out.versions)

   
    FETCH_REFERENCES.out.screen.join(MAPPING.out.coverage).multiMap{
        it ->
            screen: [it[0], it[1]]
            cov: [it[0], it[2]]
    }.set{
        input_mapping_report_ch
    }
    MAPPING_REPORT(input_mapping_report_ch.screen, input_mapping_report_ch.cov)
    //MAPPING_REPORT.out.csv.view()

    callConsensus(MAPPING.out.bam, MAPPING.out.bai, MAPPING.out.ref, MAPPING.out.ref_fai)
    ch_versions.mix(callConsensus.out.versions)
    //callConsensus.out.stats.view()
    
    RUN_TYPING(callConsensus.out.fasta, ch_typingdb)
    ch_versions.mix(RUN_TYPING.out.versions)

    //RUN_TYPING.out.tsv.view()
    callClade(RUN_TYPING.out.fasta, RUN_TYPING.out.tsv)
    //callClade.out.tsv.view()
    ch_test = callClade.out.tsv.groupTuple()//.view()
    ch_versions.mix(callClade.out.versions)
    ch_master = RUN_ILLUMINA_QC.out.raw_stats
        .join(RUN_ILLUMINA_QC.out.qc_stats)
        .join(callConsensus.out.stats)
        .join(RUN_TYPING.out.tsv)
        .join(MAPPING_REPORT.out.csv)
        .join(ch_test).view()
    
   
    //ch_master = RUN_ILLUMINA_QC.out.raw_stats.join(RUN_ILLUMINA_QC.out.qc_stats).join(callConsensus.out.stats).join(RUN_TYPING.out.tsv).join(MAPPING_REPORT.out.csv).join(callClade.out.tsv)//.view()
    //ch_master.view()
   
    //s11 has empty typing output, so there is no report produced .....
    REPORT(ch_master)
    CONCAT_STATS_SUMMARY(REPORT.out.oneline_csv.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "summary_report"], files)}, in_format, out_format )
    
    //ch_versions.view()
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    /* //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowIllumina.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowIllumina.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions) */
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
