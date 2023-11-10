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
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include {ILLUMINA_QC} from '../subworkflows/local/qc_illumina'
include {
    classifier_blast;
    classifier_nextclade;
} from '../subworkflows/local/classifier'

//include {callClade} from '../subworkflows/local/callClade'
include { ASSEMBLY_ILLUMINA } from '../subworkflows/local/assembly_illumina'
include { PREPARE_REFERENCES          } from '../subworkflows/local/prepare_references'
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
    CSVTK_CONCAT as CONCAT_NEXTCLADE;
    
} from '../modules/nf-core/csvtk/concat/main'


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
    
    // SUBWORKFLOW: prepare reference databases ...
    //
    PREPARE_REFERENCES ()
    ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.versions)
    //PREPARE_REFERENCES.out.ch_flu_db_fasta.view()
    //PREPARE_REFERENCES.out.ch_flu_db_msh.view()
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    illumina_reads = INPUT_CHECK.out.reads
    in_format = "tsv"
    out_format = "tsv"

     if(!params.skip_illumina_reads_qc){
        ILLUMINA_QC(illumina_reads, PREPARE_REFERENCES.out.ch_flu_primers)
        ch_versions = ch_versions.mix(ILLUMINA_QC.out.versions)
        //ILLUMINA_QC.out.qc_reads.view()
        
        //get rid of zero size contig file and avoid the downstream crash
        ILLUMINA_QC.out.qc_reads
            .filter {meta, reads -> reads[0].size() > 0 && reads[0].countFastq() > 0}
            .set { illumina_reads } 
        CONCAT_STATS_SHORT_RAW(
            ILLUMINA_QC.out.raw_stats
                .map { cfg, stats -> stats }.collect()
                .map { files -> tuple([id: "illumina_reads_raw.seqstats"], files)}, 
                in_format, 
                out_format
        )
        CONCAT_STATS_SHORT_QC(
            ILLUMINA_QC.out.qc_stats
                .map { cfg, stats -> stats }.collect()
                .map { files -> tuple([id: "illumina_reads_qc.seqstats"], files)}, 
                in_format, 
                out_format 
        )
    }

    ASSEMBLY_ILLUMINA(
        illumina_reads, 
        PREPARE_REFERENCES.out.ch_flu_db_msh, 
        PREPARE_REFERENCES.out.ch_flu_db_fasta
    )
    
    //classifier_blast(ASSEMBLY_ILLUMINA.out.consensus, PREPARE_REFERENCES.out.ch_typing_db)
    classifier_blast(ASSEMBLY_ILLUMINA.out.blastn_fasta, PREPARE_REFERENCES.out.ch_typing_db)
    ch_versions.mix(classifier_blast.out.versions)

   
    ASSEMBLY_ILLUMINA.out.blastn_fasta.join(classifier_blast.out.tsv).multiMap{
        it ->
            consensus: [it[0], it[1]]
            tsv: [it[0], it[2]]
    }.set{
        ch_input
    }
    classifier_nextclade(ch_input.consensus, ch_input.tsv)
    //classifier_nextclade.out.tsv.collectFile().view()

    ch_nextclade = classifier_nextclade.out.tsv.groupTuple()//.view()
    ch_versions.mix(classifier_nextclade.out.versions)
    
    //reporting
    ch_report = ILLUMINA_QC.out.raw_stats
        .join(ILLUMINA_QC.out.qc_stats)
        .join(ASSEMBLY_ILLUMINA.out.consensus_stats)
        .join(classifier_blast.out.tsv)
        .join(ASSEMBLY_ILLUMINA.out.coverage_summary)
        .join(ch_nextclade)//.view()
        .multiMap{
            it ->
                raw_stats: [it[0], it[1]]
                qc_stats: [it[0], it[2]]
                consensus_stats: [it[0], it[3]]
                blastn_outfmt6: [it[0], it[4]]
                coverage_summary: [it[0], it[5]]
                nextclade_csv: it[6].join(',')
        } 
    
    
    REPORT(
        ch_report.raw_stats,
        ch_report.qc_stats,
        ch_report.consensus_stats,
        ch_report.blastn_outfmt6,
        ch_report.coverage_summary,
        ch_report.nextclade_csv
    ) 

    CONCAT_STATS_SUMMARY(
        REPORT.out.oneline_csv.map { cfg, stats -> stats }.collect()
            .map { files -> tuple([id: "summary_report"], files)}, 
        in_format, 
        out_format 
    )
    
    //ch_versions.view()
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

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
