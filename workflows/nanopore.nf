/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def valid_params = [
    nanopore_reads_mapping_tools : ['minimap2']
]
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowNanopore.initialise(params, log, valid_params)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
//def checkPathParamList = [ params.input]
def checkPathParamList = [ 
    params.input, 
    params.hostile_human_ref_minimap2, 
    params.flu_primers, 
    params.typing_db, 
    params.flu_db_msh, 
    params.flu_db_fasta,
    params.clair3_variant_model,
    params.nextclade_dataset_base
]
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
include {
    INPUT_CHECK
} from '../subworkflows/local/input_check'

include {
    PREPARE_REFERENCES
} from '../subworkflows/local/prepare_references'

include {
    QC_NANOPORE
} from '../subworkflows/local/qc_nanopore'

include {
    GET_REF_BY_MASH
} from '../subworkflows/local/get_references'

include {
    MAPPING_NANOPORE
} from '../subworkflows/local/mapping_nanopore'

include {
    PREPROCESS_BAM
} from '../subworkflows/local/preprocess_bam'

include {
    MAPPING_SUMMARY;
} from '../modules/local/misc'

include {
    CLASSIFIER_BLAST;
    CLASSIFIER_NEXTCLADE;
} from '../subworkflows/local/classifier'

include {
    BAM2LOW_DEPTH_BED;
} from '../subworkflows/local/mask'

include {
    VARIANTS_NANOPORE;
} from '../subworkflows/local/variants_nanopore'

include {
    CONSENSUS
} from '../subworkflows/local/consensus'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include {
    CUSTOM_DUMPSOFTWAREVERSIONS
} from '../modules/nf-core/custom/dumpsoftwareversions/main'

include {
    CSVTK_CONCAT as CONCAT_STATS_SHORT_RAW;
    CSVTK_CONCAT as CONCAT_STATS_SHORT_QC;
    CSVTK_CONCAT as CONCAT_STATS_SUMMARY;
    CSVTK_CONCAT as CONCAT_NEXTCLADE;
} from '../modules/nf-core/csvtk/concat/main'

include {
    FREEBAYES
} from '../modules/nf-core/freebayes/main.nf'

include {
    REPORT
} from '../modules/local/report'

include {
    SAMTOOLS_FAIDX
} from '../modules/nf-core/samtools/faidx/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow NANOPORE {

    ch_versions = Channel.empty()

    // SUBWORKFLOW: prepare reference databases ...
    //
    PREPARE_REFERENCES ()
    ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.versions)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (ch_input)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    nanopore_reads = INPUT_CHECK.out.longreads
    //nanopore_reads.view()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        sequence quality control
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    QC_NANOPORE(
        nanopore_reads,
        PREPARE_REFERENCES.out.ch_hostile_ref_minimap2
    )
    QC_NANOPORE.out.qc_reads
        .filter {meta, reads -> reads.size() > 0 && reads.countFastq() > 0}
        .set { nanopore_reads }
    ch_versions = ch_versions.mix(QC_NANOPORE.out.versions)

     /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Identify the closely related rerference through mash
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    GET_REF_BY_MASH(
        nanopore_reads, 
        PREPARE_REFERENCES.out.ch_flu_db_msh, 
        PREPARE_REFERENCES.out.ch_flu_db_fasta
    )
    ch_versions.mix(GET_REF_BY_MASH.out.versions)
    ch_screen = GET_REF_BY_MASH.out.screen
    ch_fasta_fai = GET_REF_BY_MASH.out.fasta_fai

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Mapping to the closely related rerference
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    nanopore_reads.join(ch_fasta_fai).multiMap{
        it ->
            nanopore_reads: [it[0], it[1]]
            fasta: [it[0], it[2]]
    }.set{
        ch_input
    }

    MAPPING_NANOPORE(
        ch_input.nanopore_reads,
        ch_input.fasta
    )
    ch_versions.mix(MAPPING_NANOPORE.out.versions)
    

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        filter out secondary, supplementary, duplicates, and keep soft/hard clipped alignments
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    
    MAPPING_NANOPORE.out.bam_bai.join(ch_fasta_fai).multiMap{
        it ->
            bam_bai: [it[0], it[1], it[2]]
            fasta_fai: [it[0], it[3], it[4]]
    }.set{
        ch_input
    }
    PREPROCESS_BAM(ch_input.bam_bai, ch_input.fasta_fai)
    ch_versions.mix(PREPROCESS_BAM.out.versions)
    ch_screen.join(PREPROCESS_BAM.out.coverage).multiMap {
        it ->
            screen: [it[0], it[1]]
            coverage: [it[0], it[2]]
    }.set{ 
        ch_input 
    }
    MAPPING_SUMMARY(ch_input.screen, ch_input.coverage)
    ch_versions.mix(MAPPING_SUMMARY.out.versions)
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        variant calling
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    PREPROCESS_BAM.out.bam_bai.join(ch_fasta_fai).multiMap{
            it ->
                    bam_bai: [it[0], it[1], it[2]]
                    fasta_fai: [it[0], it[3], it[4]]
            }.set { ch_input }

    
    VARIANTS_NANOPORE( ch_input.bam_bai, ch_input.fasta_fai)

  

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    consensus generrating
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    BAM2LOW_DEPTH_BED(PREPROCESS_BAM.out.bam_bai)
    ch_versions.mix(BAM2LOW_DEPTH_BED.out.versions)
    
    VARIANTS_NANOPORE.out.vcf.join(ch_fasta_fai)
        .join(BAM2LOW_DEPTH_BED.out.bed)
        .multiMap{
            it ->
                vcf: [it[0], it[1]]
                fasta: [it[0], it[2]]
                mask: [it[0], it[4]]
        }.set{
            ch_input
        }

    //ch_input.vcf.view()
    CONSENSUS(ch_input.vcf, ch_input.fasta, ch_input.mask)
    ch_versions.mix(CONSENSUS.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Typing consensus
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
     /* ########## typing contigs ########## */
    CLASSIFIER_BLAST(CONSENSUS.out.fasta, PREPARE_REFERENCES.out.ch_typing_db)
    ch_versions.mix(CLASSIFIER_BLAST.out.versions)

    CONSENSUS.out.fasta.join(CLASSIFIER_BLAST.out.tsv).multiMap{
        it ->
            consensus: [it[0], it[1]]
            tsv: [it[0], it[2]]
    }.set{
        ch_input
    }

    CLASSIFIER_NEXTCLADE(ch_input.consensus, ch_input.tsv)
    ch_versions.mix(CLASSIFIER_NEXTCLADE.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    producce analysis report
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_report = QC_NANOPORE.out.input_stats.ifEmpty([])
        .join(QC_NANOPORE.out.qc_stats.ifEmpty([]))
        .join(CONSENSUS.out.stats.ifEmpty([]))
        .join(CLASSIFIER_BLAST.out.tsv.ifEmpty([]))
        .join(MAPPING_SUMMARY.out.mapping_summary.ifEmpty([]))
        .join(CLASSIFIER_NEXTCLADE.out.tsv.groupTuple())//.view()
        .multiMap{
            it ->
                raw_stats: [it[0], it[1]]
                qc_stats: [it[0], it[2]]
                consensus_stats: [it[0], it[3]]
                blastn_outfmt6: [it[0], it[4]]
                mapping_summary: [it[0], it[5]]
                nextclade_csv: it[6].join(',')
        }

    REPORT(
        ch_report.raw_stats,
        ch_report.qc_stats,
        ch_report.consensus_stats,
        ch_report.blastn_outfmt6,
        ch_report.mapping_summary,
        ch_report.nextclade_csv
    )

    //report_name = "summary_report." + params.nanopore_reads_mapping_tool + "-" + params.nanopore_variant_caller
    CONCAT_STATS_SUMMARY(
        REPORT.out.oneline_csv.map {
            cfg, stats -> stats }.collect()
            //.map { files -> tuple([id: "summary_report"], files)
            .map { files -> tuple([id: "report.consensus_" + params.nanopore_reads_mapping_tool + "-" + params.nanopore_variant_caller], files)

        },
        "tsv",
        "tsv"
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
