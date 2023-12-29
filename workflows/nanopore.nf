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
    MAPPING_NANOPORE
} from '../subworkflows/local/mapping_nanopore'

include {
    classifier_blast;
    classifier_nextclade;
} from '../subworkflows/local/classifier'

include {
    maskfasta
} from '../subworkflows/local/maskfasta'

include {
    variants_clair3;
} from '../subworkflows/local/variants_nanopore' 

include {  
    consensus
} from '../subworkflows/local/consensus'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include {   CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include {
    CSVTK_CONCAT as CONCAT_STATS_SHORT_RAW;
    CSVTK_CONCAT as CONCAT_STATS_SHORT_QC;
    CSVTK_CONCAT as CONCAT_STATS_SUMMARY;
    CSVTK_CONCAT as CONCAT_NEXTCLADE;
} from '../modules/nf-core/csvtk/concat/main'
include {   FREEBAYES   } from '../modules/nf-core/freebayes/main.nf'


include {   REPORT      } from '../modules/local/report'
include {   CLAIR3      } from '../modules/local/clair3/main'

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

    nanopore_reads = INPUT_CHECK.out.reads
    //nanopore_reads.view()
    
    QC_NANOPORE(nanopore_reads, PREPARE_REFERENCES.out.ch_hostile_ref_minimap2)
    QC_NANOPORE.out.qc_reads
        .filter {meta, reads -> reads.size() > 0 && reads.countFastq() > 0}
        .set { nanopore_reads }
    ch_versions = ch_versions.mix(QC_NANOPORE.out.versions) 
    
    MAPPING_NANOPORE(
        nanopore_reads, 
        PREPARE_REFERENCES.out.ch_flu_db_msh, 
        PREPARE_REFERENCES.out.ch_flu_db_fasta
    )
    ch_versions.mix(MAPPING_NANOPORE.out.versions)

     // mask low depth region of the fasta reference 
    maskfasta(
        MAPPING_NANOPORE.out.bam_bai, 
        MAPPING_NANOPORE.out.fasta
    )
    fasta_fai = maskfasta.out.fasta_fai
    ch_versions.mix(maskfasta.out.versions)


    MAPPING_NANOPORE.out.bam_bai.join(maskfasta.out.fasta_fai).multiMap {
        it ->
            bam_bai: [it[0], it[1], it[2]]
            fasta_fai: [it[0], it[3], it[4]]
    }
    .set{
        ch_input
    }
    variants_clair3(
        ch_input.bam_bai, 
        ch_input.fasta_fai, 
        PREPARE_REFERENCES.out.ch_clair3_variant_model
    )
   // vcf_tbi = variants_clair3.out.vcf_tbi

    variants_clair3.out.vcf_tbi.join(variants_clair3.out.fasta_fai).multiMap{
        it ->
            vcf_tbi: [it[0], it[1], it[2]]
            fasta: [it[0], it[3]]
    }.set{
        ch_input
    } 
    consensus(ch_input.vcf_tbi, ch_input.fasta)
    
     /* ########## typing contigs ########## */
    classifier_blast(consensus.out.fasta, PREPARE_REFERENCES.out.ch_typing_db)
    ch_versions.mix(classifier_blast.out.versions)

    consensus.out.fasta.join(classifier_blast.out.tsv).multiMap{
        it ->
            consensus: [it[0], it[1]]
            tsv: [it[0], it[2]]
    }.set{
        ch_input
    }

    classifier_nextclade(ch_input.consensus, ch_input.tsv)
    //classifier_nextclade.out.tsv.collectFile().view()
    ch_versions.mix(classifier_nextclade.out.versions)
   
    //reporting
    ch_report = QC_NANOPORE.out.input_stats.ifEmpty([])
        .join(QC_NANOPORE.out.qc_stats.ifEmpty([]))
        .join(consensus.out.stats.ifEmpty([]))
        .join(classifier_blast.out.tsv.ifEmpty([]))
        .join(MAPPING_NANOPORE.out.mapping_summary.ifEmpty([]))
        .join(classifier_nextclade.out.tsv.groupTuple())//.view()
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

    CONCAT_STATS_SUMMARY(
        REPORT.out.oneline_csv.map { 
            cfg, stats -> stats }.collect()
            .map { files -> tuple([id: "summary_report"], files)
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
