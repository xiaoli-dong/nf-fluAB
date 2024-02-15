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
    QC_ILLUMINA
} from '../subworkflows/local/qc_illumina'

include {
    MAPPING_ILLUMINA
} from '../subworkflows/local/mapping_illumina'

include {
    consensus
} from '../subworkflows/local/consensus'

include {
    classifier_blast;
    classifier_nextclade;
} from '../subworkflows/local/classifier'

include {
    maskfasta
} from '../subworkflows/local/maskfasta'

include {
    variants_bcftools;
    variants_freebayes;

} from '../subworkflows/local/variants_illumina'

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

include {FREEBAYES
} from '../modules/nf-core/freebayes/main.nf'

include {
    BCFTOOLS_MPILEUP
} from '../modules/nf-core/bcftools/mpileup/main.nf'

include {
    REPORT
} from '../modules/local/report'

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
    PREPARE_REFERENCES ()
    ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.versions)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (ch_input)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    illumina_reads = INPUT_CHECK.out.shortreads


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        sequence quality control
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */


    if(!params.skip_illumina_reads_qc){
        QC_ILLUMINA(
            illumina_reads,
            PREPARE_REFERENCES.out.ch_flu_primers,
            PREPARE_REFERENCES.out.ch_hostile_ref_bowtie2
        )
        QC_ILLUMINA.out.qc_stats.view()
        ch_versions = ch_versions.mix(QC_ILLUMINA.out.versions)

        //get rid of zero size read file and avoid the downstream crash
        QC_ILLUMINA.out.qc_reads
            .filter {meta, reads -> reads[0].size() > 0 && reads[0].countFastq() > 0}
            .set { illumina_reads }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Identify the closely related rerference through mapping
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    MAPPING_ILLUMINA(
        illumina_reads,
        PREPARE_REFERENCES.out.ch_flu_db_msh,
        PREPARE_REFERENCES.out.ch_flu_db_fasta
    )
    ch_versions.mix(MAPPING_ILLUMINA.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        variant calling
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // mask low depth region of the fasta reference
    maskfasta(
        MAPPING_ILLUMINA.out.bam_bai,
        MAPPING_ILLUMINA.out.fasta
    )
    fasta_fai = maskfasta.out.fasta_fai
    ch_versions.mix(maskfasta.out.versions)
    if(params.illumina_variant_caller == 'freebayes'){
        MAPPING_ILLUMINA.out.bam_bai.join(fasta_fai).multiMap {
            it ->
                bam_bai: [it[0], it[1], it[2]]
                fasta_fai: [it[0], it[3], it[4]]
        }
        .set{
            ch_input
        }
        variants_freebayes(ch_input.bam_bai, ch_input.fasta_fai)
        vcf_tbi = variants_freebayes.out.vcf_tbi
        fasta_fai = variants_freebayes.out.fasta_fai
        ch_versions.mix(variants_freebayes.out.versions)
    }
    else if(params.illumina_variant_caller == 'bcftools'){
        MAPPING_ILLUMINA.out.bam_bai.join(fasta_fai).multiMap {
            it ->
                bam_bai: [it[0], it[1], it[2]]
                fasta_fai: [it[0], it[3], it[4]]
        }
        .set{
            ch_input
        }
        variants_bcftools(ch_input.bam_bai, ch_input.fasta_fai)
        vcf_tbi = variants_bcftools.out.vcf_tbi
        fasta_fai = variants_bcftools.out.fasta_fai
        ch_versions.mix(variants_bcftools.out.versions)
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    consensus generrating
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    vcf_tbi.join(fasta_fai).multiMap{
        it ->
            //vcf_tbi: [it[0], it[1], it[2]]
            vcf: [it[0], it[1]]
            fasta: [it[0], it[3]]
    }.set{
        ch_input
    }
    consensus(ch_input.vcf, ch_input.fasta)
    ch_versions.mix(consensus.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Typing consensus
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

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
    ch_versions.mix(classifier_nextclade.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    producce analysis report
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_report = QC_ILLUMINA.out.input_stats.ifEmpty([])
        .join(QC_ILLUMINA.out.qc_stats.ifEmpty([]))
        .join(consensus.out.stats.ifEmpty([]))
        .join(classifier_blast.out.tsv.ifEmpty([]))
        .join(MAPPING_ILLUMINA.out.mapping_summary.ifEmpty([]))
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
        REPORT.out.oneline_csv.map { cfg, stats -> stats }.collect()
            .map{
                files ->
                    tuple(
                        [id: "report.consensus_" + params.illumina_reads_mapping_tool + "-" + params.illumina_variant_caller],
                        files
                    )
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
