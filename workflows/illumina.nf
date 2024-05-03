/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    illumina_reads_qc_tools        : ['bbduk', 'fastp'],
    illumina_variant_callers   : ['freebayes', 'bcftools'],
    illumina_reads_mapping_tools : ['bwa', 'minimap2']
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowIllumina.initialise(params, log, valid_params)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ 
    params.input, 
    params.hostile_human_ref_minimap2, 
    params.hostile_human_ref_bowtie2,
    params.flu_primers, 
    params.typing_db, 
    params.flu_db_msh, 
    params.flu_db_fasta,
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
    QC_ILLUMINA
} from '../subworkflows/local/qc_illumina'

include {
    GETREF_BY_MASH
} from '../subworkflows/local/get_references'

include {
    MAPPING_ILLUMINA
} from '../subworkflows/local/mapping_illumina'

include {
    PREPROCESS_BAM
} from '../subworkflows/local/preprocess_bam'

include {
    CONSENSUS
} from '../subworkflows/local/consensus'

include {
    CLASSIFIER_BLAST;
    CLASSIFIER_NEXTCLADE;
} from '../subworkflows/local/classifier'


include {
    VARIANTS_ILLUMINA;

} from '../subworkflows/local/variants_illumina'

include { 
    RUN_POLYPOLISH; 
    RUN_POLCA;    
} from '../subworkflows/local/polisher_illumina'

//
// MODULE: Installed directly from nf-core/modules
//
include {
    CUSTOM_DUMPSOFTWAREVERSIONS
} from '../modules/nf-core/custom/dumpsoftwareversions/main'

include {
    CSVTK_CONCAT as CONCAT_CONSENSU_REPORT;
} from '../modules/nf-core/csvtk/concat/main'

include {
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_LOWDEPTH;
} from '../modules/nf-core/bedtools/genomecov/main.nf'

//
// MODULE: developed locally 
//
include {
    CONSENSUS_REPORT
} from '../modules/local/consensus/report/main'


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
    //illumina_reads.view()

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
        //QC_ILLUMINA.out.qc_stats.view()
        ch_versions = ch_versions.mix(QC_ILLUMINA.out.versions)

        //get rid of zero size read file and avoid the downstream crash
        QC_ILLUMINA.out.qc_reads
            .filter {meta, reads -> reads[0].size() > 0 && reads[0].countFastq() > 0}
            .set { illumina_reads }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Identify the closely related rerference through mash
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    GETREF_BY_MASH(
        illumina_reads, 
        PREPARE_REFERENCES.out.ch_flu_db_msh, 
        PREPARE_REFERENCES.out.ch_flu_db_fasta
    )
    ch_versions.mix(GETREF_BY_MASH.out.versions)
    ch_screen = GETREF_BY_MASH.out.screen
    ch_fasta_fai = GETREF_BY_MASH.out.fasta_fai
    ch_header = GETREF_BY_MASH.out.header
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Map reads to the identified references
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    illumina_reads.join(ch_fasta_fai).multiMap{
        it ->
            illumina_reads: [it[0], it[1]]
            fasta: [it[0], it[2]]
    }.set{
        ch_input
    }

    MAPPING_ILLUMINA(
        ch_input.illumina_reads,
        ch_input.fasta
    )
    ch_versions.mix(MAPPING_ILLUMINA.out.versions)
    

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        filter out secondary, supplementary, duplicates, and soft/hard clipped reads
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    MAPPING_ILLUMINA.out.bam_bai.join(ch_fasta_fai).join(ch_header).multiMap{
        it ->
            bam_bai: [it[0], it[1], it[2]]
            fasta_fai: [it[0], it[3], it[4]]
            header: [it[0], it[5]]
    }.set{
        ch_input
    }
    PREPROCESS_BAM(ch_input.bam_bai, ch_input.fasta_fai, ch_input.header)
    ch_versions.mix(PREPROCESS_BAM.out.versions)
   
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
    VARIANTS_ILLUMINA(ch_input.bam_bai, ch_input.fasta_fai)
    ch_versions.mix(VARIANTS_ILLUMINA.out.versions)

    

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    consensus generrating
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //produce bed file for the low depth region
    ch_input = PREPROCESS_BAM.out.bam_bai
        .map{
            it -> [it[0], it[1], 1] //[meta, bam, scale]
        }

        //bed file of the low depth regions < params.mindepth
    BEDTOOLS_GENOMECOV_LOWDEPTH(ch_input, [], "bed")
    
    

    VARIANTS_ILLUMINA.out.vcf_tbi.join(ch_fasta_fai)
        .join(BEDTOOLS_GENOMECOV_LOWDEPTH.out.genomecov)
        .multiMap{
            it ->
                vcf_tbi: [it[0], it[1], it[2]]
                //vcf: [it[0], it[1]]
                fasta: [it[0], it[3]]
                mask: [it[0], it[5]]
        }.set{
            ch_input
        }
    CONSENSUS(ch_input.vcf_tbi, ch_input.fasta, ch_input.mask)
    ch_versions.mix(CONSENSUS.out.versions)
    //CONSENSUS.out.stats.view()
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Typing consensus
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    CLASSIFIER_BLAST(CONSENSUS.out.fasta, PREPARE_REFERENCES.out.ch_typing_db)
    ch_versions.mix(CLASSIFIER_BLAST.out.versions)

    CONSENSUS.out.fasta.join(CLASSIFIER_BLAST.out.tsv).multiMap{
        it ->
            consensus: [it[0], it[1]]
            tsv: [it[0], it[2]]
    }.set{
        ch_input
    }
    //ch_input.consensus.view()
    CLASSIFIER_NEXTCLADE(ch_input.consensus, ch_input.tsv)
    ch_versions.mix(CLASSIFIER_NEXTCLADE.out.versions)
    
     /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    producce consensus report
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_nextclade_dbname = CLASSIFIER_NEXTCLADE.out.dbname.groupTuple().view()
    //CLASSIFIER_NEXTCLADE.out.tsv.view()
    ch_nextclade_tsv = CLASSIFIER_NEXTCLADE.out.tsv.map{
        meta, tsv ->
            meta.remove("seqid")
            [meta, tsv]
    }.groupTuple()//.view() 
    //CLASSIFIER_NEXTCLADE.out.dbname.view()
    CONSENSUS.out.stats
        .join(PREPROCESS_BAM.out.coverage_rmdup)
        .join(CLASSIFIER_BLAST.out.tsv)
        .join(ch_nextclade_tsv)
        .join(ch_nextclade_dbname)
        .join(ch_screen)
        .multiMap{
            it -> 
                stats: [it[0], it[1]]
                cov: [it[0], it[2]]
                typing: [it[0], it[3]]
                nextclade_tsv: [it[0], it[4].join(',')]
                nextclade_dbname: [it[0], it[5].join(',')]
                screen: [it[0], it[6]]
        }.set{
            ch_input
        }
    ch_input.nextclade_dbname.view()
    CONSENSUS_REPORT(ch_input.stats, ch_input.cov, ch_input.typing, ch_input.nextclade_tsv, ch_input.nextclade_dbname, ch_input.screen, Channel.of(params.flu_db_ver) )
    CONCAT_CONSENSU_REPORT(
        CONSENSUS_REPORT.out.csv.map { cfg, stats -> stats }.collect()
            .map{
                files ->
                    tuple(
                        [id: "consensus_report"],
                        files
                    )
        },
        "csv",
        "csv"
    )  
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
