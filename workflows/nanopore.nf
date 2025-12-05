/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    nanopore_variant_callers   : ['clair3'],
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
    //params.hostile_human_ref_minimap2,
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
    SEEK_REFERENCES
} from '../subworkflows/local/seek_references'

include {
    MAPPING_NANOPORE
} from '../subworkflows/local/mapping_nanopore'

include {
    PREPROCESS_BAM
} from '../subworkflows/local/preprocess_bam'

include {
    CLASSIFIER_BLAST;
    CLASSIFIER_NEXTCLADE;
} from '../subworkflows/local/classifier'


include {
    VARIANTS_NANOPORE;
} from '../subworkflows/local/variants_nanopore'
include {
    PREPROCESS_VCF;

} from '../subworkflows/local/preprocess_vcf'
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
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//~~~~~~~~~~
    */
    if(!params.skip_nanopore_reads_qc){
        QC_NANOPORE(
            nanopore_reads,
            //PREPARE_REFERENCES.out.ch_hostile_ref_minimap2
        )
        ch_versions = ch_versions.mix(QC_NANOPORE.out.versions)
        QC_NANOPORE.out.qc_reads
            .filter {meta, reads -> reads.size() > 0 && reads.countFastq() > 0}
            .set { nanopore_reads }

    }

     /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Identify the closely related rerference through mash
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_screen = Channel.empty()
    ch_fasta_fai = Channel.empty()
    ch_seq_header = Channel.empty()

    SEEK_REFERENCES(
        nanopore_reads,
        PREPARE_REFERENCES.out.ch_flu_db_msh,
        PREPARE_REFERENCES.out.ch_flu_db_fasta
    )
    ch_versions = ch_versions.mix(SEEK_REFERENCES.out.versions)
    SEEK_REFERENCES.out.screen
        .filter{ it[1] != null}
        .filter{meta, tsv -> tsv.size() > 0 && tsv.countLines() > 0}
        .set{ch_screen}

    SEEK_REFERENCES.out.fasta_fai
        .filter{ meta, fasta, fai -> fasta.size() > 0 && fasta.countFasta() > 0}
        .set{ch_fasta_fai}

    SEEK_REFERENCES.out.seq_header
        .filter{meta, seq_header -> seq_header.size() > 0 && seq_header.countLines() > 0}
        .set{ch_seq_header}


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         Map reads to the identified references
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
    ch_versions = ch_versions.mix(MAPPING_NANOPORE.out.versions)


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        filter out secondary, supplementary, duplicates, and keep soft/hard clipped alignments
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    MAPPING_NANOPORE.out.bam_bai
        .join(ch_fasta_fai).join(ch_seq_header)
        .multiMap{
            it ->
                bam_bai: [it[0], it[1], it[2]]
                fasta: [it[0], it[3]]
                fai: [it[0], it[4]]
                seq_header: [it[0], it[5]]
        }
        .set{ch_input}

    ch_coverage = Channel.empty()

    PREPROCESS_BAM(ch_input.bam_bai, ch_input.fasta, ch_input.fai, ch_input.seq_header)
    ch_versions = ch_versions.mix(PREPROCESS_BAM.out.versions)
    PREPROCESS_BAM.out.coverage
        .filter{meta, tsv -> tsv.size() > 0 && tsv.countLines() > 0}
        .set{ch_coverage}


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
    VARIANTS_NANOPORE(ch_input.bam_bai, ch_input.fasta_fai)
    ch_versions = ch_versions.mix(VARIANTS_NANOPORE.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       preprocess vcf before consensus: norm, filter
       todo: snpeff has problems ........................
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    //produce bed file for the low depth region
    ch_input = PREPROCESS_BAM.out.bam_bai
        .map{
            it -> [it[0], it[1], 1] //[meta, bam, scale]
        }

        //bed file of the low depth regions < params.mindepth
    BEDTOOLS_GENOMECOV_LOWDEPTH(ch_input, [], "bed")

    VARIANTS_NANOPORE.out.vcf_tbi.join(ch_fasta_fai)
        .join(BEDTOOLS_GENOMECOV_LOWDEPTH.out.genomecov)
        .multiMap{
            it ->
                vcf_tbi: [it[0], it[1], it[2]]
                fasta: [it[0], it[3]]
        }.set{
            ch_input
        }

    PREPROCESS_VCF(
        ch_input.vcf_tbi,
        ch_input.fasta,
        params.snpeff_db,
        PREPARE_REFERENCES.out.ch_snpeff_config,
        PREPARE_REFERENCES.out.ch_snpeff_dataDir
    )


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    consensus generrating
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //produce bed file for the low depth region
   /*
    ch_input = PREPROCESS_BAM.out.bam_bai
        .map{
            it -> [it[0], it[1], 1] //[meta, bam, scale]
        }

        //bed file of the low depth regions < params.mindepth
    BEDTOOLS_GENOMECOV_LOWDEPTH(ch_input, [], "bed")


    VARIANTS_NANOPORE.out.vcf_tbi.join(ch_fasta_fai)
        .join(BEDTOOLS_GENOMECOV_LOWDEPTH.out.genomecov)
        .multiMap{
            it ->
                vcf_tbi: [it[0], it[1], it[2]]
                fasta: [it[0], it[3]]
                mask: [it[0], it[5]]
        }.set{
            ch_input
        }
*/
 PREPROCESS_VCF.out.vcf_tbi
        .join(ch_fasta_fai)
        .join(BEDTOOLS_GENOMECOV_LOWDEPTH.out.genomecov)
        .multiMap{
            it ->
                vcf_tbi_fasta: [it[0], it[1], it[2], it[3]]
                mask: [it[0], it[5]]
        }.set{
            ch_input
        }
    //ch_input.vcf.view()
    consensus_stats = Channel.empty()
    consensus_fasta = Channel.empty()

    CONSENSUS(ch_input.vcf_tbi_fasta, ch_input.mask)
    //CONSENSUS(ch_input.vcf_tbi, ch_input.fasta, ch_input.mask)
    ch_versions = ch_versions.mix(CONSENSUS.out.versions)
    //consensus_stats = CONSENSUS.out.stats.filter{ it != null}
    CONSENSUS.out.stats
        //.filter{ it != null}
        .filter{meta, tsv -> tsv.size() > 0 && tsv.countLines() > 0}
        .set{consensus_stats}
    //consensus_stats.view()

    CONSENSUS.out.fasta
        .filter {meta, fasta -> fasta.size() > 0 && fasta.countFasta() > 0}
        .set { consensus_fasta }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Typing consensus
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_typing = Channel.empty()
    CLASSIFIER_BLAST(consensus_fasta, PREPARE_REFERENCES.out.ch_typing_db)
    ch_versions = ch_versions.mix(CLASSIFIER_BLAST.out.versions)
    //ch_typing = CLASSIFIER_BLAST.out.tsv.filter{ it != null}
    CLASSIFIER_BLAST.out.tsv
        //.filter{ it != null}
        .filter{meta, tsv -> tsv.size() > 0 && tsv.countLines() > 1}
        .set{ch_typing}

    consensus_fasta.join(ch_typing).multiMap{
        it ->
            consensus: [it[0], it[1]]
            typing: [it[0], it[2]]
    }.set{
        ch_input
    }
    //ch_input.consensus.view()
    CLASSIFIER_NEXTCLADE(ch_input.consensus, ch_input.typing)
    ch_versions = ch_versions.mix(CLASSIFIER_NEXTCLADE.out.versions)

     /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    producce consensus report
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */


    ch_nextclade_dbs = Channel.empty()
    ch_nextclade_tsv = Channel.empty()

    ch_nextclade_dbs = CLASSIFIER_NEXTCLADE.out.dbname
        //.filter{ it != null}
        .filter{meta, tsv -> tsv.size() > 0 && tsv.countLines() > 1}
        .map{
            meta, tsv ->
                meta.remove("seqid")
                [meta, tsv]
        }.groupTuple().view()
        //.groupTuple()

    ch_nextclade_tsv = CLASSIFIER_NEXTCLADE.out.tsv
        //.filter{ it != null}
        .filter{meta, tsv -> tsv.size() > 0 && tsv.countLines() > 1}
        .map{
            meta, tsv ->
                meta.remove("seqid")
                [meta, tsv]
        }.groupTuple()


    /* ch_merge = ch_screen.join(ch_coverage)//.view()
        .join(consensus_stats)//.view()
        .join(ch_typing)//.view()
        .join(ch_nextclade_tsv)//.view()
        .join(ch_nextclade_dbs)//.view()
        .view() */

    ch_merge = ch_screen.join(ch_coverage, remainder: true)//.view()
        .join(consensus_stats, remainder: true)//.view()
        .join(ch_typing, remainder: true)//.view()
        .join(ch_nextclade_tsv, remainder: true)//.view()
        .join(ch_nextclade_dbs, remainder: true)//.view()
        .view()

    // ch_screen
    //     .join(ch_coverage, remainder: true)//.view()
    //     .join(consensus_stats, remainder: true)//.view()
    //     .join(ch_typing, remainder: true)//.view()
    //     .join(ch_nextclade_tsv, remainder: true)//.view()
    //     .join(ch_nextclade_dbs, remainder: true)//.view()
    ch_merge.multiMap{
        it ->
            screen:  it[1] != null ? [it[0], it[1]] : [[], []]
            cov:  it[2] != null ? [it[0], it[2]] : [[], []]
            stats:  it[3] != null ? [it[0], it[3]] : [[], []]
            typing:  it[4] != null ? [it[0], it[4]] : [[], []]
            nextclade_tsv: it[5] != null ? [it[0], it[5].join(',')] : [[], []]
            nextclade_dbname: it[6] != null ? [it[0], it[6].join(',')] : [[], []]
            dbver: [it[0], params.flu_db_ver]
            pipelinever: [it[0], "${workflow.manifest.version}"]
        }.set{
            ch_input
        }

    CONSENSUS_REPORT(
        ch_input.stats,
        ch_input.cov,
        ch_input.typing,
        ch_input.nextclade_tsv,
        ch_input.nextclade_dbname,
        ch_input.screen,
        ch_input.dbver,
        ch_input.pipelinever

    )

    CONCAT_CONSENSU_REPORT(
        CONSENSUS_REPORT.out.csv.map { cfg, stats -> stats }.collect()//.view()
            .map{
                files ->
                    tuple(
                        [id: "${params.mapping_tool}-${params.variant_caller}"],
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
