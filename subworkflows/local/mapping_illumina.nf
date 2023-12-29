include {   BWAMEM2_INDEX       } from '../../modules/nf-core/bwamem2/index/main'
include {   BWAMEM2_MEM         } from '../../modules/nf-core/bwamem2/mem/main'
include {   MINIMAP2_ALIGN      } from '../../modules/local/minimap2/align/main'
include {   SAMTOOLS_FAIDX      } from '../../modules/nf-core/samtools/faidx/main'
include {   SAMTOOLS_INDEX      } from '../../modules/nf-core/samtools/index/main'
include {   SAMTOOLS_COVERAGE   } from '../../modules/nf-core/samtools/coverage/main'
include {   SAMTOOLS_SORT       } from '../../modules/nf-core/samtools/sort/main'
include {   MASH_SCREEN         } from '../../modules/nf-core/mash/screen/main'
include {   SEQKIT_GREP         } from '../../modules/local/seqkit/grep/main'
include {   CSVTK_ADD_HEADER    } from '../../modules/local/csvtk/add-header/main'
include {   CSVTK_JOIN    } from '../../modules/nf-core/csvtk/join/main'

include {
    FILTERMASH;
    MAPPING_SUMMARY;
    
} from '../../modules/local/misc'

workflow MAPPING_ILLUMINA {   

    take:
        illumina_reads
        ch_flu_db_msh
        ch_flu_db_fasta

    main:
        ch_versions = Channel.empty()
        
        // ################################# fetch reference #############################

        MASH_SCREEN(illumina_reads, ch_flu_db_msh)
        ch_versions = ch_versions.mix(MASH_SCREEN.out.versions.first())

        
        //filter out empty channel
        FILTERMASH(
            MASH_SCREEN.out.screen.filter{meta, screen -> screen.countLines() > 0}
        )
        ch_versions = ch_versions.mix(FILTERMASH.out.versions.first())
         
        CSVTK_ADD_HEADER(
            FILTERMASH.out.screen, 
            "identity,shared-hashes,median-multiplicity,p-value,query-ID,query-comment"
        )
        CSVTK_ADD_HEADER.out.tsv.map{
            it ->
            def meta = it[0]
            def arr = it[1].splitCsv(skip: 1, sep: '\t')
            a = []
            arr.each{
                n ->
                a.add(n[4])
            }
            [meta, a.join(',')]
        }.set {accession_list}
        accession_list.view()
        SEQKIT_GREP(
            accession_list,
            ch_flu_db_fasta

        )
        ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions.first())
        screen = FILTERMASH.out.screen
        //reference sequences
        fasta = SEQKIT_GREP.out.fasta
        //###################################### endo fo fetch reference ##############################

        // SAMTOOLS_FAIDX(fasta)
        // fasta_fai = SAMTOOLS_FAIDX.out.fai

        /* ######################################## Mapping ############################ */
        bam = Channel.empty()
        if ( params.illumina_reads_mapping_tool == 'minimap2' ){
            bam_format = true
            cigar_paf_format = false
            cigar_bam = false
           
            illumina_reads.join(fasta).multiMap{
                it ->
                    reads: [it[0], it[1]]
                    fasta: it[2]
            }
            .set{
                ch_input
            }

            MINIMAP2_ALIGN(ch_input.reads, ch_input.fasta, bam_format, cigar_paf_format, cigar_bam)
            ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())
            bam = MINIMAP2_ALIGN.out.bam
        }
        else if(params.illumina_reads_mapping_tool == 'bwa'){
           
            BWAMEM2_INDEX(fasta)
            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())

            illumina_reads.join(BWAMEM2_INDEX.out.index).multiMap{
                it ->
                reads: [it[0], it[1]]
                index: [it[0], it[2]]
            }.set{
                ch_input
            }

            //only keep primary alignment
            //bwa mem | samtools view -h -F260
            BWAMEM2_MEM (
                ch_input.reads,
                ch_input.index,
                false
            )
            ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
            bam = BWAMEM2_MEM.out.bam
        }
        
        // ################################### process bam file ##########################
        SAMTOOLS_SORT (
            //BWAMEM2_MEM.out.bam
            bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

        //SAMTOOLS_SORT.out.bam.view()

        SAMTOOLS_INDEX (
            SAMTOOLS_SORT.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())
        //SAMTOOLS_INDEX.out.bai.view()
        
        //SAMTOOLS_DEPTH(SAMTOOLS_SORT.out.bam)
        //ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

        ch_bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)
        //input_ch.view()
        SAMTOOLS_COVERAGE(ch_bam_bai)
        coverage = SAMTOOLS_COVERAGE.out.coverage
        // ###################################### end of process bam #############

        screen.join(coverage)
            .multiMap {
                it ->
                    screen: [it[0], it[1]]
                    coverage: [it[0], it[2]]
            }.set{
                ch_input
            }
        
        ch_input_merge = coverage.join(CSVTK_ADD_HEADER.out.tsv).map{
            it -> [it[0], [it[1], it[2]]]
        }
        //CSVTK_JOIN(ch_input_merge)

        MAPPING_SUMMARY(ch_input.screen, ch_input.coverage)
        //GET_COVERAGE.out.csv.view()


        
    emit:
        bam_bai = ch_bam_bai
        fasta
        //fasta_fai
        coverage
        mapping_summary = MAPPING_SUMMARY.out.mapping_summary
        versions = ch_versions
        
}