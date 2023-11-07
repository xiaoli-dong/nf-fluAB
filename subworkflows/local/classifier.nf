include {BLAST_BLASTN } from '../../modules/local/blast/blastn'

workflow classifier_blast {   

    take:
        fasta //input contig
        typing_db // reference db
        
    main:
        ch_versions = Channel.empty()
        BLAST_BLASTN(fasta, typing_db)
        ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())
        

    emit:
        // txt = blast_out_ch
        // seq = consensus_seq_ch
        //fasta
        tsv = BLAST_BLASTN.out.tsv
        versions = ch_versions
        
}


include {NEXTCLADE_RUN } from '../../modules/local/nextclade/run/main'
include {SEGMENT2TYPEDATA } from '../../modules/local/segment2typedata'
include {
    CSVTK_CONCAT
    
} from '../../modules/nf-core/csvtk/concat/main'

workflow classifier_nextclade {   

    take:
        fasta
        tsv
    main:
       
        ch_versions = Channel.empty()
        fasta.splitFasta(by: 1, file: true)
            .set{
                ch_fasta
            }
        //ch_fasta.view()
        SEGMENT2TYPEDATA(fasta, tsv)
        SEGMENT2TYPEDATA.out.out_tsv//.view()

        SEGMENT2TYPEDATA.out.out_tsv
            .splitCsv(header: true, sep:'\t')
            .multiMap{
                it ->
                fasta: [ it[0], it[1].fasta_path]
                dataset: [it[1].typedata]
            }
            .set{nextclade_input}

        //nextclade_input.fasta.view()
        NEXTCLADE_RUN(nextclade_input.fasta, nextclade_input.dataset)
        ch_versions = ch_versions.mix(NEXTCLADE_RUN.out.versions.first())
 
    emit:
        tsv = NEXTCLADE_RUN.out.tsv
        //tsv = CSVTK_CONCAT.out.csv
        versions = ch_versions
        
}

