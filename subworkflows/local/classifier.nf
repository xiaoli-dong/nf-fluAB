include {BLAST_BLASTN } from '../../modules/local/blast/blastn'
include { CSVTK_ADD_HEADER as CSVTK_ADD_HEADER_BLASTN } from '../../modules/local/csvtk/add-header/main'

workflow CLASSIFIER_BLAST {   

    take:
        fasta //input contig
        typing_db // reference db
        
    main:
        ch_versions = Channel.empty()
        BLAST_BLASTN(fasta, typing_db)
        ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

        CSVTK_ADD_HEADER_BLASTN(
            BLAST_BLASTN.out.tsv, 
            "qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen,qcovs"
        )
        

    emit:
        tsv = CSVTK_ADD_HEADER_BLASTN.out.tsv
        versions = ch_versions
        
}


include {NEXTCLADE_RUN } from '../../modules/local/nextclade/run/main'
include {SEGMENT2TYPEDATA } from '../../modules/local/misc'
include {
    CSVTK_CONCAT
    
} from '../../modules/nf-core/csvtk/concat/main'

workflow CLASSIFIER_NEXTCLADE{   

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

