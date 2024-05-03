include {
    BLAST_BLASTN 
} from '../../modules/local/blast/blastn'

include { 
    CSVTK_ADD_HEADER as CSVTK_ADD_HEADER_BLASTN 
} from '../../modules/local/csvtk/add-header/main'

include {
    CSVTK_CONCAT as CONCAT_TYPING;
    
} from '../../modules/nf-core/csvtk/concat/main'

workflow CLASSIFIER_BLAST {   

    take:
        fasta //input contig
        typing_db // reference db
        
    main:
        ch_versions = Channel.empty()
        in_format = "tsv"
        out_format = "tsv"
        
        BLAST_BLASTN(fasta, typing_db)
        ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

        CSVTK_ADD_HEADER_BLASTN(
            BLAST_BLASTN.out.tsv, 
            "qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen,qcovs"
        )
        CONCAT_TYPING(CSVTK_ADD_HEADER_BLASTN.out.tsv.map { cfg, tsv -> tsv }.collect().map { files -> tuple([id: "typing"], files)}, in_format, out_format )


    emit:
        tsv = CSVTK_ADD_HEADER_BLASTN.out.tsv
        summary = CONCAT_TYPING.out.csv
        versions = ch_versions
        
}


include {
    NEXTCLADE_RUN 
} from '../../modules/local/nextclade/run/main'

include {
    SEGMENT2DATASET
} from '../../modules/local/misc'

include {
    CSVTK_CONCAT as CONCAT_NEXTCLADE;    
} from '../../modules/nf-core/csvtk/concat/main'

workflow CLASSIFIER_NEXTCLADE{   

    take:
        fasta
        tsv //blast output
    main:
       
        ch_versions = Channel.empty()
        in_format = "tsv"
        out_format = "tsv"
        fasta.splitFasta(by: 1, file: true)
            .set{
                ch_fasta
            }
        //ch_fasta.view()
        SEGMENT2DATASET(fasta, tsv)
        SEGMENT2DATASET.out.out_tsv.filter{meta, tsv -> tsv.size() > 0}
            .splitCsv ( header:false, sep:'\t' )
            .multiMap{
                meta, row ->
                    def new_meta = [:]
                    new_meta.id = meta.id
                    new_meta.single_end = meta.single_end
                    new_meta.basecaller_mode = meta.basecaller_mode
                    new_meta.seqid = row[0]
                fasta: [new_meta, file(row[1], checkIfExists: true)]
                dataset: [new_meta, file(row[2], checkIfExists: true)]
            }.set{
                ch_input
            }  
      
        NEXTCLADE_RUN(ch_input.fasta, ch_input.dataset)
        
        ch_versions = ch_versions.mix(NEXTCLADE_RUN.out.versions.first())
        CONCAT_NEXTCLADE(NEXTCLADE_RUN.out.tsv.map { cfg, tsv -> tsv }.collect().map { files -> tuple([id: "nextclade"], files)}, in_format, out_format )

 
    emit:
        tsv = NEXTCLADE_RUN.out.tsv
        dbname = NEXTCLADE_RUN.out.dbname
        summary = CONCAT_NEXTCLADE.out.csv
        versions = ch_versions
        
}


