include {BLAST_BLASTN } from '../../modules/local/blast/blastn'
include {BLAST_MAKEBLASTDB } from '../../modules/nf-core/blast/makeblastdb/main'

workflow RUN_TYPING {   

    take:
        fasta
        typingdb
        
    main:
        ch_versions = Channel.empty()

        BLAST_MAKEBLASTDB (typingdb)
        ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions.first())
        BLAST_BLASTN(fasta, BLAST_MAKEBLASTDB.out.db)
        ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())
        

    emit:
        // txt = blast_out_ch
        // seq = consensus_seq_ch
        fasta
        tsv = BLAST_BLASTN.out.tsv

        versions = ch_versions
        
}