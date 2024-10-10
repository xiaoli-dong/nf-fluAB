include {
    MASH_SCREEN
} from '../../modules/nf-core/mash/screen/main'

include {
   MASH_FILTER;    
} from '../../modules/local/mash/filter'

include {
    SEQKIT_GREP
} from '../../modules/local/seqkit/grep/main'

include {
    SEQKIT_SEQ
} from '../../modules/local/seqkit/seq' 

include {
    SAMTOOLS_FAIDX
} from '../../modules/nf-core/samtools/faidx/main'

workflow SEEK_REFERENCES {   
    take:
        reads
        msh_db
        fasta_db
    main:
        ch_versions = Channel.empty()
        screen_best = Channel.empty()

        MASH_SCREEN(reads, msh_db)
        ch_versions = ch_versions.mix(MASH_SCREEN.out.versions.first())

        /*
        if it has hits to multiple types, keep the rows with the best identity for each type
        for example, if the sample has hits with both H1, H3. Then keep the best records of 
        both types
        */
        MASH_FILTER(
            MASH_SCREEN.out.screen.filter{meta, screen -> screen.countLines() > 0}
        )
        ch_versions = ch_versions.mix(MASH_FILTER.out.versions.first())

        MASH_FILTER.out.screen.filter{
            meta, screen -> screen.countLines() > 0
        }.set{
            screen_best
        } 
        // [identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment]
       screen_best.map{
            it ->
                def meta = it[0]
                def arr = it[1].splitCsv(skip: 1, sep: '\t')
                a = []
                arr.each{ n -> a.add(n[4]) }
                [meta, a.join(',')]
        }.set {seqids}
        
        //seqids.view()

        SEQKIT_GREP(seqids, fasta_db)
        ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions.first())
        fasta = SEQKIT_GREP.out.fasta

        SAMTOOLS_FAIDX(fasta)
        fasta_fai = fasta.join(SAMTOOLS_FAIDX.out.fai)
        
        SEQKIT_SEQ(fasta)
        header = SEQKIT_SEQ.out.txt
        //SEQKIT_GREP.out.fasta.view()
        
    emit:
        screen = screen_best
        fasta_fai 
        header
        versions = ch_versions
        
}