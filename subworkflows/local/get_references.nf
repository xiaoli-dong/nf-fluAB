include {
    MASH_SCREEN
} from '../../modules/nf-core/mash/screen/main'

include {
    FILTERMASH;    
} from '../../modules/local/misc'

include {
    SEQKIT_GREP
} from '../../modules/local/seqkit/grep/main'

include {
    SEQKIT_HEADER
} from '../../modules/local/seqkit/header' 

include {
    SAMTOOLS_FAIDX
} from '../../modules/nf-core/samtools/faidx/main'

workflow GETREF_BY_MASH {   

    take:
        illumina_reads
        msh_db
        fasta_db

    main:
        ch_versions = Channel.empty()
        
        MASH_SCREEN(illumina_reads, msh_db)
        ch_versions = ch_versions.mix(MASH_SCREEN.out.versions.first())

        // only keep the row with the best identity
        FILTERMASH(
            MASH_SCREEN.out.screen.filter{meta, screen -> screen.countLines() > 0},
            params.min_shared_hashes
        )
        ch_versions = ch_versions.mix(FILTERMASH.out.versions.first())
          
        // [identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment]
        FILTERMASH.out.screen.map{
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
        
        SEQKIT_HEADER(fasta)
        header = SEQKIT_HEADER.out.txt
        //SEQKIT_GREP.out.fasta.view()
        
    emit:
        screen = FILTERMASH.out.screen
        fasta_fai 
        header
        versions = ch_versions
        
}