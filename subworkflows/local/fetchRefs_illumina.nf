
include { MASH_SCREEN                      } from '../../modules/nf-core/mash/screen/main'
include {FILTERMASH} from '../../modules/local/filtermash'
include { SEQKIT_GREP                      } from '../../modules/local/seqkit/grep'

workflow fetchRefs_illumina {   

    take:
        reads
        mash_db
        mash_db_fasta

    main:
        ch_versions = Channel.empty()
        //reads.view()
        //mash_db.view()
        //screen reads for containment of influenza segments
        MASH_SCREEN(reads,mash_db)
        ch_versions = ch_versions.mix(MASH_SCREEN.out.versions.first())
        
        //filter out empty channel
        FILTERMASH(
            MASH_SCREEN.out.screen.filter{meta, screen -> screen.countLines() > 0}
        )
        ch_versions = ch_versions.mix(FILTERMASH.out.versions.first())
        
        FILTERMASH.out.screen
            .filter{meta, tsv -> tsv.countLines() > 0 }
            .map{
                it ->
                def meta = it[0]
                def arr = it[1].splitCsv(sep: '\t')
                a = []
                arr.each{
                    n ->
                    a.add(n[4])
                }
                [meta, a.join(',')]
            }.set {accession_list}

        SEQKIT_GREP(
            accession_list,
            mash_db_fasta

        )
        ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions.first())

    emit:
        screen = FILTERMASH.out.screen
        fasta = SEQKIT_GREP.out.fasta
        versions = ch_versions
        
}