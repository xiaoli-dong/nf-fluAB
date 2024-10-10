
include {
    BCFTOOLS_CONSENSUS;
} from '../../modules/local/bcftools/consensus'
include {
    CONSENSUS_REFORMAT;
} from '../../modules/local/consensus/reformat/main'

include {
    SEQKIT_FX2TAB;
} from '../../modules/nf-core/seqkit/fx2tab'


workflow CONSENSUS {   

    take:
        vcf_tbi_fasta
        mask_bed_file
    main:
        ch_versions = Channel.empty()
        
        
        BCFTOOLS_CONSENSUS(vcf_tbi_fasta, mask_bed_file)
        ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)
        
        /*
        reformatted header, elimnated segments with too many
        ambiguious bases
        */  
      
        CONSENSUS_REFORMAT(BCFTOOLS_CONSENSUS.out.fasta)
        ch_versions.mix(CONSENSUS_REFORMAT.out.versions)
        fasta = CONSENSUS_REFORMAT.out.fasta

        //stats
        SEQKIT_FX2TAB(fasta)
        ch_versions = ch_versions.mix(SEQKIT_FX2TAB.out.versions)

        stats = SEQKIT_FX2TAB.out.text
    emit:
        fasta
        stats
        versions = ch_versions
        
}
