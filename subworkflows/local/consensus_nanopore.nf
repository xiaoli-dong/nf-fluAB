
include {
    BEDTOOLS_MASKFASTA
} from '../../modules/nf-core/bedtools/maskfasta/main.nf'


include {
    BEDTOOLS_GENOMECOV
} from '../../modules/nf-core/bedtools/genomecov/main.nf'

include {
    SAMTOOLS_FAIDX as INDEX_REFERENCE
} from '../../modules/nf-core/samtools/faidx/main'

include {
    TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_NORM;
    
} from '../../modules/nf-core/tabix/bgziptabix/main'

include {
    BCFTOOLS_CONSENSUS
} from '../../modules/nf-core/bcftools/consensus'

include {
    SEQKIT_FX2TAB as SEQKIT_FX2TAB_CONSENSUS
} from '../../modules/nf-core/seqkit/fx2tab'

include {
    FORMAT_CONSENSUS;
} from '../../modules/local/misc'

include {   BCFTOOLS_NORM} from '../../modules/nf-core/bcftools/norm/main.nf'

include {
    CLAIR3
} from '../../modules/local/clair3/main'

workflow consensus_clair3 {   

    take:
        bam_bai //[meta, bam, bai, [], [],[]]
        fasta_fai // [meta, fasta, fai]
        model_path
    main:
        ch_versions = Channel.empty()

        ch_input = bam_bai.map{
            it -> [it[0], it[1], 1] //[meta, bam, scale]
        }
        //mask low depth regions of the reference
        BEDTOOLS_GENOMECOV(ch_input, [], "bed")
        ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)

        BEDTOOLS_GENOMECOV.out.genomecov.join(fasta_fai).multiMap{
            it ->
                bed: [it[0], it[1]]
                fasta: it[2]
        }.set{
            ch_input
        }

        BEDTOOLS_MASKFASTA(ch_input.bed, ch_input.fasta)
        ch_versions.mix(BEDTOOLS_MASKFASTA.out.versions)

        
        INDEX_REFERENCE(BEDTOOLS_MASKFASTA.out.fasta)
        ch_versions.mix(INDEX_REFERENCE.out.versions)
        fasta_fai = BEDTOOLS_MASKFASTA.out.fasta.join(INDEX_REFERENCE.out.fai)
        //bam_bai.view()
        //fasta_fai.view()
        bam_bai.join(fasta_fai)
            .multiMap{
                it ->
                    bam_bai: [it[0], it[1], it[2], [], [], []]
                    fasta_fai: [it[0], it[3], it[4]]
            }
            .set{
                ch_input
            }

        CLAIR3(
            ch_input.bam_bai, 
            ch_input.fasta_fai, 
            model_path
        )
        ch_versions = ch_versions.mix(CLAIR3.out.versions)

        CLAIR3.out.vcf_tbi.join(fasta_fai).multiMap{
            it ->
                vcf_tbi: [it[0], it[1], it[2]]
                fasta: [it[0], it[3]]
        }.set{
            ch_input
        }

        //TODO: add filter, 

        BCFTOOLS_NORM(ch_input.vcf_tbi, ch_input.fasta)
        TABIX_BGZIPTABIX_NORM(BCFTOOLS_NORM.out.vcf)
        ch_vcf_tbi_fasta = TABIX_BGZIPTABIX_NORM.out.gz_tbi.join(fasta_fai).map{
            it ->
                [it[0], it[1], it[2], it[3]]
        }
        

        BCFTOOLS_CONSENSUS(ch_vcf_tbi_fasta)
        ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)

        FORMAT_CONSENSUS(BCFTOOLS_CONSENSUS.out.fasta)
        ch_versions.mix(FORMAT_CONSENSUS.out.versions)

        SEQKIT_FX2TAB_CONSENSUS(FORMAT_CONSENSUS.out.consensus_fasta)

        ch_versions = ch_versions.mix(SEQKIT_FX2TAB_CONSENSUS.out.versions)
       
    emit:
        fasta = FORMAT_CONSENSUS.out.consensus_fasta
        stats = SEQKIT_FX2TAB_CONSENSUS.out.text
        versions = ch_versions
        
}
