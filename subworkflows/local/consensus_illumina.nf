
include {
    BEDTOOLS_MASKFASTA
} from '../../modules/nf-core/bedtools/maskfasta/main.nf'


include {
    BEDTOOLS_GENOMECOV
} from '../../modules/nf-core/bedtools/genomecov/main.nf'


include {
    BCFTOOLS_MPILEUP
} from '../../modules/nf-core/bcftools/mpileup/main.nf'

include {
    BCFTOOLS_NORM as MPILEUP_NORM
} from '../../modules/nf-core/bcftools/norm/main.nf'

include {
    TABIX_BGZIPTABIX as BGZIPTABIX_NORM;
    TABIX_BGZIPTABIX as BGZIPTABIX_FILTER;
} from '../../modules/nf-core/tabix/bgziptabix/main'

/* include {
    BCFTOOLS_QUERY
} from '../../modules/nf-core/bcftools/query'
 */
include {
    BCFTOOLS_FILTER
} from '../../modules/nf-core/bcftools/filter'

include {
    BCFTOOLS_CONSENSUS
} from '../../modules/nf-core/bcftools/consensus'

include {
    SEQKIT_FX2TAB as SEQKIT_FX2TAB_CONSENSUS
} from '../../modules/nf-core/seqkit/fx2tab'

include {
    FORMAT_CONSENSUS;
} from '../../modules/local/misc'

workflow consensus_bcftools {   

    take:
        bam //[meta, bam]
        fasta //[meta, fasta]
    main:
        ch_versions = Channel.empty()
        ch_input = bam.map{
            it -> [it[0], it[1], 1] //[meta, bam, scale]
        }

        //mask low depth regions of the reference
        BEDTOOLS_GENOMECOV(ch_input, [], "bed")
        ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)
        BEDTOOLS_GENOMECOV.out.genomecov.join(fasta).multiMap{
            it ->
                bed: [it[0], it[1]]
                fasta: it[2]
        }.set{
            ch_input
        }

        BEDTOOLS_MASKFASTA(ch_input.bed, ch_input.fasta)
        ch_versions.mix(BEDTOOLS_MASKFASTA.out.versions)
        bam.join(BEDTOOLS_MASKFASTA.out.fasta).multiMap{
            it ->
                bam_interval: [it[0], it[1], []]
                fasta: [it[0], it[2]]
        }.set{
            ch_input
        }
        
        BCFTOOLS_MPILEUP(ch_input.bam_interval, ch_input.fasta, false)
        ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)
        
        BCFTOOLS_MPILEUP.out.vcf
            .join(BCFTOOLS_MPILEUP.out.tbi)
            .join(BEDTOOLS_MASKFASTA.out.fasta)
            .multiMap{
                it ->
                    vcf_tbi: [it[0], it[1], it[2]] //[meta, vcf.gz, vcf.gz.tbi]
                    fasta: [it[0], it[3]] 
            }
            .set{
                ch_input
            }

        MPILEUP_NORM(ch_input.vcf_tbi, ch_input.fasta)
        ch_versions.mix(MPILEUP_NORM.out.versions)

       /*  //created bed file to mask low depth region in the consensus calling
        BGZIPTABIX_NORM(MPILEUP_NORM.out.vcf)
        BCFTOOLS_QUERY(BGZIPTABIX_NORM.out.gz_tbi, [], [], [])
        ch_versions.mix(BGZIPTABIX_NORM.out.versions)
        ch_versions.mix(BCFTOOLS_QUERY.out.versions)
 */
        //filter out variants: keep variant 0.25 < maf < 0.75, dp > 10
        BCFTOOLS_FILTER(MPILEUP_NORM.out.vcf) //[meta, vcf]
        ch_versions.mix(BCFTOOLS_FILTER.out.versions)

        BGZIPTABIX_FILTER(BCFTOOLS_FILTER.out.vcf)
        ch_versions.mix(BGZIPTABIX_FILTER.out.versions)

        ch_vcf_tbi_fasta = BGZIPTABIX_FILTER.out.gz_tbi
            .join(BEDTOOLS_MASKFASTA.out.fasta)
            //.join(BCFTOOLS_QUERY.out.output)
            .map{
                it -> vcf_tbi_fasta: [it[0], it[1], it[2], it[3]]
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

include {   FREEBAYES   } from '../../modules/nf-core/freebayes/main.nf'
include {
    BCFTOOLS_NORM as FREEBAYES_NORM
} from '../../modules/nf-core/bcftools/norm/main.nf'

include {
    TABIX_TABIX
} from '../../modules/nf-core/tabix/tabix'

include {
    BCFTOOLS_SORT
} from '../../modules/nf-core/bcftools/sort'
include {
    SAMTOOLS_FAIDX as INDEX_REFERENCE
} from '../../modules/nf-core/samtools/faidx/main'

workflow consensus_freebayes {   

    take:
        bam_bai //[meta, bam, bai, [], [],[]]
        fasta_fai // [meta, fasta, fai]
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

        bam_bai.join(fasta_fai)
            .multiMap{
                it ->
                    bam_bai: [it[0], it[1], it[2], [], [], []]
                    fasta: it[3]
                    fai: [it[4]]
            }
            .set{
                ch_input
            }

        FREEBAYES (ch_input.bam_bai, ch_input.fasta, ch_input.fai, [], [], [])
        ch_versions.mix(FREEBAYES.out.versions)
       
        /*
        normalize the vcf
        */
        //without sorting it gives error: [E::hts_idx_push] Unsorted positions on seque
        BCFTOOLS_SORT(FREEBAYES.out.vcf)
        TABIX_TABIX(BCFTOOLS_SORT.out.vcf)
        BCFTOOLS_SORT.out.vcf //vcf.gz file
            .join(TABIX_TABIX.out.tbi)
            .join(fasta_fai)
            .multiMap{
                it ->
                    vcf_tbi: [it[0], it[1], it[2]] //[meta, vcf.gz, vcf.gz.tbi]
                    fasta: [it[0], it[3]] 
            }
            .set{
                ch_input
            }
        
        FREEBAYES_NORM(ch_input.vcf_tbi, ch_input.fasta)
        ch_versions.mix(FREEBAYES_NORM.out.versions)
/* 
        //created bed file to mask low depth region in the consensus calling
        BGZIPTABIX_NORM(FREEBAYES_NORM.out.vcf)
        BCFTOOLS_QUERY(BGZIPTABIX_NORM.out.gz_tbi, [], [], [])
        ch_versions.mix(BGZIPTABIX_NORM.out.versions)
        ch_versions.mix(BCFTOOLS_QUERY.out.versions) */
        //TODO: 
        //filter out variants: keep variant 0.25 < maf < 0.75, dp > 10
        BCFTOOLS_FILTER(FREEBAYES_NORM.out.vcf) //[meta, vcf]
        ch_versions.mix(BCFTOOLS_FILTER.out.versions)

        BGZIPTABIX_FILTER(BCFTOOLS_FILTER.out.vcf)
        ch_versions.mix(BGZIPTABIX_FILTER.out.versions)

        ch_vcf_tbi_fasta = BGZIPTABIX_FILTER.out.gz_tbi.join(fasta_fai).map{
            it -> [it[0], it[1], it[2], it[3]]
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
