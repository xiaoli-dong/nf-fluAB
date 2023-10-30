
include { PROCESSGVCF  } from '../../modules/local/processgvcf'

include {
    BCFTOOLS_NORM as BCFTOOLS_NORM_LOWFREQ;
    BCFTOOLS_NORM as BCFTOOLS_NORM_HIGHFREQ;
} from '../../modules/nf-core/bcftools/norm'


include {
    BCFTOOLS_CONSENSUS as BCFTOOLS_CONSENSUS_AMBIGUOUS
} from '../../modules/nf-core/bcftools/consensus'
include {
    BCFTOOLS_CONSENSUS as BCFTOOLS_CONSENSUS_FIXED
} from '../../modules/local/bcftools/consensus'

include {
    TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_FIXED;
    TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_AMBIGUOUS;
} from '../../modules/nf-core/tabix/bgziptabix/main'


include {
    SPLITVCF as SPLITVCF_AMBIGUOUS;
    SPLITVCF as SPLITVCF_FIXED;
} from '../../modules/local/splitvcf'

include {
    SEQKIT_FX2TAB as SEQKIT_FX2TAB_CONSENSUS
} from '../../modules/nf-core/seqkit/fx2tab'

include { BIOAWK  } from '../../modules/nf-core/bioawk/main'


workflow make_consensus {   

    take:
        vcf //[meta, vcf_file]
        fasta //reference fasta file for mapping

    main:
        ch_versions = Channel.empty()

        /* Process a .gvcf file to create a file of consensus variants, low-frequency variants and a coverage mask
        // make depth mask, split variants into ambiguous/consensus
        // NB: this has to happen before bcftools norm or else the depth mask misses any bases exposed during normalization
        */
        PROCESSGVCF(vcf)
        lowFreq_vcf = PROCESSGVCF.out.variants
        highFreq_vcf = PROCESSGVCF.out.consensus

        lowFreq_vcf.join(fasta).multiMap{
            it ->
                input: [it[0], it[1], []]
                fasta: [it[0], it[2]]
        }.set{
            ch_input
        }
        
        /*
        ########################## 
        normalize variant records into canonical VCF representation ##############
        */
        BCFTOOLS_NORM_LOWFREQ(
            ch_input.input,
            ch_input.fasta
        )
        ch_versions = ch_versions.mix(BCFTOOLS_NORM_LOWFREQ.out.versions.first())

        highFreq_vcf.join(fasta).multiMap{
            it ->
                input: [it[0], it[1], []]
                fasta: [it[0], it[2]]
        }.set{
            ch_input
        }
        BCFTOOLS_NORM_HIGHFREQ(
            ch_input.input,
            ch_input.fasta
        )
        ch_versions = ch_versions.mix(BCFTOOLS_NORM_HIGHFREQ.out.versions.first())
        // ###################################################


        /*
        split the consensus sites file into a set that should be IUPAC codes 
        and all other bases, using the ConsensusTag in the VCF
        */
        SPLITVCF_AMBIGUOUS(BCFTOOLS_NORM_HIGHFREQ.out.vcf)
        SPLITVCF_FIXED(BCFTOOLS_NORM_HIGHFREQ.out.vcf)

        TABIX_BGZIPTABIX_AMBIGUOUS( SPLITVCF_AMBIGUOUS.out.vcf)
        TABIX_BGZIPTABIX_FIXED( SPLITVCF_FIXED.out.vcf )

        /* 
            apply ambiguous variants first using IUPAC codes. 
            this variant set cannot contain indels or the subsequent step will break
        */
        BCFTOOLS_CONSENSUS_AMBIGUOUS(
            //tuple val(meta), path(vcf), path(tbi), path(fasta)
           TABIX_BGZIPTABIX_AMBIGUOUS.out.gz_tbi.join(fasta)
        )

        TABIX_BGZIPTABIX_FIXED.out.gz_tbi
            .join(BCFTOOLS_CONSENSUS_AMBIGUOUS.out.fasta)
            .join(PROCESSGVCF.out.mask)
            .multiMap{
                it ->
                input: [it[0], it[1], it[2], it[3]]
                mask: [it[0], it[4]]

            }.set{
                ch_input
            }   
        //# apply remaninng variants, including indels
        BCFTOOLS_CONSENSUS_FIXED(
            ch_input.input,  //tuple val(meta), path(vcf), path(tbi), path(fasta)
            ch_input.mask //tuple val(meta), path(mask)

        )
        
        //format the consensus contig name
        BIOAWK(BCFTOOLS_CONSENSUS_FIXED.out.fasta)
        ch_versions = ch_versions.mix(BIOAWK.out.versions.first())

        //sort sequence by name before stats??????
        SEQKIT_FX2TAB_CONSENSUS(BIOAWK.out.output)
        ch_versions = ch_versions.mix(SEQKIT_FX2TAB_CONSENSUS.out.versions.first())

    emit:
       
        fasta = BIOAWK.out.output
        stats = SEQKIT_FX2TAB_CONSENSUS.out.text
        versions = ch_versions
        
}