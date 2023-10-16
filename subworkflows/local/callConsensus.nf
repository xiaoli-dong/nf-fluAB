
include { BIOAWK  } from '../../modules/nf-core/bioawk/main'
include { FREEBAYES } from '../../modules/nf-core/freebayes/main.nf'
include { PROCESSGVCF  } from '../../modules/local/processgvcf'
include {BCFTOOLS_NORM as BCFTOOLS_NORM_VARIANTS} from '../../modules/nf-core/bcftools/norm'
include {BCFTOOLS_NORM as BCFTOOLS_NORM_CONSENSUS} from '../../modules/nf-core/bcftools/norm'
include {BCFTOOLS_CONSENSUS as BCFTOOLS_CONSENSUS_AMBIGUOUS} from '../../modules/nf-core/bcftools/consensus'
include {BCFTOOLS_CONSENSUS as BCFTOOLS_CONSENSUS_FIXED} from '../../modules/local/bcftools/consensus'
include {TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_VARIANTS} from '../../modules/nf-core/tabix/bgziptabix/main'
include {TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_CONSENSUS} from '../../modules/nf-core/tabix/bgziptabix/main'
include {TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_CONSENSUS_FIXED} from '../../modules/nf-core/tabix/bgziptabix/main'
include {TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_CONSENSUS_AMBIGUOUS} from '../../modules/nf-core/tabix/bgziptabix/main'
include {SPLITVCF as SPLITVCF_CONSENSUS_AMBIGUOUS} from '../../modules/local/splitvcf'
include {SPLITVCF as SPLITVCF_CONSENSUS_FIXED} from '../../modules/local/splitvcf'
include {SEQKIT_FX2TAB as SEQKIT_FX2TAB_CONSENSUS} from '../../modules/nf-core/seqkit/fx2tab'


workflow callConsensus {   

    take:
        bam
        bai
        ref_fasta/*  */
        ref_fasta_fai
    main:
        ch_versions = Channel.empty()

        input = bam.join(bai).join(ref_fasta).join(ref_fasta_fai)
        // input.view()
        input.multiMap{ it ->
            finput: [it[0], it[1], it[2], [], [], []] //meta, bam, bam.bai
            fasta: it[3] //fasta ref_fastaerence file
            fai: it[4] //index file of the fasta file
        }.set{
            input_ch
        }
       

        samples     = []
        populations = []
        cnv         = []
        FREEBAYES (input_ch.finput, input_ch.fasta, input_ch.fai, [], [], [])
        ch_versions = ch_versions.mix(FREEBAYES.out.versions.first())
        //FREEBAYES.out.vcf.view()
        
        // Process a .gvcf file to create a file of consensus variants, low-frequency variants and a coverage mask
        // make depth mask, split variants into ambiguous/consensus
        // NB: this has to happen before bcftools norm or else the depth mask misses any bases exposed during normalization
        PROCESSGVCF(FREEBAYES.out.vcf)
        
        //create vcf tbi  
        TABIX_BGZIPTABIX_VARIANTS(PROCESSGVCF.out.variants)
        TABIX_BGZIPTABIX_VARIANTS.out.gz_tbi.join(ref_fasta).multiMap{
                it ->
                input: [it[0], it[1], it[2]]
                fasta: it[3]
            }.set{
                input_bcftool_ch
            }
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_VARIANTS.out.versions.first())

        // ############# normalize variant records into canonical VCF representation ##############
        // low frequency variant
        BCFTOOLS_NORM_VARIANTS(
            input_bcftool_ch.input,
            input_bcftool_ch.fasta
        )
        ch_versions = ch_versions.mix(BCFTOOLS_NORM_VARIANTS.out.versions.first())

        TABIX_BGZIPTABIX_CONSENSUS(PROCESSGVCF.out.consensus)
        TABIX_BGZIPTABIX_CONSENSUS.out.gz_tbi.join(ref_fasta).multiMap{
                it ->
                input: [it[0], it[1], it[2]]
                fasta: it[3]
            }.set{
                input_bcftool_consensus_ch
            }
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_CONSENSUS.out.versions.first())
        
        // high-frequency variants
        BCFTOOLS_NORM_CONSENSUS(
            input_bcftool_consensus_ch.input,
            input_bcftool_consensus_ch.fasta
        )
        ch_versions = ch_versions.mix(BCFTOOLS_NORM_CONSENSUS.out.versions.first())
        // ###################################################


        //split the consensus sites file into a set that should be IUPAC codes and all other bases, using the ConsensusTag in the VCF
        SPLITVCF_CONSENSUS_AMBIGUOUS(BCFTOOLS_NORM_CONSENSUS.out.vcf)
        TABIX_BGZIPTABIX_CONSENSUS_AMBIGUOUS(SPLITVCF_CONSENSUS_AMBIGUOUS.out.vcf)
        

        input_ch =  TABIX_BGZIPTABIX_CONSENSUS_AMBIGUOUS.out.gz_tbi.join(ref_fasta)
        //input_ch.view()

        // # apply ambiguous variants first using IUPAC codes. this variant set cannot contain indels or the subsequent step will break
        // consensus
        BCFTOOLS_CONSENSUS_AMBIGUOUS(
            input_ch
        )
        
        SPLITVCF_CONSENSUS_FIXED(BCFTOOLS_NORM_CONSENSUS.out.vcf)
        TABIX_BGZIPTABIX_CONSENSUS_FIXED(SPLITVCF_CONSENSUS_FIXED.out.vcf)
        input_ch =  TABIX_BGZIPTABIX_CONSENSUS_FIXED.out.gz_tbi.join(BCFTOOLS_CONSENSUS_AMBIGUOUS.out.fasta).join(PROCESSGVCF.out.mask)
        //input_ch.view()
        
        //# apply remaninng variants, including indels
        input_ch.multiMap{
            it ->
            input: [it[0], it[1], it[2], it[3]]
            mask: [it[0], it[4]]

        }.set{
            input_bcftool_consensus_ch
        }
        BCFTOOLS_CONSENSUS_FIXED(
            input_bcftool_consensus_ch.input,
            input_bcftool_consensus_ch.mask

        )
        // BCFTOOLS_CONSENSUS.out.fasta.view()
        
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