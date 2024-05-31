include {
    BCFTOOLS_NORM;
} from '../../modules/local/bcftools/norm/main.nf'

include {
    SETGT;
} from '../../modules/local/setgt.nf'

include {
     TABIX_BGZIPTABIX
} from '../../modules/nf-core/tabix/bgziptabix'


include {
    BCFTOOLS_CONSENSUS;
} from '../../modules/local/bcftools/consensus'

include {
    SEQKIT_FX2TAB as SEQKIT_FX2TAB_CONSENSUS;
    SEQKIT_FX2TAB as SEQKIT_FX2TAB_REFORMAT;
} from '../../modules/nf-core/seqkit/fx2tab'

include {
    SEQKIT_SORT
} from '../../modules/local/seqkit/sort'

include {
    CONSENSUS_REHEADER;
} from '../../modules/local/consensus/reheader'

include {
    SEQKIT_TAB2FX 
} from '../../modules/nf-core/seqkit/tab2fx/main'   

workflow CONSENSUS {   

    take:
        vcf_tbi
        fasta //[meta, fasta]
        mask_bed_file
    main:
        ch_versions = Channel.empty()
        
        BCFTOOLS_NORM(vcf_tbi, fasta)
        ch_versions.mix(BCFTOOLS_NORM.out.versions)

        /*
            filter vcf: 
                for the biallelic variants, ignore the variants whose freq < 0.25
                for the multiallelic variants, ignore the variants whose freq < 0.20
            set genotype for each record
        */
        ch_input = BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi)

        SETGT(ch_input)
        ch_versions.mix(SETGT.out.versions)

        TABIX_BGZIPTABIX(SETGT.out.vcf)
        ch_versions.mix(TABIX_BGZIPTABIX.out.versions)
        
        TABIX_BGZIPTABIX.out.gz_tbi
            .join(fasta)
            .join(mask_bed_file)
            //.join(SETGT.out.low_depth_bed)
            .multiMap{
                it ->
                    vcf_tbi_fasta: [it[0], it[1], it[2], it[3]]
                    mask_bed: [it[0], it[4]]
            }.set{
                ch_input
            }
        BCFTOOLS_CONSENSUS(ch_input.vcf_tbi_fasta, ch_input.mask_bed)
        ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)
        
        /*
        convert to tab format, change header, filter out 
        contigs has too may Ns > 25% and then covert back to 
        fasta
        */  
        SEQKIT_FX2TAB_REFORMAT(BCFTOOLS_CONSENSUS.out.fasta)
        CONSENSUS_REHEADER(SEQKIT_FX2TAB_REFORMAT.out.text)
        SEQKIT_TAB2FX(CONSENSUS_REHEADER.out.txt)     
        ch_versions.mix(CONSENSUS_REHEADER.out.versions)
        
        SEQKIT_SORT(SEQKIT_TAB2FX.out.fastx)
        ch_versions = ch_versions.mix(SEQKIT_SORT.out.versions)
        fasta = SEQKIT_SORT.out.fasta
        
        SEQKIT_FX2TAB_CONSENSUS(fasta)
        ch_versions = ch_versions.mix(SEQKIT_FX2TAB_CONSENSUS.out.versions)
        stats = SEQKIT_FX2TAB_CONSENSUS.out.text
    emit:
        fasta
        stats
        versions = ch_versions
        
}

/*
Phased variants (heterozygous and homozygous) contain a modified GT field, using a pipe symbol (|) instead of forward slash in accordance with VCF 4.1 specifications. 
in freebayes, the variants has GT
for unphased heterozygous mutations, there is no difference between 0/1 and 1/0, by convention it is written as 0/1
when the vaf is high: GT=1, 1 for alternate allele
when the vaf is low: GT = 0, o fr reference allele
For diploid organisms, it has 0 value for reference allele and 1 for the alternate allele (non-reference allele).
0/0	the sample is a homozygous reference
0/1	the sample is heterozygous (carries both reference and alternate alleles)
1/1	the sample is a homozygous alternate
./.	No genotype called or missing genotype

deletion just apply
*/
