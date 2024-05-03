include {
    BCFTOOLS_NORM as BCFTOOLS_NORM_BIALLELIC;
    BCFTOOLS_NORM as BCFTOOLS_NORM_MULTIALLELIC;
} from '../../modules/local/bcftools/norm/main.nf'

include {
    BCFTOOLS_VIEW;
} from '../../modules/local/bcftools/view/main.nf'

include {
    BCFTOOLS_SETGT;
} from '../../modules/local/bcftools/setGT/main.nf'

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
        
        
       
        /*
            Convert multiallelic to biallelic vcf first, output vcf.gz and index it
        */
        BCFTOOLS_NORM_BIALLELIC(vcf_tbi, fasta)
        ch_versions.mix(BCFTOOLS_NORM_BIALLELIC.out.versions)
       

        /*
        Filter out:
        non-variants sites and the sites whose alt is under certain value
        keep variants whose frequency >= 0.25 and depth >= mindepth
        */
        
        ch_input = BCFTOOLS_NORM_BIALLELIC.out.vcf.join(BCFTOOLS_NORM_BIALLELIC.out.tbi)
        BCFTOOLS_VIEW(ch_input, [], [], []) //output vcf.gz
        ch_versions.mix(BCFTOOLS_VIEW.out.versions)

        /*
        convert biallelic vcf to multiallelic vcf
        SNPs and indels shuld be merged separately into two records
        */
        
       
        BCFTOOLS_VIEW.out.vcf
            .join(BCFTOOLS_VIEW.out.tbi)
            .join(fasta)
            .multiMap{
                it ->
                    vcf_tbi: [it[0], it[1], it[2]]
                    fasta: [it[0], it[3]]
            }.set{
                ch_input
            }
        BCFTOOLS_NORM_MULTIALLELIC(ch_input.vcf_tbi, ch_input.fasta)
        

        //set genotype
        BCFTOOLS_SETGT(BCFTOOLS_NORM_MULTIALLELIC.out.vcf)
        ch_versions.mix(BCFTOOLS_SETGT.out.versions)

        BCFTOOLS_SETGT.out.vcf
            .join(BCFTOOLS_SETGT.out.tbi)
            .join(fasta)
            .join(mask_bed_file)
            .multiMap{
                it ->
                    vcf_tbi_fasta: [it[0], it[1], it[2], it[3]]
                    mask_bed: [it[0], it[4]]
            }.set{
                ch_input
            }
        BCFTOOLS_CONSENSUS(ch_input.vcf_tbi_fasta, ch_input.mask_bed)
        ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)
        
        //convert to tab format, change header, filter out 
        //contigs has too may Ns > 25% and then covert back to 
        //fasta  
        SEQKIT_FX2TAB_REFORMAT(BCFTOOLS_CONSENSUS.out.fasta)
        CONSENSUS_REHEADER(SEQKIT_FX2TAB_REFORMAT.out.text)
        SEQKIT_TAB2FX(CONSENSUS_REHEADER.out.txt)     

        //CONSENSUS_REHEADER(BCFTOOLS_CONSENSUS.out.fasta)
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
in freebayes, the variants has GT 
when the vaf is high: GT=1, 1 for alternate allele
when the vaf is low: GT = 0, o fr reference allele
For diploid organisms, it has 0 value for reference allele and 1 for the alternate allele (non-reference allele).
0/0	the sample is a homozygous reference
0/1	the sample is heterozygous (carries both reference and alternate alleles)
1/1	the sample is a homozygous alternate
./.	No genotype called or missing genotype

deletion just apply
*/
