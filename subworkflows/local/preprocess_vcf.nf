include {
    BCFTOOLS_NORM;
} from '../../modules/local/bcftools/norm/main.nf'

include { 
    SNPEFF_SNPEFF 
} from '../../modules/local/snpeff/snpeff'         

include {
    BCFTOOLS_FILTER as BCFTOOLS_FILTER_LOW_QUALITY_DEPTH;
    BCFTOOLS_FILTER as BCFTOOLS_FILTER_FRAMESHIFT;
} from '../../modules/local/bcftools/filter'

workflow PREPROCESS_VCF {   

    take:
        vcf_tbi
        fasta //[meta, fasta]
        snpeff_db
        snpeff_config
        snpeff_dataDir
    main:
        ch_versions = Channel.empty()
        
        //--check-ref w -m -any --output-type z  --write-index=tbi
        BCFTOOLS_NORM(vcf_tbi, fasta)
        ch_versions.mix(BCFTOOLS_NORM.out.versions)

        //annotate vcf
        SNPEFF_SNPEFF(BCFTOOLS_NORM.out.vcf, snpeff_db, snpeff_config, snpeff_dataDir )
        
        //
        BCFTOOLS_FILTER_LOW_QUALITY_DEPTH(SNPEFF_SNPEFF.out.vcf)
        BCFTOOLS_FILTER_FRAMESHIFT(BCFTOOLS_FILTER_LOW_QUALITY_DEPTH.out.vcf)

    emit:
        vcf = BCFTOOLS_FILTER_FRAMESHIFT.out.vcf
        tbi = BCFTOOLS_FILTER_FRAMESHIFT.out.tbi
        versions = ch_versions
        
}
