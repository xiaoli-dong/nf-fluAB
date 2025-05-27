include {
    BCFTOOLS_NORM;
} from '../../modules/local/bcftools/norm/main.nf'

include {
    SNPEFF_SNPEFF
} from '../../modules/local/snpeff/snpeff'

include {
    BCFTOOLS_FILTER as BCFTOOLS_FILTER_FIX; //removes sites within 7 bases of an indel.
    BCFTOOLS_FILTER as BCFTOOLS_FILTER_LOW_QUALITY_DEPTH;
    BCFTOOLS_FILTER as BCFTOOLS_FILTER_FRAMESHIFT;
} from '../../modules/local/bcftools/filter'

include {
    SETGT;
} from '../../modules/local/setgt'

include {
    TABIX_BGZIPTABIX;
} from '../../modules/nf-core/tabix/bgziptabix'
//https://content.csbs.utah.edu/~rogers/tch/archgen/topics/weur.slr.html
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
        BCFTOOLS_FILTER_FIX(BCFTOOLS_NORM.out.vcf)
        ch_versions.mix(BCFTOOLS_FILTER_FIX.out.versions)

        BCFTOOLS_FILTER_LOW_QUALITY_DEPTH(BCFTOOLS_FILTER_FIX.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_FILTER_LOW_QUALITY_DEPTH.out.versions)
        //annotate vcf
        SNPEFF_SNPEFF(BCFTOOLS_FILTER_LOW_QUALITY_DEPTH.out.vcf, snpeff_db, snpeff_config, snpeff_dataDir )
        ch_versions = ch_versions.mix(SNPEFF_SNPEFF.out.versions)

        //
        //BCFTOOLS_FILTER_LOW_QUALITY_DEPTH(SNPEFF_SNPEFF.out.vcf)
        BCFTOOLS_FILTER_FRAMESHIFT(SNPEFF_SNPEFF.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_FILTER_FRAMESHIFT.out.versions)
        SETGT(BCFTOOLS_FILTER_FRAMESHIFT.out.vcf.join(BCFTOOLS_FILTER_FRAMESHIFT.out.tbi))
        ch_versions = ch_versions.mix(SETGT.out.versions)
        TABIX_BGZIPTABIX(SETGT.out.vcf)
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    emit:
        vcf_tbi = TABIX_BGZIPTABIX.out.gz_tbi
        versions = ch_versions

}
