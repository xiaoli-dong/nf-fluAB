include {
    BCFTOOLS_NORM as BCFTOOLS_NORM_BIALLELIC;
    BCFTOOLS_NORM as BCFTOOLS_NORM_MULTIALLELIC;
} from '../../modules/nf-core/bcftools/norm/main.nf'

include {
    BCFTOOLS_VIEW;
} from '../../modules/nf-core/bcftools/view/main.nf'

include {
    BCFTOOLS_SETGT;
} from '../../modules/local/bcftools/setGT/main.nf'

include {
    TABIX_TABIX as TABIX_TABIX_VCF;
    TABIX_TABIX as TABIX_TABIX_BIALLELIC;
    TABIX_TABIX as TABIX_TABIX_BIALLELIC_FILTERED;
    TABIX_TABIX as TABIX_TABIX_SETGT;
} from '../../modules/nf-core/tabix/tabix'

include {
    BCFTOOLS_CONSENSUS;
} from '../../modules/nf-core/bcftools/consensus'

include {
    SEQKIT_FX2TAB as SEQKIT_FX2TAB_CONSENSUS;
     SEQKIT_FX2TAB as SEQKIT_FX2TAB_REFORMAT;
} from '../../modules/nf-core/seqkit/fx2tab'

include {
    SEQKIT_SORT
} from '../../modules/nf-core/seqkit/sort'

include {
    REHEADER;
} from '../../modules/local/misc'
include { SEQKIT_TAB2FX } from '../../modules/nf-core/seqkit/tab2fx/main'    
workflow consensus {   

    take:
        vcf
        fasta //[meta, fasta]
    main:
        ch_versions = Channel.empty()
        
        /*
        handle empty vcf file
        */
        /*
        vcf.filter{
            meta, sample_vcf ->
            def line_count = 0
            def num_variants = 0
            def enough_variants = false
            // make sure that the VCF has at least 1 variant, then stop counting
            InputStream fileStream = new FileInputStream(sample_vcf.toFile());
            InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream);
            BufferedReader br = new BufferedReader(new InputStreamReader(gzipStream));
            
            while ((line = br.readLine()) != null) {
                if (! line.startsWith("#")) num_variants++
                if (num_variants > 0) {
                    enough_variants = true
                    break
                    }
                line_count++
            }   
            println ">>> Read ${line_count} lines from file ${sample_vcf} and found ${num_variants} variants"
            if(! enough_variants) println ">>> WARNING: File ${sample_vcf} does not have enough variants and will not be included"
            
            enough_variants
        }.set{
            filtered_vcf
        }
        */
        TABIX_TABIX_VCF(vcf)
        
        //vcf.join(TABIX_TABIX_VCF.out.tbi).join(fasta).view()
        vcf.join(TABIX_TABIX_VCF.out.tbi).join(fasta).multiMap{
            it ->
                vcf_tbi: [it[0], it[1], it[2]] //[meta, vcf.gz, vcf.gz.tbi]
                fasta: [it[0], it[3]] 
        }.set{
            ch_input
        } 
        
        //Convert multiallelic to biallelic vcf first, output vcf.gz and index it
        BCFTOOLS_NORM_BIALLELIC(ch_input.vcf_tbi, ch_input.fasta)
        ch_versions.mix(BCFTOOLS_NORM_BIALLELIC.out.versions)
        TABIX_TABIX_BIALLELIC(BCFTOOLS_NORM_BIALLELIC.out.vcf)
        ch_versions.mix(TABIX_TABIX_BIALLELIC.out.versions)

        /*
        Filter out:
        non-variants sites and the sites whose alt is under certain value
        keep variants whose frequency >= 0.25 and depth >= mindepth
        */
        ch_input = BCFTOOLS_NORM_BIALLELIC.out.vcf.join(TABIX_TABIX_BIALLELIC.out.tbi)
        BCFTOOLS_VIEW(ch_input, [], [], []) //output vcf.gz
        ch_versions.mix(BCFTOOLS_VIEW.out.versions)

        /*
        convert bioallelic vcf to multiallelic vcf
        SNPs and indels shuld be merged separately into two records
        */
        TABIX_TABIX_BIALLELIC_FILTERED(BCFTOOLS_VIEW.out.vcf)
        
        BCFTOOLS_VIEW.out.vcf
            .join(TABIX_TABIX_BIALLELIC_FILTERED.out.tbi)
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

        TABIX_TABIX_SETGT(BCFTOOLS_SETGT.out.vcf)
        BCFTOOLS_CONSENSUS(BCFTOOLS_SETGT.out.vcf.join(TABIX_TABIX_SETGT.out.tbi).join(fasta))
        ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)
        
        SEQKIT_FX2TAB_REFORMAT(BCFTOOLS_CONSENSUS.out.fasta)
        REHEADER(SEQKIT_FX2TAB_REFORMAT.out.text)

        //reformated and filtered fasta
        SEQKIT_TAB2FX(REHEADER.out.txt)     

        //REHEADER(BCFTOOLS_CONSENSUS.out.fasta)
        ch_versions.mix(REHEADER.out.versions)
        SEQKIT_SORT(SEQKIT_TAB2FX.out.fastx)
        //SEQKIT_SORT(REHEADER.out.consensus_fasta)
        ch_versions = ch_versions.mix(SEQKIT_SORT.out.versions)
        
        SEQKIT_FX2TAB_CONSENSUS(SEQKIT_SORT.out.fastx)
        ch_versions = ch_versions.mix(SEQKIT_FX2TAB_CONSENSUS.out.versions)
       

    emit:
        fasta = SEQKIT_SORT.out.fastx
        stats = SEQKIT_FX2TAB_CONSENSUS.out.text
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


/*
bcftools +setGT mytest.vcf.gz -- -t q -i 'GT="1" && FORMAT/ALT_FREQ < 0.9' -n 'c:0/1' | \
bcftools +setGT -o mytest.setGT.vcf.gz -- -t q -i 'GT="1" && FORMAT/ALT_FREQ >= 0.9' -n 'c:1/1'
tabix -p vcf mytest.setGT.vcf.gz
bcftools consensus -p test_ -f ${REF} --mark-del '-' -m $MASK -H I -i 'FORMAT/ALT_QUAL >= 20 & INFO/DP >= 10' mytest.setGT.vcf.gz > mytest.con.fasta
*/


/*
bcftools +setGT S10/variants/freebayes/S10.norm.vcf -- -t q -i 'strlen(ref) !=strlen(alt) && FORMAT/AD[:1]/FORMAT/DP >= 0.5' -n 'c:1/1' | bcftools +setGT  -- -t q -i 'strlen(ref) !=strlen(alt) && FORMAT/AD[:1]/FORMAT/DP < 0.5' -n 'c:0/0'  | bcftools +setGT -- -t q -i 'GT="1" && FORMAT/AD[:1]/FORMAT/DP < 0.75' -n 'c:0/1' | bcftools +setGT -o  mytest.setGT.vcf.gz -- -t q -i 'GT="1" && FORMAT/AD[:1]/FORMAT/DP >=0.75'  -n 'c:1/1' 
bcftools consensus -p test_ -f S10/variants/freebayes/S10.reference_masked.fa --mark-del "-" --mark-ins lc  --mark-snv lc -I mytest.setGT.vcf.gz > S10.test.consensus.fa

*/


/*
#Convert multiallelic to biallelic vcf first
bcftools norm --fasta-ref S10.reference_masked.fa -m- --check-ref w --output-type v S10.mpieup.vcf.gz --output test/S10.biallellic.vcf

#Filter the alternative alleles under certain value
bcftools view -e 'FORMAT/AD[:1]/FORMAT/DP < 0.25'  test/S10.biallellic.vcf > test/S10.biallellic-filtered.vcf

#Convert biallelic vcf to multiallelic vcf
bcftools norm --fasta-ref S10.reference_masked.fa -m+ --check-ref w --output-type v test/S10.biallellic-filtered.vcf --output test/S10.multiallellic-filtered.vcf

#set genotype
bcftools +setGT test/S10.multiallellic-filtered.vcf -- -t q -i 'strlen(ref) !=strlen(alt) && FORMAT/AD[:1]/FORMAT/DP >= 0.5' -n 'c:1/1' | bcftools +setGT  -- -t q -i 'strlen(ref) !=strlen(alt) && FORMAT/AD[:1]/FORMAT/DP < 0.5' -n 'c:0/0'  | bcftools +setGT -- -t q -i 'GT !~ "/"  && FORMAT/AD[:1]/FORMAT/DP < 0.75' -n 'c:0/1' | bcftools +setGT -o  test/S10.setGT.vcf.gz -- -t q -i 'GT !~ "/" && FORMAT/AD[:1]/FORMAT/DP >=0.75'  -n 'c:1/1'
tabix -p vcf test/S10.setGT.vcf.gz 
bcftools consensus -p S10_ -f  S10.reference_masked.fa --mark-del "-" --mark-ins lc  --mark-snv lc -I test/S10.setGT.vcf.gz -o test/S10.consensus.fasta

## freebayers
bcftools norm --fasta-ref S10.reference_masked.fa -m- --check-ref w --output-type v S10.freebayes.vcf.gz --output test/S10.biallellic.vcf
bcftools view -e 'FORMAT/AD[:1]/FORMAT/DP < 0.25'  test/S10.biallellic.vcf > test/S10.biallellic-filtered.vcf
bcftools norm --fasta-ref S10.reference_masked.fa -m+ --check-ref w --output-type v test/S10.biallellic-filtered.vcf --output test/S10.multiallellic-filtered.vcf
 bcftools +setGT test/S10.multiallellic-filtered.vcf -- -t q -i 'strlen(ref) !=strlen(alt) && FORMAT/AD[:1]/FORMAT/DP >= 0.5' -n 'c:1/1' | bcftools +setGT  -- -t q -i 'strlen(ref) !=strlen(alt) && FORMAT/AD[:1]/FORMAT/DP < 0.5' -n 'c:0/0'  | bcftools +setGT -- -t q -i 'GT !~ "/"  && FORMAT/AD[:1]/FORMAT/DP < 0.75' -n 'c:0/1' | bcftools +setGT -o  test/S10.setGT.vcf.gz -- -t q -i 'GT !~ "/" && FORMAT/AD[:1]/FORMAT/DP >=0.75'  -n 'c:1/1'
tabix -p vcf test/S10.setGT.vcf.gz 
bcftools consensus -p S10_ -f  S10.reference_masked.fa --mark-del "-" --mark-ins lc  --mark-snv lc -I test/S10.setGT.vcf.gz -o test/S10.consensus.fasta

*/