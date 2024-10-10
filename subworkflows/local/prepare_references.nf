//
// prepare databases
//


workflow PREPARE_REFERENCES {
    main:

    ch_versions = Channel.empty()

    //
    // Prepare fasta flu primers 
    //
    ch_flu_primers = Channel.empty()
    
    if (params.flu_primers) {
        ch_flu_primers = Channel.value(file(params.flu_primers))
    }
    else {
        log.error "Please specify a valid fasta format flu primer file"
        System.exit(1)
    }
     
     //
    // hostile reference
    //
    ch_hostile_ref_minimap2 = Channel.empty()
    if (params.hostile_human_ref_minimap2) {
        ch_hostile_ref_minimap2 = Channel.value(file(params.hostile_human_ref_minimap2))
    }
    else {
        log.error "Please specify a valid hostime human reference file for minimap2"
        System.exit(1)
    }

    ch_hostile_ref_bowtie2 = Channel.empty()
    if (params.hostile_human_ref_bowtie2) {
        ch_hostile_ref_bowtie2 = Channel.value(file(params.hostile_human_ref_bowtie2))
    }
    else {
        log.error "Please specify a valid hostime human reference file from bowtie2"
        System.exit(1)
    }
    //
    // Prepare flu mash
    //
    ch_flu_db_msh = Channel.empty()
    
    if (params.flu_db_msh) {
        ch_flu_db_msh = Channel.value(file(params.flu_db_msh))
    }
    else {
        log.error "Please specify a valid  Flu database mash file"
        System.exit(1)
    }
    

    //
    // Prepare flu fasta file
    //
    ch_flu_db_fasta = Channel.empty()

    if (params.flu_db_fasta) {
        ch_flu_db_fasta = Channel.value(file(params.flu_db_fasta))
    }
    else {
        log.error "Please specify a valid  flu database fasta file"
        System.exit(1)
    }

    //
    // Prepare snpeff.config file
    //
    ch_snpeff_config = Channel.empty()
    
    if (params.snpeff_config) {
        ch_snpeff_config = Channel.value(file(params.snpeff_config))
    }
    else {
        log.error "Please specify a valid snpeff.config file"
        System.exit(1)
    }

    //
    // Prepare snpeff dataDir
    //
    ch_snpeff_dataDir = Channel.empty()
    
    if (params.snpeff_dataDir) {
        ch_snpeff_dataDir = Channel.value(file(params.snpeff_dataDir))
    }
    else {
        log.error "Please specify a valid snpeff dataDir directory"
        System.exit(1)
    }

    //
    // Prepare typing database
    //
    ch_typing_db = Channel.empty()

    if (params.typing_db) {
        ch_typing_db = Channel.value(file(params.typing_db))
    }
    else {
        log.error "Please specify a valid typing database"
        System.exit(1)
    }

    //clair3 variant_model
    ch_clair3_variant_model = Channel.empty()
    if (params.clair3_variant_model) {
        ch_clair3_variant_model = Channel.value(file(params.clair3_variant_model))
    }
    else {
        log.error "Please specify a valid clair3 variant model path"
        System.exit(1)
    }
    
    emit:
        ch_flu_primers
        ch_flu_db_msh 
        ch_flu_db_fasta
        ch_snpeff_config
        ch_snpeff_dataDir
        ch_typing_db
        ch_clair3_variant_model
        ch_hostile_ref_bowtie2
        ch_hostile_ref_minimap2
        versions = ch_versions
        
}
