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


    
    emit:
        ch_flu_primers
        ch_flu_db_msh 
        ch_flu_db_fasta 
        ch_typing_db  
        versions = ch_versions
        
}
