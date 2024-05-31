//
// Check input samplesheet and get read channels
//

include { SAMPLESHEETCHECK } from '../../modules/local/samplesheetcheck'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEETCHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_read_channels(it) }
        .set { reads }

    //reads.view()

    reads.map { 
        meta, reads, long_fastq -> [ meta, reads ] }//.view()
        .filter{ meta, reads -> !reads[0].equals('NA') && !reads[1].equals('NA') }
        .set { shortreads }
    
    //shortreads.view() 
    

    reads.map {
        meta, reads, long_fastq -> [ meta, long_fastq ] }
        .filter{ meta, long_fastq -> !long_fastq.equals('NA') }
        .set { longreads }
       
    //longreads.view()
    
    
    reads
        .map { meta, reads, long_fastq -> meta.id }
        .set {ids}

    //ids.view()

    emit:
    reads      // channel: [ val(meta), [ reads ], long_fastq ]
    shortreads // channel: [ val(meta), [ reads ] ]
    longreads  // channel: [ val(meta), long_fastq ]
    ids
    versions = SAMPLESHEETCHECK.out.versions // channel: [ versions.yml ]
}
// Function to get list of [ meta, [ fastq_1, fastq_2 ], long_fastq ]
def create_read_channels(LinkedHashMap row) {
    
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = !(row.fastq_1 == 'NA') && !(row.fastq_2 == 'NA') ? false : true
    
    def array = []
    // check short reads
    if ( !(row.fastq_1 == 'NA') ) {
        if ( !file(row.fastq_1).exists() ) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
        }
        
        if(file(row.fastq_1).size() == 0){
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file is empty!\n${row.fastq_1}" 
        }
        
        fastq_1 = file(row.fastq_1)
    } else { fastq_1 = 'NA' }
    if ( !(row.fastq_2 == 'NA') ) {
        if ( !file(row.fastq_2).exists() ) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        
        if(file(row.fastq_2).size() == 0){
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file is empty!\n${row.fastq_2}" 
        }
        fastq_2 = file(row.fastq_2)
    } else { fastq_2 = 'NA' }

    // check long_fastq
    if ( !(row.long_fastq == 'NA') ) {
        if ( !file(row.long_fastq).exists() ) {
            exit 1, "ERROR: Please check input samplesheet -> Long FastQ file does not exist!\n${row.long_fastq}"
        }
        
        if(file(row.long_fastq).size() == 0){
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file is empty!\n${row.long_fastq}" 
        }
        long_fastq = file(row.long_fastq)
    } else { long_fastq = 'NA' }

   
    // prepare output // currently does not allow single end data!
    if ( meta.single_end ) {
        array = [ meta, [fastq_1, fastq_2] , long_fastq]
    } else {
        array = [ meta, [ fastq_1, fastq_2 ], long_fastq]
    } 
    return array 
}
