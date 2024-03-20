//
// This file holds several functions specific to the workflow/influenza.nf in the nf-core/influenza pipeline
//

import groovy.text.SimpleTemplateEngine

class WorkflowNanopore {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, valid_params) {
        
        if (!valid_params['nanopore_reads_mapping_tools'].contains(params.nanopore_reads_mapping_tool)) {
            log.error "Invalid option: '${params.nanopore_reads_mapping_tool}'. Valid options for '--nanopore_reads_mapping_tool': ${valid_params['nanopore_reads_mapping_tools'].join(', ')}."
            System.exit(1)
        }
        
    }

}
