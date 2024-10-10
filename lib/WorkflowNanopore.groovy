//
// This file holds several functions specific to the workflow/influenza.nf in the nf-core/influenza pipeline
//

import groovy.text.SimpleTemplateEngine

class WorkflowNanopore {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, valid_params) {
        
        if (!valid_params['nanopore_reads_mapping_tools'].contains(params.mapping_tool)) {
            log.error "Invalid option: '${params.mapping_tool}'. Valid options for nanopore data is '--mapping_tool': ${valid_params['nanopore_reads_mapping_tools'].join(', ')}."
            System.exit(1)
        }
        if (!valid_params['nanopore_variant_callers'].contains(params.variant_caller)) {
            log.error "Invalid option: '${params.variant_caller} for nanopore data'. Valid options for '--variant_caller': ${valid_params['nanopore_variant_callers'].join(', ')}."
            System.exit(1)
        }
        
    }

}
