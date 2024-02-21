#!/usr/bin/env nextflow

include { FOLD2GO } from './workflows/fold2go'

workflow {
    
    // dump all parameters to json
    Channel.of( groovy.json.JsonOutput.toJson(params) ).collectFile(cache: true, storeDir:"${params.OUT}/logs", name:"${workflow.runName}.json")
  
    // spawn pipeline run
    FOLD2GO()

}
