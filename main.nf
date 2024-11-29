#!/usr/bin/env nextflow

include { validateParameters } from 'plugin/nf-schema'

validateParameters()

include { FOLD2GO } from './workflows/fold2go'

workflow {
    
    FOLD2GO()

}
