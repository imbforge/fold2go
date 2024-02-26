#!/usr/bin/env nextflow

include { FOLD2GO } from './workflows/fold2go'

workflow {
    
    FOLD2GO()

}
