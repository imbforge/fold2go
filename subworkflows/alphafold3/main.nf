include { MSA; MODEL } from '../../modules/alphafold3'

workflow ALPHAFOLD3 {

    take:
        input

    main:
        input
            .map { 
                json -> [ [ id: new groovy.json.JsonSlurper().parse(json).name ], json ]
            }
            .set { json }

        MSA(json) // | INFERENCE

    emit:
        prediction = Channel.empty()
}