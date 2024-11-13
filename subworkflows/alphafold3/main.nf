include { MSA; MODEL } from '../../modules/alphafold3'

workflow ALPHAFOLD3 {

    main:
        Channel
            .fromPath( params.IN )
            .map { 
                json -> [ [ id: new groovy.json.JsonSlurper().parse(json).name ], json ]
            }
            .set { json }

        MSA(json) // | MODEL

    emit:
        prediction = Channel.empty()
        count      = json.count()
}