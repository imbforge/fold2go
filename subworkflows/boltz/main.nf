include { MSA; INFERENCE } from '../../modules/boltz'

workflow BOLTZ {

    take:
        input

    main:

        input.map { yaml ->
            [ [ id: yaml.simpleName, model: "${params.BOLTZ.MODEL_PRESET}" ], yaml ]
        }
        | MSA
        | INFERENCE

    emit:
        prediction = INFERENCE.out.prediction
        jobcount   = input.count()
}