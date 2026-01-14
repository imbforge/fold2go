include { MSA; INFERENCE } from '../../modules/boltz'

workflow BOLTZ {

    take:
        input: Channel<Path>

    main:

        input
        .map { yaml ->
            [ [ id: yaml.simpleName, model: "${params.BOLTZ.MODEL_PRESET}" ], yaml ]
        }
        | MSA
        | INFERENCE

    emit:
        prediction: Channel<Tuple<Map, Set<Path>>> = INFERENCE.out.prediction
        jobcount: Channel<Integer> = input.count()
}