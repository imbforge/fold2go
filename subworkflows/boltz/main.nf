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
        
        INFERENCE(
            MSA.out.yaml.join(MSA.out.msa, by: 0)
        )

    emit:
        msa: Channel<Tuple<Map, Path>> = MSA.out.yaml.mix(MSA.out.msa)
        prediction: Channel<Tuple<Map, Path>> = INFERENCE.out.prediction
        jobcount: Channel<Integer> = input.count()
}