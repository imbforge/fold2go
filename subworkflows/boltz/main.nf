include { MSA; INFERENCE } from '../../modules/boltz'

workflow BOLTZ {

    take:
        input: Channel<Path>

    main:

        jobdef =
            input
            .map { yaml ->
                [ [ id: yaml.simpleName, model: "${params.BOLTZ.MODEL_PRESET}" ], yaml ]
            }

        MSA(jobdef)

        INFERENCE(
            MSA.out.msa,
            workDir
        )

    emit:
        prediction: Channel<Tuple<Map, Path>> = INFERENCE.out.prediction
        jobcount: Channel<Integer> = input.count()
}