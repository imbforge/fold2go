include { MSA; INFERENCE } from '../../modules/boltz'

workflow BOLTZ {

    take:
        input

    main:

        MSA(
            input
            .map { yaml ->
                [ [ id: yaml.simpleName, model: "${params.BOLTZ.MODEL_PRESET}" ], yaml ]
            }
        )

        INFERENCE(
            MSA.out.msa,
            workDir
        )

    emit:
        prediction = INFERENCE.out.prediction
        jobcount   = input.count()
}