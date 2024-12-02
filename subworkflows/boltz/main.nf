include { MSA; INFERENCE; METRICS } from '../../modules/boltz'

workflow BOLTZ {

    take:
        input

    main:

        input.map { yaml ->
            [ [id: yaml.simpleName], yaml ]
        }
        | MSA
        | INFERENCE
        | METRICS

    emit:
        metrics    = METRICS.out.metrics
        prediction = INFERENCE.out.prediction
        jobcount   = input.count()
}