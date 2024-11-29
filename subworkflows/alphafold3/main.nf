include { MSA; INFERENCE         } from '../../modules/alphafold3'
include { AF3_METRICS as METRICS } from '../../modules/pymol'

workflow ALPHAFOLD3 {

    take:
        input

    main:
        // instantiate json slurper
        def slurper = new groovy.json.JsonSlurper()

        input.map { json ->
            def jobdef = slurper.parse(json)
            [ jobdef.name, ( jobdef instanceof List ? 'alphafoldserver' : jobdef.dialect ), json ]
        }
        .groupTuple( by: params.ALPHAFOLD3.GROUP_MSA ? 1 : [0, 1] )
        .map { id, dialect, json ->
            [ [ id: id, jobsize: ( dialect == 'alphafold3' ? json.size() : id.flatten().size() ) ], json ]
        }
        .set { jobdef }

        MSA( jobdef )

        INFERENCE(
            MSA.out.json.flatten().map { [ [ id: slurper.parse(it).name ], it ] }
        )
        
        METRICS(
            INFERENCE.out.prediction
        )

    emit:
        metrics    = METRICS.out.metrics
        prediction = INFERENCE.out.prediction
        jobcount   = jobdef.sum { meta, json -> meta.jobsize }
}
