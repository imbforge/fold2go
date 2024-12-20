include { MSA; INFERENCE } from '../../modules/alphafold3'

workflow ALPHAFOLD3 {

    take:
        input

    main:

        input.map { json ->
            def jobdef = new groovy.json.JsonSlurper().parse(json)
            [ jobdef.name, ( jobdef instanceof List ? 'alphafoldserver' : jobdef.dialect ), json ]
        }
        .groupTuple( by: params.ALPHAFOLD3.GROUP_MSA ? 1 : [0, 1] )
        .map { id, dialect, json ->
            [ [ id: id, jobsize: ( dialect == 'alphafold3' ? json.size() : id.flatten().size() ) ], json ]
        }
        .set { jobdef }

        MSA( jobdef )

        INFERENCE(
            MSA.out.json.flatten().map { [ [ id: it.name.minus('_data.json'), model: 'alphafold3' ], it ] }
        )

    emit:
        prediction = INFERENCE.out.prediction
        jobcount   = jobdef.sum { meta, json -> meta.jobsize }
}
