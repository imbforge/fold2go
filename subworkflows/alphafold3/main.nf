include { MSA; INFERENCE } from '../../modules/alphafold3'

workflow ALPHAFOLD3 {

    take:
        input: Channel<Path>

    main:

        jobdef =
            input
            .map { json ->
                def jobdef = new groovy.json.JsonSlurper().parse(json)
                [ jobdef.name, ( jobdef instanceof List ? 'alphafoldserver' : jobdef.dialect ), json ]
            }
            .groupTuple( by: params.ALPHAFOLD3.GROUP_MSA ? 1 : [0, 1] )
            .map { id, dialect, json ->
                [ [ id: id, jobsize: ( dialect == 'alphafold3' ? json.size() : id.flatten().size() ) ], json as Set<Path> ]
            }

        MSA( jobdef )

        INFERENCE(
            MSA.out.msa.transpose().map { _meta, json -> [ [ id: json.name.minus('_data.json'), model: 'alphafold3' ], json ] }
        )

    emit:
        msa: Channel<Tuple<Map, Path>> = MSA.out.msa.transpose().map { _meta, json -> [ [ id: json.name.minus('_data.json').toUpperCase(), model: 'alphafold3' ], json ] }
        prediction: Channel<Tuple<Map, Path>> = INFERENCE.out.prediction
        jobcount: Channel<Integer> = jobdef.sum { meta, _json -> meta.jobsize }
}
