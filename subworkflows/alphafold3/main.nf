include { MSA; INFERENCE } from '../../modules/alphafold3'

workflow ALPHAFOLD3 {

    take:
        input

    main:
        input
            .map { json ->
                def job = new groovy.json.JsonSlurper().parse(json)
                [
                    [
                        id: job.name,
                        chains: job.sequences.collect {
                            [
                                realm: it.entrySet().key.pop(),
                                sequence: it.entrySet().value.sequence.pop(),
                                hash: it.entrySet().value.sequence.pop().md5()
                            ]
                        }
                    ],
                    json
                ]
            }
            .tap { jobdef }
            .map { meta, json ->
                meta.chains
                    .findAll { it.realm in ['protein', 'rna'] }
                    .collect {
                        [
                            name: it.hash,
                            sequences:
                            [
                                [
                                    (it.realm): [ id: 'XYZ', sequence: it.sequence ]
                                ]
                            ],
                            modelSeeds: [1],
                            dialect: 'alphafold3',
                            version: 1
                        ]
                    }   
            }
            .flatten()
            .unique()
            .collectFile {
                [ "${it.name}.json", groovy.json.JsonOutput.toJson(it) ]
            }
            .set { json }

    MSA(
        json.map { json -> [ [id: json.baseName], json] }
    )

    jobdef
        .combine(MSA.out.json)
        .filter { meta, jobdef, hash, msa -> ( hash.id in meta.chains*.hash ) }
        .map { meta, jobdef, hash, msa -> [ groupKey( meta.id, meta.chains.size() ), jobdef, msa ] }
        .groupTuple(by: [0,1])
        .set { job }

    INFERENCE(job)

    emit:
        prediction = Channel.empty()
}