include { MSA; INFERENCE } from '../../modules/alphafold3'

workflow ALPHAFOLD3 {

    take:
        input

    main:
        input
            .map { json ->
                def job = new groovy.json.JsonSlurper().parse(json)
                [
                    [ id: job.name ],
                    job.sequences.collect {
                        // get protein and rna sequences w/o user-defined msa and calculate md5 hash over their sequence 
                        def (realm, entity) = [ it*.key.pop(), it*.value.pop() ]
                        ( realm in ['protein', 'rna'] && !('unpairedMsa' in entity*.key) ) ? [ type: realm, sequence: entity.sequence, hash: entity.sequence.md5() ] : [ type: realm ]
                    },
                    json
                ]
            }
            .set { jobdef }

    // create msa definition for each rna and protein entity to scatter msa jobs
    // use md5 hash to uniquely identify and deduplicate each msa job
    // this allows caching of individual msa across individual job definitions and pipeline runs
    // TODO: figure out msa pairing
    // TODO: consider compressing MSA output to save storage space
    jobdef
        .transpose()
        .filter { meta, entity, json -> entity.sequence }
        .map { meta, entity, json ->
            [
                name: entity.hash,
                sequences: [ [ (entity.type): [ id: 'XYZ', sequence: entity.sequence ] ] ],
                modelSeeds: [1],
                dialect: 'alphafold3',
                version: 1
            ]
        }
        .unique { it.name }
        .collectFile { [ "${it.name}.json", groovy.json.JsonOutput.toJson(it) ] }
        .map { json -> [ [ id: json.baseName ], json] }
        .set { msadef }

    MSA(msadef)

    // gather msa results and pair w/ initial job definition
    jobdef
        .combine(MSA.out.json)
        .filter { meta, entities, jobdef, hash, msa -> ( hash.id in entities*.hash ) }
        .map { meta, entities, jobdef, hash, msa -> [ groupKey( meta.id, entities.count { it.hash } ), jobdef, msa ] }
        .groupTuple( by: [0,1] )
        .set { job }

    INFERENCE(job)

    emit:
        prediction = Channel.empty()
}