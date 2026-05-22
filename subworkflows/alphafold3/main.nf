nextflow.enable.types = true

include { MSA as HMMER  } from '../../modules/alphafold3'
include { MSA as MMSEQS } from '../../modules/colabfold'
include { INFERENCE     } from '../../modules/alphafold3'

workflow ALPHAFOLD3 {

    take:
        input: Channel<Record>

    main:

        jobdef = 
            input
            .map { it -> tuple(it.model, it.input) }
            .groupBy()
            .map { model, query -> 
                record(
                    model: model,
                    query: query.toSet()
                )
            }

        msa =
            ( params.ALPHAFOLD3.USE_MMSEQS? MMSEQS( jobdef ) : HMMER( jobdef ) )
            .flatMap { rec ->
                rec.json.collect { it ->
                    record(
                        id: it.simpleName.minus('_data'),
                        json: it,
                        model: rec.model
                    )
                }
            }

        prediction = INFERENCE( msa )

    emit:
        msa: Channel<Record> = msa
        prediction: Channel<Record> = prediction
        jobcount: Value<Integer> = input.collect().map { it -> it.size() }
}