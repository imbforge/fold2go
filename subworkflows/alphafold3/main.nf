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
            ( params.ALPHAFOLD3.MSA_METHOD == 'mmseqs2' ? MMSEQS( input ) : HMMER( jobdef ) )
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
        jobcount: Value<Integer> = ( params.ALPHAFOLD3.MSA_METHOD == 'mmseqs2' ? input.map { it -> it.input.countFasta() }.collect().map { it -> it.max() ?: 0 } : input.collect().map { it -> it.size() } )
}
