nextflow.enable.types = true

include { MSA; INFERENCE } from '../../modules/alphafold3'

workflow ALPHAFOLD3 {

    take:
        input: Channel<Record>

    main:

        jobdef = 
            input
            .map { it -> tuple(it.model, it.input) }
            .groupBy()
            .map { model, json -> 
                record(
                    model: model,
                    queries: json.toSet()
                )
            }

        msa = 
            MSA(
                jobdef
            )
            .flatMap { rec ->
                rec.msa.collect { it ->
                    // extract id from filename in a convoluted way (.minus() is not something the language server accepts)
                    def id = it.name.tokenize('_data.json').first()
                    record(id: id, msa: it, model: rec.model)
                }
            }
        
        prediction = INFERENCE( msa )

    emit:
        msa: Channel<Record> = msa
        prediction: Channel<Record> = prediction
        jobcount: Value<Integer> = input.collect().map { it -> it.size() }
}