nextflow.enable.types = true

include { MSA; INFERENCE } from '../../modules/colabfold'

workflow COLABFOLD {

    take:
        input: Channel<Record>

    main:

        jobdef =
            input
            .map { it -> tuple(it.model, it.input) }
            .groupBy()
            .map { model, fasta -> 
                record(
                    model: model,
                    query: fasta.toSet()
                )
            }

        msa =
            MSA(
                jobdef
            )
            .flatMap { rec ->
                rec.msa.collect { it ->
                    record(id: it.simpleName, msa: it, model: rec.model)
                }
            }

        prediction = INFERENCE(msa)

    emit:
        msa: Channel<Record> = msa
        prediction: Channel<Record> = prediction
        jobcount: Value<Integer> = msa.collect().map { it -> it.size() }
}