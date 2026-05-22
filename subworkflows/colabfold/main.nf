nextflow.enable.types = true

include { MSA; INFERENCE } from '../../modules/colabfold'

workflow COLABFOLD {

    take:
        input: Channel<Record>

    main:

        msa =
            MSA(
                input
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
        jobcount: Value<Integer> = input.map { it -> it.input.countFasta() }.collect().map { it -> it.sum() }
}