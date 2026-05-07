nextflow.enable.types = true

include { MSA; INFERENCE } from '../../modules/boltz'

workflow BOLTZ {

    take:
        input: Channel<Record>

    main:

        msa = MSA( input )

        prediction = INFERENCE( msa )

    emit:
        msa: Channel<Record> = msa
        prediction: Channel<Record> = prediction
        jobcount: Value<Integer> = input.collect().map { it -> it.size() }
}