include { MSA; INFERENCE } from '../../modules/colabfold'

workflow COLABFOLD {

    take:
        input: Channel<Path>

    main:

        MSA(
            input.collect()
        )
        
        msa = MSA.out.msa.flatMap().map { it -> tuple([ id: it.simpleName, model: 'colabfold' ], it) }

        INFERENCE(
            msa
        )

    emit:
        msa: Channel<Tuple<Map, Path>> = msa
        prediction: Channel<Tuple<Map, Path>> = INFERENCE.out.prediction
        jobcount: Value<Integer> = msa.count()
}