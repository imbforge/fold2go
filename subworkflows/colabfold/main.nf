include { COLABFOLD_SEARCH; COLABFOLD_BATCH } from '../../modules/colabfold'

workflow COLABFOLD {

    take:
        input1: Channel<Path>
        input2: Channel<Path>

    main:

        COLABFOLD_SEARCH(input1, input2)
        
        msa = COLABFOLD_SEARCH.out.msa.flatMap().map { it -> tuple([ id: it.simpleName, model: 'colabfold' ], it) }.view()

        COLABFOLD_BATCH(
            msa
        )

    emit:
        msa: Channel<Tuple<Map, Path>> = msa
        prediction: Channel<Tuple<Map, Path>> = COLABFOLD_BATCH.out.prediction
        jobcount: Value<Integer> = msa.count()
}