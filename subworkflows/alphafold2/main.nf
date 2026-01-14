include { MSA; INFERENCE_MONOMER; INFERENCE_MULTIMER } from '../../modules/alphafold2'

workflow ALPHAFOLD2 {

    take:
        input: Channel<Path>

    main:
        Boolean multimer = ( params.ALPHAFOLD2.MODEL_PRESET == 'multimer' )
        List databases = multimer ? ['uniref90', 'mgnify', 'bfd', 'uniprot'] : ['uniref90', 'mgnify', 'bfd']

        sequences = 
            input
            .map { fasta -> [ fasta, fasta ] }
            .splitFasta ( record: [ id: true ] )
            .groupTuple ( by: ( multimer ? 1 : [ 0, 1 ] ) )
            .map { record, fasta ->
                multimer ? [ [ ('A'..'H'), record.id ].transpose().collectEntries(), fasta ] : [ [ 'A': record.id ], fasta ]
            }
            .unique { meta, _fasta -> meta }

        jobdef =
            sequences
            .splitFasta ( record: [ id: true, seqString: true ] )
            .filter { meta, record -> ( record.id in meta*.value ) }
            .unique { _meta, record -> record }
            .combine ( databases )

        MSA(jobdef)

        chains = 
            sequences
            .combine ( MSA.out.msa )
            .filter { meta, _fasta, record, _msa -> ( record in meta*.value ) }
            .map { meta, fasta, _record, msa -> [ groupKey( meta, meta*.value.unique().size() * databases.size() ), fasta, msa ] }
            .groupTuple( by: 0 )
            .map { meta, fasta, msa ->
                [ [ id: meta.getGroupTarget()*.value.join('.'), model: "alphafold2_${params.ALPHAFOLD2.MODEL_PRESET}" ], fasta.first() ] + ( multimer ? ('A'..'H').collect { chain -> msa.sort().findAll { it -> it.parent.name == meta[chain] } } : [ msa.unique().sort() ] )
            }
        
        multimer
            ? INFERENCE_MULTIMER(chains)
            : INFERENCE_MONOMER(chains)

    emit:
        prediction: Channel<Tuple<Map, Set<Path>>> = ( multimer ? INFERENCE_MULTIMER : INFERENCE_MONOMER ).out.prediction
        jobcount: Channel<Integer> = input.count()
}