include { MSA; INFERENCE_MONOMER; INFERENCE_MULTIMER } from '../../modules/alphafold2'

workflow ALPHAFOLD2 {

    take:
        input

    main:

        def ( Boolean multimer, List databases) = [ ( params.ALPHAFOLD2.MODEL_PRESET == 'multimer' ), ['uniref90', 'mgnify', 'bfd', 'uniprot'] ]

        input
            .map { fasta -> [ fasta, fasta ] }
            .splitFasta ( record: [ id: true ] )
            .groupTuple ( by: ( multimer ? 1 : [ 0, 1 ] ) )
            .map { record, fasta ->
                multimer ? [ [ ('A'..'H'), record.id ].transpose().collectEntries(), fasta ] : [ [ 'A': record.id ], fasta ]
            }
            .unique { meta, fasta -> meta }
            .set { fasta }

        MSA(
            fasta.splitFasta ( record: [ id: true, seqString: true ] ).filter { meta, record -> ( record.id in meta*.value ) }.unique { meta, record -> record },
            multimer ? databases : databases.minus('uniprot')
        )

        fasta
            .combine ( MSA.out.msa )
            .filter { meta, fasta, record, msa -> ( record in meta*.value ) }
            .map { meta, fasta, record, msa -> [ groupKey( meta, meta*.value.unique().size() * databases.size() ), fasta, msa ] }
            .groupTuple( by: 0 )
            .map { meta, fasta, msa ->
                [ [ id: meta.getGroupTarget()*.value.join('.'), model: "alphafold2_${params.ALPHAFOLD2.MODEL_PRESET}" ], fasta.first() ] + ( multimer ? ('A'..'H').collect { chain -> msa.findAll { it.parent.name == meta[chain] } } : [ msa.unique() ] )
            } | ( multimer ? INFERENCE_MULTIMER : INFERENCE_MONOMER )

    emit:
        prediction = ( multimer ? INFERENCE_MULTIMER : INFERENCE_MONOMER ).out.prediction
        jobcount   = input.count()
}