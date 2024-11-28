if ( params.ALPHAFOLD2.MODEL_PRESET == 'multimer' ) {
    include { MULTIMER as INFERENCE } from '../../modules/alphafold2'
    databases = ['uniref90', 'mgnify', 'bfd', 'uniprot']
} else {
    include { MONOMER as INFERENCE } from '../../modules/alphafold2'
    databases = ['uniref90', 'mgnify', 'bfd']
}  

include { AF2_METRICS as METRICS } from '../../modules/pymol'
include { MSA } from '../../modules/alphafold2'

workflow ALPHAFOLD2 {

    take:
        input

    main:
        input
            .map { fasta -> [ fasta, fasta ] }
            .splitFasta ( record: [ id: true ] )
            .groupTuple ( by: ( params.ALPHAFOLD2.MODEL_PRESET == 'multimer' ? 1 : [ 0, 1 ] ) )
            .map { record, fasta ->
                params.ALPHAFOLD2.MODEL_PRESET == 'multimer'
                ? [ [ ('A'..'H'), record.id ].transpose().collectEntries(), fasta ]
                : [ [ 'A': record.id ], fasta ]
            }
            .unique { meta, fasta -> meta }
            .set { fasta }

        MSA(
            fasta.splitFasta ( record: [ id: true, seqString: true ] ).filter { meta, record -> ( record.id in meta*.value ) }.unique { meta, record -> record },
            databases
        )

        fasta
            .combine ( MSA.out.msa )
            .filter { meta, fasta, record, msa -> ( record in meta*.value ) }
            .map { meta, fasta, record, msa -> [ groupKey( meta, meta*.value.unique().size() * databases.size() ), fasta, msa ] }
            .groupTuple( by: 0 )
            .map { meta, fasta, msa ->
                [ meta.getGroupTarget(), fasta.first() ] + ( params.ALPHAFOLD2.MODEL_PRESET == 'multimer' ? ('A'..'H').collect { chain -> msa.findAll { it.parent.name == meta[chain] } } : [ msa.unique() ] )
            }
            .set { msa }

        INFERENCE(msa) | METRICS

    emit:
        metrics    = METRICS.out.metrics
        prediction = INFERENCE.out.prediction
        jobcount   = input.count()
}