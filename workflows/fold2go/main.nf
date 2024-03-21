switch ( params.MODEL_PRESET ) {
    case { it.startsWith('monomer') }:
        include { MONOMER as ALPHAFOLD } from '../../modules/alphafold'
        databases = ['uniref90', 'mgnify', 'bfd']
        break
    case "multimer":
        include { MULTIMER as ALPHAFOLD } from '../../modules/alphafold'
        databases = ['uniref90', 'mgnify', 'bfd', 'uniprot']
        break
}

include { MSA   } from '../../modules/alphafold'
include { PYMOL } from '../../modules/pymol'
include { COMMS } from '../../modules/comms'

workflow FOLD2GO {

    Channel
        .fromPath( params.IN )
        .map { fasta -> [ fasta, fasta ] }
        .splitFasta ( record: [ id: true ], elem: 1 )
        .groupTuple ( by: 0 )
        .map { fasta, record -> [ [ ('A'..'H'), record.id ].transpose().collectEntries(), fasta ] }
        .set { fasta }

        COMMS(
            fasta.count(),
            params.OUT,
            workflow.runName,
            workflow.launchDir
        )

        MSA(
            fasta.splitFasta ( record: [ id: true, seqString: true ] ).unique { fasta, record -> record },
            databases
        )

        fasta
            .combine ( MSA.out.msa )
            .filter { meta, fasta, record, msa -> ( record in meta*.value ) }
            .map { meta, fasta, record, msa -> [ groupKey( meta, meta*.value.unique().size() * databases.size() ), fasta, msa ] }
            .groupTuple()
            .map { meta, fasta, msa -> [ meta.getGroupTarget(), fasta.unique() ] + ( params.MODEL_PRESET == 'multimer' ? ('A'..'H').collect { chain -> msa.findAll { it.parent.name == chain } } : [ msa ] ) }
            .set { msa }

        ALPHAFOLD(msa) | PYMOL

        PYMOL.out.metrics
            .collectFile(name: "template_indep_info.tsv", storeDir: "${params.OUT}/${workflow.runName}", keepHeader: true)
            .subscribe onComplete: {
                sendMail(
                    from: "alphafold@imb-mainz.de",
                    to: "${params.EMAIL}",
                    subject: "AlphaFold (${workflow.runName})",
                    text: "Predictions are complete!",
                    attach: "${params.OUT}/${workflow.runName}/template_indep_info.tsv"
                )
            }
}
