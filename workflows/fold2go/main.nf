switch ( params.MODEL_PRESET ) {
    case { it.startsWith('monomer') }:
        include { MONOMER as ALPHAFOLD } from '../../modules/alphafold'
        databases = ['uniref90', 'mgnify', 'bfd']
        chains    = ['A']
        break
    case "multimer":
        include { MULTIMER as ALPHAFOLD } from '../../modules/alphafold'
        databases = ['uniref90', 'mgnify', 'bfd', 'uniprot']
        chains    = ('A'..'H')
        break
}

include { MSA   } from '../../modules/alphafold'
include { PYMOL } from '../../modules/pymol'
include { SHINY } from '../../modules/shiny'

workflow FOLD2GO {

    Channel
        .fromPath( params.IN )
        .map { fasta -> [ fasta, fasta ] }
        .splitFasta ( record: [ id: true ] )
        .groupTuple ( by: 1 )
        .map { record, fasta -> [ [ chains, record.id ].transpose().collectEntries(), fasta ] }
        .set { fasta }

    SHINY(
        fasta.count(),
        params.OUT,
        workflow.runName,
        workflow.launchDir
    )

    MSA(
        fasta.splitFasta ( record: [ id: true, seqString: true ] ).filter { meta, record -> ( record.id in meta*.value ) }.unique { meta, record -> record },
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
            sendMail{
                to "${params.EMAIL}"
                from "alphafold@imb-mainz.de"
                subject "AlphaFold (${workflow.runName})"
                attach "${params.OUT}/${workflow.runName}/template_indep_info.tsv"

                """
                Dear ${workflow.userName},

                AlphaFold predictions are complete, please find a table with useful metrics attached.
                Results of this run have all been stored at ${params.OUT}/${workflow.runName}.

                ---
                Deet-doot-dot, I am a bot.
                """.stripIndent()
            }
        }
}
