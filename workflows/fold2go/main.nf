switch ( params.MODEL_PRESET ) {
    case { it ==~ /^(mono|multi)mer.*/ }:
        include { ALPHAFOLD2 as ALPHAFOLD } from '../../subworkflows/alphafold2'
        break
    default:
        include { ALPHAFOLD3 as ALPHAFOLD } from '../../subworkflows/alphafold3'
}

include { SHINY } from '../../modules/shiny'

workflow FOLD2GO {

    Channel
        .fromPath( params.IN )
        .tap { input }
        .count()
        .set { count }
    
    ALPHAFOLD(input)

    SHINY(
        count,
        params.OUT,
        workflow.runName,
        workflow.launchDir
    )

    ALPHAFOLD.out.metrics
        .collectFile(name: "template_indep_info.tsv", storeDir: "${params.OUT}/${workflow.runName}", keepHeader: true)
        .subscribe onNext: {
            sendMail {
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
