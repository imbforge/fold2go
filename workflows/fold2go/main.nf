include { ALPHAFOLD2 } from '../../subworkflows/alphafold2'
include { ALPHAFOLD3 } from '../../subworkflows/alphafold3'
include { BOLTZ      } from '../../subworkflows/boltz'
include { SHINY      } from '../../modules/shiny'
include { METRICS    } from '../../modules/metrics'

workflow FOLD2GO {

    input =
        channel
        .fromPath( params.IN )
        .branch { fname ->
            fasta: fname =~ /.(fasta|fa)$/
            json : fname =~ /.json$/
            yaml : fname =~ /.(yaml|yml)$/
        }

    input.fasta | ALPHAFOLD2
    input.json  | ALPHAFOLD3
    input.yaml  | BOLTZ

    jobcount =
        ALPHAFOLD2.out.jobcount
        .mix(ALPHAFOLD3.out.jobcount)
        .mix(BOLTZ.out.jobcount)
        .sum()
        .collectFile { njobs ->
            [
            "shiny_config.json",
            """
            {
                "njobs": ${njobs},
                "data": "${params.OUT}/${workflow.runName}",
                "log": "${workflow.launchDir}/.nextflow.log"
            }
            """
            ]
        }

    SHINY(
        params.SOCKET ?: "${workflow.workDir}/shiny.sock",
        jobcount
    )

    METRICS(
        ALPHAFOLD2.out.prediction.mix(ALPHAFOLD3.out.prediction).mix(BOLTZ.out.prediction)
    )
    
    METRICS.out.metrics
    .collectFile ( storeDir: "${params.OUT}/${workflow.runName}", keepHeader: true ) {
        model, metrics -> [ "${model}_metrics.tsv", metrics ]
    }
    .collect()
    .map { metrics ->
        if( params.EMAIL ) {
            try {
                sendMail (
                    to: "${params.EMAIL}",
                    subject: "fold2go (${workflow.runName})",
                    attach: metrics,
                    body: """
                    Dear ${workflow.userName},

                    fold2go predictions are complete, please find some useful metrics attached.
                    Results of this run have all been stored at ${params.OUT}/${workflow.runName}.

                    ---
                    Deet-doot-dot, I am a bot.
                    """.stripIndent()
                )
            }
            catch( Exception e ) {
                log.warn "Failed to send notification email to ${params.EMAIL}"
                log.warn e.message
            }
        }
    }
}
