include { ALPHAFOLD2 } from '../../subworkflows/alphafold2'
include { ALPHAFOLD3 } from '../../subworkflows/alphafold3'
include { BOLTZ      } from '../../subworkflows/boltz'
include { SHINY      } from '../../modules/shiny'

workflow FOLD2GO {

    Channel
        .fromPath( params.IN )
        .branch {
            fasta: it =~ /.(fasta|fa)$/
            json : it =~ /.json$/
            yaml : it =~ /.(yaml|yml)$/
        }
        .set { input }

    input.fasta | ALPHAFOLD2
    input.json  | ALPHAFOLD3
    input.yaml  | BOLTZ

    SHINY(
        params.PORT ?: new Random().nextInt(39000, 39200),
        ALPHAFOLD2.out.jobcount.mix(ALPHAFOLD2.out.jobcount).mix(BOLTZ.out.jobcount).sum(),
        params.OUT,
        workflow.runName,
        workflow.launchDir
    )

    ALPHAFOLD2.out.metrics.mix(ALPHAFOLD3.out.metrics).mix(BOLTZ.out.metrics)
        .collectFile ( storeDir: "${params.OUT}/${workflow.runName}", keepHeader: true ) {
            preset, metrics -> [ "${preset}_metrics.tsv", metrics ]
        }
        .collect()
        .map { metrics ->
            if( params.EMAIL ) {
                try {
                    sendMail {
                        to "${params.EMAIL}"
                        from "alphafold@imb-mainz.de"
                        subject "AlphaFold (${workflow.runName})"
                        attach metrics

                        """
                        Dear ${workflow.userName},

                        AlphaFold predictions are complete, please find some useful metrics attached.
                        Results of this run have all been stored at ${params.OUT}/${workflow.runName}.

                        ---
                        Deet-doot-dot, I am a bot.
                        """.stripIndent()
                    }
                }
                catch( Exception e ) {
                    log.warn "Failed to send notification email to ${params.EMAIL}"
                    log.error e
                }
        }
        }
}
