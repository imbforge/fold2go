include { ALPHAFOLD2 } from '../../subworkflows/alphafold2'
include { ALPHAFOLD3 } from '../../subworkflows/alphafold3'
include { SHINY      } from '../../modules/shiny'

workflow FOLD2GO {

    Channel
        .fromPath( params.IN )
        .branch {
            af2: it =~ /.(fasta|fa)$/
            af3: it.isDirectory() || it =~ /.json$/
        }
        .set { input }
    
    ALPHAFOLD2(
        input.af2
    )

    ALPHAFOLD3(
        input.af3
    )

    SHINY(
        ALPHAFOLD3.out.jobcount.mix(ALPHAFOLD2.out.jobcount).sum(),
        params.OUT,
        workflow.runName,
        workflow.launchDir
    )

    ALPHAFOLD3.out.metrics
        .mix(ALPHAFOLD2.out.metrics)
        .collectFile ( storeDir: "${params.OUT}/${workflow.runName}", keepHeader: true ) {
            model, metrics -> [ "${model}_metrics.tsv", metrics ]
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
 
workflow FOLD2GOs {
    
    ALPHAFOLD(
        Channel.fromPath( params.IN )
    )

    SHINY(
        ALPHAFOLD.out.jobcount,
        params.OUT,
        workflow.runName,
        workflow.launchDir
    )

    ALPHAFOLD.out.metrics
        .collectFile(name: "template_indep_info.tsv", storeDir: "${params.OUT}/${workflow.runName}", keepHeader: true)
        .subscribe onNext: {
            if( params.EMAIL ) {
                try {
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
                catch( Exception e ) {
                    log.warn "Failed to send notification email to ${params.EMAIL}"
                    log.error e
                }
        }
        }
}
