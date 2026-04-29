include { ALPHAFOLD2 } from '../../subworkflows/alphafold2'
include { ALPHAFOLD3 } from '../../subworkflows/alphafold3'
include { BOLTZ      } from '../../subworkflows/boltz'
include { COLABFOLD  } from '../../subworkflows/colabfold'
include { SHINY      } from '../../modules/shiny'
include { METRICS    } from '../../modules/metrics'

workflow FOLD2GO {

    main:
        input = channel.fromPath(params.IN)

        COLABFOLD(
            input.filter { it -> params.BATCH_MODE && it =~ /.(fasta|fa)$/ }
        )

        ALPHAFOLD2(
            input.filter { it -> !params.BATCH_MODE && it =~ /.(fasta|fa)$/ }
        )

        ALPHAFOLD3(
            input.filter { it -> it =~ /.json$/ }
        )

        BOLTZ(
            input.filter { it -> it =~ /.(yaml|yml)$/ }
        )

        jobcount =
            channel.empty()
            .mix(ALPHAFOLD2.out.jobcount, ALPHAFOLD3.out.jobcount, BOLTZ.out.jobcount, COLABFOLD.out.jobcount)
            .sum()
            .collectFile { njobs ->
                [
                "shiny_config.json",
                """
                {
                    "njobs": ${njobs},
                    "data": "${workflow.outputDir}/${workflow.runName}",
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
            ALPHAFOLD2.out.prediction.mix(ALPHAFOLD3.out.prediction).mix(BOLTZ.out.prediction).mix(COLABFOLD.out.prediction)
        )
        
        METRICS.out.metrics
        .collectFile ( storeDir: "${workflow.outputDir}/${workflow.runName}", keepHeader: true ) {
            meta, metrics -> [ "${meta.model}_metrics.tsv", metrics ]
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
                        Results of this run have all been stored at ${workflow.outputDir}/${workflow.runName}.

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

    emit:
        msa: Channel<Tuple<Map, Set<Path>>> = ALPHAFOLD2.out.msa.mix(ALPHAFOLD3.out.msa).mix(BOLTZ.out.msa).mix(COLABFOLD.out.msa)
        predictions: Channel<Tuple<Map, Set<Path>>> = ALPHAFOLD2.out.prediction.mix(ALPHAFOLD3.out.prediction).mix(BOLTZ.out.prediction).mix(COLABFOLD.out.prediction)
        metrics: Channel<Tuple<Map, Path>> = METRICS.out.metrics
}
