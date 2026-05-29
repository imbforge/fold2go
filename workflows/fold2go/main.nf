nextflow.enable.types = true

include { ALPHAFOLD2 } from '../../subworkflows/alphafold2'
include { ALPHAFOLD3 } from '../../subworkflows/alphafold3'
include { BOLTZ      } from '../../subworkflows/boltz'
include { COLABFOLD  } from '../../subworkflows/colabfold'
include { SHINY      } from '../../modules/shiny'
include { METRICS    } from '../../modules/metrics'

workflow FOLD2GO {

    take:
        input: Channel<Record>

    main:

        alphafold2 = ALPHAFOLD2(
            input.filter { record ->
                (record.model == 'alphafold2_multimer')
            }
        )

        alphafold3 = ALPHAFOLD3(
            input.filter { record ->
                (record.model == 'alphafold3')
            }
        )

        boltz = BOLTZ(
            input.filter { record ->
                (record.model in ['boltz1', 'boltz2'])
            }
        )

        colabfold = COLABFOLD(
            input.filter { record ->
                (record.model == 'colabfold')
            }
        )

        SHINY(
            channel.empty()
                .mix(alphafold2.jobcount).mix(alphafold3.jobcount).mix(boltz.jobcount).mix(colabfold.jobcount)
                .collect()
                .map { counts ->
                    record(
                        njobs: "${counts.sum()}",
                        data: "${workflow.outputDir}/${workflow.runName}",
                        logfile: "${workflow.launchDir}/.nextflow.log"
                    )
                }
        )

        metrics = METRICS(
            alphafold2.prediction.mix(alphafold3.prediction).mix(boltz.prediction).mix(colabfold.prediction)
        )
        
        metrics
        .collectFile ( storeDir: "${workflow.outputDir}/${workflow.runName}", keepHeader: true ) {
            record -> [ "${record.model}_metrics.tsv", record.metrics ]
        }
        .collect()
        .map { tsv ->
            if( params.EMAIL ) {
                try {
                    sendMail (
                        to: "${params.EMAIL}",
                        subject: "fold2go (${workflow.runName})",
                        attach: tsv,
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
                    log.warn e.getMessage()
                }
            }
        }

    emit:
        msa: Channel<Tuple<Map, Set<Path>>> = alphafold2.msa.mix(alphafold3.msa).mix(boltz.msa).mix(colabfold.msa)
        predictions: Channel<Record> = alphafold2.prediction.mix(alphafold3.prediction).mix(boltz.prediction).mix(colabfold.prediction)
        metrics: Channel<Record> = metrics
}
