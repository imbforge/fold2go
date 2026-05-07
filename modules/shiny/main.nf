nextflow.enable.types = true

process SHINY {
    tag "${workflow.userName}@localhost:${workflow.workDir}/shiny.sock"

    input:
    record(
        njobs: String,
        data: String,
        logfile: String
    )

    stage:
    env 'FOLD2GO_NJOBS', njobs
    env 'FOLD2GO_DATA', data
    env 'FOLD2GO_LOG', logfile

    script:
    """
    #!/usr/bin/env python

    from shiny import run_app

    run_app("${moduleDir}/resources/usr/bin/app.py:app", uds="${workflow.workDir}/shiny.sock", ws="websockets-sansio")
    """
}
