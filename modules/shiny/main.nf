nextflow.preview.types = true

process SHINY {
    tag "${workflow.userName}@localhost:${socket}"

    input:
    socket: String
    _json: Path

    when:
    params.SHINY.enabled

    script:
    """
    #!/usr/bin/env python

    from shiny import run_app

    run_app("${moduleDir}/resources/usr/bin/app.py:app", uds="${socket}")
    """
}
