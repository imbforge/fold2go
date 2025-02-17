process SHINY {
    tag "${workflow.userName}@localhost:${socket}"

    when:
        params.SHINY.enabled

    input:
        val(socket)
        path(json)

    script:
        """
        #!/usr/bin/env python

        from shiny import run_app

        run_app("${moduleDir}/resources/usr/bin/app.py:app", uds="${socket}")
        """
}
