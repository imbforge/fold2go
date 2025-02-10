process SHINY {
    tag "${workflow.userName}@localhost:${SHINY_APP_PORT}"

    when:
        params.SHINY.enabled

    input:
        val(SHINY_APP_PORT)
        path(json)

    script:
        """
        shiny run \\
            --port=${SHINY_APP_PORT} \\
            --host=127.0.0.1 \\
            --log-level=debug \\
            ${moduleDir}/resources/usr/bin/app.py
        """
}
