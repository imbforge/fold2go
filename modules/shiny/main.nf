process SHINY {
    tag "${workflow.userName}@localhost:${SHINY_APP_PORT}"

    when:
        params.SHINY.enabled

    input:
        val(SHINY_APP_PORT)
        env(SHINY_APP_NJOBS)
        env(SHINY_APP_DATA)
        env(SHINY_APP_RUN_NAME)
        env(SHINY_APP_LAUNCH_DIR)

    script:
        """
        shiny run \\
            --port=${SHINY_APP_PORT} \\
            --host=127.0.0.1 \\
            ${moduleDir}/resources/usr/bin/app.py
        """
}
