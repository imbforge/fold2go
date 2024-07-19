process SHINY {
    tag "${workflow.userName}@localhost:${params.PORT}"

    when:
        params.SHINY.enabled

    input:
        env(SHINY_APP_NJOBS)
        env(SHINY_APP_DATA)
        env(SHINY_APP_RUN_NAME)
        env(SHINY_APP_LAUNCH_DIR)

    script:
        """
        shiny run \\
            --port=${params.PORT} \\
            --host=127.0.0.1 \\
            ${moduleDir}/resources/usr/bin/app.py
        """
}
