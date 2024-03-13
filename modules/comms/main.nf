process COMMS {
    tag "localhost:${params.COMMS.port}"

    when:
        params.COMMS.enabled

    input:
        env(SHINY_APP_NJOBS)
        env(SHINY_APP_DATA)
        env(SHINY_APP_LAUNCH_DIR)

    script:
        """
        shiny run \\
            --port=${params.COMMS.port} \\
            --host=0.0.0.0 \\
            ${moduleDir}/resources/usr/bin/app.py
        """
}
