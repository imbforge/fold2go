process PYMOL {
    tag "${meta}"
 
    when:
        params.PYMOL.enabled

    input:
        tuple val(meta), path(prediction, stageAs: "chains/*"), path(fasta, stageAs: "chains.fasta")

    output:
        path("template_indep_info.tsv"), emit: metrics
        path("*_contacts.pse"), optional: true, emit: contacts

    script:
        """
        python ${moduleDir}/resources/usr/bin/calculate_template_independent_metrics.py \\
            --path_to_prediction='chains' \\
            --project_name='${workflow.runName}' \\
            --prediction_name='${meta*.value.join('.')}'
        """
}
