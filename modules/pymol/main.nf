process PYMOL {
    tag "${meta}"
 
    when:
        (params.PYMOL.enabled && meta.size() == 2)

    input:
        tuple val(meta), path(prediction, stageAs: "chains/*"), path(fasta, stageAs: "chains.fasta")

    output:
        path("template_indep_info.tsv"), emit: metrics

    script:
        """
        python ${moduleDir}/resources/usr/bin/calculate_template_independent_metrics.py \\
            -path_to_prediction='./' \\
            -project_name='${workflow.runName}' \\
            -prediction_name='${meta*.value.join('.')}' \\
            -skip_write_out_contacts
        """
}
