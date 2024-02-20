process PYMOL {
    tag "${meta}"
 
    when:
        params.PYMOL.enabled

    input:
        tuple val(meta), path(fasta), path(prediction, stageAs: "prediction/*")

    output:
        path("template_indep_info.tsv"), emit: metrics

    script:
        """
        mv prediction ${meta}
        mv ${fasta} ${meta}.fasta

        python ${moduleDir}/resources/usr/bin/calculate_template_independent_metrics.py \\
            -path_to_prediction ./ \\
            -project_name ${workflow.runName} \\
            -skip_write_out_contacts
        """
}
