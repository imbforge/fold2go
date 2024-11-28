process AF2_METRICS {
    tag "${meta}"
 
    when:
        params.METRICS.enabled

    input:
        tuple val(meta), path(prediction, stageAs: "chains/*"), path(fasta, stageAs: "chains.fasta")

    output:
        tuple val("alphafold2_${params.ALPHAFOLD2.MODEL_PRESET}"), path("template_indep_info.tsv"), emit: metrics
        path("*_contacts.pse"), optional: true, emit: contacts

    script:
        """
        python ${moduleDir}/resources/usr/bin/calculate_template_independent_metrics.py \\
            --model_preset=${params.ALPHAFOLD2.MODEL_PRESET} \\
            --prediction_dir=chains \\
            --project_name=${workflow.runName} \\
            --prediction_name=${meta*.value.join('.')}
        """
}

process AF3_METRICS {
    tag "${meta}"
 
    when:
        params.METRICS.enabled

    input:
        tuple val(meta), path(prediction)

    output:
        tuple val("alphafold3"), path("*_metrics.tsv"), emit: metrics

    script:
        """
        #!/usr/bin/env python

        import json
        import pandas
        from pathlib import Path

        models = []

        for model in Path('${prediction}').glob('**/seed-*/summary_confidences.json'):
            df = pandas.read_json(model, precise_float=True).select_dtypes(exclude=['object'])
            df.insert(0, 'project_name', '${workflow.runName}')
            df.insert(1, 'prediction_name', model.parent.parent.name)
            df.insert(2, 'model_preset', 'alphafold3')
            df.insert(3, 'model_id', model.parent.name)
            models.append(df)
        
        pandas.concat(models).to_csv('${meta.id}_metrics.tsv', sep='\t', index=False)
        """
}
