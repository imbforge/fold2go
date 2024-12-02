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

        models = {}

        for model in Path('${prediction}').glob('**/seed-*/summary_confidences.json'):
            with model.open('r') as fin:
                metrics = { metric: value for metric, value in json.load(fin).items() if not isinstance(value, dict | list) }
            info = {
                'project_name': '${workflow.runName}',
                'prediction_name': model.parent.parent.name,
                'model_preset': 'alphafold3',
                'model_id': model.parent.name
            }
            models[model.stem] = { **info, **metrics }
        
        pandas.DataFrame.from_dict(models, orient='index').to_csv('${meta.id}_metrics.tsv', sep='\t', index=False)
        """
}
