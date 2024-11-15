process MSA {
    tag "${meta}"

    when:
        params.MSA.enabled

    input:
        tuple val(meta), path(json)

    output:
        tuple val(meta), path("msa/**/*.json"), emit: json

    script:
        """
        python /app/alphafold/run_alphafold.py \\
            --norun_inference \\
            --db_dir=${params.DATABASE} \\
            --json_path=${json} \\
            --output_dir=msa
        """
}

process INFERENCE {
    tag "${meta}"
    label "gpu"

    when:
        params.INFERENCE.enabled

    input:
        tuple val(meta), path(jobdef, stageAs: 'input.json'), path(msa, stageAs: 'msa/*')

    output:
        tuple val(meta), path("*.json"), emit: jobdef
        //tuple val(meta), path("seed-*"), emit: prediction

    script:
        """
        python ${moduleDir}/resources/usr/bin/collect_msas.py --json_path=${jobdef}
        """

//      python /app/alphafold/run_alphafold.py \\
//          --json_path=${jobdef}
//          --norun_data_pipeline \\
//          --db_dir=${params.DATABASE} \\
//          --model_dir=${params.WEIGHTS}/${params.MODEL_PRESET} \\
//          --json_path=${json} \\
//          --output_dir=predictions
}