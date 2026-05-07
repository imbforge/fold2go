nextflow.enable.types = true

process MSA {
    tag "${input.simpleName}"

    input:
    record(
        model: String,
        input: Path
    )

    stage:
    stageAs input, 'input/*'

    output:
    record(
        id   : input.simpleName,
        msa  : files("${input.simpleName}*.{yaml,csv}"),
        model: model
    )

    script:
    """
    #!/usr/bin/env python

    # prefetch msas, so we don't waste gpu time on http requests

    import yaml
    from boltz.main import compute_msa
    from pathlib import Path

    with Path("${input}").open('r') as fin:
        targets = yaml.safe_load(fin)

    sequences = {}

    for entity in targets['sequences']:
        if 'protein' in entity and entity['protein'].get('msa') is None:
            chain_id, seq = entity['protein']['id'], entity['protein']['sequence']
            if isinstance(chain_id, list):
                msa_id = f"${input.simpleName}_{'_'.join(chain_id)}"
            else :
                msa_id = f"${input.simpleName}_{chain_id}"
            sequences[msa_id] = seq
            entity['protein']['msa'] = f"{msa_id}.csv"

    compute_msa(
        data=sequences,
        target_id="${input.simpleName}",
        msa_dir=Path.cwd(),
        msa_server_url="${params.BOLTZ.MSA_SERVER_URL}",
        msa_pairing_strategy="${params.BOLTZ.MSA_PAIRING_STRATEGY}"
    )

    with Path("${input.simpleName}.yaml").open("w") as fout:
        yaml.dump(targets, fout)
    """
}

process INFERENCE {
    tag "${id}"
    label "gpu"

    input:
    record(
        id   : String,
        msa  : Set<Path>,
        model: String
    )

    output:
    record(
        id        : id,
        prediction: file("boltz_results_*/predictions/${id}", type: 'dir'),
        model     : model
    )

    script:
    """
    boltz predict ${id}.yaml \\
        --write_full_pae \\
        --model=${model} \\
        --recycling_steps=${params.BOLTZ.RECYCLING_STEPS} \\
        --sampling_steps=${params.BOLTZ.SAMPLING_STEPS} \\
        --diffusion_samples=${params.BOLTZ.DIFFUSION_SAMPLES} \\
        --cache=${workflow.workDir} \\
        ${params.BOLTZ.USE_KERNELS ? '' : '--no_kernels'} \\
        --out_dir=.
    """
}
