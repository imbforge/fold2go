process MSA {
    tag "${meta}"

    when:
        params.MSA.enabled

    input:
        tuple val(meta), path(yaml, stageAs: 'input/*')

    output:
        tuple val(meta), path("*.yaml"), path("*.csv"), emit: msa

    script:
        """
        #!/usr/bin/env python

        # prefetch msas, so we don't waste gpu time on http requests

        import yaml
        from boltz.main import compute_msa
        from pathlib import Path

        with Path("${yaml}").open('r') as fin:
            targets = yaml.safe_load(fin)

        sequences = {}

        for entity in targets['sequences']:
            if 'protein' in entity and entity['protein'].get('msa') is None:
                chain_id, seq = entity['protein']['id'], entity['protein']['sequence']
                if isinstance(chain_id, list):
                    msa_id = f"${meta.id}_{'_'.join(chain_id)}"
                else :
                    msa_id = f"${meta.id}_{chain_id}"
                sequences[msa_id] = seq
                entity['protein']['msa'] = f"{msa_id}.csv"

        compute_msa(
            data=sequences,
            target_id="${meta.id}",
            msa_dir=Path.cwd(),
            msa_server_url="${params.BOLTZ.MSA_SERVER_URL}",
            msa_pairing_strategy="${params.BOLTZ.MSA_PAIRING_STRATEGY}"
        )

        with Path("${meta.id}.yaml").open("w") as fout:
            yaml.dump(targets, fout)
        """
}

process INFERENCE {
    tag "${meta}"
    label "gpu"

    when:
        params.INFERENCE.enabled

    input:
        tuple val(meta), path(yaml), path(msa)

    output:
        tuple val(meta), path("boltz_results_*/predictions/${meta.id}", type: 'dir'), emit: prediction

    script:
        """
        boltz predict ${yaml} \\
            --write_full_pae \\
            --recycling_steps=${params.BOLTZ.RECYCLING_STEPS} \\
            --sampling_steps=${params.BOLTZ.SAMPLING_STEPS} \\
            --diffusion_samples=${params.BOLTZ.DIFFUSION_SAMPLES} \\
            --cache=${workDir} \\
            --out_dir=.
        """

}
