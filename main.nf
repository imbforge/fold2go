#!/usr/bin/env nextflow

include { validateParameters } from 'plugin/nf-schema'

include { FOLD2GO } from './workflows/fold2go'

workflow {

  main:
    validateParameters()

    FOLD2GO()
    // when a run is terminated via GUI, the jupyter-server is not stopped, which leads to a 500 error
    // workaround this by manually stopping the jupyter-server, which will shut the proxy down as well
    // FIXME: once https://github.com/jupyterhub/jupyter-server-proxy/pull/395 is merged, this can be handled programmatically
    if ( System.getenv('JUPYTERHUB_SERVICE_URL') ) {
      workflow.onComplete = { "jupyter-server stop ${new URL(System.getenv('JUPYTERHUB_SERVICE_URL')).getPort()}".execute() }
    }

  publish:
    msa: Channel<Tuple<Map, Path>> = FOLD2GO.out.msa
    predictions: Channel<Tuple<Map, Path>> = FOLD2GO.out.predictions
    metrics: Channel<Tuple<Map, Path>> = FOLD2GO.out.metrics
}

output {
  msa: Channel<Tuple<Map, Path>> {
    path { meta, msa ->
      msa >> "${workflow.runName}/msa/${meta.id}/${meta.model}/${msa.name}"
    }
    mode 'copy'
    enabled params.SAVE_MSA
  }
  predictions: Channel<Tuple<Map, Path>> {
    path { meta, prediction ->
      prediction >> "${workflow.runName}/predictions/${meta.model}/${meta.id}"
    }
    mode 'copy'
  }
  metrics: Channel<Tuple<Map, Path>> {
    path { meta, metrics -> 
      metrics >> "${workflow.runName}/metrics/${meta.id}.${meta.model}_metrics.tsv"
    }
    mode 'copy'
  }
}