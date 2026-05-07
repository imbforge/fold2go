#!/usr/bin/env nextflow

nextflow.enable.types = true

include { validateParameters } from 'plugin/nf-schema'

include { FOLD2GO } from './workflows/fold2go'

workflow {

  main:
    validateParameters()

    input =
      channel.fromPath(params.IN)
      .map { input ->
        record(model: params.MODEL_PRESET, input: input)
      }

    result = FOLD2GO(input)
   
  publish:
    msa: Channel<Tuple<Map, Path>> = result.msa
    predictions: Channel<Record> = result.predictions
    metrics: Channel<Record> = result.metrics

  onError:
    // when a run is terminated via GUI, the jupyter-server is not stopped, which leads to a 500 error
    // workaround this by manually stopping the jupyter-server, which will shut the proxy down as well
    // FIXME: once https://github.com/jupyterhub/jupyter-server-proxy/pull/395 is merged, this can be handled programmatically
    if( System.getenv('JUPYTERHUB_SERVICE_URL') ) {
      "jupyter-server stop ${new URL(System.getenv('JUPYTERHUB_SERVICE_URL')).getPort()}".execute()
    }
}

output {
  msa: Channel<Record> {
    path { record ->
      record.msa >> "${workflow.runName}/msa/${record.id}/${record.model}/${record.msa.name}"
    }
    mode 'copy'
    enabled params.SAVE_MSA
  }
  predictions: Channel<Record> {
    path { record ->
      record.prediction >> "${workflow.runName}/predictions/${record.model}/${record.id}"
    }
    mode 'copy'
  }
  metrics: Channel<Record> {
    path { record -> 
      record.metrics >> "${workflow.runName}/metrics/${record.id}.${record.model}_metrics.tsv"
    }
    mode 'copy'
  }
}
