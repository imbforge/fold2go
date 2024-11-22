import json
import os
import pandas
import signal
import plotly.express as px

from ipymolstar import PDBeMolstar
from pathlib import Path
from shiny import reactive, req
from shiny.express import input, render, ui
from shinywidgets import render_plotly, render_widget

# get environment and store vars
njobs = int(os.getenv('SHINY_APP_NJOBS'))
results_dir = Path(os.getenv("SHINY_APP_DATA")) / os.getenv("SHINY_APP_RUN_NAME")
log_file = Path(os.getenv("SHINY_APP_LAUNCH_DIR")) / '.nextflow.log'

# construct paths to model structure and confidence metrics
# outputs are different for af2 and af3, so we need separate codepaths here
def _construct_paths(record: dict) -> dict:
    prediction = results_dir / 'predictions' / record.get('prediction_name')
    if 'ranking_score' in record.keys():
        # af3 output
        return {
            'structure': prediction / record.get('model_id') / 'model.cif',
            'confidences' : prediction / record.get('model_id') / 'confidences.json'
        }
    else:
        # af2 output
        return {
            'structure': prediction / f"{record.get('model_rank')}.pdb",
            'confidences': prediction / f"pae_{record.get('model_id')}.json"
        }

# list files of completed predictions
def _list_files() -> list:
    return list( (results_dir / 'metrics').glob("*_metrics.tsv") )

# parse log file to retrieve pipeline progress
# TODO: replace this with an http endpoint and use nf-weblog
def _parse_log() -> list:
    msg = []
    with open(log_file) as txt:
        for line in txt:
            if 'INFO  nextflow' in line:
                msg.append(line.split(' - ')[-1])
    return msg

# initialize reactive value to hold metrics dataframe
df = reactive.value()

# construct data frame from list of files, check for changes every 30s
@reactive.effect
@reactive.poll(_list_files, 30)
def _():
    if (files := _list_files()):
        models = pandas.concat([pandas.read_csv(f, sep='\t') for f in files])
        models.attrs['done'] = len(files)
        df.set(models)

# display dataframe with calculated metrics and selectable rows
with ui.card():
    with ui.card_header():
        with ui.tooltip(placement="right", id="metrics_tooltip"):
            "Metrics"
            "Results will appear here as they are produced"
    @render.data_frame
    def render_frame():
        progress = ui.Progress(min=0, max=njobs)
        if not df.is_set():
            # display modal with log file content until predictions become available
            progress.set(0, message="Predictions are running", detail=f"(0/{njobs} complete)")
            ui.modal_show(
                ui.modal(
                    ui.tags.pre(_parse_log()),
                    title=ui.HTML('<a href="https://gitlab.rlp.net/imbforge/fold2go">imbforge/fold2go</a> is currently running, please check back later to see some results...'),
                    size='xl',
                    footer=[ui.modal_button("Dismiss"), ui.span(class_='spinner-border')]
                )
            )
        else:
            ui.modal_remove()
            if (done := df().attrs['done']) < njobs:
                progress.set(df().attrs['done'], message="Predictions are running", detail=f"({df().attrs['done']}/{njobs} complete)")
            else:
                progress.set(done, message="Predictions are complete", detail=f"({done}/{njobs} complete)")
    
            return render.DataGrid(df(), selection_mode='row')

    @render.download(label="Download", filename="template_indep_info.csv")
    def download_metrics():
        yield df().to_csv(index=False)

# render model structure based on selection in first card
with ui.card():
    with ui.card_header():
        with ui.tooltip(placement="right", id="structure_tooltip"):
            "Model structure"
            "Select a row to display the corresponding 3D model"
    @render_widget
    def render_structure():
        prediction = req(render_frame.data_view(selected=True).to_dict(orient='records')).pop()
        structure = _construct_paths(prediction).get('structure')
        return PDBeMolstar(
            custom_data = {
                'data'  : structure.read_bytes(),
                'format': structure.suffix[1:],
                'binary': False
            },
            alphafold_view = True
        )

# plot predicted aligned error based on selection in first card
with ui.card():
    with ui.card_header():
        with ui.tooltip(placement="right", id="pae_tooltip"):
            "Predicted aligned error"
            "Select a row to display the corresponding PAE plot"
    @render_plotly
    def render_pae():
        prediction = req(render_frame.data_view(selected=True).to_dict(orient='records')).pop()
        with open(_construct_paths(prediction).get('confidences')) as fin:
            data = json.load(fin)
        pae = data[0].get('predicted_aligned_error') if isinstance(data, list) else data.get('pae')
        return px.imshow(pae, labels={'color': 'PAE'}, color_continuous_scale=px.colors.sequential.Viridis)

# display button to terminate shiny server process
with ui.card():
    ui.input_action_button("terminate", "Shutdown", class_="btn-danger")
    @reactive.effect
    @reactive.event(input.terminate)
    def _():
        return os.kill(os.getpid(), signal.SIGUSR1)
