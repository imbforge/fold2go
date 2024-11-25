import json
import os
import pandas
import signal
import plotly.express as px

from ipymolstar import PDBeMolstar
from pathlib import Path
from shiny import reactive
from shiny.express import input, render, ui
from shinywidgets import render_plotly, render_widget

## get environment and store vars

njobs = int(os.getenv('SHINY_APP_NJOBS'))
results_dir = Path(os.getenv("SHINY_APP_DATA")) / os.getenv("SHINY_APP_RUN_NAME")
log_file = Path(os.getenv("SHINY_APP_LAUNCH_DIR")) / '.nextflow.log'

## define helper functions

# load pae metric for selected prediction
def _get_pae(record: dict) -> list:
    if 'ranking_score' in record:
        with (results_dir / 'predictions' / record.get('prediction_name') / record.get('model_id') / 'confidences.json').open() as fin:
            return json.load(fin).get('pae')
    else:
        with (results_dir / 'predictions' / record.get('prediction_name') / f"pae_{record.get('model_id')}.json").open() as fin:
            return json.load(fin)[0].get('predicted_aligned_error')

# load 3D model for selected prediction
def _get_model(record: dict) -> dict:
    if 'ranking_score' in record:
        model = results_dir / 'predictions' / record.get('prediction_name') / record.get('model_id') / 'model.cif'
    else:
        model = results_dir / 'predictions' / record.get('prediction_name') / f"{record.get('model_rank')}.pdb"
    return {
        'data'  : model.read_text(),
        'format': model.suffix[1:],
        'binary': False
    }

# list files of completed predictions
def _list_files() -> list:
    return list( (results_dir / 'metrics').glob("*_metrics.tsv") )

# parse log file to retrieve pipeline progress
# TODO: replace this with an http endpoint and use nf-weblog
def _parse_log() -> list:
    msg = []
    with log_file.open() as txt:
        for line in txt:
            if 'INFO  nextflow' in line:
                msg.append(line.split(' - ')[-1])
    return msg

## define reactive values and functions

# initialize reactive values to hold metrics dataframe and selected row
df, selection = reactive.value(), reactive.value()

# construct data frame from list of files, check for changes every 30s
@reactive.effect
@reactive.poll(_list_files, 30)
def _():
    if (files := _list_files()):
        models = pandas.concat([pandas.read_csv(f, sep='\t') for f in files])
        models.attrs['done'] = len(files)
        df.set(models)

# store row selection
@reactive.effect
def _():
    selected = render_frame.data_view(selected=True).to_dict(orient='records')
    ui.update_tooltip("structure_tip", show=(not selected))
    ui.update_tooltip("pae_tip", show=(not selected))
    if selected:
        selection.set(selected.pop())

# terminate shiny server on button press
@reactive.effect
@reactive.event(input.terminate)
def _():
    return os.kill(os.getpid(), signal.SIGUSR1)

## define user interface

# display dataframe with calculated metrics and selectable rows
with ui.card():
    with ui.card_header(class_='text-center'):
        with ui.tooltip(placement="bottom", id="metrics_tip"):
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

    with ui.card_footer(class_='text-center'):
        # display button to download metrics table
        with ui.tooltip(id="download_tip"):
            @render.download(label="Download", filename="template_indep_info.csv")
            def download_metrics():
                yield df().to_csv(index=False)
            "Press this button to download metrics table as csv"

        # display button to terminate shiny server process
        with ui.tooltip(id="terminate_tip"):
            ui.input_action_button("terminate", "Shutdown", class_="btn-danger")
            "Press this button to terminate fold2go. Unfinished predictions will be lost"


with ui.layout_columns(col_widths=(4,8)):
    # plot predicted aligned error based on selection in first card
    with ui.card():
        with ui.card_header(class_='text-center'):
            with ui.tooltip(placement="bottom", id="pae_tip"):
                "Predicted Aligned Error"
                "Select row to display corresponding PAE plot"
        @render_plotly
        def render_pae():
            return px.imshow(
                _get_pae(selection()),
                labels = {'color': 'PAE'},
                color_continuous_scale = px.colors.sequential.Viridis
            )
    # render model structure based on selection in first card
    with ui.card():
        with ui.card_header(class_='text-center'):
            with ui.tooltip(placement="bottom", id="structure_tip"):
                "Predicted Structure"
                "Select row to display corresponding 3D structure"
        @render_widget
        def render_structure():
            return PDBeMolstar(
                custom_data = _get_model(selection()),
                alphafold_view = True,
                sequence_panel = True,
                hide_animation_icon = True
            )
