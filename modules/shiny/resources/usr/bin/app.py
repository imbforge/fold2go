import json
import os
import signal
import pandas
import py3Dmol
import plotly.express as px

from pathlib import Path
from shiny import reactive, req
from shiny.express import input, render, ui
from shinywidgets import render_plotly

# get environment and store vars
njobs = int(os.getenv('SHINY_APP_NJOBS'))
metrics = Path(f'{os.getenv("SHINY_APP_DATA")}/{os.getenv("SHINY_APP_RUN_NAME")}/metrics')
predictions = Path(f'{os.getenv("SHINY_APP_DATA")}/{os.getenv("SHINY_APP_RUN_NAME")}/predictions')
logfile = Path(f'{os.getenv("SHINY_APP_LAUNCH_DIR")}/.nextflow.log')

# initialize surface styles for 3Dmol
styles = ['cartoon', 'stick', 'sphere']

# initialize color maps to select from
cmaps = {
    'viridis': px.colors.sequential.Viridis,
    'aggrnyl': px.colors.sequential.Aggrnyl,
    'portland': px.colors.diverging.Portland,
    'spectral': px.colors.diverging.Spectral_r
}

# list files of completed predictions
def _list_files() -> list:
    return list(metrics.glob("*.template_indep_metrics.tsv"))

# parse log file to retrieve pipeline progress
# TODO: replace this with an http endpoint and use nf-weblog
def _parse_log() -> list:
    msg = []
    with open(logfile) as txt:
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
        metrics = pandas.concat([pandas.read_csv(f, sep='\t') for f in files])
        metrics.attrs['done'] = len(files)
        df.set(metrics)

# display sidebar for (general) app settings
with ui.sidebar(position='left', open='closed'):
    ui.input_select('style', 'Surface style', {style:style for style in styles}, selected='cartoon')
    ui.input_select('cmap', 'Color scheme', {cmap:cmap for cmap in cmaps.keys()}, selected='viridis')

with ui.layout_columns(col_widths=(12, 6, 6)):
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
            with ui.tooltip(placement="right", id="pdb_tooltip"):
                "Model structure"
                "Select a row to display the corresponding 3D model"
        @render.ui
        def render_pdb():
            prediction = req(render_frame.data_view(selected=True).to_dict(orient='records')).pop()
            with open(f'{predictions}/{prediction.get('prediction_name')}/{prediction.get('model_rank')}.pdb') as pdb:
                model = "".join([res for res in pdb])
            view = py3Dmol.view(width=800, height=600)
            view.addModelsAsFrames(model)
            view.setStyle({'model': -1}, {input.style(): {'colorscheme': {'prop':'b', 'gradient':'linear', 'colors':list(reversed(cmaps.get(input.cmap()))), 'min':50, 'max':90}}})
            view.zoomTo()
            return ui.HTML(view._make_html())

    # plot predicted aligned error based on selection in first card
    with ui.card():
        with ui.card_header():
            with ui.tooltip(placement="right", id="pae_tooltip"):
                "Predicted aligned error"
                "Select a row to display the corresponding PAE plot"
        @render_plotly
        def render_pae():
            prediction = req(render_frame.data_view(selected=True).to_dict(orient='records')).pop()
            with open(f'{predictions}/{prediction.get('prediction_name')}/pae_{prediction.get('model_id')}.json') as fin:
                pae = json.load(fin)[0].get('predicted_aligned_error')
            return px.imshow(pae, labels={'color': 'PAE'}, color_continuous_scale=cmaps.get(input.cmap()))

    # display button to terminate shiny server process
    with ui.card():
        ui.input_action_button("terminate", "Terminate", class_="btn-danger")
        @reactive.effect
        @reactive.event(input.terminate)
        def _():
            return os.kill(os.getpid(), signal.SIGUSR1)