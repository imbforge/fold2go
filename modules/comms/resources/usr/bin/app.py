import json
import os
import signal
import pandas
import py3Dmol
import pickle
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

# initialize reactive value to hold metrics dataframe
df = reactive.value()

# list files of completed predictions
def _list_files() -> list:
    return list(metrics.glob("*.template_indep_metrics.tsv"))

# construct data frame from list of files, check for changes every 60s
@reactive.effect
@reactive.poll(_list_files, 60)
def _update_frame() -> bool:
    if (metrics := _list_files()):
        data = []
        for f in metrics:
            data.append(pandas.read_csv(f, sep='\t'))
        return df.set(pandas.concat(data))

# get rowdata for a selected prediction from dataframe
def _get_prediction() -> dict:
    idx = list(req(input.render_frame_selected_rows()))[0]
    return df().iloc[idx].to_dict()

# parse log file to retrieve pipeline progress
# TODO: replace this with an http endpoint and use nf-weblog
def _parse_log() -> list:
    log = []
    with open(logfile) as txt:
        for line in txt:
            if 'INFO  nextflow' in line:
                log.append(line.split(' - ')[-1])
    return log

# check for new log messages every 60s
@reactive.poll(_parse_log, interval_secs=60)
def _render_modal() -> ui.Tag:
    return ui.modal(
        ui.tags.pre(_parse_log()),
        title=ui.HTML('<a href="https://gitlab.rlp.net/imbforge/fold2go">imbforge/fold2go</a> is currently running, please check back later to see some results...'),
        size='xl',
        footer=ui.div(class_='spinner-border')
    )

# display sidebar for (general) app settings
with ui.sidebar(position='left', open='closed'):
    ui.input_select('style', 'Surface style', {style:style for style in styles}, selected='cartoon')
    ui.input_select('cmap', 'Color scheme', {cmap:cmap for cmap in cmaps.keys()}, selected='viridis')

with ui.layout_columns(col_widths=(12, 6, 6)):
    # display dataframe with calculated metrics and selectable rows
    with ui.card():
        @render.data_frame
        def render_frame():
            if not (done := len(_list_files())):
                ui.modal_show(_render_modal())
            else:
                ui.modal_remove()
                progress = ui.Progress(min=0, max=njobs)
                progress.set(done, message="AlphaFold is running", detail=f"({done}/{njobs} complete)")

            return render.DataGrid(df(), row_selection_mode="single")

        @render.download(label="Download", filename="template_indep_info.csv")
        def download_metrics():
            yield df().to_csv(index=False)

    # render model structure based on selection in first card
    with ui.card():
        @render.ui
        def render_pdb():
            prediction = _get_prediction()
            with open(f'{predictions}/{prediction.get('prediction_name')}/{prediction.get('model_rank')}.pdb') as pdb:
                model = "".join([res for res in pdb])
            view = py3Dmol.view(width=800, height=600)
            view.addModelsAsFrames(model)
            view.setStyle({'model': -1}, {input.style(): {'colorscheme': {'prop':'b', 'gradient':'linear', 'colors':list(reversed(cmaps.get(input.cmap()))), 'min':50, 'max':90}}})
            view.zoomTo()
            return ui.HTML(view._make_html())

    # plot predicted aligned error based on selection in first card
    with ui.card():
        @render_plotly
        def render_pae():
            prediction = _get_prediction()
            with open(f'{predictions}/{prediction.get('prediction_name')}/pae_{prediction.get('model_id')}.json') as fin:
                pae = json.load(fin)[0].get('predicted_aligned_error')
            return px.imshow(pae, labels={'color': 'PAE'}, color_continuous_scale=cmaps.get(input.cmap()))

    # display button to terminate shiny server process
    with ui.card():
        ui.input_action_button("exit", "Exit", class_="btn-danger")
        @reactive.effect
        @reactive.event(input.exit)
        def _():
           return os.kill(os.getpid(), signal.SIGUSR1)
