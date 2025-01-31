import json
import os
import polars
import signal
import numpy as np
import plotly.graph_objects as go

from ipymolstar import PDBeMolstar
from pathlib import Path
from shiny import reactive
from shiny.express import input, render, ui
from shinywidgets import render_plotly, render_widget

# set some shiny page options
ui.busy_indicators.use(pulse=False)
ui.page_opts(
    fillable=True,
    fillable_mobile=True
)

## get environment and store vars
njobs = int(os.getenv('SHINY_APP_NJOBS'))
results_dir = Path(os.getenv("SHINY_APP_DATA")) / os.getenv("SHINY_APP_RUN_NAME")
log_file = Path(os.getenv("SHINY_APP_LAUNCH_DIR")) / '.nextflow.log'

## define helper functions

# load pae metric for selected prediction
def _get_pae(record: dict) -> dict:
    match record.get('model_preset').split('_')[0]:
        case 'alphafold2':
            with (results_dir / 'predictions' / 'alphafold2' / record.get('prediction_name') / f"pae_{record.get('model_id')}.json").open('r') as fin:
                return np.array(json.load(fin)[0].get('predicted_aligned_error'), dtype=np.float16)
        case 'alphafold3':
            with (results_dir / 'predictions' / 'alphafold3' / record.get('prediction_name') / record.get('model_id') / 'confidences.json').open('r') as fin:
                return np.array(json.load(fin).get('pae'), dtype=np.float16)
        case 'boltz':
            with np.load(results_dir / 'predictions' / 'boltz' / record.get('prediction_name') / f"pae_{record.get('prediction_name')}_{record.get('model_id')}.npz") as fin:
                return np.array(fin.get('pae'), dtype=np.float16)

# load 3D model for selected prediction
def _get_model(record: dict) -> dict:
    match record.get('model_preset').split('_')[0]:
        case 'alphafold2':
            model = results_dir / 'predictions' / 'alphafold2' / record.get('prediction_name') / f"{record.get('model_rank')}.cif"
        case 'alphafold3':
            model = results_dir / 'predictions' / 'alphafold3' / record.get('prediction_name') / record.get('model_id') / 'model.cif'
        case 'boltz':
            model = results_dir / 'predictions' / 'boltz' / record.get('prediction_name') /  f"{record.get('prediction_name')}_{record.get('model_id')}.cif"
    return {
        'data'  : model.read_text(),
        'format': 'cif',
        'binary': False
    }

# list files of completed predictions
def _list_files() -> list:
    return list( (results_dir / 'metrics').glob("*_metrics.tsv") )

# parse log file to retrieve pipeline progress
# TODO: replace this with an http endpoint and use nf-weblog
def _parse_log() -> list:
    msg = []
    with log_file.open('r') as txt:
        for line in txt:
            if 'INFO  nextflow.Session' in line:
                msg.append(line.split(' - ')[-1])
    return msg

## define reactive values and functions

# initialize reactive values to hold metrics dataframe and selected row
df, selection = reactive.value(), reactive.value()

# construct data frame from list of files, check for changes every 30s
@reactive.effect
@reactive.poll(_list_files, 30)
def _():
    completed = _list_files()
    ui.update_tooltip("metrics_tip", show=(not completed))
    if completed:
        metrics = polars.concat(
            [polars.read_csv(tsv, separator='\t') for tsv in completed],
            how="diagonal"
        )
        df.set(metrics)

# store row selection
@reactive.effect
def _():
    selected = render_frame.data_view(selected=True).to_dicts()
    ui.update_tooltip("structure_tip", show=(not selected))
    ui.update_tooltip("pae_tip", show=(not selected))
    if selected:
        row = selected.pop()
        # collect chain lengths from column names (chain<single letter>_length) and store them as {<single letter>: <len>}
        row['chain_info'] = {
            label.split('_')[0].removeprefix('chain'): length 
            for label, length in row.items()
            if label.endswith('_length') and length is not None
        }
        selection.set(row)

# get selected residue pairs from PAE plot and highlight them in structure
@reactive.effect
def _():
    def _get_query_param(residue, chains=selection().get('chain_info')):
        idx = 0
        for chain, length in chains.items():
            if residue <= (idx + length):
                return {
                    "struct_asym_id": chain,
                    "residue_number": residue - idx,
                }
            idx += length
        else:
            raise ValueError(f"Residue {residue} is out of range")

    def _hover_callback(trace, points, _):
        render_structure.widget.highlight = {
            "data": [_get_query_param(res + 1) for res in points.point_inds.pop()],
            "focus": True
        }
    render_pae.widget.data[0].on_hover(_hover_callback)

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
            "Results will appear here as they become available"
    @render.data_frame
    def render_frame():
        progress = ui.Progress(min=0, max=njobs)
        if not df.is_set():
            # display modal with log file content until predictions become available
            ui.modal_show(
                ui.modal(
                    ui.tags.pre(_parse_log()),
                    title=ui.HTML('<a href="https://gitlab.rlp.net/imbforge/fold2go">imbforge/fold2go</a> is currently running, please check back later to see some results...'),
                    size='xl',
                    footer=[ui.modal_button("Dismiss"), ui.span(class_='spinner-border')]
                )
            )
            progress.set(
                value=0,
                message="Predictions are running",
                detail=f"(0/{njobs} complete)"
            )
        else:
            ui.modal_remove()
            progress.set(
                value=(ncomplete := df().n_unique(subset=['prediction_name', 'model_preset'])),
                message="Predictions are running" if ( ncomplete < njobs ) else "Predictions are complete",
                detail=f"({ncomplete}/{njobs} complete)"
            )
        return render.DataGrid(df(), selection_mode='row')

    with ui.card_footer(class_='text-center'):
        # display button to download metrics table
        with ui.tooltip(id="download_tip"):
            @render.download(label="Download", filename="template_indep_info.csv")
            def download_metrics():
                yield df().write_csv()
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
            fig = go.Figure(
                data = {
                    'type': 'heatmap',
                    'colorbar': {'title': 'PAE [Å]'},
                    'colorscale': 'Greens_r',
                    'z': (pae := _get_pae(selection())),
                    'x0': 1,
                    'y0': 1,
                    'hovertemplate': 'Scored residue: %{x}<br>Aligned residue: %{y}<br>Predicted Aligned Error: %{z:.2f} Å',
                    'zmin': 0.0,
                    'zmax': 31.75,
                    'name': ''
                },
                layout = {
                    'xaxis' : {'title': 'Scored residue', 'constrain': 'domain'},
                    'yaxis' : {'title': 'Aligned residue', 'autorange': 'reversed', 'scaleanchor': 'x'},
                }
            )
            # draw chain boundaries as dashed lines
            idx = 0.5
            for chain in selection().get('chain_info').values():
                idx += chain
                fig.add_shape(
                    type = "line",
                    x0 = idx,
                    x1 = idx,
                    y0 = 0,
                    y1 = pae.shape[1],
                    line = {'color': 'red', 'dash': 'dot', 'width': 1.5}
                )
                fig.add_shape(
                    type = "line",
                    x0 = 0,
                    x1 = pae.shape[0],
                    y0 = idx,
                    y1 = idx,
                    line = {'color': 'red', 'dash': 'dot', 'width': 1.5}
                )
            return fig

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
