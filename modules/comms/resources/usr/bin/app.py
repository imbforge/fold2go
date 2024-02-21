import os
import signal
import pandas
import py3Dmol

from pathlib import Path
from shiny import reactive, req
from shiny.express import input, render, ui

# get environment and store vars
njobs = int(os.getenv('SHINY_APP_NJOBS'))
metrics = Path(f"{os.getenv('SHINY_APP_DATA')}/metrics")
predictions = Path(f"{os.getenv('SHINY_APP_DATA')}/predictions")
header = {
    'project_name':str,
    'prediction_name':str,
    'chainA_length':int,
    'chainB_length':int,
    'model_id':str,
    'model_confidence':float,
    'chainA_intf_avg_plddt':float,
    'chainB_intf_avg_plddt':float,
    'intf_avg_plddt':float,
    'pDockQ':float,
    'iPAE':float,
    'num_chainA_intf_res':int,
    'num_chainB_intf_res':int,
    'num_res_res_contact':int,
    'num_atom_atom_contact':int
}

# list files of completed predictions
def _list_files() -> list:
    return list(metrics.glob("*_template_indep_metrics.tsv"))

# construct data frame from list of files, check for changes every 60s
@reactive.poll(_list_files, 60)
def _build_frame() -> pandas.DataFrame:
    if (metrics := _list_files()):
        data = []
        for f in metrics:
            data.append(pandas.read_csv(f, sep='\t', usecols=header))
        return pandas.concat(data)
    else: # as long as no predictions are available, just return the header
        return pandas.DataFrame(columns=header)

with ui.layout_columns(col_widths=(12, 12)):

    # first card shows dataframe with calculated metrics and selectable rows
    with ui.card():
        @render.data_frame
        def render_frame():

            progress = ui.Progress(min=0, max=njobs)

            if not (done := len(_list_files())):
                progress.set(None, message="AlphaFold is running", detail=f"({done}/{njobs} complete)")
                msg = ui.modal(
                    f"This could take a while, please check back later to see some results",
                    title=f"AlphaFold is running {njobs} predictions",
                    size='m',
                    footer=ui.HTML('<div class="spinner-border"></div>')
                )
                ui.modal_show(msg)
            else:
                progress.set(done, message="AlphaFold is running", detail=f"({done}/{njobs} complete)")
                ui.modal_remove()

            return render.DataGrid(_build_frame(), row_selection_mode="single")

        @render.download(label="Download", filename="template_indep_info.csv")
        def download_metrics():
            yield _build_frame().to_csv(index=False)

    # second card displays model structure based on selection in first card
    with ui.card():
        @render.ui
        def render_pdb():
            idx = list(req(input.render_frame_selected_rows()))[0]
            row = _build_frame().iloc[idx].to_dict()
            with open(f'{predictions}/{row.get('model_id')}.pdb') as pdb:
                model = "".join([i for i in pdb])
            view = py3Dmol.view(width=1200, height=800)
            view.addModelsAsFrames(model)
            view.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
            view.zoomTo()
            return ui.HTML(view._make_html())

    # third card provides button to terminate shiny server process
    with ui.card():
        ui.input_action_button("exit", "Exit", class_="btn-danger")
        @reactive.effect
        @reactive.event(input.exit)
        def _():
           return os.kill(os.getpid(), signal.SIGUSR1)