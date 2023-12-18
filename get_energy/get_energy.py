import click
from collections import defaultdict
import csv
import glob
import pathlib
from openmmtools import multistate
from openff.units import unit
from openfe.protocols.openmm_utils.multistate_analysis import MultistateEquilFEAnalysis
from openfecli.commands.gather import load_results, is_results_json, get_type


def get_free_energy_at_time(reporter, analyzer, time="1 ns", timestep="4 fs"):
    time_per_iteration = reporter._checkpoint_interval * unit(timestep)
    tot_time = analyzer.n_iterations * time_per_iteration
    
    if unit(time) > tot_time:
        errmsg = (f"The requested production analysis time {time} is greater "
                  f"than the total simulation time {tot_time}")
        raise ValueError(errmsg)
    
    requested_iters = round(unit(time).to('fs').m /  time_per_iteration.to('fs').m)
    
    # Overwrite the current max_n_iterations to subsample the dataset
    analyzer.max_n_iterations = requested_iters
    
    u_ln = analyzer._unbiased_decorrelated_u_ln
    N_l = analyzer._unbiased_decorrelated_N_l
    
    DG, dDG = MultistateEquilFEAnalysis._get_free_energy(
        analyzer, u_ln, N_l, unit.kilocalorie_per_mole
    )
    
    return DG, dDG


@click.command
@click.option(
    '--simulation',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    default="simulation.nc",
    help="path to the simulation `.nc` file [default is `simulation.nc`]"
)
@click.option(
    '--checkpoint',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    default="simulation.nc",
    help="path to the simulation `.nc` file [default is `simulation.nc`]"
)
@click.option(
    '--time',
    type=str,
    required=True,
    help=(
        "Simulation time (with units) to analyze production energies for. "
        "Example: `1 ns` or `100 ps`."
    )
)
@click.option(
    '--timestep',
    type=str,
    default="4 fs",
    help="Simulation timestep (with units). Default to `4 fs`"
)
def get_free_energy(simulation, checkpoint, time, timestep):
    """
    Get the absolute storage path and the checkpoint path relative to it.

    """
    if not simulation.is_file() or not checkpoint.is_file():
        errmsg = "Either the simulation or checkpoint file could nto be found"
        raise ValueError(errmsg)

    storage = simulation.resolve()
    checkpoint = checkpoint.resolve().relative_to(storage.parent)

    reporter = multistate.MultiStateReporter(
        storage=storage.as_posix(),
        open_mode='r',
        checkpoint_storage=checkpoint.as_posix()
    )
    analyzer = multistate.MultiStateSamplerAnalyzer(reporter)
    DG, dDG = get_free_energy_at_time(reporter, analyzer, time, timestep)

    print(f"DG: {DG}, MBAR error: {dDG}")


if __name__ == "__main__":
    get_free_energy()
