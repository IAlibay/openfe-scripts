import click
from collections import defaultdict
import csv
import glob
import pathlib
from openmmtools import multistate
from openff.units import unit
from openfe.protocols.openmm_utils.multistate_analysis import MultistateEquilFEAnalysis
from openfecli.commands.gather import load_results, is_results_json, get_type


def get_names(result) -> tuple[str, str]:
    # Result to tuple of ligand names
    nm = list(result['unit_results'].values())[0]['name']
    toks = nm.split()
    if toks[2] == 'repeat':
        return toks[0], toks[1]
    else:
        return toks[0], toks[2]


def _get_column(val):
    import numpy as np
    if val == 0:
        return 0

    log10 = np.log10(val)

    if log10 >= 0.0:
        col = np.floor(log10 + 1)
    else:
        col = np.floor(log10)
    return int(col)


def format_estimate_uncertainty(
    est: float,
    unc: float,
    unc_prec: int = 1,
) -> tuple[str, str]:
    import numpy as np
    # get the last column needed for uncertainty
    unc_col = _get_column(unc) - (unc_prec - 1)

    if unc_col < 0:
        est_str = f"{est:.{-unc_col}f}"
        unc_str = f"{unc:.{-unc_col}f}"
    else:
        est_str = f"{np.round(est, -unc_col + 1)}"
        unc_str = f"{np.round(unc, -unc_col + 1)}"

    return est_str, unc_str


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


def _get_files(storage: pathlib.Path, checkpoint=pathlib.Path):
    """
    Get the absolute storage path and the checkpoint path relative to it.

    Parameters
    ----------

    """
    if not storage.is_file():
        storage = storage.relative_to(storage.parts[0])
        checkpoint = checkpoint.relative_to(checkpoint.parts[0])

    
    assert (storage.is_file() and checkpoint.is_file())

    storage = storage.resolve()
    checkpoint = checkpoint.resolve().relative_to(storage.parent)

    return storage, checkpoint


def get_result_energies(results, time, timestep):
    
    pus = list(results['unit_results'].values())
    retvals = []
    
    for pu in pus:
        storage_file, checkpoint_file = _get_files(
            storage=pu['outputs']['nc'],
            checkpoint=pu['outputs']['last_checkpoint']
        )

        reporter = multistate.MultiStateReporter(
            storage=storage_file,
            open_mode='r',
            checkpoint_storage=checkpoint_file
        )
        analyzer = multistate.MultiStateSamplerAnalyzer(reporter)
        DG, dDG = get_free_energy_at_time(reporter, analyzer, time, timestep)
        retvals.append((DG, dDG))
        analyzer.clear()
        reporter.close()
    
    return retvals


def write_output(legs, output, time):
    writer = csv.writer(
        output,
        delimiter="\t",
        lineterminator="\n",
    )

    # Write header
    writer.writerow(
        [f"leg ({time})", "ligand_i", "ligand_j", "repeat",
         "DG(i->j) (kcal/mol)", "MBAR uncertainty (kcal/mol)"]
    )

    for ligpair, vals in sorted(legs.items()):
        for simtype, repeats in sorted(vals.items()):
            for i, (m, u) in enumerate(repeats):
                if m is None:
                    m, u = 'NaN', 'NaN'
                else:
                    m, u = format_estimate_uncertainty(m.m, u.m)

                writer.writerow([simtype, *ligpair, i, m, u])


@click.command
@click.option(
    '--resultsdir',
    type=click.Path(dir_okay=True, file_okay=False, path_type=pathlib.Path),
    required=True,
    help="path to openfe JSON result files directory"
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
@click.option(
    'output', '-o',
    type=click.File(mode='w'),
    default='-'
)
def analyze_all_results(resultsdir, time, timestep, output):
    # Find all possible jsons
    json_fns = glob.glob(str(resultsdir) + '/*json', recursive=True)
    
    # Filter only results jsons
    results_fns = filter(is_results_json, json_fns)
    
    # pair legs of simulations together into dict of dicts
    legs = defaultdict(dict)
    
    for result_fn in results_fns:
        result = load_results(result_fn)
        if result is None:
            continue
        elif result['estimate'] is None or result['uncertainty'] is None:
            errmsg = (f"WARNING: Calculations for {result_fn} "
                      "did not finish successfully!")
            click.echo(errmsg, err=True)
        
        try:
            names = get_names(result)
        except KeyError:
            raise ValueErorr("Failed to guess names")
        
        simtype = get_type(result)
        
        legs[names][simtype] = get_result_energies(result, time, timestep)

    write_output(legs, output, time)


if __name__ == "__main__":
    analyze_all_results()
