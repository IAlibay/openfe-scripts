# Subsampling RBFEs

A script for subsampling RBFE results from a directory of openfe result JSON files

## Usage

Usage: subsample.py [OPTIONS]

Options:
  --resultsdir DIRECTORY  path to openfe JSON result files directory
                          [required]
  --time TEXT             Simulation time (with units) to analyze production
                          energies for. Example: `1 ns` or `100 ps`.
                          [required]
  --timestep TEXT         Simulation timestep (with units). Default to `4 fs`
  -o FILENAME
  --help                  Show this message and exit.

## Notes

* The timestep is used to calculate the total simulation time, this is set to the default of `4 fs`, please change if necessary.
* Units are required for `time` and `timestep`. An error will be thrown if `time` is greater than the simulated time.
* In some cases too low a `time` value may lead to a PyMBAR convergence error.


