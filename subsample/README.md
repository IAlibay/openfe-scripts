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


### Example usage is: 

`python subsample.py --resultsdir ./results --time "500 ps" --timestep "4 fs"`

### An example output is:

```
leg (500 ps) ligand_i ligand_j repeat DG(i->j) (kcal/mol) MBAR uncertainty (kcal/mol)
complex lig_ejm_31 lig_ejm_42 0 -14.7 0.1
complex lig_ejm_31 lig_ejm_42 1 -15.0 0.1
complex lig_ejm_31 lig_ejm_42 2 -14.7 0.1
complex lig_ejm_31 lig_ejm_46 0 -40.2 0.2
complex lig_ejm_31 lig_ejm_46 1 -40.4 0.1
complex lig_ejm_31 lig_ejm_46 2 -40.6 0.2
complex lig_ejm_42 lig_ejm_43 0 -19.0 0.1
complex lig_ejm_42 lig_ejm_43 1 -18.6 0.2
complex lig_ejm_42 lig_ejm_43 2 -19.0 0.1
```

## Notes

* The timestep is used to calculate the total simulation time, this is set to the default of `4 fs`, please change if necessary.
* Units are required for `time` and `timestep`. An error will be thrown if `time` is greater than the simulated time.
* In some cases too low a `time` value may lead to a PyMBAR convergence error.


