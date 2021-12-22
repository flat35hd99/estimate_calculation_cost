# Estimate calculation cost

## Calculation process

1. Prepare a structure
2. Optimize a structure
3. Heat your system
4. Get structures to use as initial structure in NVE simulation
5. Get trajectories to get flux and conductivities using CURP
6. Get flux and conductivities

This script calculate step6.

## Usage

```
$ python cli.py --pdb-code 1VII
```

You can get `pdb1vii.ent.csv` that record all residue pairs in your system without [these residues](https://github.com/flat35hd99/estimate_calculation_cost/blob/a9664a7dbc210793e4a6e36ec70c6380b325f098/cli.py#L34).
