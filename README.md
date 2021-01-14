# parody_py

Python package to process data from PARODY-JA4.3 dynamo simulations.

## Getting Started

### Prerequisites
- [Python](https://www.python.org/)

### Installing
Conda:
```
conda install -c jnywong parody_py
```

Pip:
```
pip install parody_py
```

Git:

Find the Git repo [here](https://github.com/jnywong/nondim-slurry).

### Example script: diagnostics

Sample scripts can be found within the module package `parody_py/scripts`.

1. Open `parody_py/scripts/diagnostics.py`

2. Folder structure should be of the form `<folder>/<run_ID>/Gt_*.run_ID`. Set path to simulation data by setting

```
folder = <folder>
run_ID = <run_ID>
```

3. Run `parody_py/scripts/diagnostics.py`

4. Admire the output:

### Example script: meridional snapshots

1. Open `parody_py/scripts/meridional_snapshot.py`

2. Folder structure should be of the form `<folder>/<run_ID>/Gt_*.run_ID`. Set path to simulation data by setting

```
folder = <folder>
run_ID = <run_ID>
```

3. Specify timestamp of snapshot by setting `timestamp`

4. Run `parody_py/scripts/meridional_snapshot.py`

5. Admire the output:

### Example script: surface snapshots

1. Open `parody_py/scripts/surface_snapshot.py`

2. Folder structure should be of the form `<folder>/<run_ID>/Gt_*.run_ID`. Set path to simulation data by setting

```
folder = <folder>
run_ID = <run_ID>
```

3. Specify timestamp of snapshot by setting `timestamp`

4. Run `parody_py/scripts/surface_snapshot.py`

5. Admire the output:

## Links

* [PyPI](https://test.pypi.org/project/slurpy/)

## Authors

* **Jenny Wong** - *Institut de Physique du Globe de Paris/Institut des Sciences de la Terre*

## License

This project is licensed under the MIT License - see the [license.md](LICENSE.md) file for details

## Acknowledgments

* Del Duca Foundation
* ERC SEIS

:tada:
