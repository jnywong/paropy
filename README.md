[![PyPI version](https://badge.fury.io/py/paropy.svg)](https://badge.fury.io/py/paropy)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# paropy

Welcome! Here is my Python package to process data from PARODY-JA4.3 dynamo simulations with a stably-stratified F-layer at the base of the outer core. 

## Getting Started

### Prerequisites
- [Python3](https://www.python.org/)
- [SHTns](https://bitbucket.org/nschaeff/shtns/src/master/)
- [ChaosMagPy](https://github.com/ancklo/ChaosMagPy)

### Installing
Conda:
```
conda install -c jnywong paropy
```

Pip:
```
pip install paropy
```

Git: 

Find the Git repo [here](https://github.com/jnywong/nondim-slurry).

#### Warning: IPython issue
Please note that ipython may not be compatible with jedi release 0.18.0. Please revert to version 0.17.2 for auto-complete features to work (see [here](https://github.com/ipython/ipython/issues/12740)).

## Package structure
```
paropy/
  docs/
  paropy/
    __init__.py
    data_utils.py
    plot_utils.py
    routines.py
    scripts/
      convective_power.py
      diagnostics.py
      diagnostic_parameters.py
      filter_surface_field.py
      latitude_vs_Br.py
      latitude_vs_Br_with_CHAOS.py
      meridional_snapshot.py
      meridional_timeavg.py
      rotation_rate.py
      surface_snapshot.py
      surface_timeavg.py
    data/
      CHAOS-7.7.mat
  LICENSE.md
  MANIFEST.in
  README.md
  setup.py
```

## Examples

### Diagnostics

Example scripts can be found within the module `paropy`.

1. Open `paropy/scripts/diagnostics.py`

2. Set path to simulation data by setting

```
run_ID = <run_ID> # PARODY simulation tag
directory = <path_to_data>
saveDir = <path_to_savefigs>
```

3. Run `paropy/scripts/diagnostics.py`

4. Admire the output:

![](https://raw.githubusercontent.com/jnywong/paropy/master/docs/diag1_test.png)

![](https://raw.githubusercontent.com/jnywong/paropy/master/docs/diag2_test.png)

### Meridional snapshots

1. Open `paropy/scripts/meridional_snapshot.py`

2. Set path to simulation data and snapshot time by setting

```
run_ID, timestamp = 'c-200a', '16.84707134'
directory = <path_to_data>
saveDir = <path_to_savefigs>
```

3. Run `paropy/scripts/meridional_snapshot.py`

4. Admire the output:

![](https://raw.githubusercontent.com/jnywong/paropy/master/docs/merid_test.png)

### Surface snapshots

1. Open `paropy/scripts/surface_snapshot.py`

2. Set path to simulation data and snapshot time by setting

```
run_ID, timestamp = 'c-200a', '16.84707134'
directory = <path_to_data>
saveDir = <path_to_savefigs>
```

3. Run `paropy/scripts/surface_snapshot.py`

4. Admire the output:

![](https://raw.githubusercontent.com/jnywong/paropy/master/docs/surface_test.png)

## Links

* [PyPI](https://pypi.org/project/paropy/)
* [Anaconda Cloud](https://anaconda.org/jnywong/paropy)

## Data

* CHAOS-7.7 [(Finlay et al. 2020)](http://www.spacecenter.dk/files/magnetic-models/CHAOS-7/CHAOS-7.pdf)

## Authors

* [**Jenny Wong**](https://jnywong.netlify.app/) - *Institut de Physique du Globe de Paris - Institut des Sciences de la Terre*

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Del Duca Foundation
* ERC SEIC

:tada:
