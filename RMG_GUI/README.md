# RMG GUI

Originally written by Wenchang Lu at NCSU in Python 2.

In June 2023 updated to Python 3, including some re-writes and fixes, by Jackson Burns at MIT.

## Requirements

`conda` for managing a virtual environment, and then:
 - Python 3.7 (other minor releases _may_ work)
 - `pymol`, installed via `conda` with `conda install -c conda-forge pymol-open-source`
 - `PyQt5`, which will come bundled with the above (but if things go wrong, can be installed via `pip` with `pip install PyQt5`)
 - _optional_: `PyCifRW` for reading and writing crystallographic information files (untested)

## Usage

From within this directory (the presence of the `input` file is needed):

`python RMG_GUI.py`

Navigate to the `Configuration` tab and click `Browse...` to load your molecular coordinates file in an accepted format (`xyz` or `cif`), and then use the rest of the GUI to configure your run options for RMG.

Click `Save` on the top panel to write the input file to the directory selected by `...`.
By default, this will overwrite the provided example `input` file, so I recommend selecting a directory other than `RMG_GUI`.

### Customization

Set the default prefix and suffix for pseudopotentials in `Species.py`.

Set `_NCSMURG_ADDON_PATH` in `IOControl.py`, `Misc.py`, `Setup.py`, and `Mdscf.py` if you have any local add-ons.

_Disclaimer_: I (Jackson) don't know what (if anything) the above does.

## Notes

`RMG_GUI.py` is _very_ picky about spacing in the input file.
It will only read and write `input` files with:
 - no spaces in around `=`
 - multi-line inputs (namely `pseudopotential`, `atomic_orbital_files`, and `atoms`) should be formatted like this:

```shell
pseudopotential="
 C    ../C.pp
"
atomic_orbital_files="
 C    ../C-atom/Wave/wave 
"
atoms="
 C     0.000000000000e+00    0.000000000000e+00    0.000000000000e+00      1
...
 C     1.681856921219e+00    5.045570763658e+00    5.045570763658e+00      1
"
```
