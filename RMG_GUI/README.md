# RMG GUI

Originally written by Wenchang Lu at NCSU in Python 2, updated by Jackson Burns for Python 3 in 2023.

## Requirements

 - Python >= 3.7
 - `pymol`, installed via `conda` with `conda install -c conda-forge pymol-open-source`
 - `PyQt5`, installed via `pip` with `pip install PyQt5`
 - _optional_: `PyCifRW` for reading and writing crystallographic information files

## Usage

From within this directory:

`python RMG_GUI.py`

Set `_NCSMURG_ADDON_PATH` in `IOControl.py` if you have any local add-ons.
