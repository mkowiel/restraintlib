# RestraintLib
Bond length and angle restraints for DNA and RNA oligonucleotides.

# Citing

    M.Kowiel, D.Brzezinski, M. Gilski, M.Jaskolski (2020).
    Conformation-dependent restraints for polynucleotides: The sugar moiety.
    Nucleic Acids Res. 48, 962–973. https://doi.org/10.1093/nar/gkz1122 OpenAccess

    M.Gilski, J.Zhao, M.Kowiel, D.Brzezinski, D.H.Turner, M.Jaskolski (2019).
    Accurate geometrical restraints for Watson–Crick base pairs.
    Acta Cryst. B75, 235-245. https://doi.org/10.1107/S2052520619002002 OpenAccess

    M.Kowiel, D.Brzezinski, M.Jaskolski (2016).
    Conformation-dependent restraints for polynucleotides: I. clustering of the geometry of the phosphodiester group.
    Nucleic Acids Res. 44, 8479–8489. https://doi.org/10.1093/nar/gkw717 OpenAccess

# Dependencies

Scripts were tested with

* Python 3.7
* cctbx==2020.10 [cctbx_project](https://github.com/cctbx/cctbx_project)
* scikit-learn==0.20.3
* numpy==1.15.4

# Installation

Install cctbx with use of miniconda/anacona setup (Option 1) or manually (Option 2).
It seems cctbx installation from pypi does not work at the moment (2020.11.27).

## Option 1: use anaconda setup

1. Clone repository or download sourcecode.

    ```bash
    git clone https://github.com/mkowiel/restraintlib.git
    ```

2. Download Miniconda: [minicnda](https://docs.conda.io/en/latest/miniconda.html#linux-installers).
You can download for example https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

3. Install Miniconda

    ```bash
    sh .\Miniconda3-latest-Linux-x86_64.sh
    ```

4. Create environment (in the example with `Python3.7` and the `cctbx` env name).

    ```bash
    conda create -n cctbx python=3.7
    conda install --force-reinstall -n cctbx -c conda-forge cctbx-base
    conda install --force-reinstall -n cctbx -c conda-forge pytest
    conda install --force-reinstall -n cctbx -c conda-forge "scikit-learn>1.0,<1.1"
    conda install -c conda-forge conda-build conda-verify
    ```

4. Install library

    ```bash
    # make sure you build on base environment
    conda activate base
    conda-build . --python 3.7 -c conda-forge
    conda install --force-reinstall -n cctbx -c local --offline restraintlib
    ```

## Option 2: manual installation

Follow the installation instruction in the [cctbx_project](https://github.com/cctbx/cctbx_project)

1. Download [https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py](https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py) in the directory where the cctbx and its dependencies shall be installed.
2. On Linux or Mac OS execute it: python bootstrap.py (you may want to run it with the --help option first to discover the available options).

Then install the library.

1. Clone or download the source code.
2. Install:

    ```bash
    cctbx.python setup.py install
    ```

# Run

```bash
    conda activate cctbx
    restraintlib <printer_format> <input_file> <output_file>
```

# Test

Execute `pytest` in the main directory

# Generate functional restraints models.

The scripts are compatible with python 2.7, however the models are seriazized
(pickeld) so it might be compatible with python 2.7 or different `scikit-learn` version.
It is not the big problem, it is possible to recreate the modes from raw data
(`restraintlib/data/combined_results.csv`). Install dependencies and execute:

```bash
python retrain_restraintlib.py
```

or when installed manually

```bash
cctbx.python retrain_restraintlib.py
```

