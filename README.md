
[![License MIT](https://img.shields.io/pypi/l/sc-3D.svg?color=green)](https://github.com/GuignardLab/sc3D/raw/main/LICENSE)
[![PyPI](https://img.shields.io/pypi/v/sc-3D.svg?color=green)](https://pypi.org/project/sc-3D)
[![Python Version](https://img.shields.io/pypi/pyversions/sc-3D.svg?color=green)](https://python.org)
[![tests](https://github.com/GuignardLab/sc3D/workflows/tests/badge.svg)](https://github.com/GuignardLab/sc3D/actions)
[![codecov](https://codecov.io/gh/GuignardLab/sc3D/branch/main/graph/badge.svg)](https://codecov.io/gh/GuignardLab/sc3D)

# sc3D

__sc3D__ is a Python library to handle 3D spatial transcriptomic datasets.

__To access the 3D viewer for sc3D datasets, you can go to the following link: [GuignardLab/napari-sc3D-viewer](https://github.com/GuignardLab/napari-sc3D-viewer)__

You can find it on the Guignard Lab GitHub page: [GuignardLab/sc3D](https://github.com/GuignardLab/sc3D). In there you will be able to find jupyter notebooks giving examples about how to use the datasets.

This code was developed in the context of the following study:

[**Spatial transcriptomic maps of whole mouse embryos.**](https://www.nature.com/articles/s41588-023-01435-6) *Abhishek Sampath Kumar, Luyi Tian, Adriano Bolondi et al.*

The __sc3D__ code is based on the [anndata](https://anndata.readthedocs.io/en/latest/) and [Scanpy](https://scanpy.readthedocs.io/en/stable/) libraries and allows to read, register arrays and compute 3D differential expression.

The dataset necessary to run the tests and look at the results can be downloaded [there](https://figshare.com/s/9c73df7fd39e3ca5422d) for the unregistered dataset (and test the provided algorithms) and [there](https://figshare.com/s/1c29d867bc8b90d754d2) for the already registered atlas to browse with our visualiser. You can find the visualiser [there](https://www.github.com/guignardlab/napari-sc3d-viewer).

## Description of the GitHub repository

- data: a folder containing examples for the tissue color and tissue name files

- src: a folder containing the source code

- txt: a folder containing the text describing the method (LaTeX, pdf and docx formats are available)

- README.md: this file

- notebooks/Test-embryo.ipynb: Basic read and write examples (many different ways of writing)

- notebooks/Spatial-differential-expression.ipynb: a jupyter notebook with some examples on how to perform the spatial differential expression

- notebooks/Embryo-registration.ipynb: a jupyter notebook with an example on how to do the array registration

- setup.py: Setup file to install the library

## Installation

We strongly advise to use virtual environments to install this package. For example using conda or miniconda:

```shell
conda create -n sc-3D
conda activate sc-3D
```

If necessary, install `pip`:

```shell
conda install pip
```

Then, to install the latest stable version:

```shell
pip install sc-3D
```

or to install the latest version from the GitHub repository:

```shell
git clone https://github.com/GuignardLab/sc3D.git
cd sc3D
pip install .
```

### Troubleshooting for latest M1 MacOs chips

If working with an M1 chip, it is possible that all the necessary libraries are not yet available from the usual channels.

To overcome this issue we recommand to manually install the latest, GitHub version of __sc3D__ using [miniforge](https://github.com/conda-forge/miniforge) instead of anaconda or miniconda.

Once miniforge is installed and working, you can run the following commands:

```shell
conda create -n sc-3D
conda activate sc-3D
```

to create your environment, then:

```shell
git clone https://github.com/GuignardLab/sc3D.git
cd sc3D
conda install pip scipy numpy matplotlib pandas seaborn anndata napari
pip install .
```

If the previous commands are still not working, it is possible that you need to install the `pkg-config` package. You can find some information on how to do it there: [install pkg-config](https://gist.github.com/jl/9e5ebbc9ccf44f3c804e)

## Basic usage

Once installed, the library can be called the following way:

```python
from sc3D import Embryo
```

To import some data:

**Note: at the time being, the following conventions are expected:**

- **the x-y coordinates are stored in `data.obsm['X_spatial']`**
- **the array number should be stored in `data.obs['orig.ident']` in the format `".*_[0-9]*"` where the digits after the underscore (`_`) are the id of the array**
- **the tissue type has to be stored in `data.obs['predicted.id']`**
- **the gene names have to be stored as indices or in `data.var['feature_name']`**

Since version `0.1.2`, one can specify the name of the columns where the different necessary informations are stored using the following parameters:

- `tissue_id` to inform the tissue id column
- `array_id` to inform the array/puck/slice id column
- `pos_id` to inform the position column (an `x, y` position is expected within this given column)
- `gene_name_id` to inform the gene name column
- `pos_reg_id` when to inform the registered position column (an `x, y, z` position is expected within this given column)

```python
# To read the data
embryo = Embryo('path/to/data.h5ad')

# To remove potential spatial outliers
embryo.removing_spatial_outliers(th=outlier_threshold)

# To register the arrays and compute the
# spline interpolations
embryo.reconstruct_intermediate(embryo, th_d=th_d,
                                genes=genes_of_interest)

# To save the dataset as a registered dataset (to then look at it in the 3D visualizer)
embryo.save_anndata('path/to/out/registered.h5ad')

# To compute the 3D differential expression for selected tissues
tissues_to_process = [5, 10, 12, 18, 21, 24, 30, 31, 33, 34, 39]
th_vol = .025
_ = embryo.get_3D_differential_expression(tissues_to_process, th_vol)
```

The dataset used for the project this code is from can be downloaded [there](https://cellxgene.cziscience.com/collections/d74b6979-efba-47cd-990a-9d80ccf29055/private) (under the name `mouse_embryo_E8.5_merged_data`)

Many other functions are available that can be found used in the two provided jupyter notebooks.

## Running the notebooks

Two example notebooks are provided.
To run them one wants to first install the jupyter notebook:

```shell
conda install jupyter
```

or

```shell
pip install jupyter
```

The notebooks can the be started from a terminal in the folder containing the `.ipynb` files with the following command:

```shell
jupyter notebook
```

The notebooks should be self content.

Note that the test dataset is not included in this repository put can be downloaded from [there](https://cellxgene.cziscience.com/collections/d74b6979-efba-47cd-990a-9d80ccf29055/private).
