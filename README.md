# sc3D

sc3D is a Python library to handle 3D spatial transcriptomic datasets.

You can find it on the GuignardLab GitHub page: [GuignardLab/sc3D](https://github.com/GuignardLab/sc3D). In there you will be able to find jupyter notebooks giving examples about how to use the datasets.

This code was developed in the context of the following study:

**Spatial transcriptomic maps of whole mouse embryos reveal principles of neural tube patterning.** *Abhishek Sampath Kumar, Luyi Tian, Adriano Bolondi et al.*

The whole code is based on the [anndata](https://anndata.readthedocs.io/en/latest/)/[Scanpy](https://scanpy.readthedocs.io/en/stable/) libraries and allows to read, register arrays and compute 3D differential expression.

The dataset necessary to run the tests and look at the results can be downloaded [there](https://cellxgene.cziscience.com/collections/d74b6979-efba-47cd-990a-9d80ccf29055/private) (it is under the name `mouse_embryo_E8.5_merged_data`).

## Description of the GitHub repository

- data: a folder containing examples for the tissue color and tissue name files

- src: a folder containing the source code

- txt: a folder containing the text describing the method (LaTeX, pdf and docx formats are available)

- README.md: this file

- Test-embryo.ipynb: Basic read and write examples (many different ways of writing)

- Spatial-differential-expression.ipynb: a jupyter notebook with some examples on how to perform the spatial differential expression

- spec-file.txt: File containing all the libraries necessary to run sc3D

## Installation

The current version of sc3D has only been tested with the Python 3.8  (though it is not working with Python>=3.9 because of the open3d dependency)

We strongly advise to use virtual environments to install this package. For example:

```shell
conda create -n sc-3D python==3.8
conda activate sc-3D
```

If necessary, install `pip`:
```shell
conda install pip
```

Then, the installation can be done directly from pip:
```shell
pip install sc-3D
```

or from the GitHub repository:

```shell
pip install .
```

## Basic usage

Once installed, the library can be called the following way:

```python
from sc3D import Embryo
```

To import some data:

**Note: at the time being, the following conventions are expected:**
- **the x-y coordinates are stored in `data.obsm['X_spatial']`**
- **the array number should be stored in `data.obs['orig.ident']` in the format `"{digits}_{array_id:digit}"`**
- **the tissue type has to be stored in `data.obs['predicted.id']`**
- **the gene names have to be stored as indices or in `data.var['feature_name']`**

```python
# To read the data
embryo = Embryo('path/to/data.h5ad')

# To remove potential spatial outliers
embryo.removing_spatial_outliers(th=outlier_threshold)

# To register the arrays and compute the
# spline interpolations
embryo.reconstruct_intermediate(embryo, th_d=th_d,
                                genes=genes_of_interest)

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

Note that the test dataset is not included in this repository put can be downloaded from there:
