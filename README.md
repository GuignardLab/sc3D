# sc3D

sc3D is a Python library to handle 3D spatial transcriptomic datasets.

It is based on the [anndata](https://anndata.readthedocs.io/en/latest/)/[Scanpy](https://scanpy.readthedocs.io/en/stable/) libraries and allows to read, register arrays and compute 3D differential expression.

## Description of the repository

- data: a folder containing examples for the tissue color and tissue name files

- src: a folder containing the source code

- txt: a folder containing the text describing the method (LaTeX, pdf and docx formats are available)

- README.md: this file

- Test-embryo.ipynb: Basic read and write examples (many different ways of writing)

- Spatial-differential-expression.ipynb: a jupyter notebook with some examples on how to perform the spatial differential expression

- spec-file.txt: File containing all the libraries necessary to run sc3D

## Basic usage

Once installed, the library can be called the following way:

```python
from sc3D import Embryo
```

To import some data:

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

Many other functions are available that can be found used in the two provided jupyter notebooks.

## Installation

The current version of sc3D has only been tested with the Python 3.8  (though it is not working with Python>=3.9 because of the open3d dependency)

We strongly advise to use virtual environments to install this package. For example:

```shell
conda create -n sc3D python==3.8
conda activate sc3D
```

Then, the installation can be done either from this repository:

```shell
pip install .
```

or directly from pip:

(temporary until published)
```shell
pip install -i https://test.pypi.org/simple/ sc3D --extra-index-url https://pypi.org/simple
```
(will become the following when published)
```shell
pip install sc3D
```
