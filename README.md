# scSpatial

scSpatial is a library to handle 3D spatial transcriptomic datasets.

It is based on the [anndata]([anndata - Annotated Data &mdash; anndata 0.7.9.dev5+g62089e2 documentation](https://anndata.readthedocs.io/en/latest/))/[Scanpy]([Scanpy – Single-Cell Analysis in Python &mdash; Scanpy 1.8.2 documentation](https://scanpy.readthedocs.io/en/stable/)) libraries and allows to read, register pucks and compute 3D differential expression.

## Description of the repository

- data: a folder containing examples for the tissue color and tissue name files

- src: a folder containing the source code

- txt: a folder containing the text describing the method (LaTeX, pdf and docx formats are available)

- README.md: this file

- Spatial-differential-expression.ipynb: a jupyter notebook with some examples on how to perform the spatial differential expression

- Test-embryo.ipynb: Basic read and write examples (many different ways of writing)

- spec-file.txt: File containing all the libraries necessary to run scSpatial

## Basic usage

Once installed, the library can be called the following way:

```python
from EmbryoAlignment import Embryo
```
