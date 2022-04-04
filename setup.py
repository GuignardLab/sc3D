from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='sc-3D',
    author='Leo Guignard',
    author_email='leo.guignard@univ-amu.fr',
    version='0.1.3',
    description='Array alignment and 3D differential expression for 3D sc omics',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/leoguignard/sc3D',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10'
    ],
    packages=['sc3D'],
    package_dir={ '' : 'src' },
    scripts=['src/sc3D-visualiser'],

    install_requires=['scipy', 'numpy', 'matplotlib', 'pandas',
                      'seaborn', 'scikit-learn', 'anndata',
                      'transformations', 'napari', 'PyQt5'],
)
