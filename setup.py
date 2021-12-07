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
    name='scSpatial',
    author='Leo Guignard',
    author_email='leo.guignard@univ-amu.fr',
    version='0.1.0',
    description='Puck alignment and 3D differential expression for 3D sc omics',
    long_description=long_description,
    url='https://github.com/leoguignard/scSpatial',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
    ],
    packages=['scSpatial'],
    package_dir= { '' : 'src' },

    install_requires=['scipy', 'numpy', 'matplotlib', 'pandas',
                      'seaborn', 'scikit-learn', 'open3d',
                      'anndata', 'transformations'],
)