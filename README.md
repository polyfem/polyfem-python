# polyfem-python
![Build](https://github.com/polyfem/polyfem-python/workflows/Build/badge.svg)


[![Anaconda-Server Badge](https://anaconda.org/conda-forge/polyfempy/badges/latest_release_date.svg)](https://anaconda.org/conda-forge/polyfempy)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/polyfempy/badges/downloads.svg)](https://anaconda.org/conda-forge/polyfempy)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/polyfempy/badges/platforms.svg)](https://anaconda.org/conda-forge/polyfempy)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/polyfempy/badges/installer/conda.svg)](https://conda.anaconda.org/conda-forge)

The Python bindings are in alpha. Expect a lot of API changes and possible bugs. Use at your own peril!

<br/>
To use the Python bindings, clone the current repository and use anaconda to install:

```
conda create -n diffipc python=3.9
conda activate diffipc
conda install numpy scipy conda-forge::cython pytorch::pytorch

# optional
export N_THREADS=16

cd polyfem-python/
pip install . -v
```

For full documentation see [https://polyfem.github.io/](https://polyfem.github.io/).
