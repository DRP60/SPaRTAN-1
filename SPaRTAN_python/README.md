# Run SPaRTAN in Python

## System requirements
The package was implemented in Python. Besides python standard library packages, such as argparse, functools, etc. it requires some other package from PyPI repository, such as pandas, numpy, sklearn, scipy, and matplotlib. In order to improve the running time performance, we converted some computationally intensive python functions into two Cython modules which are platform dependentent. We have built Cython extensions for WIndow, Mac, and Linux system respectively which are .so files for Mac and Linux and .pyd for Windows. You can also build Cython extensions onsite by running setup.py file, please see details [here]#Cython extension compilation.


### installation
Download the reporsitory from https://github.com/osmanbeyoglulab/SPaRTAN

Install python3.7 or later. PySpaRTAN used the following dependencies as well: pandas, numpy, scipy, sklearn, matplotlib. 

You can install python dependencies by the following commands:
```sh
pip install pandas
pip install numpy
pip install scipy
pip install -U scikit-learn
pip install matplotlib
```
Cython is not reqired to install unless the pre-built Cython extensions do not work on your system. 

Our codes have been tested on Linux, Mac, and Windows systems. Please see Prerequisites.xlsx for each version of packages we tested on every operating system.

### Cython extension compilation

In the case of the Cython extensions which are paltform-dependent binary modules are not compatible with your operating system, you need to build the extensions on your own machine. 

First, install Cython by running the command
```sh
pip install "cython>0.21"    #Not sure if it works. this is based on the discussion at https://stackoverflow.com/questions/8795617/how-to-pip-install-a-package-with-min-and-max-version-range
```
Cython requires a C compiler to be present on the system. Please check [here](https://cython.readthedocs.io/en/latest/src/quickstart/install.html) for C compiler installation on various operating system.

After installing Cython and C compiler, navigate to the directory SPaRTAN_python and type in the command:
```sh
python setup.py build_ext --inplace
```
This generates new Cython extension .so files (or .pyd files on Windows) located in folder SPaRTAN_python. The previously downloaded .so and .pyd files are renamed to "*_old.so" and "*_old.pyd" 

## Usage

To run PySPaRTAN module, simply execute the command
```sh
python run_PySPaRTAN.py
```
This genereates the results using default parameters


PySPaRTAN model has the following parameters ....(not complete yet)
```sh
dataset_D: default="Dpbmc"

dataset_P: default="Ppbmc5kn_CD8"

dataset_Y: default="Ypbmc5kn_CD8"

input_dir: default="../data/inputs"

output_dir: default="../data/outputs"

spectrumP: default=0.7

rsL2: default=0.001

lambda: default=0.001

normalization: default="l2",  no normalization if set to ""

fold: default=0 (no cross-validation)
```
....
