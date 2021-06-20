# Run PySPaRTAN

PySPaRTAN is a  python implementation of the SPaRTAN pipeline (Single-cell Proteomic and RNA based Transcription factor Activity Network) which enables users to infer transcription factor (TF) activities and link cell-surface receptors to TFs by exploiting cellular indexing of transcriptomes and epitopes by sequencing (CITE-seq) datasets with cis-regulatory information.

## Data
In this demo, we use a subset of the CITE-seq data for 5k (Nextgen) PBMC obtained from 10X Genomics website data to train SPaRTAN.
To construct the TF – target gene prior matrix, we downloaded a gene set resource containing TF target-gene interactions from DoRothEA. 

## System requirements
Besides python standard library packages (i.e. argparse, functools, etc.), PySPaRTAN requires some other package from PyPI repository, such as pandas, numpy, sklearn, scipy, and matplotlib. In order to improve the running time performance, we converted some computationally intensive python functions into two Cython modules which, are platform dependent. We have built Cython extensions for Window, Mac, and Linux system, which are .pyd for Windows and .so files for Mac and Linux. You can also build Cython extensions onsite by running setup.py file, which is also explained in this tutorial. 


### Installation
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

If the Cython extensions, which are platform-dependent binary modules, are not compatible with your operating system. You need to build the extensions on your own machine. 

First, install Cython by running the command
```sh
pip install "cython>0.21"    #please see https://stackoverflow.com/questions/8795617/how-to-pip-install-a-package-with-min-and-max-version-range
```
Cython requires a C compiler to be present on the system. Please check [here](https://cython.readthedocs.io/en/latest/src/quickstart/install.html) for C compiler installation on various operating system.

After installing Cython and C compiler, navigate to the directory SPaRTAN_python and type in the command:
```sh
python setup.py build_ext --inplace
```
This generates new Cython extension .so files (or .pyd files on Windows) located in folder SPaRTAN_python. The previously downloaded .so and .pyd files are renamed to "*_old.so" and "*_old.pyd" 

## Usage

To run PySPaRTAN module, simply execute the command to generate output using default parameters
```sh
python run_PySPaRTAN.py
```

PySPaRTAN implementation has the following default parameters 
```sh
input_dir: default="../data/inputs"
output_dir: default="../data/outputs"
spectrumP: default=0.7                                              ##Spectrum cut-off points for P. 
rsL2: default=0.001                                                 ##L_2 norm (Ridge) regularization parameter (rsL2 >=0 )
lambda: default=0.001                                               ##L_1 norm (LASSO) regularization parameter (lambda >=0)
normalization: default="l2",  no normalization if set to ""
fold: default=0 (no cross-validation).                              ##number of folds for tuning parameters (e.g. rsL2, lambda and spectrumP)
```
....
```sh
python run_PySPaRTAN.py -h
```
```sh
usage: run_PySPaRTAN.py [-h] [--input_dir INPUT_DIR] [--output_dir OUTPUT_DIR] [--dataset_D DATASET_D]
                        [--dataset_P DATASET_P] [--dataset_Y DATASET_Y] [--spectrumP SPECTRUMP] [--lamda LAMDA]
                        [--rsL2 RSL2] [--normalization NORMALIZATION] [--fold FOLD]

optional arguments:
  -h, --help            show this help message and exit
  --input_dir INPUT_DIR
                        directory of input files
  --output_dir OUTPUT_DIR
                        directory of output files
  --dataset_D DATASET_D
                        name of gene-TF matrix
  --dataset_P DATASET_P
                        name of the dataset P which will be passed in
  --dataset_Y DATASET_Y
                        name of the dataset Y which will be passed in
  --spectrumP SPECTRUMP
                        Dimension reduction coefficient on protein space
  --lamda LAMDA         LASSO regression coefficient
  --rsL2 RSL2           ridge regression coefficient
  --normalization NORMALIZATION
                        type of normalizion performed on matrices, no normalization if set to empty
  --fold FOLD           how many folds for the cross_validation. No cross_validation and using default/specified
                        parameters if set to 0
```
## References
Ma X*, Somasundaram A*, Qi Z, Singh H, Osmanbeyoglu HU, SPaRTAN, a computational framework for linking cell-surface receptors to transcriptional regulators. bioRxiv 2020.12.22.423961: doi:https://doi.org/10.1101/2020.12.22.423961

Pelossof, R., Singh, I., Yang, J.L., Weirauch, M.T., Hughes, T.R. and Leslie, C.S. (2015) Affinity regression predicts the recognition code of nucleic acid-binding proteins. Nat Biotechnol, 33, 1242-1249.

Osmanbeyoglu, H.U., Pelossof, R., Bromberg, J.F. and Leslie, C.S. (2014) Linking signaling pathways to transcriptional programs in breast cancer. Genome Res, 24, 1869-1880.

Osmanbeyoglu, H.U., Toska, E., Chan, C., Baselga, J. and Leslie, C.S. (2017) Pancancer modelling predicts the context-specific impact of somatic mutations on transcriptional programs. Nat Commun, 8, 14249.

Garcia-Alonso, L., Holland, C.H., Ibrahim, M.M., Turei, D. and Saez-Rodriguez, J. (2019) Benchmark and integration of resources for the estimation of human transcription factor activities. Genome Res, 29, 1363-1375.




