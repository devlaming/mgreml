
# MGREML (Multivariate Genomic-relatedness matrix REstricted Maximum Likelihood) `VERSION 0.01`

`mgreml` is a command-line tool for rapid estimation of SNP-based heritability and genetic correlations for many traits at once using a genomic-relatedness matrix (GRM).

## Getting Started

In order to download `mgreml`, please clone this repository using the following commands
```  
git clone https://github.com/devlaming/mgreml.git
cd mgreml
```

In order to install the Python dependencies, you will need the [Anaconda](https://store.continuum.io/cshop/anaconda/) Python distribution and package manager. After installing Anaconda, run the following commands to create an environment with `mgreml`'s dependencies:

```
conda env create --file environment.yml
source activate mgreml
```

Once the above has completed, you can run:

```
./mgreml -h
```
to print a list of all command-line options. If these commands fail with an error, something has gone wrong during the installation process. 

