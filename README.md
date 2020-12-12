
# MGREML (Multivariate GREML) `VERSION 0.01`

`mgreml` is a command-line tool for rapid estimation of SNP-based heritability and genetic correlations for (nearly) balanced data on many traits at once using a genomic-relatedness matrix (GRM). `mgreml` can easily handle estimation of the full genetic correlation matrix for up to 100 traits observed in 20,000 individuals. `mgreml` allows users to specify structural models and test hypothesis regarding nested models (e.g. no genetic correlations). In addition, the tool can handle a considerable amount of fixed-effect covariates and a very minor degree of phenotypic missingness. Finally, `mgreml` has options to return e.g. the full set of factor coefficients and their sampling covariance matrix.

## Getting Started

In order to download `mgreml`, please clone this repository using the following commands
```  
git clone https://github.com/devlaming/mgreml.git
cd mgreml
```

In order to install the Python dependencies, you will need the [Anaconda](https://www.anaconda.com/) Python distribution and package manager. After installing Anaconda, run the following commands to create an environment with `mgreml`'s dependencies:

```
conda env create --file mgreml.yml
conda activate mgreml
```

Once the above has completed, you can run:

```
python ./mgreml -h
```
to print a list of all command-line options. If these commands fail something has gone wrong during the installation process.

## Tutorial

A short tutorial describing the basic functions of `mgreml` will be described here in due course.

## Updating `mgreml`

You can update to the newest version of `mgreml` using `git`. First, navigate to your `mgreml` directory (e.g. `cd mgreml`), then run
```
git pull
```
If `mgreml` is up to date, you will see 
```
Already up to date.
```
otherwise, you will see `git` output similar to 
```
remote: Enumerating objects: 5, done.
remote: Counting objects: 100% (5/5), done.
remote: Compressing objects: 100% (1/1), done.
remote: Total 3 (delta 2), reused 3 (delta 2), pack-reused 0
Unpacking objects: 100% (3/3), 304 bytes | 13.00 KiB/s, done.
From https://github.com/devlaming/mgreml
   5b4ca9a..b18a8cc  master     -> origin/master
Updating 5b4ca9a..b18a8cc
Fast-forward
 README.md | 6 +++---
 1 file changed, 3 insertions(+), 3 deletions(-)
 ```
which tells you which files were changed.

If you have modified the `mgreml` source code yourself, `git pull` may fail with an error such as `Your local changes to the following files would be overwritten by merge`. 

In case the Python dependencies have changed, you can update the `mgreml` environment with

```
conda env update --file mgreml.yml
```



