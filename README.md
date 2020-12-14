
# MGREML (Multivariate GREML) `VERSION 0.01`

`mgreml` is a command-line tool for rapid estimation of SNP-based heritability and genetic correlations for (nearly) balanced data on many traits in a single analysis using a genomic-relatedness matrix (GRM). `mgreml` can easily handle estimation of the full genetic correlation matrix for up to 100 traits observed in 20,000 individuals. `mgreml` allows users to specify structural models and test hypotheses regarding nested models (e.g. no genetic correlations). In addition, the tool can handle a considerable amount of fixed-effect covariates and a very minor degree of phenotypic missingness. Finally, `mgreml` has options to return e.g. the full set of factor coefficients and their sampling covariance matrix.

## Installation

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

In this short tutorial we will go over the basic functions of `mgreml`. First go over the steps in Installation.

Now, having cloned the `mgreml` repository, the main directory should contain a subdirectory called `tutorial`. This directory in turn contains several files, including `pheno.txt` and `covar.txt`. Details on how this dataset has been generated using simulation can be found in the python script in `./tutorial/simulate.py`

Let's first inspect the `pheno.txt` file. This file contains data in tab-separated format on ten phenotypes observed in a set of 5,000 individuals. The first few columns and rows of this file look as follows:

|FID | IID | Some pheno 101 | Some pheno 102 | Some pheno 103 | ... |
| --- | --- | --- | --- | --- | --- |
| FID 1 | IID 5001 | -2.0716069263263246 | -3.6759090311303346 | -2.745781585038478 | ... |
| FID 2 | IID 5002 | -1.4715012582461549 | -1.4667468932126324 | -1.4061486757142365 | ... | 
| FID 3 | IID 5003 | -3.3915982459245417 | -7.2773361819872 | -0.5107286051290811 | ... |
| FID 4 | IID 5004 | 3.236244118718508 | 1.3410024876651214 | 2.642661382564801 | ... |
| ... | ... | ... | ... | ... | ... |

For the same set of individuals, you have a binary genomic-relatedness matrix (a.k.a. GRM) e.g. computed using [PLINK](https://www.cog-genomics.org/plink/) or [GCTA](https://cnsgenomics.com/software/gcta/). In this case, the set of binary GRM files comprises `data.grm.bin`, `data.grm.N.bin`, and `data.grm.id`. We refer to this set of binary GRM files by its prefix, i.e. `data`.

The simplest command for running an `mgreml` analysis on this data is as follows:

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt --out ./tutorial/nocovs
```

Upon carrying out this command, `mgreml` will first report the follow command will be carried out:

```
Call:
mgreml \
--grm ./tutorial/data \
--pheno ./tutorial/pheno.txt \
--out ./tutorial/nocovs
```

After a few hundred BFGS iterations, `mgreml` will have finished, and have written e.g. heritability estimates to `./tutorial/nocovs.HSq.out` and genetic correlation estimates to `./tutorial/nocovs.RhoG.out`.

The estimated heritabilities are as follows:

| trait | heritability | standard error |
| --- | --- | --- |
| Some pheno 101 | 0.052012142 | 0.029377156 |
| Some pheno 102 | 0.001615452 | 0.031020903 |
| Some pheno 103 | 0.015868615 | 0.029401385 |
| Some pheno 104 | 0.013322251 | 0.029342585 | 
| Some pheno 105 | 0.037358512 | 0.029760677 |
| Some pheno 106 | 0.171444942 | 0.028955814 |
| Some pheno 107 | 0.004249289 | 0.029571263 |
| Some pheno 108 | 0.001346727 | 0.029682721 |
| Some pheno 109 | 0.016323737 | 0.031897206 |
| Some pheno 110 | 0.156523063 | 0.029363308 | 

Comparing these estimates to the true values in `./tutorial/true.HSq.txt`, printed below, we see that our estimates seem to be biased downwards.

| trait | heritability |
| --- | --- |
| Some pheno 101 | 0.47230192893158857 |
| Some pheno 102 | 0.0070794955805743845 |
| Some pheno 103 | 0.034581822843263034 |
| Some pheno 104 | 0.055009566001384846 |
| Some pheno 105 | 0.4667172940834852 |
| Some pheno 106 | 0.4999339124049756 |
| Some pheno 107 | 0.02222839567470805 |
| Some pheno 108 | 0.03412011302653774 |
| Some pheno 109 | 0.16608084405020662 |
| Some pheno 110 | 0.7421233484196411 |

The reason for this bias, is that we did not control for our fixed-effect covariates, in `./tutorial/covar.txt`, which affect the traits of interest. So we need to use the `--covar` option to specify our fixed-effect covariates. This boils down to the following `mgreml` command:

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt --out ./tutorial/covs
```

When we compare the new estimates of heritability (see below) to the true values, we notice any evidence of bias in our estimates is gone.

| trait | heritability | standard error |
| --- | --- | --- |
| Some pheno 101 | 0.449949558 | 0.025379594 |
| Some pheno 102 | 0.014902772 | 0.029610792 |
| Some pheno 103 | 0.031537671 | 0.028705274 |
| Some pheno 104 | 0.043928235 | 0.029429108 |
| Some pheno 105 | 0.481418117 | 0.024593709 |
| Some pheno 106 | 0.487289667 | 0.024525843 |
| Some pheno 107 | 0.036227795 | 0.028984128 |
| Some pheno 108 | 0.016881014 | 0.028562744 |
| Some pheno 109 | 0.165447283 | 0.028834262 |
| Some pheno 110 | 0.742442814 | 0.019310151 |

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

## Support

Before contacting us, please try the following:

1. Go over the tutorial in this `README.md`; this tutorial is quite self-contained
2. The method is described in detail in the supplementary information of the paper (citation below)

If that doesn't work, you can get in touch with us via the [google group](...).

## Citation

If you use the software, please cite

[]()

## License

This project is licensed under GNU GPL v3.

## Authors

Ronald de Vlaming (Vrije Universiteit Amsterdam)

Eric Slob (University of Cambridge)


