
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

If we compare the new estimates of heritability (see below) to the true values, taking the standard errors of the estimates into account, we see that any evidence of bias in our estimates is gone.

|  | heritability | standard error |
| --- | --- | --- |
| Some pheno 101 | 0.450 | 0.025 |
| Some pheno 102 | 0.015 | 0.030 |
| Some pheno 103 | 0.032 | 0.029 |
| Some pheno 104 | 0.044 | 0.029 |
| Some pheno 105 | 0.481 | 0.025 |
| Some pheno 106 | 0.487 | 0.025 |
| Some pheno 107 | 0.036 | 0.029 |
| Some pheno 108 | 0.017 | 0.029 |
| Some pheno 109 | 0.165 | 0.029 |
| Some pheno 110 | 0.742 | 0.019 |

In addition to reporting the heritabilities and their standard errors, `mgreml` also automatically reports genetic and environment correlations, as well as their standard errors.

In case you do not care about standard errors, you can use the `--no-se` option. Especially for a large number of traits, computing the standard errors is computationally demanding, as this requires calculating the average information matrix, which has a computational complexity of the order *NT* <sup>4</sup>, where *T* denotes the number of traits and *N* the number of observations.

`mgreml` also automatically reports the fixed-effect estimates (a.k.a. GLS estimates), including the sampling covariance matrix of those estimates, and their standard errors.

Now, suppose each trait has a different set of covariates, `mgreml` can easily handle this using the `--covar-model` option. This option should be followed by a filename which contains a binary table, indicating which covariate affects which phenotype. E.g. the `tutorial` folder contains `covar_model.txt`, of which the content is shown below:

|  | intercept | my covar 301 | my covar 302 | my covar 303 | my covar 304 | my covar 305 | my covar 306 | my covar 307 | my covar 308 | my covar 309 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Some pheno 101 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Some pheno 102 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Some pheno 103 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Some pheno 104 | 1 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| Some pheno 105 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 |
| Some pheno 106 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 |
| Some pheno 107 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| Some pheno 108 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| Some pheno 109 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 |
| Some pheno 110 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |

Clearly, this file implies that the intercept is a covariate that applies to all phenotypes, whereas all other covariates all affect different traits. We can now perform `mgreml` estimation under this model for the fixed effects using the following command:

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt  --covar-model ./tutorial/covar_model.txt \
                --out ./tutorial/different_covs
```

Similarly, users can also specify which genetic factor affects which trait and which environment factor affects which trait. Such specifications can be passed to `mgreml` using the `--genetic-model` and `--environment-model` options. Note, that any such user-specified structural model must be identified. Moreover, for the factor specification of the environment, `mgreml` requires as many factors as there are traits.

For example, we could impose a factor structure, where there is only one genetic factor, and where there are *T*=10 environment factors, each affecting a different trait. Effectively, this boils down to a model with genetic correlations all equal to one and environment correlations all equal to zero. These factor structures are shown in the files `gen_model.txt` and `env_model.txt` both found in the `tutorial` folder. Both files contain a binary table, with elements equal to one, where a given factor is permitted to affect the given phenotype, and equal to zero otherwise.

To estimate this structural model, we can simply carry out the following command:

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt \
                --genetic-model ./tutorial/gen_model.txt \
                --environment-model ./tutorial/env_model.txt \
                --out ./tutorial/custom_model
```

The estimates in the resulting file, `custom_model.RhoG.out`, reveal that all genetic correlations are estimated at either zero or one, as expected under this model:

|  | Some pheno 101 | Some pheno 102 | Some pheno 103 | ...  | Some pheno 108 | Some pheno 109 | Some pheno 110 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Some pheno 101 | 1 | -1 | 1 | ...  | -1 | -1 | -1 |
| Some pheno 102 | -1 | 1 | -1 | ...  | 1 | 1 | 1 |
| Some pheno 103 | 1 | -1 | 1 | ...  | -1 | -1 | -1 |
| ...  | ...  | ...  | ...  | ...  | ...  | ...  | ... |
| Some pheno 108 | -1 | 1 | -1 | ...  | 1 | 1 | 1 |
| Some pheno 109 | -1 | 1 | -1 | ...  | 1 | 1 | 1 |
| Some pheno 110 | -1 | 1 | -1 | ...  | 1 | 1 | 1 |

Similarly, the estimate of environment correlations in `custom_model.RhoE`, reveal these are all estimated at zero, also as expected under this model.

Notice that in `mgreml`, specifying `--genetic-model` does not require you to also specify `--environment-model` (nor the other way around).

For the specific cases of genetic correlations all equal to one or all equal to zero, and environment correlations all equal to zero, `mgreml` has two custom options that can be used for such cases instead of `--genetic-model` and `--environment-model`, namely `--rho-genetic 0` or  `--rho-genetic 1` and  `--rho-environment 0`.

So, effectively, we could have also estimated the last model using the following command:

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt \
                --rho-genetic 1 \
                --rho-environment 0 \
                --out ./tutorial/rhoG1_rhoE0
```

Inspection of the log-likelihoods in `custom_model.loglik.out` and `rhoG1_rhoE0.loglik.out` indeed reveal that these models yield an identical fit to the data:

```
Log-likelihood of model = -76460.81732177232,
based on data on 10 traits and 4980 observations,
with a model consisting of 1 genetic factors and 10 environment factors,
comprising 10 free genetic factor coefficients and 10 free environment factor coefficients in turn.
Estimates converged after 37 BFGS iterations 
```

In case you have estimated a model, either according to some structural model e.g. using `--genetic-model`, or just the saturated model we started with, you can make `mgreml` report the factor coefficients (i.e. the effect of each factor on each trait) by using the `--all-coefficients` option. Using this option not only reports the estimated factor coefficients, but also the sampling covariance matrix of those estimates. This sampling covariance matrix may grow very large for large *T*.

E.g. the command

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt \
                --all-coefficients \
                --out ./tutorial/full
```

generates, amongst others, the file `full.coeff.out`, which contains 110 estimated factor coefficients in this case, of which a few lines are shown below:

| trait | factor | coefficient |
| --- | --- | --- |
| Some pheno 101 | genetic factor 0 | 0.993 |
| Some pheno 102 | genetic factor 0 | 0.081 |
| Some pheno 103 | genetic factor 0 | 0.203 |
| ... | ... | ... |
Some pheno 109 | environment factor 8 | 0.226 |
Some pheno 110 | environment factor 8 | -0.415 |
Some pheno 110 | environment factor 9 | 0.361 |

The file `full.coeff.var.out` contains a 110-by-110 matrix representing the sampling covariance matrix of those estimates. 

`mgreml` can also be used to specify two models at once, to compare them using a likelihood-ratio test, provided the null model is nested with respect to the alternative. E.g. one can use the following command to compare the saturated model to the previously considered model assuming perfect genetic correlations and no environment correlations at all:

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt \
                --restricted-rho-genetic 1 \
                --restricted-rho-environment 0 \
                --out ./tutorial/restricted_rhoG1_rhoE0
```

Inspection of `restricted_rhoG1_rhoE0.loglik.out` reveals that the saturated model fits the data significantly better than this restricted model:

```
Log-likelihood of nested model (null hypothesis) = -76460.81732177232,
based on data on 10 traits and 4980 observations,
with a model consisting of 1 genetic factors and 10 environment factors,
comprising 10 free genetic factor coefficients and 10 free environment factor coefficients in turn.
Estimates converged after 37 BFGS iterations 

Log-likelihood of parent model (alternative hypothesis) = -66849.63370313856,
based on data on 10 traits and 4980 observations,
with a model consisting of 10 genetic factors and 10 environment factors,
comprising 55 free genetic factor coefficients and 55 free environment factor coefficients in turn.
Estimates converged after 351 BFGS iterations 

Results of likelihood-ratio test with 90 degrees of freedom:
Chi-square test statistic is 19222.36723726752
with P-value = 0.0
```

Notice that `--genetic-model` and `--environment-model` also have their restricted counterparts, i.e. `--restricted-genetic-model` and `--restricted-environment-model`. This means we could have also carried out the preceding comparison of the two models using the following command:

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt \
                --restricted-genetic-model ./tutorial/gen_model.txt \
                --restricted-environment-model ./tutorial/env_model.txt \
                --out ./tutorial/restricted_custom_model
```

Other commands: `--store-iter`, `--reinitialise`, `--restricted-reinitialise`, `--grad-tol`, `--newton`, `--rel-cutoff`, `--drop-missings`, `--ignore-pcs`

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


