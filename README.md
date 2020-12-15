
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

| FID | IID | Some pheno 101 | Some pheno 102 | ...  | Some pheno 109 | Some pheno 110 |
| --- | --- | --- | --- | --- | --- | --- |
| FID 1 | IID 5001 | -2.072 | -3.676 | ...  | 4.641 | 7.931 |
| FID 2 | IID 5002 | -1.472 | -1.467 | ...  | 6.098 | 3.570 |
| FID 3 | IID 5003 | -3.392 | -7.277 | ...  | -0.832 | -5.750 |
| ...  | ...  | ...  | ...  | ...  | ...  | ...  |
| FID 4998 | IID 9998 | 2.575 | 2.740 | ...  | 3.328 | -6.982 |
| FID 4999 | IID 9999 | -3.072 | -0.306 | ...  | 2.530 | -1.255 |
| FID 5000 | IID 10000 | -4.220 | 1.117 | ...  | 2.806 | 3.159 |

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

Comparing these estimates to the true values in `./tutorial/true.HSq.txt`, printed below, we see that our estimates seem to be strongly downwards biased.

| trait | heritability |
| --- | --- |
| Some pheno 101 | 0.472 |
| Some pheno 102 | 0.007 |
| Some pheno 103 | 0.035 |
| Some pheno 104 | 0.055 |
| Some pheno 105 | 0.467 |
| Some pheno 106 | 0.500 |
| Some pheno 107 | 0.022 |
| Some pheno 108 | 0.034 |
| Some pheno 109 | 0.166 |
| Some pheno 110 | 0.742 |

The reason for this bias, is that we did not control for our fixed-effect covariates, in `./tutorial/covar.txt`, which affect the traits of interest. So we need to use the `--covar` option to specify our fixed-effect covariates. This boils down to the following `mgreml` command:

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt --out ./tutorial/covs
```

If we compare the new estimates of heritability (see below) to the true values, taking the standard errors of the estimates into account, we see the massive downwards bias is gone.

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

Importantly, the file with covariates should **NEVER** contain principal components (PCs) from your genetic data. `mgreml` removes the effects of population stratification in the so-called canonical transformation. In essence, within this transformation `mgreml` automatically corrects for the 20 leading PCs your genetic data.

In case you want to change the number of PCs you control for, do **NOT** add them manually your file with covariate data. Instead, use the `--ignore-pcs` option, followed by the total number of leading PCs you want to control for. E.g. `--ignore-pcs 20` is equivalent to the default setting, `--ignore-pcs 40` controls for the 40 leadings PCs, and `--ignore-pcs 0` controls for no PCs at all (not recommended).

For advanced users, the `--ignore-pcs` option can also be followed by a second number, indicating the number of trailing eigenvectors from your GRM to ignore. E.g. `--ignore-pcs 100 1000` controls for 100 leading eigenvectors from your GRM and 1000 trailing eigenvectors. By default no trailing eigenvectors are ignored. However, if the trailing eigenvalues are sufficiently small, a considerable number of trailing eigenvectors may be ignored, boosting CPU time without diminishing statistical efficiency of your analysis too much.

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

Similarly, the estimate of environment correlations in `custom_model.RhoE.out`, reveal these are all estimated at zero, also as expected under this model.

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

By default, `mgreml` will not store any intermediate results. However, using the `--store-iter` option, users can specify every how many iterations they want the current parameter estimates to be stored. E.g. `--store-iter 10` will cause `mgreml` to store estimates every ten iterations. The estimates will be stored in a so-called `.pkl` with a prefix a set by the `--out` option. This `.pkl` file contains the model specification as well as the estimates of that model in a given iteration.

Such a `.pkl` file can also be used to reinitialise `mgreml` e.g. if you accidentally switched off your computer halfway through an analysis. For instance

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt \
                --store-iter 50 \
                --out ./tutorial/covar
```

causes `mgreml` to store results every 50 iterations. Then, if the preceding analysis has reached e.g. just up until iteration 350 before a power outage, we could reinitialise later on using the following command:

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt \
                --reinitialise ./tutorial/covar.estimates.iter.350.bfgs.pkl \
                --out ./tutorial/covar_reinitialised
```

Notice that as such `.pkl` files already implicitly contain the full model specification, the option `--reinitialise` cannot be combined with options such as `--genetic-model`, `--rho-environment` and so on.

In case `--store-iter` is used when estimating a nested versus alternative model (i.e. in conjunction with one of the `--restricted-...` options), `--store-iter` stores two sets of `.pkl` files, namely one set with filenames containing `.estimates.` (denoting the alternative model) and the other containing `.estimates0.` (denoting the nested model).

`.pkl` files can also be used to reinitialise a restricted model, using the `--restricted-reinitialise` option. E.g. the command

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt \
                --restricted-rho-genetic 1 \
                --restricted-rho-environment 0 \
                --store-iter 10 \
                --out ./tutorial/restricted_rhoG1_rhoE0
```
causes two sets of `.pkl` files to be stored (i.e. a file for every 10 iterations) and
```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt \
                --reinitialise ./tutorial/restricted_rhoG1_rhoE0.estimates.iter.350.bfgs.pkl \
                --restricted-reinitialise ./tutorial/restricted_rhoG1_rhoE0.estimates0.iter.30.bfgs.pkl \
                --out ./tutorial/restricted_rhoG1_rhoE0_reinitialised
```
reinitialises estimation for the null and alternative model from appropriate `.pkl` files. Notice that analogous to `--reinitialise`, the `--restricted-reinitialise` option cannot be combined with options such as `--restricted-environment-model` and `--restricted-rho-genetic`, as the `.pkl` already contain the full model specification.

`mgreml` performs basic data management, e.g. in terms of figuring out for which individuals we have phenotypic as well as GRM data (and data on covariates, if applicable). In case `--covar-model` is used `mgreml` also tests if there are any covariates that affect no phenotype at all, and if so, excludes such covariates.

`mgreml` also performs relatedness pruning when using the `--rel-cutoff` option. E.g. `--rel-cutoff 0.025` selects a subset of individuals such that relatedness in the GRM is in excess of 0.025. `mgreml` follows the same algorithm for such pruning as [PLINK](https://www.cog-genomics.org/plink/). Importantly, `mgreml` does this pruning at such a stage, that sample size is maximised (e.g. for a pair of individuals with a relatedness in excess of the threshold, we do not drop the individual for whom we have phenotype data and keep the individual for whom we do not have any phenotype data at all).

In general, `mgreml` simply tries to maximise sample size at each turn. E.g. if an individual has missing values only for a subset of the phenotypes, `mgreml` keeps that individual in the data, by introducing phenotype-by-individual-specific dummies (i.e. dummies that control for individual *i* having a missing value for trait *t*). Even when a covariate is missing, sometimes parts of that observation can still be salvaged (i.e. if the missing covariate does not affect all phenotypes according to `--covar-model`).

However, introducing these dummies to control for gaps in the data can become computationally highly demanding. Controlling for fixed effect covariates has a computational complexity of the order *NT* <sup>2</sup> provided the number of unique covariates is of the order 1. However, if e.g. missingness in each trait is a proportion of sample size, then the total number of unique covariates to control for this missingness becomes of the order *NT*, and thereby the computational complexity of controling for this missingness of the order *N* <sup>2</sup> *T* <sup>3</sup>, which is prohibitively complex for large *N* and *T*.

Therefore, `mgreml` has a `--drop-missings` option, whereby all individuals are dropped that have at least one missing phenotype and/or at least one missing covariate that is relevant (either because `--covar-model` has not been used, or because the file following `--covar-model` indicates the covariate with a missing value for a given individual affects at least one trait).

Finally, `mgreml` has a few advanced option regarding the estimation algorithm. First, `--newton` forces `mgreml` to use a Newton algorithm for solving the optimisation problem instead of BFGS. Although in theory this approach requires fewer iterations that BFGS, we strongly recommend sticking to BFGS: BFGS iterations are faster and BFGS is numerically much more stable, especially for large *T*.

In addition, `mgreml` deems the model to have converged if the length divided by the number of traits, of the gradient vector of the log-likelihoods per observation, is below 10<sup>-5</sup>. The option `--grad-tol` can be used to specify a different treshold. We do **NOT** recommend deviating from 10<sup>-5</sup> by more than one order of magnitude. E.g. you could use `--grad-tol 5E-5` or `--grad-tol 1e-6`. However, we de **NOT** recommend e.g. `--grad-tol 1E-9`, as such a threshold requires a degree of convergence that is beyond numerical precision, nor do we commend e.g. `--grad-tol 0.01`, as this is simply too lenient; `mgreml` simply has not converged when that convergence condition is met.

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


