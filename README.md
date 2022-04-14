# MGREML (Multivariate GREML) `beta v1.0.1`

[![DOI](https://zenodo.org/badge/176773372.svg)](https://zenodo.org/badge/latestdoi/176773372)

`mgreml` is a command-line tool using Python 3.x for rapid estimation of SNP-based heritability and genetic correlations for (nearly) balanced data on many traits in a single analysis, using a genomic-relatedness matrix (GRM) derived from SNP data on unrelated individuals.

`mgreml` can easily handle estimation of the full genetic correlation matrix for up to 100 traits observed in 20,000 individuals. `mgreml` allows users to specify structural models and test hypotheses regarding nested models (e.g. no genetic correlations). In addition, the tool can handle a considerable amount of fixed-effect covariates and a very minor degree of phenotypic missingness.

Finally, `mgreml` has built-in options to (i) return the full set of factor coefficients and (ii) variance components, as well as (iii) the complete covariance matrix of those estimates, and (iv) estimates of a genetic mediation model for two traits as discussed by Rietveld et al. (2022).

## Installation

:warning: Before downloading `mgreml`, please make sure [Git](https://git-scm.com/downloads) and [Anaconda](https://www.anaconda.com/) with **Python 3.x** are installed.

In order to download `mgreml`, open a command-line interface by starting [Anaconda Prompt](https://docs.anaconda.com/anaconda/user-guide/getting-started/), navigate to your working directory, and clone the `mgreml` repository using the following command:

```  
git clone https://github.com/devlaming/mgreml.git
```

Now, enter the newly created `mgreml` directory using:

```
cd mgreml
```

Then run the following commands to create a custom Python environment which has all of `mgreml`'s dependencies (i.e. an environment that has packages such as `numpy` and `pandas` pre-installed):

```
conda env create --file mgreml.yml
conda activate mgreml
```

(or `activate mgreml` instead of `conda activate mgreml` on some machines).

In case you cannot create a customised conda environment (e.g. because of insufficient user rights) or simply prefer to use Anaconda Navigator or `pip` to install packages e.g. in your base environment rather than a custom environment, please note that `mgreml` only requires Python 3.x with the packages `networkx`, `numpy`, `pandas`, `psutil`, `scipy`, and `tqdm` installed. To all [SURFsara](https://userinfo.surfsara.nl/) users: please note that the default Python 3.x environment already has all necessary packages [pre-installed](https://userinfo.surfsara.nl/systems/shared/software/python).

Once the above has completed, you can now run

```
python ./mgreml.py -h
```

to print a list of all command-line options. If this command fails, something has gone wrong during installation.

:warning: **Windows users**: in case the preceding command fails, try replacing slashes (i.e. `/`) in all `mgreml` commands by backslashes (i.e. `\`), so e.g. try

```
python .\mgreml.py -h
```

## Tutorial

In this tutorial, you will learn how to use `mgreml`. Before you start using `mgreml`, please go over the steps under [Installation](#installation). The following topics will be covered in this tutorial:

1. [Tutorial data](#tutorial-data)
2. [Basic estimation](#basic-estimation)
3. [Controlling for covariates](#controlling-for-covariates)
4. [Population stratification](#population-stratification)
5. [Standard errors](#standard-errors)
6. [Different traits with different covariates](#different-traits-with-different-covariates)
7. [Specifying structural models](#specifying-structural-models)
8. [Factor coefficients and variance components](#factor-coefficients-and-variance-components)
9. [Nested models and likelihood-ratio tests](#nested-models-and-likelihood-ratio-tests)
10. [Estimation reinitialisation](#estimation-reinitialisation)
11. [Genetic mediation analysis](#genetic-mediation-analysis)
12. [Data formats and management](#data-formats-and-management)
13. [Missing data and unbalancedness](#missing-data-and-unbalancedness)
14. [Advanced options](#advanced-options)

### Tutorial data

Now that you have cloned the `mgreml` repository, and `mgreml` is up-and-running, the main directory of `mgreml` should contain a subdirectory called `tutorial`. This directory in turn contains several files, including `pheno.txt` and `covar.txt`. Details on how this dataset has been generated using simulation can be found in the python script in `./tutorial/simulate.py`.

Let us first inspect the `pheno.txt` file. This file contains data in tab-separated format on ten phenotypes observed in a set of 5,000 individuals. The first two columns list family and individual ID, followed by the phenotypes:

| FID | IID | Some pheno 101 | Some pheno 102 | ... | Some pheno 109 | Some pheno 110 |
| --- | --- | --- | --- | --- | --- | --- |
| FID 1 | IID 5001 | 3.738 | 1.447 | ... | 0.585 | 1.848 |
| FID 2 | IID 5002 | -3.667 | -1.704 | ... | 0.317 | -0.946 |
| FID 3 | IID 5003 | -2.644 | -0.737 | ... | -2.647 | 0.093 |
| ... | ... |  |  | ... |  |  |
| FID 4998 | IID 9998 | -2.487 | 0.550 | ... | 2.467 | -2.093 |
| FID 4999 | IID 9999 | -2.460 | -2.980 | ... | 6.344 | -0.201 |
| FID 5000 | IID 10000 | -2.192 | -3.691 | ... | 0.024 | -0.906 |

Although `mgreml` in principle can handle phenotype data without header (using a modifier that we discuss later on), we recommend always including headers in your data, so e.g. your phenotypes are labelled, allowing `mgreml` output to refer to specific phenotype names rather than ambiguous indices such as `1`, `2`, `3` etc.

For the same set of individuals, you have a binary genomic-relatedness matrix (a.k.a. GRM) e.g. computed using [LDAK](https://dougspeed.com/ldak/), [PLINK](https://www.cog-genomics.org/plink/), or [GCTA](https://cnsgenomics.com/software/gcta/). In this case, the set of binary GRM files comprises `data.grm.bin`, `data.grm.N.bin`, and `data.grm.id` (MGREML ignores whether or not the `.grm.N.bin` file is present). We refer to this set of binary GRM files by its prefix, i.e. `data`.

### Basic estimation

The command for running an `mgreml` analysis on this data without correcting for any covariates at all is as follows:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \ 
                   --no-intercept --out ./tutorial/nocovs
```

:warning: The phenotypic data itself may only be numerical. E.g., values such as `yes` and `no` are not permitted as phenotypic values. Instead, please use values such as `1` and `0`. Please notice that the labels of the phenotypes (i.e., in the header row) do **not** need to be numerical, of course.

When carrying out this command, `mgreml` will first show the welcome screen and directly after that summarise the input options that you specified:

```
Your call:
./mgreml.py \
--grm ./tutorial/data \
--pheno ./tutorial/pheno.txt \
--no-intercept \
--out ./tutorial/nocovs
```

After a few hundred BFGS iterations, `mgreml` will have finished, and have written e.g. heritability estimates to `./tutorial/nocovs.HSq.out` and genetic correlation estimates to `./tutorial/nocovs.RhoG.out`.

The estimated heritabilities are as follows:

| trait | heritability | standard error |
| --- | --- | --- |
| Some pheno 101 | 0.088 | 0.016 |
| Some pheno 102 | 0.063 | 0.016 |
| Some pheno 103 | 0.088 | 0.016 |
| Some pheno 104 | 0.064 | 0.016 |
| Some pheno 105 | 0.059 | 0.016 |
| Some pheno 106 | 0.063 | 0.016 |
| Some pheno 107 | 0.085 | 0.016 |
| Some pheno 108 | 0.145 | 0.017 |
| Some pheno 109 | 0.037 | 0.015 |
| Some pheno 110 | 0.154 | 0.017 |

Comparing these estimates to the true values in `./tutorial/true.HSq.txt`, printed below, we see that our estimates seem to be biased.

| trait | heritability |
| --- | --- |
| Some pheno 101 | 0.249 |
| Some pheno 102 | 0.252 |
| Some pheno 103 | 0.248 |
| Some pheno 104 | 0.249 |
| Some pheno 105 | 0.253 |
| Some pheno 106 | 0.252 |
| Some pheno 107 | 0.248 |
| Some pheno 108 | 0.250 |
| Some pheno 109 | 0.249 |
| Some pheno 110 | 0.249 |

The simple reason for this bias is that we did not control for any fixed-effect covariates.

### Controlling for covariates

By removing the `--no-intercept` option, `mgreml` automatically adds one fixed effect per phenotype, namely a fixed effect that controls for differences in mean across phenotypes:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --out ./tutorial/intercept
```

Resulting SNP heritability estimates in `./tutorial/intercept.HSq.out`, however, show our estimates are still strongly biased:

| trait | heritability | standard error |
| --- | --- | --- |
| Some pheno 101 | 0.084 | 0.016 |
| Some pheno 102 | 0.083 | 0.016 |
| Some pheno 103 | 0.088 | 0.016 |
| Some pheno 104 | 0.076 | 0.016 |
| Some pheno 105 | 0.058 | 0.016 |
| Some pheno 106 | 0.072 | 0.016 |
| Some pheno 107 | 0.088 | 0.016 |
| Some pheno 108 | 0.151 | 0.017 |
| Some pheno 109 | 0.055 | 0.016 |
| Some pheno 110 | 0.159 | 0.017 |

The reasons this bias persists is that more fixed-effect covariates are at play than just the intercept. The file `./tutorial/covar.txt` contains the other covariates that affect the traits of interest. So we need to use the `--covar` option to specify these additional fixed-effect covariates. This boils down to the following `mgreml` command:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt --out ./tutorial/covs
```

:warning: As with the phenotypic data, the actual data on the covariates may also only be numerical. E.g., values such as `female` and `male` are not permitted as values of covariates. Instead, please use values such as `1` and `0`. Please notice that the labels of the covariates (i.e., in the header row) do **not** need to be numerical, of course.

Notice that analyses including covariates are computationally slightly more demanding. E.g. in this case we have 10 covariates (i.e. the intercept + 9 additional covariates in `./tutorial/covar.txt`), each of which is allowed to have a different effect on each trait. As we have 10 traits, this means we have 100 fixed effects in total, which our model needs to take into account.

If we compare the new estimates of heritability (see below) to the true values, taking the standard errors of the estimates into account, we see the strong bias is gone.

| trait | heritability | standard error |
| --- | --- | --- |
| Some pheno 101 | 0.253 | 0.017 |
| Some pheno 102 | 0.260 | 0.017 |
| Some pheno 103 | 0.251 | 0.017 |
| Some pheno 104 | 0.258 | 0.017 |
| Some pheno 105 | 0.236 | 0.017 |
| Some pheno 106 | 0.247 | 0.017 |
| Some pheno 107 | 0.269 | 0.017 |
| Some pheno 108 | 0.246 | 0.017 |
| Some pheno 109 | 0.255 | 0.017 |
| Some pheno 110 | 0.267 | 0.017 |

Moreover, when looking at the heatmaps of the genetic correlations as estimated by `mgreml` and the true genetic correlations, we see that `mgreml` provides highly accurate estimates in this simulation:

![Comparison of `mgreml` estimates of genetic correlations and true genetic correlations](https://github.com/devlaming/mgreml/blob/master/tutorial/rhoG_estimates_true.png?raw=true)

### Population stratification

In the previous part, you may have noticed that `covar.txt` does not include any principal components (PCs) from the genetic data as fixed-effect covariates. Perhaps surprisingly, this is perfectly normal for an `mgreml` analysis. In fact, in `mgreml`, the file with your covariates should **NEVER** contain PCs from your genetic data, as `mgreml` already removes the effects of population stratification during the so-called canonical transformation. By default, `mgreml` removes the effects of 20 leading PCs from your genetic data. The effective sample size is reduced by 20 as a result of this correction for PCs.

In case you want to change the number of PCs you control for, do **NOT** add these PCs to your file with covariate data. Instead, use the `--adjust-pcs` option, followed by the total number of leading PCs you want to control for. E.g. `--adjust-pcs 20` is equivalent to the default setting, `--adjust-pcs 40` controls for the 40 leadings PCs, and `--adjust-pcs 0` controls for no PCs at all (not recommended). In these three cases, the sample size is reduced by 20, 40, and zero respectively.

For instance, the command

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --adjust-pcs 1000 --out ./tutorial/many_pcs
```

makes `mgreml` adjust for 1000 leading PCs from the genetic data.

As there is no population stratification in our data (by virtue of our simulation design), this means adjusting for so many PCs will just reduce precision of our estimates, without eliminating any bias. If we look at `many_pcs.HSq.out` we see that our estimates indeed have considerably higher standard errors:

| trait | heritability | standard error |
| --- | --- | --- |
| Some pheno 101 | 0.230 | 0.032 |
| Some pheno 102 | 0.252 | 0.032 |
| Some pheno 103 | 0.262 | 0.031 |
| Some pheno 104 | 0.298 | 0.029 |
| Some pheno 105 | 0.235 | 0.032 |
| Some pheno 106 | 0.258 | 0.032 |
| Some pheno 107 | 0.274 | 0.031 |
| Some pheno 108 | 0.245 | 0.031 |
| Some pheno 109 | 0.230 | 0.032 |
| Some pheno 110 | 0.311 | 0.030 |

For advanced users, the `--adjust-pcs` option can also be followed by a second number, indicating the number of trailing eigenvectors from your GRM to adjust for. E.g. `--adjust-pcs 100 1000` controls for 100 leading eigenvectors from your GRM and 1000 trailing eigenvectors. Doing this decreases the overall sample size by 100 + 1000 = 1100.

By default no trailing eigenvectors are adjusted for. However, if the trailing eigenvalues are sufficiently small, we may adjust for a considerable number of them, boosting CPU time (as the overall sample size becomes smaller), without diminishing statistical efficiency of the analysis too much (as effectively only less informative bits of data are ignored). Note that adjusting for trailing eigenvectors may require you to use `--no-intercept`, as the last eigenvector from the GRM tends to be highly multicollinear with the intercept for technical reasons.

### Standard errors

In addition to reporting the heritabilities, and genetic and environment correlations, `mgreml` also automatically reports the standard errors of all estimates.

In case you do not wish `mgreml` to compute standard errors, you can use the `--no-se` option. Especially for a large number of traits, computing the standard errors is computationally demanding, as this requires calculating the average information matrix, which has a computational complexity of the order *NT*<sup> 4</sup>, where *T* denotes the number of traits and *N* the number of observations.

`mgreml` also automatically reports the fixed-effect estimates (a.k.a. GLS estimates), including the covariance matrix of those estimates, and their standard errors. If the `--no-se` option is used, the estimated covariance matrix and standard errors of the GLS estimates will not be included either.

### Different traits with different covariates

Now, suppose each trait has a different set of covariates, `mgreml` can easily handle this using the `--covar-model` option. This option should be followed by a filename which contains a binary table, indicating which covariate affects which phenotype. E.g. the `tutorial` folder contains `covar_model.txt`, of which the content is shown below:

|  | my covar 301 | my covar 302 | my covar 303 | my covar 304 | my covar 305 | my covar 306 | my covar 307 | my covar 308 | my covar 309 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Some pheno 101 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Some pheno 102 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Some pheno 103 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Some pheno 104 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| Some pheno 105 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 |
| Some pheno 106 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 |
| Some pheno 107 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| Some pheno 108 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| Some pheno 109 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 |
| Some pheno 110 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |

Clearly, each covariate affects a different trait. We can now perform `mgreml` estimation under this model for the fixed effects using the following command:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt --covar-model ./tutorial/covar_model.txt \
                   --out ./tutorial/different_covs
```

Notice here that we did not use `--no-intercept`. This means `mgreml` (1) adds the intercept to the set of covariates and (2) assumes the intercept applies to all phenotypes.

So in total, we have 10 fixed effects for the intercept (i.e. one fixed effect per trait), and 9 additional fixed effects for `my covar 301` up until `my covar 309`.

Now, `different_covs.GLS.est.out`, in the folder `tutorial`, shows the fixed-effect estimates for the intercept affecting all traits and for the covariates that affect a given trait according to `covar_model.txt`:

| trait | covariate | beta hat | standard error |
| --- | --- | --- | --- |
| Some pheno 101 | intercept | -0.328 | 0.047 |
| Some pheno 102 | intercept | -1.063 | 0.045 |
| Some pheno 102 | my covar 301 | 1.646 | 0.011 |
| Some pheno 103 | intercept | -0.307 | 0.054 |
| ... | ... | ... | ... |
| Some pheno 109 | intercept | 1.184 | 0.109 |
| Some pheno 109 | my covar 308 | 5.766 | 0.013 |
| Some pheno 110 | intercept | 0.684 | 0.035 |
| Some pheno 110 | my covar 309 | 1.398 | 0.005 |

E.g. `my covar 301` does not affect `Some pheno 101` in this case.

### Specifying structural models

Analogous to `--covar-model`, users can also specify which genetic factor affects which trait and which environment factor affects which trait. Such specifications can be passed to `mgreml` using the `--genetic-model` and `--environment-model` options. Note that any such user-specified structural model must be identified. Moreover, for the factor specification of the environment, `mgreml` requires as many factors as there are traits.

For example, we could impose a factor structure, where there is only one genetic factor, and where there are *T*=10 environment factors, each affecting only a single trait, and no trait being affected by two distinct environment factors.

Effectively, this boils down to a model with genetic correlations all equal to one and environment correlations all equal to zero. These factor structures are shown in the files `gen_model.txt` and `env_model.txt` both found in the `tutorial` folder. Both files contain a binary table, with elements equal to one, where a given factor is permitted to affect the given phenotype, and equal to zero otherwise.

To estimate this structural model, we can simply carry out the following command:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --genetic-model ./tutorial/gen_model.txt \
                   --environment-model ./tutorial/env_model.txt \
                   --out ./tutorial/custom_model
```

The estimates in the resulting file, `custom_model.RhoG.out`, reveal that all genetic correlations are estimated at either zero or one, as expected under this model:

|  | Some pheno 101 | Some pheno 102 | Some pheno 103 | ... | Some pheno 108 | Some pheno 109 | Some pheno 110 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Some pheno 101 | 1 | 1 | -1 | ... | 1 | 1 | -1 |
| Some pheno 102 | 1 | 1 | -1 | ... | 1 | 1 | -1 |
| Some pheno 103 | -1 | -1 | 1 | ... | -1 | -1 | 1 |
| ... | ... | ... | ... | ... | ... | ... | ... |
| Some pheno 108 | 1 | 1 | -1 | ... | 1 | 1 | -1 |
| Some pheno 109 | 1 | 1 | -1 | ... | 1 | 1 | -1 |
| Some pheno 110 | -1 | -1 | 1 | ... | -1 | -1 | 1 |

Similarly, the estimates of environment correlations, in `custom_model.RhoE.out`, show these are all estimated at zero, also as expected under this model.

Notice that in `mgreml`, specifying `--genetic-model` does not require you to also specify `--environment-model` (nor the other way around).

For the specific cases of genetic correlations all equal to one or all equal to zero, and environment correlations all equal to zero, `mgreml` has two custom options that can be used for such cases instead of `--genetic-model` and `--environment-model`, namely `--rho-genetic 0` or `--rho-genetic 1` and `--rho-environment 0`.

So, effectively, we could have also estimated the last model using the following command:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --rho-genetic 1 \
                   --rho-environment 0 \
                   --out ./tutorial/rhoG1_rhoE0
```

Inspection of the log-likelihoods in `custom_model.loglik.out` and `rhoG1_rhoE0.loglik.out` indeed reveal that these models yield an identical fit to the data:

```
Log-likelihood of model = -103463.59460195134,
based on data on 10 traits and 4980 observations,
with a model consisting of 1 genetic factors and 10 environment factors,
comprising 10 free genetic factor coefficients and 10 free environment factor coefficients in turn.
Controlled for 100 fixed-effect covariates in total in this model.
Estimates converged after 36 BFGS iterations. 
```

Notice that the option `--rho-genetic` cannot be combined with `--genetic-model` and, similarly, that `--rho-environment` cannot be combined with `--environment-model`.

In addition, for the specific case of no genetic variance at all, `mgreml` also has the custom option `--no-var-genetic`. This enforces genetic variance to be absent for all traits in your data. E.g. the following command

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --no-var-genetic \
                   --out ./tutorial/novarG
```

yields heritability estimates all equal to zero, as expected, in `novarG.HSq.out`.

Notice that the option `--no-var-genetic` cannot be combined with `--rho-genetic` and/or `--genetic-model`.

### Factor coefficients and variance components

In case you estimate a model using `mgreml`, either according to some specific structural model (e.g. using `--genetic-model`) or the default fully saturated model we started with, `mgreml` can report the factor coefficients (i.e. the estimated effect of each factor on each trait) by using the `--factor-coefficients` option. Unless `--no-se` is used, the `--factor-coefficients` option not only reports the estimated factor coefficients, but also the complete covariance matrix of those estimates. :warning: This covariance matrix may grow very large for large *T*.

E.g. the command

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --factor-coefficients \
                   --out ./tutorial/factors
```

generates, amongst others, the file `factors.coeff.out`, which contains 110 estimated factor coefficients in this case, of which a few lines are shown below:

| trait | factor | estimate | standard error |
| --- | --- | --- | --- |
| Some pheno 101 | genetic factor 0 | 1.005 | 0.039 |
| Some pheno 102 | genetic factor 0 | -0.197 | 0.056 |
| Some pheno 103 | genetic factor 0 | 0.285 | 0.056 |
| Some pheno 104 | genetic factor 0 | -0.341 | 0.054 |
| ... | ... | ... | ... |
| Some pheno 110 | environment factor 7 | 0.034 | 0.005 |
| Some pheno 109 | environment factor 8 | 0.929 | 0.013 |
| Some pheno 110 | environment factor 8 | 0.230 | 0.004 |
| Some pheno 110 | environment factor 9 | 0.106 | 0.001 |

The file `factors.coeff.var.out` contains a 110-by-110 matrix representing the covariance matrix of those estimates. 

Similarly, `mgreml` can also return the estimated variance components (again either based on some structural model, or just the saturated model), also including the covariance matrix of those estimated variance components (unless `--no-se` is used). To get these results, use the `--variance-components` option. E.g. the command

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --variance-components \
                   --out ./tutorial/components
```

generates, amongst others, the file `components.VCs.out`, which contains 110 estimated covariance components in this case, of which a few lines are shown below:

| component type | first trait | second trait | estimate | standard error |
| --- | --- | --- | --- | --- |
| genetic covariance | Some pheno 101 | Some pheno 101 | 1.010 | 0.078 |
| genetic covariance | Some pheno 101 | Some pheno 102 | -0.198 | 0.055 |
| genetic covariance | Some pheno 101 | Some pheno 103 | 0.286 | 0.055 |
| genetic covariance | Some pheno 101 | Some pheno 104 | -0.343 | 0.055 |
| ... | ... | ... | ... | ... |
| environment covariance | Some pheno 108 | Some pheno 110 | 0.006 | 0.054 |
| environment covariance | Some pheno 109 | Some pheno 109 | 2.977 | 0.076 |
| environment covariance | Some pheno 109 | Some pheno 110 | 0.985 | 0.057 |
| environment covariance | Some pheno 110 | Some pheno 110 | 2.959 | 0.076 |

The file `components.VCs.var.out` contains a 110-by-110 matrix representing the covariance matrix of those estimates. 

### Nested models and likelihood-ratio tests

`mgreml` can also be used to specify two models at once, to compare them using a likelihood-ratio test, provided the null model is nested with respect to the alternative. E.g. one can use the following command to compare the saturated model to the previously considered model assuming perfect genetic correlations and no environment correlations at all:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --restricted-rho-genetic 1 \
                   --restricted-rho-environment 0 \
                   --out ./tutorial/restricted_rhoG1_rhoE0
```

Inspection of `restricted_rhoG1_rhoE0.loglik.out` reveals that the saturated model fits the data significantly better than this restricted model:

```
Log-likelihood of nested model (null hypothesis) = -103463.59460195134,
based on data on 10 traits and 4980 observations,
with a model consisting of 1 genetic factors and 10 environment factors,
comprising 10 free genetic factor coefficients and 10 free environment factor coefficients in turn.
Controlled for 100 fixed-effect covariates in total in this model.
Estimates converged after 36 BFGS iterations.

Log-likelihood of parent model (alternative hypothesis) = -85227.36770921224,
based on data on 10 traits and 4980 observations,
with a model consisting of 10 genetic factors and 10 environment factors,
comprising 55 free genetic factor coefficients and 55 free environment factor coefficients in turn.
Controlled for 100 fixed-effect covariates in total in this model.
Estimates converged after 53 BFGS iterations.

Results of likelihood-ratio test with 90 degrees of freedom:
Chi-square test statistic is 36472.453785478196
with P-value = 0.0
```

Notice that `--no-var-genetic`, `--genetic-model`, and `--environment-model` also have their restricted counterparts, i.e. `--restricted-no-var-genetic`, `--restricted-genetic-model`, and `--restricted-environment-model`. This means we could have also carried out the preceding comparison of the two models using the following command:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --restricted-genetic-model ./tutorial/gen_model.txt \
                   --restricted-environment-model ./tutorial/env_model.txt \
                   --out ./tutorial/restricted_custom_model
```

As expected, the file `restricted_custom_model.loglik.out` contains results that are identical to those found in `restricted_rhoG1_rhoE0.loglik.out`. 

As a further example, to test if the genetic components even improve the fit of the model in the first place, we could carry out the command

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --restricted-no-var-genetic \
                   --out ./tutorial/restricted_novarG
```

where output file `restricted_novarG.loglik.out` reveals that allowing for genetic variance significantly improves the fit of the model:

```
Log-likelihood of nested model (null hypothesis) = -94524.0137520493,
based on data on 10 traits and 4980 observations,
with a model consisting of 1 genetic factors and 10 environment factors,
comprising 0 free genetic factor coefficients and 55 free environment factor coefficients in turn.
Controlled for 100 fixed-effect covariates in total in this model.
Estimates converged after 20 BFGS iterations.

Log-likelihood of parent model (alternative hypothesis) = -85227.36770921224,
based on data on 10 traits and 4980 observations,
with a model consisting of 10 genetic factors and 10 environment factors,
comprising 55 free genetic factor coefficients and 55 free environment factor coefficients in turn.
Controlled for 100 fixed-effect covariates in total in this model.
Estimates converged after 53 BFGS iterations.

Results of likelihood-ratio test with 55 degrees of freedom:
Chi-square test statistic is 18593.292085674126
with P-value = 0.0
```

As before, `--restricted-no-var-genetic`, `--restricted-rho-genetic`, and/or `--restricted-genetic-model` cannot be combined with one another. Similarly, `--restricted-rho-environment` and `--restricted-environment-model` cannot be combined with each other.

:warning: **when using options such as `--restricted-genetic-model` and `--genetic-model`, it is your own responsibility to ensure the restricted model is nested with respect to the other model.** `mgreml` only inspects nestedness superficially. The best way to allow `mgreml` to assert nestedness is to appropriately label the factors in both models.

### Estimation reinitialisation 

By default, `mgreml` will not store any intermediate results. However, using the `--store-iter` option, users can specify every how many iterations they want the current parameter estimates to be stored. E.g. `--store-iter 10` will make `mgreml` store estimates every ten iterations. The estimates will be stored in a so-called `.pkl` with a prefix a set by the `--out` option. This `.pkl` file contains the model specification as well as the estimates of that model in a given iteration.

Such a `.pkl` file can also be used to reinitialise `mgreml` e.g. if you accidentally switched off your computer halfway through an analysis. For instance

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --store-iter 50 \
                   --out ./tutorial/covar
```

makes `mgreml` store results every 50 iterations. Then, if the preceding analysis has reached e.g. up until iteration 52 before a power outage, we could reinitialise later on using the following command:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --reinitialise ./tutorial/covar.estimates.iter.50.bfgs.pkl \
                   --out ./tutorial/covar_reinitialised
```

Notice that as such `.pkl` files already implicitly contain the full model specification, the option `--reinitialise` cannot be combined with options such as `--genetic-model`, `--rho-environment`, `--no-var-genetic`, and so on.

In case `--store-iter` is used when estimating a nested versus alternative model (i.e. in conjunction with one of the `--restricted-...` options), `--store-iter` stores two sets of `.pkl` files, namely one set with filenames containing `.estimates.` (denoting the alternative model) and the other containing `.estimates0.` (denoting the nested model).

`.pkl` files can also be used to reinitialise a restricted model, using the `--restricted-reinitialise` option. E.g. the command

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --restricted-rho-genetic 1 \
                   --restricted-rho-environment 0 \
                   --store-iter 10 \
                   --out ./tutorial/restricted_rhoG1_rhoE0
```
makes `mgreml` store two sets of `.pkl` files (i.e. a file for every ten iterations, for both the restricted and alternative model) and
```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --reinitialise ./tutorial/restricted_rhoG1_rhoE0.estimates.iter.30.bfgs.pkl \
                   --restricted-reinitialise ./tutorial/restricted_rhoG1_rhoE0.estimates0.iter.30.bfgs.pkl \
                   --out ./tutorial/restricted_rhoG1_rhoE0_reinitialised
```
reinitialises estimation for the null and alternative model from appropriate `.pkl` files. Notice that analogous to `--reinitialise`, the `--restricted-reinitialise` option cannot be combined with options such as `--restricted-environment-model` and `--restricted-rho-genetic`, as the `.pkl` file already contains the full model specification.

### Genetic mediation analysis

Rietveld et al. (2022; see [Citation](#citation)) propose a structural equations model (SEM), which can be used to answer the question to which degree the genetic variance of outcome *Y* is mediated by supposed mediator *M*.

`mgreml` has a `--mediation` option, which estimates the relevant parameters from this SEM by, first, fitting a bivariate saturated model for *M* and *Y* and, second, by transforming the estimated variance components and their sampling variance matrix to estimates and standard errors of the parameters in the SEM.

The `--mediation` option forces `mgreml` to only consider the first two phenotypes in your phenotype file, where the first phenotype is treated as mediator *M* and the second phenotype as outcome *Y*. Any subsequent phenotypes will be ignored.

In case you do not wish standard errors to be reported, you can combine `--mediation` with the `--no-se` option. Please note that `--mediation` cannot be combined with any of the following options: `--(restricted-)genetic-model`, `--(restricted-)rho-genetic`, `--(restricted-)no-var-genetic`, `--(restricted-)environment-model`, `--(restricted-)rho-environment`, and `--(restricted-)reinitialise`.

The resulting SEM estimates will be stored in an output file ending in `.mediation.out`.

As an example, consider `mediation.txt` in the subdirectory `tutorial`. This file comprises two phenotypes, labelled Mediator and Outcome. A few lines from this file are shown below:

| FID | IID | Mediator | Outcome |
| --- | --- | --- | --- |
| FID 1 | IID 5001 | -10.670 | -11.437 |
| FID 2 | IID 5002 | 1.947 | 2.618 |
| FID 3 | IID 5003 | 5.800 | 10.586 |
| ... | ... | ... | ... |
| FID 4998 | IID 9998 | -12.892 | -11.536 |
| FID 4999 | IID 9999 | -2.089 | -3.536 |
| FID 5000 | IID 10000 | -13.177 | -9.323 |

An overview of the SEM that underlies this data is shown in the figure below:

![Structural equations model used to generate phenotypes `mediation.txt`](https://github.com/devlaming/mgreml/blob/development/tutorial/sem.png?raw=true)

This SEM for *M* and *Y* is equivalent to the following two equations: (1) *M* = 3*G* + 4*G*<sup>\*</sup> + 5*E*<sup>\*</sup> and (2) *Y* = *M* + 2*G* + 4*E* = 5*G* + 4*G*<sup>\*</sup> + 5*E*<sup>\*</sup> + 4*E*. The last expression in Equation (2) is found by substituting *M* by its underlying terms.

In this model, *M* has a genetic variance of 25, of which 9 is caused by a genetic factor that also has a direct effect on the outcome *Y*. The remaining genetic variance is caused by a genetic factor that has no direct bearing on *Y*. In addition, *M* has an environment variance of 25. Thus, the SNP-based heritability of *M* is 50%. Finally, *M* is affected by the fixed-effect covariates in `covariates.txt`.

Also, under this model, *M* has a direct effect on *Y* equal to 1. Moreover, the aforementioned genetic factor *G* that directly affects both *M* and *Y*, has a direct effect of 2 on *Y*. Moreover, *Y* has an idiosyncratic environment factor, which adds 16 to its variance. The total genetic variance and environment variance of *Y* are both equal 41, putting the SNP-based heritability of *Y* also at 50%. Finally, *Y* is also affected by the fixed-effect covariates in `covariates.txt`.

Under this model, the genetic variance of *Y* that is mediated by *M* equals 25. This number effectively quantifies the so-called indirect effect that is often reported in the mediation literature. Here, this indirect effect reflects (i) the total effect genes have on *M* and (ii) the effect *M*, in turn, has on *Y*.

Moreover, under this model, if would consider *R* = *Y* - *Mb*, where *b* is the true effect of *M* on *Y* (here, *b* = 1) and, thus, *R* is the part of *Y* that remains if we would correct *Y* for *M* without any bias, then in this model *R* = 2*G* + 4*E*. In other words, *R* has an idiosyncratic genetic variance equal to four. Put differently, out of the full genetic variance of *Y* (which equals 41), only 4 is truly non-mediated. Thus, 4/41 = 9.76% of the genetic variance of *Y* is non-mediated.

Bearing these considerations, let's run the following `mgreml` command:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/mediation.txt \
                   --covar ./tutorial/covar.txt --mediation \
                   --out ./tutorial/try_mediation
```

Now, let's have a look at the output file `try_mediation.mediation.out`:

```
Mediation analysis in line with Rietveld et al. (2022):
Mediator M = Mediator; Outcome Y = Outcome
Effect M on Y (S.E.; Wald test statistic; asymptotic P value)* = 0.9852870634198768 (0.015112109436740342; 4250.8460961267; 0.0)
Total genetic variance of M (S.E.) = 24.84236880237618 (1.1967344214880011)
Total genetic variance of Y (S.E.) = 42.45446095524771 (2.005691645227429)
Indirect genetic effect = Genetic variance Y mediated by M (S.E.; Wald test statistic; asymptotic P value)* = 24.116738049063056 (1.3575913481004387; 315.57239745323005; 0.0)
Direct genetic effect = Genetic variance Y not mediated by M (S.E.; Wald test statistic; asymptotic P value)* = 4.46526039455577 (0.4541257330245469; 96.68104957436566; 0.0)
Proportion of genetic variance Y not mediated by M (S.E.) = 0.10517764904052629 (0.00920912144463799)
* Wald test statistic and P value under null hypothesis that parameter of interest = 0
Log-likelihood restricted model with no genetic variance of mediator = -31083.968230838338
Log-likelihood restricted model with no effect mediator on outcome = -32233.66196561707
Supremum restricted models is achieved under model with no genetic variance of mediator: test has 2 degrees of freedom 
Log-likelihood unrestricted model with genetic mediation = -30531.931576735766
Chi-square test statistic for presence of genetic mediation = 1104.0733082051447
with P-value = 0.0
```

Estimates are all less than two standard errors away from the true parameters of the structural model. Moreover, estimates in `try_mediation.HSq.out` also show that the estimated heritabilities are less than two standard errors removed from the true value (50% for both). In addition, based on Wald tests, observe that the estimated effect of *M* on *Y* is significant, the indirect effect of genes on *Y* via *M* is significant, and the direct effect of genes on *Y* is also significant. Finally, observe that the more reliable likelihood-ratio test also finds the estimated indirect genetic effect on *Y* (i.e., mediated by *M*) to be highly significant.

### Data formats and management

The input files that follow the options `--pheno`, `--covar`, `--covar-model`, `--genetic-model`, `--environment-model`, `--restricted-genetic-model`, `--restricted-environment-model` can be comma-, tab-, or space-separated. Just make sure to be completely consistent within each file. :warning: Please make sure your labels (e.g. phenotype labels, covariate labels, etc.) do not contain commas.

The options `--pheno`, `--covar`, `--covar-model`, `--genetic-model`, `--environment-model`, `--restricted-genetic-model`, `--restricted-environment-model` have modifiers `nolabelpheno`, `nolabelcovar`, and `nolabelfactor` to indicate when headers (e.g. phenotype labels) are absent. E.g. the following options are possible in `mgreml`: `--pheno myphen.txt nolabelpheno` and `--covar-model mycovmodel.txt nolabelpheno nolabelcovar`). However, as `mgreml` is a multivariate method, we strongly recommend always providing headers to `mgreml`, so everything is labelled in terms of input as well as output.

Naturally, `mgreml` performs basic data management, e.g. in terms of figuring out for which individuals we have phenotypic data as well as GRM data (and data on covariates, if applicable). In case `--covar-model` is used, `mgreml` also tests if there are any covariates that affect no phenotype at all, and if so, excludes such covariates.

Finally, `mgreml` can also perform relatedness pruning when using the `--grm-cutoff` option. E.g. the following command selects a subset of individuals such that relatedness in the GRM is nowhere in excess of 0.1:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --grm-cutoff 0.1 \
                   --out ./tutorial/pruned
```

which, in this case, causes the sample size to decrease by 44, as in shown in part of the log-file below:

```
2. CLEANING YOUR DATA
[...]
APPLYING RELATEDNESS CUTOFF
Current memory usage is 284MB
Removing individuals such that there is no relatedness in excess of 0.1 in GRM
First pass: dropped 40 individuals with only one relatedness value in excess of 0.1
Second pass: dropped 4 individuals with relatedness values in excess of 0.1
Relatedness cutoff has been applied
Remaining sample size is 4956
[...]
```

Here, `mgreml` follows the greedy algorithm developed by Boppana and Halldórsson (1992); [doi:10.1007/BF01994876](https://link.springer.com/article/10.1007/BF01994876). Importantly, `mgreml` does this pruning at such a stage and in such a manner that sample size is maximised. E.g. for a pair of individuals with a relatedness in excess of the threshold, we try to keep the observation with the lowest missingness in terms of phenotypes and covariates.

### Missing data and unbalancedness

In general, `mgreml` simply tries to preserve sample size at each turn. E.g. if an individual has missing values only for a subset of the phenotypes, `mgreml` keeps that individual in the data, by introducing phenotype-by-individual-specific dummies (i.e. dummies that control for individual *i* having a missing value for trait *t*). Even when a covariate is missing, sometimes parts of that observation can still be salvaged (i.e. if the missing covariate does not affect all phenotypes according to `--covar-model`).

Introducing dummies to control for gaps in the data can become computationally highly demanding. Controlling for fixed-effect covariates has a computational complexity of the order *NT*<sup> 2</sup> provided the number of unique covariates is of the order 1. However, if the missingness in each trait is on average proportional to sample size (*N*), then the total number of unique covariates needed to control for this missingness becomes of the order *NT*, and thereby the computational complexity of controlling for this missingness of the order *N*<sup> 3</sup>*T*<sup> 2</sup>, which is prohibitively complex for large *N* and *T*.

Therefore, `mgreml` has a `--drop-missings` option, whereby all individuals are dropped that have at least one missing phenotype and/or at least one missing covariate that is relevant (either because `--covar-model` has not been used, or because the file following `--covar-model` indicates the covariate with a missing value for a given individual affects at least one trait).

MGREML can handle missing values in the phenotype and/or covariate file when encoded using one of the following formats (comma-separated list of formats, with each format between single quotation marks):

‘’, ‘-999’, ‘#N/A’, ‘#N/A N/A’, ‘#NA’, ‘-1.#IND’, ‘-1.#QNAN’, ‘-NaN’, ‘-nan’, ‘1.#IND’, ‘1.#QNAN’, ‘’, ‘N/A’, ‘NA’, ‘NULL’, ‘NaN’, ‘n/a’, ‘nan’, ‘null’.

This list is based on https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html, which is the reference manual of the pandas function that MGREML uses to read in phenotype and covariate data.

:warning: Your phenotypes and covariates should not contain non-numeric stuff other than missing values!

### Advanced options

Note that `mgreml` has a few advanced options regarding the estimation algorithm. First, `--newton` forces `mgreml` to use a Newton algorithm for solving the optimisation problem instead of BFGS. Although in theory this approach requires fewer iterations than BFGS to converge, we recommend using BFGS: especially for large *T*, single BFGS iterations are so much faster than single Newton iterations, that the overall runtime of the BFGS algorithm is much lower than the runtime of the Newton algorithm.

In addition, `mgreml` considers the estimates to have converged if the root mean squared sum of the gradient vector (taking scale of traits and sample size into account) is below 10<sup>&minus;5</sup>. For conceptual ease, we refer to this criterion as length of the gradient. The option `--grad-tol` can be used to specify a different treshold.

:warning: We do **NOT** recommend deviating from 10<sup>&minus;5</sup> by more than one order of magnitude. E.g. you could use `--grad-tol 5E-5` or `--grad-tol 1e-6`. However, we de **NOT** recommend e.g. `--grad-tol 1E-9`, as such a threshold requires a degree of convergence that is often beyond numerical precision, nor do we recommend e.g. `--grad-tol 0.01`, as this is too lenient; the optimisation procedure simply has not converged when the latter convergence criterion is met.

When `mgreml` finds there to be strong multicollinearity in your phenotypes when controlling for the fixed-effect covariates, `mgreml` will return an error. In general, models with very high collinearity between phenotypes are poorly identified. Thus, in case `mgreml` returns this error, we recommend you consider a modified set of phenotypes (e.g. a strict subset of your phenotypes, or at least one phenotype defined slightly differently). In case you want to override this error, please use the `--ignore-collinearity` option. However, when using this option, `mgreml` may fail to converge at all, and should `mgreml` converge, the resulting estimates should be interpreted with great caution. Therfore, this option should only be used as a last resort.

Finally, `mgreml` also provides the option to perform pairwise bivariate estimation when analysing three or more traits. Such pairwise estimation is performed by using the `--pairwise` option. In case of strong multicollinearity between your phenotypes, this option can sometimes be helpful. However, in general, we would recommend figuring out what exactly causes the multicollinearity and whether it can be resolved in some other fashion (e.g. by analysing a slightly different subset of phenotypes), rather than switchting from multivariate estimation to pairwise bivariate estimation. Also, note that `--pairwise` cannot be combined with the following options: `--mediation`, `--(restricted-)genetic-model`, `--(restricted-)no-var-genetic`, `--(restricted-)environment-model`, `--factor-coefficients`, `--variance-components`, `--store-iter`, and/or `--(restricted-)reinitialise`. :warning: The `--pairwise` option generates *T*×(*T*−1)/2 sets of output files (one set for each unique combination of two traits). Thus, **the number of output files can be very large**.

## Overview of commands

An overview of all `mgreml` commands is listed below:

| Command | Usage |
| --- | --- |
| `-h`, `--help` | show help message and exit |
| `--grm PREFIX` | prefix of binary GRM |
| `--grm-cutoff THRESHOLD` | option to drop individuals using a greedy algorithm, such that there is no relatedness in GRM in excess of threshold for remaining individuals |
| `--adjust-pcs INTEGER [INTEGER]` | option to specify for how many leading principal components (PCs) from genetic data to adjust (to control for population stratification) and for how many trailing PCs to adjust (for computational efficiency); if just one non-negative integer is specified this is taken as the number of leading PCs to adjust for |
| `--pheno FILENAME [nolabelpheno]` | phenotype file: should be comma-, space-, or tab-separated, with one row per individual, with FID and IID as first two fields, followed by a field per phenotype; can be followed by optional flag `nolabelpheno`, e.g. `--pheno` `mypheno.txt nolabelpheno`, but we recommend to label phenotypes |
| `--mediation` | option to perform a genetic mediation analysis, in line with the structural equations model proposed by Rietveld et al. (2022) and based on estimates from a saturated bivariate model; the first phenotype in the phenotype file is assumed to act as mediator for the genetic component of the second phenotype in the phenotype file; all further phenotypes are ignored; cannot be combined with `--(restricted-)genetic-model`, `--(restricted-)rho-genetic`, `--(restricted-)no-var-genetic`, `--(restricted-)environment-model`, `--(restricted-)rho-environment`, and `--(restricted-)reinitialise` |
| `--drop-missings` | option to drop all observations from data with at least one missing phenotype or at least one missing covariate |
| `--no-intercept` | option to indicate an intercept should not be included automatically as covariate |
| `--covar FILENAME [nolabelcovar]` | optional covariate file: should be comma-, space-, or tab- separated, with one row per individual, with FID and IID as first two fields, followed by a field per covariate; can be followed by optional flag `nolabelcovar`, e.g. `--covar mycovar.txt nolabelcovar`, but we recommend to label covariates; :warning: do not include principal components from genetic data as covariates, use `--adjust-pcs` instead |
| `--covar-model FILENAME [nolabelpheno] [nolabelcovar]` | optional covariate model file: should be comma-, space-, or tab-separated, with one row per phenotype and one column per covariate; can be followed by optional flags `nolabelpheno` and/or `nolabelcovar`, but we recommend to label phenotypes and covariates; without `--covar-model`, all covariates are assumed to apply to all traits |
| `--genetic-model FILENAME [nolabelpheno] [nolabelfactor]` |  optional genetic model file: should be comma-, space-, or tab- separated, with one row per phenotype and one column per genetic factor; can be followed by optional flags `nolabelpheno` and/or `nolabelfactor`, but we recommend to label phenotypes and genetic factors |
| `--rho-genetic 0`, `1`  | option followed by `0` or `1`, forcing all genetic correlations to take on the specified value; this flag cannot be combined with `--genetic-model` |
| `--no-var-genetic` | option to force all genetic variances to equal zero; this flag cannot be combined with `--genetic-model` and/or `--rho-genetic` |
| `--restricted-genetic-model FILENAME [nolabelpheno] [nolabelfactor]` | optional restricted genetic model file: should be comma-, space-, or tab-separated, with one row per phenotype and one column per genetic factor; can be followed by optional flags `nolabelpheno` and/or `nolabelfactor`, but we recommend to label phenotypes and genetic factors |
| `--restricted-rho-genetic 0`, `1` | option followed by `0` or `1`, forcing all genetic correlations in the restricted model to take on the specified value; this flag cannot be combined with `--restricted-genetic-model` |
| `--restricted-no-var-genetic` | option to force all genetic variances in the restricted model to equal zero; this flag cannot be combined with `--restricted-genetic-model` and/or `--restricted-rho-genetic` |
| `--environment-model FILENAME [nolabelpheno] [nolabelfactor]` |  optional environment model file: should be comma-, space-, or tab-separated, with one row per phenotype and one column per environment factor; can be followed by optional flags `nolabelpheno` and/or `nolabelfactor`, but we recommend to label phenotypes and environment factors |
| `--rho-environment 0`  | option followed by `0`, forcing all environment correlations to zero; this flag cannot be combined with `--environment-model` |
| `--restricted-environment-model FILENAME [nolabelpheno] [nolabelfactor]` | optional restricted environment model file: should be comma-, space-, or tab-separated, with one row per phenotype and one column per environment factor; can be followed by optional flags `nolabelpheno` and/or `nolabelfactor`, but we recommend to label phenotypes and environment factors |
| `--restricted-rho-environment 0` | option followed by `0`, forcing all environment correlations in the restricted model to zero; this flag cannot be combined with `--restricted-environment-model` |
| `--no-se` | option to skip calculation of standard errors and covariance matrix of estimates |
| `--factor-coefficients` | option to report estimated factor coefficients |
| `--variance-components` | option to report estimated variance components |
| `--newton` | option to use Newton method instead of BFGS; not recommended, unless the model is well-defined, starting values are of good quality, and the number of traits is small |
| `--grad-tol THRESHOLD` | option to set convergence threshold on the length of the gradient vector per parameter, per observation, different from the default value of `1E-5` |
| `--store-iter INTEGER` | option to specify every how many iterations you want to store results |
| `--reinitialise FILENAME` | option to reinitialise `mgreml` for a model and its estimates from a `.pkl` file stored by `--store-iter` |
| `--restricted-reinitialise FILENAME` | option to reinitialise `mgreml` for a restricted model and its estimates from a `.pkl` file generated by `--store-iter` |
| `--pairwise` | option to perform pairwise bivariate estimation instead of multivariate estimation; cannot be combined with `--mediation`, `--(restricted-)genetic-model`, `--(restricted-)no-var-genetic`, `--(restricted-)environment-model`, `--factor-coefficients`, `--variance-components`, `--store-iter`, and/or `--(restricted-)reinitialise`; :warning: the number of output files can be very large |
| `--ignore-collinearity` | option to ignore multicollinearity between phenotypes; please use this option only as a last resort; model may be poorly identified when your phenotype data is perfectly collinear; preferred route to solve collinearity is to consider e.g. a subset of phenotypes |
| `--out PREFIX` |  prefix of output files |

## Updating `mgreml`

You can update to the newest version of `mgreml` using `git`. First, navigate to your `mgreml` directory (e.g. `cd mgreml`), then run
```
git pull
```
If `mgreml` is up to date, you will see 
```
Already up-to-date.
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

If you have modified the `mgreml` source code yourself, `git pull` may fail with an error such as `error: Your local changes [...] would be overwritten by merge`. 

In case the Python dependencies have changed, you can update the `mgreml` environment with

```
conda env update --file mgreml.yml
```

## Support

Before contacting us, please try the following:

1. Go over the tutorial in this `README.md` file
2. Go over the method, described in detail in the supplementary information of the paper (see [Citation](#citation))

### Contact

In case you have a question that is not resolved by going over the preceding two steps, or in case you have encountered a bug, please send an e-mail to r\[dot\]devlaming\[at\]vu\[dot\]nl.

## Citation

In general, if you use the software, please cite

[R. de Vlaming, E.A.W. Slob, P.R. Jansen, A. Dagher, P.D. Koellinger, P.J.F. Groenen, and C.A. Rietveld (2021). Multivariate analysis reveals shared genetic architecture of brain morphology and human behavior. *Commun Biol* **4**, 1180](https://doi.org/10.1038/s42003-021-02712-y)

In addition, if you use the `--mediation` option, please also cite

[C.A. Rietveld, R. de Vlaming, E.A.W. Slob (2022). *tba*]

## Derivations

For full details on the derivation of the MGREML method, see the [Supplementary Information](https://www.biorxiv.org/content/biorxiv/early/2021/04/19/2021.04.19.440478/DC1/embed/media-1.pdf), available on bioRxiv.

For derivations on the structural model used in the mediation analysis, see *tba*.

## License

This project is licensed under GNU GPL v3.

## Authors

Ronald de Vlaming (Vrije Universiteit Amsterdam)

Eric Slob (University of Cambridge)
