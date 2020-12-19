# MGREML (Multivariate GREML) `BETA VERSION 0.01`

`mgreml` is a command-line tool for rapid estimation of SNP-based heritability and genetic correlations for (nearly) balanced data on many traits in a single analysis using a genomic-relatedness matrix (GRM) using Python 3.x.

`mgreml` can easily handle estimation of the full genetic correlation matrix for up to 100 traits observed in 20,000 individuals. `mgreml` allows users to specify structural models and test hypotheses regarding nested models (e.g. no genetic correlations). In addition, the tool can handle a considerable amount of fixed-effect covariates and a very minor degree of phenotypic missingness.

Finally, `mgreml` has options to return the full set of factor coefficients and variance components, as well as the complete covariance matrix of those estimates.

Please note that this is still a beta version.

## Installation

Before downloading `mgreml`, please make sure [Git](https://git-scm.com/downloads) and [Anaconda](https://www.anaconda.com/) with **Python 3.x** are installed.

In order to download `mgreml`, open a command-line interface by starting [Anaconda Prompt](https://docs.anaconda.com/anaconda/user-guide/getting-started/), navigate to your working directory, and clone the `mgreml` repository using the following command:

```  
git clone https://github.com/devlaming/mgreml.git
```

Then enter the newly created `mgreml` directory using:

```
cd mgreml
```

Then run the following commands to create a custom Python environment which has all of `mgreml`'s dependencies (i.e. an environment that has packages such as `numpy` and `pandas` pre-installed):

```
conda env create --file mgreml.yml
conda activate mgreml
```

(or `activate mgreml` instead of `conda activate mgreml` on some machines).

In case you cannot create a customised conda environment (e.g. because of insufficient user rights) or simply prefer to use Anaconda Navigator or `pip` to install packages e.g. in your base environment rather than a custom environment, please note that `mgreml` only requires Python 3.x with the packages `numpy`, `pandas`, `scipy`, and `tqdm` installed.

Once the above has completed, you can now run

```
python ./mgreml -h
```

to print a list of all command-line options. If this command fails, something has gone wrong during installation.

*Windows users*: in case the preceding command fails, try replacing slashes (i.e. `/`) in all your `mgreml` commands by backslashes (i.e. `\`), so e.g. try

```
python .\mgreml -h
```

## Tutorial

In this short tutorial we will go over the basic functions of `mgreml`. First, go over the steps in Installation.

Now that you have cloned the `mgreml` repository, and `mgreml` is up-and-running, the main directory of `mgreml` should contain a subdirectory called `tutorial`. This directory in turn contains several files, including `pheno.txt` and `covar.txt`. Details on how this dataset has been generated using simulation can be found in the python script in `./tutorial/simulate.py`

Let's first inspect the `pheno.txt` file. This file contains data in tab-separated format on ten phenotypes observed in a set of 5,000 individuals. The first two columns list family and individual ID, followed by the phenotypes:

| FID | IID | Some pheno 101 | Some pheno 102 | ...  | Some pheno 109 | Some pheno 110 |
| --- | --- | --- | --- | --- | --- | --- |
| FID 1 | IID 5001 | -2.072 | -3.676 | ...  | 4.641 | 7.931 |
| FID 2 | IID 5002 | -1.472 | -1.467 | ...  | 6.098 | 3.570 |
| FID 3 | IID 5003 | -3.392 | -7.277 | ...  | -0.832 | -5.750 |
| ...  | ...  | ...  | ...  | ...  | ...  | ...  |
| FID 4998 | IID 9998 | 2.575 | 2.740 | ...  | 3.328 | -6.982 |
| FID 4999 | IID 9999 | -3.072 | -0.306 | ...  | 2.530 | -1.255 |
| FID 5000 | IID 10000 | -4.220 | 1.117 | ...  | 2.806 | 3.159 |

Although `mgreml` can handle phenotype data without header, using a modifier that we discuss later, we recommend always including headers in your data, so the `mgreml` output refers e.g. to specific phenotype names rather than ambiguous indices such as `1`, `2`, etc.

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
| Some pheno 101 | 0.052 | 0.029 |
| Some pheno 102 | 0.002 | 0.031 |
| Some pheno 103 | 0.016 | 0.029 |
| Some pheno 104 | 0.013 | 0.029 |
| Some pheno 105 | 0.037 | 0.030 |
| Some pheno 106 | 0.171 | 0.029 |
| Some pheno 107 | 0.004 | 0.030 |
| Some pheno 108 | 0.001 | 0.030 |
| Some pheno 109 | 0.016 | 0.032 |
| Some pheno 110 | 0.157 | 0.029 |

Comparing these estimates to the true values in `./tutorial/true.HSq.txt`, printed below, we see that our estimates seem to be biased.

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

The simple reason for this bias is that we did not control for our fixed-effect covariates, in `./tutorial/covar.txt`, which affect the traits of interest. So we need to use the `--covar` option to specify our fixed-effect covariates. This boils down to the following `mgreml` command:

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt --out ./tutorial/covs
```

If we compare the new estimates of heritability (see below) to the true values, taking the standard errors of the estimates into account, we see the strong bias is gone.

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

Notice upon inspection of `covar.txt` that it does include the intercept but does not include any principal components (PCs) from the genetic data.

First, if you want to control for the intercept, the intercept **MUST** be included as a separate covariate in your covariate file (i.e. as a vector of ones), as `mgreml` assumes the intercept is absent by default (opposed e.g. to `gcta`).

Second, the file with covariates should **NEVER** contain PCs from your genetic data, as `mgreml` already removes the effects of population stratification in the so-called canonical transformation. By default, `mgreml` removes the effects of 20 leading PCs from your genetic data. The effective sample size is reduced by 20 as a result of this correction for PCs.

In case you want to change the number of PCs you control for, do **NOT** add these PCs to your file with covariate data. Instead, use the `--adjust-pcs` option, followed by the total number of leading PCs you want to control for. E.g. `--adjust-pcs 20` is equivalent to the default setting, `--adjust-pcs 40` controls for the 40 leadings PCs, and `--adjust-pcs 0` controls for no PCs at all (not recommended). In these three cases, the sample size is reduced by 20, 40, and zero respectively.

For advanced users, the `--adjust-pcs` option can also be followed by a second number, indicating the number of trailing eigenvectors from your GRM to adjust for. E.g. `--adjust-pcs 100 1000` controls for 100 leading eigenvectors from your GRM and 1000 trailing eigenvectors. Doing this decreases the overall sample size by 100 + 1000 = 1100. By default no trailing eigenvectors are adjusted for. However, if the trailing eigenvalues are sufficiently small, a considerable number of trailing eigenvectors may be adjusted for, boosting CPU time (as the sample size becomes smaller) without diminishing statistical efficiency of your analysis too much (as effectively only less informative bits of data are ignored).

In addition to reporting the heritabilities and their standard errors, `mgreml` also automatically reports genetic and environment correlations, as well as their standard errors.

In case you care neither about standard errors nor the covariance matrix of estimates, you can use the `--no-se` option. Especially for a large number of traits, computing the standard errors is computationally demanding, as this requires calculating the average information matrix, which has a computational complexity of the order *NT*<sup> 4</sup>, where *T* denotes the number of traits and *N* the number of observations.

`mgreml` also automatically reports the fixed-effect estimates (a.k.a. GLS estimates), including the covariance matrix of those estimates, and their standard errors.

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

Now, `different_covs.GLS.est.out`, in the folder `tutorial`, only shows fixed-effect estimates for covariates that affect the given trait according to `covar_model.txt`:

| trait | covariate | beta hat | standard error |
| --- | --- | --- | --- |
| Some pheno 101 | intercept | -1.413 | 0.041 |
| Some pheno 102 | intercept | -1.777 | 0.031 |
| Some pheno 102 | my covar 301 | -1.432 | 0.017 |
| ...  | ...  | ...  | ...  |
| Some pheno 109 | intercept | 2.269 | 0.037 |
| Some pheno 109 | my covar 308 | -1.282 | 0.014 |
| Some pheno 110 | intercept | -0.397 | 0.057 |
| Some pheno 110 | my covar 309 | 1.172 | 0.031 |

E.g. `my covar 301` does not affect `Some pheno 101` in this case.

Analogous to `--covar-model`, users can also specify which genetic factor affects which trait and which environment factor affects which trait. Such specifications can be passed to `mgreml` using the `--genetic-model` and `--environment-model` options. Note, that any such user-specified structural model must be identified. Moreover, for the factor specification of the environment, `mgreml` requires as many factors as there are traits.

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

Similarly, the estimates of environment correlations, in `custom_model.RhoE.out`, show these are all estimated at zero, also as expected under this model.

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
Log-likelihood of model = -76460.81732177256,
based on data on 10 traits and 4980 observations,
with a model consisting of 1 genetic factors and 10 environment factors,
comprising 10 free genetic factor coefficients and 10 free environment factor coefficients in turn.
Estimates converged after 37 BFGS iterations 
```

Regarding specific factor models, `mgreml` also allows users to force all genetic variances to zero using `--no-var-genetic` and doing the same in the restricted model using `--restricted-no-var-genetic`. The option `--no-var-genetic` cannot be combined with `--rho-genetic` and/or  `--genetic-model`. Similarly, `--restricted-no-var-genetic` cannot be combined with `--restricted-rho-genetic` and/or  `--restricted-genetic-model`.

In case you have estimated a model, either according to some structural model e.g. using `--genetic-model`, or just the saturated model we started with, you can make `mgreml` report the factor coefficients (i.e. the effect of each factor on each trait) by using the `--all-coefficients` option. Using this option not only reports the estimated factor coefficients, but also the covariance matrix of those estimates. This covariance matrix may grow very large for large *T*.

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
| Some pheno 109 | environment factor 8 | 0.227 |
| Some pheno 110 | environment factor 8 | -0.416 |
| Some pheno 110 | environment factor 9 | 0.360 |

The file `full.coeff.var.out` contains a 110-by-110 matrix representing the covariance matrix of those estimates. 

Similarly, `mgreml` can also return the estimated variance components (again either based on some structural model, or just the saturated model), including the covariance matrix of those estimated variance components. To get these results, use the `--variance-components` option. E.g. the command

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt \
                --variance-components \
                --out ./tutorial/vcs
```

generates, amongst others, the file `vcs.VCs.out`, which contains 110 estimated covariance components in this case, of which a few lines are shown below:

| component | first trait | second trait | estimate |
| --- | --- | --- | --- |
| genetic covariance | Some pheno 101 | Some pheno 101 | 0.987 |
| genetic covariance | Some pheno 101 | Some pheno 102 | 0.081 |
| genetic covariance | Some pheno 101 | Some pheno 103 | 0.202 |
| ...  | ...  | ...  | ...  |
| environment covariance | Some pheno 109 | Some pheno 109 | 0.646 |
| environment covariance | Some pheno 109 | Some pheno 110 | 0.164 |
| environment covariance | Some pheno 110 | Some pheno 110 | 0.960 |

The file `vcs.VCs.var.out` contains a 110-by-110 matrix representing the covariance matrix of those estimates. 

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
Log-likelihood of nested model (null hypothesis) = -76460.81732177256,
based on data on 10 traits and 4980 observations,
with a model consisting of 1 genetic factors and 10 environment factors,
comprising 10 free genetic factor coefficients and 10 free environment factor coefficients in turn.
Estimates converged after 37 BFGS iterations 

Log-likelihood of parent model (alternative hypothesis) = -66849.63485188212,
based on data on 10 traits and 4980 observations,
with a model consisting of 10 genetic factors and 10 environment factors,
comprising 55 free genetic factor coefficients and 55 free environment factor coefficients in turn.
Estimates converged after 268 BFGS iterations 

Results of likelihood-ratio test with 90 degrees of freedom:
Chi-square test statistic is 19222.364939780877
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

causes `mgreml` to store results every 50 iterations. Then, if the preceding analysis has reached e.g. up until iteration 250 before a power outage, we could reinitialise later on using the following command:

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt \
                --reinitialise ./tutorial/covar.estimates.iter.250.bfgs.pkl \
                --out ./tutorial/covar_reinitialised
```

Notice that as such `.pkl` files already implicitly contain the full model specification, the option `--reinitialise` cannot be combined with options such as `--genetic-model`, `--rho-environment`, `--no-var-genetic`, and so on.

In case `--store-iter` is used when estimating a nested versus alternative model (i.e. in conjunction with one of the `--restricted-...` options), `--store-iter` stores two sets of `.pkl` files, namely one set with filenames containing `.estimates.` (denoting the alternative model) and the other containing `.estimates0.` (denoting the nested model).

`.pkl` files can also be used to reinitialise a restricted model, using the `--restricted-reinitialise` option. E.g. the command

```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt \
                --restricted-rho-genetic 1 \
                --restricted-rho-environment 0 \
                --store-iter 5 \
                --out ./tutorial/restricted_rhoG1_rhoE0
```
causes two sets of `.pkl` files to be stored (i.e. a file for every 5 iterations of both the restricted and alternative model) and
```
python ./mgreml --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                --covar ./tutorial/covar.txt \
                --reinitialise ./tutorial/restricted_rhoG1_rhoE0.estimates.iter.265.bfgs.pkl \
                --restricted-reinitialise ./tutorial/restricted_rhoG1_rhoE0.estimates0.iter.25.bfgs.pkl \
                --out ./tutorial/restricted_rhoG1_rhoE0_reinitialised
```
reinitialises estimation for the null and alternative model from appropriate `.pkl` files. Notice that analogous to `--reinitialise`, the `--restricted-reinitialise` option cannot be combined with options such as `--restricted-environment-model` and `--restricted-rho-genetic`, as the `.pkl` file already contains the full model specification.

`mgreml` performs basic data management, e.g. in terms of figuring out for which individuals we have phenotypic as well as GRM data (and data on covariates, if applicable). In case `--covar-model` is used `mgreml` also tests if there are any covariates that affect no phenotype at all, and if so, excludes such covariates.

*Warning: the following option still needs to be implemented:* `mgreml` also performs relatedness pruning when using the `--rel-cutoff` option. E.g. `--rel-cutoff 0.025` selects a subset of individuals such that relatedness in the GRM is nowhere in excess of 0.025. `mgreml` follows the same algorithm for such a cutoff as [PLINK](https://www.cog-genomics.org/plink/1.9/distance#rel_cutoff). Importantly, `mgreml` does this pruning at such a stage that sample size is maximised (e.g. for a pair of individuals with a relatedness in excess of the threshold, we do not drop the individual for whom we have phenotype data and keep the individual for whom we do not have any phenotype data at all).

In general, `mgreml` simply tries to maximise sample size at each turn. E.g. if an individual has missing values only for a subset of the phenotypes, `mgreml` keeps that individual in the data, by introducing phenotype-by-individual-specific dummies (i.e. dummies that control for individual *i* having a missing value for trait *t*). Even when a covariate is missing, sometimes parts of that observation can still be salvaged (i.e. if the missing covariate does not affect all phenotypes according to `--covar-model`).

However, introducing these dummies to control for gaps in the data can become computationally highly demanding. Controlling for fixed effect covariates has a computational complexity of the order *NT* <sup>2</sup> provided the number of unique covariates is of the order 1. However, if e.g. missingness in each trait is a proportion of sample size, then the total number of unique covariates to control for this missingness becomes of the order *NT*, and thereby the computational complexity of controling for this missingness of the order *N* <sup>2</sup> *T* <sup>3</sup>, which is prohibitively complex for large *N* and *T*.

Therefore, `mgreml` has a `--drop-missings` option, whereby all individuals are dropped that have at least one missing phenotype and/or at least one missing covariate that is relevant (either because `--covar-model` has not been used, or because the file following `--covar-model` indicates the covariate with a missing value for a given individual affects at least one trait).

Note that `mgreml` has a few advanced options regarding the estimation algorithm. First, `--newton` forces `mgreml` to use a Newton algorithm for solving the optimisation problem instead of BFGS. Although in theory this approach requires fewer iterations that BFGS, we strongly recommend sticking to BFGS: BFGS iterations are faster and BFGS is numerically much more stable, especially for large *T*.

In addition, `mgreml` considers the estimates to have converged if the root mean squared sum of the gradient vector (taking scale of traits and sample size into account) is below 10<sup>&minus;5</sup>. For conceptual ease, we refer to this criterion as length of the gradient. The option `--grad-tol` can be used to specify a different treshold. We do **NOT** recommend deviating from 10<sup>&minus;5</sup> by more than one order of magnitude. E.g. you could use `--grad-tol 5E-5` or `--grad-tol 1e-6`. However, we de **NOT** recommend e.g. `--grad-tol 1E-9`, as such a threshold requires a degree of convergence that is beyond numerical precision, nor do we commend e.g. `--grad-tol 0.01`, as this is too lenient; the optimisation procedure simply has not converged when the latter convergence criterion is met.

Finally, input files following the options `--pheno`,  `--covar`,  `--covar-model`, `--genetic-model`, `--environment-model`, `--restricted-genetic-model`, `--restricted-environment-model` can be comma, tab, or space-separated. Just make sure to be consistent within each file. Also, (1) when your labels (e.g. phenotype labels or FID and IIDs) contain spaces, make sure the file is not space-separated and (2) when those labels contain commas, make sure the file is not comma-separated. As a final remark, note that the command `python ./mgreml -h` shows that these options all have modifiers, such as `nolabelpheno`, to indicate when headers (e.g. phenotype labels) are absent using (e.g. using the option `--pheno my_pheno_file.txt nolabelpheno`). However, we recommend always providing headers to `mgreml`, so everything is labelled in terms of input as well as output.

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

If you have modified the `mgreml` source code yourself, `git pull` may fail with an error such as `error: Your local changes [...] would be overwritten by merge`. 

In case the Python dependencies have changed, you can update the `mgreml` environment with

```
conda env update --file mgreml.yml
```

## Support

Before contacting us, please try the following:

1. Go over the tutorial in this `README.md` file
2. Go over the method, described in detail in the supplementary information of the paper (citation below)

If that doesn't work, you can get in touch with us via the [e.g. google group; to be decided](...).

## Citation

If you use the software, please cite

[]()

## License

This project is licensed under GNU GPL v3.

## Authors

Ronald de Vlaming (Vrije Universiteit Amsterdam)

Eric Slob (University of Cambridge)


