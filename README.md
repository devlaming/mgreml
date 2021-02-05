# MGREML (Multivariate GREML) `BETA VERSION 0.01`

`mgreml` is a command-line tool for rapid estimation of SNP-based heritability and genetic correlations for (nearly) balanced data on many traits in a single analysis using a genomic-relatedness matrix (GRM) using Python 3.x.

`mgreml` can easily handle estimation of the full genetic correlation matrix for up to 100 traits observed in 20,000 individuals. `mgreml` allows users to specify structural models and test hypotheses regarding nested models (e.g. no genetic correlations). In addition, the tool can handle a considerable amount of fixed-effect covariates and a very minor degree of phenotypic missingness.

Finally, `mgreml` has options to return the full set of factor coefficients and variance components, as well as the complete covariance matrix of those estimates.

Please note that this is still a beta version.

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

In case you cannot create a customised conda environment (e.g. because of insufficient user rights) or simply prefer to use Anaconda Navigator or `pip` to install packages e.g. in your base environment rather than a custom environment, please note that `mgreml` only requires Python 3.x with the packages `numpy`, `pandas`, `scipy`, and `tqdm` installed.

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

In this tutorial we will go over the basic functions of `mgreml`. First, go over the steps in Installation.

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

Although `mgreml` in principle can handle phenotype data without header (using a modifier that we discuss later on), we recommend always including headers in your data, so e.g. your phenotypes are labelled, allowing `mgreml` output to refer to specific phenotype names rather than ambiguous indices such as `1`, `2`, `3` etc.

For the same set of individuals, you have a binary genomic-relatedness matrix (a.k.a. GRM) e.g. computed using [PLINK](https://www.cog-genomics.org/plink/) or [GCTA](https://cnsgenomics.com/software/gcta/). In this case, the set of binary GRM files comprises `data.grm.bin`, `data.grm.N.bin`, and `data.grm.id`. We refer to this set of binary GRM files by its prefix, i.e. `data`.

The command for running an `mgreml` analysis on this data without correcting for any covariates at all is as follows:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \ 
                   --no-intercept --out ./tutorial/nocovs
```

Upon carrying out this command, `mgreml` will first report the follow command will be carried out:

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
| Some pheno 101 | 0.052 | 0.029 |
| Some pheno 102 | 0.002 | 0.031 |
| Some pheno 103 | 0.017 | 0.029 |
| Some pheno 104 | 0.013 | 0.029 |
| Some pheno 105 | 0.038 | 0.030 |
| Some pheno 106 | 0.171 | 0.029 |
| Some pheno 107 | 0.005 | 0.030 |
| Some pheno 108 | 0.001 | 0.030 |
| Some pheno 109 | 0.016 | 0.032 |
| Some pheno 110 | 0.156 | 0.029 |

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

The simple reason for this bias is that we did not control for any fixed-effect covariates. By removing the `--no-intercept` option, `mgreml` automatically adds one fixed effect per phenotype, namely a fixed effect that controls for differences in mean across phenotypes:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --out ./tutorial/intercept
```

Resulting SNP heritability estimates in `./tutorial/intercept.HSq.out`, however, show our estimates are still strongly biased:

| trait | heritability | standard error |
| --- | --- | --- |
| Some pheno 101 | 0.070 | 0.030 |
| Some pheno 102 | 0.004 | 0.030 |
| Some pheno 103 | 0.016 | 0.029 |
| Some pheno 104 | 0.013 | 0.029 |
| Some pheno 105 | 0.041 | 0.030 |
| Some pheno 106 | 0.190 | 0.029 |
| Some pheno 107 | 0.008 | 0.029 |
| Some pheno 108 | 0.001 | 0.030 |
| Some pheno 109 | 0.027 | 0.030 |
| Some pheno 110 | 0.157 | 0.029 |

The reasons this bias persists is that more fixed-effect covariates are at play than just the intercept. The file `./tutorial/covar.txt` contains the other covariates that affect the traits of interest. So we need to use the `--covar` option to specify these additional fixed-effect covariates. This boils down to the following `mgreml` command:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt --out ./tutorial/covs
```

Notice that analyses including covariates are computationally slightly more demanding. E.g. in this case we have 10 covariates (i.e. the intercept + 9 additional covariates in `./tutorial/covar.txt`), each of which is allowed to have a different effect on each trait. As we have 10 traits, this means we have 100 fixed effects in total, which our model needs to take into account.

If we compare the new estimates of heritability (see below) to the true values, taking the standard errors of the estimates into account, we see the strong bias is gone.

| trait | heritability | standard error |
| --- | --- | --- |
| Some pheno 101 | 0.450 | 0.025 |
| Some pheno 102 | 0.015 | 0.030 |
| Some pheno 103 | 0.031 | 0.029 |
| Some pheno 104 | 0.044 | 0.029 |
| Some pheno 105 | 0.482 | 0.025 |
| Some pheno 106 | 0.487 | 0.025 |
| Some pheno 107 | 0.037 | 0.029 |
| Some pheno 108 | 0.017 | 0.029 |
| Some pheno 109 | 0.166 | 0.029 |
| Some pheno 110 | 0.742 | 0.019 |

Notice upon inspection of `covar.txt` that it does not include any principal components (PCs) from the genetic data.

In fact, in `mgreml`, the file with your covariates should **NEVER** contain PCs from your genetic data, as the tool already removes the effects of population stratification during the so-called canonical transformation. By default, `mgreml` removes the effects of 20 leading PCs from your genetic data. The effective sample size is reduced by 20 as a result of this correction for PCs.

In case you want to change the number of PCs you control for, do **NOT** add these PCs to your file with covariate data. Instead, use the `--adjust-pcs` option, followed by the total number of leading PCs you want to control for. E.g. `--adjust-pcs 20` is equivalent to the default setting, `--adjust-pcs 40` controls for the 40 leadings PCs, and `--adjust-pcs 0` controls for no PCs at all (not recommended). In these three cases, the sample size is reduced by 20, 40, and zero respectively.

For instance, the command

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --adjust-pcs 1000 --out ./tutorial/many_pcs
```

causes `mgreml` to adjust for 1000 leading PCs from the genetic data.

As there is no population stratification in our data (by virtue of our simulation design), this means adjusting for so many PCs will just reduce precision of our estimates, without eliminating any bias. If we look at `many_pcs.HSq.out` we see that our estimates indeed have considerably higher standard errors:

| trait | heritability | standard error |
| --- | --- | --- |
| Some pheno 101 | 0.459 | 0.035 |
| Some pheno 102 | 0.070 | 0.050 |
| Some pheno 103 | 0.029 | 0.050 |
| Some pheno 104 | 0.099 | 0.049 |
| Some pheno 105 | 0.468 | 0.035 |
| Some pheno 106 | 0.499 | 0.034 |
| Some pheno 107 | 0.034 | 0.051 |
| Some pheno 108 | 0.022 | 0.050 |
| Some pheno 109 | 0.187 | 0.046 |
| Some pheno 110 | 0.724 | 0.025 |

For advanced users, the `--adjust-pcs` option can also be followed by a second number, indicating the number of trailing eigenvectors from your GRM to adjust for. E.g. `--adjust-pcs 100 1000` controls for 100 leading eigenvectors from your GRM and 1000 trailing eigenvectors. Doing this decreases the overall sample size by 100 + 1000 = 1100.

By default no trailing eigenvectors are adjusted for. However, if the trailing eigenvalues are sufficiently small, a considerable number of trailing eigenvectors may be adjusted for, boosting CPU time (as the overall sample size becomes smaller) without diminishing statistical efficiency of your analysis too much (as effectively only less informative bits of data are ignored). Note that adjusting for trailing eigenvectors may require you to use `--no-intercept`, as the last eigenvector from the GRM tends to be highly multicollinear with the intercept.

In addition to reporting the heritabilities and their standard errors, `mgreml` also automatically reports genetic and environment correlations, as well as their standard errors.

In case you care neither about standard errors nor the covariance matrix of estimates, you can use the `--no-se` option. Especially for a large number of traits, computing the standard errors is computationally demanding, as this requires calculating the average information matrix, which has a computational complexity of the order *NT*<sup> 4</sup>, where *T* denotes the number of traits and *N* the number of observations.

`mgreml` also automatically reports the fixed-effect estimates (a.k.a. GLS estimates), including the covariance matrix of those estimates, and their standard errors. If the `--no-se` option is used, the GLS estimates will not include their standard errors and covariance matrix.

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
| Some pheno 101 | intercept | -1.413 | 0.041 |
| Some pheno 102 | intercept | -1.777 | 0.031 |
| Some pheno 102 | my covar 301 | -1.432 | 0.017 |
| ... | ... | ... | ... |
| Some pheno 109 | intercept | 2.269 | 0.037 |
| Some pheno 109 | my covar 308 | -1.282 | 0.014 |
| Some pheno 110 | intercept | -0.397 | 0.057 |
| Some pheno 110 | my covar 309 | 1.172 | 0.031 |

E.g. `my covar 301` does not affect `Some pheno 101` in this case.

Analogous to `--covar-model`, users can also specify which genetic factor affects which trait and which environment factor affects which trait. Such specifications can be passed to `mgreml` using the `--genetic-model` and `--environment-model` options. Note that any such user-specified structural model must be identified. Moreover, for the factor specification of the environment, `mgreml` requires as many factors as there are traits.

For example, we could impose a factor structure, where there is only one genetic factor, and where there are *T*=10 environment factors, each affecting a different trait. Effectively, this boils down to a model with genetic correlations all equal to one and environment correlations all equal to zero. These factor structures are shown in the files `gen_model.txt` and `env_model.txt` both found in the `tutorial` folder. Both files contain a binary table, with elements equal to one, where a given factor is permitted to affect the given phenotype, and equal to zero otherwise.

To estimate this structural model, we can simply carry out the following command:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
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
Log-likelihood of model = -76460.81732259798,
based on data on 10 traits and 4980 observations,
with a model consisting of 1 genetic factors and 10 environment factors,
comprising 10 free genetic factor coefficients and 10 free environment factor coefficients in turn.
Estimates converged after 25 BFGS iterations 
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

In case you have estimated a model, either according to some structural model e.g. using `--genetic-model`, or just the saturated model we started with, `mgreml` can report the factor coefficients (i.e. the estimated effect of each factor on each trait) by using the `--factor-coefficients` option. Using this option not only reports the estimated factor coefficients, but also the covariance matrix of those estimates (unless `--no-se` is used). This covariance matrix may grow very large for large *T*.

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
| Some pheno 101 | genetic factor 0 | 0.993 | 0.032 |
| Some pheno 102 | genetic factor 0 | 0.081 | 0.038 |
| Some pheno 103 | genetic factor 0 | 0.203 | 0.041 |
| ... | ... | ... | ... |
| Some pheno 109 | environment factor 8 | 0.227 | 0.077 |
| Some pheno 110 | environment factor 8 | -0.414 | 0.143 |
| Some pheno 110 | environment factor 9 | 0.362 | 0.139 |

The file `factors.coeff.var.out` contains a 110-by-110 matrix representing the covariance matrix of those estimates. 

Similarly, `mgreml` can also return the estimated variance components (again either based on some structural model, or just the saturated model), including the covariance matrix of those estimated variance components (unless `--no-se` is used). To get these results, use the `--variance-components` option. E.g. the command

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --variance-components \
                   --out ./tutorial/components
```

generates, amongst others, the file `components.VCs.out`, which contains 110 estimated covariance components in this case, of which a few lines are shown below:

| component | first trait | second trait | estimate | standard error |
| --- | --- | --- | --- | --- |
| genetic covariance | Some pheno 101 | Some pheno 101 | 0.986 | 0.064 |
| genetic covariance | Some pheno 101 | Some pheno 102 | 0.081 | 0.036 |
| genetic covariance | Some pheno 101 | Some pheno 103 | 0.202 | 0.043 |
| ... | ... | ... | ... | ... |
| environment covariance | Some pheno 109 | Some pheno 109 | 0.646 | 0.024 |
| environment covariance | Some pheno 109 | Some pheno 110 | 0.164 | 0.029 |
| environment covariance | Some pheno 110 | Some pheno 110 | 0.960 | 0.066 |

The file `components.VCs.var.out` contains a 110-by-110 matrix representing the covariance matrix of those estimates. 

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
Log-likelihood of nested model (null hypothesis) = -76460.81732259798,
based on data on 10 traits and 4980 observations,
with a model consisting of 1 genetic factors and 10 environment factors,
comprising 10 free genetic factor coefficients and 10 free environment factor coefficients in turn.
Estimates converged after 25 BFGS iterations 

Log-likelihood of parent model (alternative hypothesis) = -66849.63346504091,
based on data on 10 traits and 4980 observations,
with a model consisting of 10 genetic factors and 10 environment factors,
comprising 55 free genetic factor coefficients and 55 free environment factor coefficients in turn.
Estimates converged after 269 BFGS iterations 

Results of likelihood-ratio test with 90 degrees of freedom:
Chi-square test statistic is 19222.367715114146
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

Also, to test whether the genetic components improve the fit of the model significantly, we could carry out the command

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --restricted-no-var-genetic \
                   --out ./tutorial/restricted_novarG
```

where output file `restricted_novarG.loglik.out` reveals that allowing for genetic variance significantly improves the fit of the model:

```
Log-likelihood of nested model (null hypothesis) = -68327.76083331279,
based on data on 10 traits and 4980 observations,
with a model consisting of 1 genetic factors and 10 environment factors,
comprising 0 free genetic factor coefficients and 55 free environment factor coefficients in turn.
Estimates converged after 23 BFGS iterations 

Log-likelihood of parent model (alternative hypothesis) = -66849.63346504091,
based on data on 10 traits and 4980 observations,
with a model consisting of 10 genetic factors and 10 environment factors,
comprising 55 free genetic factor coefficients and 55 free environment factor coefficients in turn.
Estimates converged after 269 BFGS iterations 

Results of likelihood-ratio test with 55 degrees of freedom:
Chi-square test statistic is 2956.2547365437495
with P-value = 0.0
```

As before, `--restricted-no-var-genetic`, `--restricted-rho-genetic`, and/or `--restricted-genetic-model` cannot be combined with one another. Similarly, `--restricted-rho-environment` and `--restricted-environment-model` cannot be combined with each other.

:warning: **when using options such as `--restricted-genetic-model` and `--genetic-model`, it is your own responsibility to ensure the restricted model is nested with respect to the other model.** `mgreml` only inspects nestedness superficially. The best way to allow `mgreml` to assert nestedness is to appropriately label the factors in both models.

By default, `mgreml` will not store any intermediate results. However, using the `--store-iter` option, users can specify every how many iterations they want the current parameter estimates to be stored. E.g. `--store-iter 10` will cause `mgreml` to store estimates every ten iterations. The estimates will be stored in a so-called `.pkl` with a prefix a set by the `--out` option. This `.pkl` file contains the model specification as well as the estimates of that model in a given iteration.

Such a `.pkl` file can also be used to reinitialise `mgreml` e.g. if you accidentally switched off your computer halfway through an analysis. For instance

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --store-iter 50 \
                   --out ./tutorial/covar
```

causes `mgreml` to store results every 50 iterations. Then, if the preceding analysis has reached e.g. up until iteration 200 before a power outage, we could reinitialise later on using the following command:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --reinitialise ./tutorial/covar.estimates.iter.200.bfgs.pkl \
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
causes two sets of `.pkl` files to be stored (i.e. a file for every ten iterations, for both the restricted and alternative model) and
```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --reinitialise ./tutorial/restricted_rhoG1_rhoE0.estimates.iter.200.bfgs.pkl \
                   --restricted-reinitialise ./tutorial/restricted_rhoG1_rhoE0.estimates0.iter.20.bfgs.pkl \
                   --out ./tutorial/restricted_rhoG1_rhoE0_reinitialised
```
reinitialises estimation for the null and alternative model from appropriate `.pkl` files. Notice that analogous to `--reinitialise`, the `--restricted-reinitialise` option cannot be combined with options such as `--restricted-environment-model` and `--restricted-rho-genetic`, as the `.pkl` file already contains the full model specification.

`mgreml` performs basic data management, e.g. in terms of figuring out for which individuals we have phenotypic as well as GRM data (and data on covariates, if applicable). In case `--covar-model` is used `mgreml` also tests if there are any covariates that affect no phenotype at all, and if so, excludes such covariates.

`mgreml` also performs relatedness pruning when using the `--grm-cutoff` option. E.g. the following command selects a subset of individuals such that relatedness in the GRM is nowhere in excess of 0.05:

```
python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
                   --covar ./tutorial/covar.txt \
                   --grm-cutoff 0.05 \
                   --out ./tutorial/pruned
```

which in this case causes the sample size to decrease by 13, as in shown in part of the log-file below:

```
2. CLEANING YOUR DATA
[...]
APPLYING RELATEDNESS CUTOFF
Removing individuals such that there is no relatedness in excess of 0.05 in GRM
First pass: dropped 13 individuals with only one relatedness value in excess of 0.05
Second pass not required, as there are no individuals with relatedness in excess of 0.05 left after first pass
Relatedness cutoff has been applied
Remaining sample size is 4987
```

`mgreml` follows the greedy algorithm developed by Boppana and Halldórsson (1992); [doi:10.1007/BF01994876](https://link.springer.com/article/10.1007/BF01994876). Importantly, `mgreml` does this pruning at such a stage that sample size is maximised. E.g. for a pair of individuals with a relatedness in excess of the threshold, we try to keep the observation with the lowest missingness.

In general, `mgreml` simply tries to preserve sample size at each turn. E.g. if an individual has missing values only for a subset of the phenotypes, `mgreml` keeps that individual in the data, by introducing phenotype-by-individual-specific dummies (i.e. dummies that control for individual *i* having a missing value for trait *t*). Even when a covariate is missing, sometimes parts of that observation can still be salvaged (i.e. if the missing covariate does not affect all phenotypes according to `--covar-model`).

Introducing dummies to control for gaps in the data can become computationally highly demanding. Controlling for fixed-effect covariates has a computational complexity of the order *NT*<sup> 2</sup> provided the number of unique covariates is of the order 1. However, if e.g. missingness in each trait is a proportion of sample size, then the total number of unique covariates needed to control for this missingness becomes of the order *NT*, and thereby the computational complexity of controlling for this missingness of the order *N*<sup> 2</sup>*T*<sup> 3</sup>, which is prohibitively complex for large *N* and *T*.

Therefore, `mgreml` has a `--drop-missings` option, whereby all individuals are dropped that have at least one missing phenotype and/or at least one missing covariate that is relevant (either because `--covar-model` has not been used, or because the file following `--covar-model` indicates the covariate with a missing value for a given individual affects at least one trait).

Note that `mgreml` has a few advanced options regarding the estimation algorithm. First, `--newton` forces `mgreml` to use a Newton algorithm for solving the optimisation problem instead of BFGS. Although in theory this approach requires fewer iterations that BFGS, we strongly recommend sticking to BFGS: BFGS iterations are faster and BFGS is numerically much more stable, especially for large *T*.

In addition, `mgreml` considers the estimates to have converged if the root mean squared sum of the gradient vector (taking scale of traits and sample size into account) is below 10<sup>&minus;5</sup>. For conceptual ease, we refer to this criterion as length of the gradient. The option `--grad-tol` can be used to specify a different treshold. We do **NOT** recommend deviating from 10<sup>&minus;5</sup> by more than one order of magnitude. E.g. you could use `--grad-tol 5E-5` or `--grad-tol 1e-6`. However, we de **NOT** recommend e.g. `--grad-tol 1E-9`, as such a threshold requires a degree of convergence that is often beyond numerical precision, nor do we commend e.g. `--grad-tol 0.01`, as this is too lenient; the optimisation procedure simply has not converged when the latter convergence criterion is met.

Input files following the options `--pheno`, `--covar`, `--covar-model`, `--genetic-model`, `--environment-model`, `--restricted-genetic-model`, `--restricted-environment-model` can be comma, tab, or space-separated. Just make sure to be consistent within each file. Also, (1) when your labels (e.g. phenotype labels or FID and IIDs) contain spaces, make sure the file is not space-separated and (2) when those labels contain commas, make sure the file is not comma-separated.

Finally, that the options `--pheno`, `--covar`, `--covar-model`, `--genetic-model`, `--environment-model`, `--restricted-genetic-model`, `--restricted-environment-model` have modifiers `nolabelpheno`, `nolabelcovar`, and `nolabelfactor` to indicate when headers (e.g. phenotype labels) are absent. E.g. the following options are possible in `mgreml`: `--pheno myphen.txt nolabelpheno` and `--covar-model mycovmodel.txt nolabelpheno nolabelcovar`). However, as `mgreml` is a multivariate method, we strongly recommend always providing headers to `mgreml`, so everything is labelled in terms of input as well as output.

## Overview of commands

An overview of all `mgreml` commands is listed below:

| Command | Usage |
| --- | --- |
| -h, --help | show help message and exit |
| --grm PREFIX | prefix of binary GRM |
| --grm-cutoff THRESHOLD | option to drop individuals using greedy algorithm, <br> such that there is no relatedness in GRM in excess of <br> threshold for remaining individuals |
| --adjust-pcs INTEGER [INTEGER] | option to specify for how <br> many leading principal components (PCs) from genetic <br> data to adjust (to control for population <br> stratification) and for how many trailing PCs to <br> adjust (for computational efficiency); if just one <br> non-negative integer is specified this is taken as the <br> number of leading PCs to adjust for |
| --pheno FILENAME [nolabelpheno]  | phenotype file: <br> should be comma-, space-, or tab-separated, with one <br> row per individual, with FID and IID as first two <br> fields, followed by a field per phenotype; can be <br> followed by optional flag nolabelpheno, e.g. --pheno <br> mypheno.txt nolabelpheno, but we recommend to label <br> phenotypes |
| --drop-missings | option to drop all observations from data with at <br> least one missing phenotype or at least one missing <br> covariate |
| --no-intercept | option to indicate an intercept should not be included <br> automatically as covariate |
| --covar FILENAME [nolabelcovar]  | optional covariate file: should be comma-, space-, or tab- <br> separated, with one row per individual, with FID and <br> IID as first two fields, followed by a field per <br> covariate; can be followed by optional flag <br> nolabelcovar, e.g. --covar mycovar.txt nolabelcovar, <br> but we recommend to label covariates; WARNING: do not <br> include principal components from genetic data as <br> covariates, use --adjust-pcs instead |
| --covar-model FILENAME <br> [nolabelpheno] [nolabelcovar] | optional covariate model file: should be comma-, space-, or <br> tab-separated, with one row per phenotype and one <br> column per covariate; can be followed by optional <br> flags nolabelpheno and/or nolabelcovar, but we <br> recommend to label phenotypes and covariates; without <br> --covar-model, all covariates are assumed to apply to <br> all traits |
| --genetic-model FILENAME <br> [nolabelpheno] [nolabelfactor] |  optional genetic model file: should be comma-, space-, or tab- <br> separated, with one row per phenotype and one column <br> per genetic factor; can be followed by optional flags <br> nolabelpheno and/or nolabelfactor, but we recommend to <br> label phenotypes and genetic factors |
| --rho-genetic 0 or 1  | option followed by 0 or 1, forcing all genetic <br> correlations to take on the specified value; this flag <br> cannot be combined with --genetic-model |
| --no-var-genetic | option to force all genetic variances to equal zero; <br> this flag cannot be combined with --genetic-model <br> and/or --rho-genetic |
| --restricted-genetic-model FILENAME <br> [nolabelpheno] [nolabelfactor] |  optional restricted genetic model file: should be comma-, <br> space-, or tab-separated, with one row per phenotype <br> and one column per genetic factor; can be followed by <br> optional flags nolabelpheno and/or nolabelfactor, but <br> we recommend to label phenotypes and genetic factors |
| --restricted-rho-genetic 0 or 1 | option followed by 0 or 1, forcing all genetic <br> correlations in the restricted model to take on the <br> specified value; this flag cannot be combined with <br> --restricted-genetic-model |
| --restricted-no-var-genetic | option to force all genetic variances in the <br> restricted model to equal zero; this flag cannot be <br> combined with --restricted-genetic-model and/or <br> --restricted-rho-genetic |
| --environment-model FILENAME <br> [nolabelpheno] [nolabelfactor] |  optional environment model file: should be comma-, space-, or <br> tab-separated, with one row per phenotype and one <br> column per environment factor; can be followed by <br> optional flags nolabelpheno and/or nolabelfactor, but <br> we recommend to label phenotypes and environment <br> factors |
| --rho-environment 0  | option followed by 0, forcing all environment <br> correlations to zero; this flag cannot be combined <br> with --environment-model |
| --restricted-environment-model FILENAME <br> [nolabelpheno] [nolabelfactor] |  optional restricted environment model file: should be comma-, <br> space-, or tab-separated, with one row per phenotype <br> and one column per environment factor; can be followed <br> by optional flags nolabelpheno and/or nolabelfactor, <br> but we recommend to label phenotypes and environment <br> factors |
| --restricted-rho-environment 0 | option followed by 0, forcing all environment <br> correlations in the restricted model to zero; this <br> flag cannot be combined with --restricted-environment- <br> model |
| --no-se | option to skip calculation of standard errors and <br> covariance matrix of estimates |
| --factor-coefficients | option to report estimated factor coefficients |
| --variance-components | option to report estimated variance components |
| --newton | option to use Newton method instead of BFGS; not <br> recommended, unless the model is well-defined, <br> starting values are of good quality, and the number of <br> traits is small |
| --grad-tol THRESHOLD | option to set convergence threshold on the length of <br> the gradient vector per parameter, per observation, <br> different from the default value of 1E-5 |
| --store-iter INTEGER | option to specify every how many iterations you want <br> to store results |
| --reinitialise FILENAME | option to reinitialise mgreml for a model and its <br> estimates from a .pkl file stored by --store-iter |
| --restricted-reinitialise FILENAME | option to reinitialise mgreml for a restricted model <br> and its estimates from a .pkl file generated by <br> --store-iter |
| --out PREFIX |  prefix of output files |

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

## Citation

If you use the software, please cite

[]()

## License

This project is licensed under GNU GPL v3.

## Authors

Ronald de Vlaming (Vrije Universiteit Amsterdam)

Eric Slob (University of Cambridge)


