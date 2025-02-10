Mutation Accumulation Experiments
================

- [Installation](#installation)
- [Data format](#data-format)
- [Fitting the GCM model](#fitting-the-gcm-model)
  - [Accessing posterior samples](#accessing-posterior-samples)
  - [Various MCMC and model checks](#various-mcmc-and-model-checks)
- [Fitting the model with MMR
  saturation](#fitting-the-model-with-mmr-saturation)
  - [Accessing posterior samples](#accessing-posterior-samples-1)
  - [Various MCMC and model checks](#various-mcmc-and-model-checks-1)

[![DOI](https://zenodo.org/badge/923498484.svg)](https://doi.org/10.5281/zenodo.14847205)

<!-- README.md is generated from README.Rmd. Please edit that file -->

This is a package to analyse data generated from Mutation Accumulation
Experiments.

# Installation

Clone the repository on your computer using

    git clone git@forgemia.inra.fr:konkam/MutAccExperiments.git

or

    git clone https://forgemia.inra.fr/konkam/MutAccExperiments.git

You should then install the present package on your computer, using a
command such a the following, from inside the package folder:

    Rscript -e "devtools::install_local('.')"

Alternatively, you may open the folder using the Rstudio editor and
press `Ctrl + Alt + B`.

- Dependencies should be:
  - R package devtools
  - JAGS + the R package RJags (On ubuntu, jags is on the canonical
    repository, available by apt install jags)
  - R package tidyverse
  - R package coda

# Data format

``` r
library(MutAccExperiments)
library(tidyverse)
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.5
#> ✔ forcats   1.0.0     ✔ stringr   1.5.1
#> ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
#> ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
#> ✔ purrr     1.0.2     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

``` r
data(example_input_data)

input_data_onestrain <- example_input_data %>%
  filter(strain == first(strain))
```

The minimum data required to run the analysis is a data frame with the
following columns:

``` r
minimal_input_data_onestrain <- example_input_data %>%
  filter(strain == first(strain)) %>%
  select(mutation_id, m, n, t)
minimal_input_data_onestrain
#>   mutation_id   m       n      t
#> 1           1 118 1661176 251000
#> 2           2  23 1661176 251000
#> 3           3   8 1661176 251000
#> 4           4 120 2133850 251000
#> 5           5  25 2133850 251000
#> 6           6  25 2133850 251000
```

``` r
input_data_onestrain %>%
  check_input_format_GCM()
#> [1] "All good!"
```

# Fitting the GCM model

The GCM model is defined as follows:

$$
\begin{align}
m_i &\sim \text{Poisson}(n_i \mu_i t_i) \quad \forall i = 1, \ldots, I\\
\log_{10}(\mu_i) &\sim \mathcal{N}(\nu, \sigma) \\
\mu_0 &\sim \mathcal{N}(-8, 3) \\
\sigma &\sim \text{t}^+(0, 3, 5)\\
\end{align}
$$

where $m_i$ is the number of mutations observed of the $i$-th type,
$n_i$ is the number of possible mutation sites of type $i$ in the
genome, $t_i$ is the number of elapsed generations and $\mu_i$ is the
mutation rate for mutations of type $i$. $n_i, t_i$ are introduced so
that $\mu_i$ can be interpreted as a mutation rate per site and per
generation, which will be comparable over mutation types. $\nu$ is the
mean of the $\log10$ mutation rate, and $\sigma$ is the standard
deviation of the $\log10$ mutation rate. $\text{t}^+$ denotes the t
distribution truncated to positive values.

This is a hierarchical structure, where $\nu$ represents the average
tendency of the $\log10$ mutation rates (over mutation types), and
$\sigma$ represents the variability of the $\log10$ mutation rates
between mutation types. These two values caracterise the average
tendency and the dispersion of the mutation rates.

``` r
fit_GCM_model <- input_data_onestrain %>%
  EstimateMusGCM_onestrain()
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 6
#>    Unobserved stochastic nodes: 17
#>    Total graph size: 75
#> 
#> Initializing model
#> Warning: No initial value blocks found and n.chains not specified: 2 chains
#> were used
#> Warning: No initial values were provided - JAGS will use the same initial
#> values for all chains
#> 
#> Auto-run JAGS
#> 
#> Running a pilot chain...
#> Compiling rjags model...
#> Calling the simulation using the rjags method...
#> Adapting the model for 1000 iterations...
#> Burning in the model for 4000 iterations...
#> Running the model for 10000 iterations...
#> Simulation complete
#> Finished running the simulation
#> 
#> Calculating the Gelman-Rubin statistic for 18 variables....
#> The Gelman-Rubin statistic is below 1.05 for all parameters
#> 
#> Calculating the necessary sample length based on the Raftery and
#> Lewis's diagnostic...
#> Indicated sample length achieved
#> Auto-run JAGS complete
```

## Accessing posterior samples

``` r
extract_posterior_samples(fit_GCM_model, type = "mu")
#> # A tibble: 8,000 × 8
#>    `log10_mu[1]` `log10_mu[2]` `log10_mu[3]` `log10_mu[4]` `log10_mu[5]`
#>            <dbl>         <dbl>         <dbl>         <dbl>         <dbl>
#>  1         -9.57         -10.3         -10.2         -9.61         -10.4
#>  2         -9.55         -10.2         -10.8         -9.60         -10.4
#>  3         -9.52         -10.2         -10.6         -9.59         -10.3
#>  4         -9.56         -10.1         -10.4         -9.70         -10.3
#>  5         -9.51         -10.2         -10.4         -9.68         -10.4
#>  6         -9.53         -10.3         -10.4         -9.68         -10.2
#>  7         -9.56         -10.2         -10.8         -9.74         -10.2
#>  8         -9.57         -10.2         -10.6         -9.66         -10.3
#>  9         -9.56         -10.3         -10.6         -9.65         -10.2
#> 10         -9.53         -10.2         -10.7         -9.70         -10.4
#> # ℹ 7,990 more rows
#> # ℹ 3 more variables: `log10_mu[6]` <dbl>, iteration <int>, chain_id <int>
```

``` r
extract_posterior_samples(fit_GCM_model, type = "hyperparameters")
#> # A tibble: 8,000 × 5
#>    log10_sigma log10_mean loglikelihood iteration chain_id
#>          <dbl>      <dbl>         <dbl>     <int>    <int>
#>  1       0.522     -10.0          -25.2         1        1
#>  2       0.348      -9.83         -17.7         2        1
#>  3       0.432     -10.2          -19.3         3        1
#>  4       0.420      -9.83         -22.2         4        1
#>  5       0.576     -10.2          -20.2         5        1
#>  6       0.272     -10.0          -20.6         6        1
#>  7       0.423      -9.99         -23.3         7        1
#>  8       0.491      -9.86         -19.4         8        1
#>  9       0.407     -10.1          -18.6         9        1
#> 10       0.498      -9.99         -17.7        10        1
#> # ℹ 7,990 more rows
```

``` r
extract_posterior_samples(fit_GCM_model, type = "predictive")
#> # A tibble: 8,000 × 8
#>    `m_pred[1]` `m_pred[2]` `m_pred[3]` `m_pred[4]` `m_pred[5]` `m_pred[6]`
#>          <dbl>       <dbl>       <dbl>       <dbl>       <dbl>       <dbl>
#>  1         121          26          27         121          33          16
#>  2         118          23          12         131          21          23
#>  3         140          22           9         144          27          28
#>  4         111          42          17          94          30          29
#>  5         140          17          12         106          20          24
#>  6         133          23          19         111          24          15
#>  7         105          20          10          95          26          18
#>  8         106          30          11         118          24          29
#>  9         113          30          16         126          25          41
#> 10          98          24          10         109          21          30
#> # ℹ 7,990 more rows
#> # ℹ 2 more variables: iteration <int>, chain_id <int>
```

``` r
extract_posterior_samples(fit_GCM_model, type = "prior")
#> # A tibble: 8,000 × 5
#>    log10_mean_prior log10_sigma_prior log10_mu_prior iteration chain_id
#>               <dbl>             <dbl>          <dbl>     <int>    <int>
#>  1            -6.80            0.0977          -7.39         1        1
#>  2           -12.1             0.657          -11.8          2        1
#>  3            -9.62            0.319           -9.07         3        1
#>  4            -7.79            0.787           -7.05         4        1
#>  5           -11.1             0.306          -11.8          5        1
#>  6            -8.74            0.434           -8.57         6        1
#>  7            -9.94            0.223          -10.1          7        1
#>  8            -8.11            0.532           -7.91         8        1
#>  9            -7.12            0.791           -6.77         9        1
#> 10           -10.2             0.182           -9.48        10        1
#> # ℹ 7,990 more rows
```

Note that the m_pred variables are samples from the posterior predictive
distribution evaluated at the observed data points.

The variables with a suffix “\_prior” are samples from the prior
distribution.

You can also access the posterior samples at once:

``` r
extract_posterior_samples(fit_GCM_model)
#> # A tibble: 8,000 × 20
#>    iteration `m_pred[1]` `m_pred[2]` `m_pred[3]` `m_pred[4]` `m_pred[5]`
#>        <int>       <dbl>       <dbl>       <dbl>       <dbl>       <dbl>
#>  1         1         121          26          27         121          33
#>  2         2         118          23          12         131          21
#>  3         3         140          22           9         144          27
#>  4         4         111          42          17          94          30
#>  5         5         140          17          12         106          20
#>  6         6         133          23          19         111          24
#>  7         7         105          20          10          95          26
#>  8         8         106          30          11         118          24
#>  9         9         113          30          16         126          25
#> 10        10          98          24          10         109          21
#> # ℹ 7,990 more rows
#> # ℹ 14 more variables: `m_pred[6]` <dbl>, `log10_mu[1]` <dbl>,
#> #   `log10_mu[2]` <dbl>, `log10_mu[3]` <dbl>, `log10_mu[4]` <dbl>,
#> #   `log10_mu[5]` <dbl>, `log10_mu[6]` <dbl>, log10_mean <dbl>,
#> #   log10_sigma <dbl>, log10_mean_prior <dbl>, log10_sigma_prior <dbl>,
#> #   log10_mu_prior <dbl>, loglikelihood <dbl>, chain_id <int>
```

## Various MCMC and model checks

``` r
fit_GCM_model %>%
  traceplot()
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

``` r
fit_GCM_model %>%
  summary()
#>                         Lower95      Median     Upper95        Mean          SD
#> m_pred[1]          8.700000e+01 116.0000000 146.0000000 116.7613750 15.28428553
#> m_pred[2]          1.000000e+01  23.0000000  36.0000000  23.3825000  6.80310437
#> m_pred[3]          2.000000e+00   9.0000000  17.0000000   9.1460000  4.20635182
#> m_pred[4]          8.900000e+01 119.0000000 149.0000000 119.1158750 15.66496479
#> m_pred[5]          1.200000e+01  25.0000000  39.0000000  25.3816250  7.04387602
#> m_pred[6]          1.100000e+01  25.0000000  38.0000000  25.3207500  7.01961743
#> log10_mu[1]       -9.632650e+00  -9.5547774  -9.4754345  -9.5552007  0.04032195
#> log10_mu[2]       -1.043628e+01 -10.2589299 -10.0896339 -10.2619523  0.08898681
#> log10_mu[3]       -1.095701e+01 -10.6739188 -10.3958138 -10.6792255  0.14303977
#> log10_mu[4]       -9.734471e+00  -9.6542691  -9.5768387  -9.6551589  0.04011236
#> log10_mu[5]       -1.050082e+01 -10.3296148 -10.1662107 -10.3315513  0.08584925
#> log10_mu[6]       -1.051194e+01 -10.3295424 -10.1721960 -10.3331596  0.08613524
#> log10_mean        -1.056414e+01 -10.1216082  -9.6673078 -10.1214369  0.22596775
#> log10_sigma        2.480212e-01   0.4810637   0.8972228   0.5178723  0.18728680
#> log10_mean_prior  -1.366545e+01  -7.9796206  -2.0127306  -7.9810054  2.99045633
#> log10_sigma_prior  2.732636e-05   0.4193661   1.4807098   0.5462398  0.50145213
#> log10_mu_prior    -1.382192e+01  -7.9923695  -2.0611397  -7.9837055  3.03259120
#> loglikelihood     -2.259161e+01 -18.8171228 -16.4429524 -19.1273173  1.73656566
#>                   Mode        MCerr MC%ofSD SSeff        AC.20      psrf
#> m_pred[1]          115 0.1799253035     1.2  7216  0.006225207 0.9999415
#> m_pred[2]           22 0.0761119400     1.1  7989 -0.008541070 1.0004899
#> m_pred[3]            7 0.0486753186     1.2  7468 -0.020809646 1.0005447
#> m_pred[4]          120 0.1797516731     1.1  7595 -0.022659007 1.0005138
#> m_pred[5]           23 0.0799794947     1.1  7757  0.017771930 1.0006961
#> m_pred[6]           25 0.0814758101     1.2  7423  0.006154575 0.9998849
#> log10_mu[1]         NA 0.0004956962     1.2  6617  0.032270961 1.0005197
#> log10_mu[2]         NA 0.0010530311     1.2  7141 -0.014162985 0.9999156
#> log10_mu[3]         NA 0.0017582409     1.2  6618 -0.022070576 1.0003380
#> log10_mu[4]         NA 0.0004862265     1.2  6806 -0.023951466 1.0011036
#> log10_mu[5]         NA 0.0010344316     1.2  6888  0.002502742 1.0002648
#> log10_mu[6]         NA 0.0010515408     1.2  6710  0.004561343 1.0000277
#> log10_mean          NA 0.0024734151     1.1  8346 -0.003319086 1.0010346
#> log10_sigma         NA 0.0031833830     1.7  3461 -0.015268920 1.0018809
#> log10_mean_prior    NA 0.0334343182     1.1  8000  0.004211017 1.0004387
#> log10_sigma_prior   NA 0.0055372357     1.1  8201 -0.021312352 1.0007093
#> log10_mu_prior      NA 0.0334759109     1.1  8207  0.008075057 1.0004445
#> loglikelihood       NA 0.0227060206     1.3  5849  0.013536138 1.0000777
```

``` r
fit_GCM_model %>%
  plot_prior_posterior()
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

``` r
posterior_predictive_one_strain(fit_GCM_model)
#> Joining with `by = join_by(mutation_id)`
#> # A tibble: 6 × 17
#>   mutation_id m_pred_mean m_pred_median m_pred_infCI m_pred_supCI
#>         <int>       <dbl>         <dbl>        <dbl>        <dbl>
#> 1           1      117.             116           89          148
#> 2           2       23.4             23           12           38
#> 3           3        9.15             9            2           19
#> 4           4      119.             119           90          151
#> 5           5       25.4             25           13           41
#> 6           6       25.3             25           13           40
#> # ℹ 12 more variables: m_pred_infquart <dbl>, m_pred_supquart <dbl>,
#> #   genotype <chr>, mutation_label <chr>, nposinref <int>, ngeninMA <dbl>,
#> #   bps.n <int>, strain <chr>, context <chr>, m <int>, n <int>, t <dbl>
```

``` r
plot_posterior_predictive_onestrain(fit_GCM_model)
#> Joining with `by = join_by(mutation_id)`
#> Warning in geom_segment(aes(y = m_pred_infCI/(n * t), yend = m_pred_supCI/(n *
#> : Ignoring unknown parameters: `width`
#> Warning in geom_segment(aes(y = m_pred_infquart/(n * t), yend =
#> m_pred_supquart/(n * : Ignoring unknown parameters: `width`
```

<img src="man/figures/README-unnamed-chunk-16-1.png" width="100%" />

# Fitting the model with MMR saturation

Compared to the previous unconstrained model, the MMR saturation model
assumes that the mutation rates of the different strains are related,
i.e. constrained by each other. Four strains are considered, wild-type
(wt), MMR-deactivated (MMR-), proofreading-deficient (proof-) and
MMR-deactivated and proofreading-deficient (MMR-\_proof-). The model
includes an explicit parameter for the saturation mechanism, which would
be equal to 0 in the absence of saturation. Therefore, this parameter
quantifies both the evidence for the existence of a saturation mechanism
and its magnitude.

n practice, for each mutation type we assume that there is a baseline
error rate $\gamma_i$. These mutations may be corrected by the
Proofreading mechanism with a probability $1-q_{\text{proofreading},i}$,
and further corrected by the MMR mechanism with a probability
$1-q_{\text{MMR},i}$. Therefore, in the wild-type strain, the number of
observed mutation is impacted by the error rate and the two repair
mechanisms. This error rate is observed directly for the mutant
MMR-\_proof- for which both MMR and Proofreading are inactivated. The
strain MMR- has the baseline error rate moderated by the Proofreading
mechanism, the strain C\* has the baseline error moderated by the MMR
mechanism which can saturate. Saturation is modelled as follows: a
proportion $1-\theta$ of mutations are corrected by the MMR mechanism
with probability $1-q_{\text{MMR},i}$, until saturation may occur. At
this point, the remaining proportion of mutations remain uncorrected. As
mentioned above, if $\theta=0$, there is no saturation, and an estimate
of $\theta$ significatively different from 0 provides evidence for
saturation. This results in the following system of equations for the
rate of mutation of type $i$ per site per generation in the four
strains:

$$
\begin{align}
\mu_{\text{MMR-\_proof-},i} &= \gamma_i \\
\mu_{\text{MMR-},i} &= \gamma_i q_{\text{proofreading},i} \\
\mu_{\text{proof-},i} &= \gamma_i \left(\theta + (1-\theta) q_{\text{MMR},i} \right) \\
\mu_{\text{wt},i} &= \gamma_i q_{\text{proofreading},i}q_{\text{MMR},i} \\
\end{align}
$$

We choose the following generative model:

$$
\begin{align}
m_{s, i} &\sim \text{Poisson}(n_{s, i} \mu_{s, i} t_{s, i}) \quad \text{for } i = 1, \ldots, I, s \in \left\{\text{wt}, \text{MMR-}, \text{proof-}, \text{MMR-proof-}\right\}\\
\log10 \mu_{s, i} &\sim \mathcal{N}(\nu_{s}, \sigma_{s})  \quad \text{for } i = 1, \ldots, I, s \in \left\{\text{wt}, \text{MMR-}, \text{MMR-proof-}\right\} \\
\nu_s &\sim \mathcal{N}(-8, 3) \quad \text{for } s \in \left\{\text{wt}, \text{MMR-}, \text{MMR-proof-}\right\}\\
\sigma_s &\sim \text{t}^+(0, 3, 5) \quad \text{for } s \in \left\{\text{wt}, \text{MMR-}, \text{MMR-proof-}\right\}\\
\theta &\sim \text{Beta}(1, 1)\\
\end{align}
$$

``` r
minimal_input_data <- example_input_data %>%
  select(strain, mutation_id, m, n, t)

minimal_input_data
#>           strain mutation_id    m       n         t
#> 1             wt           1  118 1661176 251000.00
#> 2             wt           2   23 1661176 251000.00
#> 3             wt           3    8 1661176 251000.00
#> 4             wt           4  120 2133850 251000.00
#> 5             wt           5   25 2133850 251000.00
#> 6             wt           6   25 2133850 251000.00
#> 7       MMRminus           1 2419 1661176  38000.00
#> 8       MMRminus           2   14 1661176  38000.00
#> 9       MMRminus           3   21 1661176  38000.00
#> 10      MMRminus           4 2292 2133850  38000.00
#> 11      MMRminus           5   58 2133850  38000.00
#> 12      MMRminus           6   40 2133850  38000.00
#> 13    proofminus           1  210 1661176   1895.14
#> 14    proofminus           2   19 1661176   1895.14
#> 15    proofminus           3    2 1661176   1895.14
#> 16    proofminus           4  138 2133850   1895.14
#> 17    proofminus           5   13 2133850   1895.14
#> 18    proofminus           6   13 2133850   1895.14
#> 19 MMRproofminus           1  274 1661176    230.49
#> 20 MMRproofminus           2    2 1661176    230.49
#> 21 MMRproofminus           3    2 1661176    230.49
#> 22 MMRproofminus           4  210 2133850    230.49
#> 23 MMRproofminus           5   12 2133850    230.49
#> 24 MMRproofminus           6    2 2133850    230.49
```

``` r
fit_MMRsaturation_model <- EstimateMusMMRsaturation(example_input_data)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 24
#>    Unobserved stochastic nodes: 53
#>    Total graph size: 232
#> 
#> Initializing model
#> Warning: No initial value blocks found and n.chains not specified: 2 chains
#> were used
#> Warning: No initial values were provided - JAGS will use the same initial
#> values for all chains
#> 
#> Auto-run JAGS
#> 
#> Running a pilot chain...
#> Compiling rjags model...
#> Calling the simulation using the rjags method...
#> Adapting the model for 1000 iterations...
#> Burning in the model for 4000 iterations...
#> Running the model for 10000 iterations...
#> Simulation complete
#> Finished running the simulation
#> 
#> Calculating the Gelman-Rubin statistic for 60 variables....
#> The Gelman-Rubin statistic is below 1.05 for all parameters
#> 
#> Calculating the necessary sample length based on the Raftery and
#> Lewis's diagnostic...
#> The model will need to be run for a further 2728 updates.  This will
#> take approximately 0.2 seconds.
#> 
#> Calling the simulation using the rjags method...
#> Note: the model did not require adaptation
#> Running the model for 2728 iterations...
#> Simulation complete
#> Finished running the simulation
#> Indicated sample length achieved
#> Note: Summary statistics were not produced as there are >50 monitored
#> variables
#> [To override this behaviour see ?add.summary and ?runjags.options]
#> FALSEAuto-run JAGS complete
```

## Accessing posterior samples

``` r
extract_posterior_samples(fit_MMRsaturation_model, type = "mu")
#> # A tibble: 8,000 × 26
#>    `mu_wt[1]` `mu_wt[2]` `mu_wt[3]` `mu_wt[4]` `mu_wt[5]` `mu_wt[6]`
#>         <dbl>      <dbl>      <dbl>      <dbl>      <dbl>      <dbl>
#>  1   2.50e-10   6.12e-11   1.50e-11   2.21e-10   4.84e-11   5.17e-11
#>  2   3.03e-10   5.38e-11   3.02e-11   2.16e-10   4.13e-11   5.08e-11
#>  3   2.86e-10   6.37e-11   2.48e-11   2.04e-10   4.08e-11   4.82e-11
#>  4   2.89e-10   5.53e-11   2.15e-11   2.17e-10   4.28e-11   3.64e-11
#>  5   2.80e-10   4.55e-11   2.38e-11   2.03e-10   4.36e-11   5.09e-11
#>  6   2.89e-10   5.40e-11   3.37e-11   2.27e-10   6.77e-11   4.77e-11
#>  7   2.48e-10   6.08e-11   3.08e-11   2.27e-10   6.48e-11   5.03e-11
#>  8   2.76e-10   7.16e-11   3.73e-11   2.33e-10   4.03e-11   5.56e-11
#>  9   2.77e-10   5.27e-11   3.05e-11   1.88e-10   5.90e-11   4.43e-11
#> 10   2.62e-10   5.51e-11   2.15e-11   2.35e-10   5.91e-11   3.87e-11
#> # ℹ 7,990 more rows
#> # ℹ 20 more variables: `mu_proofminus[1]` <dbl>, `mu_proofminus[2]` <dbl>,
#> #   `mu_proofminus[3]` <dbl>, `mu_proofminus[4]` <dbl>,
#> #   `mu_proofminus[5]` <dbl>, `mu_proofminus[6]` <dbl>, `mu_MMRminus[1]` <dbl>,
#> #   `mu_MMRminus[2]` <dbl>, `mu_MMRminus[3]` <dbl>, `mu_MMRminus[4]` <dbl>,
#> #   `mu_MMRminus[5]` <dbl>, `mu_MMRminus[6]` <dbl>,
#> #   `mu_MMRproofminus[1]` <dbl>, `mu_MMRproofminus[2]` <dbl>, …
```

``` r
extract_posterior_samples(fit_MMRsaturation_model, type = "hyperparameters")
#> # A tibble: 8,000 × 10
#>    log10_sigma_wt log10_sigma_MMRminus log10_sigma_MMRproofminus log10_mean_wt
#>             <dbl>                <dbl>                     <dbl>         <dbl>
#>  1          0.403                1.23                      0.486        -10.2 
#>  2          0.378                0.886                     0.980        -10.2 
#>  3          0.336                1.09                      0.702        -10.4 
#>  4          0.467                1.35                      0.766        -10.1 
#>  5          0.312                1.18                      0.796        -10.3 
#>  6          0.420                0.907                     0.798         -9.95
#>  7          0.424                1.19                      0.684        -10.0 
#>  8          0.355                0.897                     0.651        -10.2 
#>  9          0.442                1.20                      0.850         -9.96
#> 10          0.506                1.00                      0.957        -10.2 
#> # ℹ 7,990 more rows
#> # ℹ 6 more variables: log10_mean_MMRminus <dbl>,
#> #   log10_mean_MMRproofminus <dbl>, theta4 <dbl>, loglikelihood <dbl>,
#> #   iteration <int>, chain_id <int>
```

``` r
extract_posterior_samples(fit_MMRsaturation_model, type = "predictive")
#> # A tibble: 8,000 × 26
#>    `m_wt_pred[1]` `m_wt_pred[2]` `m_wt_pred[3]` `m_wt_pred[4]` `m_wt_pred[5]`
#>             <dbl>          <dbl>          <dbl>          <dbl>          <dbl>
#>  1            118             29              9            145             28
#>  2            117             30             14            104             19
#>  3            128             35             10            110             31
#>  4            116             30              8            140             23
#>  5            105             22             12            103             17
#>  6            128             24             12            137             29
#>  7            116             33             20            124             35
#>  8            135             25             22            122             20
#>  9            140             22              9            120             40
#> 10             95             19             11            110             37
#> # ℹ 7,990 more rows
#> # ℹ 21 more variables: `m_wt_pred[6]` <dbl>, `m_proofminus_pred[1]` <dbl>,
#> #   `m_proofminus_pred[2]` <dbl>, `m_proofminus_pred[3]` <dbl>,
#> #   `m_proofminus_pred[4]` <dbl>, `m_proofminus_pred[5]` <dbl>,
#> #   `m_proofminus_pred[6]` <dbl>, `m_MMRminus_pred[1]` <dbl>,
#> #   `m_MMRminus_pred[2]` <dbl>, `m_MMRminus_pred[3]` <dbl>,
#> #   `m_MMRminus_pred[4]` <dbl>, `m_MMRminus_pred[5]` <dbl>, …
```

Note that the m_pred variables are samples from the posterior predictive
distribution evaluated at the observed data points.

The variables with a suffix “\_prior” are samples from the prior
distribution.

You can also access the posterior samples at once:

``` r
extract_posterior_samples(fit_MMRsaturation_model)
#> # A tibble: 8,000 × 62
#>    iteration `mu_wt[1]` `mu_wt[2]` `mu_wt[3]` `mu_wt[4]` `mu_wt[5]` `mu_wt[6]`
#>        <int>      <dbl>      <dbl>      <dbl>      <dbl>      <dbl>      <dbl>
#>  1         1   2.50e-10   6.12e-11   1.50e-11   2.21e-10   4.84e-11   5.17e-11
#>  2         2   3.03e-10   5.38e-11   3.02e-11   2.16e-10   4.13e-11   5.08e-11
#>  3         3   2.86e-10   6.37e-11   2.48e-11   2.04e-10   4.08e-11   4.82e-11
#>  4         4   2.89e-10   5.53e-11   2.15e-11   2.17e-10   4.28e-11   3.64e-11
#>  5         5   2.80e-10   4.55e-11   2.38e-11   2.03e-10   4.36e-11   5.09e-11
#>  6         6   2.89e-10   5.40e-11   3.37e-11   2.27e-10   6.77e-11   4.77e-11
#>  7         7   2.48e-10   6.08e-11   3.08e-11   2.27e-10   6.48e-11   5.03e-11
#>  8         8   2.76e-10   7.16e-11   3.73e-11   2.33e-10   4.03e-11   5.56e-11
#>  9         9   2.77e-10   5.27e-11   3.05e-11   1.88e-10   5.90e-11   4.43e-11
#> 10        10   2.62e-10   5.51e-11   2.15e-11   2.35e-10   5.91e-11   3.87e-11
#> # ℹ 7,990 more rows
#> # ℹ 55 more variables: `mu_proofminus[1]` <dbl>, `mu_proofminus[2]` <dbl>,
#> #   `mu_proofminus[3]` <dbl>, `mu_proofminus[4]` <dbl>,
#> #   `mu_proofminus[5]` <dbl>, `mu_proofminus[6]` <dbl>, `mu_MMRminus[1]` <dbl>,
#> #   `mu_MMRminus[2]` <dbl>, `mu_MMRminus[3]` <dbl>, `mu_MMRminus[4]` <dbl>,
#> #   `mu_MMRminus[5]` <dbl>, `mu_MMRminus[6]` <dbl>,
#> #   `mu_MMRproofminus[1]` <dbl>, `mu_MMRproofminus[2]` <dbl>, …
```

## Various MCMC and model checks

``` r
fit_MMRsaturation_model %>%
  traceplot()
```

<img src="man/figures/README-unnamed-chunk-23-1.png" width="100%" />

``` r
fit_MMRsaturation_model %>%
  summary()
#> Calculating summary statistics...
#> Calculating the Gelman-Rubin statistic for 60 variables....
#>                                 Lower95        Median       Upper95
#> mu_wt[1]                   2.316088e-10  2.801983e-10  3.327877e-10
#> mu_wt[2]                   3.898391e-11  6.017714e-11  8.333290e-11
#> mu_wt[3]                   8.335881e-12  2.094019e-11  3.491409e-11
#> mu_wt[4]                   1.839339e-10  2.202054e-10  2.632560e-10
#> mu_wt[5]                   3.085319e-11  4.665977e-11  6.523821e-11
#> mu_wt[6]                   3.318521e-11  5.053250e-11  7.049253e-11
#> mu_proofminus[1]           5.703821e-08  6.467118e-08  7.233851e-08
#> mu_proofminus[2]           3.036368e-09  5.003104e-09  7.492941e-09
#> mu_proofminus[3]           1.795143e-10  7.200352e-10  1.508074e-09
#> mu_proofminus[4]           3.180374e-08  3.630132e-08  4.147641e-08
#> mu_proofminus[5]           2.092919e-09  3.304847e-09  4.697447e-09
#> mu_proofminus[6]           1.168381e-09  2.216073e-09  3.429724e-09
#> mu_MMRminus[1]             3.675916e-08  3.828836e-08  3.977777e-08
#> mu_MMRminus[2]             1.004744e-10  1.897259e-10  2.984462e-10
#> mu_MMRminus[3]             2.097424e-10  3.364727e-10  4.878875e-10
#> mu_MMRminus[4]             2.718252e-08  2.826871e-08  2.944247e-08
#> mu_MMRminus[5]             5.501376e-10  7.174993e-10  9.123567e-10
#> mu_MMRminus[6]             3.297614e-10  4.681641e-10  6.227534e-10
#> mu_MMRproofminus[1]        6.510188e-07  7.283817e-07  8.025818e-07
#> mu_MMRproofminus[2]        6.133720e-09  1.338109e-08  2.100686e-08
#> mu_MMRproofminus[3]        1.276385e-09  5.119222e-09  1.068294e-08
#> mu_MMRproofminus[4]        3.588760e-07  4.069356e-07  4.572664e-07
#> mu_MMRproofminus[5]        1.472130e-08  2.306665e-08  3.338705e-08
#> mu_MMRproofminus[6]        6.226581e-09  1.217048e-08  1.903521e-08
#> log10_mean_wt             -1.058113e+01 -1.011224e+01 -9.688167e+00
#> log10_mean_MMRminus       -9.668356e+00 -8.764848e+00 -7.919508e+00
#> log10_mean_MMRproofminus  -8.178133e+00 -7.385191e+00 -6.624328e+00
#> log10_sigma_wt             2.461381e-01  4.774903e-01  8.881818e-01
#> log10_sigma_MMRminus       5.271417e-01  9.621457e-01  1.586757e+00
#> log10_sigma_MMRproofminus  4.911006e-01  8.677894e-01  1.454264e+00
#> theta4                     7.008305e-02  8.202657e-02  9.447116e-02
#> log10_mean_prior          -1.392868e+01 -7.908641e+00 -2.208318e+00
#> log10_sigma_prior          1.291041e-04  4.115669e-01  1.451936e+00
#> log10_mu_prior            -1.407854e+01 -7.910198e+00 -1.790458e+00
#> theta4_prior               3.913924e-02  4.990964e-01  9.877131e-01
#> m_wt_pred[1]               8.800000e+01  1.160000e+02  1.470000e+02
#> m_wt_pred[2]               1.300000e+01  2.500000e+01  3.900000e+01
#> m_wt_pred[3]               2.000000e+00  8.000000e+00  1.700000e+01
#> m_wt_pred[4]               8.600000e+01  1.180000e+02  1.460000e+02
#> m_wt_pred[5]               1.100000e+01  2.500000e+01  3.800000e+01
#> m_wt_pred[6]               1.300000e+01  2.700000e+01  4.100000e+01
#> m_proofminus_pred[1]       1.640000e+02  2.030000e+02  2.380000e+02
#> m_proofminus_pred[2]       5.000000e+00  1.600000e+01  2.600000e+01
#> m_proofminus_pred[3]       0.000000e+00  2.000000e+00  6.000000e+00
#> m_proofminus_pred[4]       1.170000e+02  1.460000e+02  1.770000e+02
#> m_proofminus_pred[5]       5.000000e+00  1.300000e+01  2.200000e+01
#> m_proofminus_pred[6]       2.000000e+00  9.000000e+00  1.600000e+01
#> m_MMRminus_pred[1]         2.280000e+03  2.417000e+03  2.546000e+03
#> m_MMRminus_pred[2]         3.000000e+00  1.200000e+01  2.100000e+01
#> m_MMRminus_pred[3]         1.000000e+01  2.100000e+01  3.400000e+01
#> m_MMRminus_pred[4]         2.162000e+03  2.293000e+03  2.427000e+03
#> m_MMRminus_pred[5]         3.600000e+01  5.800000e+01  7.800000e+01
#> m_MMRminus_pred[6]         2.200000e+01  3.800000e+01  5.500000e+01
#> m_MMRproofminus_pred[1]    2.330000e+02  2.790000e+02  3.200000e+02
#> m_MMRproofminus_pred[2]    0.000000e+00  5.000000e+00  1.000000e+01
#> m_MMRproofminus_pred[3]    0.000000e+00  2.000000e+00  5.000000e+00
#> m_MMRproofminus_pred[4]    1.630000e+02  2.000000e+02  2.350000e+02
#> m_MMRproofminus_pred[5]    3.000000e+00  1.100000e+01  1.900000e+01
#> m_MMRproofminus_pred[6]    1.000000e+00  6.000000e+00  1.200000e+01
#> loglikelihood             -1.113369e+01 -1.001761e+01 -9.133831e+00
#>                                    Mean           SD         Mode        MCerr
#> mu_wt[1]                   2.811626e-10 2.596957e-11 1.904657e-10          Inf
#> mu_wt[2]                   6.092464e-11 1.151680e-11 2.720260e-11          Inf
#> mu_wt[3]                   2.158322e-11 7.020899e-12 5.346114e-12          Inf
#> mu_wt[4]                   2.208087e-10 2.027553e-11 1.550353e-10          Inf
#> mu_wt[5]                   4.722159e-11 8.985195e-12 2.233567e-11          Inf
#> mu_wt[6]                   5.108709e-11 9.693198e-12 2.274875e-11          Inf
#> mu_proofminus[1]           6.475059e-08 3.964205e-09           NA          Inf
#> mu_proofminus[2]           5.091549e-09 1.156863e-09 2.016002e-09          Inf
#> mu_proofminus[3]           7.761157e-10 3.636149e-10 5.990096e-11          Inf
#> mu_proofminus[4]           3.636128e-08 2.445917e-09           NA          Inf
#> mu_proofminus[5]           3.344112e-09 6.814655e-10 1.443905e-09          Inf
#> mu_proofminus[6]           2.265673e-09 5.897461e-10 7.029891e-10          Inf
#> mu_MMRminus[1]             3.829764e-08 7.777372e-10           NA          Inf
#> mu_MMRminus[2]             1.951918e-10 5.231363e-11 6.416635e-11          Inf
#> mu_MMRminus[3]             3.418277e-10 7.166295e-11 1.472561e-10          Inf
#> mu_MMRminus[4]             2.827443e-08 5.788359e-10           NA          Inf
#> mu_MMRminus[5]             7.208175e-10 9.355002e-11 4.204445e-10          Inf
#> mu_MMRminus[6]             4.727200e-10 7.649947e-11 2.337053e-10          Inf
#> mu_MMRproofminus[1]        7.292387e-07 3.876777e-08           NA 5.591735e-10
#> mu_MMRproofminus[2]        1.370819e-08 3.896984e-09           NA          Inf
#> mu_MMRproofminus[3]        5.514409e-09 2.589636e-09           NA          Inf
#> mu_MMRproofminus[4]        4.075846e-07 2.525548e-08           NA 3.378651e-10
#> mu_MMRproofminus[5]        2.346389e-08 4.833182e-09           NA          Inf
#> mu_MMRproofminus[6]        1.246873e-08 3.383658e-09           NA          Inf
#> log10_mean_wt             -1.011299e+01 2.247981e-01           NA 2.522805e-03
#> log10_mean_MMRminus       -8.764206e+00 4.341214e-01           NA 4.705666e-03
#> log10_mean_MMRproofminus  -7.390719e+00 3.885012e-01           NA 4.343575e-03
#> log10_sigma_wt             5.158166e-01 1.829345e-01           NA 2.545003e-03
#> log10_sigma_MMRminus       1.015121e+00 2.926597e-01           NA 4.030640e-03
#> log10_sigma_MMRproofminus  9.165807e-01 2.763054e-01           NA 3.778979e-03
#> theta4                     8.222722e-02 6.270731e-03           NA 8.867577e-05
#> log10_mean_prior          -7.950648e+00 3.021411e+00           NA 3.378041e-02
#> log10_sigma_prior          5.358693e-01 5.005591e-01           NA 5.689241e-03
#> log10_mu_prior            -7.956450e+00 3.118144e+00           NA 3.486191e-02
#> theta4_prior               4.993872e-01 2.903662e-01           NA 3.273800e-03
#> m_wt_pred[1]               1.172157e+02 1.537047e+01 1.120000e+02 1.718470e-01
#> m_wt_pred[2]               2.540737e+01 6.927151e+00 2.500000e+01 8.295390e-02
#> m_wt_pred[3]               8.979375e+00 4.144285e+00 7.000000e+00 4.960976e-02
#> m_wt_pred[4]               1.183599e+02 1.542310e+01 1.180000e+02 1.607268e-01
#> m_wt_pred[5]               2.536100e+01 7.046480e+00 2.200000e+01 7.878205e-02
#> m_wt_pred[6]               2.744012e+01 7.384543e+00 2.700000e+01 8.137997e-02
#> m_proofminus_pred[1]       2.035948e+02 1.890116e+01 2.000000e+02 2.113213e-01
#> m_proofminus_pred[2]       1.600425e+01 5.422329e+00 1.600000e+01 6.062349e-02
#> m_proofminus_pred[3]       2.473750e+00 1.954044e+00 1.000000e+00 2.184688e-02
#> m_proofminus_pred[4]       1.468876e+02 1.560449e+01 1.430000e+02 1.744635e-01
#> m_proofminus_pred[5]       1.354700e+01 4.591288e+00 1.200000e+01 5.017839e-02
#> m_proofminus_pred[6]       9.175125e+00 3.843443e+00 9.000000e+00 4.297100e-02
#> m_MMRminus_pred[1]         2.416889e+03 6.943973e+01 2.417000e+03 7.763597e-01
#> m_MMRminus_pred[2]         1.226125e+01 4.801159e+00 1.200000e+01 5.775676e-02
#> m_MMRminus_pred[3]         2.159287e+01 6.493277e+00 2.300000e+01 7.259704e-02
#> m_MMRminus_pred[4]         2.293229e+03 6.775671e+01 2.303000e+03 7.575430e-01
#> m_MMRminus_pred[5]         5.857450e+01 1.084004e+01 5.900000e+01 1.228366e-01
#> m_MMRminus_pred[6]         3.834375e+01 8.761060e+00 3.500000e+01 9.981655e-02
#> m_MMRproofminus_pred[1]    2.792479e+02 2.239846e+01 2.760000e+02 2.738882e-01
#> m_MMRproofminus_pred[2]    5.218125e+00 2.731662e+00 5.000000e+00 3.288208e-02
#> m_MMRproofminus_pred[3]    2.098625e+00 1.764309e+00 1.000000e+00 2.015765e-02
#> m_MMRproofminus_pred[4]    2.004787e+02 1.875507e+01 1.980000e+02 2.262012e-01
#> m_MMRproofminus_pred[5]    1.151538e+01 4.155980e+00 1.100000e+01 4.646527e-02
#> m_MMRproofminus_pred[6]    6.129375e+00 2.954273e+00 6.000000e+00 3.435409e-02
#> loglikelihood             -1.007444e+01 5.211115e-01           NA 6.256127e-03
#>                           MC%ofSD SSeff         AC.30      psrf
#> mu_wt[1]                      Inf     0  0.0091310799 1.0001852
#> mu_wt[2]                      Inf     0 -0.0014507455 1.0009218
#> mu_wt[3]                      Inf     0  0.0006945708 1.0002506
#> mu_wt[4]                      Inf     0 -0.0129736621 1.0007723
#> mu_wt[5]                      Inf     0  0.0094102515 0.9999170
#> mu_wt[6]                      Inf     0  0.0149162562 0.9999584
#> mu_proofminus[1]              Inf     0  0.0052484494 1.0005388
#> mu_proofminus[2]              Inf     0  0.0035031193 1.0006089
#> mu_proofminus[3]              Inf     0 -0.0151714812 1.0005327
#> mu_proofminus[4]              Inf     0  0.0031127759 1.0001394
#> mu_proofminus[5]              Inf     0  0.0059667455 1.0001627
#> mu_proofminus[6]              Inf     0 -0.0034571894 1.0004701
#> mu_MMRminus[1]                Inf     0 -0.0037430708 0.9999483
#> mu_MMRminus[2]                Inf     0  0.0315998799 1.0007269
#> mu_MMRminus[3]                Inf     0 -0.0270617555 0.9998799
#> mu_MMRminus[4]                Inf     0 -0.0057899435 1.0002554
#> mu_MMRminus[5]                Inf     0  0.0125623565 1.0001456
#> mu_MMRminus[6]                Inf     0  0.0103605818 1.0001869
#> mu_MMRproofminus[1]           1.4  4807  0.0157584154 1.0001210
#> mu_MMRproofminus[2]           Inf     0  0.0349208031 1.0007020
#> mu_MMRproofminus[3]           Inf     0 -0.0063451462 1.0002215
#> mu_MMRproofminus[4]           1.3  5588  0.0234519690 1.0017567
#> mu_MMRproofminus[5]           Inf     0  0.0082710082 0.9999679
#> mu_MMRproofminus[6]           Inf     0  0.0020787535 1.0002439
#> log10_mean_wt                 1.1  7940 -0.0076176149 0.9999511
#> log10_mean_MMRminus           1.1  8511 -0.0078785186 1.0000166
#> log10_mean_MMRproofminus      1.1  8000  0.0053304740 0.9998866
#> log10_sigma_wt                1.4  5167 -0.0216240597 1.0016999
#> log10_sigma_MMRminus          1.4  5272 -0.0101332910 1.0000294
#> log10_sigma_MMRproofminus     1.4  5346  0.0049126841 1.0000627
#> theta4                        1.4  5001  0.0068933585 1.0008303
#> log10_mean_prior              1.1  8000 -0.0116836903 0.9998998
#> log10_sigma_prior             1.1  7741  0.0115503007 0.9999157
#> log10_mu_prior                1.1  8000 -0.0138759284 0.9998792
#> theta4_prior                  1.1  7867  0.0031666844 1.0001105
#> m_wt_pred[1]                  1.1  8000  0.0040211845 0.9998811
#> m_wt_pred[2]                  1.2  6973 -0.0077195012 1.0009552
#> m_wt_pred[3]                  1.2  6979 -0.0040193247 1.0000022
#> m_wt_pred[4]                  1.0  9208 -0.0164073442 1.0005304
#> m_wt_pred[5]                  1.1  8000  0.0098121386 0.9999135
#> m_wt_pred[6]                  1.1  8234 -0.0018382057 1.0005319
#> m_proofminus_pred[1]          1.1  8000 -0.0174230996 0.9999204
#> m_proofminus_pred[2]          1.1  8000  0.0103787923 1.0000379
#> m_proofminus_pred[3]          1.1  8000 -0.0209070146 1.0008440
#> m_proofminus_pred[4]          1.1  8000 -0.0016446548 0.9999191
#> m_proofminus_pred[5]          1.1  8372 -0.0113564949 1.0015305
#> m_proofminus_pred[6]          1.1  8000 -0.0100479145 1.0002220
#> m_MMRminus_pred[1]            1.1  8000 -0.0097128848 1.0000149
#> m_MMRminus_pred[2]            1.2  6910  0.0158663900 1.0002873
#> m_MMRminus_pred[3]            1.1  8000 -0.0121568774 1.0000158
#> m_MMRminus_pred[4]            1.1  8000  0.0004401771 1.0005308
#> m_MMRminus_pred[5]            1.1  7788  0.0019972728 1.0009042
#> m_MMRminus_pred[6]            1.1  7704  0.0150440256 0.9998828
#> m_MMRproofminus_pred[1]       1.2  6688  0.0046887302 0.9999511
#> m_MMRproofminus_pred[2]       1.2  6901  0.0107137758 1.0002882
#> m_MMRproofminus_pred[3]       1.1  7661  0.0070651982 1.0001058
#> m_MMRproofminus_pred[4]       1.2  6875 -0.0017587473 1.0019364
#> m_MMRproofminus_pred[5]       1.1  8000 -0.0051705079 1.0000921
#> m_MMRproofminus_pred[6]       1.2  7395 -0.0006199068 0.9999879
#> loglikelihood                 1.2  6938  0.0146270865 1.0000509
```

``` r
fit_MMRsaturation_model %>%
  plot_prior_posterior()
```

<img src="man/figures/README-unnamed-chunk-25-1.png" width="100%" />

``` r
fit_MMRsaturation_model %>%
  posterior_predictive()
#> Joining with `by = join_by(strain, mutation_id)`
#> # A tibble: 24 × 18
#>    m_pred_mean m_pred_median m_pred_infCI m_pred_supCI m_pred_infquart
#>          <dbl>         <dbl>        <dbl>        <dbl>           <dbl>
#>  1      117.             116           88          149             107
#>  2       25.4             25           13           40              21
#>  3        8.98             8            2           18               6
#>  4      118.             118           90          150             108
#>  5       25.4             25           13           40              20
#>  6       27.4             27           14           43              22
#>  7      204.             203          167          241             191
#>  8       16.0             16            7           28              12
#>  9        2.47             2            0            7               1
#> 10      147.             146          117          178             136
#> # ℹ 14 more rows
#> # ℹ 13 more variables: m_pred_supquart <dbl>, colnam <chr>, strain <chr>,
#> #   mutation_id <dbl>, genotype <chr>, mutation_label <chr>, nposinref <int>,
#> #   ngeninMA <dbl>, bps.n <int>, context <chr>, m <int>, n <int>, t <dbl>
```

``` r
fit_MMRsaturation_model %>%
  plot_posterior_predictive()
#> Joining with `by = join_by(strain, mutation_id)`
```

<img src="man/figures/README-unnamed-chunk-27-1.png" width="100%" />
