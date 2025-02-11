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

    git clone git@github.com:konkam/MutAccExperiments.git

or

    git clone https://github.com/konkam/MutAccExperiments.git

You should then install the present package on your computer, using a
command such a the following, from inside the package folder:

    Rscript -e "devtools::install_local('.')"

Alternatively, you may open the folder using the Rstudio editor and
press `Ctrl + Alt + B`.

- Dependencies should be:
  - R package devtools
  - JAGS + the R package RJags (On ubuntu, jags is on the canonical
    repository, available by `apt install jags`)
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
#>  1         -9.52         -10.2         -10.4         -9.64         -10.3
#>  2         -9.60         -10.3         -10.8         -9.69         -10.4
#>  3         -9.65         -10.3         -10.7         -9.59         -10.2
#>  4         -9.54         -10.2         -10.6         -9.65         -10.1
#>  5         -9.54         -10.3         -10.8         -9.66         -10.4
#>  6         -9.56         -10.2         -10.6         -9.65         -10.3
#>  7         -9.52         -10.4         -10.7         -9.69         -10.5
#>  8         -9.54         -10.1         -10.6         -9.68         -10.3
#>  9         -9.54         -10.3         -10.6         -9.71         -10.3
#> 10         -9.62         -10.4         -10.7         -9.69         -10.4
#> # ℹ 7,990 more rows
#> # ℹ 3 more variables: `log10_mu[6]` <dbl>, iteration <int>, chain_id <int>
```

``` r
extract_posterior_samples(fit_GCM_model, type = "hyperparameters")
#> # A tibble: 8,000 × 5
#>    log10_sigma log10_mean loglikelihood iteration chain_id
#>          <dbl>      <dbl>         <dbl>     <int>    <int>
#>  1       0.400     -10.0          -20.3         1        1
#>  2       0.376      -9.75         -18.6         2        1
#>  3       0.338     -10.2          -21.6         3        1
#>  4       0.285     -10.1          -19.3         4        1
#>  5       0.492     -10.3          -17.4         5        1
#>  6       0.480      -9.97         -18.0         6        1
#>  7       0.902     -10.3          -19.9         7        1
#>  8       0.761     -10.1          -19.2         8        1
#>  9       0.655      -9.91         -17.7         9        1
#> 10       0.760      -9.67         -19.9        10        1
#> # ℹ 7,990 more rows
```

``` r
extract_posterior_samples(fit_GCM_model, type = "predictive")
#> # A tibble: 8,000 × 8
#>    `m_pred[1]` `m_pred[2]` `m_pred[3]` `m_pred[4]` `m_pred[5]` `m_pred[6]`
#>          <dbl>       <dbl>       <dbl>       <dbl>       <dbl>       <dbl>
#>  1         118          16          22         135          32          30
#>  2         103          19           5         106          21          31
#>  3          84          20           6         125          28          21
#>  4         100          16           9         133          39          25
#>  5         115          19           4         121          16          29
#>  6         119          29          11         118          22          21
#>  7         138          12           6         120          19          15
#>  8         106          44          10         121          24          28
#>  9         109          20          12         106          30          26
#> 10          95          19           8         107          19          30
#> # ℹ 7,990 more rows
#> # ℹ 2 more variables: iteration <int>, chain_id <int>
```

``` r
extract_posterior_samples(fit_GCM_model, type = "prior")
#> # A tibble: 8,000 × 5
#>    log10_mean_prior log10_sigma_prior log10_mu_prior iteration chain_id
#>               <dbl>             <dbl>          <dbl>     <int>    <int>
#>  1           -10.5           0.000588         -10.6          1        1
#>  2            -3.94          0.803             -3.94         2        1
#>  3            -7.70          0.150             -7.53         3        1
#>  4           -11.6           0.380            -11.6          4        1
#>  5           -12.6           0.506            -11.8          5        1
#>  6           -11.2           2.63             -10.8          6        1
#>  7           -11.2           1.19             -10.8          7        1
#>  8           -14.1           0.154            -14.8          8        1
#>  9            -8.71          0.0511            -9.57         9        1
#> 10            -4.08          0.284             -3.77        10        1
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
#>  1         1         118          16          22         135          32
#>  2         2         103          19           5         106          21
#>  3         3          84          20           6         125          28
#>  4         4         100          16           9         133          39
#>  5         5         115          19           4         121          16
#>  6         6         119          29          11         118          22
#>  7         7         138          12           6         120          19
#>  8         8         106          44          10         121          24
#>  9         9         109          20          12         106          30
#> 10        10          95          19           8         107          19
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
#> m_pred[1]          8.800000e+01 116.0000000 147.0000000 116.5787500 15.28949867
#> m_pred[2]          1.000000e+01  23.0000000  36.0000000  23.3861250  6.75919694
#> m_pred[3]          2.000000e+00   9.0000000  17.0000000   9.1478750  4.28891680
#> m_pred[4]          8.900000e+01 118.0000000 148.0000000 118.9848750 15.21816662
#> m_pred[5]          1.200000e+01  25.0000000  39.0000000  25.4488750  7.05312729
#> m_pred[6]          1.200000e+01  25.0000000  38.0000000  25.3811250  6.98492787
#> log10_mu[1]       -9.632294e+00  -9.5539425  -9.4758132  -9.5550358  0.04019080
#> log10_mu[2]       -1.045075e+01 -10.2579999 -10.0973513 -10.2619441  0.08997468
#> log10_mu[3]       -1.096935e+01 -10.6722170 -10.4067876 -10.6800667  0.14571068
#> log10_mu[4]       -9.732272e+00  -9.6554109  -9.5775713  -9.6553494  0.03936527
#> log10_mu[5]       -1.049594e+01 -10.3271315 -10.1685507 -10.3310126  0.08452829
#> log10_mu[6]       -1.049927e+01 -10.3294252 -10.1678181 -10.3328598  0.08460154
#> log10_mean        -1.058880e+01 -10.1226905  -9.6810239 -10.1231112  0.22938509
#> log10_sigma        2.406210e-01   0.4796592   0.9039064   0.5218049  0.19252703
#> log10_mean_prior  -1.378902e+01  -8.0016290  -2.1642006  -8.0119910  2.96632955
#> log10_sigma_prior  1.059375e-04   0.4227572   1.5290972   0.5556562  0.51829465
#> log10_mu_prior    -1.375097e+01  -8.0055901  -2.0012493  -8.0135195  3.00630091
#> loglikelihood     -2.241103e+01 -18.7703165 -16.4554838 -19.0977112  1.71096327
#>                   Mode        MCerr MC%ofSD SSeff        AC.20      psrf
#> m_pred[1]          117 0.1762085143     1.2  7529  0.001944458 1.0001292
#> m_pred[2]           21 0.0732649031     1.1  8511 -0.004048939 1.0003128
#> m_pred[3]            8 0.0499708671     1.2  7367  0.004988326 1.0000809
#> m_pred[4]          117 0.1657291382     1.1  8432 -0.009479317 1.0003799
#> m_pred[5]           25 0.0804265269     1.1  7691 -0.005390145 1.0005324
#> m_pred[6]           23 0.0808836901     1.2  7458  0.019391190 1.0005947
#> log10_mu[1]         NA 0.0004841613     1.2  6891 -0.003165741 0.9999582
#> log10_mu[2]         NA 0.0010801284     1.2  6939  0.010159673 0.9999725
#> log10_mu[3]         NA 0.0017942551     1.2  6595  0.018369357 1.0010661
#> log10_mu[4]         NA 0.0004612202     1.2  7285 -0.012331715 0.9999286
#> log10_mu[5]         NA 0.0009972663     1.2  7184  0.002629354 1.0006181
#> log10_mu[6]         NA 0.0010052541     1.2  7083  0.021934281 0.9999742
#> log10_mean          NA 0.0026017867     1.1  7773 -0.002533986 1.0003134
#> log10_sigma         NA 0.0033794040     1.8  3246 -0.005283180 1.0005609
#> log10_mean_prior    NA 0.0331645725     1.1  8000 -0.007203961 1.0008385
#> log10_sigma_prior   NA 0.0058777534     1.1  7776  0.006558803 0.9999282
#> log10_mu_prior      NA 0.0336114659     1.1  8000 -0.004364755 1.0010271
#> loglikelihood       NA 0.0222607390     1.3  5907  0.005350306 1.0002975
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
#> 1           1      117.             116           88          148
#> 2           2       23.4             23           12           38
#> 3           3        9.15             9            2           19
#> 4           4      119.             118           91          151
#> 5           5       25.4             25           13           41
#> 6           6       25.4             25           13           40
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
MMR-deactivated and proofreading-deficient (MMR-proof-). The model
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
MMR-proof- for which both MMR and Proofreading are inactivated. The
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
\mu_{\text{MMR-proof-},i} &= \gamma_i \\
\mu_{\text{MMR-},i} &= \gamma_i q_{\text{proofreading},i} \\
\mu_{\text{proof-},i} &= \gamma_i \left(\theta + (1-\theta) q_{\text{MMR},i} \right) \\
\mu_{\text{wt},i} &= \gamma_i q_{\text{proofreading},i}q_{\text{MMR},i} \\
\end{align}
$$

We choose the following generative model:

$$
\begin{align}
m_{s, i} &\sim \text{Poisson}(n_{s, i} \mu_{s, i} t_{s, i}) \qquad \text{for } i = 1, \ldots, I, s = \text{wt}, \text{MMR-}, \text{proof-}, \text{MMR-proof-}\\
\log10 \mu_{s, i} &\sim \mathcal{N}(\nu_{s}, \sigma_{s})  \qquad \text{for } i = 1, \ldots, I, s = \text{wt}, \text{MMR-}, \text{MMR-proof-} \\
\nu_s &\sim \mathcal{N}(-8, 3) \qquad \text{for } s = \text{wt}, \text{MMR-}, \text{MMR-proof-}\\
\sigma_s &\sim \text{t}^+(0, 3, 5) \qquad \text{for } s = \text{wt}, \text{MMR-}, \text{MMR-proof-}\\
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
#> The model will need to be run for a further 2111 updates.  This will
#> take approximately 0.1 seconds.
#> 
#> Calling the simulation using the rjags method...
#> Note: the model did not require adaptation
#> Running the model for 2111 iterations...
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
#>  1   2.51e-10   6.06e-11   1.15e-11   1.86e-10   6.20e-11   6.43e-11
#>  2   3.41e-10   5.33e-11   2.57e-11   2.89e-10   4.25e-11   5.38e-11
#>  3   2.74e-10   6.07e-11   2.56e-11   2.36e-10   3.78e-11   5.41e-11
#>  4   2.96e-10   5.35e-11   2.37e-11   2.59e-10   4.07e-11   5.83e-11
#>  5   3.10e-10   8.47e-11   1.31e-11   2.01e-10   6.82e-11   5.45e-11
#>  6   2.80e-10   6.08e-11   2.13e-11   2.02e-10   3.16e-11   4.38e-11
#>  7   2.84e-10   7.97e-11   2.00e-11   2.14e-10   4.97e-11   5.00e-11
#>  8   2.77e-10   5.11e-11   3.06e-11   2.18e-10   3.61e-11   5.01e-11
#>  9   2.69e-10   5.99e-11   1.69e-11   2.36e-10   2.98e-11   4.68e-11
#> 10   2.88e-10   3.65e-11   8.56e-12   1.94e-10   3.66e-11   6.01e-11
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
#>  1          1.02                 0.974                     0.796        -10.2 
#>  2          0.478                1.13                      1.07          -9.95
#>  3          0.384                0.945                     0.870         -9.97
#>  4          0.686                0.674                     1.22         -10.5 
#>  5          0.428                0.741                     1.37          -9.87
#>  6          0.593                0.797                     1.21         -10.4 
#>  7          0.491                1.13                      0.859        -10.1 
#>  8          0.406                0.661                     1.04         -10.1 
#>  9          0.793                1.00                      0.442        -10.2 
#> 10          0.453                0.993                     0.744        -10.4 
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
#>  1            106             25              9             81             25
#>  2            141             26             16            153             21
#>  3            117             30             14            131             18
#>  4            128             32             14            133             18
#>  5            118             28             10            120             45
#>  6            116             25             12             99             18
#>  7            133             33             10            137             25
#>  8            112             18             21            110             23
#>  9             90             34              7            112             17
#> 10            126             21              4            113             13
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
#>  1         1   2.51e-10   6.06e-11   1.15e-11   1.86e-10   6.20e-11   6.43e-11
#>  2         2   3.41e-10   5.33e-11   2.57e-11   2.89e-10   4.25e-11   5.38e-11
#>  3         3   2.74e-10   6.07e-11   2.56e-11   2.36e-10   3.78e-11   5.41e-11
#>  4         4   2.96e-10   5.35e-11   2.37e-11   2.59e-10   4.07e-11   5.83e-11
#>  5         5   3.10e-10   8.47e-11   1.31e-11   2.01e-10   6.82e-11   5.45e-11
#>  6         6   2.80e-10   6.08e-11   2.13e-11   2.02e-10   3.16e-11   4.38e-11
#>  7         7   2.84e-10   7.97e-11   2.00e-11   2.14e-10   4.97e-11   5.00e-11
#>  8         8   2.77e-10   5.11e-11   3.06e-11   2.18e-10   3.61e-11   5.01e-11
#>  9         9   2.69e-10   5.99e-11   1.69e-11   2.36e-10   2.98e-11   4.68e-11
#> 10        10   2.88e-10   3.65e-11   8.56e-12   1.94e-10   3.66e-11   6.01e-11
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
#> mu_wt[1]                   2.327839e-10  2.798985e-10  3.341728e-10
#> mu_wt[2]                   3.892566e-11  6.014513e-11  8.442610e-11
#> mu_wt[3]                   9.634712e-12  2.085711e-11  3.585598e-11
#> mu_wt[4]                   1.810848e-10  2.203752e-10  2.592399e-10
#> mu_wt[5]                   3.085932e-11  4.645664e-11  6.612574e-11
#> mu_wt[6]                   3.343753e-11  5.082250e-11  7.007582e-11
#> mu_proofminus[1]           5.701508e-08  6.460470e-08  7.244388e-08
#> mu_proofminus[2]           3.063153e-09  5.029541e-09  7.428369e-09
#> mu_proofminus[3]           1.644918e-10  7.237817e-10  1.479493e-09
#> mu_proofminus[4]           3.183394e-08  3.623452e-08  4.135694e-08
#> mu_proofminus[5]           2.113965e-09  3.299763e-09  4.756957e-09
#> mu_proofminus[6]           1.211615e-09  2.210797e-09  3.458214e-09
#> mu_MMRminus[1]             3.678521e-08  3.829738e-08  3.981578e-08
#> mu_MMRminus[2]             1.010496e-10  1.886429e-10  2.973472e-10
#> mu_MMRminus[3]             1.986521e-10  3.367688e-10  4.848518e-10
#> mu_MMRminus[4]             2.713162e-08  2.826840e-08  2.940014e-08
#> mu_MMRminus[5]             5.401802e-10  7.183516e-10  9.051220e-10
#> mu_MMRminus[6]             3.317789e-10  4.677523e-10  6.238417e-10
#> mu_MMRproofminus[1]        6.572755e-07  7.292252e-07  8.098724e-07
#> mu_MMRproofminus[2]        6.874178e-09  1.327500e-08  2.143996e-08
#> mu_MMRproofminus[3]        1.133555e-09  5.117094e-09  1.052482e-08
#> mu_MMRproofminus[4]        3.574918e-07  4.066247e-07  4.547359e-07
#> mu_MMRproofminus[5]        1.457417e-08  2.318122e-08  3.312995e-08
#> mu_MMRproofminus[6]        6.605727e-09  1.216298e-08  1.914941e-08
#> log10_mean_wt             -1.057290e+01 -1.011072e+01 -9.665970e+00
#> log10_mean_MMRminus       -9.585402e+00 -8.760153e+00 -7.907844e+00
#> log10_mean_MMRproofminus  -8.129701e+00 -7.377963e+00 -6.605429e+00
#> log10_sigma_wt             2.449086e-01  4.777782e-01  8.700775e-01
#> log10_sigma_MMRminus       5.409720e-01  9.590875e-01  1.575047e+00
#> log10_sigma_MMRproofminus  4.849970e-01  8.637623e-01  1.448814e+00
#> theta4                     6.973326e-02  8.187160e-02  9.448267e-02
#> log10_mean_prior          -1.363009e+01 -7.982266e+00 -1.904545e+00
#> log10_sigma_prior          8.837835e-05  4.253404e-01  1.496424e+00
#> log10_mu_prior            -1.377523e+01 -7.980894e+00 -1.620060e+00
#> theta4_prior               1.598565e-02  5.034570e-01  9.616849e-01
#> m_wt_pred[1]               8.800000e+01  1.170000e+02  1.470000e+02
#> m_wt_pred[2]               1.100000e+01  2.500000e+01  3.800000e+01
#> m_wt_pred[3]               2.000000e+00  9.000000e+00  1.700000e+01
#> m_wt_pred[4]               8.800000e+01  1.180000e+02  1.470000e+02
#> m_wt_pred[5]               1.100000e+01  2.500000e+01  3.800000e+01
#> m_wt_pred[6]               1.400000e+01  2.700000e+01  4.100000e+01
#> m_proofminus_pred[1]       1.650000e+02  2.030000e+02  2.390000e+02
#> m_proofminus_pred[2]       6.000000e+00  1.600000e+01  2.600000e+01
#> m_proofminus_pred[3]       0.000000e+00  2.000000e+00  6.000000e+00
#> m_proofminus_pred[4]       1.140000e+02  1.470000e+02  1.760000e+02
#> m_proofminus_pred[5]       5.000000e+00  1.300000e+01  2.200000e+01
#> m_proofminus_pred[6]       2.000000e+00  9.000000e+00  1.600000e+01
#> m_MMRminus_pred[1]         2.279000e+03  2.416000e+03  2.552000e+03
#> m_MMRminus_pred[2]         3.000000e+00  1.200000e+01  2.100000e+01
#> m_MMRminus_pred[3]         8.000000e+00  2.100000e+01  3.300000e+01
#> m_MMRminus_pred[4]         2.159000e+03  2.291000e+03  2.425000e+03
#> m_MMRminus_pred[5]         3.800000e+01  5.800000e+01  7.900000e+01
#> m_MMRminus_pred[6]         2.100000e+01  3.800000e+01  5.400000e+01
#> m_MMRproofminus_pred[1]    2.370000e+02  2.790000e+02  3.240000e+02
#> m_MMRproofminus_pred[2]    0.000000e+00  5.000000e+00  1.000000e+01
#> m_MMRproofminus_pred[3]    0.000000e+00  2.000000e+00  5.000000e+00
#> m_MMRproofminus_pred[4]    1.610000e+02  2.000000e+02  2.350000e+02
#> m_MMRproofminus_pred[5]    4.000000e+00  1.100000e+01  1.900000e+01
#> m_MMRproofminus_pred[6]    0.000000e+00  6.000000e+00  1.100000e+01
#> loglikelihood             -1.113607e+01 -1.001488e+01 -9.159670e+00
#>                                    Mean           SD         Mode        MCerr
#> mu_wt[1]                   2.811217e-10 2.574373e-11 1.997443e-10          Inf
#> mu_wt[2]                   6.102108e-11 1.176555e-11 2.475396e-11          Inf
#> mu_wt[3]                   2.153485e-11 6.974056e-12 3.372108e-12          Inf
#> mu_wt[4]                   2.208314e-10 2.020393e-11 1.477559e-10          Inf
#> mu_wt[5]                   4.713274e-11 9.113829e-12 1.854151e-11          Inf
#> mu_wt[6]                   5.134210e-11 9.492030e-12 2.392044e-11          Inf
#> mu_proofminus[1]           6.467672e-08 3.967187e-09           NA          Inf
#> mu_proofminus[2]           5.122546e-09 1.143418e-09 1.918807e-09          Inf
#> mu_proofminus[3]           7.797698e-10 3.627272e-10 7.144058e-11          Inf
#> mu_proofminus[4]           3.631105e-08 2.452031e-09           NA          Inf
#> mu_proofminus[5]           3.342654e-09 6.819923e-10 1.221175e-09          Inf
#> mu_proofminus[6]           2.269260e-09 5.895548e-10 7.211594e-10          Inf
#> mu_MMRminus[1]             3.830643e-08 7.815281e-10           NA          Inf
#> mu_MMRminus[2]             1.934225e-10 5.132651e-11 5.067978e-11          Inf
#> mu_MMRminus[3]             3.419416e-10 7.366526e-11 1.464420e-10          Inf
#> mu_MMRminus[4]             2.827725e-08 5.860791e-10           NA          Inf
#> mu_MMRminus[5]             7.228444e-10 9.389197e-11 4.451510e-10          Inf
#> mu_MMRminus[6]             4.726093e-10 7.586726e-11 2.498110e-10          Inf
#> mu_MMRproofminus[1]        7.294306e-07 3.905606e-08           NA 5.183197e-10
#> mu_MMRproofminus[2]        1.369284e-08 3.891838e-09           NA          Inf
#> mu_MMRproofminus[3]        5.548305e-09 2.600412e-09           NA          Inf
#> mu_MMRproofminus[4]        4.075423e-07 2.492300e-08           NA 3.290904e-10
#> mu_MMRproofminus[5]        2.352464e-08 4.850475e-09           NA          Inf
#> mu_MMRproofminus[6]        1.244061e-08 3.307501e-09           NA          Inf
#> log10_mean_wt             -1.011018e+01 2.233758e-01           NA 2.497418e-03
#> log10_mean_MMRminus       -8.756444e+00 4.239059e-01           NA 4.739412e-03
#> log10_mean_MMRproofminus  -7.382478e+00 3.872854e-01           NA 4.280934e-03
#> log10_sigma_wt             5.143702e-01 1.797435e-01           NA 2.441226e-03
#> log10_sigma_MMRminus       1.013248e+00 2.956250e-01           NA 4.076501e-03
#> log10_sigma_MMRproofminus  9.122984e-01 2.673245e-01           NA 3.494114e-03
#> theta4                     8.210845e-02 6.314115e-03           NA 9.127651e-05
#> log10_mean_prior          -7.971547e+00 3.011641e+00           NA 3.310299e-02
#> log10_sigma_prior          5.499588e-01 5.025581e-01           NA 5.551708e-03
#> log10_mu_prior            -7.987467e+00 3.109614e+00           NA 3.476654e-02
#> theta4_prior               5.012467e-01 2.879129e-01           NA 3.223324e-03
#> m_wt_pred[1]               1.173259e+02 1.526558e+01 1.130000e+02 1.706743e-01
#> m_wt_pred[2]               2.546088e+01 7.061530e+00 2.300000e+01 8.146265e-02
#> m_wt_pred[3]               8.977375e+00 4.148194e+00 7.000000e+00 4.690009e-02
#> m_wt_pred[4]               1.182097e+02 1.540610e+01 1.190000e+02 1.767465e-01
#> m_wt_pred[5]               2.526663e+01 7.061499e+00 2.500000e+01 7.782326e-02
#> m_wt_pred[6]               2.748438e+01 7.215886e+00 2.600000e+01 8.090840e-02
#> m_proofminus_pred[1]       2.036267e+02 1.880738e+01 2.040000e+02 2.126367e-01
#> m_proofminus_pred[2]       1.615800e+01 5.456511e+00 1.500000e+01 5.919989e-02
#> m_proofminus_pred[3]       2.455250e+00 1.912643e+00 1.000000e+00 2.144559e-02
#> m_proofminus_pred[4]       1.467046e+02 1.587335e+01 1.480000e+02 1.836681e-01
#> m_proofminus_pred[5]       1.350150e+01 4.636802e+00 1.200000e+01 5.184103e-02
#> m_proofminus_pred[6]       9.171125e+00 3.865030e+00 8.000000e+00 4.373695e-02
#> m_MMRminus_pred[1]         2.417119e+03 7.003943e+01 2.394000e+03 8.036214e-01
#> m_MMRminus_pred[2]         1.224687e+01 4.775409e+00 1.100000e+01 5.680905e-02
#> m_MMRminus_pred[3]         2.157587e+01 6.578727e+00 2.100000e+01 7.284547e-02
#> m_MMRminus_pred[4]         2.292427e+03 6.736446e+01 2.287000e+03 7.649078e-01
#> m_MMRminus_pred[5]         5.861800e+01 1.083418e+01 5.500000e+01 1.211298e-01
#> m_MMRminus_pred[6]         3.829212e+01 8.666656e+00 3.600000e+01 9.689616e-02
#> m_MMRproofminus_pred[1]    2.791751e+02 2.226015e+01 2.770000e+02 2.671946e-01
#> m_MMRproofminus_pred[2]    5.246875e+00 2.756515e+00 5.000000e+00 3.355613e-02
#> m_MMRproofminus_pred[3]    2.136500e+00 1.774545e+00 1.000000e+00 1.984002e-02
#> m_MMRproofminus_pred[4]    2.005677e+02 1.891208e+01 2.000000e+02 2.269337e-01
#> m_MMRproofminus_pred[5]    1.150950e+01 4.118165e+00 1.000000e+01 4.624103e-02
#> m_MMRproofminus_pred[6]    6.133125e+00 2.944314e+00 5.000000e+00 3.291843e-02
#> loglikelihood             -1.007130e+01 5.164096e-01           NA 5.984995e-03
#>                           MC%ofSD SSeff         AC.30      psrf
#> mu_wt[1]                      Inf     0  0.0002666779 0.9998872
#> mu_wt[2]                      Inf     0 -0.0062990835 1.0004683
#> mu_wt[3]                      Inf     0 -0.0036710141 0.9998863
#> mu_wt[4]                      Inf     0 -0.0088135536 1.0000105
#> mu_wt[5]                      Inf     0 -0.0002916232 0.9999735
#> mu_wt[6]                      Inf     0  0.0041456800 0.9999593
#> mu_proofminus[1]              Inf     0  0.0012840763 1.0000249
#> mu_proofminus[2]              Inf     0 -0.0199248923 0.9999506
#> mu_proofminus[3]              Inf     0  0.0083741161 1.0007460
#> mu_proofminus[4]              Inf     0 -0.0089782740 1.0000823
#> mu_proofminus[5]              Inf     0  0.0039062567 0.9999277
#> mu_proofminus[6]              Inf     0 -0.0090499987 0.9999658
#> mu_MMRminus[1]                Inf     0 -0.0034945463 1.0000320
#> mu_MMRminus[2]                Inf     0  0.0013181808 1.0006709
#> mu_MMRminus[3]                Inf     0 -0.0003107580 1.0001585
#> mu_MMRminus[4]                Inf     0  0.0102051207 0.9998908
#> mu_MMRminus[5]                Inf     0 -0.0145885696 1.0005639
#> mu_MMRminus[6]                Inf     0 -0.0093704349 1.0009287
#> mu_MMRproofminus[1]           1.3  5678  0.0084176836 1.0001481
#> mu_MMRproofminus[2]           Inf     0 -0.0006416608 1.0010367
#> mu_MMRproofminus[3]           Inf     0 -0.0015565942 1.0007018
#> mu_MMRproofminus[4]           1.3  5735 -0.0024392970 0.9999334
#> mu_MMRproofminus[5]           Inf     0  0.0132802603 0.9999539
#> mu_MMRproofminus[6]           Inf     0 -0.0106195547 1.0000266
#> log10_mean_wt                 1.1  8000 -0.0016557657 1.0000640
#> log10_mean_MMRminus           1.1  8000  0.0122104234 1.0003004
#> log10_mean_MMRproofminus      1.1  8184  0.0072226424 1.0000332
#> log10_sigma_wt                1.4  5421 -0.0017874287 0.9999114
#> log10_sigma_MMRminus          1.4  5259 -0.0080359727 1.0006801
#> log10_sigma_MMRproofminus     1.3  5853  0.0026706622 0.9999200
#> theta4                        1.4  4785  0.0048460046 0.9999427
#> log10_mean_prior              1.1  8277 -0.0109989996 0.9998825
#> log10_sigma_prior             1.1  8194  0.0043285419 1.0000194
#> log10_mu_prior                1.1  8000 -0.0082849098 0.9999043
#> theta4_prior                  1.1  7978 -0.0097648100 1.0000670
#> m_wt_pred[1]                  1.1  8000  0.0004906000 0.9999935
#> m_wt_pred[2]                  1.2  7514 -0.0104719685 1.0007704
#> m_wt_pred[3]                  1.1  7823 -0.0009212799 1.0002639
#> m_wt_pred[4]                  1.1  7598  0.0013301167 1.0001187
#> m_wt_pred[5]                  1.1  8233 -0.0069745648 0.9998754
#> m_wt_pred[6]                  1.1  7954 -0.0010703466 1.0001595
#> m_proofminus_pred[1]          1.1  7823 -0.0117973431 1.0000351
#> m_proofminus_pred[2]          1.1  8495  0.0059753034 1.0001348
#> m_proofminus_pred[3]          1.1  7954  0.0049936821 1.0002101
#> m_proofminus_pred[4]          1.2  7469  0.0006156780 0.9998898
#> m_proofminus_pred[5]          1.1  8000  0.0039134175 0.9999348
#> m_proofminus_pred[6]          1.1  7809  0.0057142936 0.9998875
#> m_MMRminus_pred[1]            1.1  7596 -0.0006406876 0.9999412
#> m_MMRminus_pred[2]            1.2  7066  0.0119573913 1.0003053
#> m_MMRminus_pred[3]            1.1  8156  0.0064757903 0.9999198
#> m_MMRminus_pred[4]            1.1  7756  0.0166506333 0.9999801
#> m_MMRminus_pred[5]            1.1  8000 -0.0207418408 0.9999403
#> m_MMRminus_pred[6]            1.1  8000 -0.0066298886 1.0001356
#> m_MMRproofminus_pred[1]       1.2  6941 -0.0042902622 1.0003107
#> m_MMRproofminus_pred[2]       1.2  6748  0.0039150690 1.0000233
#> m_MMRproofminus_pred[3]       1.1  8000 -0.0009595961 1.0000889
#> m_MMRproofminus_pred[4]       1.2  6945 -0.0024534009 1.0002174
#> m_MMRproofminus_pred[5]       1.1  7931  0.0029576522 1.0001751
#> m_MMRproofminus_pred[6]       1.1  8000 -0.0096677466 1.0003121
#> loglikelihood                 1.2  7445  0.0000227049 1.0011400
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
#>  1      117.             117           89         149              107
#>  2       25.5             25           13          41               20
#>  3        8.98             9            2          18                6
#>  4      118.             118           90         149.             107
#>  5       25.3             25           13          40               20
#>  6       27.5             27           15          43               22
#>  7      204.             203          168         242              191
#>  8       16.2             16            7          28               12
#>  9        2.46             2            0           7                1
#> 10      147.             147          117         180              136
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
