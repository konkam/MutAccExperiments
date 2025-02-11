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
#>  1         -9.49         -10.3         -10.7         -9.62         -10.2
#>  2         -9.51         -10.3         -10.7         -9.71         -10.3
#>  3         -9.56         -10.3         -10.7         -9.66         -10.3
#>  4         -9.52         -10.3         -10.5         -9.66         -10.2
#>  5         -9.53         -10.2         -10.6         -9.66         -10.1
#>  6         -9.56         -10.2         -10.8         -9.72         -10.2
#>  7         -9.52         -10.2         -10.6         -9.61         -10.5
#>  8         -9.55         -10.3         -10.8         -9.63         -10.4
#>  9         -9.65         -10.3         -10.9         -9.62         -10.5
#> 10         -9.63         -10.3         -10.4         -9.69         -10.2
#> # ℹ 7,990 more rows
#> # ℹ 3 more variables: `log10_mu[6]` <dbl>, iteration <int>, chain_id <int>
```

``` r
extract_posterior_samples(fit_GCM_model, type = "hyperparameters")
#> # A tibble: 8,000 × 5
#>    log10_sigma log10_mean loglikelihood iteration chain_id
#>          <dbl>      <dbl>         <dbl>     <int>    <int>
#>  1       0.551     -10.3          -18.5         1        1
#>  2       0.354     -10.2          -18.7         2        1
#>  3       0.544     -10.2          -17.5         3        1
#>  4       0.488     -10.0          -18.8         4        1
#>  5       0.401     -10.1          -20.2         5        1
#>  6       0.507      -9.91         -19.7         6        1
#>  7       0.392     -10.0          -18.8         7        1
#>  8       0.409     -10.2          -17.4         8        1
#>  9       0.470     -10.1          -22.0         9        1
#> 10       0.679      -9.83         -22.4        10        1
#> # ℹ 7,990 more rows
```

``` r
extract_posterior_samples(fit_GCM_model, type = "predictive")
#> # A tibble: 8,000 × 8
#>    `m_pred[1]` `m_pred[2]` `m_pred[3]` `m_pred[4]` `m_pred[5]` `m_pred[6]`
#>          <dbl>       <dbl>       <dbl>       <dbl>       <dbl>       <dbl>
#>  1         140          32           6         144          35          32
#>  2         124          27           4          86          25          22
#>  3         106          23          12         133          20          40
#>  4         121          25          12         109          30          29
#>  5         126          31           9         108          41          33
#>  6         116          20           8         109          28          37
#>  7         131          22           7         111          20          26
#>  8         118          17          10         135          21          16
#>  9          91          20           3         146          25          24
#> 10          97          27          15         102          34          41
#> # ℹ 7,990 more rows
#> # ℹ 2 more variables: iteration <int>, chain_id <int>
```

``` r
extract_posterior_samples(fit_GCM_model, type = "prior")
#> # A tibble: 8,000 × 5
#>    log10_mean_prior log10_sigma_prior log10_mu_prior iteration chain_id
#>               <dbl>             <dbl>          <dbl>     <int>    <int>
#>  1           -10.1             0.228          -10.3          1        1
#>  2            -7.37            0.0478          -7.18         2        1
#>  3            -6.28            1.37            -7.02         3        1
#>  4            -9.88            0.783           -9.53         4        1
#>  5            -8.30            0.171           -8.70         5        1
#>  6           -11.5             1.56           -12.3          6        1
#>  7           -14.3             0.0272         -13.6          7        1
#>  8            -8.57            0.669           -8.78         8        1
#>  9            -8.01            0.0512          -8.55         9        1
#> 10            -9.78            0.328           -9.88        10        1
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
#>  1         1         140          32           6         144          35
#>  2         2         124          27           4          86          25
#>  3         3         106          23          12         133          20
#>  4         4         121          25          12         109          30
#>  5         5         126          31           9         108          41
#>  6         6         116          20           8         109          28
#>  7         7         131          22           7         111          20
#>  8         8         118          17          10         135          21
#>  9         9          91          20           3         146          25
#> 10        10          97          27          15         102          34
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
#> m_pred[1]          8.500000e+01 116.0000000 145.0000000 116.8507500 15.29772784
#> m_pred[2]          9.000000e+00  23.0000000  35.0000000  23.2358750  6.67111692
#> m_pred[3]          1.000000e+00   9.0000000  17.0000000   9.1831250  4.25333417
#> m_pred[4]          8.700000e+01 118.0000000 147.0000000 118.9542500 15.44015728
#> m_pred[5]          1.300000e+01  25.0000000  39.0000000  25.4637500  6.96322459
#> m_pred[6]          1.100000e+01  25.0000000  38.0000000  25.3215000  7.07128971
#> log10_mu[1]       -9.631653e+00  -9.5538081  -9.4742681  -9.5544327  0.04061876
#> log10_mu[2]       -1.044171e+01 -10.2596422 -10.0928162 -10.2636564  0.08975740
#> log10_mu[3]       -1.096679e+01 -10.6698842 -10.4039231 -10.6794690  0.14540333
#> log10_mu[4]       -9.737088e+00  -9.6552554  -9.5802550  -9.6559433  0.03999487
#> log10_mu[5]       -1.049184e+01 -10.3268616 -10.1658659 -10.3299029  0.08442310
#> log10_mu[6]       -1.050677e+01 -10.3313432 -10.1697478 -10.3346788  0.08641319
#> log10_mean        -1.058157e+01 -10.1245873  -9.6845702 -10.1223498  0.22440396
#> log10_sigma        2.369463e-01   0.4838387   0.8626661   0.5155386  0.17598247
#> log10_mean_prior  -1.396346e+01  -8.0197667  -2.3445246  -8.0071474  2.99217560
#> log10_sigma_prior  8.634133e-05   0.4158693   1.4519655   0.5413715  0.49911097
#> log10_mu_prior    -1.370850e+01  -8.0419005  -1.8884892  -8.0174145  3.02930521
#> loglikelihood     -2.255556e+01 -18.7929203 -16.4778883 -19.1358460  1.73339548
#>                   Mode        MCerr MC%ofSD SSeff         AC.20      psrf
#> m_pred[1]          112 0.1680013522     1.1  8291  0.0070851022 1.0000327
#> m_pred[2]           23 0.0785262559     1.2  7217  0.0038756070 1.0009275
#> m_pred[3]            7 0.0532463729     1.3  6381  0.0087538439 0.9998800
#> m_pred[4]          118 0.1815593106     1.2  7232 -0.0087333234 0.9998832
#> m_pred[5]           25 0.0787992000     1.1  7809  0.0054157602 1.0003916
#> m_pred[6]           25 0.0831556921     1.2  7231 -0.0092646689 0.9999902
#> log10_mu[1]         NA 0.0004731776     1.2  7369 -0.0031226355 0.9999676
#> log10_mu[2]         NA 0.0010782207     1.2  6930  0.0040584368 1.0011683
#> log10_mu[3]         NA 0.0018440814     1.3  6217  0.0154728179 1.0000783
#> log10_mu[4]         NA 0.0004875253     1.2  6730 -0.0003053877 1.0002503
#> log10_mu[5]         NA 0.0010037881     1.2  7074 -0.0114055701 1.0006017
#> log10_mu[6]         NA 0.0010454930     1.2  6832 -0.0110968640 0.9999722
#> log10_mean          NA 0.0025089126     1.1  8000  0.0254882617 1.0004358
#> log10_sigma         NA 0.0026553236     1.5  4392 -0.0142779188 1.0000693
#> log10_mean_prior    NA 0.0334323444     1.1  8010 -0.0023991193 1.0003183
#> log10_sigma_prior   NA 0.0056433769     1.1  7822  0.0070936329 1.0001844
#> log10_mu_prior      NA 0.0342827786     1.1  7808 -0.0042830191 1.0004383
#> loglikelihood       NA 0.0223584205     1.3  6011  0.0127706561 0.9999209
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
#> 1           1      117.             116           88         149 
#> 2           2       23.2             23           11          37 
#> 3           3        9.18             9            2          19 
#> 4           4      119.             118           90         150.
#> 5           5       25.5             25           13          40 
#> 6           6       25.3             25           13          40 
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
#> The model will need to be run for a further 2415 updates.  This will
#> take approximately 0.3 seconds.
#> 
#> Calling the simulation using the rjags method...
#> Note: the model did not require adaptation
#> Running the model for 2415 iterations...
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
#>  1   2.93e-10   8.68e-11   3.68e-11   2.02e-10   4.77e-11   6.13e-11
#>  2   2.86e-10   5.40e-11   1.99e-11   1.63e-10   4.84e-11   5.99e-11
#>  3   3.06e-10   6.55e-11   2.00e-11   2.26e-10   4.46e-11   3.77e-11
#>  4   2.54e-10   7.32e-11   3.16e-11   1.91e-10   4.61e-11   6.05e-11
#>  5   2.45e-10   6.34e-11   3.61e-11   2.33e-10   4.59e-11   5.17e-11
#>  6   2.67e-10   5.18e-11   2.98e-11   2.30e-10   4.29e-11   4.07e-11
#>  7   2.45e-10   6.73e-11   1.83e-11   2.26e-10   4.39e-11   4.39e-11
#>  8   2.85e-10   5.77e-11   1.54e-11   1.88e-10   3.64e-11   6.46e-11
#>  9   3.01e-10   5.75e-11   1.32e-11   2.47e-10   4.85e-11   6.56e-11
#> 10   2.82e-10   6.62e-11   2.95e-11   2.23e-10   4.60e-11   6.16e-11
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
#>  1          0.332                1.22                      0.839        -10.0 
#>  2          0.613                0.757                     0.862        -10.0 
#>  3          0.337                1.46                      0.815        -10.0 
#>  4          0.578                0.621                     0.791        -10.2 
#>  5          0.479                0.702                     0.771         -9.76
#>  6          0.730                1.14                      0.548        -10.2 
#>  7          0.478                1.32                      0.562        -10.4 
#>  8          0.416                1.12                      0.581        -10.0 
#>  9          0.487                1.27                      0.747        -10.1 
#> 10          0.662                0.896                     0.634         -9.65
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
#>  1            129             26             15            101             18
#>  2            124             21             14             82             15
#>  3            131             27             16            115             15
#>  4             95             31              8             87             25
#>  5             83             31             21            121             23
#>  6            114             20             11            104             30
#>  7            116             37              8            110             29
#>  8            109             40              9            115             19
#>  9            123             29              2            134             34
#> 10            110             15             16            110             25
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
#>  1         1   2.93e-10   8.68e-11   3.68e-11   2.02e-10   4.77e-11   6.13e-11
#>  2         2   2.86e-10   5.40e-11   1.99e-11   1.63e-10   4.84e-11   5.99e-11
#>  3         3   3.06e-10   6.55e-11   2.00e-11   2.26e-10   4.46e-11   3.77e-11
#>  4         4   2.54e-10   7.32e-11   3.16e-11   1.91e-10   4.61e-11   6.05e-11
#>  5         5   2.45e-10   6.34e-11   3.61e-11   2.33e-10   4.59e-11   5.17e-11
#>  6         6   2.67e-10   5.18e-11   2.98e-11   2.30e-10   4.29e-11   4.07e-11
#>  7         7   2.45e-10   6.73e-11   1.83e-11   2.26e-10   4.39e-11   4.39e-11
#>  8         8   2.85e-10   5.77e-11   1.54e-11   1.88e-10   3.64e-11   6.46e-11
#>  9         9   3.01e-10   5.75e-11   1.32e-11   2.47e-10   4.85e-11   6.56e-11
#> 10        10   2.82e-10   6.62e-11   2.95e-11   2.23e-10   4.60e-11   6.16e-11
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
#> mu_wt[1]                   2.327951e-10  2.802698e-10  3.311966e-10
#> mu_wt[2]                   3.878450e-11  6.023438e-11  8.454468e-11
#> mu_wt[3]                   8.553859e-12  2.078016e-11  3.505977e-11
#> mu_wt[4]                   1.834168e-10  2.203658e-10  2.620052e-10
#> mu_wt[5]                   3.030778e-11  4.644774e-11  6.447045e-11
#> mu_wt[6]                   3.271482e-11  5.072579e-11  6.975354e-11
#> mu_proofminus[1]           5.738331e-08  6.461374e-08  7.322929e-08
#> mu_proofminus[2]           2.967323e-09  5.017188e-09  7.378170e-09
#> mu_proofminus[3]           1.622340e-10  7.300670e-10  1.505963e-09
#> mu_proofminus[4]           3.150991e-08  3.626280e-08  4.110381e-08
#> mu_proofminus[5]           1.997442e-09  3.316200e-09  4.673927e-09
#> mu_proofminus[6]           1.154593e-09  2.222003e-09  3.408469e-09
#> mu_MMRminus[1]             3.683481e-08  3.830242e-08  3.980518e-08
#> mu_MMRminus[2]             9.581847e-11  1.878703e-10  2.883081e-10
#> mu_MMRminus[3]             2.109516e-10  3.365824e-10  4.952567e-10
#> mu_MMRminus[4]             2.712172e-08  2.824642e-08  2.940736e-08
#> mu_MMRminus[5]             5.407363e-10  7.151529e-10  9.028115e-10
#> mu_MMRminus[6]             3.285426e-10  4.673548e-10  6.212174e-10
#> mu_MMRproofminus[1]        6.524792e-07  7.284840e-07  8.088735e-07
#> mu_MMRproofminus[2]        6.742319e-09  1.319025e-08  2.131176e-08
#> mu_MMRproofminus[3]        1.223218e-09  5.227506e-09  1.086695e-08
#> mu_MMRproofminus[4]        3.582062e-07  4.069891e-07  4.549029e-07
#> mu_MMRproofminus[5]        1.485456e-08  2.323592e-08  3.373737e-08
#> mu_MMRproofminus[6]        6.079979e-09  1.218727e-08  1.872308e-08
#> log10_mean_wt             -1.056377e+01 -1.011476e+01 -9.655947e+00
#> log10_mean_MMRminus       -9.613273e+00 -8.754237e+00 -7.933441e+00
#> log10_mean_MMRproofminus  -8.147945e+00 -7.378490e+00 -6.631994e+00
#> log10_sigma_wt             2.451438e-01  4.805677e-01  8.888567e-01
#> log10_sigma_MMRminus       5.442971e-01  9.643834e-01  1.619155e+00
#> log10_sigma_MMRproofminus  5.069622e-01  8.657077e-01  1.473773e+00
#> theta4                     7.020775e-02  8.186484e-02  9.440094e-02
#> log10_mean_prior          -1.393770e+01 -7.977687e+00 -2.333797e+00
#> log10_sigma_prior          8.707254e-05  4.151186e-01  1.477721e+00
#> log10_mu_prior            -1.400765e+01 -7.980584e+00 -2.102963e+00
#> theta4_prior               2.538266e-02  4.978741e-01  9.728945e-01
#> m_wt_pred[1]               8.500000e+01  1.170000e+02  1.440000e+02
#> m_wt_pred[2]               1.100000e+01  2.500000e+01  3.800000e+01
#> m_wt_pred[3]               2.000000e+00  9.000000e+00  1.700000e+01
#> m_wt_pred[4]               8.700000e+01  1.180000e+02  1.460000e+02
#> m_wt_pred[5]               1.200000e+01  2.500000e+01  3.800000e+01
#> m_wt_pred[6]               1.300000e+01  2.700000e+01  4.100000e+01
#> m_proofminus_pred[1]       1.660000e+02  2.030000e+02  2.400000e+02
#> m_proofminus_pred[2]       6.000000e+00  1.600000e+01  2.600000e+01
#> m_proofminus_pred[3]       0.000000e+00  2.000000e+00  6.000000e+00
#> m_proofminus_pred[4]       1.150000e+02  1.460000e+02  1.760000e+02
#> m_proofminus_pred[5]       5.000000e+00  1.300000e+01  2.200000e+01
#> m_proofminus_pred[6]       2.000000e+00  9.000000e+00  1.600000e+01
#> m_MMRminus_pred[1]         2.283000e+03  2.417000e+03  2.552000e+03
#> m_MMRminus_pred[2]         3.000000e+00  1.200000e+01  2.100000e+01
#> m_MMRminus_pred[3]         8.000000e+00  2.100000e+01  3.300000e+01
#> m_MMRminus_pred[4]         2.155000e+03  2.290000e+03  2.417000e+03
#> m_MMRminus_pred[5]         3.800000e+01  5.800000e+01  7.900000e+01
#> m_MMRminus_pred[6]         2.100000e+01  3.800000e+01  5.400000e+01
#> m_MMRproofminus_pred[1]    2.340000e+02  2.790000e+02  3.220000e+02
#> m_MMRproofminus_pred[2]    0.000000e+00  5.000000e+00  1.000000e+01
#> m_MMRproofminus_pred[3]    0.000000e+00  2.000000e+00  6.000000e+00
#> m_MMRproofminus_pred[4]    1.630000e+02  2.000000e+02  2.350000e+02
#> m_MMRproofminus_pred[5]    3.000000e+00  1.100000e+01  1.900000e+01
#> m_MMRproofminus_pred[6]    1.000000e+00  6.000000e+00  1.200000e+01
#> loglikelihood             -1.111924e+01 -1.001953e+01 -9.140585e+00
#>                                    Mean           SD         Mode        MCerr
#> mu_wt[1]                   2.811662e-10 2.559439e-11 1.963550e-10          Inf
#> mu_wt[2]                   6.102106e-11 1.182950e-11 2.404215e-11          Inf
#> mu_wt[3]                   2.163441e-11 6.973825e-12 4.390707e-12          Inf
#> mu_wt[4]                   2.207529e-10 2.020440e-11 1.563142e-10          Inf
#> mu_wt[5]                   4.703840e-11 8.796767e-12 2.034990e-11          Inf
#> mu_wt[6]                   5.117721e-11 9.558054e-12 2.393531e-11          Inf
#> mu_proofminus[1]           6.471076e-08 4.014246e-09           NA          Inf
#> mu_proofminus[2]           5.100460e-09 1.148052e-09 2.033352e-09          Inf
#> mu_proofminus[3]           7.899690e-10 3.734857e-10 4.361637e-11          Inf
#> mu_proofminus[4]           3.631126e-08 2.452277e-09           NA          Inf
#> mu_proofminus[5]           3.359837e-09 6.836617e-10 1.315701e-09          Inf
#> mu_proofminus[6]           2.278051e-09 5.972362e-10 6.717950e-10          Inf
#> mu_MMRminus[1]             3.830499e-08 7.654536e-10           NA          Inf
#> mu_MMRminus[2]             1.928288e-10 5.070796e-11 5.773543e-11          Inf
#> mu_MMRminus[3]             3.428882e-10 7.368787e-11 1.300001e-10          Inf
#> mu_MMRminus[4]             2.825227e-08 5.839610e-10           NA          Inf
#> mu_MMRminus[5]             7.200479e-10 9.313302e-11 4.184361e-10          Inf
#> mu_MMRminus[6]             4.712380e-10 7.529038e-11 2.466588e-10          Inf
#> mu_MMRproofminus[1]        7.300237e-07 3.984512e-08           NA 5.431307e-10
#> mu_MMRproofminus[2]        1.360602e-08 3.860490e-09           NA          Inf
#> mu_MMRproofminus[3]        5.617113e-09 2.655363e-09           NA          Inf
#> mu_MMRproofminus[4]        4.076550e-07 2.502621e-08           NA 3.293856e-10
#> mu_MMRproofminus[5]        2.362706e-08 4.873306e-09           NA          Inf
#> mu_MMRproofminus[6]        1.248955e-08 3.323500e-09           NA          Inf
#> log10_mean_wt             -1.011683e+01 2.268273e-01           NA 2.536007e-03
#> log10_mean_MMRminus       -8.750278e+00 4.217767e-01           NA 4.773414e-03
#> log10_mean_MMRproofminus  -7.382480e+00 3.812310e-01           NA 4.195541e-03
#> log10_sigma_wt             5.175604e-01 1.851305e-01           NA 2.563279e-03
#> log10_sigma_MMRminus       1.020221e+00 3.036281e-01           NA 4.017757e-03
#> log10_sigma_MMRproofminus  9.121044e-01 2.685323e-01           NA 3.452433e-03
#> theta4                     8.208294e-02 6.318175e-03           NA 7.984890e-05
#> log10_mean_prior          -8.008181e+00 2.993814e+00           NA 3.256495e-02
#> log10_sigma_prior          5.424683e-01 5.145400e-01           NA 5.752732e-03
#> log10_mu_prior            -8.011770e+00 3.083522e+00           NA 3.350059e-02
#> theta4_prior               4.989018e-01 2.872846e-01           NA 3.211939e-03
#> m_wt_pred[1]               1.171835e+02 1.513076e+01 1.190000e+02 1.710729e-01
#> m_wt_pred[2]               2.541525e+01 6.979857e+00 2.200000e+01 8.228545e-02
#> m_wt_pred[3]               9.060375e+00 4.216480e+00 8.000000e+00 4.715195e-02
#> m_wt_pred[4]               1.182043e+02 1.536275e+01 1.190000e+02 1.760338e-01
#> m_wt_pred[5]               2.519375e+01 6.887354e+00 2.300000e+01 7.388653e-02
#> m_wt_pred[6]               2.742650e+01 7.370338e+00 2.800000e+01 8.365249e-02
#> m_proofminus_pred[1]       2.038060e+02 1.920542e+01 1.960000e+02 2.143858e-01
#> m_proofminus_pred[2]       1.605588e+01 5.382019e+00 1.500000e+01 6.017280e-02
#> m_proofminus_pred[3]       2.461250e+00 1.919885e+00 1.000000e+00 2.146497e-02
#> m_proofminus_pred[4]       1.467281e+02 1.567650e+01 1.400000e+02 1.752686e-01
#> m_proofminus_pred[5]       1.363875e+01 4.641033e+00 1.200000e+01 5.099714e-02
#> m_proofminus_pred[6]       9.197125e+00 3.866411e+00 9.000000e+00 4.322779e-02
#> m_MMRminus_pred[1]         2.417946e+03 6.902582e+01 2.443000e+03 7.717322e-01
#> m_MMRminus_pred[2]         1.212262e+01 4.824469e+00 1.100000e+01 5.798598e-02
#> m_MMRminus_pred[3]         2.167250e+01 6.554320e+00 1.900000e+01 7.278639e-02
#> m_MMRminus_pred[4]         2.290294e+03 6.748818e+01 2.281000e+03 7.545407e-01
#> m_MMRminus_pred[5]         5.829600e+01 1.076361e+01 5.900000e+01 1.218167e-01
#> m_MMRminus_pred[6]         3.812150e+01 8.662195e+00 4.000000e+01 9.684628e-02
#> m_MMRproofminus_pred[1]    2.796675e+02 2.259138e+01 2.770000e+02 2.753690e-01
#> m_MMRproofminus_pred[2]    5.197625e+00 2.734731e+00 5.000000e+00 3.170613e-02
#> m_MMRproofminus_pred[3]    2.179250e+00 1.790536e+00 1.000000e+00 2.024978e-02
#> m_MMRproofminus_pred[4]    2.005220e+02 1.874260e+01 2.020000e+02 2.206713e-01
#> m_MMRproofminus_pred[5]    1.163588e+01 4.187658e+00 1.100000e+01 4.465622e-02
#> m_MMRproofminus_pred[6]    6.144625e+00 2.989787e+00 5.000000e+00 3.342683e-02
#> loglikelihood             -1.006969e+01 5.196376e-01           NA 6.248291e-03
#>                           MC%ofSD SSeff         AC.30      psrf
#> mu_wt[1]                      Inf     0 -0.0299206349 1.0009091
#> mu_wt[2]                      Inf     0 -0.0067555435 1.0000397
#> mu_wt[3]                      Inf     0  0.0102495069 1.0012426
#> mu_wt[4]                      Inf     0 -0.0089279049 0.9999538
#> mu_wt[5]                      Inf     0 -0.0033088434 0.9999276
#> mu_wt[6]                      Inf     0  0.0055710863 1.0006566
#> mu_proofminus[1]              Inf     0 -0.0046076357 1.0006555
#> mu_proofminus[2]              Inf     0  0.0301897031 0.9999085
#> mu_proofminus[3]              Inf     0 -0.0099551549 0.9998763
#> mu_proofminus[4]              Inf     0 -0.0170502140 1.0001638
#> mu_proofminus[5]              Inf     0 -0.0157022389 1.0002884
#> mu_proofminus[6]              Inf     0  0.0025930327 1.0000412
#> mu_MMRminus[1]                Inf     0 -0.0040191319 1.0000185
#> mu_MMRminus[2]                Inf     0  0.0122256516 0.9999147
#> mu_MMRminus[3]                Inf     0 -0.0072941686 1.0000672
#> mu_MMRminus[4]                Inf     0  0.0139515971 1.0002459
#> mu_MMRminus[5]                Inf     0 -0.0159355788 1.0003884
#> mu_MMRminus[6]                Inf     0  0.0082192825 0.9999351
#> mu_MMRproofminus[1]           1.4  5382 -0.0186183687 0.9999480
#> mu_MMRproofminus[2]           Inf     0 -0.0195780954 0.9999554
#> mu_MMRproofminus[3]           Inf     0 -0.0040938886 1.0000076
#> mu_MMRproofminus[4]           1.3  5773 -0.0057739919 0.9999355
#> mu_MMRproofminus[5]           Inf     0 -0.0149633728 1.0002123
#> mu_MMRproofminus[6]           Inf     0  0.0025052949 0.9999785
#> log10_mean_wt                 1.1  8000  0.0183731247 0.9999888
#> log10_mean_MMRminus           1.1  7807 -0.0005059905 1.0001584
#> log10_mean_MMRproofminus      1.1  8257  0.0038325116 1.0002916
#> log10_sigma_wt                1.4  5216 -0.0027790154 1.0001651
#> log10_sigma_MMRminus          1.3  5711  0.0137332211 1.0005588
#> log10_sigma_MMRproofminus     1.3  6050 -0.0017546800 1.0003655
#> theta4                        1.3  6261 -0.0161741802 1.0002607
#> log10_mean_prior              1.1  8452 -0.0150341807 1.0000549
#> log10_sigma_prior             1.1  8000  0.0056843592 0.9999635
#> log10_mu_prior                1.1  8472 -0.0157812433 0.9999814
#> theta4_prior                  1.1  8000 -0.0033222214 1.0003683
#> m_wt_pred[1]                  1.1  7823 -0.0299419749 0.9999361
#> m_wt_pred[2]                  1.2  7195 -0.0092572633 0.9998845
#> m_wt_pred[3]                  1.1  7997  0.0043540660 1.0012286
#> m_wt_pred[4]                  1.1  7616 -0.0201016401 1.0006514
#> m_wt_pred[5]                  1.1  8689 -0.0031991078 0.9999434
#> m_wt_pred[6]                  1.1  7763  0.0042534341 1.0009673
#> m_proofminus_pred[1]          1.1  8025 -0.0117740989 1.0011675
#> m_proofminus_pred[2]          1.1  8000  0.0072555415 0.9999745
#> m_proofminus_pred[3]          1.1  8000  0.0030874913 0.9999205
#> m_proofminus_pred[4]          1.1  8000 -0.0206284764 1.0008646
#> m_proofminus_pred[5]          1.1  8282 -0.0083700450 0.9999399
#> m_proofminus_pred[6]          1.1  8000  0.0017812940 0.9999223
#> m_MMRminus_pred[1]            1.1  8000  0.0056875314 1.0001707
#> m_MMRminus_pred[2]            1.2  6922  0.0066786750 1.0000139
#> m_MMRminus_pred[3]            1.1  8109  0.0005668369 1.0005649
#> m_MMRminus_pred[4]            1.1  8000  0.0245971647 1.0000094
#> m_MMRminus_pred[5]            1.1  7807 -0.0116870111 1.0008903
#> m_MMRminus_pred[6]            1.1  8000 -0.0074941831 1.0000579
#> m_MMRproofminus_pred[1]       1.2  6731  0.0108754189 0.9999381
#> m_MMRproofminus_pred[2]       1.2  7439 -0.0269064626 0.9999640
#> m_MMRproofminus_pred[3]       1.1  7819  0.0097472287 1.0002654
#> m_MMRproofminus_pred[4]       1.2  7214 -0.0155682847 1.0000283
#> m_MMRproofminus_pred[5]       1.1  8794 -0.0076853121 0.9999616
#> m_MMRproofminus_pred[6]       1.1  8000 -0.0009421045 0.9999737
#> loglikelihood                 1.2  6916  0.0059638308 1.0003360
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
#>  1      117.             117          89          148             106 
#>  2       25.4             25          13           40              21 
#>  3        9.06             9           2           19               6 
#>  4      118.             118          89          149             108.
#>  5       25.2             25          13           40              20 
#>  6       27.4             27          15           43              22 
#>  7      204.             203         168          242.            191 
#>  8       16.1             16           7           28              12 
#>  9        2.46             2           0            7               1 
#> 10      147.             146         118.         179             136 
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
