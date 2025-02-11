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
#>  1         -9.57         -10.4         -10.7         -9.65         -10.3
#>  2         -9.57         -10.2         -10.6         -9.62         -10.3
#>  3         -9.47         -10.3         -10.6         -9.66         -10.2
#>  4         -9.55         -10.3         -10.7         -9.62         -10.2
#>  5         -9.60         -10.3         -10.6         -9.72         -10.2
#>  6         -9.56         -10.4         -10.9         -9.63         -10.2
#>  7         -9.56         -10.3         -10.6         -9.67         -10.4
#>  8         -9.56         -10.3         -10.7         -9.71         -10.4
#>  9         -9.58         -10.3         -10.9         -9.64         -10.3
#> 10         -9.53         -10.3         -10.7         -9.67         -10.5
#> # ℹ 7,990 more rows
#> # ℹ 3 more variables: `log10_mu[6]` <dbl>, iteration <int>, chain_id <int>
```

``` r
extract_posterior_samples(fit_GCM_model, type = "hyperparameters")
#> # A tibble: 8,000 × 5
#>    log10_sigma log10_mean loglikelihood iteration chain_id
#>          <dbl>      <dbl>         <dbl>     <int>    <int>
#>  1       0.380     -10.2          -17.4         1        1
#>  2       0.417     -10.4          -16.9         2        1
#>  3       0.389      -9.95         -20.6         3        1
#>  4       0.570     -10.1          -17.7         4        1
#>  5       0.570     -10.0          -20.6         5        1
#>  6       0.438     -10.2          -18.8         6        1
#>  7       0.440     -10.0          -16.8         7        1
#>  8       0.677     -10.1          -18.1         8        1
#>  9       0.686      -9.78         -17.1         9        1
#> 10       0.464     -10.0          -18.0        10        1
#> # ℹ 7,990 more rows
```

``` r
extract_posterior_samples(fit_GCM_model, type = "predictive")
#> # A tibble: 8,000 × 8
#>    `m_pred[1]` `m_pred[2]` `m_pred[3]` `m_pred[4]` `m_pred[5]` `m_pred[6]`
#>          <dbl>       <dbl>       <dbl>       <dbl>       <dbl>       <dbl>
#>  1         125          16           4         115          29          19
#>  2         126          24           9         144          29          38
#>  3         144          17           6          99          32          39
#>  4         119          25          13         112          29          33
#>  5          98          18           5         107          38          20
#>  6         120          21           6         144          32          20
#>  7         124          19          10         115          20          20
#>  8         121          21           8         108          17          26
#>  9         121          18           6         137          24          35
#> 10         117          23          12         110          19          20
#> # ℹ 7,990 more rows
#> # ℹ 2 more variables: iteration <int>, chain_id <int>
```

``` r
extract_posterior_samples(fit_GCM_model, type = "prior")
#> # A tibble: 8,000 × 5
#>    log10_mean_prior log10_sigma_prior log10_mu_prior iteration chain_id
#>               <dbl>             <dbl>          <dbl>     <int>    <int>
#>  1           -11.0             0.245          -11.1          1        1
#>  2           -11.4             1.38           -11.8          2        1
#>  3           -12.1             0.106          -12.1          3        1
#>  4            -6.41            0.855           -6.18         4        1
#>  5            -4.23            0.0393          -3.71         5        1
#>  6           -11.8             0.0371         -11.3          6        1
#>  7            -9.49            0.608           -9.33         7        1
#>  8            -4.77            0.690           -5.38         8        1
#>  9            -8.87            1.55            -8.59         9        1
#> 10           -15.1             0.894          -15.3         10        1
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
#>  1         1         125          16           4         115          29
#>  2         2         126          24           9         144          29
#>  3         3         144          17           6          99          32
#>  4         4         119          25          13         112          29
#>  5         5          98          18           5         107          38
#>  6         6         120          21           6         144          32
#>  7         7         124          19          10         115          20
#>  8         8         121          21           8         108          17
#>  9         9         121          18           6         137          24
#> 10        10         117          23          12         110          19
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
#> m_pred[1]          8.700000e+01 116.0000000 146.0000000 116.8916250 15.22047429
#> m_pred[2]          1.000000e+01  23.0000000  36.0000000  23.2863750  6.80846229
#> m_pred[3]          1.000000e+00   9.0000000  17.0000000   9.1955000  4.23178665
#> m_pred[4]          8.700000e+01 118.0000000 147.0000000 118.9963750 15.32342856
#> m_pred[5]          1.100000e+01  25.0000000  38.0000000  25.4275000  6.99727547
#> m_pred[6]          1.100000e+01  25.0000000  38.0000000  25.3717500  7.01433903
#> log10_mu[1]       -9.628216e+00  -9.5536892  -9.4694184  -9.5547107  0.04011571
#> log10_mu[2]       -1.042877e+01 -10.2570560 -10.0844221 -10.2611115  0.08898758
#> log10_mu[3]       -1.096065e+01 -10.6727184 -10.3979179 -10.6795646  0.14533935
#> log10_mu[4]       -9.732548e+00  -9.6547187  -9.5782344  -9.6553231  0.03959221
#> log10_mu[5]       -1.049503e+01 -10.3281620 -10.1602316 -10.3303734  0.08432389
#> log10_mu[6]       -1.049982e+01 -10.3292575 -10.1717682 -10.3313439  0.08399191
#> log10_mean        -1.056697e+01 -10.1215766  -9.6377657 -10.1204384  0.23165872
#> log10_sigma        2.483487e-01   0.4800730   0.8768879   0.5169073  0.18262777
#> log10_mean_prior  -1.392607e+01  -7.9417473  -2.2169174  -7.9544715  2.98762169
#> log10_sigma_prior  5.792443e-05   0.4209298   1.4752635   0.5491363  0.54390311
#> log10_mu_prior    -1.385244e+01  -7.9529628  -1.8467165  -7.9575286  3.04446014
#> loglikelihood     -2.243717e+01 -18.7531720 -16.4565098 -19.0818619  1.73333013
#>                   Mode        MCerr MC%ofSD SSeff        AC.20      psrf
#> m_pred[1]          113 0.1757826459     1.2  7497  0.007806499 0.9999011
#> m_pred[2]           21 0.0761209225     1.1  8000  0.006241835 0.9999664
#> m_pred[3]            8 0.0495724050     1.2  7287  0.009799433 0.9999091
#> m_pred[4]          115 0.1763384603     1.2  7551 -0.006756269 0.9999074
#> m_pred[5]           23 0.0805118304     1.2  7553  0.013324914 0.9999894
#> m_pred[6]           24 0.0804978429     1.1  7593  0.002764834 1.0007107
#> log10_mu[1]         NA 0.0004769692     1.2  7074 -0.005332135 0.9999020
#> log10_mu[2]         NA 0.0010506090     1.2  7174 -0.002306907 0.9999471
#> log10_mu[3]         NA 0.0018113556     1.2  6438 -0.011889332 1.0004182
#> log10_mu[4]         NA 0.0004721892     1.2  7031 -0.018740107 1.0000858
#> log10_mu[5]         NA 0.0010368465     1.2  6614 -0.004678526 1.0017140
#> log10_mu[6]         NA 0.0009984654     1.2  7076 -0.004543747 1.0011281
#> log10_mean          NA 0.0025900232     1.1  8000 -0.005447146 1.0014890
#> log10_sigma         NA 0.0029484372     1.6  3837 -0.009408964 1.0007030
#> log10_mean_prior    NA 0.0334026260     1.1  8000 -0.008227313 1.0007235
#> log10_sigma_prior   NA 0.0064076506     1.2  7205  0.004532133 1.0000246
#> log10_mu_prior      NA 0.0340380991     1.1  8000 -0.005201422 1.0007120
#> loglikelihood       NA 0.0227806141     1.3  5789  0.007169189 1.0001586
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
#> 2           2       23.3             23           11           38
#> 3           3        9.20             9            2           19
#> 4           4      119.             118           90          150
#> 5           5       25.4             25           13           40
#> 6           6       25.4             25           13           41
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
m_{s, i} &\sim \text{Poisson}(n_{s, i} \mu_{s, i} t_{s, i}) \quad \text{for } i = 1, \ldots, I, s \in \{\text{wt}, \text{MMR-}, \text{proof-}, \text{MMR-proof-}\}\\
\log10 \mu_{s, i} &\sim \mathcal{N}(\nu_{s}, \sigma_{s})  \quad \text{for } i = 1, \ldots, I, s \in \{\text{wt}, \text{MMR-}, \text{MMR-proof-}\} \\
\nu_s &\sim \mathcal{N}(-8, 3) \quad \text{for } s \in \{\text{wt}, \text{MMR-}, \text{MMR-proof-}\}\\
\sigma_s &\sim \text{t}^+(0, 3, 5) \quad \text{for } s \in \{\text{wt}, \text{MMR-}, \text{MMR-proof-}\}\\
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
#> The model will need to be run for a further 2323 updates.  This will
#> take approximately 0.1 seconds.
#> 
#> Calling the simulation using the rjags method...
#> Note: the model did not require adaptation
#> Running the model for 2323 iterations...
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
#>  1   2.88e-10   7.09e-11   1.50e-11   2.16e-10   4.19e-11   5.33e-11
#>  2   2.75e-10   5.42e-11   2.15e-11   1.65e-10   6.79e-11   7.21e-11
#>  3   3.28e-10   5.56e-11   2.47e-11   2.15e-10   7.30e-11   4.79e-11
#>  4   2.51e-10   6.71e-11   2.33e-11   2.17e-10   5.78e-11   4.93e-11
#>  5   2.98e-10   4.58e-11   2.07e-11   1.89e-10   5.34e-11   6.01e-11
#>  6   3.03e-10   5.29e-11   3.15e-11   1.99e-10   5.93e-11   4.95e-11
#>  7   2.76e-10   6.80e-11   2.02e-11   2.21e-10   3.87e-11   5.50e-11
#>  8   3.35e-10   7.45e-11   1.31e-11   2.16e-10   4.96e-11   4.06e-11
#>  9   2.68e-10   8.05e-11   1.63e-11   2.21e-10   4.08e-11   5.77e-11
#> 10   2.92e-10   5.24e-11   2.09e-11   2.02e-10   3.94e-11   3.54e-11
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
#>  1          0.884                1.29                      1.48          -10.6
#>  2          0.306                0.769                     1.26          -10.0
#>  3          0.570                0.871                     0.791         -10.1
#>  4          0.744                0.859                     1.21          -10.8
#>  5          0.450                0.759                     0.856         -10.4
#>  6          0.383                0.823                     0.686         -10.0
#>  7          0.354                0.971                     1.10          -10.3
#>  8          0.534                1.08                      0.717         -10.1
#>  9          0.541                1.92                      0.864         -10.4
#> 10          0.651                0.729                     0.904         -10.8
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
#>  1            127             29              7            104             19
#>  2            117             26              4            100             37
#>  3            150             15             12            122             39
#>  4            120             34              7            116             31
#>  5            132             24              8             98             25
#>  6            118             23             14            118             34
#>  7            107             42              7            105             15
#>  8            129             29              9            104             27
#>  9            107             35              5            125             26
#> 10            134             17              9            101             14
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
#>  1         1   2.88e-10   7.09e-11   1.50e-11   2.16e-10   4.19e-11   5.33e-11
#>  2         2   2.75e-10   5.42e-11   2.15e-11   1.65e-10   6.79e-11   7.21e-11
#>  3         3   3.28e-10   5.56e-11   2.47e-11   2.15e-10   7.30e-11   4.79e-11
#>  4         4   2.51e-10   6.71e-11   2.33e-11   2.17e-10   5.78e-11   4.93e-11
#>  5         5   2.98e-10   4.58e-11   2.07e-11   1.89e-10   5.34e-11   6.01e-11
#>  6         6   3.03e-10   5.29e-11   3.15e-11   1.99e-10   5.93e-11   4.95e-11
#>  7         7   2.76e-10   6.80e-11   2.02e-11   2.21e-10   3.87e-11   5.50e-11
#>  8         8   3.35e-10   7.45e-11   1.31e-11   2.16e-10   4.96e-11   4.06e-11
#>  9         9   2.68e-10   8.05e-11   1.63e-11   2.21e-10   4.08e-11   5.77e-11
#> 10        10   2.92e-10   5.24e-11   2.09e-11   2.02e-10   3.94e-11   3.54e-11
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
#> mu_wt[1]                   2.346732e-10  2.802420e-10  3.347770e-10
#> mu_wt[2]                   3.990693e-11  6.028459e-11  8.410767e-11
#> mu_wt[3]                   9.002371e-12  2.096186e-11  3.532989e-11
#> mu_wt[4]                   1.837094e-10  2.205730e-10  2.612236e-10
#> mu_wt[5]                   2.995582e-11  4.649660e-11  6.452962e-11
#> mu_wt[6]                   3.325154e-11  5.054758e-11  7.009606e-11
#> mu_proofminus[1]           5.743922e-08  6.453268e-08  7.306685e-08
#> mu_proofminus[2]           3.011400e-09  4.993006e-09  7.399738e-09
#> mu_proofminus[3]           1.785233e-10  7.226669e-10  1.484263e-09
#> mu_proofminus[4]           3.153756e-08  3.622161e-08  4.103700e-08
#> mu_proofminus[5]           2.105891e-09  3.290076e-09  4.714938e-09
#> mu_proofminus[6]           1.158123e-09  2.221478e-09  3.457355e-09
#> mu_MMRminus[1]             3.679103e-08  3.828998e-08  3.986927e-08
#> mu_MMRminus[2]             9.793107e-11  1.902041e-10  2.974671e-10
#> mu_MMRminus[3]             2.026470e-10  3.368017e-10  4.859358e-10
#> mu_MMRminus[4]             2.710569e-08  2.826797e-08  2.943475e-08
#> mu_MMRminus[5]             5.460320e-10  7.166602e-10  9.065256e-10
#> mu_MMRminus[6]             3.305588e-10  4.667129e-10  6.209079e-10
#> mu_MMRproofminus[1]        6.513876e-07  7.302575e-07  8.073133e-07
#> mu_MMRproofminus[2]        6.886659e-09  1.332352e-08  2.135654e-08
#> mu_MMRproofminus[3]        1.425708e-09  5.133347e-09  1.074243e-08
#> mu_MMRproofminus[4]        3.615616e-07  4.071825e-07  4.579919e-07
#> mu_MMRproofminus[5]        1.431806e-08  2.327643e-08  3.322611e-08
#> mu_MMRproofminus[6]        6.211514e-09  1.214355e-08  1.909955e-08
#> log10_mean_wt             -1.056596e+01 -1.011625e+01 -9.651590e+00
#> log10_mean_MMRminus       -9.598510e+00 -8.756848e+00 -7.911217e+00
#> log10_mean_MMRproofminus  -8.229006e+00 -7.388547e+00 -6.638753e+00
#> log10_sigma_wt             2.358002e-01  4.839379e-01  8.912410e-01
#> log10_sigma_MMRminus       5.459948e-01  9.634439e-01  1.587935e+00
#> log10_sigma_MMRproofminus  4.769053e-01  8.662879e-01  1.444495e+00
#> theta4                     6.932249e-02  8.164766e-02  9.400779e-02
#> log10_mean_prior          -1.393593e+01 -8.013697e+00 -2.148096e+00
#> log10_sigma_prior          1.076715e-04  4.177497e-01  1.502070e+00
#> log10_mu_prior            -1.431773e+01 -8.018451e+00 -2.088910e+00
#> theta4_prior               4.656475e-02  4.984014e-01  9.942811e-01
#> m_wt_pred[1]               8.700000e+01  1.170000e+02  1.460000e+02
#> m_wt_pred[2]               1.200000e+01  2.500000e+01  3.800000e+01
#> m_wt_pred[3]               1.000000e+00  8.000000e+00  1.600000e+01
#> m_wt_pred[4]               8.700000e+01  1.180000e+02  1.460000e+02
#> m_wt_pred[5]               1.300000e+01  2.500000e+01  3.900000e+01
#> m_wt_pred[6]               1.200000e+01  2.700000e+01  4.000000e+01
#> m_proofminus_pred[1]       1.690000e+02  2.030000e+02  2.420000e+02
#> m_proofminus_pred[2]       6.000000e+00  1.600000e+01  2.600000e+01
#> m_proofminus_pred[3]       0.000000e+00  2.000000e+00  6.000000e+00
#> m_proofminus_pred[4]       1.160000e+02  1.470000e+02  1.760000e+02
#> m_proofminus_pred[5]       5.000000e+00  1.300000e+01  2.200000e+01
#> m_proofminus_pred[6]       3.000000e+00  9.000000e+00  1.700000e+01
#> m_MMRminus_pred[1]         2.277000e+03  2.416000e+03  2.551000e+03
#> m_MMRminus_pred[2]         3.000000e+00  1.200000e+01  2.100000e+01
#> m_MMRminus_pred[3]         1.000000e+01  2.100000e+01  3.400000e+01
#> m_MMRminus_pred[4]         2.154000e+03  2.292000e+03  2.422000e+03
#> m_MMRminus_pred[5]         3.800000e+01  5.800000e+01  7.900000e+01
#> m_MMRminus_pred[6]         2.000000e+01  3.800000e+01  5.300000e+01
#> m_MMRproofminus_pred[1]    2.340000e+02  2.790000e+02  3.220000e+02
#> m_MMRproofminus_pred[2]    0.000000e+00  5.000000e+00  1.000000e+01
#> m_MMRproofminus_pred[3]    0.000000e+00  2.000000e+00  5.000000e+00
#> m_MMRproofminus_pred[4]    1.630000e+02  2.000000e+02  2.350000e+02
#> m_MMRproofminus_pred[5]    4.000000e+00  1.100000e+01  1.900000e+01
#> m_MMRproofminus_pred[6]    1.000000e+00  6.000000e+00  1.200000e+01
#> loglikelihood             -1.118248e+01 -1.001530e+01 -9.187250e+00
#>                                    Mean           SD         Mode        MCerr
#> mu_wt[1]                   2.812772e-10 2.577545e-11 2.017908e-10          Inf
#> mu_wt[2]                   6.091528e-11 1.145548e-11 2.387638e-11          Inf
#> mu_wt[3]                   2.158648e-11 6.999346e-12 3.933457e-12          Inf
#> mu_wt[4]                   2.208796e-10 2.001593e-11 1.486362e-10          Inf
#> mu_wt[5]                   4.699105e-11 8.964928e-12 1.976067e-11          Inf
#> mu_wt[6]                   5.110926e-11 9.540608e-12 2.365798e-11          Inf
#> mu_proofminus[1]           6.465736e-08 3.978274e-09           NA          Inf
#> mu_proofminus[2]           5.099006e-09 1.142334e-09 1.800858e-09          Inf
#> mu_proofminus[3]           7.794495e-10 3.575841e-10 6.721562e-11          Inf
#> mu_proofminus[4]           3.627262e-08 2.431600e-09           NA          Inf
#> mu_proofminus[5]           3.353608e-09 6.832621e-10 1.352331e-09          Inf
#> mu_proofminus[6]           2.273806e-09 5.918910e-10 6.756633e-10          Inf
#> mu_MMRminus[1]             3.828988e-08 7.827476e-10           NA          Inf
#> mu_MMRminus[2]             1.953392e-10 5.231576e-11 5.652844e-11          Inf
#> mu_MMRminus[3]             3.419484e-10 7.342530e-11 1.303270e-10          Inf
#> mu_MMRminus[4]             2.827041e-08 5.979387e-10           NA          Inf
#> mu_MMRminus[5]             7.210499e-10 9.302143e-11 4.293625e-10          Inf
#> mu_MMRminus[6]             4.708238e-10 7.545652e-11 2.641203e-10          Inf
#> mu_MMRproofminus[1]        7.306589e-07 4.038236e-08           NA 5.296771e-10
#> mu_MMRproofminus[2]        1.373455e-08 3.859558e-09           NA          Inf
#> mu_MMRproofminus[3]        5.549758e-09 2.565012e-09           NA          Inf
#> mu_MMRproofminus[4]        4.079119e-07 2.498331e-08           NA 3.177329e-10
#> mu_MMRproofminus[5]        2.364307e-08 4.905587e-09           NA          Inf
#> mu_MMRproofminus[6]        1.248525e-08 3.346616e-09           NA          Inf
#> log10_mean_wt             -1.011646e+01 2.266868e-01           NA 2.567195e-03
#> log10_mean_MMRminus       -8.759980e+00 4.246250e-01           NA 4.659154e-03
#> log10_mean_MMRproofminus  -7.391457e+00 3.968703e-01           NA 4.248721e-03
#> log10_sigma_wt             5.224393e-01 1.895571e-01           NA 2.635706e-03
#> log10_sigma_MMRminus       1.017511e+00 2.943200e-01           NA 3.849087e-03
#> log10_sigma_MMRproofminus  9.177937e-01 2.782809e-01           NA 3.798859e-03
#> theta4                     8.193506e-02 6.351782e-03           NA 8.758803e-05
#> log10_mean_prior          -8.022766e+00 3.003108e+00           NA 3.427575e-02
#> log10_sigma_prior          5.535563e-01 5.283986e-01           NA 5.843978e-03
#> log10_mu_prior            -8.031352e+00 3.111015e+00           NA 3.539430e-02
#> theta4_prior               4.990707e-01 2.885651e-01           NA 3.226256e-03
#> m_wt_pred[1]               1.171098e+02 1.516050e+01 1.180000e+02 1.694995e-01
#> m_wt_pred[2]               2.537112e+01 6.902126e+00 2.200000e+01 7.942418e-02
#> m_wt_pred[3]               8.943750e+00 4.147739e+00 7.000000e+00 4.729264e-02
#> m_wt_pred[4]               1.182765e+02 1.515090e+01 1.170000e+02 1.714767e-01
#> m_wt_pred[5]               2.517162e+01 6.960189e+00 2.500000e+01 7.997520e-02
#> m_wt_pred[6]               2.728900e+01 7.301174e+00 2.800000e+01 7.845608e-02
#> m_proofminus_pred[1]       2.036507e+02 1.884830e+01 2.050000e+02 2.063910e-01
#> m_proofminus_pred[2]       1.604725e+01 5.411555e+00 1.400000e+01 6.050302e-02
#> m_proofminus_pred[3]       2.477750e+00 1.924699e+00 1.000000e+00 2.189546e-02
#> m_proofminus_pred[4]       1.466678e+02 1.562297e+01 1.500000e+02 1.547348e-01
#> m_proofminus_pred[5]       1.351250e+01 4.557405e+00 1.200000e+01 5.095334e-02
#> m_proofminus_pred[6]       9.208000e+00 3.873708e+00 8.000000e+00 4.386633e-02
#> m_MMRminus_pred[1]         2.416615e+03 6.985288e+01 2.396000e+03 7.816610e-01
#> m_MMRminus_pred[2]         1.232125e+01 4.823791e+00 1.100000e+01 5.696549e-02
#> m_MMRminus_pred[3]         2.155262e+01 6.532577e+00 2.000000e+01 7.303643e-02
#> m_MMRminus_pred[4]         2.292565e+03 6.894234e+01 2.302000e+03 7.561185e-01
#> m_MMRminus_pred[5]         5.858237e+01 1.075118e+01 5.800000e+01 1.202019e-01
#> m_MMRminus_pred[6]         3.815525e+01 8.638228e+00 3.600000e+01 9.657832e-02
#> m_MMRproofminus_pred[1]    2.795113e+02 2.268610e+01 2.780000e+02 2.754245e-01
#> m_MMRproofminus_pred[2]    5.233875e+00 2.717462e+00 4.000000e+00 3.198744e-02
#> m_MMRproofminus_pred[3]    2.121625e+00 1.748275e+00 1.000000e+00 1.954631e-02
#> m_MMRproofminus_pred[4]    2.007508e+02 1.868031e+01 1.970000e+02 2.230841e-01
#> m_MMRproofminus_pred[5]    1.165488e+01 4.157604e+00 1.000000e+01 4.648343e-02
#> m_MMRproofminus_pred[6]    6.147250e+00 2.958743e+00 5.000000e+00 3.307976e-02
#> loglikelihood             -1.007639e+01 5.240954e-01           NA 6.286835e-03
#>                           MC%ofSD SSeff         AC.30      psrf
#> mu_wt[1]                      Inf     0 -2.687797e-03 1.0001725
#> mu_wt[2]                      Inf     0  1.006787e-03 1.0000740
#> mu_wt[3]                      Inf     0 -2.167658e-02 1.0004088
#> mu_wt[4]                      Inf     0  1.146405e-02 0.9999415
#> mu_wt[5]                      Inf     0 -2.711701e-03 1.0003869
#> mu_wt[6]                      Inf     0 -1.049941e-02 0.9998868
#> mu_proofminus[1]              Inf     0  2.244955e-02 1.0003085
#> mu_proofminus[2]              Inf     0 -1.236688e-03 1.0001456
#> mu_proofminus[3]              Inf     0 -1.599685e-02 1.0006681
#> mu_proofminus[4]              Inf     0 -7.479278e-03 1.0007780
#> mu_proofminus[5]              Inf     0 -1.523009e-02 1.0000576
#> mu_proofminus[6]              Inf     0 -1.341972e-02 1.0003531
#> mu_MMRminus[1]                Inf     0  2.791568e-03 1.0000355
#> mu_MMRminus[2]                Inf     0 -2.330483e-02 0.9999305
#> mu_MMRminus[3]                Inf     0 -9.778845e-03 1.0036895
#> mu_MMRminus[4]                Inf     0 -6.084406e-03 1.0005551
#> mu_MMRminus[5]                Inf     0  2.458489e-02 0.9999729
#> mu_MMRminus[6]                Inf     0  1.264988e-02 0.9998807
#> mu_MMRproofminus[1]           1.3  5812 -5.006192e-04 1.0006569
#> mu_MMRproofminus[2]           Inf     0  3.738864e-05 1.0000988
#> mu_MMRproofminus[3]           Inf     0 -1.126982e-02 1.0000134
#> mu_MMRproofminus[4]           1.3  6183  3.337419e-03 1.0005053
#> mu_MMRproofminus[5]           Inf     0 -1.095660e-03 0.9998931
#> mu_MMRproofminus[6]           Inf     0  4.909044e-03 1.0001033
#> log10_mean_wt                 1.1  7797  2.302772e-02 1.0015948
#> log10_mean_MMRminus           1.1  8306 -1.364432e-02 1.0004608
#> log10_mean_MMRproofminus      1.1  8725  4.044564e-03 0.9999127
#> log10_sigma_wt                1.4  5172  1.168536e-02 0.9999028
#> log10_sigma_MMRminus          1.3  5847  1.318297e-02 1.0005603
#> log10_sigma_MMRproofminus     1.4  5366  1.117450e-02 0.9999292
#> theta4                        1.4  5259  9.632220e-03 1.0010432
#> log10_mean_prior              1.1  7677 -1.687378e-02 1.0005964
#> log10_sigma_prior             1.1  8175  1.265504e-02 0.9998868
#> log10_mu_prior                1.1  7726 -1.399251e-02 1.0006042
#> theta4_prior                  1.1  8000  2.296586e-02 1.0003344
#> m_wt_pred[1]                  1.1  8000 -7.805156e-03 1.0005602
#> m_wt_pred[2]                  1.2  7552 -4.512795e-03 1.0009856
#> m_wt_pred[3]                  1.1  7692 -1.346787e-02 1.0001116
#> m_wt_pred[4]                  1.1  7807  8.330883e-03 1.0001671
#> m_wt_pred[5]                  1.1  7574 -7.617474e-03 1.0005443
#> m_wt_pred[6]                  1.1  8660 -1.059827e-02 0.9999457
#> m_proofminus_pred[1]          1.1  8340  1.296861e-02 0.9999261
#> m_proofminus_pred[2]          1.1  8000  4.853988e-03 1.0001337
#> m_proofminus_pred[3]          1.1  7727 -5.209716e-03 1.0007264
#> m_proofminus_pred[4]          1.0 10194 -2.066958e-02 0.9999794
#> m_proofminus_pred[5]          1.1  8000 -5.898708e-03 0.9999255
#> m_proofminus_pred[6]          1.1  7798 -1.417664e-02 1.0013876
#> m_MMRminus_pred[1]            1.1  7986 -4.552238e-03 1.0001284
#> m_MMRminus_pred[2]            1.2  7171 -1.859888e-02 1.0000904
#> m_MMRminus_pred[3]            1.1  8000  1.411063e-04 1.0008881
#> m_MMRminus_pred[4]            1.1  8314  1.192716e-02 1.0004729
#> m_MMRminus_pred[5]            1.1  8000  2.503522e-02 1.0003116
#> m_MMRminus_pred[6]            1.1  8000  4.364141e-03 0.9999336
#> m_MMRproofminus_pred[1]       1.2  6784 -8.193784e-03 1.0000442
#> m_MMRproofminus_pred[2]       1.2  7217 -2.587453e-03 1.0004707
#> m_MMRproofminus_pred[3]       1.1  8000 -6.461518e-03 0.9999992
#> m_MMRproofminus_pred[4]       1.2  7012  1.125915e-02 0.9999810
#> m_MMRproofminus_pred[5]       1.1  8000 -3.260838e-03 0.9999211
#> m_MMRproofminus_pred[6]       1.1  8000 -3.461880e-03 1.0004650
#> loglikelihood                 1.2  6950 -1.905349e-03 1.0000618
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
#>  1      117.             117           89          148             107
#>  2       25.4             25           13           40              20
#>  3        8.94             8            2           18               6
#>  4      118.             118           90          150             108
#>  5       25.2             25           13           40              20
#>  6       27.3             27           15           43              22
#>  7      204.             203          168          242             191
#>  8       16.0             16            7           28              12
#>  9        2.48             2            0            7               1
#> 10      147.             147          117          178             136
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
