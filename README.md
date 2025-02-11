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
#>  1         -9.48         -10.4         -10.7         -9.70         -10.4
#>  2         -9.54         -10.3         -10.7         -9.68         -10.2
#>  3         -9.50         -10.3         -10.7         -9.65         -10.3
#>  4         -9.60         -10.3         -10.6         -9.73         -10.4
#>  5         -9.50         -10.3         -11.0         -9.61         -10.4
#>  6         -9.51         -10.4         -10.6         -9.70         -10.5
#>  7         -9.57         -10.2         -10.7         -9.72         -10.3
#>  8         -9.46         -10.3         -10.7         -9.67         -10.2
#>  9         -9.56         -10.3         -10.8         -9.65         -10.4
#> 10         -9.53         -10.1         -10.8         -9.67         -10.2
#> # ℹ 7,990 more rows
#> # ℹ 3 more variables: `log10_mu[6]` <dbl>, iteration <int>, chain_id <int>
```

``` r
extract_posterior_samples(fit_GCM_model, type = "hyperparameters")
#> # A tibble: 8,000 × 5
#>    log10_sigma log10_mean loglikelihood iteration chain_id
#>          <dbl>      <dbl>         <dbl>     <int>    <int>
#>  1       0.429     -10.1          -23.2         1        1
#>  2       0.393     -10.0          -17.3         2        1
#>  3       0.500     -10.0          -17.5         3        1
#>  4       0.351     -10.2          -20.5         4        1
#>  5       0.650     -10.2          -20.5         5        1
#>  6       0.560     -10.3          -20.1         6        1
#>  7       0.522     -10.3          -18.7         7        1
#>  8       0.541      -9.78         -20.6         8        1
#>  9       0.421     -10.1          -20.6         9        1
#> 10       0.543     -10.1          -18.5        10        1
#> # ℹ 7,990 more rows
```

``` r
extract_posterior_samples(fit_GCM_model, type = "predictive")
#> # A tibble: 8,000 × 8
#>    `m_pred[1]` `m_pred[2]` `m_pred[3]` `m_pred[4]` `m_pred[5]` `m_pred[6]`
#>          <dbl>       <dbl>       <dbl>       <dbl>       <dbl>       <dbl>
#>  1         156          24          13         104          33          52
#>  2         120          28           4         101          34          19
#>  3         130          16          11         118          32          26
#>  4         114          25          11          94          26          22
#>  5         126          23           5         128          22          42
#>  6         117          17          13         122          16          20
#>  7         108          18           5          96          16          35
#>  8         142          27           6         115          31          38
#>  9         101          17          11         115          20          41
#> 10         115          34           8         121          23          15
#> # ℹ 7,990 more rows
#> # ℹ 2 more variables: iteration <int>, chain_id <int>
```

``` r
extract_posterior_samples(fit_GCM_model, type = "prior")
#> # A tibble: 8,000 × 5
#>    log10_mean_prior log10_sigma_prior log10_mu_prior iteration chain_id
#>               <dbl>             <dbl>          <dbl>     <int>    <int>
#>  1            -9.39            0.689           -9.40         1        1
#>  2            -6.46            1.32            -6.63         2        1
#>  3           -10.9             1.13           -11.1          3        1
#>  4            -2.12            0.363           -1.78         4        1
#>  5            -7.12            0.215           -6.20         5        1
#>  6            -2.03            0.380           -1.92         6        1
#>  7            -7.49            0.0704          -7.09         7        1
#>  8            -7.71            0.316           -6.47         8        1
#>  9            -6.44            0.573           -6.99         9        1
#> 10            -9.24            1.01           -10.1         10        1
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
#>  1         1         156          24          13         104          33
#>  2         2         120          28           4         101          34
#>  3         3         130          16          11         118          32
#>  4         4         114          25          11          94          26
#>  5         5         126          23           5         128          22
#>  6         6         117          17          13         122          16
#>  7         7         108          18           5          96          16
#>  8         8         142          27           6         115          31
#>  9         9         101          17          11         115          20
#> 10        10         115          34           8         121          23
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
#> m_pred[1]          8.400000e+01 116.0000000 145.0000000 116.3262500 15.51380500
#> m_pred[2]          1.000000e+01  23.0000000  36.0000000  23.2427500  6.75159392
#> m_pred[3]          2.000000e+00   9.0000000  17.0000000   9.2725000  4.27004367
#> m_pred[4]          8.700000e+01 119.0000000 148.0000000 118.9926250 15.58571253
#> m_pred[5]          1.100000e+01  25.0000000  38.0000000  25.5791250  7.10544159
#> m_pred[6]          1.200000e+01  25.0000000  39.0000000  25.5393750  7.08738164
#> log10_mu[1]       -9.633936e+00  -9.5551013  -9.4731188  -9.5555620  0.04079907
#> log10_mu[2]       -1.043895e+01 -10.2596272 -10.0945929 -10.2627064  0.08844526
#> log10_mu[3]       -1.096206e+01 -10.6692109 -10.3957522 -10.6771988  0.14426978
#> log10_mu[4]       -9.732965e+00  -9.6542550  -9.5754254  -9.6551632  0.04037843
#> log10_mu[5]       -1.050346e+01 -10.3266799 -10.1622614 -10.3295729  0.08624568
#> log10_mu[6]       -1.049874e+01 -10.3293383 -10.1598902 -10.3311989  0.08617821
#> log10_mean        -1.059662e+01 -10.1198251  -9.6967623 -10.1239584  0.22569251
#> log10_sigma        2.475757e-01   0.4778156   0.8871072   0.5144408  0.17893496
#> log10_mean_prior  -1.387150e+01  -7.9960672  -2.1917508  -8.0053404  3.00955037
#> log10_sigma_prior  4.324250e-05   0.4197039   1.4834137   0.5488072  0.50117939
#> log10_mu_prior    -1.412475e+01  -7.9866476  -2.1466059  -8.0075285  3.06703190
#> loglikelihood     -2.264979e+01 -18.8503358 -16.4560548 -19.1659594  1.75053665
#>                   Mode        MCerr MC%ofSD SSeff         AC.20      psrf
#> m_pred[1]          113 0.1762052403     1.1  7752  0.0003942117 1.0004895
#> m_pred[2]           24 0.0754851148     1.1  8000  0.0062567154 0.9999371
#> m_pred[3]            7 0.0500401272     1.2  7282  0.0066770783 0.9999455
#> m_pred[4]          124 0.1752028213     1.1  7914  0.0139238030 1.0002275
#> m_pred[5]           22 0.0812187723     1.1  7654 -0.0128402131 1.0000546
#> m_pred[6]           24 0.0828365841     1.2  7320 -0.0002559283 0.9999426
#> log10_mu[1]         NA 0.0004733226     1.2  7430  0.0103179746 1.0000303
#> log10_mu[2]         NA 0.0010551321     1.2  7026  0.0048564099 1.0001461
#> log10_mu[3]         NA 0.0016467306     1.1  7675 -0.0140916043 1.0001867
#> log10_mu[4]         NA 0.0004653807     1.2  7528  0.0073987177 0.9999055
#> log10_mu[5]         NA 0.0010591990     1.2  6630 -0.0204885733 0.9999035
#> log10_mu[6]         NA 0.0010429496     1.2  6828  0.0028261199 1.0002037
#> log10_mean          NA 0.0025707997     1.1  7707  0.0105541997 1.0012217
#> log10_sigma         NA 0.0028683305     1.6  3892 -0.0090322414 1.0013394
#> log10_mean_prior    NA 0.0336477961     1.1  8000  0.0159578228 0.9999937
#> log10_sigma_prior   NA 0.0058415092     1.2  7361 -0.0129673344 1.0000805
#> log10_mu_prior      NA 0.0342904591     1.1  8000  0.0163794039 0.9999650
#> loglikelihood       NA 0.0227618386     1.3  5915 -0.0032295656 1.0001307
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
#> 1           1      116.             116           87          149
#> 2           2       23.2             23           11           38
#> 3           3        9.27             9            2           19
#> 4           4      119.             119           90          151
#> 5           5       25.6             25           13           41
#> 6           6       25.5             25           13           41
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
#> The model will need to be run for a further 2262 updates.  This will
#> take approximately 0.3 seconds.
#> 
#> Calling the simulation using the rjags method...
#> Note: the model did not require adaptation
#> Running the model for 2262 iterations...
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
#>  1   2.90e-10   6.93e-11   2.23e-11   2.61e-10   5.15e-11   6.69e-11
#>  2   2.67e-10   7.15e-11   1.57e-11   2.16e-10   4.95e-11   5.20e-11
#>  3   2.55e-10   7.00e-11   2.25e-11   1.77e-10   5.52e-11   4.53e-11
#>  4   2.58e-10   5.92e-11   2.66e-11   2.56e-10   4.07e-11   4.43e-11
#>  5   2.72e-10   5.74e-11   6.32e-12   2.06e-10   4.38e-11   5.95e-11
#>  6   2.29e-10   5.96e-11   1.38e-11   2.30e-10   3.76e-11   6.78e-11
#>  7   2.43e-10   5.83e-11   1.81e-11   2.17e-10   4.79e-11   6.98e-11
#>  8   2.48e-10   4.37e-11   1.49e-11   1.93e-10   4.53e-11   4.69e-11
#>  9   2.32e-10   7.30e-11   2.26e-11   2.30e-10   5.75e-11   6.66e-11
#> 10   3.12e-10   6.12e-11   2.17e-11   2.22e-10   4.04e-11   3.91e-11
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
#>  1          0.494                1.22                      0.753         -9.91
#>  2          0.423                0.946                     0.717        -10.3 
#>  3          0.510                1.20                      0.947         -9.99
#>  4          0.517                0.940                     0.912        -10.2 
#>  5          0.709                0.664                     1.09         -10.2 
#>  6          0.706                0.478                     0.769        -10.1 
#>  7          0.478                1.11                      0.585         -9.87
#>  8          0.452                1.05                      0.893        -10.4 
#>  9          0.294                1.01                      0.993        -10.0 
#> 10          0.533                0.993                     1.06         -10.3 
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
#>  1            131             25             12            143             32
#>  2            106             24              8            125             34
#>  3            120             25             13             83             37
#>  4            100             17              7            140             29
#>  5            101             23              3            112             25
#>  6            111             25              3            129             10
#>  7            101             29              7            120             20
#>  8             96             21              6            118             28
#>  9             85             25             10            125             38
#> 10            133             32              8            140             21
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
#>  1         1   2.90e-10   6.93e-11   2.23e-11   2.61e-10   5.15e-11   6.69e-11
#>  2         2   2.67e-10   7.15e-11   1.57e-11   2.16e-10   4.95e-11   5.20e-11
#>  3         3   2.55e-10   7.00e-11   2.25e-11   1.77e-10   5.52e-11   4.53e-11
#>  4         4   2.58e-10   5.92e-11   2.66e-11   2.56e-10   4.07e-11   4.43e-11
#>  5         5   2.72e-10   5.74e-11   6.32e-12   2.06e-10   4.38e-11   5.95e-11
#>  6         6   2.29e-10   5.96e-11   1.38e-11   2.30e-10   3.76e-11   6.78e-11
#>  7         7   2.43e-10   5.83e-11   1.81e-11   2.17e-10   4.79e-11   6.98e-11
#>  8         8   2.48e-10   4.37e-11   1.49e-11   1.93e-10   4.53e-11   4.69e-11
#>  9         9   2.32e-10   7.30e-11   2.26e-11   2.30e-10   5.75e-11   6.66e-11
#> 10        10   3.12e-10   6.12e-11   2.17e-11   2.22e-10   4.04e-11   3.91e-11
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
#> mu_wt[1]                   2.309943e-10  2.805059e-10  3.329021e-10
#> mu_wt[2]                   3.897531e-11  6.048328e-11  8.436733e-11
#> mu_wt[3]                   9.191965e-12  2.086995e-11  3.529250e-11
#> mu_wt[4]                   1.826677e-10  2.206570e-10  2.613945e-10
#> mu_wt[5]                   3.056929e-11  4.644606e-11  6.569831e-11
#> mu_wt[6]                   3.328631e-11  5.035869e-11  7.009627e-11
#> mu_proofminus[1]           5.709543e-08  6.466295e-08  7.264147e-08
#> mu_proofminus[2]           2.830018e-09  5.011139e-09  7.291215e-09
#> mu_proofminus[3]           1.825496e-10  7.271745e-10  1.527759e-09
#> mu_proofminus[4]           3.167575e-08  3.632244e-08  4.120231e-08
#> mu_proofminus[5]           2.046311e-09  3.288081e-09  4.659237e-09
#> mu_proofminus[6]           1.180330e-09  2.220343e-09  3.436939e-09
#> mu_MMRminus[1]             3.680053e-08  3.829363e-08  3.989757e-08
#> mu_MMRminus[2]             9.773440e-11  1.894193e-10  2.949948e-10
#> mu_MMRminus[3]             2.041180e-10  3.357053e-10  4.882021e-10
#> mu_MMRminus[4]             2.711783e-08  2.826939e-08  2.942962e-08
#> mu_MMRminus[5]             5.402958e-10  7.166859e-10  9.041290e-10
#> mu_MMRminus[6]             3.296964e-10  4.694978e-10  6.190088e-10
#> mu_MMRproofminus[1]        6.509692e-07  7.288979e-07  8.068831e-07
#> mu_MMRproofminus[2]        6.751637e-09  1.331252e-08  2.152598e-08
#> mu_MMRproofminus[3]        1.157064e-09  5.148699e-09  1.058671e-08
#> mu_MMRproofminus[4]        3.597541e-07  4.072072e-07  4.582039e-07
#> mu_MMRproofminus[5]        1.416244e-08  2.311879e-08  3.289411e-08
#> mu_MMRproofminus[6]        6.635907e-09  1.219343e-08  1.923059e-08
#> log10_mean_wt             -1.058463e+01 -1.010945e+01 -9.684179e+00
#> log10_mean_MMRminus       -9.588177e+00 -8.761531e+00 -7.862296e+00
#> log10_mean_MMRproofminus  -8.146931e+00 -7.387789e+00 -6.644603e+00
#> log10_sigma_wt             2.381231e-01  4.795321e-01  8.684787e-01
#> log10_sigma_MMRminus       5.402337e-01  9.639503e-01  1.568893e+00
#> log10_sigma_MMRproofminus  4.899533e-01  8.628168e-01  1.455412e+00
#> theta4                     7.005629e-02  8.200440e-02  9.461223e-02
#> log10_mean_prior          -1.402098e+01 -7.889006e+00 -2.245692e+00
#> log10_sigma_prior          1.117846e-05  4.214700e-01  1.495824e+00
#> log10_mu_prior            -1.400222e+01 -7.877955e+00 -1.809920e+00
#> theta4_prior               1.409309e-02  4.928566e-01  9.616617e-01
#> m_wt_pred[1]               8.900000e+01  1.170000e+02  1.470000e+02
#> m_wt_pred[2]               1.100000e+01  2.500000e+01  3.800000e+01
#> m_wt_pred[3]               2.000000e+00  9.000000e+00  1.700000e+01
#> m_wt_pred[4]               8.700000e+01  1.180000e+02  1.470000e+02
#> m_wt_pred[5]               1.100000e+01  2.500000e+01  3.800000e+01
#> m_wt_pred[6]               1.300000e+01  2.700000e+01  4.100000e+01
#> m_proofminus_pred[1]       1.670000e+02  2.040000e+02  2.400000e+02
#> m_proofminus_pred[2]       6.000000e+00  1.600000e+01  2.600000e+01
#> m_proofminus_pred[3]       0.000000e+00  2.000000e+00  6.000000e+00
#> m_proofminus_pred[4]       1.170000e+02  1.460000e+02  1.780000e+02
#> m_proofminus_pred[5]       5.000000e+00  1.300000e+01  2.200000e+01
#> m_proofminus_pred[6]       2.000000e+00  9.000000e+00  1.600000e+01
#> m_MMRminus_pred[1]         2.280000e+03  2.417000e+03  2.549000e+03
#> m_MMRminus_pred[2]         3.000000e+00  1.200000e+01  2.100000e+01
#> m_MMRminus_pred[3]         9.000000e+00  2.100000e+01  3.400000e+01
#> m_MMRminus_pred[4]         2.156000e+03  2.292000e+03  2.421000e+03
#> m_MMRminus_pred[5]         3.700000e+01  5.800000e+01  7.800000e+01
#> m_MMRminus_pred[6]         2.000000e+01  3.800000e+01  5.400000e+01
#> m_MMRproofminus_pred[1]    2.320000e+02  2.790000e+02  3.200000e+02
#> m_MMRproofminus_pred[2]    0.000000e+00  5.000000e+00  1.000000e+01
#> m_MMRproofminus_pred[3]    0.000000e+00  2.000000e+00  5.000000e+00
#> m_MMRproofminus_pred[4]    1.620000e+02  2.000000e+02  2.360000e+02
#> m_MMRproofminus_pred[5]    3.000000e+00  1.100000e+01  1.900000e+01
#> m_MMRproofminus_pred[6]    0.000000e+00  6.000000e+00  1.100000e+01
#> loglikelihood             -1.114472e+01 -1.002294e+01 -9.115662e+00
#>                                    Mean           SD         Mode        MCerr
#> mu_wt[1]                   2.815859e-10 2.616874e-11 2.001872e-10          Inf
#> mu_wt[2]                   6.125202e-11 1.177185e-11 2.397534e-11          Inf
#> mu_wt[3]                   2.162194e-11 6.871799e-12 5.275226e-12          Inf
#> mu_wt[4]                   2.210649e-10 2.026063e-11 1.457461e-10          Inf
#> mu_wt[5]                   4.700987e-11 9.037119e-12 2.096982e-11          Inf
#> mu_wt[6]                   5.105577e-11 9.654325e-12 2.148515e-11          Inf
#> mu_proofminus[1]           6.476676e-08 3.967211e-09           NA          Inf
#> mu_proofminus[2]           5.109620e-09 1.163694e-09 1.817134e-09          Inf
#> mu_proofminus[3]           7.866943e-10 3.732443e-10 9.409103e-11          Inf
#> mu_proofminus[4]           3.636549e-08 2.447466e-09           NA          Inf
#> mu_proofminus[5]           3.342300e-09 6.791346e-10 1.190747e-09          Inf
#> mu_proofminus[6]           2.268498e-09 5.942629e-10 6.851348e-10          Inf
#> mu_MMRminus[1]             3.829828e-08 7.839010e-10           NA          Inf
#> mu_MMRminus[2]             1.943095e-10 5.155628e-11 4.537015e-11          Inf
#> mu_MMRminus[3]             3.412454e-10 7.334341e-11 1.537163e-10          Inf
#> mu_MMRminus[4]             2.826851e-08 5.939342e-10           NA          Inf
#> mu_MMRminus[5]             7.212265e-10 9.315607e-11 4.233517e-10          Inf
#> mu_MMRminus[6]             4.720628e-10 7.540831e-11 2.323957e-10          Inf
#> mu_MMRproofminus[1]        7.292538e-07 4.007250e-08           NA 4.996235e-10
#> mu_MMRproofminus[2]        1.366219e-08 3.898626e-09           NA          Inf
#> mu_MMRproofminus[3]        5.571841e-09 2.627062e-09           NA          Inf
#> mu_MMRproofminus[4]        4.075044e-07 2.524732e-08           NA 3.223669e-10
#> mu_MMRproofminus[5]        2.350401e-08 4.833592e-09           NA          Inf
#> mu_MMRproofminus[6]        1.246062e-08 3.332777e-09           NA          Inf
#> log10_mean_wt             -1.010900e+01 2.226481e-01           NA 2.520035e-03
#> log10_mean_MMRminus       -8.752981e+00 4.299264e-01           NA 4.656639e-03
#> log10_mean_MMRproofminus  -7.393362e+00 3.795588e-01           NA 4.292235e-03
#> log10_sigma_wt             5.141040e-01 1.825862e-01           NA 2.696311e-03
#> log10_sigma_MMRminus       1.013846e+00 2.873829e-01           NA 3.712449e-03
#> log10_sigma_MMRproofminus  9.154800e-01 2.777162e-01           NA 4.084346e-03
#> theta4                     8.225287e-02 6.376853e-03           NA 9.008896e-05
#> log10_mean_prior          -7.944172e+00 3.002813e+00           NA 3.357247e-02
#> log10_sigma_prior          5.468959e-01 4.957421e-01           NA 5.508573e-03
#> log10_mu_prior            -7.944234e+00 3.099529e+00           NA 3.416662e-02
#> theta4_prior               4.956480e-01 2.886095e-01           NA 3.226752e-03
#> m_wt_pred[1]               1.173451e+02 1.529623e+01 1.130000e+02 1.706066e-01
#> m_wt_pred[2]               2.559487e+01 7.057636e+00 2.400000e+01 7.918795e-02
#> m_wt_pred[3]               9.040375e+00 4.184846e+00 7.000000e+00 4.678800e-02
#> m_wt_pred[4]               1.184486e+02 1.528543e+01 1.160000e+02 1.708963e-01
#> m_wt_pred[5]               2.518525e+01 7.043251e+00 2.300000e+01 7.874594e-02
#> m_wt_pred[6]               2.723525e+01 7.370088e+00 2.700000e+01 8.353917e-02
#> m_proofminus_pred[1]       2.042259e+02 1.884248e+01 2.020000e+02 2.148914e-01
#> m_proofminus_pred[2]       1.600250e+01 5.359182e+00 1.400000e+01 5.991748e-02
#> m_proofminus_pred[3]       2.504125e+00 1.953192e+00 1.000000e+00 2.201841e-02
#> m_proofminus_pred[4]       1.470419e+02 1.566589e+01 1.410000e+02 1.808370e-01
#> m_proofminus_pred[5]       1.349013e+01 4.594022e+00 1.100000e+01 5.136273e-02
#> m_proofminus_pred[6]       9.216500e+00 3.846359e+00 9.000000e+00 4.300360e-02
#> m_MMRminus_pred[1]         2.418256e+03 6.921071e+01 2.436000e+03 7.737992e-01
#> m_MMRminus_pred[2]         1.221575e+01 4.714788e+00 1.200000e+01 6.260775e-02
#> m_MMRminus_pred[3]         2.155300e+01 6.493340e+00 1.900000e+01 7.259775e-02
#> m_MMRminus_pred[4]         2.291424e+03 6.812901e+01 2.292000e+03 7.869651e-01
#> m_MMRminus_pred[5]         5.848613e+01 1.063776e+01 6.000000e+01 1.242786e-01
#> m_MMRminus_pred[6]         3.824650e+01 8.821959e+00 3.600000e+01 9.898615e-02
#> m_MMRproofminus_pred[1]    2.794525e+02 2.256442e+01 2.770000e+02 2.743445e-01
#> m_MMRproofminus_pred[2]    5.297000e+00 2.730425e+00 4.000000e+00 3.168405e-02
#> m_MMRproofminus_pred[3]    2.095500e+00 1.767773e+00 1.000000e+00 1.976430e-02
#> m_MMRproofminus_pred[4]    2.004311e+02 1.893138e+01 1.930000e+02 2.264571e-01
#> m_MMRproofminus_pred[5]    1.164213e+01 4.182522e+00 1.200000e+01 4.676201e-02
#> m_MMRproofminus_pred[6]    6.136000e+00 2.971508e+00 5.000000e+00 3.322247e-02
#> loglikelihood             -1.007988e+01 5.256296e-01           NA 6.117975e-03
#>                           MC%ofSD SSeff         AC.30      psrf
#> mu_wt[1]                      Inf     0  0.0046212852 1.0000406
#> mu_wt[2]                      Inf     0 -0.0056275496 0.9999663
#> mu_wt[3]                      Inf     0 -0.0161750759 1.0004612
#> mu_wt[4]                      Inf     0  0.0098513548 1.0000878
#> mu_wt[5]                      Inf     0 -0.0013303571 1.0005317
#> mu_wt[6]                      Inf     0  0.0038057734 0.9999585
#> mu_proofminus[1]              Inf     0  0.0044417409 1.0010225
#> mu_proofminus[2]              Inf     0  0.0038151509 1.0008234
#> mu_proofminus[3]              Inf     0  0.0073011409 0.9999727
#> mu_proofminus[4]              Inf     0  0.0065036921 1.0007172
#> mu_proofminus[5]              Inf     0  0.0010145543 1.0003267
#> mu_proofminus[6]              Inf     0  0.0121476278 0.9999685
#> mu_MMRminus[1]                Inf     0 -0.0070796212 0.9999577
#> mu_MMRminus[2]                Inf     0  0.0032082541 0.9999974
#> mu_MMRminus[3]                Inf     0  0.0102162951 0.9999026
#> mu_MMRminus[4]                Inf     0  0.0121604582 1.0000462
#> mu_MMRminus[5]                Inf     0  0.0121105795 0.9998965
#> mu_MMRminus[6]                Inf     0  0.0229107941 1.0001090
#> mu_MMRproofminus[1]           1.2  6433  0.0042041134 1.0003128
#> mu_MMRproofminus[2]           Inf     0 -0.0038490104 1.0000459
#> mu_MMRproofminus[3]           Inf     0  0.0068485546 0.9999788
#> mu_MMRproofminus[4]           1.3  6134 -0.0165126148 0.9999614
#> mu_MMRproofminus[5]           Inf     0 -0.0102419837 1.0002414
#> mu_MMRproofminus[6]           Inf     0 -0.0019180205 0.9999926
#> log10_mean_wt                 1.1  7806 -0.0036110306 1.0017593
#> log10_mean_MMRminus           1.1  8524  0.0107059243 0.9999840
#> log10_mean_MMRproofminus      1.1  7820  0.0052681647 0.9998982
#> log10_sigma_wt                1.5  4586  0.0048622282 1.0009198
#> log10_sigma_MMRminus          1.3  5992 -0.0257435483 1.0002227
#> log10_sigma_MMRproofminus     1.5  4623 -0.0089777455 1.0002161
#> theta4                        1.4  5010  0.0055618974 1.0010929
#> log10_mean_prior              1.1  8000  0.0155880609 0.9998902
#> log10_sigma_prior             1.1  8099 -0.0153546047 0.9998872
#> log10_mu_prior                1.1  8230  0.0124435154 1.0000316
#> theta4_prior                  1.1  8000  0.0075210777 1.0008356
#> m_wt_pred[1]                  1.1  8039 -0.0057424156 1.0003010
#> m_wt_pred[2]                  1.1  7943  0.0001019550 1.0000775
#> m_wt_pred[3]                  1.1  8000  0.0008452583 0.9999272
#> m_wt_pred[4]                  1.1  8000  0.0014227227 1.0005124
#> m_wt_pred[5]                  1.1  8000 -0.0151809878 0.9999365
#> m_wt_pred[6]                  1.1  7783  0.0135063448 1.0002766
#> m_proofminus_pred[1]          1.1  7688 -0.0065790523 0.9999992
#> m_proofminus_pred[2]          1.1  8000 -0.0054937285 1.0000833
#> m_proofminus_pred[3]          1.1  7869  0.0205573848 1.0000878
#> m_proofminus_pred[4]          1.2  7505  0.0155526710 1.0007597
#> m_proofminus_pred[5]          1.1  8000  0.0001407080 1.0005880
#> m_proofminus_pred[6]          1.1  8000  0.0111215212 1.0006923
#> m_MMRminus_pred[1]            1.1  8000  0.0063529066 1.0000050
#> m_MMRminus_pred[2]            1.3  5671  0.0061060370 1.0002145
#> m_MMRminus_pred[3]            1.1  8000 -0.0067180216 0.9999997
#> m_MMRminus_pred[4]            1.2  7495 -0.0006107319 1.0004200
#> m_MMRminus_pred[5]            1.2  7327 -0.0086648256 0.9999114
#> m_MMRminus_pred[6]            1.1  7943 -0.0020774680 1.0009344
#> m_MMRproofminus_pred[1]       1.2  6765 -0.0012994499 1.0000252
#> m_MMRproofminus_pred[2]       1.2  7426 -0.0114327361 1.0005053
#> m_MMRproofminus_pred[3]       1.1  8000  0.0229920946 1.0002798
#> m_MMRproofminus_pred[4]       1.2  6989 -0.0158274452 1.0000386
#> m_MMRproofminus_pred[5]       1.1  8000 -0.0130330135 0.9999657
#> m_MMRproofminus_pred[6]       1.1  8000 -0.0034351568 1.0001782
#> loglikelihood                 1.2  7381  0.0149903998 1.0000089
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
#>  1      117.             117         89            148           107  
#>  2       25.6             25         13             41            20.8
#>  3        9.04             9          2             19             6  
#>  4      118.             118         90.0          150           108  
#>  5       25.2             25         13             40            20  
#>  6       27.2             27         14             43            22  
#>  7      204.             204        168            243           191  
#>  8       16.0             16          7             28            12  
#>  9        2.50             2          0              7             1  
#> 10      147.             146        118            179           136  
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
