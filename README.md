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
  - A functioning latex installation
  - cmake

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
input_data = read.delim(file = "data-raw/SAT_mcmc_chain_labmut_R3610MMR-3610_inputfile.csv", sep = " ") %>% 
  rename(mutation_id = context_id,
         m = m_sc, 
         n = n_c,
         t = t_s, 
         mutation_label = labmut) %>% 
  mutate(strain = gsub("MMR-", "MMRminus", strain))  %>%
  mutate(strain = gsub("WT3610", "wt", strain)) %>% 
  mutate(strain = gsub("polC_mutL", "MMRproofminus", strain)) %>% 
  mutate(strain = gsub("polC", "proofminus", strain))

input_data_onestrain = input_data %>% 
  filter(strain == first(strain)) 
```

The minimum data required to run the analysis is a data frame with the
following columns:

``` r
minimal_input_data_onestrain = input_data %>% 
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

``` r
fit_GCM_model = input_data_onestrain %>% 
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
#>  1         -9.52         -10.1         -10.7         -9.67         -10.4
#>  2         -9.54         -10.3         -10.6         -9.65         -10.2
#>  3         -9.55         -10.3         -10.7         -9.63         -10.3
#>  4         -9.57         -10.4         -10.9         -9.70         -10.3
#>  5         -9.54         -10.2         -10.6         -9.64         -10.2
#>  6         -9.58         -10.3         -10.6         -9.69         -10.3
#>  7         -9.60         -10.3         -10.6         -9.68         -10.2
#>  8         -9.59         -10.2         -10.6         -9.63         -10.4
#>  9         -9.48         -10.3         -10.7         -9.68         -10.3
#> 10         -9.55         -10.3         -10.8         -9.62         -10.4
#> # ℹ 7,990 more rows
#> # ℹ 3 more variables: `log10_mu[6]` <dbl>, iteration <int>, chain_id <int>
```

``` r
extract_posterior_samples(fit_GCM_model, type = "hyperparameters")
#> # A tibble: 8,000 × 5
#>    log10_sigma log10_mean loglikelihood iteration chain_id
#>          <dbl>      <dbl>         <dbl>     <int>    <int>
#>  1       0.503     -10.2          -18.4         1        1
#>  2       0.428     -10.2          -18.3         2        1
#>  3       0.461     -10.3          -17.0         3        1
#>  4       0.570     -10.1          -19.3         4        1
#>  5       0.751      -9.53         -17.3         5        1
#>  6       0.481     -10.2          -18.0         6        1
#>  7       0.492     -10.2          -18.1         7        1
#>  8       0.328     -10.1          -17.2         8        1
#>  9       0.453     -10.4          -18.2         9        1
#> 10       0.377     -10.3          -17.7        10        1
#> # ℹ 7,990 more rows
```

``` r
extract_posterior_samples(fit_GCM_model, type = "predictive")
#> # A tibble: 8,000 × 8
#>    `m_pred[1]` `m_pred[2]` `m_pred[3]` `m_pred[4]` `m_pred[5]` `m_pred[6]`
#>          <dbl>       <dbl>       <dbl>       <dbl>       <dbl>       <dbl>
#>  1         119          42           6         120          17          26
#>  2         102          20           7         113          32          25
#>  3         104          24           8         136          19          22
#>  4         109          25           2         113          29          28
#>  5         120          32           7         119          41          23
#>  6         117          18          10         108          25          39
#>  7          97          18          12         124          23          22
#>  8         126          26          12         117          17          23
#>  9         161          19           7         114          36          24
#> 10         124          12           8         117          16          26
#> # ℹ 7,990 more rows
#> # ℹ 2 more variables: iteration <int>, chain_id <int>
```

``` r
extract_posterior_samples(fit_GCM_model, type = "prior")
#> # A tibble: 8,000 × 5
#>    log10_mean_prior log10_sigma_prior log10_mu_prior iteration chain_id
#>               <dbl>             <dbl>          <dbl>     <int>    <int>
#>  1            -6.32            0.0568          -6.98         1        1
#>  2            -5.83            0.345           -6.75         2        1
#>  3            -7.43            0.313           -7.53         3        1
#>  4           -13.9             1.38           -13.5          4        1
#>  5           -10.4             2.69           -10.5          5        1
#>  6            -9.01            0.259           -8.56         6        1
#>  7            -7.97            0.258           -8.07         7        1
#>  8           -11.6             0.195          -11.9          8        1
#>  9            -2.47            0.155           -3.23         9        1
#> 10           -10.3             0.134          -11.0         10        1
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
#>  1         1         119          42           6         120          17
#>  2         2         102          20           7         113          32
#>  3         3         104          24           8         136          19
#>  4         4         109          25           2         113          29
#>  5         5         120          32           7         119          41
#>  6         6         117          18          10         108          25
#>  7         7          97          18          12         124          23
#>  8         8         126          26          12         117          17
#>  9         9         161          19           7         114          36
#> 10        10         124          12           8         117          16
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
  traceplot
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

``` r
fit_GCM_model %>% 
  summary
#>                         Lower95      Median     Upper95        Mean          SD
#> m_pred[1]          8.600000e+01 116.0000000 146.0000000 116.8637500 15.42694671
#> m_pred[2]          1.000000e+01  23.0000000  36.0000000  23.1215000  6.76348326
#> m_pred[3]          2.000000e+00   9.0000000  17.0000000   9.1456250  4.24684269
#> m_pred[4]          8.800000e+01 119.0000000 147.0000000 119.2938750 15.24458304
#> m_pred[5]          1.200000e+01  25.0000000  39.0000000  25.4291250  7.15130194
#> m_pred[6]          1.100000e+01  25.0000000  38.0000000  25.4785000  7.06797583
#> log10_mu[1]       -9.631881e+00  -9.5541419  -9.4739386  -9.5544863  0.04057941
#> log10_mu[2]       -1.044534e+01 -10.2605788 -10.0950172 -10.2634517  0.08952878
#> log10_mu[3]       -1.096282e+01 -10.6707332 -10.4010667 -10.6798175  0.14535859
#> log10_mu[4]       -9.732319e+00  -9.6542977  -9.5782116  -9.6543513  0.03960798
#> log10_mu[5]       -1.050012e+01 -10.3297875 -10.1646802 -10.3326401  0.08605969
#> log10_mu[6]       -1.050092e+01 -10.3284469 -10.1651612 -10.3306282  0.08498760
#> log10_mean        -1.057998e+01 -10.1219731  -9.6794782 -10.1222264  0.22690121
#> log10_sigma        2.384310e-01   0.4825386   0.8880191   0.5211335  0.19077198
#> log10_mean_prior  -1.361292e+01  -8.0130327  -2.0651334  -7.9872935  2.98271723
#> log10_sigma_prior  8.920633e-05   0.4077116   1.4853073   0.5426845  0.51806345
#> log10_mu_prior    -1.409692e+01  -8.0403123  -2.2062202  -7.9926409  3.03183015
#> loglikelihood     -2.242250e+01 -18.8235532 -16.4458200 -19.1260931  1.72083624
#>                   Mode        MCerr MC%ofSD SSeff        AC.20      psrf
#> m_pred[1]          113 0.1773113463     1.1  7570  0.025608936 1.0001351
#> m_pred[2]           20 0.0756180417     1.1  8000  0.013380226 0.9999760
#> m_pred[3]            7 0.0502354585     1.2  7147  0.016255714 1.0001684
#> m_pred[4]          118 0.1745416916     1.1  7628 -0.019456397 1.0000509
#> m_pred[5]           23 0.0840972386     1.2  7231  0.010315396 1.0000333
#> m_pred[6]           24 0.0827498435     1.2  7296 -0.001397931 1.0002942
#> log10_mu[1]         NA 0.0004732438     1.2  7353  0.008381707 0.9999805
#> log10_mu[2]         NA 0.0011282623     1.3  6297  0.003464174 1.0000125
#> log10_mu[3]         NA 0.0018129078     1.2  6429 -0.000977797 1.0000954
#> log10_mu[4]         NA 0.0004484250     1.1  7802 -0.006816225 1.0000127
#> log10_mu[5]         NA 0.0010591795     1.2  6602 -0.009233665 0.9999905
#> log10_mu[6]         NA 0.0010450737     1.2  6613 -0.021217756 1.0002867
#> log10_mean          NA 0.0025368327     1.1  8000  0.018040875 1.0002007
#> log10_sigma         NA 0.0034715395     1.8  3020  0.006664236 1.0002023
#> log10_mean_prior    NA 0.0333477924     1.1  8000  0.003574647 1.0001390
#> log10_sigma_prior   NA 0.0057491111     1.1  8120 -0.009531894 1.0002215
#> log10_mu_prior      NA 0.0341461606     1.1  7884  0.005351398 1.0002839
#> loglikelihood       NA 0.0219309641     1.3  6157  0.001864577 1.0000794
```

``` r
fit_GCM_model %>% 
  plot_prior_posterior
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

``` r
posterior_predictive_one_strain(fit_GCM_model)
#> Joining with `by = join_by(mutation_id)`
#> # A tibble: 6 × 17
#>   mutation_id m_pred_mean m_pred_median m_pred_infCI m_pred_supCI
#>         <int>       <dbl>         <dbl>        <dbl>        <dbl>
#> 1           1      117.             116         89            149
#> 2           2       23.1             23         12.0           38
#> 3           3        9.15             9          2             19
#> 4           4      119.             119         91            150
#> 5           5       25.4             25         13             41
#> 6           6       25.5             25         13             40
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

``` r
minimal_input_data = input_data %>% 
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
fit_MMRsaturation_model = EstimateMusMMRsaturation(input_data)
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
#> The model will need to be run for a further 2605 updates.  This will
#> take approximately 0.2 seconds.
#> 
#> Calling the simulation using the rjags method...
#> Note: the model did not require adaptation
#> Running the model for 2605 iterations...
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
#>  1   2.71e-10   5.10e-11   1.44e-11   2.06e-10   4.15e-11   4.99e-11
#>  2   3.47e-10   6.55e-11   1.78e-11   2.48e-10   3.97e-11   4.29e-11
#>  3   2.79e-10   6.15e-11   2.11e-11   2.10e-10   5.46e-11   4.85e-11
#>  4   2.83e-10   5.04e-11   2.61e-11   2.62e-10   4.20e-11   6.90e-11
#>  5   3.00e-10   8.35e-11   2.32e-11   2.18e-10   4.99e-11   5.07e-11
#>  6   2.63e-10   4.13e-11   3.49e-11   2.17e-10   3.25e-11   4.51e-11
#>  7   2.77e-10   7.41e-11   1.86e-11   2.23e-10   3.92e-11   5.79e-11
#>  8   2.58e-10   6.19e-11   2.66e-11   2.58e-10   5.37e-11   5.14e-11
#>  9   2.95e-10   8.04e-11   2.57e-11   1.96e-10   5.53e-11   5.49e-11
#> 10   3.01e-10   7.42e-11   1.60e-11   2.41e-10   4.54e-11   5.25e-11
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
#>  1          0.543                0.889                     0.913        -10.4 
#>  2          0.849                1.26                      0.888         -9.95
#>  3          0.731                1.17                      0.654        -10.1 
#>  4          0.680                1.02                      0.719         -9.79
#>  5          0.307                0.684                     0.639        -10.2 
#>  6          0.678                0.577                     0.879        -10.2 
#>  7          0.601                0.652                     0.860        -10.1 
#>  8          0.401                1.19                      0.672        -10.1 
#>  9          0.223                1.24                      0.678         -9.99
#> 10          0.305                1.27                      0.623        -10.3 
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
#>  1            141             29              5            102             27
#>  2            173             24              4            128             24
#>  3            118             40              6            107             29
#>  4            117             21              8            143             23
#>  5            118             37              7            103             33
#>  6            115             20             11            101             11
#>  7            134             25              5            124             22
#>  8            111             25             14            139             33
#>  9            115             46              8             91             24
#> 10            144             29              3            125             28
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
#>  1         1   2.71e-10   5.10e-11   1.44e-11   2.06e-10   4.15e-11   4.99e-11
#>  2         2   3.47e-10   6.55e-11   1.78e-11   2.48e-10   3.97e-11   4.29e-11
#>  3         3   2.79e-10   6.15e-11   2.11e-11   2.10e-10   5.46e-11   4.85e-11
#>  4         4   2.83e-10   5.04e-11   2.61e-11   2.62e-10   4.20e-11   6.90e-11
#>  5         5   3.00e-10   8.35e-11   2.32e-11   2.18e-10   4.99e-11   5.07e-11
#>  6         6   2.63e-10   4.13e-11   3.49e-11   2.17e-10   3.25e-11   4.51e-11
#>  7         7   2.77e-10   7.41e-11   1.86e-11   2.23e-10   3.92e-11   5.79e-11
#>  8         8   2.58e-10   6.19e-11   2.66e-11   2.58e-10   5.37e-11   5.14e-11
#>  9         9   2.95e-10   8.04e-11   2.57e-11   1.96e-10   5.53e-11   5.49e-11
#> 10        10   3.01e-10   7.42e-11   1.60e-11   2.41e-10   4.54e-11   5.25e-11
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
  traceplot
```

<img src="man/figures/README-unnamed-chunk-23-1.png" width="100%" />

``` r
fit_MMRsaturation_model %>% 
  summary
#> Calculating summary statistics...
#> Calculating the Gelman-Rubin statistic for 60 variables....
#>                                 Lower95        Median       Upper95
#> mu_wt[1]                   2.319737e-10  2.805517e-10  3.333281e-10
#> mu_wt[2]                   4.020673e-11  6.013787e-11  8.464137e-11
#> mu_wt[3]                   9.523458e-12  2.096973e-11  3.582148e-11
#> mu_wt[4]                   1.818733e-10  2.201120e-10  2.610554e-10
#> mu_wt[5]                   3.048891e-11  4.626423e-11  6.514124e-11
#> mu_wt[6]                   3.391805e-11  5.045271e-11  7.016785e-11
#> mu_proofminus[1]           5.689379e-08  6.466949e-08  7.224961e-08
#> mu_proofminus[2]           3.025031e-09  5.012935e-09  7.471035e-09
#> mu_proofminus[3]           1.722407e-10  7.318250e-10  1.519398e-09
#> mu_proofminus[4]           3.167510e-08  3.626778e-08  4.115815e-08
#> mu_proofminus[5]           2.125790e-09  3.320537e-09  4.752809e-09
#> mu_proofminus[6]           1.130549e-09  2.212426e-09  3.430021e-09
#> mu_MMRminus[1]             3.681826e-08  3.829187e-08  3.987386e-08
#> mu_MMRminus[2]             1.007386e-10  1.890770e-10  2.960109e-10
#> mu_MMRminus[3]             2.041492e-10  3.369981e-10  4.856705e-10
#> mu_MMRminus[4]             2.709532e-08  2.826071e-08  2.939754e-08
#> mu_MMRminus[5]             5.395248e-10  7.172403e-10  8.985671e-10
#> mu_MMRminus[6]             3.309799e-10  4.640555e-10  6.220455e-10
#> mu_MMRproofminus[1]        6.571874e-07  7.290655e-07  8.091372e-07
#> mu_MMRproofminus[2]        6.967082e-09  1.328221e-08  2.158305e-08
#> mu_MMRproofminus[3]        1.210026e-09  5.235682e-09  1.077328e-08
#> mu_MMRproofminus[4]        3.625000e-07  4.074955e-07  4.586586e-07
#> mu_MMRproofminus[5]        1.471051e-08  2.325813e-08  3.370842e-08
#> mu_MMRproofminus[6]        6.281081e-09  1.206100e-08  1.921639e-08
#> log10_mean_wt             -1.058437e+01 -1.011088e+01 -9.682178e+00
#> log10_mean_MMRminus       -9.569483e+00 -8.758980e+00 -7.927872e+00
#> log10_mean_MMRproofminus  -8.134109e+00 -7.382121e+00 -6.580994e+00
#> log10_sigma_wt             2.332630e-01  4.821714e-01  8.808577e-01
#> log10_sigma_MMRminus       5.262165e-01  9.590418e-01  1.569814e+00
#> log10_sigma_MMRproofminus  4.826015e-01  8.611181e-01  1.454626e+00
#> theta4                     7.003505e-02  8.186476e-02  9.514560e-02
#> log10_mean_prior          -1.404947e+01 -8.110199e+00 -2.286707e+00
#> log10_sigma_prior          2.050898e-04  4.178389e-01  1.486792e+00
#> log10_mu_prior            -1.413478e+01 -8.112506e+00 -2.052095e+00
#> theta4_prior               5.483254e-02  5.062389e-01  9.999854e-01
#> m_wt_pred[1]               8.500000e+01  1.170000e+02  1.450000e+02
#> m_wt_pred[2]               1.100000e+01  2.500000e+01  3.800000e+01
#> m_wt_pred[3]               2.000000e+00  8.000000e+00  1.700000e+01
#> m_wt_pred[4]               8.700000e+01  1.180000e+02  1.470000e+02
#> m_wt_pred[5]               1.200000e+01  2.500000e+01  3.800000e+01
#> m_wt_pred[6]               1.200000e+01  2.700000e+01  4.000000e+01
#> m_proofminus_pred[1]       1.650000e+02  2.030000e+02  2.380000e+02
#> m_proofminus_pred[2]       5.000000e+00  1.600000e+01  2.600000e+01
#> m_proofminus_pred[3]       0.000000e+00  2.000000e+00  6.000000e+00
#> m_proofminus_pred[4]       1.180000e+02  1.460000e+02  1.790000e+02
#> m_proofminus_pred[5]       4.000000e+00  1.300000e+01  2.200000e+01
#> m_proofminus_pred[6]       2.000000e+00  9.000000e+00  1.600000e+01
#> m_MMRminus_pred[1]         2.283000e+03  2.417000e+03  2.550000e+03
#> m_MMRminus_pred[2]         3.000000e+00  1.200000e+01  2.100000e+01
#> m_MMRminus_pred[3]         9.000000e+00  2.100000e+01  3.400000e+01
#> m_MMRminus_pred[4]         2.166000e+03  2.292000e+03  2.430000e+03
#> m_MMRminus_pred[5]         3.600000e+01  5.800000e+01  7.800000e+01
#> m_MMRminus_pred[6]         2.200000e+01  3.800000e+01  5.500000e+01
#> m_MMRproofminus_pred[1]    2.360000e+02  2.800000e+02  3.230000e+02
#> m_MMRproofminus_pred[2]    0.000000e+00  5.000000e+00  1.000000e+01
#> m_MMRproofminus_pred[3]    0.000000e+00  2.000000e+00  5.000000e+00
#> m_MMRproofminus_pred[4]    1.630000e+02  2.000000e+02  2.360000e+02
#> m_MMRproofminus_pred[5]    3.000000e+00  1.100000e+01  1.900000e+01
#> m_MMRproofminus_pred[6]    1.000000e+00  6.000000e+00  1.200000e+01
#> loglikelihood             -1.110589e+01 -1.001143e+01 -9.099817e+00
#>                                    Mean           SD         Mode        MCerr
#> mu_wt[1]                   2.812666e-10 2.609427e-11 1.955085e-10          Inf
#> mu_wt[2]                   6.094489e-11 1.152010e-11 2.690342e-11          Inf
#> mu_wt[3]                   2.165710e-11 6.928706e-12 4.880494e-12          Inf
#> mu_wt[4]                   2.208716e-10 2.038386e-11 1.520201e-10          Inf
#> mu_wt[5]                   4.704917e-11 9.014261e-12 1.815064e-11          Inf
#> mu_wt[6]                   5.118714e-11 9.465765e-12 2.178262e-11          Inf
#> mu_proofminus[1]           6.475921e-08 3.976790e-09           NA          Inf
#> mu_proofminus[2]           5.112709e-09 1.156008e-09 1.881662e-09          Inf
#> mu_proofminus[3]           7.880761e-10 3.664566e-10 3.921636e-11          Inf
#> mu_proofminus[4]           3.634283e-08 2.443337e-09           NA          Inf
#> mu_proofminus[5]           3.361434e-09 6.885851e-10 1.319377e-09          Inf
#> mu_proofminus[6]           2.269283e-09 6.002067e-10 6.900236e-10          Inf
#> mu_MMRminus[1]             3.829306e-08 7.761509e-10           NA          Inf
#> mu_MMRminus[2]             1.941415e-10 5.125821e-11 6.649447e-11          Inf
#> mu_MMRminus[3]             3.424925e-10 7.336915e-11 1.305927e-10          Inf
#> mu_MMRminus[4]             2.826877e-08 5.926066e-10           NA          Inf
#> mu_MMRminus[5]             7.206391e-10 9.223960e-11 4.012494e-10          Inf
#> mu_MMRminus[6]             4.690869e-10 7.581991e-11 2.020792e-10          Inf
#> mu_MMRproofminus[1]        7.300599e-07 3.906134e-08           NA 5.054774e-10
#> mu_MMRproofminus[2]        1.371709e-08 3.902208e-09           NA          Inf
#> mu_MMRproofminus[3]        5.595958e-09 2.616311e-09           NA          Inf
#> mu_MMRproofminus[4]        4.077563e-07 2.489885e-08           NA 3.201302e-10
#> mu_MMRproofminus[5]        2.364332e-08 4.914758e-09           NA          Inf
#> mu_MMRproofminus[6]        1.241422e-08 3.388701e-09           NA          Inf
#> log10_mean_wt             -1.011003e+01 2.222820e-01           NA 2.408500e-03
#> log10_mean_MMRminus       -8.759222e+00 4.139085e-01           NA 4.562723e-03
#> log10_mean_MMRproofminus  -7.387007e+00 3.900132e-01           NA 4.360480e-03
#> log10_sigma_wt             5.192794e-01 1.852319e-01           NA 2.660443e-03
#> log10_sigma_MMRminus       1.009511e+00 2.869501e-01           NA 3.627894e-03
#> log10_sigma_MMRproofminus  9.141557e-01 2.765872e-01           NA 3.944370e-03
#> theta4                     8.214723e-02 6.441805e-03           NA 9.056633e-05
#> log10_mean_prior          -8.086714e+00 3.015130e+00           NA 3.169297e-02
#> log10_sigma_prior          5.446705e-01 5.093683e-01           NA 5.694910e-03
#> log10_mu_prior            -8.085269e+00 3.088948e+00           NA 3.412793e-02
#> theta4_prior               5.050175e-01 2.862870e-01           NA 3.200786e-03
#> m_wt_pred[1]               1.173326e+02 1.532692e+01 1.150000e+02 1.766590e-01
#> m_wt_pred[2]               2.537425e+01 6.943808e+00 2.400000e+01 8.101092e-02
#> m_wt_pred[3]               8.956375e+00 4.206654e+00 7.000000e+00 4.703182e-02
#> m_wt_pred[4]               1.184005e+02 1.541359e+01 1.170000e+02 1.842355e-01
#> m_wt_pred[5]               2.523075e+01 6.885850e+00 2.500000e+01 7.698614e-02
#> m_wt_pred[6]               2.741737e+01 7.282242e+00 2.600000e+01 8.287663e-02
#> m_proofminus_pred[1]       2.037340e+02 1.905284e+01 2.020000e+02 2.209837e-01
#> m_proofminus_pred[2]       1.603637e+01 5.437095e+00 1.300000e+01 6.219545e-02
#> m_proofminus_pred[3]       2.477375e+00 1.954153e+00 1.000000e+00 2.184810e-02
#> m_proofminus_pred[4]       1.469670e+02 1.574149e+01 1.440000e+02 1.759953e-01
#> m_proofminus_pred[5]       1.361250e+01 4.661283e+00 1.300000e+01 5.553725e-02
#> m_proofminus_pred[6]       9.177000e+00 3.890284e+00 9.000000e+00 4.349470e-02
#> m_MMRminus_pred[1]         2.417387e+03 6.902064e+01 2.411000e+03 7.716743e-01
#> m_MMRminus_pred[2]         1.221612e+01 4.754444e+00 1.000000e+01 5.727315e-02
#> m_MMRminus_pred[3]         2.164075e+01 6.583757e+00 1.900000e+01 7.360864e-02
#> m_MMRminus_pred[4]         2.292356e+03 6.786253e+01 2.285000e+03 7.587261e-01
#> m_MMRminus_pred[5]         5.837237e+01 1.074481e+01 5.700000e+01 1.201306e-01
#> m_MMRminus_pred[6]         3.808212e+01 8.779110e+00 3.500000e+01 1.008800e-01
#> m_MMRproofminus_pred[1]    2.797666e+02 2.250055e+01 2.830000e+02 2.677558e-01
#> m_MMRproofminus_pred[2]    5.238875e+00 2.743220e+00 5.000000e+00 3.255478e-02
#> m_MMRproofminus_pred[3]    2.157875e+00 1.752736e+00 1.000000e+00 1.959619e-02
#> m_MMRproofminus_pred[4]    2.005884e+02 1.880124e+01 2.000000e+02 2.248316e-01
#> m_MMRproofminus_pred[5]    1.167775e+01 4.194711e+00 1.000000e+01 4.796742e-02
#> m_MMRproofminus_pred[6]    6.152375e+00 2.977753e+00 5.000000e+00 3.329229e-02
#> loglikelihood             -1.007070e+01 5.251627e-01           NA 6.236941e-03
#>                           MC%ofSD SSeff         AC.30      psrf
#> mu_wt[1]                      Inf     0 -8.331133e-03 0.9999356
#> mu_wt[2]                      Inf     0 -2.957342e-03 1.0003352
#> mu_wt[3]                      Inf     0  1.650429e-02 0.9999776
#> mu_wt[4]                      Inf     0  1.950419e-02 1.0004964
#> mu_wt[5]                      Inf     0  2.864284e-02 1.0001427
#> mu_wt[6]                      Inf     0  1.056009e-02 1.0001458
#> mu_proofminus[1]              Inf     0 -1.368480e-02 1.0002003
#> mu_proofminus[2]              Inf     0  3.129713e-03 1.0008406
#> mu_proofminus[3]              Inf     0 -8.760059e-03 1.0012076
#> mu_proofminus[4]              Inf     0 -6.824937e-03 1.0000452
#> mu_proofminus[5]              Inf     0  2.965672e-03 0.9999302
#> mu_proofminus[6]              Inf     0  3.471410e-03 1.0007780
#> mu_MMRminus[1]                Inf     0  1.027597e-02 1.0004235
#> mu_MMRminus[2]                Inf     0  5.988439e-03 1.0002984
#> mu_MMRminus[3]                Inf     0  7.590804e-03 1.0001931
#> mu_MMRminus[4]                Inf     0  1.857384e-02 0.9998771
#> mu_MMRminus[5]                Inf     0 -1.835841e-02 0.9999426
#> mu_MMRminus[6]                Inf     0  3.717450e-02 1.0004058
#> mu_MMRproofminus[1]           1.3  5972 -1.765419e-02 1.0004080
#> mu_MMRproofminus[2]           Inf     0 -1.186330e-02 1.0004088
#> mu_MMRproofminus[3]           Inf     0 -1.013736e-02 1.0011583
#> mu_MMRproofminus[4]           1.3  6049  8.219139e-03 0.9999857
#> mu_MMRproofminus[5]           Inf     0  4.138820e-03 0.9999020
#> mu_MMRproofminus[6]           Inf     0  1.864982e-02 0.9999336
#> log10_mean_wt                 1.1  8518 -2.385154e-02 1.0008587
#> log10_mean_MMRminus           1.1  8229  4.115650e-03 0.9999396
#> log10_mean_MMRproofminus      1.1  8000 -6.770418e-03 0.9998762
#> log10_sigma_wt                1.4  4848  1.469630e-02 0.9999581
#> log10_sigma_MMRminus          1.3  6256  2.167384e-02 1.0003239
#> log10_sigma_MMRproofminus     1.4  4917 -2.764774e-02 1.0010621
#> theta4                        1.4  5059 -1.298522e-02 1.0003120
#> log10_mean_prior              1.1  9051  3.302569e-03 1.0000848
#> log10_sigma_prior             1.1  8000 -1.110158e-02 1.0004010
#> log10_mu_prior                1.1  8192  3.844898e-03 1.0000139
#> theta4_prior                  1.1  8000 -1.100546e-02 1.0001324
#> m_wt_pred[1]                  1.2  7527  4.070673e-03 1.0009727
#> m_wt_pred[2]                  1.2  7347 -8.470873e-03 1.0003681
#> m_wt_pred[3]                  1.1  8000  9.537739e-03 0.9999535
#> m_wt_pred[4]                  1.2  6999  1.344584e-02 1.0001003
#> m_wt_pred[5]                  1.1  8000  1.915474e-02 1.0004049
#> m_wt_pred[6]                  1.1  7721  1.286802e-02 1.0002480
#> m_proofminus_pred[1]          1.2  7434 -2.286425e-02 1.0004527
#> m_proofminus_pred[2]          1.1  7642  5.517018e-03 1.0002499
#> m_proofminus_pred[3]          1.1  8000 -7.112717e-03 1.0011270
#> m_proofminus_pred[4]          1.1  8000 -9.203639e-03 1.0000249
#> m_proofminus_pred[5]          1.2  7044  1.555872e-02 1.0003228
#> m_proofminus_pred[6]          1.1  8000 -1.504450e-03 1.0002646
#> m_MMRminus_pred[1]            1.1  8000  1.463674e-02 1.0003327
#> m_MMRminus_pred[2]            1.2  6891  7.807436e-03 0.9999776
#> m_MMRminus_pred[3]            1.1  8000  1.683589e-03 0.9999445
#> m_MMRminus_pred[4]            1.1  8000  5.574090e-03 0.9999582
#> m_MMRminus_pred[5]            1.1  8000 -4.117237e-03 1.0001322
#> m_MMRminus_pred[6]            1.1  7573  8.607904e-03 1.0005060
#> m_MMRproofminus_pred[1]       1.2  7062  1.660871e-04 1.0000284
#> m_MMRproofminus_pred[2]       1.2  7101 -2.387523e-05 1.0000524
#> m_MMRproofminus_pred[3]       1.1  8000 -8.834158e-03 1.0009183
#> m_MMRproofminus_pred[4]       1.2  6993  1.181100e-02 1.0007702
#> m_MMRproofminus_pred[5]       1.1  7647  1.079666e-02 1.0007339
#> m_MMRproofminus_pred[6]       1.1  8000 -2.334791e-03 1.0002756
#> loglikelihood                 1.2  7090  1.252002e-03 1.0003206
```

``` r
fit_MMRsaturation_model %>% 
  plot_prior_posterior
```

<img src="man/figures/README-unnamed-chunk-25-1.png" width="100%" />

``` r
fit_MMRsaturation_model %>% 
  posterior_predictive
#> Joining with `by = join_by(strain, mutation_id)`
#> # A tibble: 24 × 18
#>    m_pred_mean m_pred_median m_pred_infCI m_pred_supCI m_pred_infquart
#>          <dbl>         <dbl>        <dbl>        <dbl>           <dbl>
#>  1      117.             117           89        149               107
#>  2       25.4             25           13         40                20
#>  3        8.96             8            2         18.0               6
#>  4      118.             118           89        150               108
#>  5       25.2             25           13         40                20
#>  6       27.4             27           15         43                22
#>  7      204.             203          168        242               190
#>  8       16.0             16            7         28                12
#>  9        2.48             2            0          7                 1
#> 10      147.             146          118        179               136
#> # ℹ 14 more rows
#> # ℹ 13 more variables: m_pred_supquart <dbl>, colnam <chr>, strain <chr>,
#> #   mutation_id <dbl>, genotype <chr>, mutation_label <chr>, nposinref <int>,
#> #   ngeninMA <dbl>, bps.n <int>, context <chr>, m <int>, n <int>, t <dbl>
```

``` r
fit_MMRsaturation_model %>% 
  plot_posterior_predictive
#> Joining with `by = join_by(strain, mutation_id)`
```

<img src="man/figures/README-unnamed-chunk-27-1.png" width="100%" />
