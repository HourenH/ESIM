## Overview

This directory contains necessary scripts and data to reproduce the
results and figures presented in *“Robust Regression for Spherical Data
in an Extrinsic Semiparametric Framework”*. Since the simulation studies
are computationally intensive, we provide the simulated data and
estimation results in the corresponding .RData files, which will be used
to produce relevant figures.

## Files in the Directory

-   **data**: This folder contains **gemas\_merge.csv** for real data
    application and .RData files that store simulated data and
    estimation results.

-   **plot.R**: this script produces the graphical illustrations of the
    asymptotic relative efficiency and standardized gross error
    sensitivity in Section 5. Graphical illustrations of the contours of
    the distributions on the sphere and shape of link function are also
    included.

-   **ESLFns.R**: this script provides estimation methods for all
    regression models considered in the numerical studies, including the
    extrinsic single index model, Fréchet single-index model with
    extrinsic distance and single-index quantile regression.

-   **ESIM\_contaminate.R**: this script reproduces simulations for
    Section 6.1. The relevant simulated data and estimation results are
    saved in **simu\_contaminate\_XXX.RData**, which can be used to
    generate the corresponding figures.

-   **ESIM\_shape.R**: this script reproduces simulations and relevant
    figures for Section 6.2. The relevant simulated data and estimation
    results are saved in **simu\_shape.RData**.

-   **ESIM\_power.R**: this script reproduces simulations for Section
    6.3.

-   **GEMAS\_sand.R**: this script reproduces the data application to
    **gemas\_merge.csv** in Section 6.4, with the estimation results
    saved in **GEMAS\_sand.RData**.

-   **ESIM\_cov.R** and **heavy\_tailed.R**: these scripts reproduces
    additional numerical results in the supplementary material, with the
    estimation results saved in **simu\_cov\_XXX.RData** and
    **simu\_heavy.RData**, respectively.

## Key Functions

We provide a brief overview of the main functions in **ESLFns.R** that
implement the proposed methodology. A detailed introduction to their
implementations can be found in Section S3 of the supplementary
material.

-   **ls\_esl**: This function performs the estimation for the extrinsic
    single-index model with the least squares.

-   **esl\_est**: This function performs the estimation for the
    extrinsic single-index model with the exponential squared loss. It
    involves iteratively reweighted least squares, and includes
    estimation of the tuning parameter *λ*. The function
    **cost\_lambda\_median** is used to estimate *λ* during this
    iterative process.

-   **siqr\_est**: This function performs the estimation for the
    single-index quantile regression with a quantile level of 0.5. It
    adapts the estimation approach proposed by Wu et al. (2010) and
    extends to multivariate response.

-   **fsim\_est**: This function performs the estimation for the
    single-index Fréchet regression proposed by Bhattacharjee and Müller
    (2023). For a detailed explanation of the methodology, please refer
    to Section S3.2 of the supplementary material. The code is adapted
    from <https://github.com/functionaldata/tFrechet>.

-   **simProject**: This function generate random samples following the
    scaled von Mises Fisher distribution. See Scealy and Wood (2019) for
    details.

## Required R packages:

-   Directional (6.9)
-   KernSmooth (2.23-24)
-   quantreg (5.99)
-   frechet (0.3.0)
-   MASS (7.3-60.2)

## Reference

-   Bhattacharjee, S. and Müller, H.-G. (2023). Single-index Fréchet
    regression. *The Annals of Statistics* **51**.(4), pp. 1770-1798.
-   Scealy, J. L. and Wood, A. T. A. (2019). Scaled von Mises-Fisher
    distributions and regression models for paleomagnetic directional
    data. *Journal of the American Statistical Association*
    **114**.(528), pp. 1547-1560.
-   Wu, T. Z., Yu, K., and Yu, Y.(2010). Single-index quantile
    regression. *Journal of Multivariate Analysis* **101**.(7),
    pp. 1607-1621.
