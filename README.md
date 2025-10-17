<!--
SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>

SPDX-License-Identifier: AGPL-3.0-only
-->


<img align="right" src="https://avatars.githubusercontent.com/u/177750891?v=4" width="110">

[![DOI](https://img.shields.io/badge/DOI-10.48550%2FarXiv.2412.16317-blue.svg)](https://doi.org/10.48550/arXiv.2509.26274)
# Incomplete Bessel Function on top of EpsteinLib

This is a fork of [EpsteinLib](https://github.com/epsteinlib/epsteinlib) that includes an implementation of the incomplete Bessel function with the algorithm introduced in [1].

For $\nu\in\mathbb R$ and two $d$-dimensional vectors $\boldsymbol{k},\boldsymbol{r}\in\mathbb R^d$ the incomplete Bessel function is defined as
$$
G_{\nu}(\boldsymbol{k},\boldsymbol{r})=2\int_0^1 t^{-\nu} e^{-\pi \boldsymbol{k}^2/t^2} e^{-\pi \boldsymbol{r}^2 t^2}\, \frac{\mathrm d t}{t}.
$$

In this library, the function is implemented as
```c
double incomplete_bessel_g(
    double nu,
    unsigned int dim,
    const double *k,
    const double *r
);
```

## Examples
A Mathematica notebook replicating the results in [1] is included in `examples/mathematica/MicromagneticsReplicated.wls`.

## Installation
The library can be installed following the [EpsteinLib installation instructions](https://github.com/epsteinlib/epsteinlib/blob/main/README.md#installation).
## References

[1] Andreas A. Buchheit, Jonathan K. Busse, Torsten Keßler, and Filipp N. Rybakov.  “Zeta expansion for long-range interactions under periodic boundary conditions with applications to micromagnetics”.  arXiv preprint arXiv:2509.26274.
