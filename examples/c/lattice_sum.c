// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file lattice_sum.c
 * @brief Calculates a function value of the incomplete Bessel function.
 *
 * Minimal working example for the Epstein Zeta Library.
 * If the library is installed, compile with `gcc -o lattice_sum lattice_sum.c
 * -lm -lepsteinZeta` If the library is not installed, compile with `gcc -o
 * lattice_sum lattice_sum.c -lm -L/path/to/library -lepsteinZeta
 * -I/path/to/include`
 * @author Andreas Buchheit, Jonathan Busse and Ruben Gutendorf.
 * @date 06/06/2024
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
// If library is not installed, use #include "/path/to/epsteinZeta.h"
#include "epsteinZeta.h"

/**
 * @brief calculate a function value of the incomplete Bessel function and compare to
 * precomputed value.
 *
 * @return true, if the difference to precomputed value is smaller than
 * 10**(-14)
 */
int main() {

    int dim = 1;
    double nu = 2.1;
    double x[] = {1.2};
    double y[] = {1.3};

    double besselRef =
        0.00003616792891719726; // Reference value obtained by integration
    double bessel = creal(incomplete_bessel_g(nu, dim, x, y));

    printf("C library Bessel value:\t\t %.16lf\n", creal(bessel));
    printf("Integration reference value:\t %.16lf\n", besselRef);
    printf("Relative error:\t\t\t +%.2e\n",
           fabs(besselRef - bessel) / fabs(besselRef));

    return fabs(bessel - besselRef) > pow(10, -14);
}
