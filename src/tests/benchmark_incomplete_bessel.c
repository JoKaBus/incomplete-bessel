// SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan.busse@dlr.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file benchmark_incomplete_bessel.c
 * @author Jonathan Busse
 * @date 06/06/2025
 * @section Description: Generate epsteinZeta(Reg) reference values
 */

#include <complex.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "epsteinZeta.h"

#define BASE_PATH "src/tests/csv"

#ifndef MAX_PATH_LENGTH
#define MAX_PATH_LENGTH 1024
#endif

/**
 * @brief Opens a file.
 * @param path Path to the file.
 * @param mode 'r' to read or 'w' to write.
 * @return FILE* Pointer to the opened file.
 * @note Exits the program if the file cannot be opened.
 */
FILE *open(char *path, char *mode) {
    FILE *file = fopen(path, mode);
    if (file) {
        return file;
    }
    printf("File '%s' does not exist.\n", path);
    exit(1); // NOLINT
}

/**
 * @brief Sorts an array of doubles in ascending order using bubble sort.
 * @param arr: The array to be sorted.
 * @param size: The size of the array.
 */
void sort(double *arr, int size) {
    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - i - 1; j++) {
            if (*(arr + j) > *(arr + j + 1)) {
                double temp = *(arr + j);
                *(arr + j) = *(arr + j + 1);
                *(arr + j + 1) = temp;
            }
        }
    }
}

/**
 * @brief Benchmarks the Epstein Zeta and regularized Epstein Zeta functions.
 * @param dim Dimension of the lattice.
 * @param a Array representing the lattice matrix.
 * @param x Array of x coordinates.
 * @param y Array of y coordinates.
 * @param zetaDataString Filename for storing Epstein Zeta function data.
 * @param zetaRegDataString Filename for storing regularized Epstein Zeta function
 * data.
 * @return  0 on successful execution.
 */
int benchmark(double nu, char zetaDataString[]) {
    FILE *zetaData = open(zetaDataString, "w");
    if (zetaData == NULL) {
        printf("%s\n", strerror(errno)); // NOLINT
        return 1;
    }
    int dim = 1;
    double x[1];
    double y[1];

    double complex zetaReg;
    double offset = 0.1;
    double elapsedTime;
    int iterations = 25;
    double *elapsedTimes = malloc(iterations * sizeof(double));
    clock_t timeStart;
    clock_t timeEnd;
    for (int i = 0; i < 10 + 1; i++) {
        for (int j = 0; j < 10 + 1; j++) {
            x[0] = i * 0.1 + offset;
            y[0] = j * 0.1 + offset;
            // zeta
            for (int n = 0; n < iterations; n++) {
                timeStart = clock();
                zetaReg = incomplete_bessel_g(nu, dim, x, y);
                timeEnd = clock();
                elapsedTimes[n] = ((double)(timeEnd - timeStart)) / CLOCKS_PER_SEC;
            }
            sort(elapsedTimes, iterations);
            elapsedTime = elapsedTimes[iterations / 2];
            fprintf(zetaData, "%.16lf,%.16lf,%.16lf,%.16lf,%.16lf\n", x[0], // NOLINT
                    y[0], creal(zetaReg), cimag(zetaReg), elapsedTime);
        }
        printf("k:\t %.16lf\t", x[0]);
        printf("r:\t %.16lf\t", y[0]);
        printf("bessel:\t%.16lf %+.16lf I, \t t: %.8f seconds\n", creal(zetaReg),
               cimag(zetaReg), elapsedTime);
    }
    free(elapsedTimes);
    if (fclose(zetaData) != 0) {
        return fprintf(stderr, "Error closing file: %d\n", errno);
    }
    return 0;
}

/*!
 * @brief Evaluating the incomplete bessel function over a grid of x and y points for
 * a constant nu indicated by the function name.
 * @return  0 on successful execution.
 */
int num4() {
    double nu = -4.;
    char zetaDataString[MAX_PATH_LENGTH];

    if (snprintf(zetaDataString, MAX_PATH_LENGTH, "%s/incomplete_bessel_%s.csv",
                 BASE_PATH, __func__) >= MAX_PATH_LENGTH) {
        return fprintf(stderr, "Error: filename too long\n");
    }
    printf("\n========== Benchmarking %s() ==========\n", __func__);
    return benchmark(nu, zetaDataString);
}

/*!
 * @brief Evaluating the incomplete bessel function over a grid of x and y points for
 * a constant nu indicated by the function name.
 * @return  0 on successful execution.
 */
int num2() {
    double nu = -4.;
    char zetaDataString[MAX_PATH_LENGTH];

    if (snprintf(zetaDataString, MAX_PATH_LENGTH, "%s/incomplete_bessel_%s.csv",
                 BASE_PATH, __func__) >= MAX_PATH_LENGTH) {
        return fprintf(stderr, "Error: filename too long\n");
    }
    printf("\n========== Benchmarking %s() ==========\n", __func__);
    return benchmark(nu, zetaDataString);
}

/*!
 * @brief Evaluating the incomplete bessel function over a grid of x and y points for
 * a constant nu indicated by the function name.
 * @return  0 on successful execution.
 */
int nu0() {
    double nu = -4.;
    char zetaDataString[MAX_PATH_LENGTH];

    if (snprintf(zetaDataString, MAX_PATH_LENGTH, "%s/incomplete_bessel_%s.csv",
                 BASE_PATH, __func__) >= MAX_PATH_LENGTH) {
        return fprintf(stderr, "Error: filename too long\n");
    }
    printf("\n========== Benchmarking %s() ==========\n", __func__);
    return benchmark(nu, zetaDataString);
}

/*!
 * @brief Evaluating the incomplete bessel function over a grid of x and y points for
 * a constant nu indicated by the function name.
 * @return  0 on successful execution.
 */
int nu2() {
    double nu = -4.;
    char zetaDataString[MAX_PATH_LENGTH];

    if (snprintf(zetaDataString, MAX_PATH_LENGTH, "%s/incomplete_bessel_%s.csv",
                 BASE_PATH, __func__) >= MAX_PATH_LENGTH) {
        return fprintf(stderr, "Error: filename too long\n");
    }
    printf("\n========== Benchmarking %s() ==========\n", __func__);
    return benchmark(nu, zetaDataString);
}

/**
 * @brief Main function to run all benchmark tests.
 * @return number of failed executions.
 */
int main() {
    int failed = 0;
    failed += num4();
    failed += num2();
    failed += nu0();
    failed += nu2();
    return failed;
}
