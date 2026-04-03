#ifndef CCUS_MATH_H
#define CCUS_MATH_H

#include <stdlib.h>

static inline void spsolve(int n, double *a, double *b, double *c, double *d, double *x) {
    double *c_star = (double *)malloc(n * sizeof(double));
    double *d_star = (double *)malloc(n * sizeof(double));

    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (int i = 1; i < n; i++) {
        double m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
        c_star[i] = c[i] * m;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
    }

    x[n - 1] = d_star[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = d_star[i] - c_star[i] * x[i + 1];
    }

    free(c_star);
    free(d_star);
}

#endif
