#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void takeStep(double alpha, double v, double x, double dx, double *v_next, double *alpha_next) {
    double vp = (alpha * (x - v) - 2.0 / x) * (x - v) / ((x - v) * (x - v) - 1);
    double alphap = (alpha - 2.0 / x * (x - v)) * (x - v) * alpha / ((x - v) * (x - v) - 1);
    double vt = v + vp * dx;
    double alphat = alpha + alphap * dx;
    double xt = x + dx;
    double vp1 = (alphat * (xt - vt) - 2.0 / xt) * (xt - vt) / ((xt - vt) * (xt - vt) - 1);
    double alphap1 = (alphat - 2.0 / xt * (xt - vt)) * (xt - vt) * alphat / ((xt - vt) * (xt - vt) - 1);
    vp = 0.5 * (vp + vp1);
    alphap = 0.5 * (alphap + alphap1);

    *v_next = v + vp * dx;
    *alpha_next = alpha + alphap * dx;
}

void getSolution(double A, int Nsteps, double **xvec, double **alphavec, double **vvec) {
    *xvec = (double *)malloc(Nsteps * sizeof(double));
    *alphavec = (double *)malloc(Nsteps * sizeof(double));
    *vvec = (double *)malloc(Nsteps * sizeof(double));

    double log_start = log10(2);
    double log_end = -12.0;
    double log_step = (log_end - log_start) / (Nsteps - 1);

    for (int i = 0; i < Nsteps; i++) {
        (*xvec)[i] = pow(10, log_start + i * log_step);
    }

    (*vvec)[0] = (-(A - 2) / (*xvec)[0]) - ((1 - A / 6) * (A - 2)) / ((*xvec)[0] * (*xvec)[0] * (*xvec)[0]);
    (*alphavec)[0] = A / ((*xvec)[0] * (*xvec)[0]) - A * (A - 2) / (2 * (*xvec)[0] * (*xvec)[0] * (*xvec)[0] * (*xvec)[0]);

    for (int i = 0; i < Nsteps - 1; i++) {
        double alpha = (*alphavec)[i];
        double x = (*xvec)[i];
        double v = (*vvec)[i];
        double dx = (*xvec)[i + 1] - (*xvec)[i];

        takeStep(alpha, v, x, dx, &(*vvec)[i + 1], &(*alphavec)[i + 1]);
    }
}

void saveToBinaryFile(const char *filename, double *xvec, double *alphavec, double *vvec, int Nsteps) {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        printf("Error opening file %s for writing.\n", filename);
        exit(1);
    }

    fwrite(&Nsteps, sizeof(int), 1, file);

    fwrite(xvec, sizeof(double), Nsteps, file);
    fwrite(alphavec, sizeof(double), Nsteps, file);
    fwrite(vvec, sizeof(double), Nsteps, file);

    fclose(file);
}

int main() {
    int Nsteps = 100;
    double A = 2.0000001;
    double *xvec, *alphavec, *vvec;

    getSolution(A, Nsteps, &xvec, &alphavec, &vvec);

    saveToBinaryFile("shusolution.bin", xvec, alphavec, vvec, Nsteps);

    free(xvec);
    free(alphavec);
    free(vvec);

    return 0;
}

