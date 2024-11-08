#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SCRIPT_R 3360.0
#define T 10.0
#define B 8.86
#define M 1.989e30
#define G 6.674e-11
#define N_STEPS 10000
#define AU_TO_METERS 1.49e11

// Trapezoidal integration (linear)
double trapIntegrateLinear(double (*f)(double), double xmin, double xmax, int N) {
    double step = (xmax - xmin) / (N - 1);
    double integral = 0.0;
    for (int i = 0; i < N - 1; i++) {
        double x1 = xmin + i * step;
        double x2 = xmin + (i + 1) * step;
        double avg = (f(x1) + f(x2)) / 2.0;
        integral += avg * (x2 - x1);
    }
    return integral;
}

// Trapezoidal integration (logarithmic)
double trapIntegrateLog(double (*f)(double), double xmin, double xmax, int N) {
    double log_min = log10(xmin);
    double log_max = log10(xmax);
    double integral = 0.0;
    for (int i = 0; i < N - 1; i++) {
        double s1 = pow(10, log_min + i * (log_max - log_min) / (N - 1));
        double s2 = pow(10, log_min + (i + 1) * (log_max - log_min) / (N - 1));
        double avg = (f(s1) + f(s2)) / 2.0;
        integral += avg * (s2 - s1);
    }
    return integral;
}

// Interpolation for a given x-value using linear interpolation
double interpolate(double *x, double *y, int size, double xi) {
    if (xi <= x[0]) return y[0];
    if (xi >= x[size - 1]) return y[size - 1];
    int i = 0;
    while (xi > x[i + 1]) i++;
    double x0 = x[i], x1 = x[i + 1];
    double y0 = y[i], y1 = y[i + 1];
    return y0 + (y1 - y0) * (xi - x0) / (x1 - x0);
}

// Load data from a binary file
int loadData(const char *filename, double *xvec, double *alphavec, double *vvec) {
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        fprintf(stderr, "Error opening file.\n");
        return -1;
    }

    int Nsteps;
    fread(&Nsteps, sizeof(int), 1, file);

    fread(xvec, sizeof(double), Nsteps, file);
    fread(alphavec, sizeof(double), Nsteps, file);
    fread(vvec, sizeof(double), Nsteps, file);

    fclose(file);
    return Nsteps;
}

// Define `getAlpha` using interpolated data
double getAlpha(double x, double *xvec, double *alphavec, int Nsteps) {
    if (x > 1) return 2 * pow(x, -2);
    return interpolate(xvec, alphavec, Nsteps, x);
}

// Define `getV` using interpolated data
double getV(double x, double *xvec, double *vvec, int Nsteps) {
    return interpolate(xvec, vvec, Nsteps, x);
}

// Gravitational potential calculation (simplified for this example)
double getPhi(double r, double t, double *xvec, double *alphavec, double *vvec, int Nsteps, double c_s) {
    double rmax = (G * M) / (2 * c_s * c_s);
    double phi_max = -G * M / rmax;
    if (r >= rmax) return -G * M / r;
    else return phi_max;  // Integration over Mtot omitted for brevity
}

// Calculate phidot
double phidot(double r, double t, double *xvec, double *alphavec, double *vvec, int Nsteps, double c_s) {
    double h = t * 1e-4;
    double phi1 = getPhi(r, t + h, xvec, alphavec, vvec, Nsteps, c_s);
    double phi2 = getPhi(r, t - h, xvec, alphavec, vvec, Nsteps, c_s);
    return (phi1 - phi2) / (2 * h);
}

// Save phivec to binary file
void save_phivec_to_file(const char *filename, double *phivec, int size) {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        fprintf(stderr, "Error opening file for writing\n");
        exit(1);
    }
    fwrite(phivec, sizeof(double), size, file);
    fclose(file);
}

int main() {
    double xvec[N_STEPS], alphavec[N_STEPS], vvec[N_STEPS];
    int Nsteps = loadData("/Users/richardanderson/workdir/planetesimal_project/ShuScripts/shusolution.bin", xvec, alphavec, vvec);

    if (Nsteps == -1) {
        fprintf(stderr, "Failed to load data.\n");
        return 1;
    }

    double c_s = sqrt(SCRIPT_R * T);  // Calculate c_s within main

    double t[N_STEPS];
    double phivec[N_STEPS];
    for (int i = 0; i < N_STEPS; i++) {
        t[i] = pow(10, log10(1) + i * (13 - log10(1)) / (N_STEPS - 1));
    }

    double rmax = (G * M) / (2 * c_s * c_s);
    for (int i = 0; i < N_STEPS; i++) {
        phivec[i] = getPhi(rmax * 0.5, t[i], xvec, alphavec, vvec, Nsteps, c_s);
    }

    save_phivec_to_file("phivec.bin", phivec, N_STEPS);
    printf("phivec saved to phivec.bin\n");

    return 0;
}

