#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "rebound.h"

// constants
const double script_R = 3.36e3;               // m^2/s^2/K
const double T = 10.0;                        // Kelvin
const double B = 8.86;
const double M = 1.989e30;                    // kg
const double G = 6.674e-11;                   // m^3 kg^-1 s^-2
const double a = 0.975471932310942752606143078377;
const double c_s = 1.8330302780e2; 
double global_t = 0.0;
double global_val2 = 0.0;
double tmax = 40.;
#define N_POINTS 10000
#define DATA_FILE "/Users/richardanderson/workdir/planetesimal_project/clusterscripts/silsbeescripts/shuInt.txt"
#define AU_IN_METERS 1.49e11
#define SECONDS_IN_YEAR 31557600.0

gsl_interp_accel *acc;
gsl_spline *spline;
double spline_min_x = 0.0;
double spline_max_x = 0.0;

void additional_forces(struct reb_simulation* const r);
void heartbeat(struct reb_simulation* const r);

double au_to_meters(double au) {
    return au * AU_IN_METERS;
}

double meters_to_au(double meters) {
    return meters / AU_IN_METERS;
}

double years_to_seconds(double years) {
    return years * SECONDS_IN_YEAR;
}

double seconds_to_years(double seconds) {
    return seconds / SECONDS_IN_YEAR;
}

double meters3_to_au3(double cubic_meters) {
    double au_cubed = pow(AU_IN_METERS, 3);
    return cubic_meters / au_cubed;
}

typedef double (*func_ptr)(double);

int load_data(const char *filename, double *y_arr, int n_points) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return 0;
    }

    for (int i = 0; i < n_points; i++) {
        if (fscanf(file, "%lf", &y_arr[i]) != 1) {
            fprintf(stderr, "Error reading data at line %d\n", i+1);
            fclose(file);
            return 0;
        }
    }
    fclose(file);
    return 1;
}

double trapIntegrateLog(func_ptr f, double xmin, double xmax, int N) {
    double *s = malloc(N * sizeof(double));
    double *fvec = malloc(N * sizeof(double));
    double logxmin = log10(xmin);
    double logxmax = log10(xmax);
    double m = 0.0;

    for (int i = 0; i < N; i++) {
        double exponent = logxmin + (logxmax - logxmin) * i / (N - 1);
        s[i] = pow(10.0, exponent);
        fvec[i] = f(s[i]);
    }

    for (int i = 0; i < N - 1; i++) {
        double deltax = s[i + 1] - s[i];
        double av = (fvec[i + 1] + fvec[i]) / 2.0;
        m += deltax * av;
    }

    free(s);
    free(fvec);
    return m;
}

double trapIntegrateLinear(func_ptr f, double xmin, double xmax, int N) {
    double *s = malloc(N * sizeof(double));
    double *fvec = malloc(N * sizeof(double));
    double m = 0.0;

    for (int i = 0; i < N; i++) {
        s[i] = xmin + (xmax - xmin) * i / (N - 1);
        fvec[i] = f(s[i]);
    }

    for (int i = 0; i < N - 1; i++) {
        double deltax = s[i + 1] - s[i];
        double av = (fvec[i + 1] + fvec[i]) / 2.0;
        m += deltax * av;
    }

    free(s);
    free(fvec);
    return m;
}

void logspace(double *array, double log_min, double log_max, int num_points) {
    double delta = (log_max - log_min) / (num_points - 1);
    for (int i = 0; i < num_points; i++) {
        array[i] = pow(10.0, log_min + delta * i);
    }
}

double get_x_integral(double x, double val2) {
    static int printed_spline_range = 0;
    /*
    if (!printed_spline_range && spline != NULL) {
        printf("Spline valid range: [%e, %e]\n", spline_min_x, spline_max_x);
        fflush(stdout);
        printed_spline_range = 1;
    }
    */

    if (spline == NULL || acc == NULL) {
        printf("Error: Spline or accelerator is NULL!\n");
        fflush(stdout);
        return 0.0;
    }

    if (!isfinite(x)) {
        printf("Warning: x is not finite: %e\n", x);
        fflush(stdout);
        return 0.0;
    }

    if (x < spline_min_x) {
        printf("Warning: x = %e < spline min = %e. Returning 0.\n", x, spline_min_x);
        fflush(stdout);
        return 0.0;
    }

    if (x > spline_max_x) {
        // printf("Warning: x = %e > spline max = %e. Using extrapolation.\n", x, spline_max_x);
        // fflush(stdout);
        return val2 + 2.0 * (x - 2.0);
    }

    double result = gsl_spline_eval(spline, x, acc);

    if (!isfinite(result)) {
        printf("Warning: gsl_spline_eval returned non-finite for x = %e\n", x);
        fflush(stdout);
        return 0.0;
    }

    return result;
}

double Mtot(double r, double t) {
    double r0 = c_s * t;
    double centralMass = 0.975502 * c_s * c_s * r0 / G;
    double Mcalc = centralMass + pow(r0, 3) / (G * t * t) * get_x_integral(r / r0, global_val2);
    return (Mcalc < M) ? Mcalc : M;
}

double Mtot_wrapper(double rp) {
    double r0 = c_s * global_t;
    double centralMass = 0.975502 * c_s * c_s * r0 / G;
    double Mcalc = centralMass + pow(r0, 3) / (G * global_t * global_t) * get_x_integral(rp / r0, global_val2);
    return Mcalc < M ? Mcalc / (rp * rp) : M / (rp * rp);
}

double trapIntegrateLog_wrapper(double r, double t, double rmax) {
    global_t = t;
    return trapIntegrateLog(Mtot_wrapper, r, rmax, 10000);
}

double getPhi(double r, double t, double val2) {
    double rmax = (G * M) / (2.0 * c_s * c_s);
    double phi_max = -G * M / rmax;

    if (r >= rmax) {
        return -G * M / r;
    } else {
        global_t = t;
        global_val2 = val2;

        double integral = trapIntegrateLog(Mtot_wrapper, r, rmax, 10000);
        return phi_max - G * integral;
    }
}

double get_bmax(double R, double v0, double t) {
    double phi = getPhi(R, t, global_val2);
    return R * sqrt(1.0 - 2.0 * phi / (v0 * v0));
}

void get_params(double b, double R, double v0, double t, double *vr, double *vazimuthal) {
    double potential = getPhi(R, t, global_val2);
    double v = sqrt(v0 * v0 - 2.0 * potential);
    *vazimuthal = v0 * b / R;
    *vr = sqrt(v * v - (*vazimuthal) * (*vazimuthal));
}

typedef struct {
    double x, y, vx, vy;
} Particle;

double getEnergy(Particle p, double t) {
    double r = sqrt(p.x * p.x + p.y * p.y);
    double v = sqrt(p.vx * p.vx + p.vy * p.vy);
    return 0.5 * v * v + getPhi(r, t, global_val2);
}

void simulate_particle(double v0, double tf) {
    struct reb_simulation* r = reb_simulation_create();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, time(NULL));
    r->ri_ias15.min_dt = 1e-6;
    // r->dt             = 1000.; 
    r->integrator     = REB_INTEGRATOR_IAS15;
    r->gravity        = REB_GRAVITY_NONE;
    r->G = G;
    r->additional_forces = additional_forces;
    r->force_is_velocity_dependent = 0;

    r->t = 0.0;

    int Nsteps = 100;
    int pid = 0;
    double* tvec = malloc((Nsteps + 1) * sizeof(double));
    double tmin = tf / 1e9;
    double tmax_local = tf;
    logspace(tvec, log10(tmin), log10(tmax_local), Nsteps + 1);
    /*
    char pos_filename[100];
    snprintf(pos_filename, sizeof(pos_filename), "positions_%.2f.txt", v0);
    FILE *fout_pos = fopen(pos_filename, "w");
    
    if (!fout_pos) {
        fprintf(stderr, "Error: Could not open %s for writing.\n", pos_filename);
    }
    */
    for (int i = 0; i < Nsteps; i++) {
        double t_current = tvec[i];
        double r0 = c_s * t_current;
        double bmiss = get_bmax(r0, v0, t_current);
        double rand_val = ((double)rand()) / RAND_MAX;
        double b = sqrt(bmiss * bmiss * rand_val);

        double N_exp = 0.0;
        if (bmiss > b) {
            double delta_t = tvec[i+1] - tvec[i];
            N_exp = (10.0 / tf) * delta_t;
        }

        int N_static = gsl_ran_poisson(rng, N_exp); 

        for (int j = 0; j < N_static; j++) {
            double rand_val2 = ((double)rand()) / RAND_MAX;
            double b2 = sqrt(bmiss * bmiss * rand_val2);
            if (b2 < bmiss) {
                double vr, vtheta;
                get_params(b2, r0, v0, t_current, &vr, &vtheta);

                struct reb_particle p = {0};
                p.m  = 0.0;
                p.x  = -r0;       
                p.y  = 0;
                p.vx = vr;
                p.vy = vtheta;
                reb_simulation_add(r, p);
                pid+=1;
            }
        }

        if (r->N > 0) {
            reb_simulation_integrate(r, tvec[i+1]);
        } else {
            r->t = tvec[i+1];
        }

        // Write positions of all particles at current time step
        /*
        if (fout_pos) {
            const struct reb_particle* const particles = r->particles;
            for (int k = 0; k < r->N; k++) {
                fprintf(fout_pos, "%.6e %.6e %.6e %d\n", r->t, particles[k].x, particles[k].y, pid);
            }
        }
        */
    }

    gsl_rng_free(rng);
    FILE *fout = fopen("energies.txt", "a");
    if (!fout) {
        fprintf(stderr, "Error: Could not open energies.txt for writing.\n");
    }
    const struct reb_particle* const particles = r->particles;

    for (int i = 0; i < r->N; i++) {
        double rpos = sqrt(particles[i].x * particles[i].x + particles[i].y * particles[i].y);
        double vmag = sqrt(particles[i].vx * particles[i].vx + particles[i].vy * particles[i].vy);

        Particle p = { particles[i].x, particles[i].y, particles[i].vx, particles[i].vy };
        double energy = getEnergy(p, r->t);

        printf("Particle %d: r = %.3e m, v = %.3e m/s, E = %.6e J\n", i, rpos, vmag, energy);

        if (fout) {
            // fprintf(fout, "%.6e\n",energy);
            fprintf(fout, "%.6e %.6e\n", v0, energy);
        }
    }

    if (fout) fclose(fout);

    free(tvec);
    reb_simulation_free(r);
}

void additional_forces(struct reb_simulation* const sim) {
    const double t = sim->t;
    const double G_local = sim->G;

    struct reb_particle* const particles = sim->particles;
    const int N = sim->N;

    for (int i = 0; i < N; i++) { 
        double x = particles[i].x;
        double y = particles[i].y;
        double r = sqrt(x*x + y*y);

        if (r > 1e4) {
            double M_enc = Mtot(r, t);  
            double a_mag = G_local * M_enc / (r * r);

            particles[i].ax = -x / r * a_mag;
            particles[i].ay = -y / r * a_mag;
        }
    }
}

int main(int argc, char *argv[]) {
    double c_s = sqrt(script_R * T);
    if (argc < 3) {
        printf("Usage: %s v0 t_f\n", argv[0]);
        return 1;
    }

    double* v_x = malloc(N_POINTS * sizeof(double));
    double* v_integral = malloc(N_POINTS * sizeof(double));
    if (!v_x || !v_integral) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 1;
    }

    logspace(v_x, log10(1e-12), log10(2.0), N_POINTS);

    if (!load_data(DATA_FILE, v_integral, N_POINTS)) {
        free(v_x);
        free(v_integral);
        return 1;
    }

    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, N_POINTS);
    gsl_spline_init(spline, v_x, v_integral, N_POINTS);

    spline_min_x = v_x[0];
    spline_max_x = v_x[N_POINTS - 1];

    double v0 = atof(argv[1]);
    double t_f = atof(argv[2])*SECONDS_IN_YEAR;

    printf("Running Simulation");

    simulate_particle(v0, t_f);

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    free(v_x);
    free(v_integral);

    return 0;
}

