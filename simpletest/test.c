#include <stdio.h>
#include <math.h>
#include <rebound.h>

const double G_SI = 6.67430e-11;
const double SUN_MASS = 1.989e30;
const double EARTH_MASS = 5.972e24;
const double AU = 1.496e11;
const double SECONDS_PER_YEAR = 3.154e7;
const double TF = 350000; // in years

struct reb_simulation* sim2;
void custom_sun_force(struct reb_simulation* r) {
    // r is the raw pointer passed from REBOUND, but we use the global sim2 to match your Python code style
    for (int i = 0; i < sim2->N; i++) {
        struct reb_particle* p = &sim2->particles[i];
        double dx = -p->x;
        double dy = -p->y;
        double dz = -p->z;
        double r3 = pow(dx*dx + dy*dy + dz*dz, 1.5);
        if (r3 == 0) continue;
        double a = G_SI * SUN_MASS / r3;
        p->ax += a * dx;
        p->ay += a * dy;
        p->az += a * dz;
    }
}

int main() {
    double tf_seconds = TF * SECONDS_PER_YEAR;

    // ---------------------------
    // Simulation 1: Real Sun
    struct reb_simulation* sim = reb_create_simulation();
    sim->G = G_SI;

    reb_add(sim, (struct reb_particle){.m=SUN_MASS});
    double v_circ = sqrt(G_SI * SUN_MASS / AU);
    reb_add(sim, (struct reb_particle){.m=EARTH_MASS, .x=AU, .vy=v_circ});

    reb_integrate(sim, tf_seconds);

    struct reb_particle earth1 = sim->particles[1];

    // ---------------------------
    // Simulation 2: Custom Force
    sim2 = reb_create_simulation();
    sim2->G = G_SI;
    reb_add(sim2, (struct reb_particle){.m=EARTH_MASS, .x=AU, .vy=v_circ});
    sim2->additional_forces = custom_sun_force;

    reb_integrate(sim2, tf_seconds);

    struct reb_particle earth2 = sim2->particles[0];

    // ---------------------------
    // Print final positions
    printf("Final positions:\n");
    printf("  Real:   (%.3e, %.3e)\n", earth1.x, earth1.y);
    printf("  Custom: (%.3e, %.3e)\n", earth2.x, earth2.y);

    reb_free_simulation(sim);
    reb_free_simulation(sim2);
    return 0;
}

