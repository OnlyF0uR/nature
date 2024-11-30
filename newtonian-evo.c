#include <stdio.h>
#include <math.h>

#define C 299792458  // Speed of light in m/s
#define G 6.67430e-11  // Gravitational constant in m^3 kg^-1 s^-2

// Schwarzschild metric
void schwarzschild_metric(double r, double M, double *g_tt, double *g_rr, double *g_theta_theta, double *g_phi_phi) {
    *g_tt = -(1 - 2 * G * M / (r * C * C));
    *g_rr = 1 / (1 - 2 * G * M / (r * C * C));
    *g_theta_theta = r * r;
    *g_phi_phi = r * r * sin(M_PI / 2) * sin(M_PI / 2); // Assume theta = pi/2 for simplicity
}

// Christoffel symbols (simplified)
void christoffel_symbols(double r, double M, double *Gamma_ttr, double *Gamma_rtt, double *Gamma_rrr, double *Gamma_rtheta_theta,
                          double *Gamma_rphi_phi, double *Gamma_theta_r_theta, double *Gamma_phi_r_phi, double *Gamma_phi_theta_phi) {
    *Gamma_ttr = G * M / (r * r * (1 - 2 * G * M / (r * C * C)));
    *Gamma_rtt = (G * M * C * C / r / r) * (1 - 2 * G * M / (r * C * C));
    *Gamma_rrr = -G * M / (r * r * (1 - 2 * G * M / (r * C * C)));
    *Gamma_rtheta_theta = -r * (1 - 2 * G * M / (r * C * C));
    *Gamma_rphi_phi = -r * sin(M_PI / 2) * sin(M_PI / 2) * (1 - 2 * G * M / (r * C * C));
    *Gamma_theta_r_theta = 1 / r;
    *Gamma_phi_r_phi = 1 / r;
    *Gamma_phi_theta_phi = 1 / tan(M_PI / 2); // Assume theta = pi/2 for simplicity
}

// Simple time evolution using Euler's method
void simple_evolution(double initial_r, double initial_v, double M, double dt, int steps) {
    double r = initial_r;
    double v = initial_v;
    double time = 0;

    printf("Time(s)\tRadius(m)\tVelocity(m/s)\tVelocity Delta(m/s)");
    printf("-------------------------------------------------------------------------------\n");

    double prev_v = v;

    for (int i = 0; i < steps; i++) {
        // Calculate acceleration using simplified geodesic equation
        double a = -G * M / (r * r);

        // Update position and velocity
        r += v * dt;
        v += a * dt;
        time += dt;

        double v_delta = v - prev_v;

        printf("%.2f\t%.5f\t%.5f\t%.5f\n", time, r, v, v_delta);

        // Check for collision with the center (r <= 0)
        if (r <= 0) {
            printf("The object 'collides' with the center after: %.2f seconds\n", time);
            return;
        }

        prev_v = v;
    }

    printf("The object never 'collides' with the center within the simulated time.\n");
}

int main() {
    // Variables for Schwarzschild metric and Christoffel symbols
    double r = 1e8;  // Initial radius in meters
    double M = 1.989e30;  // Mass of the Sun in kg
    double g_tt, g_rr, g_theta_theta, g_phi_phi;

    // Schwarzschild metric calculation
    schwarzschild_metric(r, M, &g_tt, &g_rr, &g_theta_theta, &g_phi_phi);
    printf("Schwarzschild Metric at r = %.2f m:\n", r);
    printf("g_tt = %.5f, g_rr = %.5f, g_theta_theta = %.5f, g_phi_phi = %.5f\n\n", g_tt, g_rr, g_theta_theta, g_phi_phi);

    // Christoffel symbols calculation
    double Gamma_ttr, Gamma_rtt, Gamma_rrr, Gamma_rtheta_theta, Gamma_rphi_phi, Gamma_theta_r_theta, Gamma_phi_r_phi, Gamma_phi_theta_phi;
    christoffel_symbols(r, M, &Gamma_ttr, &Gamma_rtt, &Gamma_rrr, &Gamma_rtheta_theta, &Gamma_rphi_phi, &Gamma_theta_r_theta, &Gamma_phi_r_phi, &Gamma_phi_theta_phi);
    printf("Christoffel Symbols at r = %.2f m:\n", r);
    printf("Gamma_ttr = %.5f, Gamma_rtt = %.5f, Gamma_rrr = %.5f, Gamma_rtheta_theta = %.5f, Gamma_rphi_phi = %.5f\n", 
           Gamma_ttr, Gamma_rtt, Gamma_rrr, Gamma_rtheta_theta, Gamma_rphi_phi);
    printf("Gamma_theta_r_theta = %.5f, Gamma_phi_r_phi = %.5f, Gamma_phi_theta_phi = %.5f\n\n", 
           Gamma_theta_r_theta, Gamma_phi_r_phi, Gamma_phi_theta_phi);

    // Basic time evolution
    double initial_v = 1000;  // Initial velocity in m/s
    double dt = 0.1;  // Time step in seconds
    int steps = 1000;  // Number of time steps

    simple_evolution(r, initial_v, M, dt, steps);

    return 0;
}
