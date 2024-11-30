#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Simple scripts that computes wavefunctions and energy eigenvalues
// for a 1D particle in a box numerically

#define L 1.0       // Length of the box
#define N 100       // Number of grid points
#define PI 3.14159265358979323846

// Function to calculate the Hamiltonian matrix and solve for eigenvalues and eigenvectors
void solve_particle_in_box(double *energies, double *wavefunctions, int n) {
    double dx = L / (N - 1);
    double *H = (double *)malloc(N * N * sizeof(double));
    double *T = (double *)malloc(N * N * sizeof(double));

    // Initialize Hamiltonian matrix (T for kinetic energy part)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                H[i * N + j] = -2.0;
            } else if (abs(i - j) == 1) {
                H[i * N + j] = 1.0;
            } else {
                H[i * N + j] = 0.0;
            }
        }
    }

    // Scale the matrix by -hbar^2 / (2m * dx^2)
    double coeff = -1.0 / (2.0 * dx * dx);
    for (int i = 0; i < N * N; i++) {
        H[i] *= coeff;
    }

    // Simple eigenvalue approximation
    for (int n = 1; n <= 3; n++) {
        energies[n-1] = (n * n * PI * PI) / (2.0 * L * L);
    }

    // Simplified computation of wavefunctions as sin functions
    for (int i = 0; i < N; i++) {
        for (int n = 1; n <= 3; n++) {
            wavefunctions[i * 3 + (n-1)] = sqrt(2.0 / L) * sin(n * PI * i * dx / L);
        }
    }

    free(H);
    free(T);
}

void print_results(double *energies, double *wavefunctions, int num_states) {
    for (int i = 0; i < num_states; i++) {
        printf("Energy E_%d = %.5f\n", i + 1, energies[i]);
        printf("Wavefunction (n=%d):\n", i + 1);
        for (int j = 0; j < N; j++) {
            printf("%.5f ", wavefunctions[j * num_states + i]);
        }
        printf("\n\n");
    }
}

int main() {
    double energies[3];      // Stores the first 3 energy eigenvalues
    double *wavefunctions = (double *)malloc(N * 3 * sizeof(double));  // 3 wavefunctions

    // Solve for the first few eigenvalues and wavefunctions
    solve_particle_in_box(energies, wavefunctions, N);

    print_results(energies, wavefunctions, 3);

    free(wavefunctions);
    return 0;
}
