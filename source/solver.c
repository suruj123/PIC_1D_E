#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

int main() {
    int N = 512; // Number of grid points
    double L = 50.0; // Length of the domain
    double dx = L / N;
    double *rho = (double*) fftw_malloc(sizeof(double) * N);
    fftw_complex *rho_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *phi_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    double *phi = (double*) fftw_malloc(sizeof(double) * N);

    int i, k;

    // Initialize charge density rho(x)
    for ( i = 0; i < N; i++) {
        double x = i * dx;
        rho[i] = sin(2 * M_PI * x / L); // Example: sin(2 * pi * x / L)
    }

    // Create FFTW plans
    fftw_plan forward_plan = fftw_plan_dft_r2c_1d(N, rho, rho_k, FFTW_ESTIMATE);
    fftw_plan backward_plan = fftw_plan_dft_c2r_1d(N, phi_k, phi, FFTW_ESTIMATE);

    // Perform forward FFT on rho
    fftw_execute(forward_plan);

    // Solve for phi in Fourier space
    phi_k[0][0] = 0.0; // k = 0 mode
    phi_k[0][1] = 0.0;
    for ( k = 1; k < N / 2 + 1; k++) {
        double k_val = 2 * M_PI * k / L;
        phi_k[k][0] = -rho_k[k][0] / (k_val * k_val); // Real part
        phi_k[k][1] = -rho_k[k][1] / (k_val * k_val); // Imaginary part
    }

    // Perform inverse FFT to get phi(x)
    fftw_execute(backward_plan);

    // Normalize the inverse FFT output
    for ( i = 0; i < N; i++) {
        phi[i] /= N;
    }

    // Output the result
    for ( i = 0; i < N; i++) {
        double x = i * dx;
        printf("%f %f\n", x, phi[i]);
    }

    // Cleanup
    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(backward_plan);
    fftw_free(rho);
    fftw_free(rho_k);
    fftw_free(phi_k);
    fftw_free(phi);

    return 0;
}

