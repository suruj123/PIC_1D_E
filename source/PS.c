#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "global.h"

void PS(){

    int i, k;


    // Create FFTW plans
    fftw_plan forward_plan = fftw_plan_dft_r2c_1d(Nx, rho, rho_k, FFTW_ESTIMATE);
    fftw_plan backward_plan = fftw_plan_dft_c2r_1d(Nx, phi_k, phi, FFTW_ESTIMATE);

    // Perform forward FFT on rho
    fftw_execute(forward_plan);

    // Solve for phi in Fourier space
    phi_k[0][0] = 0.0; // k = 0 mode
    phi_k[0][1] = 0.0;
    for ( k = 1; k < Nx / 2 + 1; k++) {
        double k_val = 2 * M_PI * k / region[1];
        phi_k[k][0] = -rho_k[k][0] / (k_val * k_val); // Real part
        phi_k[k][1] = -rho_k[k][1] / (k_val * k_val); // Imaginary part
    }

    // Perform inverse FFT to get phi(x)
    fftw_execute(backward_plan);

    // Normalize the inverse FFT output
    for ( i = 0; i < Nx; i++) {
        phi[i] /= Nx;
    }

    for ( i = 0; i < Nx; i++) {

    psi[i+1] = phi[i];

    }

    for ( i = 1; i <= Nx; i++) {

    //printf("%lf\n", psi[i]);

    }

    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(backward_plan);
}

