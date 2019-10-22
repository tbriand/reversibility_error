/* SPDX-License-Identifier: GPL-2.0+
 *
 * Thibaud Briand <briand.thibaud@gmail.com>
 *
 * Copyright (c) 2018-2019, Thibaud Briand
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * You should have received a copy of the GNU General Pulic License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <complex.h>
#include <fftw3.h>

#define FFTW_NTHREADS // comment to disable multithreaded FFT

// Start threaded FFTW if FFTW_NTHREADS is defined
void init_fftw(void) {
    #ifdef FFTW_NTHREADS
    fftw_init_threads();
    #ifdef _OPENMP
    fftw_plan_with_nthreads(omp_get_max_threads());
    #endif
    #endif
}

// Clean FFTW
void clean_fftw(void) {
    fftw_cleanup();
    #ifdef FFTW_NTHREADS
    fftw_cleanup_threads(); 
    #endif
}

// Compute the DFT of a real-valued image
void do_fft_real(fftw_complex *out, const double *in, int nx, int ny, int nz)
{
    // memory allocation
    fftw_complex *in_plan = (fftw_complex *) malloc(nx*ny*sizeof(fftw_complex));
    fftw_complex *out_plan = (fftw_complex *) malloc(nx*ny*sizeof(fftw_complex));
    fftw_plan plan = fftw_plan_dft_2d(ny, nx, in_plan, out_plan, FFTW_FORWARD, FFTW_ESTIMATE);

    // loop over the channels
    for (int l = 0; l < nz; l++) {
        // Real --> complex
        for(int i = 0; i < nx*ny; i++)
            in_plan[i] = (double complex) in[i + l*nx*ny];

        // compute fft
        fftw_execute(plan);

        // copy to output
        memcpy(out + l*nx*ny, out_plan, nx*ny*sizeof(fftw_complex));
    }

    // free
    fftw_destroy_plan(plan);
    fftw_free(in_plan);
    fftw_free(out_plan);
}

// Compute the real part of the iDFT of a complex-valued image
void do_ifft_real(double *out, const fftw_complex *in, int nx, int ny, int nz)
{
    // memory allocation
    fftw_complex *in_plan = (fftw_complex *) malloc(nx*ny*sizeof(fftw_complex));
    fftw_complex *out_plan = (fftw_complex *) malloc(nx*ny*sizeof(fftw_complex));
    fftw_plan plan = fftw_plan_dft_2d (ny, nx, in_plan, out_plan, FFTW_BACKWARD, FFTW_ESTIMATE);

    // normalization constant
    double norm = 1.0/(nx*ny);

    // loop over the channels
    for (int l = 0; l < nz; l++) {
        // copy to input
        memcpy(in_plan, in + l*nx*ny, nx*ny*sizeof(fftw_complex));

        // compute ifft
        fftw_execute(plan);

        // complex to real + normalization
        for(int i = 0; i < nx*ny; i++)
            out[i + l*nx*ny] = creal(out_plan[i])*norm;
    }

    // free
    fftw_destroy_plan(plan);
    fftw_free(in_plan);
    fftw_free(out_plan);
}

// Compute the fftshift of a complex-valued image
void fftshift(fftw_complex *fshift, fftw_complex *fhat, int nx, int ny, int nz) {
    int nx2 = nx/2;
    int ny2 = ny/2;
    int cx = (nx+1)/2;
    int cy = (ny+1)/2;
    int i, j, l, i2, j2;

    for(j = 0; j < ny; j++) {
        j2 = (j < cy) ? j + ny2 : j - cy;
        for(i = 0; i < nx; i++) {
            i2 = (i < cx) ? i + nx2 : i - cx;
            for(l = 0; l < nz; l++)
                fshift[i2 + j2*nx + l*nx*ny] = fhat[i + j*nx + l*nx*ny];
        }
    }
}

// Compute the DFT coefficients of the up-sampled image
// See https://www.ipol.im/pub/art/2019/273/ (Line 3 of Algorithm 3 using Proposition 11)
void upsampling_fourier(fftw_complex *out, fftw_complex *in,
                               int nxin, int nyin, int nxout, int nyout, int nz, int interp)
{
    int i, j, l, i2, j2;
    
    // normalization constant
    double norm = nxout*nyout*1.0/(nxin*nyin);
    
    // indices for the fftshift
    int nx2 = (nxin+1)/2; 
    int ny2 = (nyin+1)/2;

    // fill the output dft with zeros
    for (i = 0; i < nxout*nyout*nz; i++)
        out[i] = 0.0;

    // fill the corners with the values
    for(j = 0; j < nyin; j++) {
        j2 = (j < ny2) ? j : j + nyout-nyin;
        for(i = 0; i < nxin; i++) {
            i2 = (i < nx2) ? i : i + nxout-nxin;
            for (l = 0; l < nz; l++)
                out[i2 + j2*nxout + l*nxout*nyout] = norm*in[i + j*nxin + l*nxin*nyin];
        }
    }

    // real part
        // useless in practice if the real part of the image is taken afterwards
        if ( !(nxin%2) && nxout>nxin) {
            i = nx2; // positive in output and negative in input
            i2 = nx2 + nxout-nxin; // negative in output (already initialized)
            for(j = 0; j < nyin; j++) {
                j2 = (j < ny2) ? j : j + nyout-nyin;
                for (l = 0; l < nz; l++) {
                    out[i2 + j2*nxout + l*nxout*nyout] *= 0.5;
                    out[i + j2*nxout + l*nxout*nyout] = out[i2 + j2*nxout + l*nxout*nyout];
                }
            }
        }

        if ( !(nyin%2) && nyout>nyin) {
            j = ny2; // positive in output and negative in input
            j2 = ny2 + nyout-nyin; // negative in output (already initialized)
            for(i = 0; i < nxin; i++) {
                i2 = (i < nx2) ? i : i + nxout-nxin;
                for (l = 0; l < nz; l++) {
                    out[i2 + j2*nxout + l*nxout*nyout] *= 0.5;
                    out[i2 + j*nxout + l*nxout*nyout] = out[i2 + j2*nxout + l*nxout*nyout];
                }
            }
        }
        
        if ( !interp && !(nxin%2) && !(nyin%2) && nxout>nxin && nyout>nyin) {
            i = nx2; // positive in output and negative in input
            i2 = nx2 + nxout-nxin; // negative in output
            j = ny2; // positive in output and negative in input
            j2 = ny2 + nyout-nyin; // negative in output
            for (l = 0; l < nz; l++) {
                double complex hf = norm*0.5*in[i + j*nxin + l*nxin*nyin];
                out[i + j*nxout + l*nxout*nyout] = out[i2 + j2*nxout + l*nxout*nyout] = hf;
            }
        }
        
    // real convention
    if ( interp && !(nxin%2) && !(nyin%2) && nxout>nxin && nyout>nyin) {
        i = nx2; // positive in output and negative in input
        i2 = nx2 + nxout-nxin; // negative in output
        j = ny2; // positive in output and negative in input
        j2 = ny2 + nyout-nyin; // negative in output
        for (l = 0; l < nz; l++) {
            double complex hf = norm*0.25*in[i + j*nxin + l*nxin*nyin];
            out[i + j*nxout + l*nxout*nyout] = out[i2 + j*nxout + l*nxout*nyout] = out[i + j2*nxout + l*nxout*nyout] = out[i2 + j2*nxout + l*nxout*nyout] = hf;
        }
    }
}

// Up-sampling of an image using TPI
// See https://www.ipol.im/pub/art/2019/273/ (Algorithm 3)
void upsampling(double *out, double *in, int nxin, int nyin, int nxout, int nyout, int nz, int interp) 
{
    // allocate memory for fourier transform
    fftw_complex *inhat = fftw_malloc(nxin*nyin*nz*sizeof*inhat);
    fftw_complex *outhat = fftw_malloc(nxout*nyout*nz*sizeof*outhat);

    // compute DFT of the input
    do_fft_real(inhat, in, nxin, nyin, nz);

    // phase shift (complex convention)
    upsampling_fourier(outhat, inhat, nxin, nyin, nxout, nyout, nz, interp);

    // compute iDFT of the output
    do_ifft_real(out, outhat, nxout, nyout, nz);

    // free memory
    fftw_free(inhat);
    fftw_free(outhat);
}

// Spectrum clipping of an image in the Fourier domain (Equation 8)
static void spectrum_clipping_fourier(fftw_complex *outhat, fftw_complex *inhat, int nx, int ny, int nz, double r)
{
    // indices for the fftshift
    int nx2 = (nx+1)/2; 
    int ny2 = (ny+1)/2;
    
    int i2, j2, factori, factorj;
    
    for(int j = 0; j < ny; j++) {
        j2 = (j < ny2) ? j : j - ny;
        factorj = ( 2*fabs(j2) > (1 - r)*ny ) ? 0 : 1;
        for(int i = 0; i < nx; i++) {
            i2 = (i < nx2) ? i : i - nx;
            factori = ( 2*fabs(i2) > (1 - r)*nx ) ? 0 : 1;
            for(int l = 0; l < nz; l++)
                outhat[i+j*nx+l*nx*ny] = factori*factorj*inhat[i+j*nx+l*nx*ny];
        }
    }
}

// Spectrum clipping of an image (Definition 6)
void spectrum_clipping(double *out, double *in, int nx, int ny, int nz, double r)
{
    // allocate memory for fourier transform
    fftw_complex *inhat = fftw_malloc(nx*ny*nz*sizeof*inhat);
    
    // compute DFT of the input
    do_fft_real(inhat, in, nx, ny, nz);
    
    // phase shift (complex convention)
    spectrum_clipping_fourier(inhat, inhat, nx, ny, nz, r);
    
    // compute iDFT of the output
    do_ifft_real(out, inhat, nx, ny, nz);
    
    // free memory
    fftw_free(inhat);
}
