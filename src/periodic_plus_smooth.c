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
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include "fft_core.h"

/* M_PI is a POSIX definition */
#ifndef M_PI
/** macro definition for Pi */
#define M_PI 3.14159265358979323846
#endif                          /* !M_PI */

// Compute the jumps at the boundary of the image
static void jumps(double *out, const double *in, int w, int h, int pd)
{
    // initialization
    for (int i = 0; i < w*h*pd; i++)
        out[i] = 0;
    
    // loop over the channels
    for (int l = 0; l < pd; l++) {
        // horizontal jumps
        for (int j = 0; j < h; j++) {
                out[j*w + l*w*h] = in[j*w + l*w*h] - in[j*w + w-1 + l*w*h];  
                out[j*w + w-1 + l*w*h] -= in[j*w + l*w*h] - in[j*w + w-1 + l*w*h];
        }
        // vertical jumps    
        for (int i = 0; i < w; i++) {
            out[i + l*w*h] += in[i + l*w*h] - in[(h-1)*w + i + l*w*h];    
            out[(h-1)*w + i + l*w*h] -= in[i + l*w*h] - in[(h-1)*w + i + l*w*h]; 
        }
    }
}

// Compute the smooth component of an image using Fourier computations
static void compute_smooth_component(fftw_complex *shat, const double *in, int w, int h, int pd)
{
    // allocate memory
    double *v = (double *) malloc(w*h*pd*sizeof*v);
    
    // compute jumps
    jumps(v, in, w, h, pd);
    
    // compute the fft of jumps
    do_fft_real(shat, v, w, h, pd);

    double tmp;
    double factorh = 2*M_PI/h;
    double factorw = 2*M_PI/w;
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {
            tmp = 1.0/(4-2*cos(j*factorh)-2*cos(i*factorw));
            for (int l = 0; l < pd; l++)
                shat[j*w+i+l*w*h] = shat[j*w+i+l*w*h]*tmp;
        }
    
    // set the mean to 0
    for (int l = 0; l < pd; l++)
                shat[l*w*h] = 0.0;
    
    // free memory
    free(v);
}

// Compute the periodic component
static void compute_periodic_component(double *periodic, const double *smooth,
                                       const double *in, int w, int h, int pd)
{
    // difference of inhat and shat
    for (int i = 0; i < pd*w*h; i++)
        periodic[i] = in[i] - smooth[i];
}

// Compute the periodic plus smooth decomposition of an image
// The periodic component is possibly zoomed by TPI
void periodic_plus_smooth_decomposition(double *periodic, double *smooth, const double *in,
                                        int w, int h, int pd, int zoom)
{
    // out sizes
    int hout = zoom*h;
    int wout = zoom*w;

    // memory allocation
    fftw_complex *shat = malloc(w*h*pd*sizeof*shat);
    fftw_complex *phat = malloc(w*h*pd*sizeof*phat);
    fftw_complex *phat_zoom = malloc(wout*hout*pd*sizeof*phat_zoom);
    
    // compute smooth component
    compute_smooth_component(shat, in, w, h, pd);
    do_ifft_real(smooth, shat, w, h, pd);

    // compute the zoomed periodic component
    // 1) image - sComponent
    compute_periodic_component(periodic, smooth, in, w, h, pd);
    // 2) fft
    do_fft_real(phat, periodic, w, h, pd);
    // 3) zero-padding
    upsampling_fourier(phat_zoom, phat, w, h, wout, hout, pd, 1);
    // 4) fft inverse
    do_ifft_real(periodic, phat_zoom, wout, hout, pd);

    // free memory
    fftw_free(shat);
    fftw_free(phat);
    fftw_free(phat_zoom);
}
