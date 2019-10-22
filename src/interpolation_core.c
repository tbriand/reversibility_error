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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "interpolation_core.h"
#include "homography_core.h"
#include "periodic_plus_smooth.h"
#include "splinter.h"
#include "bicubic.h"
#include "tpi.h"
#include "fft_core.h"

// Read boundary extension
BoundaryExt read_ext(const char* boundary) {
    if(0 == strncmp(boundary, "constant", strlen(boundary)))
        return BOUNDARY_CONSTANT;
    if(0 == strncmp(boundary, "periodic", strlen(boundary)))
        return BOUNDARY_PERIODIC;
    if(0 == strncmp(boundary, "hsymmetric", strlen(boundary)))
        return BOUNDARY_HSYMMETRIC;
    if(0 == strncmp(boundary, "wsymmetric", strlen(boundary)))
        return BOUNDARY_WSYMMETRIC;
    fprintf(stderr,"Unknown boundary condition %s\n",boundary);
    exit(EXIT_FAILURE);
}

// Compare end of string (useful for reading the interpolation method)
static int EndsWith(const char *str, const char *suffix)
{
    if (!str || !suffix)
        return 0;
    size_t lenstr = strlen(str);
    size_t lensuffix = strlen(suffix);
    if (lensuffix >  lenstr)
        return 0;
    return strncmp(str + lenstr - lensuffix, suffix, lensuffix) == 0;
}

// Resampling of an image at given locations (x,y) using B-spline interpolation
static void splinter_at(double *out, double *in, int w, int h, int pd, int n,
                        BoundaryExt bc, double precision, int larger, double *x,
                        double *y, int numPixels) {
    // init plan (prefiltering)
    splinter_plan_t plan = splinter_plan(in, w, h, pd, n, bc, precision, larger);
    
    // computation of the pixel locations
    double *outp = malloc(pd*sizeof*outp);
    for(int i = 0; i < numPixels; i++) {
            splinter(outp, x[i], y[i], plan);
            for(int k = 0; k < pd; k++)
                out[k*numPixels] = outp[k];
            ++out;
    }
    
    free(outp);
    splinter_destroy_plan(plan);
}

// Resampling of an image at given locations (x,y)
// using a base interpolation method
static void interpolate_at(double *out, double *in, int w, int h, int pd,
                           char *interp, BoundaryExt bc, double *x,
                           double *y, int numPixels) {
    if (0 == strncmp(interp, "bic", 3))
        interpolate_bicubic(out, in, w, h, pd, bc, x, y, numPixels);
    else if (0 == strncmp(interp, "tpi", 3))
        interpolate_at_locations_nfft(out, in, w, h, pd, x, y, numPixels, 1);
    else if (0 == strncmp(interp, "spline", 6)) {
        // order
        int order = -1;
        sscanf(interp, "spline%d", &order);
        
        if ( order < 0 ) {
            order = 0;
            printf("Negative order in B-spline interpolation...");
            printf(" Switching to order 0\n");
        }
        else if ( order > 16 ) {
            order = 16;
            printf("Maximal order is 16...\n");
        }
        
        // precision
        double precision = 1e-12;
        
        // larger algorithm or not
        int larger = 0;
        if ( bc == BOUNDARY_CONSTANT )
            larger = 1;
        
        // interpolate at locations using B-spline interpolation
        splinter_at(out, in, w, h, pd, order,
                    bc, precision, larger, x, y, numPixels);
    }
    else 
        printf("Unknown interpolation method...\n");
}

// Resampling of an image at given locations (x,y)
// using an interpolation method (base, zoomed or p+s)
// For the zoomed version this corresponds to Algorithm 3
// For the p+s version this corresponds to Algorithm 4
static void interpolate_image_at_method(double *out, double *in, int w, int h,
                                        int pd, char *interp, BoundaryExt bc,
                                        double *x, double *y, int numPixels) {
    if (0 == strncmp(interp, "p+s", 3)) {
        // periodic plus smooth version (Algorithm 4)
        int zoom = 2;
        int wper = w*zoom;
        int hper = h*zoom;
        
        // periodic plus smooth decomposition
        double *pComp_in_zoomed = malloc(wper*hper*pd*sizeof(double));
        double *sComp_in = malloc(w*h*pd*sizeof(double));
        periodic_plus_smooth_decomposition(pComp_in_zoomed, sComp_in,
                                           in, w, h, pd, zoom);
        
        // extract interpolation method for each component
        char *interp_perio  = strchr(interp, '-') + 1;
        char *interp_smooth = strrchr(interp, '-') + 1;
        
        // interpolate smooth
        interpolate_at(out, sComp_in, w, h, pd, interp_smooth, bc,
                       x, y, numPixels);
        
        // create pixel locations for the zoomed periodic version
        for (int i = 0; i < numPixels; i++) {
                x[i] *= zoom;
                y[i] *= zoom;
        }
        
        // interpolate periodic component
        double *pComp = malloc(numPixels*pd*sizeof(double));
        interpolate_at(pComp, pComp_in_zoomed, wper, hper, pd,
                       interp_perio, BOUNDARY_PERIODIC, x, y, numPixels);
        
        // sum
        for (int k = 0; k < numPixels*pd; k++)
            out[k] += pComp[k];
        
        // free memory
        free(pComp_in_zoomed);
        free(sComp_in);
        free(pComp);
    }
    else if ( EndsWith(interp,"-z2") ) { // zoomed version (Algorithm 3)
        int zoom = 2;
        int w2 = w*zoom;
        int h2 = h*zoom;

        // create pixel locations for the zoomed version
        for (int i = 0; i < numPixels; i++) {
                x[i] *= zoom;
                y[i] *= zoom;
        }
        
        // up-sample the input image
        double *in_zoomed = malloc(w2*h2*pd*sizeof(double));
        upsampling(in_zoomed, in, w, h, w2, h2, pd, 1);
        
        // interpolation at locations
        interpolate_at(out, in_zoomed, w2, h2, pd, interp, bc, x, y, numPixels);
        
        // free memory
        free(in_zoomed);
    }
    else
        interpolate_at(out, in, w, h, pd, interp, bc, x, y, numPixels);
}

// Geometric transformation of an image (by an homography)
// using an interpolation method
void interpolate_image_homography(double *out, double *in, int w, int h, int pd,
                                  double H[9], char *interp, BoundaryExt bc,
                                  float zoom) {
    // output sizes
    int wout = w/zoom;
    int hout = h/zoom;
    int numPixels = wout*hout;
    
    // create pixel locations
    double iH[9];
    invert_homography(iH, H);
    double *x = malloc(numPixels*sizeof*x);
    double *y = malloc(numPixels*sizeof*y);
    double p[2], q[2];
    for (int j = 0; j < hout; j++) {
        p[1] = j*zoom;
        for (int i = 0; i < wout; i++) {
            p[0] = i*zoom;
            apply_homography(q, p, iH);
            x[j*wout+i] = q[0];
            y[j*wout+i] = q[1];
        }
    }
    
    // interpolation at the locations using the interpolation method
    interpolate_image_at_method(out, in, w, h, pd, interp, bc, x, y, numPixels);
    
    // free memory
    free(x);
    free(y);
}
