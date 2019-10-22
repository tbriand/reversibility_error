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

#ifndef BOUNDARY_DEFINITION
#define BOUNDARY_DEFINITION
typedef enum
{
    BOUNDARY_CONSTANT = 0,
    BOUNDARY_HSYMMETRIC = 1,
    BOUNDARY_WSYMMETRIC = 2,
    BOUNDARY_PERIODIC = 3
} BoundaryExt;
#endif

typedef double (*getsample_operator)(double*,int,int,int,int);

// Modulus with correct values
static int good_modulus(int n, int p)
{
    if (!p) return 0;
    if (p < 1) return good_modulus(n, -p);

    int r;
    if (n >= 0)
        r = n % p;
    else {
        r = p - (-n) % p;
        if (r == p)
            r = 0;
    }
    return r;
}

// Reflexion for hsym
inline static int positive_reflex(int i, int N) {
    while(1) {
        if(i < 0)
            i = -1-i;
        else if(i >= N)
            i = (2*N-1)-i;
        else
            return i;
    }
}

// Reflexion for wsym
inline static int positive_reflex2(int i, int N) {
    while(1) {
        if(i < 0)
            i = -i;
        else if(i >= N)
            i = (2*N-2)-i;
        else
            return i;
    }
}

// Extrapolate by reflection (hsym)
inline
static double getsample_hsym(double *x, int w, int h, int i, int j)
{
    i = positive_reflex(i, w);
    j = positive_reflex(j, h);
    return x[i+j*w];
}

// Extrapolate by reflection (wsym)
inline
static double getsample_wsym(double *x, int w, int h, int i, int j)
{
    i = positive_reflex2(i, w);
    j = positive_reflex2(j, h);
    return x[i+j*w];
}

// Extrapolate by periodicity
inline
static double getsample_per(double *x, int w, int h, int i, int j)
{
    i = good_modulus(i, w);
    j = good_modulus(j, h);
    return x[i+j*w];
}

// Extrapolate by constant
inline static
double getsample_constant(double *x, int w, int h, int i, int j)
{
    if (i < 0)
        i = 0;
    if (i >= w)
        i = w - 1;
    if (j < 0)
        j = 0;
    if (j >= h)
        j = h - 1;
    return x[i+j*w];
}

// Cubic interpolation
static double cubic_interpolation(double v[4], float x)
{
    return v[1] + 0.5*x*(v[2] - v[0]
        + x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3] //x^2 term
        + x*(3.0*(v[1] - v[2]) + v[3] - v[0])));   //x^3 term
}

// Bicubic interpolation
static double bicubic_interpolation_cell(double p[4][4], float x, float y)
{
    double v[4];
    v[0] = cubic_interpolation(p[0], y);
    v[1] = cubic_interpolation(p[1], y);
    v[2] = cubic_interpolation(p[2], y);
    v[3] = cubic_interpolation(p[3], y);
    return cubic_interpolation(v, x);
}

// Resampling of an image at locations (xpos,ypos) using bicubic interpolation
void interpolate_bicubic(double *out, double *in, int w, int h, int pd,
                         BoundaryExt bc, double *xpos, double *ypos,
                         int numPixels) {
    int ix, iy;
    double x, y, c[4][4];
    
    // boundary handling
    getsample_operator p;
    if ( bc == BOUNDARY_PERIODIC )
        p = getsample_per;
    if ( bc == BOUNDARY_CONSTANT )
        p = getsample_constant;
    if ( bc == BOUNDARY_HSYMMETRIC )
        p = getsample_hsym;
    if ( bc == BOUNDARY_WSYMMETRIC )
        p = getsample_wsym;

    // loop over the locations
    for (int k = 0; k < numPixels; k++) {
        x = xpos[k] - 1;
        y = ypos[k] - 1;

        ix = floor(x);
        iy = floor(y);
        for (int l = 0; l < pd; l ++) {
            for (int j = 0; j < 4; j++)
                for (int i = 0; i < 4; i++)
                    c[i][j] = p(in + l*w*h, w, h, ix + i, iy + j);
            out[k + l*numPixels] = bicubic_interpolation_cell(c, x - ix, y - iy);
        }
    }
}
