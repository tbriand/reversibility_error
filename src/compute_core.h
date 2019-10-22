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

#ifndef COMPUTE_CORE
#define COMPUTE_CORE

#include <stdlib.h>
#include <math.h>

// Projection of x onto [min, max]
static int clip(int x, int min, int max)
{
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

// Cropping of an image
static void crop(float *out, int *cw, int *ch, float *in, int w, int h, int pd,
                 int x0, int y0, int xf, int yf)
{
    // if non-positive update
    if (xf <= 0)
        xf = w + xf;
    if (yf <= 0)
        yf = h + yf;
    
    // adjust bounds
    x0 = clip(x0, 0, w);
    xf = clip(xf, 0, w);
    y0 = clip(y0, 0, h);
    yf = clip(yf, 0, h);

    // output size
    *cw = xf - x0;
    *ch = yf - y0;

    // compute crop
    for (int l = 0; l < pd; l++)
        for (int j = 0; j < *ch; j++)
            for (int i = 0; i < *cw; i++)
                out[i + j*(*cw) + l*(*cw)*(*ch)] = in[i+x0 + (j+y0)*w + l*w*h];
}

// Compute the root mean square error (RMSE) of an array
static double rmse(double *in, int N) {
    
    double val = 0.0;
    for (int i=0; i<N; i++)
        val += in[i]*in[i];
    val /= (double) N;
    val = sqrt(val);
    
    return val;
}

#endif