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

#ifndef INTERPOLATION_CORE_H
#define INTERPOLATION_CORE_H

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

// Read boundary extension
BoundaryExt read_ext(const char* boundary);
// Geometric transformation of an image (by an homography) using an interpolation method
void interpolate_image_homography(double *out, double *in, int w, int h, int pd, double H[9], 
                                  char *interp, BoundaryExt boundaryExt, float zoom);

#endif