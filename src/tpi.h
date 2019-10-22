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

#ifndef TPI_H
#define TPI_H

// Transformation of an image using trigonometric polynomial interpolation
void interpolate_at_locations_nfft(double *out, const double *in, int nx, int ny, int nz,
                                   double *x, double *y, int numPixels, int interp);

#endif
