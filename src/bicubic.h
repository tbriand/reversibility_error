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

#ifndef BICUBIC_H
#define BICUBIC_H

// Resampling of an image at locations (xpos,ypos) using bicubic interpolation 
void interpolate_bicubic(double *out, double *in, int w, int h, int pd,
                         BoundaryExt bc, double *xpos, double *ypos,
                         int numPixels);

#endif
