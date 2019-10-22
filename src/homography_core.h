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

#ifndef HOMOGRAPHY_H
#define HOMOGRAPHY_H

// Compute the inverse of an homography
void invert_homography(double iH[9], const double H[9]);
// Compute the image y of the vector x by the homography H
void apply_homography(double y[2], const double x[2], const double H[9]);
// translate homography in order to compute the result of a crop
void translate_homography(double *h_out, double *h_in, double tx, double ty);
// zoom homography in order to be compatible with a zoom of an image
void zoom_homography(double *h_out, double *h_in, double zx, double zy);
// Draw a random translation
void create_random_translation(double H[9], double L);
// Draw a random homography
void create_random_homography(double H[9], int w, int h, double L);
// Draw a random homography of a given type
void create_random_transformation(double H[9], double L, int w, int h, int type);

#endif