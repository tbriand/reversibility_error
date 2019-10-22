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
#include "iio.h"
#include "compute_core.h"

// Main function for computing the cropped version of an image
int main(int c, char *v[])
{
    // read parameters
    if (c != 7) {
            fprintf(stderr, "usage:\n\t%s x0 y0 xf yf in out\n", *v);
            //                          0 1  2  3  4  5  6
            return EXIT_FAILURE;
    }
    int x0 = atoi(v[1]);
    int y0 = atoi(v[2]);
    int xf = atoi(v[3]);
    int yf = atoi(v[4]);
    char *filename_in = v[5];
    char *filename_out = v[6];

    // read image
    int w, h, pd;
    float *image_in = iio_read_image_float_split(filename_in, &w, &h, &pd);
    float *image_out = malloc(w*h*pd*sizeof*image_out);

    // crop the image
    int cw, ch;
    crop(image_out, &cw, &ch, image_in, w, h, pd,
                    x0, y0, xf, yf);

    // write output
    iio_write_image_float_split(filename_out, image_out, cw, ch, pd);

    return EXIT_SUCCESS;
}

