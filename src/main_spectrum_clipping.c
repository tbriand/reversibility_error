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

#include "iio.h"
#include "fft_core.h"

#define PAR_DEFAULT_RATIO 0.01

// display help usage
void print_help(char *name)
{
    printf("\n<Usage>: %s input1 input2 [-r ratio]\n\n", name);
    printf("The optional parameter is:\n");
    printf("-r, \t Specify the ratio of clipped high-frequencies (by default %lf)\n", PAR_DEFAULT_RATIO);
}

// read command line parameters
static int read_parameters(int argc, char *argv[], char **infile, char **outfile, double *ratio)
{
    // display usage
    if (argc < 3) {
        print_help(argv[0]);
        return 0;
    }
    else {
        int i = 1;
        *infile = argv[i++];
        *outfile = argv[i++];

        // "default" value initialization
        *ratio   = PAR_DEFAULT_RATIO;
        
        //read each parameter from the command line
        while(i < argc) {
            if(strcmp(argv[i],"-r")==0)
                if(i < argc-1)
                    *ratio = atof(argv[++i]);

            i++;
        }
        
        // sanity check
        *ratio = (*ratio >= 0 && *ratio <= 1) ? *ratio : PAR_DEFAULT_RATIO;
        
        return 1;
    }
}

// Main function for computing the spectrum clipping
int main(int c, char *v[])
{
    char *filename_in, *filename_out;
    double ratio;
    
    int result = read_parameters(c, v, &filename_in, &filename_out, &ratio);

    if ( result ) {
        // initialize FFTW
        init_fftw();
        
        // read image
        int w, h, pd;
        double *in = iio_read_image_double_split(filename_in, &w, &h, &pd);
        
        // memory allocation
        double *out = malloc(w*h*pd*sizeof*out);
        
        // spectrum clipping of the input
        spectrum_clipping(out, in, w, h, pd, ratio);

        // write output image
        iio_write_image_double_split(filename_out, out, w, h, pd);
        
        // free memory
        free(in);
        free(out);
        clean_fftw();
    }

    return EXIT_SUCCESS;
}
