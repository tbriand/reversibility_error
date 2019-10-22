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
#include "compute_core.h"

#define PAR_DEFAULT_RATIO 0.01

// display help usage
void print_help(char *name)
{
    printf("\n<Usage>: %s input1 input2 clipped -r ratio\n\n", name);
    printf("\t Set clipped to 1 to compute the clipped reversibility error\n");
    printf("The optional parameter is:\n");
    printf("-r, \t Specify the ratio of clipped high-frequencies (by default %lf)\n", PAR_DEFAULT_RATIO);
}

// read command line parameters
static int read_parameters(int argc, char *argv[], char **infile1, char **infile2,
                           int *clipped, double *ratio)
{
    // display usage
    if (argc < 4) {
        print_help(argv[0]);
        return 0;
    }
    else {
        int i = 1;
        *infile1 = argv[i++];
        *infile2 = argv[i++];
        *clipped = atoi(argv[i++]);

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

// Main function for computing the reversibility error (Definition 7)
int main(int c, char *v[])
{
    char *filename_in, *filename_in2;
    int clipped;
    double ratio;
    
    int result = read_parameters(c, v, &filename_in, &filename_in2, &clipped, &ratio);

    if ( result ) {
        // initialize FFTW
        init_fftw();
        
        // read image
        int w, h, pd, w2, h2, pd2;
        double *in = iio_read_image_double_split(filename_in, &w, &h, &pd);
        double *in2 = iio_read_image_double_split(filename_in2, &w2, &h2, &pd2);
        
        // check sizes
        if( w != w2 || h != h2 || pd != pd2 ) {
            fprintf(stderr, "Images must have the same size (reversibility error)\n");
            return EXIT_FAILURE;
        }
        
        // compute the difference image
        for(int i = 0; i < w*h*pd; i++)
            in[i] -= in2[i];
        
        // optionnal spectrum clipping 
        if ( clipped )
            spectrum_clipping(in, in, w, h, pd, ratio);
        
        // compute the error
        double reversibility_error = rmse(in, w*h*pd);

        // print error
        printf("%1.14lg\n", reversibility_error);
        
        // free memory
        free(in);
        free(in2);
        clean_fftw();
    }

    return EXIT_SUCCESS;
}
