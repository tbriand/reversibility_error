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

#include "random.h"
#include "iio.h"
#include "xmtime.h"
#include "interpolation_core.h"
#include "homography_core.h"
#include "fft_core.h"

#define PAR_DEFAULT_L 3
#define PAR_DEFAULT_TYPE 8
#define PAR_DEFAULT_ZOOM 1
#define PAR_DEFAULT_CROP 0
#define PAR_DEFAULT_SIGMA 0
#define PAR_DEFAULT_SEED 0

// display help usage
void print_help(char *name)
{
    printf("\nusage:\n\t%s input base number [OPTIONS]\n\n", name);
    printf("\t Output images are written as base_%%i.tiff\n");
    printf("\t Output cropped images are written as base_crop_%%i.tiff\n");
    printf("\t Homographies are written as base_%%i.hom\n");
    printf("\t Homographies after crop are written as base_crop_%%i.hom\n");
    printf("\t The first image of the generated sequence is the reference (identity)\n");
    printf("\nThe optional parameters are:\n");
    printf("-c, \t Specify the crop size (by default %i)\n", PAR_DEFAULT_CROP);
    printf("-L, \t Specify the displacement of the corners (by default %i)\n", PAR_DEFAULT_L);
    printf("-i, \t Specify the interpolation method (by default p+s-spline11-spline1)\n");
    printf("-b, \t Specify the boundary condition between hsym, wsym, per and constant (by default hsym)\n");
    printf("-t, \t Specify the type of transformation (2 translation, 3 euclidean, 6 affinity, 8 homography (default))\n");
    printf("-z, \t Specify the down-sampling factor (by default %i)\n", PAR_DEFAULT_ZOOM);
    printf("-n, \t Specify the standard deviation of the noise (by default %i)\n", PAR_DEFAULT_SIGMA);
    printf("-s, \t Specify the seed of the random generator (by default %i)\n", PAR_DEFAULT_SEED);
}

// read command line parameters
static int read_parameters(int argc, char *argv[], char **infile, char **outfile,
                           int *n, char **interp, char **boundary, double *L, int *type,
                           double *zoom, int *crop, double *sigma, unsigned long *seed)
{
    // display usage
    if (argc < 4) {
        print_help(argv[0]);
        return 0;
    }
    else {
        int i = 1;
        *infile  = argv[i++];
        *outfile = argv[i++];
        *n = atoi(argv[i++]);

        // "default" value initialization
        *interp   = "p+s-spline11-spline1";
        *boundary = "hsym";
        *L         = PAR_DEFAULT_L;
        *type      = PAR_DEFAULT_TYPE;
        *zoom      = PAR_DEFAULT_ZOOM;
        *sigma     = PAR_DEFAULT_SIGMA;
        *crop      = PAR_DEFAULT_CROP;
        *seed      = PAR_DEFAULT_SEED;
        
        //read each parameter from the command line
        while(i < argc) {
            if(strcmp(argv[i],"-i")==0)
                if(i < argc-1)
                    *interp = argv[++i];

            if(strcmp(argv[i],"-b")==0)
                if(i < argc-1)
                    *boundary = argv[++i];

            if(strcmp(argv[i],"-L")==0)
                if(i < argc-1)
                    *L = atof(argv[++i]);
                
            if(strcmp(argv[i],"-t")==0)
                if(i < argc-1)
                    *type = atoi(argv[++i]);
                
            if(strcmp(argv[i],"-z")==0)
                if(i < argc-1)
                    *zoom = atof(argv[++i]);

            if(strcmp(argv[i],"-n")==0)
                if(i < argc-1)
                    *sigma = atof(argv[++i]);
                
            if(strcmp(argv[i],"-c")==0)
                if(i < argc-1)
                    *crop = atoi(argv[++i]);
                
            if(strcmp(argv[i],"-s")==0)
                if(i < argc-1)
                    *seed = atoi(argv[++i]);

            i++;
        }
        
        // sanity check
        *L = (*L >= 0) ? *L : PAR_DEFAULT_L;
        *zoom = (*zoom > 0) ? *zoom : PAR_DEFAULT_ZOOM;
        *sigma = (*sigma >= 0) ? *sigma : PAR_DEFAULT_SIGMA;
        *crop = (*crop >= 0) ? *crop : PAR_DEFAULT_CROP;

        return 1;
    }
}

int main(int c, char *v[])
{
    char *filename_in, *base_out, *interp, *boundary;
    int n, type, crop;
    unsigned long seed;
    double L, zoom, sigma;
    
    int result = read_parameters(c, v, &filename_in, &base_out, &n, &interp, &boundary,
                                 &L, &type, &zoom, &crop, &sigma, &seed);

    if ( result ) {
        // initialize FFTW
        init_fftw();
        
        // read data
        int w, h, pd;
        double *in = iio_read_image_double_split(filename_in, &w, &h, &pd);

        // initialize time
        unsigned long t1 = xmtime();
        
        // rand initialization
        xsrand(seed);

        // create all the homographies
        // it is done in the beginning so that for a given seed we always have the same homographies
        double *homographies = malloc(9*n*sizeof(double));
        double H[9], H2[9];
        char filename_homo[500];

        // randomly compute homographies
        for (int j = 0; j < n; j++) {
            if ( j > 0 )
                    create_random_transformation(H, L, w, h, type);
            else { // first transformation is the identity
                H[0] = 1.0;
                H[1] = 0.0;
                H[2] = 0.0;
                H[3] = 0.0;
                H[4] = 1.0;
                H[5] = 0.0;
                H[6] = 0.0;
                H[7] = 0.0;
                H[8] = 1.0;
            }

            // save the homography
            for (int i = 0; i < 9; i++)
                homographies[9*j+i] = H[i];

            // write it in a file by taking into account the eventual zoom
            zoom_homography(H2, H, 1.0/zoom, 1.0/zoom);
            sprintf(filename_homo, "%s_%i.hom", base_out, j+1);
            FILE *f = fopen(filename_homo, "w");
            for (int i = 0; i < 9; i++) {
                fprintf(f,"%1.16lg%c", H2[i], i==8 ? '\n' : ' ');
            }
            fclose(f);

            // crop case
            if ( crop ) {
                translate_homography(H, H2, crop, crop);

                // save homography and write it in a file
                sprintf(filename_homo, "%s_crop_%i.hom", base_out, j+1);
                FILE *g = fopen(filename_homo, "w");
                for (int i = 0; i < 9; i++)
                    fprintf(g,"%1.16lg%c", H[i], i==8 ? '\n' : ' ');
                fclose(g);
            }
        }

        // Boudary condition 
        BoundaryExt boundaryExt = read_ext(boundary);  

        // compute images
        int wout = w/zoom;
        int hout = h/zoom;
        char filename_out[500];
        double *out = malloc(wout*hout*pd*sizeof*out);
        
        for (int j = 0; j < n; j++) {
            for(int i = 0; i < 9; i++)
                H[i] = homographies[9*j+i];
            
            // apply geometric transformation to the input
            interpolate_image_homography(out, in, w, h, pd, H, interp, boundaryExt, zoom);
            
            // add noise
            if ( sigma > 0 )
                for(int i = 0; i < wout*hout*pd; i++)
                    out[i] += sigma*random_normal();
            
            // write transformed image
            sprintf(filename_out, "%s_%i.tiff", base_out, j+1);
            iio_write_image_double_split(filename_out, out, wout, hout, pd);
            
            // crop case
            if ( crop ) {
                int wcrop = wout - 2*crop;
                int hcrop = hout - 2*crop;
                double *out_crop = malloc(wcrop*hcrop*pd*sizeof(double));
                
                for(int l = 0; l < pd; l++)
                    for(int q = 0; q < hcrop; q++)
                        for(int p = 0; p < wcrop; p++)
                            out_crop[p + q*wcrop + l*wcrop*hcrop] = out[p + crop + (q+crop)*wout + l*wout*hout];
                
                sprintf(filename_out, "%s_crop_%i.tiff", base_out, j+1);
                iio_write_image_double_split(filename_out, out_crop, wcrop, hcrop, pd);
                free(out_crop);
            }
            
        }

        // final time and print time
        unsigned long t2 = xmtime();
        printf("Burst created in %.3f seconds \n", (float) (t2-t1)/1000);

        // free memory
        free(in);
        free(homographies);
        free(out);
        clean_fftw();
    }
    
    return EXIT_SUCCESS;
}
