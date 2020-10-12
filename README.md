# Reversibility Error of Image Interpolation Methods: Definition and Improvements #

## Summary ##
This repository contains programs for the evaluation of the reversibility error
of interpolation methods.
It is part of an [IPOL publication](https://doi.org/10.5201/ipol.2019.277)

## Authors ##

* Thibaud Briand <briand.thibaud@gmail.com>

Laboratoire d'Informatique Gaspard Monge (LIGM)/ Ecole des Ponts ParisTech
Centre de mathématiques et de leurs applications (CMLA)/ ENS Paris-Saclay

## Versions ##

Version 1.0,  released on 2019-10-16 : Original version.

Version 1.01, released on 2020-10-12 : Fix requirements for the IPOL demo.

## License ##

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2018-2019, Thibaud Briand <briand.thibaud@gmail.com>

All rights reserved.

## Build ##

Required environment: Any unix-like system with a standard compilation
environment (make and C compiler) and [Cmake](https://cmake.org/).

Required libraries:
[libpng](http://libpng.org/pub/png/libpng.html),
[lipjpeg](http://ijg.org/),
[libtiff](http://simplesystems.org/libtiff/)
[libfftw3](http://www.fftw.org/)

Optional libraries:
[libopenmp](https://www.openmp.org/)
[GSL](https://www.gnu.org/software/gsl/)

Build instructions:

     mkdir build
     cd build
     cmake -DCMAKE_BUILD_TYPE=Release ..
     make

It produces programs "create_burst", "crop", "interpolation", "reversibility_error" and "spectrum_clipping".

## Usage of create_burst ##

The program reads an input image, a number of images, optionnally takes some parameters and
produces a burst of images.

    <Usage>: ./create_burst input base number [OPTIONS]

	 Output images are written as base_%i.tiff
	 Output cropped images are written as base_crop_%i.tiff
	 Homographies are written as base_%i.hom
	 Homographies after crop are written as base_crop_%i.hom
	 The first image of the generated sequence is the reference (identity)

The optional parameters are:
-c,      Specify the crop size (by default 0)
-L,      Specify the displacement of the corners (by default 3)
-i,      Specify the interpolation method (by default p+s-spline11-spline1)
-b,      Specify the boundary condition between hsym, wsym, per and constant (by default hsym)
-t,      Specify the type of transformation (2 translation, 3 euclidean, 6 affinity, 8 homography (default))
-z,      Specify the down-sampling factor (by default 1)
-n,      Specify the standard deviation of the noise (by default 0)
-s,      Specify the seed of the random generator (by default 0)

Execution examples:

  1.  Burst of 10 + 1 images:

       ./create_burst input.png base 11

  2.  Burst of 100 + 1 images with a crop of 20 pixels using the zoomed bicubic interpolation method:

       ./create_burst input.png base 101 -c 20 -i bic-z2

## Usage of crop ##

The program reads an input image and the four bounds of the image domain,
and performs the crop of the image.

   <Usage>: ./crop x0 y0 xf yf in out

Execution examples:

  1.  Crop in a band of 20 pixels:

	./crop 10 10 -10 -10 input.png output.png

  2.  Remove the last column of an image:

	./crop 0 0 -1 0 input.png output.png

## Usage of interpolation ##

The program reads an input image, an homography, optionnally takes some parameters
and produces the corresponding transformed image.
This may use Algorithm 3 or Algorithm 4.

   <Usage>: ./interpolation input output "h11 h12 h13 h21 h22 h23 h31 h32 h33" [OPTIONS]

The optional parameters are:
-i,      Specify the interpolation method (by default p+s-spline11-spline1)
-b,      Specify the boundary condition between hsym, wsym, per and constant (by default hsym)
-t,      Set to 1 to apply the inverse transform (by default 0)

Execution examples:

  1.  Translation of (1.5,-2.3):

       ./interpolation input.png output.tiff "1 0 1.5 0 1 -2.3 0 0 1"

  2.  Apply the inverse of an homography using bicubic interpolation with periodic extension:

       ./interpolation input.png output.tiff "h11 h12 h13 h21 h22 h23 h31 h32 h33" -i bic -b periodic -t 1

## Usage of reversibility_error ##

The program reads two input images and computes the reversibility error (or clipped reversibility error).

   <Usage>: ./reversibility_error input1 input2 clipped -r ratio

	 Set clipped to 1 to compute the clipped reversibility error
The optional parameter is:
-r,      Specify the ratio of clipped high-frequencies (by default 0.050000)

Execution examples:

  1.  Reversibility error:

       ./reversibility_error input1.tiff input2.tiff 0

  2.  Clipped reversibility error with ratio 10%:

       ./reversibility_error input1.tiff input2.tiff 1 -r 0.1

## Usage of spectrum_clipping ##

The program reads an input image and its spectrum clipped version.

   <Usage>: ./spectrum_clipping input1 input2 [-r ratio]

The optional parameter is:
-r,     Specify the ratio of clipped high-frequencies (by default 0.050000)

Execution examples:

  1.  Spectrum clipping with ratio 10%:

       ./spectrum_clipping input.png output.tiff -r 0.1

## Usage of the demo script run.sh

The script reads an image, the displacement of the image four corners,
the interpolation method and the reversibility error parameters.
It computes the transformed image, the reversibility errors, the difference images
and its DFT.

   <Usage>: bin=path_to_demo_directory ./run.sh in baseout hx1 hy1 hx2 hy2 hx3 hy3 hx4 hy4 crop interp boundary ratio

Execution example:

       bin=demo/ demo/run.sh input.png baseout 1 1 -1 -1 0 0 1 1 20 spline11 hsym 0.01

Note that the execution of the script requires the programs "crop", "interpolation" and
"reversibility_error". The user should add these programs to the Path. For instance
this can be done using the command

    PATH=path_to_build_directory:$PATH

## List of files ##

* CMakeList.txt          : CMake file for the compilation
* FindGSL.cmake		 : CMAke file for finding GNU scientific library GSL
* LICENSE		 : License file
* README.txt             : This file

In the data/ directory:

* rubberwhale_gray.png   : Test image (grayscale)
* baboon.png		 : Test image (color)

In the demo/ directory:

* difference_and_fft.py       : Program to display the difference image and its spectrum in the demo
* hom4p.py                    : Program to compute a homography from four corresponding points
* requirements.txt            : Package used by the demo
* run.sh                      : Main script of the demo. This corresponds to Algorithm 2 with Ntransf=1
* translate_homography.py     : Program to translate a homography

In the src/ directory:

* bicubic.[hc]                : Functions to perform bicubic interpolation
* compute_core.h	      : Utility functions for the crop and the rmse
* fft_core.[hc]               : Functions related to the Fourier computations
* homography_core.[hc]	      : Functions related to homographies (contains Algorithm 1)
* interpolation_core.[hc]     : Functions to perform a geometric transformation using interpolation (contains Algorithm 3 and Algorithm 4)
* main_create_burst.c         : Main program for creating a burst from an image (in particular it generates random homographies)
* main_crop.c                 : Main program for cropping an image
* main_interpolation.c        : Main program for input/ouput (Algorithm 3 and Algorithm 4)
* main_reversibility_error.c  : Main program for computing the reversibility error
* main_spectrum_clipping.c    : Main program for computing the spectrum clipping
* periodic_plus_smooth.[hc]   : Functions to compute the periodic plus smooth decomposition of an image
* tpi.[hc]                    : Functions to perform trigonometric polynomial interpolation

Additional files are provided in the external/ directory:

* cmphomod.h             : Functions to compute an homography from four correspondences
* iio.[hc]               : Functions for opening images in any format
* random.h               : Functions to generate random numbers
* xmtime.h               : Clock with millisecond precision
* bspline                : Functions to perform B-spline interpolation
* nfft-3.5.0             : NFFT3 library
