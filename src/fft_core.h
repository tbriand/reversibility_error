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

#ifndef FFT_CORE_H
#define FFT_CORE_H

#include <complex.h>
#include <fftw3.h>

// Start threaded FFTW if FFTW_NTHREADS is defined
void init_fftw(void);
// Clean FFTW
void clean_fftw(void);
// Compute the DFT of a real-valued image
void do_fft_real(fftw_complex *out, const double *in, int nx, int ny, int nz);
// Compute the real part of the iDFT of a complex-valued image
void do_ifft_real(double *out, const fftw_complex *in, int nx, int ny, int nz);
// Compute the fftshift of a complex-valued image
void fftshift(fftw_complex *fshift, fftw_complex *fhat, int nx, int ny, int nz);
// Compute the DFT coefficients of the up-sampled image
void upsampling_fourier(fftw_complex *out, fftw_complex *in,
                        int nxin, int nyin, int nxout, int nyout, int nz, int interp);
// Up-sampling of an image using TPI
void upsampling(double *out, double *in, int nxin, int nyin, int nxout, int nyout, int nz, int interp);
// Spectrum clipping of an image
void spectrum_clipping(double *out, double *in, int nx, int ny, int nz, double r);

#endif
