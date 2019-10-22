#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-2.0+
#
# Thibaud Briand <briand.thibaud@gmail.com>
#
# Copyright (c) 2018-2019, Thibaud Briand
# All rights reserved.
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# You should have received a copy of the GNU General Pulic License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# This file contains the function for displaying the difference image
# and the corresponding spectrum

import sys
import numpy as np
from skimage import io

# read images
image1 = io.imread(sys.argv[1])
image2 = io.imread(sys.argv[2])

# difference image
diff = image1 - image2

# fft in log scale
fft = np.fft.fft2(diff,axes=(0,1))
fft = np.fft.fftshift(fft,(0,1))
fft = np.log(1 + np.abs(fft))

fft = 15*fft
fft = np.array(fft).astype('uint8')
diff = 100*diff + 128
diff = np.array(diff).astype('uint8')

# save images
io.imsave(sys.argv[3],diff)
io.imsave(sys.argv[4],fft)
