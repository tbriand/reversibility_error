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

# This file contains the function for translating a homography
# by a shift of t (useful to handle crops)

# Translate a homography by t
def translate_homography(h0, h1, h2, h3, h4, h5, h6, h7, h8, t):
    R = [[0]*3, [0]*3, [0]*3]
    R[2][2] = h6*t + h7*t + h8
    R[0][0] = h0 - t*h6
    R[0][1] = h1 - t*h7
    R[0][2] = (h0*t + h1*t + h2) - t*R[2][2]
    R[1][0] = h3 - t*h6
    R[1][1] = h4 - t*h7
    R[1][2] = (h3*t + h4*t + h5) - t*R[2][2]
    R[2][0] = h6
    R[2][1] = h7
    return R

import sys

h0 = float(sys.argv[1])
h1 = float(sys.argv[2])
h2 = float(sys.argv[3])
h3 = float(sys.argv[4])
h4 = float(sys.argv[5])
h5 = float(sys.argv[6])
h6 = float(sys.argv[7])
h7 = float(sys.argv[8])
h8 = float(sys.argv[9])
t  = float(sys.argv[10])

R = translate_homography(h0, h1, h2, h3, h4, h5, h6, h7, h8, t)

print(R[0][0], R[0][1], R[0][2], R[1][0], R[1][1], R[1][2], R[2][0], R[2][1], R[2][2])
