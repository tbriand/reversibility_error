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

# This file contains the functions for computing a homography from
# four corresponding points

# Compute the homography sending [0,0] , [0,1], [1,1] and [1,0] to x, y, z, w.
def homography_from_4pt(x, y, z, w):
    t1 = x[0]
    t2 = z[0]
    t4 = y[1]
    t5 = t1 * t2 * t4
    t6 = w[1]
    t7 = t1 * t6
    t8 = t2 * t7
    t9 = z[1]
    t10 = t1 * t9
    t11 = y[0]
    t14 = x[1]
    t15 = w[0]
    t16 = t14 * t15
    t18 = t16 * t11
    t20 = t15 * t11 * t9
    t21 = t15 * t4
    t24 = t15 * t9
    t25 = t2 * t4
    t26 = t6 * t2
    t27 = t6 * t11
    t28 = t9 * t11
    t30 = 0.1e1 / (-t24 + t21 - t25 + t26 - t27 + t28)
    t32 = t1 * t15
    t35 = t14 * t11
    t41 = t4 * t1
    t42 = t6 * t41
    t43 = t14 * t2
    t46 = t16 * t9
    t48 = t14 * t9 * t11
    t51 = t4 * t6 * t2
    t55 = t6 * t14
     
    R = [[0]*3, [0]*3, [0]*3]
    R[0][0] = -(-t5 + t8 + t10 * t11 - t11 * t7 - t16 * t2 + t18 - t20 + t21 * t2) * t30
    R[0][1] = (t5 - t8 - t32 * t4 + t32 * t9 + t18 - t2 * t35 + t27 * t2 - t20) * t30
    R[0][2] = t1
    R[1][0] = (-t9 * t7 + t42 + t43 * t4 - t16 * t4 + t46 - t48 + t27 * t9 - t51) * t30
    R[1][1] = (-t42 + t41 * t9 - t55 * t2 + t46 - t48 + t55 * t11 + t51 - t21 * t9) * t30
    R[1][2] = t14
    R[2][0] = (-t10 + t41 + t43 - t35 + t24 - t21 - t26 + t27) * t30
    R[2][1] = (-t7 + t10 + t16 - t43 + t27 - t28 - t21 + t25) * t30
    R[2][2] = 1
    
    return R


# Compute homogaphy from 4 corresponding points.
def homography_from_4corresp(a, b, c, d, x, y, z, w):
    Hr = homography_from_4pt(a, b, c, d)
    Hl = homography_from_4pt(x, y, z, w)

    # the following code computes R = Hl * inverse Hr
    t2 = Hr[1][1]-Hr[2][1]*Hr[1][2]
    t4 = Hr[0][0]*Hr[1][1]
    t5 = Hr[0][0]*Hr[1][2]
    t7 = Hr[1][0]*Hr[0][1]
    t8 = Hr[0][2]*Hr[1][0]
    t10 = Hr[0][1]*Hr[2][0]
    t12 = Hr[0][2]*Hr[2][0]
    t15 = 1/(t4-t5*Hr[2][1]-t7+t8*Hr[2][1]+t10*Hr[1][2]-t12*Hr[1][1])
    t18 = -Hr[1][0]+Hr[1][2]*Hr[2][0]
    t23 = -Hr[1][0]*Hr[2][1]+Hr[1][1]*Hr[2][0]
    t28 = -Hr[0][1]+Hr[0][2]*Hr[2][1]
    t31 = Hr[0][0]-t12
    t35 = Hr[0][0]*Hr[2][1]-t10
    t41 = -Hr[0][1]*Hr[1][2]+Hr[0][2]*Hr[1][1]
    t44 = t5-t8
    t47 = t4-t7
    t48 = t2*t15
    t49 = t28*t15
    t50 = t41*t15
    
    R = [[0]*3, [0]*3, [0]*3]
    R[0][0] = Hl[0][0]*t48+Hl[0][1]*(t18*t15)-Hl[0][2]*(t23*t15)
    R[0][1] = Hl[0][0]*t49+Hl[0][1]*(t31*t15)-Hl[0][2]*(t35*t15)
    R[0][2] = -Hl[0][0]*t50-Hl[0][1]*(t44*t15)+Hl[0][2]*(t47*t15)
    R[1][0] = Hl[1][0]*t48+Hl[1][1]*(t18*t15)-Hl[1][2]*(t23*t15)
    R[1][1] = Hl[1][0]*t49+Hl[1][1]*(t31*t15)-Hl[1][2]*(t35*t15)
    R[1][2] = -Hl[1][0]*t50-Hl[1][1]*(t44*t15)+Hl[1][2]*(t47*t15)
    R[2][0] = Hl[2][0]*t48+Hl[2][1]*(t18*t15)-t23*t15
    R[2][1] = Hl[2][0]*t49+Hl[2][1]*(t31*t15)-t35*t15
    R[2][2] = -Hl[2][0]*t50-Hl[2][1]*(t44*t15)+t47*t15

    return R


import sys

x1 = float(sys.argv[1])
y1 = float(sys.argv[2])
hx1 = float(sys.argv[3])
hy1 = float(sys.argv[4])

x2 = float(sys.argv[5])
y2 = float(sys.argv[6])
hx2 = float(sys.argv[7])
hy2 = float(sys.argv[8])

x3 = float(sys.argv[9])
y3 = float(sys.argv[10])
hx3 = float(sys.argv[11])
hy3 = float(sys.argv[12])

x4 = float(sys.argv[13])
y4 = float(sys.argv[14])
hx4 = float(sys.argv[15])
hy4 = float(sys.argv[16])

R = homography_from_4corresp([x1, y1], [x2, y2], [x3, y3], [x4, y4], [hx1, hy1], [hx2, hy2], [hx3, hy3], [hx4, hy4])

print(R[0][0], R[0][1], R[0][2], R[1][0], R[1][1], R[1][2], R[2][0], R[2][1], R[2][2])
