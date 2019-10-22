#!/bin/bash

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

# This program computes the reversibility error for an homography
# given by the displacement of the four corners
# This corresponds to Algorithm 2 with Ntransf=1
# It also computes the difference image and its DFT

# read parameters
if [ "$#" -lt "14" ]; then
    echo "usage:\n\t$0 in baseout hx1 hy1 hx2 hy2 hx3 hy3 hx4 hy4 crop interp boundary ratio"
    exit 1
fi

IN=$1
BASEOUT=$2
hx1=$3
hy1=$4
hx2=$5
hy2=$6
hx3=$7
hy3=$8
hx4=$9
hy4=${10}
crop=${11}
interp=${12}
boundary=${13}
ratio=${14}

outtiff=$BASEOUT.tiff
outpng=$BASEOUT.png
outdiff=${BASEOUT}_diff.png
outfft=${BASEOUT}_fft.png

# set the bin variable if it does not exist
if [ -z ${bin+x} ]; then
    bin="."
fi

# crop input
incrop=crop.png
crop $crop $crop -$crop -$crop $IN $incrop
crop $crop $crop -$crop -$crop $incrop $incrop

# corners shifting
w=`identify -format %w $IN`
h=`identify -format %h $IN`
w=$(($w-1))
h=$(($h-1))
hx2=`echo "$w + $hx2" | bc`
hy3=`echo "$h + $hy3" | bc`
hx4=`echo "$w + $hx4" | bc`
hy4=`echo "$h + $hy4" | bc`

# compute the corresponding homography
transform=`python3 ${bin}/hom4p.py 0 0 $hx1 $hy1 $w 0 $hx2 $hy2 0 $h $hx3 $hy3 $w $h $hx4 $hy4`

echo "Transformation by the homography: $transform"

# direct interpolation (Line 2)
interpolation $IN $outtiff "$transform" -i $interp -b $boundary -t 0 > /dev/null
crop 0 0 0 0 $outtiff $outpng # trick to save png file

# crop (Line 3)
crop $crop $crop -$crop -$crop $outtiff $outtiff

# inverse transformation (Line 4)
transform2=`python3 ${bin}/translate_homography.py $transform $crop`
interpolation $outtiff $outtiff "$transform2" -i $interp -b $boundary -t 1 > /dev/null

# crop output
crop $crop $crop -$crop -$crop $outtiff $outtiff

# compute reversibility errors (Line 5 to 7)
error=`reversibility_error $incrop $outtiff 0 -r $ratio`
echo "Reversibility error: $error"
error=`reversibility_error $incrop $outtiff 1 -r $ratio`
echo "Clipped reversibility error: $error"

# compute difference and fft
python3 ${bin}/difference_and_fft.py $incrop $outtiff $outdiff $outfft