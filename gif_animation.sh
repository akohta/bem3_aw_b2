#!/bin/bash

# settings
Delay=10
Loop=0

if [ $# -ne 1 ]; then
 echo "Usage : $0 image_directory"
 exit
fi

if [ ! -d $1 ]; then
 echo "$1 does not exist or is not a directory."
 exit
fi

convert -delay $Delay -loop $Loop $1/xy_p_*.png xy_p.gif &
convert -delay $Delay -loop $Loop $1/yz_p_*.png yz_p.gif &
convert -delay $Delay -loop $Loop $1/xz_p_*.png xz_p.gif &

wait
