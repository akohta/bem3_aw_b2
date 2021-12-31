# Gnuplot script file
set terminal postscript eps color enhanced "Arial" 20 size 7in,4in
set output "I_example4.eps"

set encoding iso

set angle degree
set polar
set grid polar
set border 0
set view equal xy
set size square
set noxtics
set noytics
set ttics format "%g{\260}"
set ttics 0,30,330
set rrange [0.0 : 1]
unset raxis
set key outside

set title "normalized sound pressure intensity"
plot "fsIp_yz.txt" using 1:3 w l title "x=0 plane, 0{\260}: z-axis, 270{\260}: y-axis", "fsIp_xz.txt" using 1:3 w l title "y=0 plane, 0{\260}: z-axis,   90{\260}: x-axis", "fsIp_xy.txt" using 1:3 w l title "z=0 plane, 0{\260}: x-axis,   90{\260}: y-axis"
