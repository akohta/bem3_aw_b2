# Gnuplot script file
set terminal postscript eps color enhanced "Arial" 20 size 7in,7in
set output "I_example2.eps"

set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')

set pm3d map

set size square
set multiplot layout 2,2

set xrange[-0.04 : 0.04]
set yrange[-0.04 : 0.04]
set xtics -0.04, 0.02, 0.04
set ytics -0.04, 0.02, 0.04

set title "sound pressure intensity on y=0 plane"
set xlabel "{/Arial-Italic x}"
set ylabel "{/Arial-Italic z}"
splot "Ip_xz.txt"

set title "sound pressure intensity on x=0 plane"
set xlabel "{/Arial-Italic y}"
splot "Ip_yz.txt"

set title "sound pressure intensity on z=0 plane"
set xlabel "{/Arial-Italic x}"
set ylabel "{/Arial-Italic y}"
splot "Ip_xy.txt"

unset multiplot
set terminal x11
