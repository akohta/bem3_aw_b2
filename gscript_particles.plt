# Gnuplot script file for visualization of particles file.
# usage on command line : gnuplot -e 'nodes="particles_filename"' gscript_particles.plg

set terminal postscript eps color enhanced "Arial" 14 size 3.5in,3.5in
set output "particles.eps"

if ( !exists("nodes") ) nodes="ex.particles"

set view equal xyz
set noborder
set noxtics
set noytics
set noztics

splot nodes with points pt 0 lc 0 notitle 
