#!/usr/bin/env gnuplot 

set term png enhanced size 600,300 

set size ratio -1

set output "vx_map.png"
set view map
set title "vx"
splot "vx.dat" using 1:2:3 with pm3d

set output "vy_map.png"
set view map
set title "vy"
splot "vy.dat" using 1:2:3 with pm3d