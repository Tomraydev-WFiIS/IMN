#!/usr/bin/env gnuplot 

set term png enhanced size 600,300 

set size ratio -1

set output "psi_map_-1000.png"
set cbr [-55:-50]
set view map
splot "Q_-1000.txt" using 1:2:3 with pm3d notitle
