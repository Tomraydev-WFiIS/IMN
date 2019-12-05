#!/usr/bin/env gnuplot 

set term png enhanced size 600,300 

set size ratio -1

set output "zeta_map_-1000.png"
set cbr [-200:300]
set view map
splot "Q_-1000.txt" using 1:2:4 with pm3d notitle
