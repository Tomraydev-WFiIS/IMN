#!/usr/bin/env gnuplot 

set term png enhanced size 600,300 

set size ratio -1

# ******** -1000 ********
set output "v_-1000.png"
set view map
splot "Q_-1000.txt" using 1:2:6 with pm3d notitle

# ******** -4000 ********
set output "v_-4000.png"
set view map
splot "Q_-4000.txt" using 1:2:6 with pm3d notitle