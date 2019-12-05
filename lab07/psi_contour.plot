#!/usr/bin/env gnuplot 

set term png enhanced size 600,300 

set size ratio -1

set contour
set cntrparam levels 40 # lub ponizsze:
unset surface
set view map
unset key

# ******** -1000 ********
set output "psi_contour_-1000.png"
set cbr [-55:-50]
set cntrparam levels increment -55,0.2,50
splot "Q_-1000.txt" using 1:2:3:3 with lines lt -1 palette  notitle

# ******** -4000 ********
set output "psi_contour_-4000.png"
set cbr [-218:-202]
set cntrparam levels increment -218,1,-202
splot "Q_-4000.txt" using 1:2:3:3 with lines lt -1 palette  notitle

# ******** +4000 ********
set output "psi_contour_+4000.png"
set cbr [218:202]
set cntrparam levels increment 218,1,202
splot "Q_+4000.txt" using 1:2:3:3 with lines lt -1 palette  notitle