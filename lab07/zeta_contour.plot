#!/usr/bin/env gnuplot 

set term png enhanced size 600,300 

set size ratio -1

set contour
set cntrparam levels 40 # lub ponizsze:
unset surface
set view map
unset key

# ******** -1000 ********
set output "zeta_contour_-1000.png"
set cbr [-200:300]
set cntrparam levels increment -200,10,300
splot "Q_-1000.txt" using 1:2:4:4 with lines lt -1 palette  notitle 

# ******** -4000 ********
set output "zeta_contour_-1000.png"
set cbr [-500:1000]
set cntrparam levels increment -500,50,1000
splot "Q_-4000.txt" using 1:2:4:4 with lines lt -1 palette  notitle 
