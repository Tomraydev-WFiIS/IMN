#!/usr/bin/gnuplot 

set term png enhanced size 600,300 

set size ratio -1

set o "psi-1000.png"
set contour
#set cntrparam levels 40 # lub ponizsze:
set cbr [-200:350]
set cntrparam levels increment -200,10,350
unset surface
set view map
unset key
sp "zad_-1000.0.txt" u 1:2:3:3 w l lt -1 palette  t '' 
