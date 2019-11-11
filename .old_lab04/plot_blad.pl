#!/usr/bin/gnuplot
set term png

set output "blad1.0.png"
set xlabel "x"
set ylabel "y"
set title "Błąd rozwiazania delta(x,y, relaksacja globalna 1.0"

set pm3d map
set palette rgbformulae 33,13,10
set size ratio -1

splot [0:15][0:10] "blad1.0_glob.dat" i 0 u 1:2:3


set output "blad0.6.png"
set xlabel "x"
set ylabel "y"
set title "Błąd rozwiazania delta(x,y), relaksacja globalna 0.6"

set pm3d map
set palette rgbformulae 33,13,10
set size ratio -1

splot [0:15][0:10] "blad0.6_glob.dat" i 0 u 1:2:3
