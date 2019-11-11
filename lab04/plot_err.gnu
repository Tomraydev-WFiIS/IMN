#!/usr/bin/env gnuplot
set terminal svg

# 0.6
set output "err_0.6.svg"
set xlabel "x"
set ylabel "y"
set title "Relaksacja globalna ω = 0.6"

set pm3d map
set palette rgbformulae 33,13,10
set size ratio -1

splot [0:15][0:10] "err_glob_0.6.dat" using 1:2:3

# 1.0
set output "err_1.0.svg"
set xlabel "x"
set ylabel "y"
set title "Relaksacja globalna ω = 1.0"

set pm3d map
set palette rgbformulae 33,13,10
set size ratio -1

splot [0:15][0:10] "err_glob_1.0.dat" using 1:2:3
