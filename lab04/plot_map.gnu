#!/usr/bin/env gnuplot
set terminal svg

# Charge density
set output "charge.svg"
set xlabel "x"
set ylabel "y"
set title "Gęstość ładunku"

set pm3d map
set palette defined (0 "blue", 5 "white", 10 "red")
set size ratio -1

splot [0:15][0:10] "charge.dat" using 1:2:3

# 0.6
set output "map_0.6.svg"
set xlabel "x"
set ylabel "y"
set title "Mapa potencjału V(x,y). Relaksacja globalna ω = 0.6"

set pm3d map
set palette defined (0 "blue", 5 "white", 10 "red")
set size ratio -1

splot [0:15][0:10] "map_glob_0.6.dat" using 1:2:3

# 1.0
set output "map_1.0.svg"
set xlabel "x"
set ylabel "y"
set title "Mapa potencjału V(x,y). Relaksacja globalna ω = 1.0"

set pm3d map
set palette defined (0 "blue", 5 "white", 10 "red")
set size ratio -1

splot [0:15][0:10] "map_glob_1.0.dat" using 1:2:3
