#!/usr/bin/gnuplot
set term png

set output "mapa1.0.png"
set xlabel "x"
set ylabel "y"
set title "Rozwiązanie: potencjał V(x,y), relaksacja globalna 1.0"

set pm3d map
set palette defined (0 "blue", 5 "white", 10 "red")
set size ratio -1

splot [0:15][0:10] "mapa1.0_glob.dat" i 0 u 1:2:3


set output "mapa0.6.png"
set xlabel "x"
set ylabel "y"
set title "Rozwiązanie: potencjał V(x,y), relaksacja globalna 0.6"

set pm3d map
set palette defined (0 "blue", 5 "white", 10 "red")
set size ratio -1

splot [0:15][0:10] "mapa0.6_glob.dat" i 0 u 1:2:3

