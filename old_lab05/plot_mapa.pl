#!/usr/bin/gnuplot
set term png

set output "mapa16.png"
set xlabel "x"
set ylabel "y"
set title "Rozwiązanie: potencjał V(x,y), k=16"
reset
set pm3d map
set palette defined (-1 "blue", 0 "white", 1 "red")
set size ratio -1

splot [0:30][0:30] "mapy.dat" i 0 u 1:2:3

reset
set output "mapa8.png"
set xlabel "x"
set ylabel "y"
set title "Rozwiązanie: potencjał V(x,y), k=8"

set pm3d map
set palette defined (0 "blue", 5 "white", 10 "red")
set size ratio -1

splot [0:30][0:30] "mapy.dat" i 1 u 1:2:3


reset
set output "mapa4.png"
set xlabel "x"
set ylabel "y"
set title "Rozwiązanie: potencjał V(x,y), k=4"

set pm3d map
set palette defined (0 "blue", 5 "white", 10 "red")
set size ratio -1

splot [0:30][0:30] "mapy.dat" i 2 u 1:2:3


reset
set output "mapa2.png"
set xlabel "x"
set ylabel "y"
set title "Rozwiązanie: potencjał V(x,y), k=2"

set pm3d map
set palette defined (0 "blue", 5 "white", 10 "red")
set size ratio -1

splot [0:30][0:30] "mapy.dat" i 3 u 1:2:3


reset
set output "mapa1.png"
set xlabel "x"
set ylabel "y"
set title "Rozwiązanie: potencjał V(x,y), k=1"

set pm3d map
set palette defined (0 "blue", 5 "white", 10 "red")
set size ratio -1

splot [0:30][0:30] "mapy.dat" i 4 u 1:2:3