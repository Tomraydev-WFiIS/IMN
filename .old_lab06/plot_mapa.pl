#!/usr/bin/gnuplot
set term png

set output "mapa50.png"
set xlabel "x"
set ylabel "y"
set title "nx=ny=50"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:5][0:5] "mapa50.dat" i 0 u 2:1:3

reset

set output "mapa100.png"
set xlabel "x"
set ylabel "y"
set title "nx=ny=100"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:10][0:10] "mapa100.dat" i 0 u 2:1:3

reset

set output "mapa200.png"
set xlabel "x"
set ylabel "y"
set title "nx=ny=200"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:20][0:20] "mapa200.dat" i 0 u 2:1:3

reset

set output "mapa2a.png"
set xlabel "x"
set ylabel "y"
set title "2a) eps1=eps2 = 1"
set pm3d map
set palette defined (-1 "blue", 0 "white", 1 "red")
set size ratio -1

splot [0:10][0:10][-0.8:0.8] "mapa2a.dat" i 0 u 2:1:3


reset

set output "mapa2b.png"
set xlabel "x"
set ylabel "y"
set title "2b) eps1=1 eps2 = 2"
set pm3d map
set palette defined (-1 "blue", 0 "white", 1 "red")
set size ratio -1
set cbrange [-0.8:0.8]

splot [0:10][0:10][-0.8:0.8] "mapa2b.dat" i 0 u 2:1:3


reset

set output "mapa2c.png"
set xlabel "x"
set ylabel "y"
set title "2c) eps1=1 eps2 = 10"
set pm3d map
set palette defined (-1 "blue", 0 "white", 1 "red")
set size ratio -1
set cbrange [-0.8:0.8]

splot [0:10][0:10][-0.8:0.8] "mapa2c.dat" i 0 u 2:1:3
