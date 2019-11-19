#!/usr/bin/env gnuplot

set terminal svg

reset
set output "map_16.svg"
set xlabel "x"
set ylabel "y"
set title "k=16"

set pm3d map
set palette defined ( 0 "blue", 4 "cyan", 5 "green", 6 "yellow", 10 "red" )
set size ratio -1

splot [0:25.6][0:25.6] "maps.dat" index 0 using 1:2:3 notitle

reset
set output "map_8.svg"
set xlabel "x"
set ylabel "y"
set title "k=8"

set pm3d map
set palette defined ( 0 "blue", 4 "cyan", 5 "green", 6 "yellow", 10 "red" )
set size ratio -1

splot [0:25.6][0:25.6] "maps.dat" index 1 using 1:2:3 notitle


reset
set output "map_4.svg"
set xlabel "x"
set ylabel "y"
set title "k=4"

set pm3d map
set palette defined ( 0 "blue", 4 "cyan", 5 "green", 6 "yellow", 10 "red" )
set size ratio -1

splot [0:25.6][0:25.6] "maps.dat" index 2 using 1:2:3 notitle


reset
set output "map_2.svg"
set xlabel "x"
set ylabel "y"
set title "k=2"

set pm3d map
set palette defined ( 0 "blue", 4 "cyan", 5 "green", 6 "yellow", 10 "red" )
set size ratio -1

splot [0:25.6][0:25.6] "maps.dat" index 3 using 1:2:3 notitle


reset
set output "map_1.svg"
set xlabel "x"
set ylabel "y"
set title "k=1"

set pm3d map
set palette defined ( 0 "blue", 4 "cyan", 5 "green", 6 "yellow", 10 "red" )
set size ratio -1

splot [0:25.6][0:25.6] "maps.dat" index 4 using 1:2:3 notitle
