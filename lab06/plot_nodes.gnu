#!/usr/bin/env gnuplot

set xl offset 0, 1
set yl offset 1
set title offset 0,-0.6
set ytics 1 offset 0.5
set xtics 1 offset 0, 0.4

set term png enhanced size 350,350 font "Serif, 11" 

set size ratio -1

set title "wezly-numeracja"
set xl "i"
set yl "j"

# kolejne kolumny w pliku zawierajÄ: indeks l, indeks i, indeks j
set output "wezly.png"
plot "zad1.txt" using 2:3 with points pointtype 7 pointsize 1.2 notitle , \
     "" u ($2+0.2):($3+0.2):1 with labels notitle 
