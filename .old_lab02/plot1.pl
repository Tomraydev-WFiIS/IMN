#!/usr/bin/gnuplot


set term png
set grid

set output "zad1_wykres.png"
set xlabel "t"
set ylabel "u(t), z(t)"
set title "Niejawna metoda trapezów z iteracją Picarda max it = 5"

plot "zad1_rozw.dat" u 1:2 w l lw 2 t "u(t)", "" u 1:3 w l lw 2 t "z(t)"