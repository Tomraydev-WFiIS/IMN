#!/usr/bin/gnuplot


set term png
set grid

set output "zad3_wykres.png"
set xlabel "t"
set ylabel "u(t), z(t)"
set title "Niejawna metoda RK2 max it = 3"

plot "zad3_rozw.dat" u 1:2 w l lw 2 t "u(t)", "" u 1:3 w l lw 2 t "z(t)"