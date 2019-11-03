#!/usr/bin/gnuplot

set grid
set term png


# 
set output "RK2_x(t).png"
set xlabel "t"
set ylabel "x(t)"
set title "Metoda RK2"

plot "RK2.dat" i 0 u 1:3 w l lw 1 t "TOL 10^{-2}", \
    "" i 1 u 1:3 w l lw 1 t "TOL 10^{-5}"


set output "RK2_v(t).png"
set xlabel "t"
set ylabel "v(t)"
set title "Metoda RK2"

plot "RK2.dat" i 0 u 1:4 w l lw 1 t "TOL 10^{-2}", \
    "" i 1 u 1:4 w l lw 1 t "TOL 10^{-5}"


set output "RK2_dt(t).png"
set xlabel "t"
set ylabel "dt(t)"
set title "Metoda RK2"

plot "RK2.dat" i 0 u 1:2 w l lw 1 t "TOL 10^{-2}", \
    "" i 1 u 1:2 w l lw 1 t "TOL 10^{-5}"


set output "RK2_v(x).png"
set xlabel "x"
set ylabel "v(x)"
set title "Metoda RK2"

plot "RK2.dat" i 0 u 3:4 w l lw 1 t "TOL 10^{-2}", \
    "" i 1 u 3:4 w l lw 1 t "TOL 10^{-5}"