#!/usr/bin/gnuplot

set term png
set grid

set output "trapezy_x(t).png"
set xlabel "t"
set ylabel "x(t)"
set title "Niejawna metoda trapezów x(t)"

plot "trapezy.dat" i 0 u 1:3 w l lw 1 t "TOL 10^{-2}", "" i 1 u 1:3 w l lw 1 t "TOL 10^{-5}"


set output "trapezy_v(t).png"
set xlabel "t"
set ylabel "v(t)"
set title "Niejawna metoda trapezów v(t)"

plot "trapezy.dat" i 0 u 1:4 w l lw 1 t "TOL 10^{-2}", "" i 1 u 1:4 w l lw 1 t "TOL 10^{-5}"


set output "trapezy_dt(t).png"
set xlabel "t"
set ylabel "dt(t)"
set title "Niejawna metoda trapezów dt(t)"

plot "trapezy.dat" i 0 u 1:2 w l lw 1 t "TOL 10^{-2}", "" i 1 u 1:2 w l lw 1 t "TOL 10^{-5}"


set output "trapezy_v(x).png"
set xlabel "x"
set ylabel "v(x)"
set title "Niejawna metoda trapezów v(x)"

plot "trapezy.dat" i 0 u 3:4 w l lw 1 t "TOL 10^{-2}", "" i 1 u 3:4 w l lw 1 t "TOL 10^{-5}"