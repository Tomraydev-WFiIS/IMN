#!/usr/bin/gnuplot

set grid

set style line 1 lc rgb '#f25900' linetype 1 linewidth 2 pointtype 7 pi 1 pointsize 1
set style line 2 lc rgb '#5060d0' linetype 1 linewidth 2 pointtype 1 pi 2 pointsize 1

# x(t)
set terminal png
set output "trapezy_x(t).png"
set xlabel "t"
set ylabel "x(t)"
set title "metoda trapezów"

plot "trapezy.dat" i 0 using 1:3 with lines linewidth 1 title "TOL 10^{-2}", \
    "" i 1 using 1:3 with lines linewidth 1 title "TOL 10^{-5}"

# set terminal x11 0
# set output
# replot

# v(t)
set terminal png
set output "trapezy_v(t).png"
set xlabel "t"
set ylabel "v(t)"
set title "metoda trapezów"

plot "trapezy.dat" i 0 using 1:4 with lines linewidth 1 title "TOL 10^{-2}", \
    "" i 1 using 1:4 with lines linewidth 1 title "TOL 10^{-5}"

# set terminal x11 1
# set output
# replot

# dt(t)
set terminal png
set output "trapezy_dt(t).png"
set xlabel "t"
set ylabel "dt(t)"
set title "metoda trapezów"

plot "trapezy.dat" i 0 using 1:2 with linespoints linestyle 1 title "TOL 10^{-2}", \
    "" i 1 using 1:2 with linespoints linestyle 2 title "TOL 10^{-5}"

# set terminal x11 2
# set output
# replot

# v(x)
set terminal png
set output "trapezy_v(x).png"
set xlabel "x"
set ylabel "v(x)"
set title "metoda trapezów"

plot "trapezy.dat" i 0 using 3:4 with lines linewidth 1 title "TOL 10^{-2}", \
    "" i 1 using 3:4 with lines linewidth 1 title "TOL 10^{-5}"

# set terminal x11 3
# set output
# replot
# pause -1