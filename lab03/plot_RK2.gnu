#!/usr/bin/gnuplot

set grid

set style line 1 lc rgb '#f25900' linetype 1 linewidth 2 pointtype 7 pi 1 pointsize 1
set style line 2 lc rgb '#5060d0' linetype 1 linewidth 2 pointtype 1 pi 2 pointsize 1

# x(t)
set terminal svg
set output "RK2_x(t).svg"
set xlabel "t"
set ylabel "x(t)"
set title "Metoda RK2"

plot "RK2.dat" index 0 using 1:3 with lines linewidth 2 title "TOL 10^{-2}" , \
    "" index 1 using 1:3 with lines linewidth 2 title "TOL 10^{-5}"

# set terminal x11 0
# set output
# replot


# v(t)
set terminal svg
set output "RK2_v(t).svg"
set xlabel "t"
set ylabel "v(t)"
set title "Metoda RK2"

plot "RK2.dat" index 0 using 1:4 with lines linewidth 2 title "TOL 10^{-2}", \
    "" index 1 using 1:4 with lines linewidth 2 title "TOL 10^{-5}"

# set terminal x11 1
# set output
# replot

# dt(t)
set terminal svg
set output "RK2_dt(t).svg"
set xlabel "t"
set ylabel "dt(t)"
set title "Metoda RK2"

plot "RK2.dat" index 0 using 1:2 with linespoints linestyle 1 title "TOL 10^{-2}", \
    "" index 1 using 1:2 with linespoints linestyle 2 title "TOL 10^{-5}"

# set terminal x11 2
# set output
# replot

# v(x)
set terminal svg
set output "RK2_v(x).svg"
set xlabel "x"
set ylabel "v(x)"
set title "Metoda RK2"

plot "RK2.dat" index 0 using 3:4 with lines linewidth 2 title "TOL 10^{-2}", \
    "" index 1 using 3:4 with lines linewidth 2 title "TOL 10^{-5}"

# set terminal x11 3
# set output
# replot
# pause -1