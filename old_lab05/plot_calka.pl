#!/usr/bin/gnuplot
set term png

set output "calka.png"
set xlabel "it"
set ylabel "S(it)"
set title "Całka funkcjonalna S(it)"

plot "calka.dat" u 2:3 w l t ""

