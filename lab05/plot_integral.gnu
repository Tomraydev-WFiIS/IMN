#!/usr/bin/env gnuplot
set terminal svg

set output "integral.svg"
set xlabel "it"
set ylabel "S"
set title "S(it)"

plot "integral.dat" using 2:3 with lines
