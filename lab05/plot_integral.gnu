#!/usr/bin/env gnuplot
set terminal svg

set output "integral.svg"
set xlabel "it"
set ylabel "S"
set title "S(it)"

plot "integral.dat" index 0 using 2:3 with lines title "k=16", \
     "integral.dat" index 1 using 2:3 with lines title "k=8", \
     "integral.dat" index 2 using 2:3 with lines title "k=4", \
     "integral.dat" index 3 using 2:3 with lines title "k=2", \
     "integral.dat" index 4 using 2:3 with lines title "k=1"
