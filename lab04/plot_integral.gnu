#!/usr/bin/env gnuplot
set terminal svg

set output "integral_glob.svg"
set logscale x
set xlabel "nr iteracji"
set ylabel "S"
set title "Relaksacja globalna"

plot "integral_glob_1.0.dat" using 1:2 with lines linewidth 2 title "ω = 1.0 it = 23096" , \
     "integral_glob_0.6.dat" using 1:2 with lines linewidth 2 title "ω = 0.6 it = 36886"

set output "integral_loc.svg"
set logscale x
set xlabel "nr iteracji"
set ylabel "S"
set title "Relaksacja lokalna"

plot "integral_loc_1.0.dat" using 1:2 with lines linewidth 2 title "ω = 1.0, 12193 it" , \
     "integral_loc_1.4.dat" using 1:2 with lines linewidth 2 title "ω = 1.4, 5547 it", \
     "integral_loc_1.8.dat" using 1:2 with lines linewidth 2 title "ω = 1.8, 1552 it", \
     "integral_loc_1.9.dat" using 1:2 with lines linewidth 2 title "ω = 1.9, 750 it"
