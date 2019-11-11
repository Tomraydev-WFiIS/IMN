#!/usr/bin/gnuplot
set term png

set output "glob_calka.png"
set logscale x
set xlabel "it"
set ylabel "S(it)"
set title "Całka funkcjonalna S(it), relaksacja globalna"

plot "calka1.0_glob.dat" u 1:2 w l t "omega 1.0 maxit = 23097" , "calka0.6_glob.dat" u 1:2 w l lw 2 t "omega 0.6 maxit=36887"


set output "loc_calka.png"
set logscale x
set logscale x
set xlabel "it"
set ylabel "S(it)"
set title "Całka funkcjonalna S(it), relaksacja lokalna"

plot "calka1.0_loc.dat" u 1:2 w l t "omega 1.0 itmax = 12194" , "calka1.4_loc.dat" u 1:2 w l lw 2 t "omega 1.4 itmax = 5548", "calka1.8_loc.dat" u 1:2 w l lw 2 t "omega 1.8 itmax = 1553", "calka1.9_loc.dat" u 1:2 w l lw 2 t "omega 1.9 itmax=751"
