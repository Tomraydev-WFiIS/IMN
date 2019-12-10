# ------------------------------------------------------------------
# Rysowanie gifow
# ------------------------------------------------------------------

reset
set term gif size 800,300 animate delay 10
set output "anim.gif"
n=49    #liczba klatek
set view map # widok z gory
set size ratio -1
set cbr [0:]

do for [i=0:n] {
  file = sprintf("out/zad5_it=%i.txt",i)
  splot file u 1:2:3 w pm3d  title sprintf("t=%i",i)
} 