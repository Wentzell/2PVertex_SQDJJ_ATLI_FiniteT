set term postscript eps enhanced 35 color solid
set xlabel ".XLABEL"
set ylabel ".YLABEL"

set output "plots/.FILENAME.eps"
plot "data/.FILENAME.dat" u 1:2 t ""