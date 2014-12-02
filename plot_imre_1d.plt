set xlabel ".XLABEL"
set ylabel ".YLABEL"
set xtics out nomirror
set ytics out nomirror

set term postscript eps enhanced 20 color solid

set output "plots/.FILENAME_IMRE.eps"

set multiplot layout 1,2

plot "data/.FILENAME_RE.dat" u 1:2 t "" with lines lw 3

plot "data/.FILENAME_IM.dat" u 1:2 t "" with lines lw 3

unset multiplot
