set pm3d map

set size square
set xlabel ".XLABEL"
set ylabel ".YLABEL"
set cblabel ".FILENAME"
set xtics out nomirror
set ytics out nomirror
set term x11
set palette rgbformulae 22,13,-31
splot "data/.FILENAME.dat" u 1:2:3 t ""

set palette defined (GPVAL_Z_MIN "black",(3*GPVAL_Z_MIN/4 + 3.5/4) "blue",(2*GPVAL_Z_MIN/4 + 2*3.5/4)"white", (GPVAL_Z_MIN/4 + 3*3.5/4) "red", 3.5 "dark-red") #(-0.003 "blue", 0 "white", 0.012 "red") #
set term postscript eps enhanced 35 color solid
set output "plots/.FILENAME.eps"
splot "data/.FILENAME.dat" u 1:2:3 t ""