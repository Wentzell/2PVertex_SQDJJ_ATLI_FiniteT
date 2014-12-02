set pm3d map

set size square
set xlabel "{/Symbol w}_1"
set xlabel "{/Symbol w}_2"
set xtics out nomirror
set ytics out nomirror
#set term x11
set palette rgbformulae 22,13,-31

set cbrange[-5:5]

splot "data/.FILENAME_RE.dat" u 1:2:3 t ""

MIN_RE = GPVAL_Z_MIN
MAX_RE = GPVAL_Z_MAX #3.5

splot "data/.FILENAME_IM.dat" u 1:2:3 t ""

MIN_IM = GPVAL_Z_MIN
MAX_IM = GPVAL_Z_MAX #3.5

set term postscript eps enhanced 13 color solid
set output "plots/.FILENAME_IMRE.eps"

set multiplot layout 1,2

#set size 0.5,0.5
set cblabel "Re {/Symbol g}"
set palette defined (MIN_RE "black",(3*MIN_RE/4 + MAX_RE/4) "blue",(2*MIN_RE/4 + 2*MAX_RE/4)"white", (MIN_RE/4 + 3*MAX_RE/4) "red", MAX_RE "dark-red") #(-0.003 "blue", 0 "white", 0.012 "red") #
#set palette defined (MIN_RE "black", (MIN_RE/2) "blue",0"white", MAX_RE/2 "red", MAX_RE "dark-red") #(-0.003 "blue", 0 "white", 0.012 "red") #
splot "data/.FILENAME_RE.dat" u 1:2:3 t ""

#set size 0.5,0.5
set cblabel "Im {/Symbol g}"
set palette defined (MIN_IM "black",(3*MIN_IM/4 + MAX_IM/4) "blue",(2*MIN_IM/4 + 2*MAX_IM/4)"white", (MIN_IM/4 + 3*MAX_IM/4) "red", MAX_IM "dark-red") #(-0.003 "blue", 0 "white", 0.012 "red") #
#set palette defined (MIN_IM "black", (MIN_IM/2) "blue",0"white", MAX_IM/2 "red", MAX_IM "dark-red") #(-0.003 "blue", 0 "white", 0.012 "red") #
splot "data/.FILENAME_IM.dat" u 1:2:3 t ""

unset multiplot
