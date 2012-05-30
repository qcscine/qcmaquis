# !/bin/bash

gnuplot << EOF
set terminal postscript eps color 
set key left top
set output "InnerProduct.eps"
set xlabel "Order"
set ylabel "Time [s]"
set xtics 1
set title "InnerProduct $1"

#red family
set style line 1  lt 1 lw 3 pt 1 lc rgb "#FF0000" #red
set style line 2  lt 1 lw 3 pt 1 lc rgb "#FF1493" #deep pink
set style line 3  lt 1 lw 3 pt 1 lc rgb "#FF6347" #tomato
set style line 4  lt 1 lw 3 pt 1 lc rgb "#FFC8CB" #pink
                  
#blue family      
set style line 5  lt 1 lw 3 pt 1 lc rgb "#0000FF" #blue
set style line 6  lt 1 lw 3 pt 1 lc rgb "#1E90FF" #dodger blue
set style line 7  lt 1 lw 3 pt 1 lc rgb "#4682B4" #steelblue
set style line 8  lt 1 lw 3 pt 1 lc rgb "#E6E6FA" #lavender

#green family
set style line 9  lt 1 lw 3 pt 1 lc rgb "#008000" #green
set style line 10 lt 1 lw 3 pt 1 lc rgb "#3CB371" #medium seagreen
set style line 11 lt 1 lw 3 pt 1 lc rgb "#7FFFd4" #Aquamarine
set style line 12 lt 1 lw 3 pt 1 lc rgb "#7FFF00" #chartreuse:w

set out "InnerProduct.eps";plot   "data.gnuplot.2"  using 2:3 t "IP_128_256" w l ls 1
set out "InnerProduct.eps";replot "data.gnuplot.3"  using 2:3 t "IP_192_384" w l ls 5
set out "InnerProduct.eps";replot "data.gnuplot.4"  using 2:3 t "IP_256_512" w l ls 9
set output
EOF
