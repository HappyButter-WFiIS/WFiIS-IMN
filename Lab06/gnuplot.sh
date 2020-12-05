#!/usr/bin/gnuplot
set term png

set out "zad2_50.png"
set title "nx=ny=50 - V(x,y)"
set xlabel "x"
set ylabel "y"
set xrange [0:5]
set yrange [0:5]
set size ratio -1
plot "zad2_50.dat" u 1:2:3 with image

reset

set out "zad2_100.png"
set title "nx=ny=100 - V(x,y)"
set xlabel "x"
set ylabel "y"
set xrange [0:10]
set yrange [0:10]
set size ratio -1
plot "zad2_100.dat" u 1:2:3 with image

reset

set out "zad2_200.png"
set title "nx=ny=200 - V(x,y)"
set xlabel "x"
set ylabel "y"
set xrange [0:20]
set yrange [0:20]
set size ratio -1
plot "zad2_200.dat" u 1:2:3 with image

reset

set out "zad3_a.png"
set title "epsilon1=epsilon2=1 - V(x,y)"
set xlabel "x"
set ylabel "y"
set xrange [0:10]
set yrange [0:10]
set zrange [-0.8:0.8]
set size ratio -1
plot "zad3_a.dat" u 1:2:3 with image

reset

set out "zad3_b.png"
set title "epsilon1=1; epsilon2=2 - V(x,y)"
set xlabel "x"
set ylabel "y"
set xrange [0:10]
set yrange [0:10]
set cbrange [-0.8:0.8]
set size ratio -1
plot "zad3_b.dat" u 1:2:3 with image

reset

set out "zad3_c.png"
set title "epsilon1=1; epsilon2=10 - V(x,y)"
set xlabel "x"
set ylabel "y"
set xrange [0:10]
set yrange [0:10]
set cbrange [-0.8:0.8]
set size ratio -1
plot "zad3_c.dat" u 1:2:3 with image