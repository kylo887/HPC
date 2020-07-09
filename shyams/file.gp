set title "Visualize the particles at t=0"
unset hidden3d
set autoscale
set parametric
set key box
set xlabel "x-axis"
set ylabel "y-axis"
set zlabel "z-axis"
set xrange [-200:200]
set yrange [-200:200]
set zrange [-200:200]



splot "time_0.txt" u 1:2:3 with points pointtype 7

