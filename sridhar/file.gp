set title "Particles at the end of simulation"

set xlabel "x -line"
set ylabel "y_line"
set zlabel "z_line"
set xrange [-2000:2000]
set yrange [-2000:2000]
set zrange [-2000:2000]



splot "coord_end.txt" u 1:2:3 with points pointtype 9 lt rgb "red"

