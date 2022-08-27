# plot.plt
set nokey
set grid
set xlabel "x"
set ylabel "y"
m="numerovdata.dat"
plot m using 1:2 with linespoints
