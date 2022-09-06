# plot.plt
set nokey
set grid
set xlabel "x"
set ylabel "y"
m="solutiondata.dat"
plot m u 1:2, "" u 1:3
