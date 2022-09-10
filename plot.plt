# plot.plt
set nokey
set grid
set xlabel "r"
set ylabel "prob dens"
#set xrange [0:10]
m="solutiondata.dat"
plot m u 1:($2**2), "" u 1:3
