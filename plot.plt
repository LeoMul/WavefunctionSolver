# plot.plt
set nokey
set grid
set xlabel "r"
set ylabel "prob dens"
#set xrange [0:10]
m="solutiondata.dat"
f(x) = 4*(x**2)*exp(-2*x)
plot m u 1:($2**2),f(x)
