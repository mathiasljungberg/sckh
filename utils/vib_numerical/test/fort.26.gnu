set terminal postscript enhanced color solid
#set key right off
#set size square 1,1
#set xrange [-0.4:0.4]
#set yrange [-0.4:0.4]
#set style data points
set output "fort.26.ps"

plot "fort.26"  u 1:2 w l

plot "fort.27"  u 1:2 w l

plot "fort.28"  u 1:2 w l






