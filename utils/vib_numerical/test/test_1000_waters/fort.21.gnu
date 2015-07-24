set terminal postscript enhanced color solid
#set key right off
#set size square 1,1
#set xrange [-0.4:0.4]
#set yrange [-0.4:0.4]
#set style data points
set output "fort.21.ps"

plot\
 "fort.21"  u 1:4 w l,\
 "fort.22"  u 2:4 w l,\
 "fort.23"  u 3:4 w l


