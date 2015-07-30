set terminal postscript enhanced color solid
#set key right off
#set size square 1,1
#set xrange [-0.4:0.4]
#set yrange [-0.4:0.4]
#set style data points
set output "fort.19.ps"

plot "fort.19"  u 1:5 w l\

plot "fort.19"  u 1:2 w l\
,"fort.19"  u 1:3 w l\
,"fort.19"  u 1:4 w l\
