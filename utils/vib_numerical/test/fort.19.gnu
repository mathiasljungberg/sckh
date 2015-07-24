set terminal postscript enhanced color solid
#set key right off
#set size square 1,1
#set xrange [-0.4:0.4]
#set yrange [-0.4:0.4]
#set style data points
set output "fort.19.ps"

plot "fort.19"  u 1:2 w l

plot [][0:5]"fort.19"  u 1:2 w l

plot [0.7:1.3][] "fort.19"  u 1:2 w l

plot [0.8:1.0][] "fort.19"  u 1:2 w l

plot [0.94:0.97][] "fort.19"  u 1:2 w l


plot [104:106][] "fort.20"  u 1:2 w l