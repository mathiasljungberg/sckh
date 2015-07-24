set terminal postscript enhanced color solid
#set key right off
#set size square 1,1
#set xrange [-0.4:0.4]
#set yrange [-0.4:0.4]
#set style data points
set output "fort.16.ps"

plot [][-.1:.1] "fort.16"  u 1:2 w p,\
"fort.19"  u 1:2 w l

plot [-0.8:0.8][] "fort.16"  u 1:2 w p,\
"fort.19"  u 1:2 w l

plot [-0.2:0.2][] "fort.16"  u 1:2 w p,\
"fort.19"  u 1:2 w l

plot "fort.17"  u 1:2 w p,\
"fort.20"  u 1:2 w l

plot [-0.8:0.8][] "fort.17"  u 1:2 w p,\
"fort.20"  u 1:2 w l

plot [-0.2:0.2][] "fort.17"  u 1:2 w p,\
"fort.20"  u 1:2 w l

plot "fort.18"  u 1:2 w p,\
"fort.21"  u 1:2 w l

plot [-0.8:0.8][] "fort.18"  u 1:2 w p,\
"fort.21"  u 1:2 w l

plot [-0.2:0.2][] "fort.18"  u 1:2 w p,\
"fort.21"  u 1:2 w l
