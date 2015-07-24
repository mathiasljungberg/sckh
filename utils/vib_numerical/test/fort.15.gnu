set terminal postscript enhanced color solid
#set key right off
#set size square 1,1
#set xrange [-0.4:0.4]
#set yrange [-0.4:0.4]
#set style data points
set output "fort.15.ps"

plot "fort.15"  u 1:2 w l,\
	"fort.15"  u 1:3 w l,\
	"fort.15"  u 1:4 w l	

plot "fort.15"  u 1:2 w l
plot "fort.15"  u 1:3 w l
plot "fort.15"  u 1:4 w l

plot "fort.16"  u 1:2 w p
plot "fort.17"  u 1:2 w p
plot "fort.18"  u 1:2 w p


