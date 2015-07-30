set terminal postscript enhanced color solid
#set key right off
#set size square 1,1
#set xrange [-0.4:0.4]
#set yrange [-0.4:0.4]
#set style data points
set output "outputfile_time_h_0_mom.dat.ps"

plot "outputfile_time_h_0_mom.dat"  u 1:2 w l

