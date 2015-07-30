set terminal postscript enhanced color solid
#set key right off
#set size square 1,1
#set xrange [-0.4:0.4]
#set yrange [-0.4:0.4]
#set style data points
set output "outputfile_sigma.dat.ps"

plot [518:528] "outputfile_sigma.dat"  u 1:2 w l

plot "outputfile_sigma.dat"  u 1:2 w l

