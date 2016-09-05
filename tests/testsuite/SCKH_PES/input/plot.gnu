set terminal postscript enhanced color solid
#set key right off
#set size square 1,1
#set xrange [-0.4:0.4]
#set yrange [-0.4:0.4]
#set style data points
set output "plot.ps"

plot "outputfile_time_h_0.dat"  u 1:2 w l lt 1 lw 4\
,"outputfile_time_h_10.dat"  u 1:2 w l lt 2 lw 4\
,"outputfile_time_h_20.dat"  u 1:2 w l lt 3 lw 4\
,"outputfile_time_h_30.dat"  u 1:2 w l lt 4 lw 4\

#,"outputfile_time_h_40.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_50.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_60.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_70.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_80.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_90.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_100.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_110.dat"  u 1:2 w l lt 1 lw 4\

#,"outputfile_time_h_12.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_13.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_14.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_15.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_16.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_17.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_18.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_19.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_20.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_21.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_22.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_23.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_24.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_25.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_26.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_27.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_28.dat"  u 1:2 w l lt 1 lw 4\
#,"outputfile_time_h_29.dat"  u 1:2 w l lt 1 lw 4\
