set terminal postscript enhanced color solid
#set key right off
#set size square 1,1
#set xrange [-0.4:0.4]
#set yrange [-0.4:0.4]
#set style data points
set output "fort.12.ps"

#set mapping cylindrical
#set dgrid3d 100,100,10
#set grid polar
#set cntrparam levels auto 100
#set contour surface
set view 0, 0, 1, 1
#unset surface
#set hidden3d
#set colorbox horiz user origin .1,.1 size .8,.04
#set parametric
#set pm3d


splot "fort.12"  u 2:3:4 w p

