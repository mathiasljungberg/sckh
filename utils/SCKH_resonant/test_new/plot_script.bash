ls -1 spectrum_resonant_sc_*.*.dat > tmp
#cat fort.31 > tmp

n=(`wc tmp`)
n=${n[0]}

echo $n

rm tmp.gnu

a=`sed -n 1p tmp` 
echo "\
#set format x \"%e\"; set format y \"%e\"
set terminal postscript enhanced color
set output \"tmp.ps\"
d=0.0 #1e11
#set nokey
plot [180:192][] \"$a\" u 1:(\$2 +1*d) w l\\" >> tmp.gnu

for ((i=2;i<=n; i++ ))
do

a=`sed -n ${i}p tmp` 
#echo ",\"$a\" u 1:(\$2 +$i*d) w l\\" >> tmp.gnu
echo " 

plot [180:192][] \"$a\" u 1:(\$2 +$i*d) w l\\" >> tmp.gnu

done
