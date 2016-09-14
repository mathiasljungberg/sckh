#!/usr/bin/python3
import matplotlib.pyplot as plt
import sys
x=[]
y=[]
for i in range(1,len(sys.argv)):
	for line in open(sys.argv[i]):
		line=line.rsplit()
		x.append(float(line[0]))
		y.append(float(line[1]))
	plt.plot(x,y,label=sys.argv[i])
	plt.legend()
plt.show()
#plt.savefig("TEST.pdf",dpi=600)
