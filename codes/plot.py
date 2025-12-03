import numpy as np
import matplotlib.pyplot as plt

fname='../data/ljw_dz0.1.dat'
data=np.zeros((200,200))
k=1
with open(fname,"r") as fl:
	
	for line in fl:
		if(k>=40000):
			break
		dat=line.split(" ")
		i=int(dat[0])
		j=int(dat[1])
		val=float(dat[2])
		data[i][j]=val
		k+=1

plt.imshow(data)
plt.colorbar()
plt.show()

y=data[100][:]
x=np.arange(0,20.,.1)
plt.plot(x,y)
plt.show()
