import numpy as np
import matplotlib.pyplot as plt
from numba import njit

pi=np.pi
Lz=50
R=0.5

dz=.01

Nz=int(Lz/dz)
NiR=int(R/dz)

rhob=.304665
mu=0.6045291013
muex=mu-np.log(rhob)
alpha=0.1



@njit
def omega3(i,j):
	Dz=abs(i-j)*dz
	dsds= Dz*Dz
	
	if dsds>R*R:
		return 0
	if abs(i-j)==NiR:
		return 3./8*pi*(R*R-dsds)
	if abs(i-j)==NiR-1:
		return 7./6*pi*(R*R-dsds)
	if abs(i-j)==NiR-2:
		return 23./24*pi*(R*R-dsds)
		
	return pi*(R*R-dsds)

@njit
def omega2(i,j):
	Dz=abs(i-j)*dz
	dsds= Dz*Dz
	
	if dsds>R*R:
		return 0
	
	if abs(i-j)==NiR:
		return 3./8*2*pi*R
	if abs(i-j)==NiR-1:
		return 7./6*2*pi*R
	if abs(i-j)==NiR-2:
		return 23./24*2*pi*R
	
	return 2*pi*R

@njit
def omega1(i,j):
	return omega2(i,j)/4/pi/R

@njit
def omega0(i,j):
	return omega2(i,j)/4/pi/R/R



@njit
def omega2v(i,j):
	Dz=(i-j)*dz
	dsds= Dz*Dz
	
	if dsds>R*R:
		return 0
	
	if abs(i-j)==NiR:
		return 3./8*2*pi*Dz
	if abs(i-j)==NiR-1:
		return 7./6*2*pi*Dz
	if abs(i-j)==NiR-2:
		return 23./24*2*pi*Dz
	
	return 2*pi*Dz

@njit
def omega1v(i,j):
	return omega2v(i,j)/4/pi/R


@njit
def getn(rho,n0,n1,n2,n3,n1v,n2v):		
	for i in range(Nz):
		for j in range(Nz):
			n0[i]+=(rho[j]*omega0(i,j))*dz
			n1[i]+=(rho[j]*omega1(i,j))*dz
			n2[i]+=(rho[j]*omega2(i,j))*dz
			n3[i]+=(rho[j]*omega3(i,j))*dz
			n1v[i]+=(rho[j]*omega1v(i,j))*dz
			n2v[i]+=(rho[j]*omega2v(i,j))*dz
	
			

@njit
def getc1_fmt(c1,n0,n1,n2,n3,n1v,n2v):
	dphidn0 = -np.log(1-n3)
	dphidn1 =  n2/(1-n3)
	dphidn2 = n1/(1-n3) + (3*n2*n2 - 3*n2v*n2v)/(24*pi*(1-n3)*(1-n3))
	dphidn3 = n0/(1-n3) + (n1*n2 - n1v*n2v)/(1-n3)/(1-n3) + (n2*n2*n2 - 3*n2*n2v*n2v)/12/pi/(1-n3)/(1-n3)/(1-n3)
	dphidn1v = -n2v/(1-n3)
	dphidn2v = -n1v/(1-n3) - 3*n2*n2v/(12*pi*(1-n3)*(1-n3)) 
	
	
	for i in range(Nz):
		for j in range(Nz):
			c1[i]+= -(dphidn0[j]*omega0(i,j) + dphidn1[j]*omega1(i,j) +dphidn2[j]*omega2(i,j) +dphidn3[j]*omega3(i,j)  - dphidn1v[j]*omega1v(i,j) - dphidn2v[j]*omega2v(i,j))*dz

@njit
def filterc(c1):
	c1temp=c1[Nz-NiR*10]
	
	for i in range(Nz-NiR*10+1,Nz):
		c1[i]=c1temp

@njit
def iterate(rho,Vext):
	for i in range(NiR-1):
		rho[i] = 0	
	n0,n1,n2,n3,n1v,n2v = np.zeros(Nz),np.zeros(Nz),np.zeros(Nz),np.zeros(Nz),np.zeros(Nz),np.zeros(Nz)
	c1 = np.zeros(Nz)
	getn(rho,n0,n1,n2,n3,n1v,n2v)
	getc1_fmt(c1,n0,n1,n2,n3,n1v,n2v)
	filterc(c1)
	
	
	#return (c1[int(Nz/2)])
	rhonew = rhob * np.exp(-Vext + c1 + muex)
	
	rhonew=alpha*rhonew + (1-alpha)*rho
	
	return rhonew

@njit
def getVext():
	Vext=np.zeros(Nz)
	for i in range(NiR-1):
		Vext[i] = 1000.
	return Vext 

@njit
def init_rho():
	Vext=getVext()
	rho = np.ones(Nz)*rhob
	for i in range(Nz):
		rho[i]=rho[i]*np.exp(-Vext[i])
	return rho
	
@njit
def getdata():
	rho=init_rho()
	Vext=getVext()
	
	#h=iterate(rho,Vext)
	for _ in range(100):
		rho=iterate(rho,Vext)
	
	return rho#h


if __name__=="__main__":
	rho=getdata()
	#print(rho)
	z=np.arange(0,Lz,dz)
	for i in range(Nz):
		print(f"{z[i]} {rho[i]}")
	#plt.imshow(rho)
	
	#plt.colorbar()
	#plt.show()	
		
		
		
	
			
			
			
			
			
			
			
			
			
			
			
	
