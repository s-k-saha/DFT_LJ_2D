import numpy as np
import matplotlib.pyplot as plt
from numba import njit

pi=np.pi
Lz=10
Lx=10
R=.5
sigma=2*R
dz=dx=.1

Nx=int(Lx/dx)
Nz=int(Lz/dz)
NiR=int(R/dz)

rhob=.304665/sigma/sigma/sigma #Nigel Wilding's densities are mentioned in units of sigma^3
mu=0.6045291013
muex=mu-np.log(rhob)
alpha=0.01

#@njit
def init_rho():
	return np.ones((Nx,Nz))*rhob

#@njit
def dist(i1,j1,i2,j2):
	Dx=abs(i2-i1)*dx
	Dx=min(Dx,Lx-Dx)
	Dz=abs(j2-j1)*dz
	return Dx*Dx+Dz*Dz

#@njit
def omega3(i1,j1,i2,j2):
	dsds=dist(i1,j1,i2,j2)
	
	if dsds>=R*R:
		return 0
	return 2*np.sqrt(R*R-dsds)

#@njit
def omega2(i1,j1,i2,j2):
	dsds=dist(i1,j1,i2,j2)
	
	if dsds>=R*R:
		return 0
	return 2*R/np.sqrt(R*R-dsds)
	
#@njit
def omega1(i1,j1,i2,j2):
	dsds=dist(i1,j1,i2,j2)
	
	if dsds>=R*R:
		return 0
	return 1/2/pi/np.sqrt(R*R-dsds)

#@njit
def omega0(i1,j1,i2,j2):
	dsds=dist(i1,j1,i2,j2)
	
	if dsds>=R*R:
		return 0
	return 1/2/pi/R/np.sqrt(R*R-dsds)

#@njit
def omega1v(i1,j1,i2,j2):
	dsds=dist(i1,j1,i2,j2)
	ds=np.sqrt(dsds)
	if dsds>=R*R:
		return 0
	return ds/2/pi/R/np.sqrt(R*R-dsds)

#@njit
def omega2v(i1,j1,i2,j2):
	dsds=dist(i1,j1,i2,j2)
	ds=np.sqrt(dsds)
	if dsds>=R*R:
		return 0
	return 2*ds/np.sqrt(R*R-dsds)

#@njit
def getn(rho,n0,n1,n2,n3,n1v,n2v):		
	for i1 in range(Nx):
		for j1 in range(Nz):
			for i2 in range(Nx):
				for j2 in range(Nz):
					n0[i1][j1]+=(rho[i2][j2]*omega0(i1,j1,i2,j2))*dx*dz
					n1[i1][j1]+=(rho[i2][j2]*omega1(i1,j1,i2,j2))*dx*dz
					n2[i1][j1]+=(rho[i2][j2]*omega2(i1,j1,i2,j2))*dx*dz
					n3[i1][j1]+=(rho[i2][j2]*omega3(i1,j1,i2,j2))*dx*dz
					n1v[i1][j1]+=(rho[i2][j2]*omega1v(i1,j1,i2,j2))*dx*dz
					n2v[i1][j1]+=(rho[i2][j2]*omega2v(i1,j1,i2,j2))*dx*dz
					

#@njit
def getc1_fmt(c1,n0,n1,n2,n3,n1v,n2v):
	dphidn0 = -np.log(1-n3)
	dphidn1 =  n2/(1-n3)
	dphidn2 = n1/(1-n3) + (3*n2*n2 - 3*n2v*n2v)/(24*pi*(1-n3)*(1-n3))
	dphidn3 = n0/(1-n3) + (n1*n2 - n1v*n2v)/(1-n3)/(1-n3) + (n2*n2*n2 - 3*n2*n2v*n2v)/12/pi/(1-n3)/(1-n3)/(1-n3)
	dphidn1v = -n2v/(1-n3)
	dphidn2v = -n1v/(1-n3) - 3*n2*n2v/(12*pi*(1-n3)*(1-n3)) 
	
	
	for i1 in range(Nx):
		for j1 in range(Nz):
			for i2 in range(Nx):
				for j2 in range(Nz):
					c1[i1][j1]+= -(dphidn0[i2][j2]*omega0(i2,j2,i1,j1) + dphidn1[i2][j2]*omega1(i2,j2,i1,j1) +dphidn2[i2][j2]*omega2(i2,j2,i1,j1) +dphidn3[i2][j2]*omega3(i2,j2,i1,j1)  - dphidn1v[i2][j2]*omega1v(i2,j2,i1,j1) - dphidn2v[i2][j2]*omega2v(i2,j2,i1,j1))*dx*dz


def pn(ni):
	plt.imshow(ni)
	plt.colorbar()
	plt.grid()
	plt.show()	
	
	
#@njit
def iterate(rho,Vext):
	for i in range(Nx):
		for j in range(NiR-1):
			rho[i][j] = 0	
	
	n0,n1,n2,n3,n1v,n2v = np.zeros((Nx,Nz)),np.zeros((Nx,Nz)),np.zeros((Nx,Nz)),np.zeros((Nx,Nz)),np.zeros((Nx,Nz)),np.zeros((Nx,Nz))
	c1 = np.zeros((Nx,Nz))
	getn(rho,n0,n1,n2,n3,n1v,n2v)
	
	pn(n3)
	
	
	getc1_fmt(c1,n0,n1,n2,n3,n1v,n2v)
	
	rhonew = rhob * np.exp(-Vext + c1 + muex)
	
	rhonew=alpha*rhonew + (1-alpha)*rho
	
	return rhonew

#@njit
def getVext():
	Vext=np.zeros((Nx,Nz))
	for i in range(Nx):
		for j in range(NiR-1):
			Vext[i][j] = 1000.
	return Vext 

#@njit
def getdata():
	rho=init_rho()
	Vext=getVext()
	
	for _ in range(1):
		rho=iterate(rho,Vext)
	
	return rho


if __name__=="__main__":
	rho=getdata()
	
	
	plt.imshow(rho)
	plt.colorbar()
	plt.grid()
	plt.show()	
		
		
		
	
			
			
			
			
			
			
			
			
			
			
			
	
