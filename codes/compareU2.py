import numpy as np
import math

rc=2.5
rmin=2.0**(1./6)
eps=1.0

def U3D(r):
	if r>rc:
		return 0
	if r>rmin:
		return 4*eps*(r**(-12)-r**(-6)) 
	return -eps

def U1D(r):
	if r>rc:
		return 0
	if r>rmin:
		return .4*eps*(r**(-10)) - eps*(r**(-4)) - .4*eps*(rc**(-10))+eps*(rc**(-4))
	return eps/2*(r*r-rmin*rmin) + .4*eps*(rmin**(-10)) - eps*(rmin**(-4))-.4*eps*(rc**(-10))+eps*(rc**(-4)) 



def aux_J11(r,a):
	temp=0
	if a==r:
		temp=0
	else:
		temp=math.atan((np.sqrt(a*a/r/r -1)))
		
	fl=(63.*np.pi/512)
	#print(f"fl={fl}")
	t1=r*math.sqrt(a*a-r*r)*(315./1280) * (a**(-2)) * (r**(0))+(315./1280*temp)* (r**(0))
	#print(f"t1={t1}")
	#t1=(r*math.sqrt(a*a-r*r)*t1)+(315./1280)*temp)
	t1=t1-fl
	print(f"t1={t1}")
	#print(f"t1={t1}")
	t2=(210./1280) * (a**(-4))*(r**(2))
	t3=(168./1280) * (a**(-6))*(r**(4))
	t4=(144./1280) * (a**(-8))*(r**(6))
	t5=(1./10) * (a**(-10))*(r**(8))
	
	
	
	
	return (r*math.sqrt(a*a-r*r)*(t2+t3+t4+t5)+t1)


def aux_J5(r,a):
	temp=0
	if a==r:
		temp=0.
	else:	
		temp=math.atan(np.sqrt(a*a/r/r -1))
	
	t1=(3./8)*(a**(-2))*(r**(-5))
	t2=(1./4)*(a**(-4))*(r**(-3))
		
	return (r*np.sqrt(a*a-r*r)*(t1+t2)) + (3./8)*(r**(-5))*temp
	
	
def J11(r,a,b):
	if r>=rmin:
		return aux_J11(r,b)-aux_J11(r,a)
	return (aux_J11(r,b)-aux_J11(r,a))*(r**(-11))
	'''
	t1=(1./11) * (-b**(-11)+a**(-11))
	t2=(1./26) * (-b**(-13)+a**(-13)) * (r**2)
	t3= (1./40) * (-b**(-15)+a**(-15)) * (r**4)
	t4=	(5./272) * (-b**(-17)+a**(-17)) * (r**6)
	t5= (35./2432) * (-b**(-19)+a**(-19)) * (r**8)
	
	return (t1+t2+t3+t4+t5)
	'''
	
def J5(r,a,b):
	if r>rmin:
		return aux_J5(r,b)-aux_J5(r,a)
	t1=(1./5) * (-b**(-5)+a**(-5))
	t2=(1./14) * (-b**(-7)+a**(-7)) * (r**2)
	t3= (1./24) * (-b**(-9)+a**(-9)) * (r**4)
	t4=	(5./176) * (-b**(-11)+a**(-11)) * (r**6)
	t5= (35./1664) * (-b**(-13)+a**(-13)) * (r**8)
	return (t1+t2+t3+t4+t5)


def U2D(r):
	if r>rc:
		return 0
	if r>rmin:
		return 8.*eps*(J11(r,r,rc)-J5(r,r,rc))
	return -2*eps*math.sqrt(rmin*rmin-r*r)+8*eps*( J11(r,rmin,rc)-J5(r,rmin,rc))#1./11*(rmin**(-11)-rc**(-11))-.2*(rmin**(-5)-rc**(-5)))
	#return -2*eps*math.sqrt(rmin*rmin-r*r) + 8*eps*(J11(r,rmin,rc)-J5(r,rmin,rc))

#print(J11(1.1225,1.1225,rc))
#print(J11(1.131,1.131,rc))
print(aux_J11(0.01,rc))
print(aux_J11(0.01,rmin))
#print(J11(0.01,rmin,rc))
'''
r_arr=np.linspace(0.0,10,1000)
U3=[U3D(r) for r in r_arr]
U1=[U1D(r) for r in r_arr]
U2=[U2D(r) for r in r_arr]

for i in range(1000):
	print(f"{r_arr[i]} {U1[i]} {U2[i]} {U3[i]}")
'''
