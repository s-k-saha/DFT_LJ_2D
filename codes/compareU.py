import numpy as np
import math

PI=np.pi
R=0.5
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
	t1=(315./1280) * (a**(-2))*(r**(-11))
	t2=(210./1280) * (a**(-4))*(r**(-9))
	t3=(168./1280) * (a**(-6))*(r**(-7))
	t4=(144./1280) * (a**(-8))*(r**(-5))
	t5=(1./10) * (a**(-10))*(r**(-3))
	
	temp=0
	if a==r:
		temp=0
	else:
		temp=math.atan((np.sqrt(a*a/r/r -1)))
	
	
	return (r*math.sqrt(a*a-r*r)*(t1+t2+t3+t4+t5)) + (315./1280)*temp*(r**(-11))


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
	if r>1.0:
		return aux_J11(r,b)-aux_J11(r,a)
	t1=(1./11) * (-b**(-11)+a**(-11))
	t2=(1./26) * (-b**(-13)+a**(-13)) * (r**2)
	t3= (1./40) * (-b**(-15)+a**(-15)) * (r**4)
	t4=	(5./272) * (-b**(-17)+a**(-17)) * (r**6)
	t5= (35./2432) * (-b**(-19)+a**(-19)) * (r**8)
	t6= (3./256)* (-b**(-21)+a**(-21))*(r**10)
	t7=(231./23552)* (-b**(-23)+a**(-23))*(r**12)
	t8= (429./51200)* (-b**(-25)+a**(-25)) * (r**14)
	
	return (t1+t2+t3+t4+t5+t6+t7+t8)
	
def J5(r,a,b):
	if r>1.0:
		return aux_J5(r,b)-aux_J5(r,a)
	t1=(1./5) * (-b**(-5)+a**(-5))
	t2=(1./14) * (-b**(-7)+a**(-7)) * (r**2)
	t3= (1./24) * (-b**(-9)+a**(-9)) * (r**4)
	t4=	(5./176) * (-b**(-11)+a**(-11)) * (r**6)
	t5= (35./1664) * (-b**(-13)+a**(-13)) * (r**8)
	t6= (21./1280)* (-b**(-15)+a**(-15)) * (r**10)
	t7= (231./17408)* (-b**(-17)+a**(-17)) * (r**12)
	t8= (429./38912)* (-b**(-19)+a**(-19)) * (r**14)
	
	return (t1+t2+t3+t4+t5+t6+t7+t8)


def U2D(r):
	if r>rc:
		return 0
	if r>rmin:
		return 8.*eps*(J11(r,r,rc)-J5(r,r,rc))
	return -2*eps*math.sqrt(rmin*rmin-r*r)+8*eps*( J11(r,rmin,rc)-J5(r,rmin,rc))#1./11*(rmin**(-11)-rc**(-11))-.2*(rmin**(-5)-rc**(-5)))
	#return -2*eps*math.sqrt(rmin*rmin-r*r) + 8*eps*(J11(r,rmin,rc)-J5(r,rmin,rc))

#print(J11(1.1225,1.1225,rc))
#print(J11(1.131,1.131,rc))
#print(aux_J11(1.1225,1.1225))

r_arr=np.linspace(0.0,10,1000)
#U3=[U3D(r) for r in r_arr]
#U1=[U1D(r) for r in r_arr]
U2=[U2D(r) for r in r_arr]

#for i in range(1000):
#	print(f"{r_arr[i]} {U1[i]} {U2[i]} {U3[i]}")

def calculate():
	rhob=.304665
	n3=rhob*4./3*PI*R*R*R
	n0=rhob
	n1=rhob*R
	n2=rhob*4*PI*R*R
	
	return -np.log(1-n3)+n2/(1-n3)*R+(n1/(1-n3)+1./24/PI/(1-n3)/(1-n3)*3*n2*n2)*4*PI*R*R+(n0/(1-n3)+n1*n2/(1-n3)/(1-n3)+1./12/PI/(1-n3)/(1-n3)*n2*n2*n2)*4./3*PI*R*R*R

'''
integral=0
for i in range(1000):
	r=.01*i
	integral+=r*.01*U2D(r)
integral=2*np.pi*integral*.597846
print(integral)
'''
print(calculate())
