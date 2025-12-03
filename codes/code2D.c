#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define MIN(x,y) ((x<y)? x : y)
#define MAX(x,y) ((x>y)? x : y)

#define PI 3.14159265358979
#define R 0.5
#define sigma 1
#define dz 0.025
#define dx 0.025
#define Lx 10
#define Lz 10
#define mu 0.6045291013
#define rhob 0.304665

#define Nx 400
#define Nz 400
#define NiR 20
#define R2 0.25

#define N 160000 //Nx*Nz
#define NiR2 1600 //4*NiR*NiR

/*
const int Nx = (int)(Lx/dx);
const int Nz = (int)(Lz/dz);
const int NiR = (int) (R/dz);
const double R2=R*R;
*/
//double rhob=.304665/sigma/sigma/sigma;
//double muex=mu-log(rhob);
const double alpha=0.01;

double rho[N],rhonew[N],Vext[N];
double n0[N],n1[N],n2[N],n3[N],n1v[N],n2v[N];
double dphidn0[N],dphidn1[N],dphidn2[N],dphidn3[N],dphidn1v[N],dphidn2v[N];
double c1[N];

int nbrarr[N][NiR2];
double distarr[N][NiR2];

void initnbr();
double dist(int,int,int,int);
double omega3(double);
double omega2(double);
double omega1(double);
double omega0(double);
double omega1v(double);
double omega2v(double);
void getn(), getc1_fmt(), getVext(), rhoinit(), iterate();

double dist(int i1,int j1,int i2,int j2)
{
	double Dx=fabs(i2-i1)*dx;
	double Dz=fabs(j2-j1)*dz;
	return Dx*Dx+Dz*Dz;
}

void initnbr()
{
	int k;
	double dsds;
	int i1,j1,i2,j2;
	for(int i=0;i<N;i++)
	{
		i1=i/Nz;
		j1=i%Nz;
		k=0;
		for(i2=MAX(0,i1-NiR);i2<=MIN(i1+NiR,Nx-1);i2++)
			for(j2=MAX(0,j1-NiR);j2<=MIN(j1+NiR,Nz-1);j2++)
				{
					dsds=dist(i1,j1,i2,j2);
					if(dsds<R2)
						{
							nbrarr[i][k]=i2*Nz+j2;
							distarr[i][k]=dsds;
							k+=1;
						}
				}
	}
}
	
	

double omega3(double dsds)
{
	if(dsds>=R2)
	return 0;
	return 2*sqrt(R2-dsds);
}

double omega2(double dsds)
{
	if(dsds>=R2)
	return 0;
	return 2*R/sqrt(R2-dsds);
}

double omega1(double dsds)
{
	
	if(dsds>=R2)
	return 0;
	return 1/2/PI/sqrt(R2-dsds);
}

double omega0(double dsds)
{
	
	if(dsds>=R2)
	return 0;
	
	return 1/2/PI/R/sqrt(R2-dsds);
}

double omega1v(double dsds)
{
	
	if(dsds>=R2)
	return 0;
	double ds=sqrt(dsds);
	
	return ds/2/PI/R/sqrt(R2-dsds);
}

double omega2v(double dsds)
{
	
	if(dsds>=R2)
	return 0;
	double ds=sqrt(dsds);
	
	return 2*ds/sqrt(R*R-dsds);
}


void getn()
{

	for(int i=0;i<N;i++)
	{
		n0[i]=0;
		n1[i]=0;
		n2[i]=0;
		n3[i]=0;
		n1v[i]=0;
		n2v[i]=0;
	}
	
	for(int i=0;i<N;i++)
	{
		for(int k=0;k<NiR2;k++)
		{
			int j=nbrarr[i][k];
			
			if(j==0)
			continue;
			
			double dsds=distarr[i][k];
			n0[i]+=(rho[j]*omega0(dsds))*dx*dz;
			n1[i]+=(rho[j]*omega1(dsds))*dx*dz;
			n2[i]+=(rho[j]*omega2(dsds))*dx*dz;
			n3[i]+=(rho[j]*omega3(dsds))*dx*dz;
			n1v[i]+=(rho[j]*omega1v(dsds))*dx*dz;
			n2v[i]+=(rho[j]*omega2v(dsds))*dx*dz;
			
			
		}
	}
}

void getc1_fmt()
{
	for(int i=0;i<N;i++)
	{
		dphidn0[i] = -log(1-n3[i]);
		dphidn1[i] =  n2[i]/(1-n3[i]);
		dphidn2[i] = n1[i]/(1-n3[i]) + (3*n2[i]*n2[i] - 3*n2v[i]*n2v[i])/(24*PI*(1-n3[i])*(1-n3[i]));
		dphidn3[i] = n0[i]/(1-n3[i]) + (n1[i]*n2[i] - n1v[i]*n2v[i])/(1-n3[i])/(1-n3[i]) + (n2[i]*n2[i]*n2[i] - 3*n2[i]*n2v[i]*n2v[i])/12/PI/(1-n3[i])/(1-n3[i])/(1-n3[i]);
		dphidn1v[i] = -n2v[i]/(1-n3[i]);
		dphidn2v[i] = -n1v[i]/(1-n3[i]) - 3*n2[i]*n2v[i]/(12*PI*(1-n3[i])*(1-n3[i]));
	}	
	
	
	
	for(int i=0;i<N;i++)
	c1[i]=0;
	
	
	for(int i=0;i<N;i++)
	{
		for(int k=0;k<NiR2;k++)
		{
			int j=nbrarr[i][k];
			
			if(j==0)
			continue;
			
			double dsds=distarr[i][k];
			
			c1[i]+= -(dphidn0[j]*omega0(dsds) + dphidn1[j]*omega1(dsds) +dphidn2[j]*omega2(dsds) +dphidn3[j]*omega3(dsds)  - dphidn1v[j]*omega1v(dsds) - dphidn2v[j]*omega2v(dsds))*dx*dz;
	
		}
	}
				
}	


void getVext(){
	for(int i=0;i<Nx;i++)
		for(int j=0;j<NiR;j++)
			Vext[i*Nz+j] = 1000.;
}

void rhoinit()
{
	for(int i=0;i<N;i++)
		rho[i]=rhob*exp(-Vext[i]);
}

void iterate(){	
	getn();
	getc1_fmt();
	double muex=mu-log(rhob);
	for(int i=0;i<N;i++)
		rhonew[i]=rhob*exp(-Vext[i]+c1[i]+muex);
	
	for(int i=0;i<N;i++)
		rho[i]=alpha*rhonew[i] + (1-alpha)*rho[i];
}




void main()
{
	clock_t start = clock();
	getVext();
	initnbr();
	rhoinit();
	
	for(int i=1;i<=100;++i)
		iterate();
	clock_t end = clock();
	double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
	
	for(int i=0;i<N;i++)
	{
		int i2=i/Nz;
		int j2=i%Nz;
		printf("%d %d %f\n",i2,j2,rho[i]);
		
	}
	printf("------------------\ntime: %f s\n",elapsed);
}
