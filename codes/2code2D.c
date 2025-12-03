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
#define dz 0.05
#define dx 0.05
#define Lx 10
#define Lz 10
#define mu 0.6045291013
#define rhob 0.304665

#define Nx 200
#define Nz 200
#define NiR 10
#define NiLJ 50//rc/dz
#define R2 0.25

#define N 40000 //Nx*Nz
#define NiR2 400 //4*NiR*NiR


/*
const int Nx = (int)(Lx/dx);
const int Nz = (int)(Lz/dz);
const int NiR = (int) (R/dz);
const double R2=R*R;
*/
//double rhob=.304665/sigma/sigma/sigma;
//double muex=mu-log(rhob);
const double alpha=0.01;
const double rmin=pow(2.,1./6);
const double rc=2.5;
const double eps=1.0;

double rho[Nx][Nz],rhonew[Nx][Nz],Vext[Nx][Nz];
double n0[Nx][Nz],n1[Nx][Nz],n2[Nx][Nz],n3[Nx][Nz],n1vx[Nx][Nz],n1vz[Nx][Nz],n2vx[Nx][Nz],n2vz[Nx][Nz];
double dphidn0[Nx][Nz],dphidn1[Nx][Nz],dphidn2[Nx][Nz],dphidn3[Nx][Nz],dphidn1vx[Nx][Nz],dphidn1vz[Nx][Nz],dphidn2vx[Nx][Nz],dphidn2vz[Nx][Nz];
double c1[Nx][Nz];

double dist(int,int,int,int);
double omega3(double);
double omega2(double);
double omega1(double);
double omega0(double);
double omega1vx(int,int,int,int);
double omega1vz(int,int,int,int);
double omega2vx(int,int,int,int);
double omega2vz(int,int,int,int);
void getn(), getc1_fmt(),getc1_LJ(),filterc1(), getVext(), rhoinit(), iterate();

double J11(double ,double , double),J5(double ,double , double),ULJ(double);

double dist(int i1,int j1,int i2,int j2)
{
	double Dx=fabs(i2-i1)*dx;
	double Dz=fabs(j2-j1)*dz;
	return Dx*Dx+Dz*Dz;
}


double omega3(double dsds)
{
	if(dsds>=R2)
	return 0.;
	return 2.*sqrt(R2-dsds);
}

double omega2(double dsds)
{
	if(dsds>=R2)
	return 0.;
	return 2.*R/sqrt(R2-dsds);
}

double omega1(double dsds)
{
	
	if(dsds>=R2)
	return 0.;
	return 1./2/PI/sqrt(R2-dsds);
}

double omega0(double dsds)
{
	
	if(dsds>=R2)
	return 0;
	
	return 1./2/PI/R/sqrt(R2-dsds);
}

double omega1vx(int i1,int j1,int i2, int j2)
{
	double Dx=(i1-i2);
	double Dz=(j1-j2);
	double dsds=Dx*Dx+Dz*Dz;
	if(dsds>=R2)
	return 0;
	double ds=sqrt(dsds);
	
	return Dx/2/PI/R/sqrt(R2-dsds);
}

double omega1vz(int i1,int j1,int i2, int j2)
{
	double Dx=(i1-i2);
	double Dz=(j1-j2);
	double dsds=Dx*Dx+Dz*Dz;
	if(dsds>=R2)
	return 0;
	double ds=sqrt(dsds);
	
	return Dz/2/PI/R/sqrt(R2-dsds);
}



double omega2vx(int i1,int j1,int i2, int j2)
{
	double Dx=(i1-i2);
	double Dz=(j1-j2);
	double dsds=Dx*Dx+Dz*Dz;
	if(dsds>=R2)
	return 0;
	
	return 2.*Dx/sqrt(R2-dsds);
}

double omega2vz(int i1,int j1,int i2, int j2)
{
	double Dx=(i1-i2);
	double Dz=(j1-j2);
	double dsds=Dx*Dx+Dz*Dz;
	if(dsds>=R2)
	return 0;
	
	return 2.*Dz/sqrt(R2-dsds);
}

double J11(double r, double a, double b)
{

    double tha = (fabs(r - a) < 1e-12 ? M_PI/2 : asin(r/a));
    double thb = asin(r/b);
    double integral =
        (1./512.)*(252.*(thb-tha)
        -105.*(sin(2.*thb)-sin(2.*tha))
        +52.5*(sin(4.*thb)-sin(4.*tha))
        -(28./3.)*(sin(6.*thb)-sin(6.*tha))
        +(7./8.)*(sin(8.*thb)-sin(8.*tha)));
    return pow(r,-11)*integral;
}

double J5(double r, double a, double b)
{

    double tha = (fabs(r - a) < 1e-12 ? M_PI/2 : asin(r/a));
    double thb = asin(r/b);
    double integral =
        (3./8.)*(thb-tha)
        -(1./4.)*(sin(2.*thb)-sin(2.*tha))
        +(1./32.)*(sin(4.*thb)-sin(4.*tha));
    return pow(r,-5)*integral;
}


double ULJ(double dsds)
{
	double r=sqrt(dsds);
	if(r>rc)
	return 0.;
	
	if(r>rmin)
	{
	double temp=8.*eps*(J11(r,r,rc)-J5(r,r,rc));
	return temp;
	}
	if(r>0)
	return -2.*eps*sqrt(rmin*rmin-dsds) + 8.*eps*(J11(r,rmin,rc)-J5(r,rmin,rc));
	
	return -2.*eps*rmin+8.*eps*((pow(rmin,-11.)-pow(rc,-11.))/11. - (pow(rmin,-5.)-pow(rc,-5.))/5. );
		
}

void getn()
{

	for(int i1=0;i1<Nx;i1++)
		for(int j1=0;j1<Nz;j1++)
			{
				n0[i1][j1]=0;
				n1[i1][j1]=0;
				n2[i1][j1]=0;
				n3[i1][j1]=0;
				n1vx[i1][j1]=0;
				n1vz[i1][j1]=0;
				n2vx[i1][j1]=0;
				n2vz[i1][j1]=0;
			}


	for(int i1=0;i1<Nx;i1++)
		for(int j1=0;j1<Nz;j1++)
			for(int i2=MAX(0,i1-NiR);i2<=MIN(i1+NiR,Nx-1);i2++)
				for(int j2=MAX(0,j1-NiR);j2<=MIN(j1+NiR,Nz-1);j2++)
					{
						double dsds=dist(i1,j1,i2,j2);
						
						n0[i1][j1]+=(rho[i2][j2]*omega0(dsds))*dx*dz;
						n1[i1][j1]+=(rho[i2][j2]*omega1(dsds))*dx*dz;
						n2[i1][j1]+=(rho[i2][j2]*omega2(dsds))*dx*dz;
						n3[i1][j1]+=(rho[i2][j2]*omega3(dsds))*dx*dz;
						n1vx[i1][j1]+=(rho[i2][j2]*omega1vx(i2,j2,i1,j1))*dx*dz;
						n1vz[i1][j1]+=(rho[i2][j2]*omega1vz(i2,j2,i1,j1))*dx*dz;
						n2vx[i1][j1]+=(rho[i2][j2]*omega2vx(i2,j2,i1,j1))*dx*dz;
						n2vz[i1][j1]+=(rho[i2][j2]*omega2vz(i2,j2,i1,j1))*dx*dz;
						
					}

}

void getc1_fmt()
{
	for(int i=0;i<Nx;i++)
		for(int j=0;j<Nz;j++)
			{
				dphidn0[i][j] = -log(1-n3[i][j]);
				dphidn1[i][j] =  n2[i][j]/(1-n3[i][j]);
				dphidn2[i][j] = n1[i][j]/(1-n3[i][j]) + (3*n2[i][j]*n2[i][j] - 3*(n2vx[i][j]*n2vx[i][j] + n2vz[i][j]*n2vz[i][j]))/(24*PI*(1-n3[i][j])*(1-n3[i][j]));
				dphidn3[i][j] = n0[i][j]/(1-n3[i][j]) + (n1[i][j]*n2[i][j] - n1vx[i][j]*n2vx[i][j] - n1vz[i][j]*n2vz[i][j])/(1-n3[i][j])/(1-n3[i][j]) + (n2[i][j]*n2[i][j]*n2[i][j] - 3*n2[i][j]*(n2vx[i][j]*n2vx[i][j] + n2vz[i][j]*n2vz[i][j]))/12/PI/(1-n3[i][j])/(1-n3[i][j])/(1-n3[i][j]);
				dphidn1vx[i][j] = -n2vx[i][j]/(1-n3[i][j]);
				dphidn1vz[i][j] = -n2vz[i][j]/(1-n3[i][j]);
				dphidn2vx[i][j] = -n1vx[i][j]/(1-n3[i][j]) - 3*n2[i][j]*n2vx[i][j]/(12*PI*(1-n3[i][j])*(1-n3[i][j]));
				dphidn2vz[i][j] = -n1vz[i][j]/(1-n3[i][j]) - 3*n2[i][j]*n2vz[i][j]/(12*PI*(1-n3[i][j])*(1-n3[i][j]));
				
			}
	
	for(int i1=0;i1<Nx;i1++)
		for(int j1=0;j1<Nz;j1++)
			{
				c1[i1][j1]=0;
			}
	
	
	for(int i1=0;i1<Nx;i1++)
		for(int j1=0;j1<Nz;j1++)
			for(int i2=MAX(0,i1-NiR);i2<=MIN(i1+NiR,Nx-1);i2++)
				for(int j2=MAX(0,j1-NiR);j2<=MIN(j1+NiR,Nz-1);j2++)
					{
					double dsds=dist(i1,j1,i2,j2);
					c1[i1][j1]+= -(dphidn0[i2][j2]*omega0(dsds) + dphidn1[i2][j2]*omega1(dsds) +dphidn2[i2][j2]*omega2(dsds) +dphidn3[i2][j2]*omega3(dsds)+ dphidn1vx[i2][j2]*omega1vx(i2,j2,i1,j1) + dphidn1vz[i2][j2]*omega1vz(i2,j2,i1,j1) + dphidn2vx[i2][j2]*omega2vx(i2,j2,i1,j1) + dphidn2vz[i2][j2]*omega2vz(i2,j2,i1,j1))*dx*dz;
					}
				
}

void getc1_LJ()
{
	
	for(int i1=0;i1<Nx;i1++)
		for(int j1=0;j1<Nz;j1++)
			for(int i2=MAX(0,i1-NiLJ);i2<=MIN(i1+NiLJ,Nx-1);i2++)
				for(int j2=MAX(0,j1-NiLJ);j2<=MIN(j1+NiLJ,Nz-1);j2++)
					{
					double dsds=dist(i1,j1,i2,j2);
					c1[i1][j1]+= -(rho[i2][j2]*ULJ(dsds))*dx*dz;
					
					}
}

void filterc1(int NiW)
{
	double temp=c1[Nx/2][Nz/2];
	for(int i=0;i<Nx;i++)
	{
		for(int j=Nz-(2*NiW+1);j<Nz;j++)
		{
			c1[i][j]=temp;
		}
	}
	
	
	for(int j=0;j<Nz;j++)
	{
		temp=c1[Nx/2][j];
		for(int i=0;i<=2*NiW+1;i++)
		{
			c1[i][j]=temp;
			c1[Nx-i][j]=temp;
		}
	}
}

void getVext(){
	for(int i=0;i<Nx;i++)
		for(int j=0;j<NiR;j++)
			Vext[i][j] = 1000.;
}

void rhoinit()
{
	for(int i=0;i<Nx;i++)
		for(int j=0;j<Nz;j++)
			rho[i][j]=rhob*exp(-Vext[i][j]);
}

void iterate(){	
	getn();
	getc1_fmt();
	//filterc1(NiR);
	getc1_LJ();
	//filterc1(NiLJ);
	
	double muex=-c1[Nx/2][Nz/2];//mu-log(rhob);
	for(int i=0;i<Nx;i++)
		for(int j=0;j<Nz;j++)
			rhonew[i][j]=rhob*exp(-Vext[i][j]+c1[i][j]+muex);
	
	for(int i=0;i<Nx;i++)
		for(int j=0;j<Nz;j++)
			rho[i][j]=alpha*rhonew[i][j] + (1-alpha)*rho[i][j];
	
	//printf("%f\n",rho[Nx/2][Nz-1]);
	//printf("%f %f %f %f \n",n0[Nx/2][Nz/2],n1[Nx/2][Nz/2],n2[Nx/2][Nz/2],n3[Nx/2][Nz/2]);
	double rhotemp=rho[Nx/2][Nz-1];
	for(int i=0;i<Nx;i++)
		for(int j=0;j<Nz;j++)
			rho[i][j]=rho[i][j]/rhotemp*rhob;
	
}




void main()
{
	clock_t start = clock();
	getVext();
	rhoinit();
	/*
	for(int i=1;i<=1;++i)
		iterate();
	clock_t end = clock();
	double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
	
	for(int i=0;i<Nx;i++)
		for(int j=0;j<Nz;j++)
			printf("%d %d %f\n",i,j,rho[i][j]);
		
	printf("------------------\ntime: %f s\n",elapsed);
	*/
	
	/*
	for(int i=0;i<1000;i++)
	{
		double ds=0.01*i;
		printf("%f %f\n",ds,ULJ(ds*ds));
	}*/
	printf("%f\n",J11(rmin,rmin,rc) );
}
