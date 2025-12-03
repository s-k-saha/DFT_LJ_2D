#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>
#include<limits.h>


#define MIN(x,y) ((x<y)? x : y)
#define MAX(x,y) ((x>y)? x : y)

#define PI 3.14159265358979
#define R 0.5
#define sigma 1
#define dz 0.01
#define dx 0.01
#define mu 0.6045291013
#define rhob 0.304665

#define rc 2.5



#define Lx 20.0
#define Lz 20.0

#define Nx ((int)(Lx/dx))
#define Nz ((int)(Lz/dz))

#define K (Nz/2+1)//Nx/2+1 
#define NiR (R/dz)
#define NiLJ ((int)(rc/dx))//rc/dz
#define R2 R*R

#define N (Nx*Nz) //Nx*Nz
#define NiR2 (4*NiR*NiR) //4*NiR*NiR


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
//const double rc=2.5;
const double eps=1.0;
double ew;

double rho[Nx*Nz],rhonew[Nx*Nz],rhocopy[Nx*Nz],Vext[Nx*Nz],Ufilter[Nx*Nz];


double n0[Nx*Nz],n1[Nx*Nz],n2[Nx*Nz],n3[Nx*Nz],n1vx[Nx*Nz],n1vz[Nx*Nz],n2vx[Nx*Nz],n2vz[Nx*Nz];
double dphidn0[Nx*Nz],dphidn1[Nx*Nz],dphidn2[Nx*Nz],dphidn3[Nx*Nz],dphidn1vx[Nx*Nz],dphidn1vz[Nx*Nz],dphidn2vx[Nx*Nz],dphidn2vz[Nx*Nz];
double c1[Nx*Nz];
double c1_temp[Nx*Nz];

fftw_complex  omega3[Nx*K];
fftw_complex  omega2[Nx*K];
fftw_complex  omega1[Nx*K];
fftw_complex  omega0[Nx*K];
//double U[Nx][Nz];

fftw_complex omega1vx[Nx*Nz];
fftw_complex omega1vz[Nx*Nz];

fftw_complex omega2vx[Nx*Nz];
fftw_complex omega2vz[Nx*Nz];

double dist(int,int,int,int);

void getn(), getc1_fmt(),getc1_LJ(),filterc1(), getVext(), rhoinit(), iterate();
void getomega3(),getomega2(),getomega1(),getomega0(),getomega1v(),getomega2v(),rhocpy(),filterrho(),write_rho(double,int),conv_FFT2D_2(double*,fftw_complex *,double*);

double aux_J5(double, double);
double aux_J11(double, double);
double J5(double, double,double);
double J11(double, double,double);
double ULJ(double);

int getDx(int),getDz(int);
void add_c(double*,double*);
void conv_FFT2D(double*, double*,double*);



double dist(int i1,int j1,int i2,int j2)
{
	double Dx=dx*abs(i2-i1);
	double Dz=dz*abs(j2-j1);
	Dx=MIN(Dx,Lx-Dx);
	Dz=MIN(Dz,Lz-Dz);
	return Dx*Dx+Dz*Dz;
}

double aux_J11(double r, double a)
{
	double t1=(315./1280) * (pow(a,-2))*(pow(r,-11));
	double t2=(210./1280) * (pow(a,-4))*(pow(r,-9));
	double t3=(168./1280) * (pow(a,-6))*(pow(r,-7));
	double t4=(144./1280) * (pow(a,-8))*(pow(r,-5));
	double t5=(1./10) * (pow(a,-10))*(pow(r,-3));
	
	double temp=0.;
	if(a==r)
		temp=0.;
	else
		temp=atan(sqrt(a*a/r/r -1));
	
	
	return (r*sqrt(a*a-r*r)*(t1+t2+t3+t4+t5)) + (315./1280)*temp*(pow(r,-11));
}

double aux_J5(double r, double a)
{
	double temp=0;
	if(a==r)
		temp=0.;
	else	
		temp=atan(sqrt(a*a/r/r -1));
	
	double t1=(3./8)*(pow(a,-2))*(pow(r,-5));
	double t2=(1./4)*(pow(a,-4))*(pow(r,-3));
		
	return (r*sqrt(a*a-r*r)*(t1+t2)) + (3./8)*(pow(r,-5))*temp;
	
}


double J11(double r, double a, double b)
{
	if(r>=1.0)
	return aux_J11(r,b)-aux_J11(r,a);
	
	
	double t1=(1./11) * (-pow(b,-11)+pow(a,-11));
	double t2=(1./26) * (-pow(b,-13)+pow(a,-13)) * pow(r,2);
	double t3= (1./40) * (-pow(b,-15)+pow(a,-15)) * pow(r,4);
	double t4=	(5./272) * (-pow(b,-17)+pow(a,-17)) * pow(r,6);
	double t5= (35./2432) * (-pow(b,-19)+pow(a,-19)) * pow(r,8);
	double t6= (3./256)* (-pow(b,-21)+pow(a,-21)) * pow(r,10);
	double t7=(231./23552)* (-pow(b,-23)+pow(a,-23)) * pow(r,12);
	double t8= (429./51200)* (-pow(b,-25)+pow(a,-25)) * pow(r,14);
  
  return (t1+t2+t3+t4+t5+t6+t7+t8);
}

double J5(double r, double a, double b)
{
	if(r>=1.0)
		return aux_J5(r,b)-aux_J5(r,a);
		
	double t1=(1./5) * (-pow(b,-5)+pow(a,-5));
	double t2=(1./14) * (-pow(b,-7)+pow(a,-7)) * pow(r,2);
	double t3= (1./24) * (-pow(b,-9)+pow(a,-9)) * pow(r,4);
	double t4=	(5./176) *  (-pow(b,-11)+pow(a,-11)) * pow(r,6);
	double t5= (35./1664) *  (-pow(b,-13)+pow(a,-13)) * pow(r,8);
	double t6= (21./1280)* (-pow(b,-15)+pow(a,-15)) * pow(r,10);
	double t7= (231./17408)* (-pow(b,-17)+pow(a,-17)) * pow(r,12);
	double t8= (429./38912)* (-pow(b,-19)+pow(a,-19)) * pow(r,14);
	
	return (t1+t2+t3+t4+t5+t6+t7+t8);
    
}


double ULJ(double dsds)
{
	double r=sqrt(dsds);
	if(r>=rc)
	return 0.;
	
	if(r>rmin)
	{
	double temp=8.*eps*(J11(r,r,rc)-J5(r,r,rc));
	return temp;
	}
	return -2.*eps*sqrt(rmin*rmin-dsds) + 8.*eps*(J11(r,rmin,rc)-J5(r,rmin,rc));
	
	//return -2.*eps*rmin+8.*eps*((pow(rmin,-11.)-pow(rc,-11.))/11. - (pow(rmin,-5.)-pow(rc,-5.))/5. );
		
}

void getUfilter()
{
	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Nz;j++)
		{
			double d2=dist(i,j,0,0);
			if(d2>=rc*rc)
			Ufilter[i*Nz+j]=0.0;
			else
			{
				Ufilter[i*Nz+j]=ULJ(d2);
				Ufilter[i*Nz+j]=(Ufilter[i*Nz+j]>100.)?0.0:Ufilter[i*Nz+j];
			}
		}
	}
}

void getomega3()
{
	double kx=0.0,kz=0.0;
  double dkx=2.*PI/Lx,dkz=2.*PI/Lz;
  double k;
  for(int i=0;i<Nx;i++)
  {
  	if(i<=Nx/2)
  		kx=dkx*i;
  	else
  		kx=-dkx*(Nx-i);
  	for(int j=0;j<K;j++)
  	{
  		kz=dkz*j;
			k=sqrt(kx*kx+kz*kz);
			if(i==0 && j==0)
			{
				omega3[i*K+j][0]=4.*PI*R*R*R/3;
				omega3[i*K+j][1]=0.0;
			}
			else
			{
				omega3[i*K+j][0]=8./3*PI*R*R*(j1(R*k)/(k));
				omega3[i*K+j][1]=0.0;
  		}
  	}
  }
	
}

void getomega2()
{
  double kx=0.0,kz=0.0;
  double dkx=2.*PI/Lx,dkz=2.*PI/Lz;
  double k;
  for(int i=0;i<Nx;i++)
  {
  	if(i<=Nx/2)
  		kx=dkx*i;
  	else
  		kx=-dkx*(Nx-i);
  	for(int j=0;j<K;j++)
  	{
  		kz=dkz*j;
			k=sqrt(kx*kx+kz*kz);
			if(i==0 && j==0)
			{
				omega2[i*K+j][0]=4.*PI*R*R;
				omega2[i*K+j][1]=0.0;
			}
			else
			{
				omega2[i*K+j][0]=8.*PI*R*R*(j1(R*k)/(R*k));
				omega2[i*K+j][1]=0.0;
  		}
  	}
  }
}

void getomega1()
{
		for(int i=0;i<Nx;i++)
		for(int j=0;j<K;j++)
		{
		omega1[i*K+j][0]=omega2[i*Nz+j][0]*(1./4./PI/R);
		omega1[i*K+j][1]=omega2[i*Nz+j][1]*(1./4./PI/R);
		}
}
void getomega0()
{
		for(int i=0;i<Nx;i++)
		for(int j=0;j<K;j++)
		{
		omega0[i*K+j][0]=omega2[i*Nz+j][0]*(1./4./PI/R/R);
		omega0[i*K+j][1]=omega2[i*Nz+j][1]*(1./4./PI/R/R);
		}
}

void getomega2v()
{
	double kx=0.0,kz=0.0;
	double dkx=2.*PI/Lx,dkz=2.*PI/Lz;
	double k;
	for(int i=0;i<Nx;i++)
	{
		if(i<=Nx/2)
  		kx=dkx*i;
  	else
  		kx=-dkx*(Nx-i);
  	for(int j=0;j<K;j++)
  	{
  		kz=dkz*j;
			k=sqrt(kx*kx+kz*kz);
		
			omega2vx[i*K+j][1]=-kx*omega3[i*K+j][0];
			omega2vx[i*K+j][0]=0.0;
			omega2vz[i*K+j][1]=-kz*omega3[i*K+j][0];
			omega2vz[i*K+j][0]=0.0;	
		}
	
	}
}


void getomega1v()
{
	for(int i=0;i<Nx;i++)
	for(int j=0;j<K;j++)
	{
		omega1vx[i*K+j][0]=omega2vx[i*K+j][0]/(4.*PI*R);
		omega1vx[i*K+j][1]=omega2vx[i*K+j][1]/(4.*PI*R);
		
		omega1vz[i*K+j][0]=omega2vz[i*Nz+j][0]/(4.*PI*R);
		omega1vz[i*K+j][1]=omega2vz[i*Nz+j][1]/(4.*PI*R);
		
	}
}



void rhoinit()
{
	for(int i=0;i<Nx;i++)
		for(int j=0;j<Nz;j++)
			rho[i*Nz+j]=rhob*exp(-Vext[i*Nz+j]);
	
	for(int i=0;i<Nx;i++) //Setting walls left 
		for(int j=0;j<NiR;j++)
			rho[i*Nz+j]=0.0;
	
	/*
	for(int i=0;i<Nx;i++)
		for(int j=Nz-NiR+1;j<Nz;j++) //Setting walls right
			rho[i*Nz+j]=0.0;
	*/
}

void getn()
{

	for(int i1=0;i1<Nx;i1++)
		for(int j1=0;j1<Nz;j1++)
			{
				n0[i1*Nz+j1]=0;
				n1[i1*Nz+j1]=0;
				n2[i1*Nz+j1]=0;
				n3[i1*Nz+j1]=0;
				n1vx[i1*Nz+j1]=0;
				n1vz[i1*Nz+j1]=0;
				n2vx[i1*Nz+j1]=0;
				n2vz[i1*Nz+j1]=0;
			}
	
	conv_FFT2D_2(rhocopy, omega0, n0);
	conv_FFT2D_2(rhocopy, omega1, n1);
	conv_FFT2D_2(rhocopy, omega2,n2);
	conv_FFT2D_2(rhocopy, omega3, n3);
	conv_FFT2D_2(rhocopy, omega1vx, n1vx);
	conv_FFT2D_2(rhocopy, omega1vz, n1vz);
	conv_FFT2D_2(rhocopy, omega2vx, n2vx);
	conv_FFT2D_2(rhocopy, omega2vz, n2vz);
	
}



void conv_FFT2D(double *f, double *g, double *h2)
{
	fftw_complex *fft_f=fftw_alloc_complex((size_t)Nx * K);
	fftw_complex *fft_g=fftw_alloc_complex((size_t)Nx * K);
	fftw_complex *fft_h=fftw_alloc_complex((size_t)Nx * K);
	
	
	fftw_plan ff = fftw_plan_dft_r2c_2d(Nx,Nz, f, fft_f, FFTW_ESTIMATE);
  fftw_execute(ff);
	
	fftw_plan fg = fftw_plan_dft_r2c_2d(Nx,Nz, g, fft_g, FFTW_ESTIMATE);
  fftw_execute(fg);
  
  
  for (size_t i = 0; i < (size_t)Nx * K; ++i) {
        double a = fft_f[i][0], b = fft_f[i][1];
        double c = fft_g[i][0], d = fft_g[i][1]; // d==0
        fft_h[i][0] = a*c - b*d;
        fft_h[i][1] = a*d + b*c;
    }
  
  
	fftw_plan fh= fftw_plan_dft_c2r_2d(Nx,Nz,fft_h,h2,FFTW_ESTIMATE);
	fftw_execute(fh);
	
	
	fftw_destroy_plan(ff);
	fftw_destroy_plan(fg);
	fftw_destroy_plan(fh);
	
	
	for(int i=0;i<Nx*Nz;i++)
	h2[i]=(h2[i]/N)*dx*dz;
	
	fftw_free(fft_f);
  fftw_free(fft_g);
  fftw_free(fft_h);
}

void conv_FFT2D_2(double *f, fftw_complex *fft_g,double *h2)
{
	fftw_complex *fft_f=fftw_alloc_complex((size_t)Nx * K);
	fftw_complex *fft_h=fftw_alloc_complex((size_t)Nx * K);
	
	
	fftw_plan ff = fftw_plan_dft_r2c_2d(Nx,Nz, f, fft_f, FFTW_ESTIMATE);
  fftw_execute(ff);
	
  
  
  for (int i = 0; i < Nx * K; ++i) {
        double a = fft_f[i][0], b = fft_f[i][1];
        double c = fft_g[i][0], d = fft_g[i][1]; // d==0
        fft_h[i][0] = a*c - b*d;
        fft_h[i][1] = a*d + b*c;
    }
  
  
	fftw_plan fh= fftw_plan_dft_c2r_2d(Nx,Nz,fft_h,h2,FFTW_ESTIMATE);
	fftw_execute(fh);
	
	
	fftw_destroy_plan(ff);
	fftw_destroy_plan(fh);
	
	
	for(int i=0;i<Nx*Nz;i++)
	h2[i]=(h2[i]/N);
	
	fftw_free(fft_f);
  fftw_free(fft_h);
}


void add_c(double *c1, double *c1_temp)
{
	for(int i=0;i<Nx*Nz;i++)
	{
		c1[i]-=c1_temp[i];
	}
}


void getc1_fmt()
{
	for(int i=0;i<Nx;i++)
		for(int j=0;j<Nz;j++)
			{
				dphidn0[i*Nz+j] = -log(1-n3[i*Nz+j]);
				dphidn1[i*Nz+j] =  n2[i*Nz+j]/(1-n3[i*Nz+j]);
				dphidn2[i*Nz+j] = n1[i*Nz+j]/(1-n3[i*Nz+j]) + (3*n2[i*Nz+j]*n2[i*Nz+j] - 3*(n2vx[i*Nz+j]*n2vx[i*Nz+j] + n2vz[i*Nz+j]*n2vz[i*Nz+j]))/(24*PI*(1-n3[i*Nz+j])*(1-n3[i*Nz+j]));
				dphidn3[i*Nz+j] = n0[i*Nz+j]/(1-n3[i*Nz+j]) + (n1[i*Nz+j]*n2[i*Nz+j] - n1vx[i*Nz+j]*n2vx[i*Nz+j] - n1vz[i*Nz+j]*n2vz[i*Nz+j])/(1-n3[i*Nz+j])/(1-n3[i*Nz+j]) + (n2[i*Nz+j]*n2[i*Nz+j]*n2[i*Nz+j] - 3*n2[i*Nz+j]*(n2vx[i*Nz+j]*n2vx[i*Nz+j] + n2vz[i*Nz+j]*n2vz[i*Nz+j]))/12/PI/(1-n3[i*Nz+j])/(1-n3[i*Nz+j])/(1-n3[i*Nz+j]);
				dphidn1vx[i*Nz+j] = -n2vx[i*Nz+j]/(1-n3[i*Nz+j]);
				dphidn1vz[i*Nz+j] = -n2vz[i*Nz+j]/(1-n3[i*Nz+j]);
				dphidn2vx[i*Nz+j] = -n1vx[i*Nz+j]/(1-n3[i*Nz+j]) - 3*n2[i*Nz+j]*n2vx[i*Nz+j]/(12*PI*(1-n3[i*Nz+j])*(1-n3[i*Nz+j]));
				dphidn2vz[i*Nz+j] = -n1vz[i*Nz+j]/(1-n3[i*Nz+j]) - 3*n2[i*Nz+j]*n2vz[i*Nz+j]/(12*PI*(1-n3[i*Nz+j])*(1-n3[i*Nz+j]));
				
			}
	
	for(int i1=0;i1<Nx;i1++)
		for(int j1=0;j1<Nz;j1++)
			{
				c1[i1*Nz+j1]=0.;
				//c1_temp[i*Nz+j1]=0;
			}
	
	conv_FFT2D_2(dphidn0, omega0, c1_temp);
	add_c(c1,c1_temp);
	conv_FFT2D_2(dphidn1, omega1, c1_temp);
	add_c(c1,c1_temp);
	conv_FFT2D_2(dphidn2, omega2, c1_temp);
	add_c(c1,c1_temp);
	conv_FFT2D_2(dphidn3, omega3, c1_temp);
	add_c(c1,c1_temp);
	conv_FFT2D_2(dphidn1vx, omega1vx, c1_temp);
	add_c(c1,c1_temp);
	conv_FFT2D_2(dphidn1vz, omega1vz, c1_temp);
	add_c(c1,c1_temp);
	conv_FFT2D_2(dphidn2vx, omega2vx, c1_temp);
	add_c(c1,c1_temp);
	conv_FFT2D_2(dphidn2vz, omega2vz, c1_temp);
	add_c(c1,c1_temp);
	
				
}

void getc1_LJ()
{
	for(int i=0;i<Nx;i++)
		for(int j=Nz-NiLJ+1;j<Nz;j++)
			rhocopy[i*Nz+j]=0.0;
			
	conv_FFT2D(rhocopy,Ufilter, c1_temp);
	add_c(c1,c1_temp);
}


void getVext(){
	for(int i=0;i<Nx;i++)
		for(int j=0;j<NiR-1;j++)
			Vext[i*Nz+j] = 1000.;
	
	
	for(int i=0;i<Nx;i++)
	for(int j=NiR-1;j<Nz;j++)
			Vext[i*Nz+j] = eps*ew*(2./15*pow(dz*j,-9)-pow(dz*j,-3));
	
}

void rhocpy()
{
	for(int i=0;i<N;i++)
		rhocopy[i]=rho[i];
	
	for(int i=0;i<Nx;i++)
		for(int j=Nz-NiR+1;j<Nz;j++)
			rhocopy[i*Nz+j]=0.0;
}
void filterrho()
{
	for(int i=0;i<Nx;i++)
		for(int j=Nz-4*NiR;j<Nz;j++)
			rho[i*Nz+j]=rhob;
}

void iterate(){	
	rhocpy();
	getn();
	getc1_fmt();
	//filterc1(NiR);
	//getc1_LJ();
	//filterc1(NiLJ);

	double muex=mu-log(rhob);//-c1[(Nx/2 * Nz ) + (Nz/2)];//mu-log(rhob);
	//printf("%f\n",muex);
	/*
	for(int i=0;i<N;i++)
	c1[i]=c1[i]/(muex)*(-3.2872526167196);
	muex=-3.2872526167196;
	*/
	
	//printf("%f\n",muex);
	for(int i=0;i<Nx;i++)
		for(int j=0;j<Nz;j++)
			rhonew[i*Nz+j]=rhob*exp(-Vext[i*Nz+j]+c1[i*Nz+j]+muex);
	
	for(int i=0;i<Nx;i++)
		for(int j=0;j<Nz;j++)
			rho[i*Nz+j]=alpha*rhonew[i*Nz+j] + (1-alpha)*rho[i*Nz+j];
	
	
	filterrho();
	
	/*
	double rhotemp=rho[(Nx/2)*Nx+(Nz-1)];
	for(int i=0;i<Nx;i++)
		for(int j=0;j<Nz;j++)
			rho[i][j]=rho[i][j]/rhotemp*rhob;
	
	*/
}


void write_rho(double elapsed,int count_iter)
{
	char fname[100];
	sprintf(fname,"../data/rho2Deps%few%fdx%frhob%f.dat",eps,ew,dx,rhob);
	FILE *F=fopen(fname,"w");
	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Nz;j++)
			fprintf(F,"%d %d %f\n",i,j,rho[i*Nz+j]);
		fprintf(F,"\n");
	}
	fprintf(F,"------------------\n100 x %d cycles time: %f s\n",count_iter,elapsed);
	fclose(F);
	
	sprintf(fname,"../data/rho1Deps%few%fdx%frhob%f.dat",eps,ew,dx,rhob);
	F=fopen(fname,"w");
	
	for(int j=0;j<Nz;j++)
	fprintf(F,"%f %f\n",j*dz,rho[(Nx/2)*Nz+j]);
	fprintf(F,"------------------\n100 x %d cycles time: %f s\n",count_iter,elapsed);
	fclose(F);
}

void main(int argc,char *argv[])
{
	
	ew=0.0;//atof(argv[1]);
	
	clock_t start = clock();clock_t end;
	double elapsed;
	
	getVext();
	rhoinit();
	getUfilter();
	getomega3();
	getomega2();
	
	getomega1();
	getomega0();
	getomega2v();
	getomega1v();
	
	//iterate();	
	
	
	for(int i=1;i<=INT_MAX;++i)
	{	
		for(int j=1;j<=100;j++)
			iterate();
		end = clock();
		elapsed = (double)(end - start) / CLOCKS_PER_SEC;
		write_rho(elapsed,i);
	}
	

	/*
	int i=(Nx/2);
	int j=(Nz/2);
	printf("%d %d %f(1.0) %f(0.5) %f(3.14) %f(0.5236) %f(0) %f(0) \n",i,j,n0[i*Nz+j],n1[i*Nz+j],n2[i*Nz+j],n3[i*Nz+j],n2vx[i*Nz+j],n2vz[i*Nz+j]);
	*/
	
	/*
	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Nz;j++)
			printf("%d %d %f\n",i,j,c1[i*Nz+j]);
		printf("\n");
	}*/
	
	/*
	for(int i=0;i<2000;i++)
	{
		double ds=0.01*i;
		printf("%f %f\n",ds,ULJ(ds*ds));
	}*/
	
	//printf("%f %f\n",aux_J5(1.16,1.16),aux_J5(1.16,rc));
}
