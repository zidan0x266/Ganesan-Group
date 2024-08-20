#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <cstring>
#include <iostream>

#define maxnumtypes 100  //number of types of atoms that occur in dump files
#define maxnummassnum 100 //maximum number of massnumbers that could occur in scattlength.txt
#define numiglin 9  //number of lines for headers
#define maxnumatoms 1000000
#define pi 3.14159265
#define boxlmax 10000.0
#define maxnumq 1000000

int massta[maxnumtypes+1]={0};
int typa[maxnumatoms+1]={0};
int typearray[1000000]={0};  //array of atom types ""
double xa[maxnumatoms+1]={0};
double ya[maxnumatoms+1]={0};
double za[maxnumatoms+1]={0};
double sq[maxnumq+1]={0};
double q[maxnumq+1]={0};


//double scattfact[11]={0.0,scattc,scattc,scatth,scatth,scattc,scatth,scattc,scatth,scattc,scattc};
double scattfact[maxnumtypes+1]={0};
double scattfactmn[maxnummassnum+1]={0};

int main(int argc, char *argv[])
{
	FILE *fPin;
	FILE *fmass;
	FILE *fPin3;
	FILE *fPout;
	
	char lin[256];
	int g,i,j,k,m,n,p;      //counter variables
	int aind;	 //atom index
	int atyp;	 //atom type
	double ax,ay,az;  //atom coordinates
	int axs,ays,azs;    //atom "box shifts"
	int c;	      // dummy variable
	double mass=0.0;    //mass of atom
	int massint=0;
	int fstep=0;	//value of 1st timestep
	int sstep=0;	//value of 2nd timestep
	int lstep=0;	//value of last timestep
	int stepsize=0;     //value of time stepsize
	int numsamples=0;   //number of time samples
	long int numlinfile=0;   //number of lines in file
	int numatoms=0;     //number of atoms in system
	double xl,xh,yl,yh,zl,zh;  //limits of x,y,z in real coordinates
	int skipsamples;
	int numtypes=0;

	double boxlx,boxly,boxlz;
	double boxlxmin=boxlmax;
	double boxlymin=boxlmax;
	double boxlzmin=boxlmax;
	int numqvec,numqval;
	double minql,maxql,minq,maxq;
	int igsamples=0,usesamples=588;
	int qcount1,qcount2;
	double qval;
	double resum,imsum,totsumsq;
	double sumscattfact;
	double rand1,rand2,rand3;
	double qx,qy,qz,qr;
	double phi;
	double sintheta,costheta;
	int fatom=1,latom=420480;
	static char s1[200]={0};
	int numscattlength=0;
	int ijunk;
	int massnumber;
	double scattfactor;
	

	int bond_length=3,numq;
	int nx, ny, nz, n_max=100; 
	double q2;
	double q_width;
	float q_int;
	int Nint=0;
	double sq_tot[400]={0};
	long Nint_cnt[400]={0} ;
	double last_totsumsq;
	int N2;	
	int n_point=100000;
	int numatype;
	int maxn_point; 
	double q_value; 
	int a; 
	int N[maxnumq+1]={0};

	fPin=fopen(argv[1],"r");
	fPout=fopen(argv[2],"w");
	maxq=2*pi/bond_length; // edited part
	minq=2*pi/200.0; //(xh-xl);  
    printf ("Characters: %f %f\n",minq, maxq);   
	for(qcount2=1;qcount2<=n_point;qcount2++){
	 	resum=imsum=0;
		totsumsq=0;
		rand1=(double)rand()/(double)RAND_MAX;
		rand2=(double)rand()/(double)RAND_MAX;
		rand3=(double)rand()/(double)RAND_MAX;
			 
		nx =(n_max+1)*rand1;
		ny =((2*n_max)+1)*rand2-n_max;
		nz =((2*n_max)+1)*rand3-n_max;
			
		N2=(nx*nx+ny*ny+nz*nz);
	    if (N2<=(n_max*n_max)){
		    qx=nx*minq;
		    qy=ny*minq;
			qz=nz*minq;
		}
		else{continue;}
		q2=qx*qx+qy*qy+qz*qz;
		q[qcount2]=sqrt(q2);
        a = int(q[qcount2]/0.1);
        N[a] = N[a] +1;
		//printf ("Characters: %d, %f, %d\n",qcount2, sqrt(q2), a);
        }
    for(i=1;i<=400;i++){
		if(N[i]!=0){
			q_value =double(i) * 0.1;
 		      	printf("%lf\n",q_value);
		}
	}
	return 0;
}
