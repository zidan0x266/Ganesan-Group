//calculate densities from dump files
#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <cstring>
#include <iostream>

#define maxnumtypes 10 //number of types of atoms that occur in dump files
#define maxnummassnum 10 //maximum number of massnumbers that could occur in scattlength.txt
#define numiglin 9  //number of lines for headers
#define maxnumatoms 430000
#define pi 3.14159265
#define maxnumframes 5000
#define maxnumq 2000

int massta[maxnumtypes+1]={0};
int typa[maxnumatoms+1]={0};
int typearray[1000000]={0};  //array of atom types ""
static double xa[maxnumatoms+1][maxnumframes+1]={0};
static double ya[maxnumatoms+1][maxnumframes+1]={0};
static double za[maxnumatoms+1][maxnumframes+1]={0};
int tstep[maxnumframes+1]={0};
static double cossum[maxnumframes+1]={0};
static double sinsum[maxnumframes+1]={0};

//static double sqt[maxnumq+1][maxnumframes+1]={0};

//double scattfact[11]={0.0,scattc,scattc,scatth,scatth,scattc,scatth,scattc,scatth,scattc,scattc};
double scattfact[maxnumtypes+1]={0};
double scattfactmn[maxnummassnum+1]={0};

int main(int argc, char *argv[]){
	printf("Code for Sq(t) analysis of membranes\n");
	FILE *fPin,*fmass,*fscat;
	FILE *fq,*fPout;
	char lin[256];
	int g,i,j,k,l,m,n,p;      //counter variables
	int aind;	  //atom index
	int atyp;	  //atom type
	double ax,ay,az;  //atom coordinates
	int axs,ays,azs;    //atom "box shifts"
	int c;	      // dummy variable
	double mass=0.0;    //mass of atom
	int massint=0;
	int fstep=0;	//value of 1st timestep
	int sstep=0;	    //value of 2nd timestep
	int lstep=0;	//value of last timestep
	int stepsize=0;     //value of time stepsize
	int numsamples=0;   //number of time samples
	int numlinfile=0;   //number of lines in file
	int numatoms=0;     //number of atoms in system
	double xl,xh,yl,yh,zl,zh;  //limits of x,y,z in real coordinates
	int numtypes;
	double boxlx,boxly,boxlz;
	double minq;
	int igsamples=0,usesamples;
	int qcount1;
	double qval;
	double sumscattfact=0;
	double qx,qy,qz,qr,q2;
	double q[maxnumq]={0};
	int numscattlength=0;
	int massnumber;
	double scattfactor;
	int cscaled;
	int incimage;
	int nx, ny, nz, n_max=25;
	int a,q1;
	double sqtval;
	double sqt[maxnumq+1][maxnumframes+1]={0};
	double sqt1[maxnumq+1][maxnumframes+1]={0};
	double sqt2[maxnumq+1][maxnumframes+1]={0};
	int N[maxnumq]={0};
	int N1[maxnumq]={0};
	int N2[maxnumq]={0};
	double n_vec[1000][3]={0};
	int count1,count2,count3;
	count1=count2=count3=2;
	fPin=fopen(argv[1],"r");
	fPout=fopen(argv[2],"w");
	fq=fopen("qvec.txt","w");
	incimage = 1;
	cscaled = 1;
	while((c = getc(fPin)) != EOF){if(c == '\n'){break;}}  //reads to end of 1st line

	fscanf(fPin,"%d",&fstep);			       //reads first step in file from 2nd line
	while((c = getc(fPin)) != EOF){if(c == '\n'){break;}}  //reads to end of 2nd line

	while((c = getc(fPin)) != EOF){if(c == '\n'){break;}}  //reads to end of 3rd line

	fscanf(fPin,"%d",&numatoms);			       //reads number of atoms from 4th line
	while((c = getc(fPin)) != EOF){if(c == '\n'){break;}}  //reads to end of 4th line

	while((c = getc(fPin)) != EOF){if(c == '\n'){break;}}  //reads to end of 5th line

	fscanf(fPin,"%lf %lf",&xl,&xh);			       //reads x limits
	while((c = getc(fPin)) != EOF){if(c == '\n'){break;}}  //reads to end of 6th line

	fscanf(fPin,"%lf %lf",&yl,&yh);			       //reads y limits
	while((c = getc(fPin)) != EOF){if(c == '\n'){break;}}  //reads to end of 7th line

	fscanf(fPin,"%lf %lf",&zl,&zh);			       //reads z limits
	while((c = getc(fPin)) != EOF){if(c == '\n'){break;}}  //reads to end of 8th line
	numlinfile=8;
	while((c = getc(fPin)) != EOF){if(c == '\n'){numlinfile++;}}  //increments number of lines in file for each \n

	numsamples=numlinfile/(numatoms+numiglin); //finds number of time samples from above data
	usesamples=numsamples-igsamples;
	printf(" Number of samples are %d\n",numsamples);
	rewind(fPin);
	boxlx=xh-xl;
	boxly=yh-yl;
	boxlz=zh-zl;
	minq=2*pi/boxlx; // edited part
	printf("Box size and min q are %lf and %lf\n",boxlx,minq);
	printf("-------------------\n");
	printf("Generating q values\n");
	printf("-------------------\n");
	j=0;
	int remain;
	for(nx=-21;nx<=21;nx++){  //why 21?
		for(ny=-21;ny<=21;ny++){
			for(nz=-21;nz<=21;nz++){
				if((nx==0)&&(ny==0)&&(nz==0)){continue;}
				if(abs(nx)>6){count1=2;}else{count1=1;}
				if(abs(ny)>6){count2=2;}else{count2=1;}
				if(abs(nz)>6){count3=2;}else{count3=1;}
				qx = nx*minq;qy = ny*minq;qz = nz*minq;
				q2 = sqrt(qx*qx+qy*qy+qz*qz);
				if(q2>1.05){continue;}  //why 1.05? seems to be the largest q
				remain = ((int)(q2*1000))%100;  //why multiple by 1000 and then divide by 100
				if((remain >=90)||(remain <=10)){
					a=(int)(round(q2*10));
					if(N2[a]<40){
						j++;
						n_vec[j][1]=qx;
						n_vec[j][2]=qy;
						n_vec[j][3]=qz;
						q[j]=q2;
						fprintf(fq,"%d %d %d %d %lf\n",j,nx,ny,nz,q2);
						N2[a]++;
					}
				}
			}
		}
	}
	for(i=1;i<=10;i=i+1){
		fprintf(fq,"For q value of %lf, q vectors found are %d\n",((double)(i))/10,N[i]);
	}
	fclose(fq);
	int max_n=j;
	printf("Total q vectors generated are %d\n",max_n);
	fmass=fopen("mass.txt","rt");
	while(fgets(lin, 256, fmass) != NULL){numtypes++;}
	rewind(fmass);
	for(i=1;i<=numtypes;i++){
		fscanf(fmass,"%d %lf",&atyp,&mass);
		massint=(int)mass;
		if(mass-(double)massint>0.5){massint++;}
		massta[atyp]=massint;
	}
	printf("Read masses\n");

	fscat=fopen("scattlength.txt","r");
	while((c = getc(fscat)) != EOF){if(c == '\n'){numscattlength++;}}  //increments number of lines in file for each \n
	rewind(fscat);
	for(i=1;i<=numscattlength;i++){
		fscanf(fscat,"%d %d %lf\n",&j,&massnumber,&scattfactor);
		scattfactmn[massnumber]=scattfactor;
	}
	fclose(fscat);
	printf("Read scattering lengths\n");
	for(i=1;i<=numtypes;i++){
		scattfact[i]=scattfactmn[massta[i]];
	}
	for(k=1;k<=igsamples;k++){
		for(m=1;m<=numatoms+9;m++){
			while((c=getc(fPin))!='\n'){;}
		}
	}
	for(k=1;k<=usesamples;k++){
		fgets(lin, 256, fPin);
		fgets(lin, 256, fPin);
		tstep[k]=atoi(lin);
		fgets(lin, 256, fPin);
		fgets(lin, 256, fPin);
		fgets(lin, 256, fPin);
		fgets(lin, 256, fPin);
		fgets(lin, 256, fPin);
		fgets(lin, 256, fPin);
		fgets(lin, 256, fPin);
//		for(m=1;m<=9;m++){while((c=getc(fPin))!='\n'){;}}  //read to end of 1st-5th lines
		for(m=1;m<=numatoms;m++){
			if(incimage){fscanf(fPin,"%d %d %lf %lf %lf %d %d %d",&aind,&atyp,&ax,&ay,&az,&axs,&ays,&azs);}
			else{fscanf(fPin,"%d %d %lf %lf %lf",&aind,&atyp,&ax,&ay,&az);}
			if(cscaled){
				xa[aind][k]=xl+boxlx*(ax+axs);
				ya[aind][k]=yl+boxly*(ay+ays);
				za[aind][k]=zl+boxlz*(az+azs);
			}
			if(k==1){
				typa[aind]=atyp;
				sumscattfact+=scattfact[typa[aind]];  // what is it used for?
			}
		}
		while((c=getc(fPin))!='\n'){;}  //read to end of last atom line
		if(k==1){printf("sumscattfact %f\n",sumscattfact);}
	}
	printf("Read all coordinated\n");
	for(qcount1=1;qcount1<=max_n;qcount1++){
		qx=n_vec[qcount1][1];
		qy=n_vec[qcount1][2];
		qz=n_vec[qcount1][3];
		q2=qx*qx+qy*qy+qz*qz;
		q[qcount1]=sqrt(q2);
		for(i=1;i<=usesamples;i++){cossum[i]=sinsum[i]=0;}  //initialize sums
		for(i=1;i<=usesamples;i++){
			for(m=1;m<=numatoms;m++){
				ax=xa[m][i];
				ay=ya[m][i];
				az=za[m][i];
				qr=qx*ax+qy*ay+qz*az;
				cossum[i]+=(scattfact[typa[m]]*cos(qr));
				sinsum[i]+=(scattfact[typa[m]]*sin(qr));
			}
		}
//		printf(" %d qcount1, %lf q vector\n",qcount1, sqrt(q2));
		a = int(q[qcount1]*100);  // why multiple by 100?
		for(i=0;i<=usesamples-1;i++){  //dframe size
			for(j=1;j<=usesamples-i;j++){  //starting frame
				k=j+i;  //ending frame
				sqt[qcount1][i]+=(2.0*(cossum[j]*cossum[k]+sinsum[j]*sinsum[k]));
//				if(i==0){fprintf(fPout2,"framesize %d sframe %d eframe %d sqt[%d][%d] = %lf\n",i,j,k,qcount1,i,sqt[qcount1][i]);}
			}
		}
		N[a]++;
		for(i=0;i<=usesamples;i++){
			sqt[qcount1][i]/=(2.0*(usesamples-i));
			sqt1[a][i]+=sqt[qcount1][i];  // what are sqt1 and sqt2 used for?
		}  //normalize sqt by number of qvectors used and number of frame combinations for each time division - if division is of size i, each q vector will have 2.0*(useframe-i) - forward and backward combinations
		if((qcount1%100)==0){printf("%d q vectors out of %d qvectors done\n",qcount1,max_n);}
	}
	for(i=5;i<=295;i=i+10){
		for(j=0;j<10;j++){
			k=i+j;
			q1 = i+5;
			if(N[k]!=0){
				N1[q1]+=N[k];
				for(l=0;l<=usesamples-1;l++){
					sqt2[q1][l]+=sqt1[k][l];
				}
			}
		}
		for(l=0;l<=usesamples-1;l++){
			sqtval=sqt2[q1][l];
			fprintf(fPout,"%lf %d %d %lf\n",(((double)(q1))/100),(l*(tstep[2]-tstep[1])),N1[q1],sqtval/sqt2[q1][0]);
		}
	}
/*
	for(qcount1=1;qcount1<=numqval;qcount1++){
		//qval=2*pi*(1/minql-1/maxql)*(qcount1-1)/(numqval-1)+2*pi/maxql;
//		qval=minq*pow(maxq/minq,((double)qcount1-1)/((double)numqval-1));
		qval=q[qcount1];
		if(qval!=0){
			for(i=0;i<=usesamples-1;i++){
				sqtval=sqt[qcount1][i];
				fprintf(fPout,"%d %.12f %d %.12f %.12f %.12f %.12f %.12f\n",qcount1,qval,i,(double)(i*stepsize),sqtval,sqtval/(sumscattfact),sqtval/(sumscattfact*sumscattfact),sqtval/sqt[qcount1][0]);
			}
		}
	}
*/
	fcloseall;
	return 0;
}
