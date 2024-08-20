//
// Created by Zhang Zidan on 2019/12/25.
//
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <cstring>
#include <iostream>

#define nframes 5
#define pi 3.1415926
#define maxnumq 1000
#define nmax 42

int main(int argc, char *argv[]) {
    int i, j, k;
    int nx, ny, nz;
    double qx,qy,qz,qr,q2;
    int a,q1;
    double minq;
    int count1,count2,count3;
    count1=count2=count3=2;
    int N[maxnumq]={0};
    int N1[maxnumq]={0};
    int N2[maxnumq]={0};
    double q[maxnumq]={0};
    double n_vec[1000][3]={0};

    minq=2 * pi / 100.0;

    printf("Code for Sq(t) analysis of polyILs\n");
    for(i=0;i<=nframes-1;i++){
        for(j=1;j<=nframes-i;j++){
            k=j+i;
            printf("%d %d %d\n", i, j, k);
        }
    }

    j=0;
    int remain;
    for(nx=-nmax;nx<=nmax;nx++){
        for(ny=-nmax;ny<=nmax;ny++){
            for(nz=-nmax;nz<=nmax;nz++){
                if((nx==0)&&(ny==0)&&(nz==0)){continue;}
                if(abs(nx)>6){count1=2;}else{count1=1;}
                if(abs(ny)>6){count2=2;}else{count2=1;}
                if(abs(nz)>6){count3=2;}else{count3=1;}
                qx = nx*minq;qy = ny*minq;qz = nz*minq;
                q2 = sqrt(qx*qx+qy*qy+qz*qz);
                if(q2>2.0){continue;}
                remain = ((int)(q2*1000))%100;
                if((remain >=90)||(remain <=10)){
                    a=(int)(round(q2*10));
                    if(N2[a]<40){
                        j++;
                        n_vec[j][0]=qx;
                        n_vec[j][1]=qy;
                        n_vec[j][2]=qz;
                        q[j]=q2;
                        printf("%d %d %d %d %lf\n",j,nx,ny,nz,q2);
                        N2[a]++;
                    }
                }
            }
        }
    }
}
