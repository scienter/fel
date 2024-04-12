#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"

double complex ***complexMemory3Asign(int harmony,int nz,int nx,int ny);
void doubleDeleteField3(double ***field,int harmony,int subNz);

void selfSeed_Field(Domain *D)
{
   int h,i,j,n,rank,N,startI,endI,numHarmony,minI,maxI,nn;
	int dataNum,start,sliceN,subN,*minmax;
	double delayT,dt,val,tmp,tau,ctau,arg,first,coef,sinTh,realV,imagV;
	double k0,shiftT,d,extincL,chi0,result,*sendData,*recvData;
	double complex compVal,J,***U;
   int myrank, nTasks;

   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	numHarmony=D->numHarmony;
   startI=1;  endI=D->subSliceN+1;
	minI=D->minI; maxI=D->maxI;
	sliceN=D->sliceN;
   N=D->nx*D->ny;

   dt = D->lambda0*D->numLambdaU/velocityC;
	delayT = D->chi_delay;
	sinTh = sin(D->bragTh*M_PI/180.0);
	d=D->chi_d;
	extincL=D->extincL;
	chi0=D->chi0;
	k0=D->ks;
   coef=M_PI*M_PI*D->chi_d/(sinTh*D->extincL*D->extincL);	


	minmax=D->minmax;
	dataNum=D->numHarmony*N*(minmax[1]-minmax[0])*2;
   sendData=(double *)malloc(dataNum*sizeof(double ));
   recvData=(double *)malloc(dataNum*sizeof(double ));

   U=complexMemory3Asign(D->numHarmony,D->subSliceN+2,D->nx,D->ny);

	for (rank=0; rank<nTasks; rank++) {

     if(myrank==rank) {
	    start=0;
       for(h=0; h<numHarmony; h++) 
         for(i=startI; i<endI; i++) 
           for(j=0; j<N; j++) {
             val=D->U[h][i][j];
             sendData[start]=creal(val);
             sendData[start+1]=cimag(val);
				 start += 2;
           }
	    
       for(n=rank+1; n<nTasks; n++) 
         MPI_Send(sendData,dataNum,MPI_DOUBLE,n,myrank,MPI_COMM_WORLD);
	  }
     for(n=rank+1; n<nTasks; n++) {
       if(myrank==n) {
         MPI_Recv(recvData,dataNum,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
			/*
			subN=minmax[rank+1]-minmax[rank];
         start=0;
			//for(h=0; h<numHarmony; h++) 
			for(h=0; h<1; h++) 
           for(i=0; i<subN; i++) 
             for(j=0; j<N; j++) {
               realV=recvData[start];
               imagV=recvData[start+1];
               compVal=realV + I*imagV;
					start+=2;

					shiftT=(minI+i-sliceN*0.5)*dt + delayT;
               for(nn=startI; nn<endI; nn++) {
	              tau = (nn-1+minI)*dt + shiftT;
				     ctau=velocityC*tau;

                 tmp=ctau*(2.0*d/sinTh+ctau/sinTh/sinTh);
		           if(tmp==0.0) J=0.5;
				     else if(tmp<0) {
					    arg=M_PI/extincL*sqrt(fabs(tmp));
						 J=gsl_sf_bessel_I1(arg); J/=arg;
						 J*=-I;
					  } else {
					    arg=M_PI/extincL*sqrt(tmp);
					    J=gsl_sf_bessel_J1(arg); J/=arg;
					  }

                 tmp=chi0*k0*(d+ctau/sinTh)/2.0/sinTh;
                 first=cexp(I*tmp);

                 result=coef*first*J*compVal;

	              U[h][nn][j]+=result*dt*velocityC;
               }
			    }
				*/
		 }  
	  }    //End for(receve n)

     MPI_Barrier(MPI_COMM_WORLD);
   }   // End for (rank)

   for(h=0; h<numHarmony; h++) 
     for(i=startI; i<endI; i++) 
       for(j=0; j<N; j++)
         D->U[h][i][j]=U[h][i][j];

   complexDeleteField3(U,D->numHarmony,D->subSliceN+2);
	free(sendData);
	free(recvData);
}



