#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <complex.h>
#include "mesh.h"
#include "constants.h"

void updateTotalEnergy(Domain *D,int iteration)
{
   double total,amp;
   int maxStep,i,j,h,N,numHarmony,startI,endI,step,start;
   double *recvData,*sendData,z,dz,minZ;
   double area,coef,coef2,dt;
   FILE *out;
   int myrank, nTasks;
   MPI_Status status; 

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   numHarmony=D->numHarmony;
   startI=1;  endI=D->subSliceN+1;
   dt=D->numSlice*D->lambda0/velocityC;
   N=D->nx*D->ny;
   dz=D->dz;
   minZ=D->minZ;

   for(h=0; h<numHarmony; h++) {
     total=0.0;
     for(i=startI; i<endI; i++) {
       for(j=0; j<N; j++) {
         amp=cabs(D->U[h][i][j]);
	      total+=amp*amp;
       }
     }
     D->totalEnergy[iteration][h]=total;
   }
   MPI_Barrier(MPI_COMM_WORLD);

   area=0.5*M_PI*D->spotSigR*D->spotSigR;
//   area=2.0*M_PI*LL->sigX*LL->sigY;
   coef=eMass*velocityC*velocityC*D->ks/eCharge;
   if(D->dimension==1)
     coef2=coef*coef/Z0*area;
   else 
     coef2=coef*coef/Z0*D->dx*D->dy;
   if(D->mode==Time_Dependent) coef2*=dt; else ;

   sendData=(double *)malloc(numHarmony*sizeof(double ));
   recvData=(double *)malloc(numHarmony*sizeof(double ));
   for(h=0; h<numHarmony; h++) 
      sendData[h] = D->totalEnergy[iteration][h];

   for(i=1; i<nTasks; i++) {
     if(myrank==i)  {
       MPI_Send(sendData,numHarmony,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
     }  else ;
   }

   if(myrank==0) {
     for(i=1; i<nTasks; i++) {
       MPI_Recv(recvData,numHarmony,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
       for(h=0; h<numHarmony; h++) 
	       D->totalEnergy[iteration][h] += recvData[h];
     }

     out=fopen("totalEnergy","a+");
     z=iteration*dz+minZ;
     fprintf(out,"%.15g",z);
     for(h=0; h<numHarmony; h++) 
        fprintf(out," %g",D->totalEnergy[iteration][h]*coef2);
     fprintf(out,"\n");
     fclose(out);
   } else ;

   free(sendData);
   free(recvData);
	
}

void saveTotalEnergy(Domain *D)
{
   int maxStep,i,h,N,numHarmony,startI,endI,step,start;
   double *recvData,*sendData,z,dz,minZ;
   double area,coef,coef2,dt;
   LoadList *LL;
   FILE *out;
   int myrank, nTasks;
   MPI_Status status; 

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   N=D->nx*D->ny;
   numHarmony=D->numHarmony;
   startI=1;  endI=D->subSliceN+1;
   maxStep=D->maxStep;
   dz=D->dz;
   minZ=D->minZ;
   dt=D->numSlice*D->lambda0/velocityC;

   LL=D->loadList;
   area=0.5*M_PI*D->spotSigR*D->spotSigR;
//   area=2.0*M_PI*LL->sigX*LL->sigY;
   coef=eMass*velocityC*velocityC*D->ks/eCharge;
   if(D->dimension==1)
     coef2=coef*coef/Z0*area;
   else 
     coef2=coef*coef/Z0*D->dx*D->dy;
   if(D->mode==Time_Dependent) coef2*=dt; else ;

   sendData=(double *)malloc(maxStep*numHarmony*sizeof(double ));
	start=0;
   for(step=0; step<maxStep; step++) {
     for(h=0; h<numHarmony; h++) 
       sendData[start+h] = D->totalEnergy[step][h];
     start += numHarmony;   
   }

   recvData=(double *)malloc(maxStep*numHarmony*sizeof(double ));
   for(i=1; i<nTasks; i++) {
     if(myrank==i)  {
       MPI_Send(sendData,maxStep*numHarmony,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
     }  else ;
   }

   if(myrank==0) {
     for(i=1; i<nTasks; i++) {
       MPI_Recv(recvData,maxStep*numHarmony,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
		 start=0;
       for(step=0; step<maxStep; step++) {
         for(h=0; h<numHarmony; h++) 
	        D->totalEnergy[step][h] += recvData[start+h];
         start += numHarmony;   
       }
     }

     out=fopen("totalEnergy","w");
     for(step=0; step<maxStep; step++) {
       z=step*dz+minZ;
       fprintf(out,"%.15g",z);
       for(h=0; h<numHarmony; h++) {
         fprintf(out," %g",D->totalEnergy[step][h]*coef2);
//			printf("step=%d, totalEnergy=%g\n",step,D->totalEnergy[step][h]*coef2);
       }
       fprintf(out,"\n");
     }
     fclose(out);
     printf("totalEnergy is made.\n");
   } else ;

   free(sendData);
   free(recvData);

   MPI_Barrier(MPI_COMM_WORLD);
}

