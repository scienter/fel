#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"
#include <gsl/gsl_sf_bessel.h>

void saveBFactor(Domain *D,int iteration)
{
   int i,j,s,startI,endI,minI,sliceI,sliceN,numSlice,rank,cnt;
	double bucketZ,z,theta,*data,*recvData;
	double complex sum;
   LoadList *LL;
	ptclList *p;
   char name[100];
   FILE *out;

   int myrank, nTasks;    
   MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   startI=1;  endI=D->subSliceN+1;
	minI=D->minI;

	numSlice=D->sliceN+2;
	recvData=(double *)malloc(numSlice*sizeof(double ));
	data=(double *)malloc(numSlice*sizeof(double ));
	for(i=0; i<numSlice; i++) data[i]=0.0;
	
   for(i=startI; i<endI; i++)
   {
	  sliceI=i+minI;
	  cnt=0;
	  sum=0.0+I*0.0;
     for(s=0; s<D->nSpecies; s++)  {
       p=D->particle[i].head[s]->pt;
       while(p) {
         theta=p->theta;
			sum+=cexp(I*theta);
			cnt++;
			p=p->next;
		 }
	  }
     if(cnt>0) data[sliceI]=cabs(sum)/(cnt*1.0);
	}

   if(myrank==0)  {
     for(rank=1; rank<nTasks; rank++) {
	    MPI_Recv(recvData,numSlice,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
       for(i=0; i<numSlice; i++)
         data[i]+=recvData[i];
     }
	} else {
	  for(rank=1; rank<nTasks; rank++) {
	    if(myrank==rank)
		   MPI_Send(data,numSlice,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
     }
   }

   bucketZ=D->numSlice*D->lambda0;
   if(myrank==0) {
     sprintf(name,"bFact%d",iteration);
     out = fopen(name,"w");

     for(i=1; i<numSlice+1; i++) {
       z=(i-numSlice*0.5)*bucketZ;
       fprintf(out,"%g %g\n",z,data[i]);
     }
	  fclose(out);
     printf("%s is made.\n",name);
   } else ;

}

