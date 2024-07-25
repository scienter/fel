#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <mpi.h>
#include <complex.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf_bessel.h>


double complex **complexMemoryAsign(int harmony,int nx,int ny);
double complex ***complexMemory3Asign(int harmony,int nz,int nx,int ny);
double ***doubleMemory3Asign(int harmony,int nz,int nx,int ny);

void boundary(Domain *D)
{
   FILE *out;
   int i,j,h,s,n,b,totalCnt,rank,nx,ny,N,nn,nz,m;
   ptclList *New;
   LoadList *LL;
   UndulatorList *UL;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   double *number,minZ,n0,sumDouble,tmpDouble,dZ,posZ,xi;
   int sliceN,max,l,subCnt,remain,tmpN,minN,maxN;

   dZ=D->lambda0*D->numSlice;
   minZ=D->minZ;
   sliceN=D->sliceN;

   // finding minI, maxI, and minmax
   sumDouble=0;
   N=D->nx*D->ny*D->numHarmony;
   number=(double *)malloc(sliceN*sizeof(double ));
   D->minmax=(int *)malloc((nTasks+1)*sizeof(int ));
   D->minmax[0]=0;
   D->minmax[nTasks]=sliceN;
   for(i=0; i<sliceN; i++) number[i]=0;   
     
   LL=D->loadList;
   while(LL->next) {
      totalCnt=LL->numBeamlet*LL->numInBeamlet;
      for(i=0; i<sliceN; i++) {
         posZ=i*dZ+minZ;
         for(l=0; l<LL->znodes-1; l++) {
            if(posZ>=LL->zpoint[l] && posZ<LL->zpoint[l+1]) {
               n0=((LL->zn[l+1]-LL->zn[l])/(LL->zpoint[l+1]-LL->zpoint[l])*(posZ-LL->zpoint[l])+LL->zn[l]);
               tmpDouble=n0*totalCnt;
	            number[i]+=tmpDouble;
               sumDouble+=tmpDouble;
	         } else ;
  	      }
      }
      LL->totalCnt=totalCnt;
      LL=LL->next;
   }

   sumDouble+=sliceN*N*2;
   for(i=0; i<sliceN; i++) number[i]+=N*2;
   if(myrank==0) {
      tmpDouble=sumDouble/(1.0*nTasks);
      max=0; rank=1;
      sumDouble=0;
      for(i=0; i<sliceN; i++) {
         sumDouble+=number[i];
         if(sumDouble>tmpDouble) {
            sumDouble=0.0;
            D->minmax[rank]=i;
            rank++;
         } else ;
      }
   } else ;
   MPI_Bcast(D->minmax,nTasks+1,MPI_INT,0,MPI_COMM_WORLD);
   D->minI=D->minmax[myrank];
   D->maxI=D->minmax[myrank+1];
   D->subSliceN=D->maxI-D->minI;
	D->startI = 1;
	D->endI = D->subSliceN+1;
   free(number);
   MPI_Barrier(MPI_COMM_WORLD);
   printf("myrank=%d,minI=%d,maxI=%d,subSliceN=%d\n",myrank,D->minI,D->maxI,D->subSliceN);
  
   // Field memory setting
   D->U=complexMemory3Asign(D->numHarmony,D->subSliceN+2,D->nx,D->ny);
   D->Uc=complexMemory3Asign(D->numHarmony,D->subSliceN+2,D->nx,D->ny);
   D->ScU=complexMemory3Asign(D->numHarmony,D->subSliceN+2,D->nx,D->ny);
   D->totalEnergy=(double **)malloc(D->maxStep*sizeof(double *));
   D->totalBunch=(double **)malloc(D->maxStep*sizeof(double *));
   for(i=0; i<D->maxStep; i++) {
     D->totalEnergy[i]=(double *)malloc(D->numHarmony*sizeof(double ));
     D->totalBunch[i]=(double *)malloc(3*sizeof(double ));
   }

   // space charge
   nz = D->subSliceN+2;
   D->Ez=(double complex ****)malloc(nz*sizeof(double complex ***));
   for(i=0; i<nz; i++) {
     D->Ez[i]=(double complex ***)malloc(D->nr*sizeof(double complex **));
     for(j=0; j<D->nr; j++) {
       D->Ez[i][j]=(double complex **)malloc(D->SCLmode*sizeof(double complex *));
       for(l=0; l<D->SCLmode; l++)
         D->Ez[i][j][l]=(double complex *)malloc(D->SCFmode*sizeof(double complex ));
     }
   }
   for(i=0; i<nz; i++)
     for(j=0; j<D->nr; j++)
       for(l=0; l<D->SCLmode; l++)
         for(m=0; m<D->SCFmode; m++)
           D->Ez[i][j][l][m]=0.0+I*0.0;



   // Memory for shift
   D->shift=0.0;

   // setting up particle's pointer
   LL=D->loadList;
   s=0; 
   while(LL->next) {
     totalCnt=LL->numBeamlet*LL->numInBeamlet;
     LL->totalCnt=totalCnt;
     LL=LL->next;
     s++;
   }

   D->particle=(Particle *)malloc((D->subSliceN+2)*sizeof(Particle ));
   for(i=0; i<D->subSliceN+2; i++) {
     D->particle[i].head = (ptclHead **)malloc(sizeof(ptclHead *));
     for(s=0; s<D->nSpecies; s++) {
       D->particle[i].head[s] = (ptclHead *)malloc(sizeof(ptclHead ));
       D->particle[i].head[s]->pt=NULL;
     }
   }
   D->avePx=0.0;
   D->avePy=0.0;
   D->aveGam=D->gamma0;
   D->totalCnt=1.0;

   // initialize prev B field
   D->prevK=0.0;

   // twiss parameter
   D->twsBX = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsGX = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsAX = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsEmitX = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsBY = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsGY = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsAY = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsEmitY = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsG = (double *)malloc((D->maxStep+1)*sizeof(double ));

   // Bessel Table
   D->BesselMax = D->harmony[D->numHarmony-1]*2*M_PI;
   D->dBessel = D->BesselMax/(1.0*D->bn);
	D->BesselMaxOrder = (D->harmony[D->numHarmony-1]+1)/2+1;
   D->BesselJ = (double **)malloc(D->bn*sizeof(double *));
	for(i=0; i<D->bn; i++) 
     D->BesselJ[i] = (double *)malloc(D->BesselMaxOrder*sizeof(double ));
   for(i=0; i<D->bn; i++) {
     xi = D->dBessel*i;
     for(j=0; j<D->BesselMaxOrder; j++) 
	    D->BesselJ[i][j]=gsl_sf_bessel_Jn(j,xi);
   }

   // define Matrix Ma
   gsl_complex g;
   int stat;
   nx=D->nx;   ny=D->ny;

   // wake field
   D->den=(double *)malloc(sliceN*sizeof(double ));   
   D->wakeF=(double *)malloc(sliceN*sizeof(double ));   
   D->wakeE=(double *)malloc(sliceN*sizeof(double ));   
   for(i=0; i<sliceN; i++) {
     D->den[i]=0.0;
     D->wakeF[i]=0.0;
     D->wakeE[i]=0.0;
   }
}

double complex **complexMemoryAsign(int harmony,int nx,int ny)
{
   int h,i,N;
   double complex **field;

   N=nx*ny;
   field = (double complex **)malloc(harmony*sizeof(double complex *));
   for(h=0; h<harmony; h++) 
     field[h] = (double complex *)malloc(N*sizeof(double complex ));
   
   for(h=0; h<harmony; h++)
     for(i=0; i<N; i++)
       field[h][i]=0.0+I*0.0;

   return field;
}

double complex ***complexMemory3Asign(int harmony,int nz,int nx,int ny)
{
   int n,h,i,j,N;
   double complex ***field;

   N=nx*ny;
   field = (double complex ***)malloc(harmony*sizeof(double complex **));
   for(h=0; h<harmony; h++) {
      field[h] = (double complex **)malloc(nz*sizeof(double complex *));
      for(i=0; i<nz; i++) 
         field[h][i] = (double complex *)malloc(N*sizeof(double complex ));
   }
   
   for(h=0; h<harmony; h++)
     for(i=0; i<nz; i++)
       for(j=0; j<N; j++)
         field[h][i][j]=0.0+I*0.0;

   return field;
}

double ***doubleMemory3Asign(int harmony,int nz,int nx,int ny)
{
   int n,h,i,j,N;
   double ***field;

   N=nx*ny;
   field = (double ***)malloc(harmony*sizeof(double **));
   for(h=0; h<harmony; h++) {
      field[h] = (double **)malloc(nz*sizeof(double *));
      for(i=0; i<nz; i++) 
         field[h][i] = (double *)malloc(N*sizeof(double ));
   }
   
   for(h=0; h<harmony; h++)
     for(i=0; i<nz; i++)
       for(j=0; j<N; j++)
         field[h][i][j]=0.0;

   return field;
}
