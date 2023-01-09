#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include <gsl/gsl_sf_bessel.h>
#include <complex.h>
#include <string.h>

void shiftField_1D(Domain *D,int iteration);
void shiftField_3D(Domain *D,int iteration);
void MPI_Transfer1F_Zminus(double complex ***f1,int harmony,int N,int fromI,int toI);
void MPI_Transfer1F_Zplus(double complex ***f1,int harmony,int N,int fromI,int toI);

void shiftField(Domain D,int iteration)
{
  int startI,endI,N;
  int myrank, nTasks;

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  N=D.nx*D.ny;
  startI=1; endI=D.subSliceN+1;
  switch(D.dimension) {

  //1D field
  case 1:
    //  shift minus
    MPI_Transfer1F_Zminus(D.U,D.numHarmony,N,startI,endI);

    shiftField_1D(&D,iteration);

    //  shift plus
    MPI_Transfer1F_Zplus(D.U,D.numHarmony,N,endI,1);
    break;
  case 3:
    //  shift minus
    MPI_Transfer1F_Zminus(D.U,D.numHarmony,N,startI,endI);

    shiftField_3D(&D,iteration);

    //  shift plus
    MPI_Transfer1F_Zplus(D.U,D.numHarmony,N,endI,1);
    break;

  default:
    printf("In EzSolve.c, what dimension?\n");
  }
}

void shiftField_3D(Domain *D,int iteration)
{
   double coef,aveGam,cnt,shiftZ,aveGam2,refGam2,avePx2,avePy2;
   int s,h,numHarmony,i,j,startI,endI,N;
   ptclList *p;

   numHarmony=D->numHarmony;
   startI=1;  endI=D->subSliceN+1;
   N=D->nx*D->ny;

   aveGam2=D->aveGam*D->aveGam;
   refGam2=D->gamma0*D->gamma0;
   avePx2=D->avePx*D->avePx;
   avePy2=D->avePy*D->avePy;

   shiftZ=((1+D->K0*D->K0)*refGam2/aveGam2+(avePx2+avePy2)*refGam2)/(1+D->KRef*D->KRef);
   if(shiftZ>=2.0) {
     printf("shiftZ=%g, aveGam2=%g,refGam2=%g, avePr2=%g\n",shiftZ,aveGam2,refGam2,avePx2+avePy2);
     exit(0);
   } else ;

   // calculating slope memory
   for(h=0; h<numHarmony; h++)  
     for(i=startI; i<endI; i++)
       for(j=0; j<N; j++)
	 D->slope[h][i][j]=D->U[h][i+1][j]-D->U[h][i][j];


   // field update
   for(i=endI-1; i>=startI; i--) {
/*	   
     cnt=0.0;
     aveGam=1.0;
     for(s=0; s<D->nSpecies; s++)  {
       p=D->particle[i].head[s]->pt;
       while(p) {
	 aveGam+=p->gamma;
	 cnt+=1.0;
	 p=p->next;
       }
     }
     if(cnt==0.0) cnt=1.0; else;
     aveGam/=cnt;
*/     

     for(h=0; h<numHarmony; h++)  
       for(j=0; j<N; j++)  {
//         D->U[h][i+1][j]=D->U[h][i][j]+D->slope[h][i][j]*(1.0-shiftZ);		      
         D->U[h][i+1][j]=D->U[h][i][j];
		 }
   }
}

void shiftField_1D(Domain *D,int iteration)
{
   double coef,aveGam,cnt,shiftZ;
   int s,h,numHarmony,i,startI,endI;
   ptclList *p;

   numHarmony=D->numHarmony;
   startI=1;  endI=D->subSliceN+1;

   shiftZ=(1+D->K0*D->K0)/(1+D->KRef*D->KRef);

   // calculating slope memory
   for(h=0; h<numHarmony; h++)  
     for(i=startI; i<endI; i++)
	D->slope[h][i][0]=D->U[h][i+1][0]-D->U[h][i][0];

   // field update
   for(i=endI-1; i>=startI; i--) {
/*	   
     cnt=0.0;
     aveGam=1.0;
     for(s=0; s<D->nSpecies; s++)  {
    :x
    p=D->particle[i].head[s]->pt;
       while(p) {
	 cnt+=1.0;
	 aveGam+=p->gamma;
	 p=p->next;
       }
     }
     if(cnt==0.0) cnt=1.0; else;
     aveGam/=cnt;
     shiftZ=coef/aveGam/aveGam;
*/
     for(h=0; h<numHarmony; h++)  
       D->U[h][i+1][0]=D->U[h][i][0];
//       D->U[h][i+1][0]=D->U[h][i][0]+D->slope[h][i][0]*(1.0-shiftZ);
   }
}

