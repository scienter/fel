#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"

// For now, plane undulator is applied.
void updateK_quadG(Domain *D,int iteration,double half)
{
   int n,exist,airPosition,noCurrent;
   double z,K0,g,dz,q0,q1,x0,x1;
   int myrank, nTasks;
   UndulatorList *UL;
   QuadList *QD;

   dz=D->dz;
//   z=(iteration+half)*D->dz+sliceI*D->numSlice*D->lambda0+D->minZ;
   z=(iteration+half)*D->dz;

   //-------------- update K -----------------//
   UL=D->undList;

   K0=D->prevK;
   exist=0;
   airPosition=0;
   noCurrent=OFF;
   while(UL->next) {
     for(n=0; n<UL->numbers; n++) {
       if(z>=UL->undStart[n] && z<UL->undEnd[n]) {
         K0=UL->K0[n]*(1+UL->taper*(z-UL->undStart[n]));
	      if(UL->noCurrent==ON) noCurrent=ON; else ;
	      exist=1;
       } else if(z>=UL->unitStart[n] && z<UL->unitEnd[n] && UL->air==ON) {
	      airPosition=1;
       } else ;
     }
     UL=UL->next; 
   }
   
   if(airPosition==1) {
     K0=0.0;
   } else ;
   D->K0=K0;

   // save prev B field.
   D->prevK=K0;
   if(exist==1 && noCurrent==OFF)     D->currentFlag=1;
   else if(exist==1 && noCurrent==ON) D->currentFlag=0;
   else	                              D->currentFlag=0;

   //-------------- update kx, ky -----------------//
   if(UL->undMode==Linear) { D->kx=0.0; D->ky=D->ku; }
	else if(UL->undMode==Helical) { D->kx=D->ku; D->ky=D->ku; }
	else  { D->kx=0.0; D->ky=0.0; }



   //-------------- update Quad -----------------//
   g=0;
   QD=D->qdList;
   exist=0;
   x0=z;
   x1=z+dz*0.5;
   while(QD->next && exist==0) {
     for(n=0; n<QD->numbers; n++) {
       q0=QD->qdStart[n];
       q1=QD->qdEnd[n];
       if(q1<=x0 || q0>=x1) ;
       else {
         if(x0<q0 && q1<x1)    
         { g=(q1-q0)/dz*2*QD->g[n]; exist=1; n=QD->numbers; }
         else if(x0<q0 && q0<x1)  
         { g=(x1-q0)/dz*2*QD->g[n]; exist=1; n=QD->numbers; }
         else if(x0<q1 && q1<x1) 
         { g=(q1-x0)/dz*2*QD->g[n]; exist=1; n=QD->numbers; }
         else                  
         { g=1.0*QD->g[n];        exist=1; n=QD->numbers; }
       }
     }
     QD=QD->next; 
   }
   D->g=g;
}

