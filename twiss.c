#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"

void calculate_twiss(Domain *D,int iteration) 
{
   int i,s,cenI;
   double cnt,invGam,x,y,xPrime,yPrime;
   double aveX2,aveXPrime,aveCrsX,emittanceX;
   double aveY2,aveYPrime,aveCrsY,emittanceY;
   ptclList *p;

   cenI=D->sliceN/2;
   aveX2=0.0; aveXPrime=0.0; aveCrsX=0.0;
   aveY2=0.0; aveYPrime=0.0; aveCrsY=0.0;

   if(cenI>=D->minI && cenI<D->maxI) {
     cnt=0.0;
     for(s=0; s<D->nSpecies; s++)  {
       p=D->particle[cenI+1-D->minI].head[s]->pt;
       while(p) {
         x=p->x; 
         y=p->y; 
         invGam=1.0/p->gamma;

         xPrime=p->px*invGam;
         yPrime=p->py*invGam;

         aveX2+=x*x;
         aveY2+=y*y;
         aveXPrime+=xPrime*xPrime;
         aveYPrime+=yPrime*yPrime;
         aveCrsX+=x*xPrime;
         aveCrsY+=y*yPrime;

         cnt+=1.0;
         p=p->next;
       }
     }
     emittanceX=sqrt((aveX2*aveXPrime-aveCrsX*aveCrsX)/cnt/cnt);
     D->twsBX[iteration]=(aveX2/cnt)/emittanceX;
     D->twsGX[iteration]=(aveXPrime/cnt)/emittanceX;
     D->twsAX[iteration]=-aveCrsX/cnt/emittanceX;
     D->twsEmitX[iteration]=emittanceX;
//   printf("iteration=%d, sigX=%g, stdX=%g, sigXPrime=%g, stdPx=%g\n",iteration,D->loadList->sigX,sqrt(aveX2/cnt),emittanceX/D->loadList->sigX,sqrt(aveXPrime/cnt));

     emittanceY=sqrt((aveY2*aveYPrime-aveCrsY*aveCrsY)/cnt/cnt);
     D->twsBY[iteration]=(aveY2/cnt)/emittanceY;
     D->twsGY[iteration]=(aveYPrime/cnt)/emittanceY;
     D->twsAY[iteration]=-aveCrsY/cnt/emittanceY;
     D->twsEmitY[iteration]=emittanceY;

     D->twsG[iteration]=D->g;
   } else ;
}

void save_twiss(Domain *D)
{
   int i;
   double dz;
   char name[100];
   FILE *out;

   dz=D->dz;
   sprintf(name,"twiss");       
   out = fopen(name,"w");  
   for(i=0; i<D->maxStep; i++)  { 
     fprintf(out,"%g %.10g %g %.10g %g ",i*dz,D->twsEmitX[i],D->twsBX[i],D->twsGX[i],D->twsAX[i]);
     fprintf(out,"%.10g %g %.10g %g %g\n",D->twsEmitY[i],D->twsBY[i],D->twsGY[i],D->twsAY[i],D->twsG[i]);
   }
   fclose(out);
   printf("%s is made.\n",name);                                                                                
  
}

