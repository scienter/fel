#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#include <complex.h>

void main(int argc, char *argv[])
{
   int i,j,numK,numX,numS;
   double maxK,maxX,maxS,minK,minX,minS,dK,dX,dS,k,x,s,reCal;
   double dino,sum,cal,sumFlat,sumCir;
   double complex cdino,csum,ccal;
   double *Zflat,*Zcir,*Wflat,*Wcir;
   FILE *out;

   if(argc < 2)
   {  printf("impZ maxS numS\n");  exit(0);  }

   maxS=atof(argv[1]);
   numS=atoi(argv[2]);

   maxK=20;
   numK=1000;
   maxX=2000;
   numX=50000;

   minK=0.0;
   dK=maxK/(numK*1.0);
   minX=0.0;
   dX=maxX/(numX*1.0);
   minS=0.0;
   dS=maxS/(numS*1.0);

   Zflat=(double *)malloc((numK+1)*sizeof(double ));
   Zcir=(double *)malloc((numK+1)*sizeof(double ));
   Wflat=(double *)malloc((numS+1)*sizeof(double ));
   Wcir=(double *)malloc((numS+1)*sizeof(double ));
   for(i=0; i<=numK; i++) { Zflat[i]=0.0; Zcir[i]=0.0; }
   for(i=0; i<=numS; i++) { Wflat[i]=0.0; Wcir[i]=0.0; }


   for(i=1; i<=numK; i++)  {
     k=minK+i*dK;

     // flat calculation
     csum=0.0+I*0.0;
     for(j=0; j<=numX; j++)  {
       x=minX+j*dX;
       if(x==0.0) {
         ccal=sqrt(k)/(2.0+k*k*k-2*sqrt(k)*k);
       } else if(x>320) { 	       
         ccal=0.0+I*0.0;
         reCal=creal(ccal);
         if(isnan(reCal)) { printf("cal=%g, x=%g, k=%g, dino=%g\n",reCal,x,k,creal(cdino)); exit(0); }
       } else { 	       
         cdino=2.0/(1.0-I)/sqrt(k)*cosh(x)-I*k*sinh(x)/x;
         ccal=1.0/cosh(x)/cdino;
         reCal=creal(ccal);
         if(isnan(reCal)) { printf("cal=%g, x=%g, k=%g, dino=%g\n",reCal,x,k,creal(cdino)); exit(0); }
       }
       csum+=ccal*dX;
     }
     Zflat[i]=creal(csum);

     // circular calculation
//     cdino=(1.0+I)/sqrt(k)-I*k*0.5;
//          Zcir[i]=creal(1.0/cdino);
     dino=2.0/k+0.25*k*k-sqrt(k);
     Zcir[i]=1.0/sqrt(k)/dino;
   }

  
   // wakefield
   for(i=0; i<=numS; i++)  {
     s=minS+i*dS;

     sumFlat=0.0;
     sumCir=0.0;
     for(j=0; j<=numK; j++)  {
       k=minK+j*dK;
       cal=Zflat[j]*cos(k*s);
       sumFlat+=cal*dK;
       cal=Zcir[j]*cos(k*s);
       sumCir+=cal*dK;
     }
     Wflat[i]=sumFlat;  
     Wcir[i]=sumCir; 
   } 

   out=fopen("impZ","w");   
   for(i=0; i<=numK; i++)  {
      k=minK+i*dK;
      fprintf(out,"%g %g %g\n",k,Zflat[i],Zcir[i]);
   }
   fclose(out);
   printf("impZ is made\n");

   out=fopen("thwakeZ","w");   
   for(i=0; i<=numS; i++)  {
      s=minS+i*dS;
      fprintf(out,"%g %g %g\n",s,Wflat[i],Wcir[i]);
   }
   fclose(out);
   printf("thwakeZ is made\n");

   free(Zflat);
   free(Zcir);
   free(Wflat);
   free(Wcir);
}
