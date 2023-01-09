#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#include <complex.h>

void main(int argc, char *argv[])
{
   int ii,i,j,numK,numX,numS,numZ,index,ac;
   double maxK,maxX,maxS,minK,minX,minS,dK,dX,dS,k,x,s,reCal;
   double dino,sum,cal,sumFlat,sumCir,ctau;
   double z,minZ,maxZ,dZ,minZCur,maxZCur,givenCur;
   double radius,cond_Al,s0,Z0,totalCharge,coef;
   double complex cdino,csum,ccal;
   double *Zflat,*Zcir,*Wflat,*Wcir,*Eflat,*Ecir,*charge;
   FILE *out;
   char fileName[100];

   if(argc < 6)
   {  printf("impedance ac[1:ON, 0:off] minZ maxZ numZ minZCur maxZCur\n");  exit(0);  }

   ac=atoi(argv[1]);
   minZ=atof(argv[2]);
   maxZ=atof(argv[3]);
   numZ=atoi(argv[4]);
   minZCur=atof(argv[5]);
   maxZCur=atof(argv[6]);

   dZ=(maxZ-minZ)/(numZ*1.0);

//-------------- set up current ----------------//
   givenCur=110; 	//3000;	// Ampere
   charge=(double *)malloc((numZ+1)*sizeof(double ));
   for(i=0; i<=numZ; i++) {
     z=minZ+i*dZ;

     if(z>=minZCur && z<maxZCur) charge[i]=givenCur*dZ/3e8;
     else 	 		 charge[i]=0.0;
//     if(z>=35e-6 && z<50e-6) charge[i]=givenCur*dZ/3e8*(z-35e-6)/15e-6;
//     else if(z>=50e-6 && z<250e-6) charge[i]=givenCur*dZ/3e8;
//     else if(z>=250e-6 && z<265e-6) charge[i]=givenCur*dZ/3e8*((250e-6-z)/15e-6+1);
//     else 	 		 charge[i]=0.0;
   }
   totalCharge=0.0;
   for(i=0; i<=numZ; i++) totalCharge+=charge[i];
   for(i=0; i<=numZ; i++) charge[i]/=totalCharge*dZ;


//-------------- impedance calculation ----------------//
   radius=3.35e-3;
   cond_Al=3.03e7;
   Z0=377;
   ctau=2.4e-6;		// [m]
   s0=pow(2*radius*radius/Z0/cond_Al,1.0/3.0);
//   s0=8.1e-6;
//   ctau=s0;

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
   Wflat=(double *)malloc((numZ+1)*sizeof(double ));
   Wcir=(double *)malloc((numZ+1)*sizeof(double ));
   for(i=0; i<=numK; i++) { Zflat[i]=0.0; Zcir[i]=0.0; }
   for(i=0; i<=numZ; i++) { Wflat[i]=0.0; Wcir[i]=0.0; }

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
	 if(ac==1)  cdino=2.0/(1.0-I)/csqrt(k-I*k*k*ctau/s0)*cosh(x)-I*k*sinh(x)/x;
	 else       cdino=2.0/(1.0-I)/sqrt(k)*cosh(x)-I*k*sinh(x)/x;
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
     if(ac==1) {
        cdino=sqrt(2.0*s0/ctau)/k*(I+s0*0.5/k/ctau)-I*0.5*k;
        ccal=1.0/cdino;
        reCal=creal(ccal);
        Zcir[i]=creal(ccal);
      } else {
        dino=2.0/k+0.25*k*k-sqrt(k);
        Zcir[i]=1.0/sqrt(k)/dino;
      }
   }

   if(ac==1) sprintf(fileName,"impZac");
   else      sprintf(fileName,"impZdc");
   out=fopen(fileName,"w");   
   for(i=0; i<=numK; i++)  {
      k=minK+i*dK;
      fprintf(out,"%g %g %g\n",k,Zflat[i],Zcir[i]);
   }
   fclose(out);
   printf("%s is made\n",fileName);

  
//-------------- wake function calculation ----------------//
   coef=Z0*3e8/M_PI/radius/radius/M_PI;

   for(i=0; i<=numZ; i++)  {
     s=(minZ+i*dZ)/s0;

     sumFlat=0.0;
     sumCir=0.0;
     for(j=0; j<=numK; j++)  {
       k=minK+j*dK;
       cal=Zflat[j]*cos(k*s);
       sumFlat+=cal*dK;
       cal=Zcir[j]*cos(k*s);
       sumCir+=cal*dK;
     }
     Wflat[i]+=sumFlat*coef;  
     Wcir[i]+=sumCir*coef; 

     if(i==0) {
       sumFlat=0.0;
       for(j=0; j<=numK; j++)  {
         k=minK+j*dK;
         cal=Zflat[j]*cos(k*s);
         sumFlat+=cal*dK;
       }
       printf("cal*dK=%g, sum=%g\n",cal*dK,sumFlat);
     }

   } 

   if(ac==1) sprintf(fileName,"wakeFac");
   else      sprintf(fileName,"wakeFdc");
   out=fopen(fileName,"w");   
   for(i=0; i<=numZ; i++)  {
      z=minZ+i*dZ;
      fprintf(out,"%g %g %g\n",z,Wflat[i],Wcir[i]);
   }
   fclose(out);
   printf("%s is made, s0=%g\n",fileName,s0);


//-------------- wake field calculation ----------------//

   coef=totalCharge;  //[eV/m]

   Eflat=(double *)malloc((numZ+1)*sizeof(double ));
   Ecir=(double *)malloc((numZ+1)*sizeof(double ));
   for(i=0; i<=numZ; i++) { Eflat[i]=0.0; Ecir[i]=0.0; }
/*
   for(i=0; i<=numZ; i++)  {

     sumFlat=sumCir=0.0;
     for(ii=0; ii<=i; ii++)  {
       index=i-ii;
       sumFlat+=charge[ii]*Wflat[index]*dZ;
       sumCir+=charge[ii]*Wcir[index]*dZ;
     }
     Eflat[i]=sumFlat;
     Ecir[i]=sumCir;
   } 
*/

   for(i=numZ; i>=0; i--)  {
     for(ii=i; ii>=0; ii--)  {
       index=i-ii;
       Eflat[ii]+=charge[i]*Wflat[index]*dZ*coef;
       Ecir[ii]+=charge[i]*Wcir[index]*dZ*coef;
     }
   } 

   if(ac==1) sprintf(fileName,"wakeEac");
   else      sprintf(fileName,"wakeEdc");
   out=fopen(fileName,"w");   
   for(i=0; i<=numZ; i++)  {
      z=minZ+i*dZ;
      fprintf(out,"%g %g %g %g\n",z,Eflat[i],Ecir[i],charge[i]);
   }
   fclose(out);
   printf("%s is made, totalCharge=%g\n",fileName,totalCharge);

   free(Zflat);
   free(Zcir);
   free(Wflat);
   free(Wcir);
   free(Eflat);
   free(Ecir);
   free(charge);
}
