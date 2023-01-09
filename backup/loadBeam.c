#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <gsl/gsl_qrng.h>

double randomValue(double beta);
void loadBeam1D(Domain *D,LoadList *LL,int s,int iteration);
void loadBeam3D(Domain *D,LoadList *LL,int s,int iteration);
double gaussianDist_1D(double sigma);
void random_2D(double *x,double *y,gsl_qrng *q1);


void loadBeam(Domain D,LoadList *LL,int s,int iteration)
{
  switch((LL->type-1)*3+D.dimension)  {
  case ((Polygon-1)*3+1):
    loadBeam1D(&D,LL,s,iteration);
    break;
  case ((Polygon-1)*3+2):
//    loadPolygonPlasma2D(D,LL,s,iteration); 
    break;
  case ((Polygon-1)*3+3):
    loadBeam3D(&D,LL,s,iteration);
    break;
  default:
    ;
  }
}


void loadBeam3D(Domain *D,LoadList *LL,int s,int iteration)
{
   int l,b,n,numInBeamlet,beamlets,noiseONOFF,flag;
   int i,startI,endI,minI,maxI,finalFlag;
   double posZ,current,n0,dEnergy,nEnergy,bucketZ,dPhi,div,ptclCnt;
   double macro,remacro,theta,sigGam,dGam,gam,gamma0,Ns,noise;
   double sigX,sigY,emitX,emitY,gammaX,gammaY,x,y,pz,px,py,vz,eTune;
   double sigXPrime,sigYPrime,xPrime,yPrime,delTX,delTY,distanceX,distanceY;
   double y1,y2,testNum,numbers,coef,sqrt2,tmp;
   int myrank;
   ptclList *New;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   minI=D->minI;   maxI=D->maxI;
   startI=1;       endI=1+D->subSliceN;
   
   numInBeamlet=LL->numInBeamlet;
   gamma0=LL->energy/mc2;
   dGam=LL->spread*gamma0;
   emitX=LL->emitX/gamma0/M_PI;
   emitY=LL->emitY/gamma0/M_PI;
   gammaX=(1+LL->alphaX*LL->alphaX)/LL->betaX;
   gammaY=(1+LL->alphaY*LL->alphaY)/LL->betaY;   
   sigX=sqrt(emitX/gammaX);
   sigY=sqrt(emitY/gammaY);
   sigXPrime=sqrt(emitX*gammaX);
   sigYPrime=sqrt(emitY*gammaY);

   distanceX=sqrt((LL->betaX-1.0/gammaX)/gammaX);
   distanceY=sqrt((LL->betaY-1.0/gammaY)/gammaY);
   vz=sqrt(gamma0*gamma0-1.0)/gamma0;	//normalized
   if(vz==0.0) { delTX=delTY=0.0; }
   else  {
     delTX=distanceX/vz;	//normalized
     delTY=distanceY/vz;	//normalized
   }

   printf("delTX=%g, delTY=%g, betaX=%g, betaY=%g\n",delTX,delTY,LL->betaX,LL->betaY);

   current=LL->peakCurrent;		// peak current in a cell
   bucketZ=D->lambda0*D->numSlice;	// size of a big slice
   ptclCnt=numInBeamlet*LL->numBeamlet;	

   dPhi=2.0*M_PI*D->numSlice;
   div=2.0*M_PI/(1.0*numInBeamlet);
   LL->index=0;
   noiseONOFF=LL->noiseONOFF;

   sigGam=gamma0*LL->spread;
   macro=current/eCharge/velocityC*bucketZ/ptclCnt;
   sqrt2=sqrt(2.0);

   LL->totalCnt=LL->numBeamlet*LL->numInBeamlet;
   eTune=LL->eTune;

   //position define   
//   srand(myrank);
   srand((unsigned int)time(NULL));

   double v[4],**randNum;
//   gsl_qrng *q1=gsl_qrng_alloc(gsl_qrng_halton,2);
   gsl_qrng *q1=gsl_qrng_alloc(gsl_qrng_niederreiter_2,4);
   gsl_qrng *q2=gsl_qrng_alloc(gsl_qrng_niederreiter_2,4);
//   gsl_qrng *q2=gsl_qrng_alloc(gsl_qrng_reversehalton,2);
//   gsl_qrng *q1=gsl_qrng_alloc(gsl_qrng_reversehalton,2);
//   gsl_qrng *q1=gsl_qrng_alloc(gsl_qrng_sobol,2);
//   gsl_qrng *q2=gsl_qrng_alloc(gsl_qrng_sobol,2);
   randNum=(double **)malloc(LL->numBeamlet*2*sizeof(double *));
   for(i=0; i<LL->numBeamlet*2; i++)
     randNum[i]=(double *)malloc(4*sizeof(double ));
   for(i=0; i<LL->numBeamlet; i+=1) {
     flag=0;
     while(flag==0) {
       gsl_qrng_get(q1,v);
       if(v[0]==0.0 || v[0]==1.0) flag=0; else flag=1 ;
     }
     randNum[i][0]=v[0];
     randNum[i][1]=v[1];
     randNum[i][2]=v[2];
     randNum[i][3]=v[3];
   }

   for(i=LL->numBeamlet; i<LL->numBeamlet*2; i+=1) {
     flag=0;
     while(flag==0) {
       gsl_qrng_get(q1,v);
       if(v[0]==0.0 || v[0]==1.0) flag=0; else flag=1 ;
     }
     randNum[i][0]=v[0];
     randNum[i][1]=v[1];
     randNum[i][2]=v[2];
     randNum[i][3]=v[3];
   }
   gsl_qrng_free(q1);
   gsl_qrng_free(q2);
     


//   gsl_qrng *q2=gsl_qrng_alloc(gsl_qrng_sobol,2);
   for(i=startI; i<endI; i++) {
     posZ=(i-startI+minI+0.5)*bucketZ+D->minZ;
     for(l=0; l<LL->znodes-1; l++) {
       if(posZ>=LL->zpoint[l] && posZ<LL->zpoint[l+1])
       {
         n0=((LL->zn[l+1]-LL->zn[l])/(LL->zpoint[l+1]-LL->zpoint[l])*(posZ-LL->zpoint[l])+LL->zn[l]);
         dEnergy=((LL->zenergy[l+1]-LL->zenergy[l])/(LL->zpoint[l+1]-LL->zpoint[l])*(posZ-LL->zpoint[l])+LL->zenergy[l]);
	 nEnergy=1.0+dEnergy*eTune;
         noise=sqrt(3.0/macro)*noiseONOFF;
         beamlets=(int)(LL->numBeamlet*n0);
	 remacro=macro*(1.0*LL->numBeamlet)/(1.0*beamlets)*n0;

//         gsl_qrng *q2=gsl_qrng_alloc(gsl_qrng_niederreiter_2,2);
         for(b=0; b<beamlets; b++)  {
//         theta=dPhi*randomValue(1.0);
           theta=dPhi*(b+0.5)/(1.0*beamlets);
           dGam=gaussianDist_1D(sigGam);
           gam=gamma0*nEnergy+dGam;

	   //           x=gaussianDist_1D(sigX*sqrt2);
//           xPrime=gaussianDist_1D(sigXPrime*sqrt2);
//           y=gaussianDist_1D(sigY*sqrt2);
//           yPrime=gaussianDist_1D(sigYPrime*sqrt2);

           coef=sqrt(-2.0*log(randNum[b][0]));
           x=coef*cos(2*M_PI*randNum[b][1]);
           xPrime=coef*sin(2*M_PI*randNum[b][1]);
	   x*=sigX;
	   xPrime*=sigXPrime;

           coef=sqrt(-2.0*log(randNum[b][2]));
           y=coef*cos(2*M_PI*randNum[b][3]);
           yPrime=coef*sin(2*M_PI*randNum[b][3]);
	   y*=sigY;
	   yPrime*=sigYPrime;
//           coef=sqrt(-2.0*log(randNum[beamlets+b][0]));
//           y=coef*cos(2*M_PI*randNum[beamlets+b][1]);
//           yPrime=coef*sin(2*M_PI*randNum[beamlets+b][1]);
//	   y*=sigY;
//	   yPrime*=sigYPrime;

//	   if(isnan(y) || isnan(x) || isnan(xPrime) || isnan(yPrime)) {
//printf("y=%g, x=%g, xPrime=%g, yPrime=%g\n",y,x,xPrime,yPrime);
//printf("here y,v[0]=%g, v[1]=%g\n",v[0],v[1]);
//			   }

           pz=sqrt((gam*gam-1.0)/(1.0+xPrime*xPrime+yPrime*yPrime));
           px=xPrime*pz;
           py=yPrime*pz;

	   x-=delTX*px/gam;
	   y-=delTY*py/gam;
           for(n=0; n<numInBeamlet; n++)  {		     
             New = (ptclList *)malloc(sizeof(ptclList));
             New->next = D->particle[i].head[s]->pt;
             D->particle[i].head[s]->pt = New;

             New->x = x;
             New->y = y;
//             New->theta=theta+n*div+randomValue(1.0)*noise;
             tmp=theta+n*div+randomValue(1.0)*noise;
             if(tmp>=dPhi) tmp-=dPhi; else ;
             New->theta=tmp;
	     
             New->gamma=gam;
             New->px=px;
             New->py=py;   	
             New->weight=remacro;
             New->index=LL->index;  	//index
             New->core=myrank;  	
             LL->index+=1;
           }
         }
//         gsl_qrng_free(q2);
       } else ;		//End of for(l,posZ)
     }
   }     

   for(i=0; i<LL->numBeamlet*2; i++) free(randNum[i]); free(randNum);
}


void loadBeam1D(Domain *D,LoadList *LL,int s,int iteration)
{

   int i,l,b,n,numInBeamlet,beamlets,noiseONOFF,minI,maxI,cnt;
   int startI,endI;
   double posZ,current,n0,nEnergy,dEnergy,eTune,bucketZ,dPhi,div,ptclCnt,tmp;
   double macro,theta,sigGam,dGam,gam,gamma0,Ns,noise;
   double sigX,sigY,emitX,emitY,x,y,pz,px,py;
   double sigXPrime,sigYPrime,xPrime,yPrime,sigPx,sigPy;
   double y1,y2,testNum,numbers,coef,sqrt2;
   int myrank,nTasks;
   ptclList *New;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   minI=D->minI;   maxI=D->maxI;
   startI=1;	   endI=1+D->subSliceN;

   numInBeamlet=LL->numInBeamlet;
   gamma0=LL->energy/mc2+1;
   dGam=LL->spread*gamma0;

   current=LL->peakCurrent;		// peak current in a cell
   bucketZ=D->lambda0*D->numSlice;		     	// size of a big slice
   ptclCnt=numInBeamlet*LL->numBeamlet;	

   dPhi=2.0*M_PI*D->numSlice;
   div=2.0*M_PI/(1.0*numInBeamlet);
   sqrt2=sqrt(2.0);
   LL->index=0;
   noiseONOFF=LL->noiseONOFF;

   sigGam=gamma0*LL->spread;
   macro=current/eCharge/velocityC*bucketZ/ptclCnt;

   LL->totalCnt=LL->numBeamlet*LL->numInBeamlet;
   eTune=LL->eTune;

   cnt=0;
   for(i=startI; i<endI; i++) {
     //position define     
     srand(i-startI+minI);
     posZ=(i-startI+minI)*bucketZ+D->minZ;
     for(l=0; l<LL->znodes-1; l++) {
       if(posZ>=LL->zpoint[l] && posZ<LL->zpoint[l+1])
       {
         n0=((LL->zn[l+1]-LL->zn[l])/(LL->zpoint[l+1]-LL->zpoint[l])*(posZ-LL->zpoint[l])+LL->zn[l]);
         dEnergy=((LL->zenergy[l+1]-LL->zenergy[l])/(LL->zpoint[l+1]-LL->zpoint[l])*(posZ-LL->zpoint[l])+LL->zenergy[l]);
	 nEnergy=1.0+dEnergy*eTune;
         noise=sqrt(3.0/macro)*noiseONOFF;
         beamlets=LL->numBeamlet*n0;

         for(b=0; b<beamlets; b++)  {
//         theta=dPhi*randomValue(1.0);
           theta=dPhi*(b+0.5)/(1.0*beamlets);
           dGam=gaussianDist_1D(sigGam);
           gam=gamma0*nEnergy+dGam;
           for(n=0; n<numInBeamlet; n++)  {		     
             New = (ptclList *)malloc(sizeof(ptclList));
             New->next = D->particle[i].head[s]->pt;
             D->particle[i].head[s]->pt = New;

             New->x = 0.0;
             New->y = 0.0;
             tmp=theta+n*div+randomValue(1.0)*noise;
	     if(tmp>=dPhi) tmp-=dPhi; else ;
             New->theta=tmp;
             New->gamma=gam;		//gamma
             New->px=0.0;
             New->py=0.0;       	//py
             New->weight=macro;
             New->index=LL->index;  	//index
             New->core=myrank;  	//index
             LL->index+=1;
	     cnt+=1;
           }
         }
	 l=LL->znodes;
       } else ;		//End of for(l,posZ)
     }         
   }			//End of for(i)
   D->ptclCnt=cnt;
   printf("myrank=%d, ptclCnt=%d\n",myrank,D->ptclCnt);
}

void random_2D(double *x,double *y,gsl_qrng *q1)
{
   double v[2];

      gsl_qrng_get(q1,v);
      *x=v[0];
      *y=v[1];
}

double gaussianDist_1D(double sigma)
{
   double r,prob,v,z,random;
   int intRand,randRange=1e4;

   r=1.0;
   prob=0.0;
//   gsl_qrng *q=gsl_qrng_alloc(gsl_qrng_niederreiter_2,1);
   while (r>prob)  {
      intRand = rand() % randRange;
      r = ((double)intRand)/randRange;
      intRand = rand() % randRange;
      random = ((double)intRand)/randRange;
//      gsl_qrng_get(q,&v);
      z = 4.0*(random-0.5);	//up to 3 sigma
      prob=exp(-z*z);
   }
//   gsl_qrng_free(q); 
  
   return z*sigma;
}

double randomValue(double beta)
{
   double r;
   int intRand, randRange=1000, rangeDev;

   rangeDev=(int)(randRange*(1.0-beta));
   intRand = rand() % (randRange-rangeDev);
   r = ((double)intRand)/randRange+(1.0-beta);

   return r;
}

