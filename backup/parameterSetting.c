#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <gsl/gsl_sf_bessel.h>

int findBeamLoadParameters(int rank, LoadList *LL,Domain *D,char *input);
int whatBeamType(char *str);
int whatSpecies(char *str);

double randomV()
{
   double r;
   int intRand, randRange=1000, rangeDev;

   intRand = rand() % randRange;
   r = ((double)intRand)/randRange;

   return r;
}

void parameterSetting(Domain *D,char *input)
{
   LoadList *LL,*New;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   FILE *in=NULL;
   int FindParameters();
   int whatONOFF();
   int whatSaveMode();
   char str[100],name[100],fileName[100];
   int rank,fail=0;
   double B0,tmp;

   //initially
   if(FindParameters("Domain",1,"dimension",input,str)) D->dimension=atoi(str);
   else  {
      printf("in [Domain], dimension=?  (1:1D, 2:2D, 3:3D)\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"harmony",input,str)) D->harmony=atoi(str);
   else  {
      printf("in [Domain], harmony=?  (max harmony)\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"zmode",input,str)) D->zmode=atoi(str);
   else  {
      printf("in [Domain], zmode=?  (number of Ez modes)\n");
      fail=1;
   }


   //Save parameter setting
   if(FindParameters("Save",1,"save_step",input,str)) D->saveStep=atoi(str);
   else  {
      printf("In [Save], save_step=?\n");
      fail=1;
   }
   if(FindParameters("Save",1,"save_start",input,str)) D->saveStart=atoi(str);
   else  {
      printf("In [Save], save_start=?\n");
      fail=1;
   }
   if(FindParameters("Save",1,"dump_save",input,str)) 
         D->dumpSave=whatONOFF(str);
   else  D->dumpSave=OFF;
   if(FindParameters("Save",1,"dump_start",input,str)) 
         D->dumpStart=atoi(str);
   else  D->dumpStart=D->saveStart;
   if(FindParameters("Save",1,"dump_step",input,str)) 
         D->dumpStep=atoi(str);
   else  D->dumpStep=D->saveStep;
   if(FindParameters("Save",1,"field_save",input,str)) 
         D->fieldSave=whatONOFF(str);
   else  D->fieldSave=ON;
   if(FindParameters("Save",1,"particle_save",input,str)) 
         D->particleSave=whatONOFF(str);
   else  D->particleSave=ON;
   if(FindParameters("Save",1,"density_save",input,str)) 
         D->rhoSave=whatONOFF(str);
   else  D->rhoSave=ON;
   if(FindParameters("Save",1,"field_format",input,str)) 
         D->saveFieldMode=whatSaveMode(str);
   else  D->saveFieldMode=TXT;
   if(FindParameters("Save",1,"particle_format",input,str)) 
         D->saveParticleMode=whatSaveMode(str);
   else  D->saveParticleMode=TXT;
   if(FindParameters("Save",1,"density_format",input,str)) 
         D->saveDensityMode=whatSaveMode(str);
   else  D->saveDensityMode=TXT;
   if(FindParameters("Save",1,"dump_format",input,str)) 
         D->saveDumpMode=whatSaveMode(str);
   else  D->saveDumpMode=TXT;

   //Domain parameter setting
   if(FindParameters("Domain",1,"max_time",input,str)) D->maxTime=atoi(str);
   else  D->maxTime=525600;
   if(FindParameters("Domain",1,"max_step",input,str)) D->maxStep=atoi(str);
   else  {
      printf("In [Domain], max_step=? [ea].\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"lambdaUs_in_iteration",input,str)) D->numLambdaU=atoi(str);
   else  D->numLambdaU=1;

   //Electron beam
   if(FindParameters("Domain",1,"ref_energy",input,str)) D->refEnergy=atof(str);
   else  {
      printf("In [Domain], ref_energy=? [MeV].\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"ref_duration",input,str)) D->refDuration=atof(str);
   else  {
      printf("In [Domain], ref_duration=? [s].\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"ref_radius",input,str)) D->refRadius=atof(str);
   else  {
      printf("In [Domain], ref_radius=? [m].\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"slices_in_bucket",input,str)) D->numSlice=atoi(str);
   else  {
      printf("In [Domain], slices_in_bucket=? [ea].\n");
      fail=1;
   }
   if(D->numSlice<D->numLambdaU) {   
      printf("In [Domain], check the condition 'slices_in_bucket(=%d) >= lambdaUs_in_itertation(=%d).\n",D->numSlice,D->numLambdaU);
      fail=1;
   }

   //Undulator
   if(FindParameters("Domain",1,"B0_x",input,str)) D->B0x=atof(str);
   else  {
      printf("In [Domain], B0_x=? [T].\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"B0_y",input,str)) D->B0y=atof(str);
   else  {
      printf("In [Domain], B0_y=? [T].\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"lambdaU",input,str)) D->lambdaU=atof(str);
   else  {
      printf("In [Domain], lambdaU=? [m].\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"undulator_length",input,str)) D->uLength=atof(str);
   else  {
      printf("In [Domain], undulator_length=? [m].\n");
      fail=1;
   }


   //Additional parameters.
   D->ku0=2.0*M_PI/D->lambdaU;
   B0=sqrt(D->B0x*D->B0x+D->B0y*D->B0y);
   D->K0=1.0/sqrt(2.0)*eCharge*B0/eMass/velocityC/D->ku0;
   D->gamma0=D->refEnergy/mc2;
   D->lambda0=D->lambdaU*0.5/D->gamma0/D->gamma0*(1.0+D->K0*D->K0);
   D->k0=2.0*M_PI/D->lambda0;
		  
   tmp=D->K0*D->K0*0.5/(1.0+D->K0*D->K0);
   D->modK0=gsl_sf_bessel_J0(tmp)-gsl_sf_bessel_J1(tmp);
   D->area=M_PI*D->refRadius*D->refRadius;
   D->beta0=1.0-0.5*(1.0+D->K0*D->K0)/D->gamma0/D->gamma0;

   D->sliceN=D->refDuration*velocityC/D->lambda0;
   D->sliceN/=D->numSlice;
   D->sliceN=1;

   D->iterN=D->uLength/D->lambdaU;
   D->iterN/=D->numSlice;
   D->dz=D->lambdaU*D->numLambdaU;   

   //future update
   D->kx=0;
   D->ky=0;

   //Printing basic things.
   if(myrank==0) {
     printf("iterN=%d, numSlice=%d, sliceN=%d\n",D->iterN,D->numSlice,D->sliceN);
     printf("lambda0=%g[nm], lambdaU=%g[cm], dz=%g[cm]\n",D->lambda0*1e9,D->lambdaU*100,D->dz*100);
     printf("B0=%g[T], K0=%g, modK0=%g\n",B0,D->K0,D->modK0);
   }
   else ;
   MPI_Barrier(MPI_COMM_WORLD);

   //Beam parameter setting
   D->loadList = (LoadList *)malloc(sizeof(LoadList));
   D->loadList->next = NULL;
   LL = D->loadList;
   rank = 1;
   while(findBeamLoadParameters(rank, LL, D,input)) 
   {
      New = (LoadList *)malloc(sizeof(LoadList));
      New->next = NULL;
      LL->next=New;
      LL=LL->next;
      rank ++;
   }
   D->nSpecies = rank-1;

   if(fail==1)
      exit(0);
   else	;

}


int findBeamLoadParameters(int rank, LoadList *LL,Domain *D,char *input)
{
   int FindParameters();
   LoadList *New;
   char name[100], str[100];
   int i,n,cnt,species,fail=0;
   double tmp,max,min;
   double *shareDouble;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   if(FindParameters("EBeam",rank,"type",input,name)) LL->type = whatBeamType(name);
   else  LL->type=0;

   if(LL->type>0)
   {
     if(FindParameters("EBeam",rank,"species",input,name)) 
        species = whatSpecies(name);
     else  species = 0;
     LL->species=species;
     if(FindParameters("EBeam",1,"beam_energy",input,str)) LL->bEnergy=atof(str);
     else  LL->bEnergy=D->refEnergy;
     if(FindParameters("EBeam",1,"energy_spread",input,str)) LL->spread=atof(str);
     else  {
        printf("In [EBeam], energy_spread=? .\n");
        fail=1;
     }
     if(FindParameters("EBeam",1,"beam_current",input,str)) LL->bCurrent=atof(str);
     else  {
        printf("In [EBeam], beam_current=? [A].\n");
        fail=1;
     }
     if(FindParameters("EBeam",1,"beam_duration",input,str)) LL->bDuration=atof(str);
     else  LL->bDuration=D->refDuration;
     if(FindParameters("EBeam",1,"beam_radius",input,str)) LL->bRadius=atof(str);
     else  LL->bRadius=D->refRadius;
     if(FindParameters("EBeam",1,"beamlets_in_bucket",input,str)) LL->numBeamlet=atoi(str);
     else  {
        printf("In [EBeam], beamlets_in_bucket=? [ea/bucket].\n");
        fail=1;
     }
     if(FindParameters("EBeam",1,"number_in_beamlet",input,str)) LL->numInBeamlet=atoi(str);
     else  {
        printf("In [EBeam], number_in_beamlet=? [ea/beamlet].\n");
        fail=1;
     }

     switch (LL->type)  {
     case Polygon :
       if(FindParameters("EBeam",rank,"Xnodes",input,str)) LL->xnodes=atoi(str);
       else  {
         printf("in [EBeam], Xnodes=?\n");
         printf("Each nodes indicates the point of beam density changing.\n");
         fail=1;
       }
       if(LL->xnodes>0)  {
          LL->xpoint = (double *)malloc(LL->xnodes*sizeof(double));
          LL->xn = (double *)malloc(LL->xnodes*sizeof(double));   
          for(i=0; i<LL->xnodes; i++)
          {
            sprintf(name,"X%d",i);
            if(FindParameters("EBeam",rank,name,input,str)) 
              LL->xpoint[i] = atof(str);
            else 
            { printf("X%d should be defined.\n",i);  fail=1; }

            sprintf(name,"Xn%d",i);
            if(FindParameters("EBeam",rank,name,input,str)) 
              LL->xn[i] = atof(str);
            else 
            { printf("Xn%d should be defined.\n",i);  fail=1; } 
          }
       }
       else ;
     }		//End of swith
   
   }	//end of if(LL->type>0)
   else ;

   if(fail==1)
      exit(0);

   return LL->type;
}
     
     
     
/*
int findLaserParameters(int rank, LaserList *L,Domain *D,char *input)
{
   int FindParameters();
   double positionX,positionY,positionZ;
   char name[100], str[100];
   int fail=0,polarity;

   if(FindParameters("Laser",rank,"polarity",input,str)) polarity=atoi(str);
   else  polarity=0;

   if(polarity)
   {
     if(FindParameters("Laser",rank,"mode",input,str)) 
        L->mode=atoi(str);
     else  L->mode=0;
     if(FindParameters("Laser",rank,"wavelength",input,str)) 
     {
        L->lambda=atof(str);
        L->lambda*=D->gamma*(1.0+D->beta);
     }
     else  L->lambda=D->lambda;
  
     if(FindParameters("Laser",rank,"a0",input,str)) 
        L->amplitude=atof(str);
     else  {
        printf("in [Laser], a0=??\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"rU",input,str)) L->rU=atof(str);
     else  {
        printf("in [Laser], rU=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"rD",input,str)) L->rD=atof(str);
     else  {
        printf("in [Laser], rD=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"retard",input,str)) L->retard=atof(str);
     else  {
        printf("in [Laser], retard=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"loadPositionX",input,str)) positionX=atof(str);
     else  {
        printf("in [Laser], loadPositionX=?  [m]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"loadPositionY",input,str)) positionY=atof(str);
     else  positionY=0;
     if(FindParameters("Laser",rank,"loadPositionZ",input,str)) positionZ=atof(str);
     else  positionZ=0;
     if(FindParameters("Laser",rank,"beamWaist",input,str)) L->beamWaist=atof(str);
     else  {
        printf("in [Laser], beamWaist=?  [m]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"focus",input,str)) L->focus=atof(str);
     else  {
        printf("in [Laser], focus=?  [m]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"flat",input,str)) L->flat=atof(str);
     else  L->flat=0.0;
     if(FindParameters("Laser",rank,"direction",input,str)) L->direction=atoi(str);
     else  L->direction=1;

     //additional laser parameters
     L->polarity=polarity;
     L->omega=2*pi*velocityC/L->lambda;
     L->loadPointX=((int)(positionX/D->lambda/D->dx));   
     L->loadPointY=((int)(positionY/D->lambda/D->dy));   
     L->loadPointZ=((int)(positionZ/D->lambda/D->dz));   
     L->rayleighLength=pi/(L->lambda/D->gamma/(1.0+D->beta))*L->beamWaist*L->beamWaist/D->lambda;
     L->beamWaist=L->beamWaist/D->lambda;
     L->focus=L->focus/D->lambda;
     if(fail==1)
        exit(0);
   }
   return polarity;
}
*/


int whatONOFF(char *str)
{
   if(strstr(str,"ON")) 		return ON;
   else if(strstr(str,"OFF"))   	return OFF;
   else 				return OFF;
}

int whatSaveMode(char *str)
{
   if(strstr(str,"TXT")) 		return TXT;
   else if(strstr(str,"HDF"))   	return HDF;
   else 				return TXT;
}


int whatBeamType(char *str)
{
   if(strstr(str,"Polygon"))        	return Polygon; 
   else if(strstr(str,"Gaussian"))     	return Gaussian;
   else   {
     printf("No Beam type!\n"); 
     return 0;
     exit(0);
   }
}

int whatSpecies(char *str)
{
   if(strstr(str,"Electron"))         	return Electron;
   else 				return 0;
}



