#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <complex.h>
#include "mpi.h"
#include <gsl/gsl_sf_bessel.h>


void initialFileSave(Domain *D)
{
   int h;
   FILE *out1,*out2,*out3,*out4;

   int myrank, nTasks;    
   MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   if(myrank==0) {
	   out1=fopen("totalEnergy","w");
		out2=fopen("bunching","w");
		out3=fopen("twissFile","w");
      out4=fopen("divergence","w");
		   fprintf(out4,"z[m] ");
			for(h=0; h<D->numHarmony; h++) 
            fprintf(out4,"tot_divX%d tot_divY%d ",h,h);
	      for(h=0; h<D->numHarmony; h++) 
            fprintf(out4,"cen_divX%d cen_divY%d unit[rad/um]",h,h);
         fprintf(out4,"\n");
      fclose(out1);
      fclose(out2);
      fclose(out3);
      fclose(out4);
   } else ;

	MPI_Barrier(MPI_COMM_WORLD);
}



void updateFELCharacter(Domain *D,int iteration)
{
   int i,j,h,startI,endI,minI,rank,numHarmony,cenI,sliceI,start,nx,ny,dataNum;
	double z,dx,dy,dz,minZ,angle,cnt,divX,divY,val,phase1,phase2;
	double data[D->numHarmony*2*2],sendData[D->numHarmony*2*2],recvData[D->numHarmony*2*2];
   char name[100];
   FILE *out;

   int myrank, nTasks;    
   MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   startI=1;  endI=D->subSliceN+1;
	numHarmony=D->numHarmony;
	minI=D->minI;
   dx=D->dx;  dy=D->dy;  dz=D->dz;
   minZ=D->minZ;
	cenI=D->sliceN/2;
	nx=D->nx; ny=D->ny;
   
	start=0;
   for(h=0; h<numHarmony; h++) {
	   data[start+h]=0.0;
	   data[start+h+1]=0.0;
		start+=2;
	}
   for(h=0; h<numHarmony; h++) {
	   data[start+h]=0.0;
	   data[start+h+1]=0.0;
		start+=2;
	}

   //total divergence
   start=0;
	for(h=0; h<numHarmony; h++)  {
      for(sliceI=startI; sliceI<endI; sliceI++)  {
         for(j=1; j<ny-1; j++)  
            for(i=1; i<nx-1; i++) { 
				   phase2=carg(D->U[h][sliceI][j*nx+(i+1)]);
					phase1=carg(D->U[h][sliceI][j*nx+(i-1)]);
					if(phase2-phase1>M_PI)       phase2 -= 2*M_PI;
					else if(phase2-phase1<-M_PI) phase2 += 2*M_PI;
	      	   val=(phase2-phase1)*0.5;
	      	   data[start+h]+=fabs(val);
				}
         for(i=1; i<nx-1; i++)  
            for(j=1; j<ny-1; j++)  {
	      	   phase2=carg(D->U[h][sliceI][(j+1)*nx+i]);
					phase1=carg(D->U[h][sliceI][(j-1)*nx+i]);
					if(phase2-phase1>M_PI)       phase2 -= 2*M_PI;
					else if(phase2-phase1<-M_PI) phase2 += 2*M_PI;
	      	   val=(phase2-phase1)*0.5;
	      	   data[start+h+1]+=fabs(val);
				}
		}
		start+=2;
	}   

   //center divergence
	for(h=0; h<numHarmony; h++)  {	
		if(D->minI<=cenI && cenI<D->maxI) {
         sliceI=cenI+1-minI;
      
         for(j=1; j<ny-1; j++)  
            for(i=1; i<nx-1; i++)  {
				   phase2=carg(D->U[h][sliceI][j*nx+(i+1)]);
					phase1=carg(D->U[h][sliceI][j*nx+(i-1)]);
					if(phase2-phase1>M_PI)       phase2 -= 2*M_PI;
					else if(phase2-phase1<-M_PI) phase2 += 2*M_PI;
	      	   val=(phase2-phase1)*0.5;
	      	   data[start+h]+=fabs(val);
				}
         for(i=1; i<nx-1; i++)  
            for(j=1; j<ny-1; j++)  {
	      	   phase2=carg(D->U[h][sliceI][(j+1)*nx+i]);
					phase1=carg(D->U[h][sliceI][(j-1)*nx+i]);
					if(phase2-phase1>M_PI)       phase2 -= 2*M_PI;
					else if(phase2-phase1<-M_PI) phase2 += 2*M_PI;
	      	   val=(phase2-phase1)*0.5;
	      	   data[start+h+1]+=fabs(val);
				}
		} else ;
		start+=2;
	}

   dataNum=D->numHarmony*2*2;
   for(i=0; i<dataNum; i++) sendData[i] = data[i];

   for(rank=1; rank<nTasks; rank++) {
      if(myrank==rank) MPI_Send(sendData,dataNum,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD); 
		else ;
	}

   if(myrank==0) {
      for(rank=1; rank<nTasks; rank++) {
	      MPI_Recv(recvData,dataNum,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
         for(i=0; i<dataNum; i++) data[i] += recvData[i];
      }

      out=fopen("divergence","a+");
      z=iteration*dz+minZ;
      fprintf(out,"%.15g ",z);
	   start=0;
		divX=1.0/(D->sliceN*dx*(nx-2)*(ny-2));
		divY=1.0/(D->sliceN*dy*(nx-2)*(ny-2));
      for(h=0; h<numHarmony; h++) {
         fprintf(out,"%g %g ",data[start+h]*divX*1e-6,data[start+h+1]*divY*1e-6);
		   start+=2;
	   }
		divX=1.0/(dx*(nx-2)*(ny-2));
		divY=1.0/(dy*(nx-2)*(ny-2));
      for(h=0; h<numHarmony; h++) {
         fprintf(out,"%g %g ",data[start+h]*divX,data[start+h+1]*divY);
		   start+=2;
	   }
      fprintf(out,"\n");
      fclose(out);
   } else ;

   MPI_Barrier(MPI_COMM_WORLD);
}


void updateBFactor(Domain *D,int iteration)
{
   int n,i,h,s,startI,endI,minI,sliceI,sliceN,rank,N[D->nSpecies],numInBeamlet;
	int cenI;
	double bucketZ,z,dz,minZ,theta,cnt,sumI,sumR;
	double dataB[6],sendData[6],recvData[6];
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
	sliceN=D->sliceN;
   dz=D->dz;
   minZ=D->minZ;
	cenI=sliceN/2;


   LL=D->loadList;
   s=0;
   while(LL->next) {
      N[s]=LL->numInBeamlet;
      LL=LL->next;
      s++;
   }

   // total bunching
   sum=0.0+I*0.0;
	cnt=0.0;
   for(i=startI; i<endI; i++)
   {
     for(s=0; s<D->nSpecies; s++)  {
	    numInBeamlet=N[s];
       p=D->particle[i].head[s]->pt;
       while(p) {
         for(n=0; n<numInBeamlet; n++) {
			  sum+=cexp(I*p->theta[n]);
			  cnt+=1.0;
         }
			p=p->next;
		 }
	  }
	}
   dataB[0]=creal(sum);
   dataB[1]=cimag(sum);
   dataB[2]=cnt;

   // center Bunching factor
   sum=0.0+I*0.0;
	cnt=0.0;
	if(D->minI<=cenI && cenI<D->maxI) {
      i=cenI+1-minI;
      for(s=0; s<D->nSpecies; s++)  {
	      numInBeamlet=N[s];
         p=D->particle[i].head[s]->pt;
         while(p) {
            for(n=0; n<numInBeamlet; n++) {
			      sum+=cexp(I*p->theta[n]);
			      cnt+=1.0;
            }
		  	   p=p->next;
		   }
	   }
	}
   dataB[3]=creal(sum);
   dataB[4]=cimag(sum);
   dataB[5]=cnt;

   MPI_Barrier(MPI_COMM_WORLD);

   for(h=0; h<6; h++) sendData[h] = dataB[h];

   for(i=1; i<nTasks; i++) {
      if(myrank==i) MPI_Send(sendData,6,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD); else ;
	}

   if(myrank==0) {
      for(i=1; i<nTasks; i++) {
	      MPI_Recv(recvData,6,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
         for(h=0; h<6; h++) dataB[h] += recvData[h];
      }

      out=fopen("bunching","a+");
      z=iteration*dz+minZ;
      fprintf(out,"%.15g %g %g\n",z,sqrt(dataB[0]*dataB[0]+dataB[1]*dataB[1])/dataB[2],sqrt(dataB[3]*dataB[3]+dataB[4]*dataB[4])/dataB[5]);
      fclose(out);
   } else ;

   MPI_Barrier(MPI_COMM_WORLD);
}


/*
void saveBFactor(Domain *D,int iteration)
{
   int i,j,s,startI,endI,minI,sliceI,sliceN,rank,cnt;
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

	sliceN=D->sliceN;
	recvData=(double *)malloc((sliceN+2)*sizeof(double ));
	data=(double *)malloc((sliceN+2)*sizeof(double ));
	for(i=0; i<(sliceN+2); i++) data[i]=0.0;

	
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
	    MPI_Recv(recvData,sliceN+2,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
       for(i=0; i<sliceN+2; i++)
         data[i]+=recvData[i];
     }
	} else {
	  for(rank=1; rank<nTasks; rank++) {
	    if(myrank==rank)
		   MPI_Send(data,sliceN+2,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
     }
   }

   bucketZ=D->numSlice*D->lambda0;
   if(myrank==0) {
     sprintf(name,"bFact%d",iteration);
     out = fopen(name,"w");

     for(i=1; i<sliceN+1; i++) {
       z=(i-1)*bucketZ;
       fprintf(out,"%g %g\n",z,data[i]);
     }
	  fclose(out);
     printf("%s is made.\n",name);
   } else ;

}
*/
