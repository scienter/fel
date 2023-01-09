#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include <gsl/gsl_sf_bessel.h>
#include <complex.h>
#include <string.h>

//void solveField1D(Domain *D,int sliceI,int iteration);
void solve_Sc_1D(Domain *D,int iteration);
void solve_Sc_3D(Domain *D,int iteration);
void solve_Field_U_1D(Domain *D,int iteration);
void solve_Field_U_3D(Domain *D,int iteration);
void MPI_Transfer1F_Z(double complex ***f1,int harmony,int N,int fromI,int toI);
//void doubleMPI_Transfer1F_Zplus(double ***f1,int harmony,int N,int fromI,int toI);

void solveField(Domain D,int iteration)
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
    solve_Sc_1D(&D,iteration);
    solve_Field_U_1D(&D,iteration);
    break;
  case 3:
    solve_Sc_3D(&D,iteration);
    solve_Field_U_3D(&D,iteration);
    break;

  default:
    printf("In EzSolve.c, what dimension?\n");
  }
}

void solve_Field_U_3D(Domain *D,int iteration)
{
   int h,harmony,i,j,sliceI,startI,endI;  
   int n,nx,ny;
   double ks,dx,dy,dz;
   double complex alpha,beta,later,*CC,*DD,*dd;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   startI=1;   endI=D->subSliceN+1;
   harmony=D->harmony;

   dx=D->dx;  dy=D->dy;  dz=D->dz;
   ks=D->ks;
   nx=D->nx;  ny=D->ny;

   // first step
   CC=(double complex *)malloc(nx*sizeof(double complex));
   DD=(double complex *)malloc(nx*sizeof(double complex));
   dd=(double complex *)malloc(nx*sizeof(double complex));

   for(h=0; h<harmony; h++)  {
     alpha=-I*dz*0.25/(ks*(h+1))/dx/dx;
     beta=-I*dz*0.25/(ks*(h+1))/dy/dy;
     CC[0]=alpha/(1.0-2*alpha);
     for(sliceI=startI; sliceI<endI; sliceI++) {
//       memcpy(&(D->tmpU[0]),&(D->U[h][sliceI][0]),nx*ny*sizeof(double complex ));
       for(j=1; j<ny-1; j++) {
         for(i=0; i<nx; i++)
           dd[i]=-beta*D->U[h][sliceI][(j+1)*nx+i]+(1+2*beta)*D->U[h][sliceI][j*nx+i]-beta*D->U[h][sliceI][(j-1)*nx+i]+D->ScU[h][sliceI][j*nx+i];
	 DD[0]=dd[0]/(1.0-2*alpha);
         for(i=1; i<nx; i++) {
           CC[i]=alpha/(1-2*alpha-alpha*CC[i-1]);
           DD[i]=(dd[i]-alpha*DD[i-1])/(1-2*alpha-alpha*CC[i-1]);
	 }
	 i=nx-1;
	 later=DD[i];
	 D->Uc[h][sliceI][j*nx+i]=later;
	 for(i=nx-2; i>=0; i--) {
	   later=DD[i]-CC[i]*later;
	   D->Uc[h][sliceI][j*nx+i]=later;
	 }
       }
       j=0;
         for(i=0; i<nx; i++) {
           dd[i]=-beta*D->U[h][sliceI][(j+1)*nx+i]+(1+2*beta)*D->U[h][sliceI][j*nx+i]+D->ScU[h][sliceI][j*nx+i];
//           dd[i]=D->tmpU[j*nx+i]+D->ScU[h][sliceI][j*nx+i];
	 }
	 DD[0]=dd[0]/(1.0-2*alpha);
         for(i=1; i<nx; i++) {
           CC[i]=alpha/(1-2*alpha-alpha*CC[i-1]);
           DD[i]=(dd[i]-alpha*DD[i-1])/(1-2*alpha-alpha*CC[i-1]);
	 }
	 i=nx-1;
	 later=DD[i];
	 D->Uc[h][sliceI][j*nx+i]=later;
	 for(i=nx-2; i>=0; i--) {
	   later=DD[i]-CC[i]*later;
	   D->Uc[h][sliceI][j*nx+i]=later;
	 }
       j=ny-1;
         for(i=0; i<nx; i++) {
           dd[i]=(1+2*beta)*D->U[h][sliceI][j*nx+i]-beta*D->U[h][sliceI][(j-1)*nx+i]+D->ScU[h][sliceI][j*nx+i];
//           dd[i]=D->tmpU[j*nx+i]+D->ScU[h][sliceI][j*nx+i];
         }
	 DD[0]=dd[0]/(1.0-2*alpha);
         for(i=1; i<nx; i++) {
           CC[i]=alpha/(1-2*alpha-alpha*CC[i-1]);
           DD[i]=(dd[i]-alpha*DD[i-1])/(1-2*alpha-alpha*CC[i-1]);
	 }
	 i=nx-1;
	 later=DD[i];
	 D->Uc[h][sliceI][j*nx+i]=later;
	 for(i=nx-2; i>=0; i--) {
	   later=DD[i]-CC[i]*later;
	   D->Uc[h][sliceI][j*nx+i]=later;
	 }
     }
   }
   free(CC); 
   free(DD);
   free(dd);

   // second step
   CC=(double complex *)malloc(ny*sizeof(double complex));
   DD=(double complex *)malloc(ny*sizeof(double complex));
   dd=(double complex *)malloc(ny*sizeof(double complex));

   for(h=0; h<harmony; h++)  {
     alpha=-I*dz*0.25/(ks*(h+1))/dx/dx;
     beta=-I*dz*0.25/(ks*(h+1))/dy/dy;
     CC[0]=beta/(1-2*beta);
     for(sliceI=startI; sliceI<endI; sliceI++) {
//       memcpy(&(D->tmpU[0]),&(D->U[h][sliceI][0]),nx*ny*sizeof(double complex ));
       for(i=1; i<nx-1; i++) {
         for(j=0; j<ny; j++) 
           dd[j]=-alpha*D->Uc[h][sliceI][j*nx+(i+1)]+(1+2*alpha)*D->Uc[h][sliceI][j*nx+i]-alpha*D->Uc[h][sliceI][j*nx+(i-1)]+D->ScU[h][sliceI][j*nx+i];
	 DD[0]=dd[0]/(1-2*beta);
         for(j=1; j<ny; j++) {
           CC[j]=beta/(1-2*beta-beta*CC[j-1]);
           DD[j]=(dd[j]-beta*DD[j-1])/(1-2*beta-beta*CC[j-1]);
	 }
	 j=ny-1;
	 later=DD[j];
	 D->U[h][sliceI][j*nx+i]=later;
	 for(j=ny-2; j>=0; j--) {
	   later=DD[j]-CC[j]*later;
	   D->U[h][sliceI][j*nx+i]=later;
	 }
       }
       i=0;
         for(j=0; j<ny; j++) {
           dd[j]=-alpha*D->Uc[h][sliceI][j*nx+(i+1)]+(1+2*alpha)*D->Uc[h][sliceI][j*nx+i]+D->ScU[h][sliceI][j*nx+i];
//           dd[j]=D->tmpU[j*nx+i]+D->ScU[h][sliceI][j*nx+i];
	 }
	 DD[0]=dd[0]/(1-2*beta);
         for(j=1; j<ny; j++) {
           CC[j]=beta/(1-2*beta-beta*CC[j-1]);
           DD[j]=(dd[j]-beta*DD[j-1])/(1-2*beta-beta*CC[j-1]);
	 }
	 j=ny-1;
	 later=DD[j];
	 D->U[h][sliceI][j*nx+i]=later;
	 for(j=ny-2; j>=0; j--) {
	   later=DD[j]-CC[j]*later;
	   D->U[h][sliceI][j*nx+i]=later;
	 }
       i=nx-1;
         for(j=0; j<ny; j++) {
           dd[j]=(1+2*alpha)*D->Uc[h][sliceI][j*nx+i]-alpha*D->Uc[h][sliceI][j*nx+(i-1)]+D->ScU[h][sliceI][j*nx+i];
//           dd[j]=D->tmpU[j*nx+i]+D->ScU[h][sliceI][j*nx+i];
	 }
	 DD[0]=dd[0]/(1-2*beta);
         for(j=1; j<ny; j++) {
           CC[j]=beta/(1-2*beta-beta*CC[j-1]);
           DD[j]=(dd[j]-beta*DD[j-1])/(1-2*beta-beta*CC[j-1]);
	 }
	 j=ny-1;
	 later=DD[j];
	 D->U[h][sliceI][j*nx+i]=later;
	 for(j=ny-2; j>=0; j--) {
	   later=DD[j]-CC[j]*later;
	   D->U[h][sliceI][j*nx+i]=later;
	 }
     }
   }
   free(CC); 
   free(DD);
   free(dd);

}

void solve_Sc_3D(Domain *D,int iteration)
{
   int sliceI,i,j,ii,jj,s,h,harmony,order,nx,ny,N;  
   int startI,endI;
   double coef,tmp,J1,J2,K,K0,xi,macro,invG,area; 
   double kx,ky,x,y,dx,dy,dz,theta,minX,minY,wx[2],wy[2];
   double complex macro_K_invG_expTheta_coef,tmpComp;
   double JJ[D->harmony];
   ptclList *p;
   LoadList *LL;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   startI=1; endI=D->subSliceN+1;
   harmony=D->harmony;

   minX=D->minX;	minY=D->minY;
   nx=D->nx;   ny=D->ny;
   dx=D->dx;   dy=D->dy;   dz=D->dz;
   K0=D->K0;
   kx=D->kx; ky=D->ky;

   coef=dz*eCharge*eCharge*mu0*0.25/eMass/D->ks/(D->lambda0*D->numSlice)/dx/dy;

   N=nx*ny;   
   for(h=0; h<harmony; h++)
     for(i=0; i<=endI; i++)
       for(j=0; j<N; j++) {
         D->ScU[h][i][j]=0.0+I*0.0;
         D->ScEz[h][i][j]=0.0+I*0.0;
       }


   for(sliceI=startI; sliceI<endI; sliceI++)
   {
     LL=D->loadList;
     s=0;
     while(LL->next) {
       p=D->particle[sliceI].head[s]->pt;
       while(p) {
         x=p->x;  y=p->y;   theta=p->theta;
         macro=p->weight;   invG=1.0/p->gamma;

         K=K0*(1.0+kx*kx*x*x+ky*ky*y*y);
         xi=K*K*0.5/(1+K*K);
         i=(int)((x-minX)/dx);
         j=(int)((y-minY)/dy);
         wx[1]=(x-minX)/dx-i;   wx[0]=1.0-wx[1];
         wy[1]=(y-minY)/dy-j;   wy[0]=1.0-wy[1];	  
         macro_K_invG_expTheta_coef=macro*K*coef*cexp(-I*theta)*invG;
         if(i>=0 && i<nx-1 && j>=0 && j<ny-1)  {
           for(h=1; h<=harmony; h++)  {
             if(h%2==1)  {  //odd harmony
               tmp=pow(-1.0,(h-1)/2);
  	       order=(h-1)/2;
	       J1=gsl_sf_bessel_Jn(order,h*xi);
	       order=(h+1)/2;
	       J2=gsl_sf_bessel_Jn(order,h*xi);
	       JJ[h-1]=tmp*(J1-J2);
	     } else {
	       JJ[h-1]=0.0;
             }

             tmpComp=JJ[h-1]*macro_K_invG_expTheta_coef;
             D->ScU[h-1][sliceI][j*nx+i]+=tmpComp;
//	     if(isnan(creal(tmpComp))) { printf("Sc real=%g,Sc imag=%g\n",creal(tmpComp),cimag(tmpComp)); exit(0); }
// 	     for(ii=0; ii<2; ii++)
//               for(jj=0; jj<2; jj++) {
//                 D->ScU[h-1][sliceI][(j+jj)*nx+(i+ii)]+=wx[ii]*wy[jj]*JJ[h-1]*macro_K_invG_expTheta_coef;
//               }
	   }		//End of harmony
         } else ;	//End of for(i,j)

         p=p->next;
       }

       s++;
       LL=LL->next;
     }
   }		//End of for(sliceI)
}


void solve_Field_U_1D(Domain *D,int iteration)
{
   double Kr,K0;
   int h,harmony,i,startI,endI;

   harmony=D->harmony;
   startI=1;  endI=D->subSliceN+1;

   // field update
   for(h=0; h<harmony; h++)  {
     for(i=startI; i<endI; i++) {
       D->U[h][i][0]=D->U[h][i][0]+D->ScU[h][i][0]*D->currentFlag;
       D->Ez[h][i][0]=D->ScEz[h][i][0];
     }
   }
}

void solve_Sc_1D(Domain *D,int iteration)
{
   int i,s,h,harmony,order,n,step;
   int startI,endI;  
   double coef1,coef2,tmp,J1,J2,K,Kr,K0,xi,macro; 
   double dz,theta,area;
   double complex factor1,factor2,alpha,beta,sum,tmpc;
   double JJ[D->harmony];
   ptclList *p;
   LoadList *LL;
   int myrank, nTasks,rank;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   startI=1; endI=D->subSliceN+1;
   dz=D->dz;
   harmony=D->harmony;
   coef1=dz*eCharge*eCharge*mu0*0.5/eMass/D->ks/(D->lambda0*D->numSlice);
   coef2=eCharge*mu0*velocityC*velocityC/D->ks/(D->lambda0*D->numSlice);
   K0=D->K0;
   
   for(h=0; h<harmony; h++)
     for(i=0; i<endI+1; i++) {
       D->ScU[h][i][0]=0.0+I*0.0;
       D->ScEz[h][i][0]=0.0+I*0.0;
     }


   for(i=startI; i<endI; i++)
   {
     LL=D->loadList;
     s=0;
     while(LL->next) {
       area=M_PI*LL->sigX*LL->sigY;
//       area=M_PI*D->spotSigR*D->spotSigR;
       p=D->particle[i].head[s]->pt;
       while(p) {
         theta=p->theta;
         macro=p->weight;

         K=K0;
         xi=K*K*0.5/(1+K*K);
         factor1=macro*K*coef1*cexp(-I*theta)/p->gamma/area;
         for(h=1; h<=harmony; h++)  {
           factor2=macro*coef2*cexp(-h*I*theta)/area;

           if(h%2==1)  {  //odd harmony
             tmp=pow(-1.0,(h-1)/2);
             order=(h-1)/2;
             J1=gsl_sf_bessel_Jn(order,h*xi);
             order=(h+1)/2;
             J2=gsl_sf_bessel_Jn(order,h*xi);
             JJ[h-1]=tmp*(J1-J2);
           } else {
                   JJ[h-1]=0.0;
           }

	   D->ScU[h-1][i][0]+=JJ[h-1]*factor1;
	   D->ScEz[h-1][i][0]+=-I*factor2/(h*1.0);
	 }		//End of harmony

         p=p->next;
       }	//End of while(p)

       s++;
       LL=LL->next;
     }	
   }		//End of for(i)

}

