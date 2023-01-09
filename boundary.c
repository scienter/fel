#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <mpi.h>
#include <complex.h>
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
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   double *number,minZ,n0,sumDouble,tmpDouble,dZ,posZ,xi;
   int sliceN,max,l,subCnt,remain,tmpN,minN,maxN;

   dZ=D->lambda0*D->numSlice;
   minZ=D->minZ;
   sliceN=D->sliceN;
   sumDouble=0;

   if(D->dimension==1) {
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

     free(number);

   } else if(D->dimension==3) {
	  subCnt=D->sliceN/nTasks;
     remain=D->sliceN%nTasks;
     minN=maxN=0;
     for(rank=0; rank<nTasks; rank++) {
       if(rank<remain) tmpN=subCnt+1;
       else            tmpN=subCnt;
       minN=maxN;
       maxN=minN+tmpN;
       if(myrank==rank) {
         D->minI=minN;
         D->maxI=maxN;
	 D->subSliceN=tmpN;
       } else ;
     }
   }


   MPI_Barrier(MPI_COMM_WORLD);
   printf("myrank=%d,minI=%d,maxI=%d,subSliceN=%d\n",myrank,D->minI,D->maxI,D->subSliceN);
  


   // Field memory setting
   D->U=complexMemory3Asign(D->numHarmony,D->subSliceN+2,D->nx,D->ny);
   D->Uc=complexMemory3Asign(D->numHarmony,D->subSliceN+2,D->nx,D->ny);
   D->ScU=complexMemory3Asign(D->numHarmony,D->subSliceN+2,D->nx,D->ny);
   D->totalEnergy=(double **)malloc(D->maxStep*sizeof(double *));
   for(i=0; i<D->maxStep; i++)
     D->totalEnergy[i]=(double *)malloc(D->numHarmony*sizeof(double ));

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
   D->slope=complexMemory3Asign(D->numHarmony,D->subSliceN+2,D->nx,D->ny);
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
//   tmpM=(double complex *)malloc(N*sizeof(double complex ));
/*
   D->Ma=complexMemoryAsign(D->harmony,D->nx,D->nx);
   D->Mb=complexMemoryAsign(D->harmony,D->ny,D->ny);
   D->invMa=complexMemoryAsign(D->harmony,D->nx,D->nx);
   D->invMb=complexMemoryAsign(D->harmony,D->ny,D->ny);

//   D->tmpU=(double complex *)malloc(D->nx*D->ny*sizeof(double complex ));

   if(D->dimension==3)  {
     for(h=0; h<D->harmony; h++) {

     //----------------- calculate Ma------------------//
       N=nx*nx;
       alpha[h]=-I*D->dz*0.25/(D->ks*(h+1))/D->dx/D->dx;
       for(i=0; i<N; i++) D->Ma[h][i]=0.0+I*0.0;
       for(i=1; i<nx-1; i++) {
         D->Ma[h][i*(nx+1)+0]=1.0-2*alpha[h];
         D->Ma[h][i*(nx+1)+1]=alpha[h];
         D->Ma[h][i*(nx+1)-1]=alpha[h];
       }
       i=0;
         D->Ma[h][i*(nx+1)+0]=1.0-2*alpha[h];
         D->Ma[h][i*(nx+1)+1]=alpha[h];
       i=nx-1;
         D->Ma[h][i*(nx+1)+0]=1.0-2*alpha[h];
         D->Ma[h][i*(nx+1)-1]=alpha[h];

       gsl_matrix_complex *Ma = gsl_matrix_complex_alloc(nx,nx);

       for(i=0; i<nx; i++)
         for(j=0; j<nx; j++) {
           data=D->Ma[h][i*nx+j];
           GSL_SET_COMPLEX(&g,creal(data),cimag(data));
           gsl_matrix_complex_set(Ma,i,j,g);
         }

       // Inverse Matrix Ma
       gsl_permutation *pa = gsl_permutation_alloc(nx);
       gsl_linalg_complex_LU_decomp(Ma,pa,&stat);
       gsl_matrix_complex *invA = gsl_matrix_complex_alloc(nx,nx);
       gsl_linalg_complex_LU_invert(Ma,pa,invA);

       for(i=0; i<nx; i++)
         for(j=0; j<nx; j++) {
           g=gsl_matrix_complex_get(invA,i,j);
           D->invMa[h][i*nx+j]=g.dat[0]+I*g.dat[1];
         }

       //----------------- calculate Ma------------------//
       N=ny*ny;
       beta[h]=-I*D->dz*0.25/(D->ks*(h+1))/D->dy/D->dy;
       for(i=1; i<ny-1; i++) {
         D->Mb[h][i*(ny+1)+0]=1.0-2*beta[h];
         D->Mb[h][i*(ny+1)+1]=beta[h];
         D->Mb[h][i*(ny+1)-1]=beta[h];
       }
       i=0;
         D->Mb[h][i*(ny+1)+0]=1.0-2*beta[h];
         D->Mb[h][i*(ny+1)+1]=beta[h];
       i=nx-1;
         D->Mb[h][i*(ny+1)+0]=1.0-2*beta[h];
         D->Mb[h][i*(ny+1)-1]=beta[h];

       gsl_matrix_complex *Mb = gsl_matrix_complex_alloc(ny,ny);

       for(i=0; i<ny; i++)
         for(j=0; j<ny; j++) {
           data=D->Mb[h][i*ny+j];
           GSL_SET_COMPLEX(&g,creal(data),cimag(data));
           gsl_matrix_complex_set(Mb,i,j,g);
         }

       // Inverse Matrix Mb
       gsl_permutation *pb = gsl_permutation_alloc(ny);
       gsl_linalg_complex_LU_decomp(Mb,pb,&stat);
       gsl_matrix_complex *invB = gsl_matrix_complex_alloc(ny,ny);
       gsl_linalg_complex_LU_invert(Mb,pb,invB);

       for(i=0; i<ny; i++)
         for(j=0; j<ny; j++) {
           g=gsl_matrix_complex_get(invB,i,j);
           D->invMb[h][i*ny+j]=g.dat[0]+I*g.dat[1];
         }

// test matrix     
//     for(i=0; i<nx; i++) {
//       for(j=0; j<nx; j++) {
//         sum=0.0+I*0.0;
//         for(nn=0; nn<nx; nn++)
//           sum+=D->Ma[h][i*nx+nn]*D->invMa[h][nn*nx+j];
//         tmpM[i*nx+j]=sum;
//       }
//     }
//     out=fopen("matrix","w");
//     for(i=0; i<nx; i++) {
//       for(j=0; j<nx; j++) {
//         data=tmpM[i*nx+j];
//         fprintf(out,"%4.2g ",creal(data));
//       }
//       fprintf(out,"\n");
//     }
//     fclose(out);

       gsl_permutation_free(pa);
       gsl_matrix_complex_free(Ma);
       gsl_matrix_complex_free(invA);
       gsl_permutation_free(pb);
       gsl_matrix_complex_free(Mb);
       gsl_matrix_complex_free(invB);
     }
   } else ; 	//End of if(D->dimension==3)
*/

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
