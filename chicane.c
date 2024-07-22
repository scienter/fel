#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"

void calParticleDelay(Domain *D,int iteration)
{
	double B0,ld,L1,L2,vz,gamma0,ldbyR,th0,numSlice;
	double dz1,dz2,dz3,shiftZ0,dPhi;
	int sliceI,startI,endI,s,nSpecies,subSliceN,n,*N;
	double gamma,invGam,shiftZ,ks;
   LoadList *LL;
   ptclList *p;
    int myrank, nTasks;

    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   ks=D->ks;
	B0=D->dipoleB;
	ld=D->ld;
	L1=D->L1;
	L2=D->L2;
   gamma0=D->gamma0;
	vz=velocityC*sqrt(1.0-1.0/gamma0/gamma0);
	numSlice=D->numSlice;
	dPhi=D->numSlice*2*M_PI;

	ldbyR=eCharge*ld*B0/gamma0/eMass/vz;
	th0=atan(ldbyR+0.5*ldbyR*ldbyR*ldbyR);
   dz1=gamma0*eMass*vz/eCharge/B0*th0-ld;
	dz2=(L1-ld)*(1.0/cos(th0)-1.0);
	dz3=(L2-ld)*(1.0-vz/velocityC);
	shiftZ0=(4.0*dz1+2.0*dz2+dz3);

	D->shiftSlice=(int)(shiftZ0/(D->lambda0*D->numLambdaU));
   
   startI=1;  endI=D->subSliceN+1;
	subSliceN=D->subSliceN;

   N=(int *)malloc(D->nSpecies*sizeof(int ));
   LL=D->loadList;
   s=0;
   while(LL->next) {
      N[s]=LL->numInBeamlet;
      LL=LL->next;
      s++;
   }

	nSpecies=D->nSpecies;
   for(sliceI=startI; sliceI<endI; sliceI++)
   {    
     for(s=0; s<nSpecies; s++)  {
       p=D->particle[sliceI].head[s]->pt;
       while(p) {
		   for(n=0; n<N[s]; n++) {
           gamma=p->gamma[n];  invGam=1.0/gamma;

	        vz=velocityC*sqrt(1.0-invGam*invGam);
	        ldbyR=eCharge*ld*B0*invGam/eMass/vz;
	        th0=atan(ldbyR+0.5*ldbyR*ldbyR*ldbyR);
           dz1=gamma*eMass*vz/eCharge/B0*th0-ld;
           dz2=(L1-ld)*(1.0/cos(th0)-1.0);
           dz3=(L2-ld)*(1.0-vz/velocityC);
           shiftZ=4.0*dz1+2.0*dz2+dz3;
	      
           p->theta[n]+=(shiftZ0-shiftZ)*ks;
         }
         p=p->next;
       }
     }		//End of for(s)
   }		//End of for(sliceI)
}

void rearrangeChicaneParticle(Domain *D)
{
    Particle *particle;
    particle=D->particle;

    int i,s,ii,intZ,cnt,deleteFlag,shiftFlag,domainShiftFlag,dataCnt;
    int l,rank,startI,endI,nSpecies,minI,maxI,indexI,n,*N;
    double dPhi,theta,aveTh,delTh;
	 int *sendN,*recvN,*start,*cntN;
	 double **sendData,**recvData;
    LoadList *LL;	 
    ptclList *p,*New,*prev,*tmp;
    int myrank, nTasks;
    double a,b,c,d,e,f;

    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    startI=1;  endI=D->subSliceN+1;
    minI=D->minI;  maxI=D->maxI;
    nSpecies=D->nSpecies;
    dPhi=D->numSlice*2*M_PI;

    sendN=(int *)malloc(nTasks*sizeof(int ));
    recvN=(int *)malloc(nTasks*sizeof(int ));

    N=(int *)malloc(D->nSpecies*sizeof(int ));
    LL=D->loadList;
    s=0;
    while(LL->next) {
       N[s]=LL->numInBeamlet;
       LL=LL->next;
       s++;
    }

    for(s=0; s<nSpecies; s++) 
	 {
       dataCnt=6*N[s]+4;
       for(i=0; i<nTasks; i++) { sendN[i]=0; recvN[i]=0; }

	    for(i=startI; i<endI; i++)   {
          p=particle[i].head[s]->pt;
	       while(p)  {
			    aveTh=0.0;
			    for(n=0; n<N[s]; n++) aveTh+=p->theta[n];
				 aveTh/=1.0*N[s];

      	    if(aveTh>=dPhi)  intZ=(int)(aveTh/dPhi);
         	 else if(aveTh<0) intZ=((int)(aveTh/dPhi))-1;
	          else             intZ=0;

             intZ*=-1;
   	       indexI=i-startI+minI+intZ;
      	    for(rank=0; rank<nTasks; rank++) {
         	   if(D->minmax[rank]<=indexI && indexI<D->minmax[rank+1]) {
					  if(rank!=myrank)  sendN[rank]+=1; else ;
					  rank=nTasks;
					} else ;
				 }
			    p=p->next;
	       }	//End of while(p)
       }		//End of for(s)

       for(rank=0; rank<nTasks; rank++)   {
          if(myrank!=rank)
             MPI_Send(&sendN[rank],1,MPI_INT,rank,myrank,MPI_COMM_WORLD);
	    }
	    for(rank=0; rank<nTasks; rank++)   {
	       if(myrank!=rank)    {
	          MPI_Recv(&recvN[rank],1,MPI_INT,rank,rank,MPI_COMM_WORLD,&status);
	       }  else ;
	    }
	    MPI_Barrier(MPI_COMM_WORLD);

       recvData=(double **)malloc(nTasks*sizeof(double *));
       sendData=(double **)malloc(nTasks*sizeof(double *));
       for(rank=0; rank<nTasks; rank++)   {
          sendData[rank]=(double *)malloc(sendN[rank]*dataCnt*sizeof(double ));
          recvData[rank]=(double *)malloc(recvN[rank]*dataCnt*sizeof(double ));
       }
       
       cntN=(int *)malloc(nTasks*sizeof(int ));
       for(rank=0; rank<nTasks; rank++) cntN[rank]=0; 
       for(i=startI; i<endI; i++)
       {
          cnt=1;
          p=particle[i].head[s]->pt;
          while(p)  {
             if(cnt==1)  prev=p; else ;
              
			    deleteFlag=OFF;
			    shiftFlag=OFF;
			    domainShiftFlag=OFF;

			    aveTh=0.0;
			    for(n=0; n<N[s]; n++) aveTh+=p->theta[n];
			    aveTh/=1.0*N[s];
             if(aveTh>=dPhi)  {
                intZ=(int)(aveTh/dPhi);
                delTh=dPhi*intZ;
                //theta-=dPhi*intZ;
				    shiftFlag=ON;	
             }
             else if(aveTh<0) {              
                intZ=((int)(aveTh/dPhi))-1;
                delTh=dPhi*intZ;
                //theta+=-1.0*dPhi*intZ;
				    shiftFlag=ON;	
             } 
             else  {
				    intZ=0;
					 delTh=0.0;
				 }

			    if(D->mode==Static) intZ=0; else ;

             intZ*=-1;
             indexI=i-startI+minI+intZ;
			    if(indexI<0 || indexI>D->sliceN) deleteFlag=ON; else ;

             for(rank=0; rank<nTasks; rank++) {
                if(D->minmax[rank]<=indexI && indexI<D->minmax[rank+1]) {
				       if(myrank!=rank) {
						    domainShiftFlag=ON;
				          for(n=0; n<N[s]; n++) {
                         sendData[rank][cntN[rank]*dataCnt+n*6+0]=p->theta[n]-delTh;
                         sendData[rank][cntN[rank]*dataCnt+n*6+1]=p->x[n];
                         sendData[rank][cntN[rank]*dataCnt+n*6+2]=p->y[n];
                         sendData[rank][cntN[rank]*dataCnt+n*6+3]=p->gamma[n];
                         sendData[rank][cntN[rank]*dataCnt+n*6+4]=p->px[n];
                         sendData[rank][cntN[rank]*dataCnt+n*6+5]=p->py[n];
				          }
                      sendData[rank][cntN[rank]*dataCnt+N[s]*6+0]=p->weight;
                      sendData[rank][cntN[rank]*dataCnt+N[s]*6+1]=p->index;
                      sendData[rank][cntN[rank]*dataCnt+N[s]*6+2]=p->core;
                      sendData[rank][cntN[rank]*dataCnt+N[s]*6+3]=indexI;
                      cntN[rank]+=1;
					       deleteFlag=ON;
				       } else ;
				       rank=nTasks;
				    } else ; 
			    }		//End of for(n<nTasks)
               
			    if(deleteFlag==ON || domainShiftFlag==ON) {
                if(cnt==1)  {
                   particle[i].head[s]->pt = p->next;
				       free(p->x);      free(p->y);
				       free(p->px);     free(p->py);
				       free(p->theta);  free(p->gamma);
                   p->next = NULL;
		             free(p);
                   p=particle[i].head[s]->pt;
                   cnt=1;
   	          } else {
                   prev->next = p->next;
				       free(p->x);      free(p->y);
				       free(p->px);     free(p->py);
				       free(p->theta);  free(p->gamma);
				       p->next = NULL;
				       free(p);
	                p=prev->next;
				    }
			    } else if(shiftFlag==ON && domainShiftFlag==OFF) {
                if(cnt==1)  {
					    for(n=0; n<N[s]; n++) p->theta[n]-=delTh;
                   particle[i].head[s]->pt = p->next;
	                p->next = particle[i+intZ].head[s]->pt;
	                particle[i+intZ].head[s]->pt = p;
	                p=particle[i].head[s]->pt;
	                cnt=1;
	             } else {
	                prev->next = p->next;
					    for(n=0; n<N[s]; n++) p->theta[n]-=delTh;
	                p->next = particle[i+intZ].head[s]->pt;
	                particle[i+intZ].head[s]->pt = p;
	                p=prev->next;
	             }
			    } else {
			       prev=p;
			       p=p->next;
			       cnt++;
			    }	

          }	//End of while(p)
       }		//End of for(i)
       
       // domain-shifting particles
       for(rank=0; rank<nTasks; rank++)  {
          if(myrank==rank)  {
             for(l=0; l<nTasks; l++) {
                if(rank!=l)
                   MPI_Send(sendData[l],sendN[l]*dataCnt,MPI_DOUBLE,l,myrank,MPI_COMM_WORLD);
                else ;
             }
          }
          else  
			 {
             MPI_Recv(recvData[rank],recvN[rank]*dataCnt,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
				  
             for(ii=0; ii<recvN[rank]; ii++)  {
                indexI=recvData[rank][ii*dataCnt+N[s]*6+3];
                i=indexI-minI+startI;          
                New = (ptclList *)malloc(sizeof(ptclList));
                New->next = particle[i].head[s]->pt;
                particle[i].head[s]->pt = New;

                New->weight=recvData[rank][ii*dataCnt+N[s]*6+0];
                New->index=recvData[rank][ii*dataCnt+N[s]*6+1];
                New->core=recvData[rank][ii*dataCnt+N[s]*6+2];
                New->theta=(double *)malloc(N[s]*sizeof(double ));
                New->x=(double *)malloc(N[s]*sizeof(double ));
                New->y=(double *)malloc(N[s]*sizeof(double ));
                New->px=(double *)malloc(N[s]*sizeof(double ));
                New->py=(double *)malloc(N[s]*sizeof(double ));
                New->gamma=(double *)malloc(N[s]*sizeof(double ));
					 for(n=0; n<N[s]; n++) {
                   New->theta[n]=recvData[rank][ii*dataCnt+n*6+0];
                   New->x[n]=recvData[rank][ii*dataCnt+n*6+1];
                   New->y[n]=recvData[rank][ii*dataCnt+n*6+2];
                   New->gamma[n]=recvData[rank][ii*dataCnt+n*6+3];
                   New->px[n]=recvData[rank][ii*dataCnt+n*6+4];
                   New->py[n]=recvData[rank][ii*dataCnt+n*6+5];
					 }
             }
			    	 
          }
		    MPI_Barrier(MPI_COMM_WORLD);  
       }
       
       for(rank=0; rank<nTasks; rank++) {
          free(sendData[rank]);
          free(recvData[rank]);
       }
       free(sendData);
       free(recvData);
	    free(cntN);

    }  //End of for(s)
	 free(sendN);
	 free(recvN);
	 free(N);
}


void shiftChicaneField(Domain *D)
{
   int s,h,numHarmony,i,j,n,ii,nn,startI,endI,N;
	int shiftN,max,indexI,maxI,minI,num;
	int *minmax,*sendN,*recvN,*start;
	double **sendData,**recvData,realV,imagV;
	double complex tmpComp,val;
   ptclList *p;
   int myrank, nTasks;

   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	numHarmony=D->numHarmony;
   startI=1;  endI=D->subSliceN+1;
	minI=D->minI; maxI=D->maxI;
   N=D->nx*D->ny;
	shiftN=D->shiftSlice;
	printf("myrank=%d, shiftSlice=%d, sliceN=%d\n",myrank,D->shiftSlice,shiftN);

	minmax=D->minmax;
	sendN=(int *)malloc(nTasks*sizeof(int ));
	recvN=(int *)malloc(nTasks*sizeof(int ));
	for(i=0; i<nTasks; i++) { sendN[i]=0; recvN[i]=0; }

   for(i=endI-1; i>=startI; i--) {
     indexI=i-startI+minI+shiftN;
	  if(indexI>=maxI) {
       for(n=myrank+1; n<nTasks; n++) {
         if(D->minmax[n]<=indexI && indexI<D->minmax[n+1]) {
			  sendN[n]+=1;
			} else ;
	    }
	  }
	}

     for(n=0; n<nTasks; n++)   {
       if(myrank!=n)
         MPI_Send(&sendN[n],1,MPI_INT,n,myrank,MPI_COMM_WORLD);
     }
     for(n=0; n<nTasks; n++)   {
	    if(myrank!=n)    {
         MPI_Recv(&recvN[n],1,MPI_INT,n,n,MPI_COMM_WORLD,&status);
	    }  else ;
     }
     MPI_Barrier(MPI_COMM_WORLD);	  

     num=N*numHarmony*3;
	  recvData=(double **)malloc(nTasks*sizeof(double *));
     sendData=(double **)malloc(nTasks*sizeof(double *));
     for(n=0; n<nTasks; n++)   {
       sendData[n]=(double *)malloc(sendN[n]*num*sizeof(double ));
       recvData[n]=(double *)malloc(recvN[n]*num*sizeof(double ));
     }
     start=(int *)malloc(nTasks*sizeof(int ));
     for(n=0; n<nTasks; n++) start[n]=0;

     for(i=endI-1; i>=startI; i--) {
       indexI=i-startI+minI+shiftN;
		 if(indexI>=maxI) {
         for(n=myrank+1; n<nTasks; n++) {
           if(D->minmax[n]<=indexI && indexI<D->minmax[n+1]) {
             for(h=0; h<numHarmony; h++) {
               for(j=0; j<N; j++) {
                 val=D->U[h][i][j];
	 	  	        sendData[n][start[n]+0]=creal(val);
	 	  	        sendData[n][start[n]+1]=cimag(val);
	 	  	        sendData[n][start[n]+2]=indexI;
				     start[n]+=3;
				     D->U[h][i][j]=0.0+I*0.0;
				   }
			    }
			  } else ;
         }
		 } else {
         for(h=0; h<numHarmony; h++)
           for(j=0; j<N; j++) {
             D->U[h][i+shiftN][j]=D->U[h][i][j];
             D->U[h][i][j]=0.0*I*0.0;
			  }			  				 
		 } 
	  }	//End of for(i=startI; i<endI; i++)


     for(n=0; n<nTasks; n++)  {
       if(myrank==n)  {
         for(nn=0; nn<nTasks; nn++) {
           if(n!=nn)
             MPI_Send(sendData[nn],sendN[nn]*num,MPI_DOUBLE,nn,myrank,MPI_COMM_WORLD);
			  else ;
			}
       }
       else  {
         MPI_Recv(recvData[n],recvN[n]*num,MPI_DOUBLE,n,n,MPI_COMM_WORLD,&status);
			start[n]=0;
         for(ii=0; ii<recvN[n]; ii++)  {
			  indexI=recvData[n][start[n]+2];
			  i=indexI-minI+startI;
//			  if(myrank==0) printf("indexI=%d, i=%d\n",indexI,i);
           for(h=0; h<numHarmony; h++) {
             for(j=0; j<N; j++) {
	 	  	      realV=recvData[n][start[n]+0];
	 	  	      imagV=recvData[n][start[n]+1];
					D->U[h][i][j]=realV+I*imagV;
				   start[n]+=3;
			    }
			  }
			}
		 }
       MPI_Barrier(MPI_COMM_WORLD);
	  }	//End of for(n<nTasks)



     for(n=0; n<nTasks; n++) {
	    free(sendData[n]);
	    free(recvData[n]);
	  }
	  free(sendData);
	  free(recvData);
	  free(start);
	  

   free(sendN);   
   free(recvN);   
}



void chicane_test(Domain *D,int iteration)
{
	double z0,z1,x0,x1;
   ChiList *Chi;

   z0=(iteration-1.0)*D->dz;
   z1=iteration*D->dz;

   D->chicaneFlag=OFF;
   D->calChicaneFlag=OFF;
   Chi=D->chiList;
	
   while(Chi->next) {
     x0=Chi->chiStart;
     x1=Chi->chiEnd;

     if(z0<=x1 && x1<z1) {
       D->chicaneFlag=ON;
       D->calChicaneFlag=ON;
		 D->dipoleB=Chi->B0;
		 D->ld=Chi->ld;
		 D->L1=Chi->L1;
		 D->L2=Chi->L2;

       D->chi_delay=Chi->delay;
       // belows are self-seeding option
       D->chi_SSON=Chi->selfSeedON;
       D->chi_d=Chi->d;
       D->bragTh=Chi->bragTh;
       D->extincL=Chi->extincL;
       D->chi0=Chi->chi0;
       D->chi_noiseONOFF=Chi->noiseONOFF;
       D->chi_washONOFF=Chi->washONOFF;
       D->rangeE = Chi->rangeE;
       D->shiftE = Chi->shiftE;
	  } else if(x0<=z1 && z1<x1) {
       D->chicaneFlag=ON;
	  } else ;
     Chi=Chi->next; 
   }
	
}

void set_chicane_zero(Domain *D)
{
  D->dipoleB=0.0;
  D->ld=0.0;
  D->L1=0.0;
  D->L2=0.0;
  D->shiftSlice=0;
  // belows are self-seeding option
  D->chi_washONOFF=OFF;
  D->chi_noiseONOFF=ON;
}
