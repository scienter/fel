#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"

void calParticleDelay(Domain *D,int iteration)
{
	double B0,ld,L1,L2,vz,gamma0,ldbyR,th0,numSlice;
	double dz1,dz2,dz3,shiftZ0,dXi,dPhi;
	int sliceI,startI,endI,s,nSpecies,subSliceN;
	double gamma,invGam,shiftZ,tmp;
   ptclList *p;

	B0=D->dipoleB;
	ld=D->ld;
	L1=D->L1;
	L2=D->L2;
   gamma0=D->gamma0;
	vz=velocityC*sqrt(1.0-1.0/gamma0/gamma0);
	dXi=D->numSlice*D->lambda0;
	dPhi=D->numSlice*2*M_PI;

	ldbyR=eCharge*ld*B0/gamma0/eMass/vz;
	th0=atan(ldbyR+0.5*ldbyR*ldbyR*ldbyR);
   dz1=gamma0*eMass*vz/eCharge/B0*th0-ld;
	dz2=(L1-ld)*(1.0/cos(th0)-1.0);
	dz3=(L2-ld)*(1.0-vz/velocityC);
	shiftZ0=4.0*dz1+2.0*dz2+dz3;

	D->shiftSlice=(int)(shiftZ0);
   
   startI=1;  endI=D->subSliceN+1;
	subSliceN=D->subSliceN;

	nSpecies=D->nSpecies;
   for(sliceI=startI; sliceI<endI; sliceI++)
   {    
     for(s=0; s<nSpecies; s++)  {
       p=D->particle[sliceI].head[s]->pt;
       while(p) {
         gamma=p->gamma;  invGam=1.0/gamma;

	      vz=velocityC*sqrt(1.0-invGam*invGam);
	      ldbyR=eCharge*ld*B0*invGam/eMass/vz;
	      th0=atan(ldbyR+0.5*ldbyR*ldbyR*ldbyR);
         dz1=gamma*eMass*vz/eCharge/B0*th0-ld;
         dz2=(L1-ld)*(1.0/cos(th0)-1.0);
         dz3=(L2-ld)*(1.0-vz/velocityC);
         shiftZ=4.0*dz1+2.0*dz2+dz3;
	
         tmp=(shiftZ-shiftZ0)/D->lambda0*2*M_PI;
         //printf("gamma0=%g, gamma=%g, tmp=%g\n",gamma0,gamma,tmp);
         p->theta+=tmp;
         

			//if(fabs(tmp/D->numSlice)>=subSliceN) {
			//	printf("shiftSlice=%g, subSliceN=%d it is too much.\n",tmp/D->numSlice,subSliceN);
			//   exit(0);
			//}  else ;

         p=p->next;
       }
     }		//End of for(s)
   }		//End of for(sliceI)
}

void rearrangeChicaneParticle(Domain *D)
{
    Particle *particle;
    particle=D->particle;

    int i,s,n,ii,nn,intZ,cnt,deleteFlag,shiftFlag,dataCnt=11;
    int startI,endI,nSpecies,minI,maxI,indexI;
    double dPhi,theta;
	 int *sendN,*recvN,*start,*cntN;
	 double **sendData,**recvData;
    ptclList *p,*New,*prev,*tmp;
    int myrank, nTasks;

    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    startI=1;  endI=D->subSliceN+1;
    minI=D->minI;  maxI=D->maxI;
    nSpecies=D->nSpecies;
    dPhi=D->numSlice*2*M_PI;

    sendN=(int *)malloc(nTasks*sizeof(int ));
    recvN=(int *)malloc(nTasks*sizeof(int ));
    for(i=0; i<nTasks; i++) { sendN[i]=0; recvN[i]=0; }

	 if(D->mode!=Static) {
	    for(i=startI; i<endI; i++)   {
   	   for(s=0; s<nSpecies; s++) {
      	  p=particle[i].head[s]->pt;
	        while(p)  {
   	       theta=p->theta;
      	    if(theta>=dPhi)  intZ=(int)(theta/dPhi);
         	 else if(theta<0) intZ=((int)(theta/dPhi))-1;
	          else             intZ=0;
 
   	       indexI=i-startI+minI+intZ;
      	    for(n=0; n<nTasks; n++) {
         	   if(D->minmax[n]<=indexI && indexI<D->minmax[n+1]) {
					  if(n!=myrank)  sendN[n]+=1; else ;
					  n=nTasks;
					} else ;
				 }
				 p=p->next;
	        }	//End of while(p)
   	   }		//End of for(s)
	    }		//End of for(i)
	 } else ;

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
	 

    recvData=(double **)malloc(nTasks*sizeof(double *));
    sendData=(double **)malloc(nTasks*sizeof(double *));
    for(n=0; n<nTasks; n++)   {
      sendData[n]=(double *)malloc(sendN[n]*dataCnt*sizeof(double ));
      recvData[n]=(double *)malloc(recvN[n]*dataCnt*sizeof(double ));
    }

//    for(n=0; n<nTasks; n++) printf("myrank=%d, sendN[%d]=%d, recvN[%d]=%d\n",myrank,n,sendN[n],n,recvN[n]);


    cntN=(int *)malloc(nTasks*sizeof(int ));
    for(n=0; n<nTasks; n++) cntN[n]=0; 
    for(i=startI; i<endI; i++)
    {
      for(s=0; s<nSpecies; s++) {
        cnt=1;
        p=particle[i].head[s]->pt;
        while(p)  {
          if(cnt==1)  prev=p; else ;
              
			 deleteFlag=0;
			 shiftFlag=0;
          theta=p->theta;
          if(theta>=dPhi)  {
            intZ=(int)(theta/dPhi);
            theta-=dPhi*intZ;
				shiftFlag=1;	
          }
          else if(theta<0) {              
            intZ=((int)(theta/dPhi))-1;
            theta+=-1.0*dPhi*intZ;
				shiftFlag=1;	
          } 
          else   intZ=0;

			 if(D->mode==Static) intZ=0; else ;

          indexI=i-startI+minI+intZ;
			 if(indexI<minI || indexI>=maxI) shiftFlag=0; else ;
			 if(myrank==0 && indexI<minI) deleteFlag=1; else ;
			 if(myrank==nTasks-1 && indexI>=maxI) deleteFlag=1; else ;

          for(n=0; n<nTasks; n++) {
            if(D->minmax[n]<=indexI && indexI<D->minmax[n+1]) {
				  if(myrank!=n) {
                sendData[n][cntN[n]*dataCnt+0]=theta;
                sendData[n][cntN[n]*dataCnt+1]=p->x;
                sendData[n][cntN[n]*dataCnt+2]=p->y;
                sendData[n][cntN[n]*dataCnt+3]=p->gamma;
                sendData[n][cntN[n]*dataCnt+4]=p->px;
                sendData[n][cntN[n]*dataCnt+5]=p->py;
                sendData[n][cntN[n]*dataCnt+6]=p->weight;
                sendData[n][cntN[n]*dataCnt+7]=p->index;
                sendData[n][cntN[n]*dataCnt+8]=p->core;
                sendData[n][cntN[n]*dataCnt+9]=s;
                sendData[n][cntN[n]*dataCnt+10]=indexI;

                cntN[n]+=1;
					 deleteFlag=1;
				  } else ;
				  n=nTasks;
				} else ; 
			 }		//End of for(n<nTasks)

			 if(deleteFlag==1 && shiftFlag==0) {
            if(cnt==1)  {
              particle[i].head[s]->pt = p->next;
              p->next = NULL;
		        free(p);
              p=particle[i].head[s]->pt;
              cnt=1;
   	      } else {
              prev->next = p->next;
				  p->next = NULL;
				  free(p);
	           p=prev->next;
				}
			 } else if(deleteFlag==0 && shiftFlag==1) {
            if(cnt==1)  {
              p->theta=theta;
              particle[i].head[s]->pt = p->next;
	           p->next = particle[i+intZ].head[s]->pt;
	           particle[i+intZ].head[s]->pt = p;
	           p=particle[i].head[s]->pt;
	           cnt=1;
	         } else {
	           prev->next = p->next;
	           p->theta=theta;
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
      }		//End of for(s)
	 }			//End of for(i)


    for(n=0; n<nTasks; n++)  {
      if(myrank==n)  {
        for(nn=0; nn<nTasks; nn++) {
          if(n!=nn)
             MPI_Send(sendData[nn],sendN[nn]*dataCnt,MPI_DOUBLE,nn,myrank,MPI_COMM_WORLD);
          else ;
        }
      }
      else  {
        MPI_Recv(recvData[n],recvN[n]*dataCnt,MPI_DOUBLE,n,n,MPI_COMM_WORLD,&status);
        for(ii=0; ii<recvN[n]; ii++)  {
          indexI=recvData[n][ii*dataCnt+10];
          i=indexI-minI+startI;
          s=recvData[n][ii*dataCnt+9];

          New = (ptclList *)malloc(sizeof(ptclList));
          New->next = particle[i].head[s]->pt;
          particle[i].head[s]->pt = New;
          New->theta=recvData[n][ii*dataCnt+0];
          New->x=recvData[n][ii*dataCnt+1];
          New->y=recvData[n][ii*dataCnt+2];
          New->gamma=recvData[n][ii*dataCnt+3];
          New->px=recvData[n][ii*dataCnt+4];
          New->py=recvData[n][ii*dataCnt+5];
          New->weight=recvData[n][ii*dataCnt+6];
          New->index=recvData[n][ii*dataCnt+7];
          New->core=recvData[n][ii*dataCnt+8];
        }
      }
		MPI_Barrier(MPI_COMM_WORLD);  
    }

/*
    for(i=startI; i<endI; i++)
    {
      for(s=0; s<nSpecies; s++) {
        p=particle[i].head[s]->pt;
        while(p)  {
			  theta=p->theta;
			  if(theta<0 || theta>=dPhi) printf("myrank=%d,i=%d, theta=%g\n",myrank,i,theta); 
			  p=p->next;
		  }		
		}
	 }
*/

    for(n=0; n<nTasks; n++) {
      free(sendData[n]);
      free(recvData[n]);
    }
    free(sendData);
    free(recvData);
	 free(sendN);
	 free(recvN);
	 free(cntN);


}

/*
void rearrangeChicaneParticle(Domain *D)
{
    Particle *particle;
    particle=D->particle;

    int i,s,intZ,cnt,deleteFlag=0;
    int startI,endI,nSpecies,targetI;
    double dPhi,theta;
    ptclList *p,*New,*prev,*tmp;

    startI=1;  endI=D->subSliceN+1;
    nSpecies=D->nSpecies;
    dPhi=D->numSlice*2*M_PI;

    for(i=startI; i<endI; i++)
    {
      for(s=0; s<nSpecies; s++) {
        cnt=1;
        p=particle[i].head[s]->pt;
        while(p)  {
          if(cnt==1)
            prev=p;
          deleteFlag=0;
              
          theta=p->theta;
          if(theta>=dPhi)  {
            intZ=(int)(theta/dPhi);
            theta-=dPhi*intZ;
            deleteFlag=1;
          }
          else if(theta<0) {              
            intZ=((int)(theta/dPhi))-1;
            theta+=-1.0*dPhi*intZ;
            deleteFlag=1;
          } 
          else   intZ=0;

          if(deleteFlag==1)  {
            if(cnt==1)  {
				  if(i+intZ>endI)   { targetI=endI; theta+=dPhi*(i+intZ-endI); }
				  else if(i+intZ<0) { targetI=0;    theta+=dPhi*(i+intZ); }
				  else              { targetI=i+intZ; }
              p->theta=theta;    
              particle[i].head[s]->pt = p->next;
              p->next = particle[targetI].head[s]->pt;
              particle[targetI].head[s]->pt = p;
              p=particle[i].head[s]->pt;
              cnt=1;
            } else {
				  if(i+intZ>endI)   { targetI=endI; theta+=dPhi*(i+intZ-endI); }
				  else if(i+intZ<0) { targetI=0;    theta+=dPhi*(i+intZ); }
				  else              { targetI=i+intZ; }
              prev->next = p->next;
              p->theta=theta;    
              p->next = particle[targetI].head[s]->pt;
              particle[targetI].head[s]->pt = p;
              p=prev->next;
            }
          }		//End of if(deleteFlag==1)
          else {
            prev=p;
            p=p->next;
            cnt++;
          }              
        }	//End of while(p)
      }		//End of for(s)
    }		//End of for(i)
}
*/

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

	minmax=D->minmax;
	sendN=(int *)malloc(nTasks*sizeof(int ));
	recvN=(int *)malloc(nTasks*sizeof(int ));
	for(i=0; i<nTasks; i++) { sendN[i]=0; recvN[i]=0; }

   max=D->sliceN-shiftN;
	if(max<=0) ;
	else {
     for(i=startI; i<endI; i++) {
       indexI=i-startI+D->minI+shiftN;
       for(n=myrank+1; n<nTasks; n++) {
	      if(D->minmax[n]<=indexI && indexI<D->minmax[n+1]) {
			  sendN[n]+=1;
			} else ;
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

     num=N*numHarmony*2;
	  recvData=(double **)malloc(nTasks*sizeof(double *));
     sendData=(double **)malloc(nTasks*sizeof(double *));
     for(n=0; n<nTasks; n++)   {
       sendData[n]=(double *)malloc(sendN[n]*(num+1)*sizeof(double ));
       recvData[n]=(double *)malloc(recvN[n]*(num+1)*sizeof(double ));
     }
     start=(int *)malloc(nTasks*sizeof(int ));
     for(n=0; n<nTasks; n++) start[n]=0;

     for(i=endI-1; i>=startI; i--) {
       indexI=i-startI+minI+shiftN;
		 if(indexI>=maxI) {
         for(n=myrank+1; n<nTasks; n++) {
           if(D->minmax[n]<=indexI && indexI<D->minmax[n+1]) {
	 	  	    sendData[n][start[n]]=indexI;
				 start[n]+=1;
             for(h=0; h<numHarmony; h++) {
               for(j=0; j<N; j++) {
					  val=D->U[h][i][j];
	 	  	        sendData[n][start[n]+2*j+0]=creal(val);
	 	  	        sendData[n][start[n]+2*j+1]=cimag(val);
					  D->U[h][i][j]=0.0+I*0.0;
					}
					start[n]+=N*2;
				 }
			  } else ;
         }
		 } else if(minI<=indexI && indexI<maxI) {
         for(h=0; h<numHarmony; h++)
           for(j=0; j<N; j++) {
             D->U[h][i+shiftN][j]=D->U[h][i][j];
             D->U[h][i][j]=0.0*I*0.0;
			  }			  				 
		 } else ;

	  }	//End of for(i=startI; i<endI; i++)


     for(n=0; n<nTasks; n++)  {
       if(myrank==n)  {
         for(nn=0; nn<nTasks; nn++) {
           if(n!=nn)
             MPI_Send(sendData[nn],sendN[nn]*(num+1),MPI_DOUBLE,nn,myrank,MPI_COMM_WORLD);
			  else ;
			}
       }
       else  {
         MPI_Recv(recvData[n],recvN[n]*(num+1),MPI_DOUBLE,n,n,MPI_COMM_WORLD,&status);
			start[n]=0;
         for(ii=0; ii<recvN[n]; ii++)  {
			  indexI=recvData[n][start[n]];
			  i=indexI-minI+startI;
//			  if(myrank==0) printf("indexI=%d, i=%d\n",indexI,i);
			  start[n]+=1;
           for(h=0; h<numHarmony; h++) {
             for(j=0; j<N; j++) {
	 	  	      realV=recvData[n][start[n]+2*j+0];
	 	  	      imagV=recvData[n][start[n]+2*j+1];
					D->U[h][i][j]=realV+I*imagV;
			    }
				 start[n]+=N*2;
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
	}

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

   //printf("in test, chiON=%d\n",Chi->chiON);   //lala
	
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
		 D->chi_SSON=Chi->selfSeedON;
		 D->chi_d=Chi->d;
       D->bragTh=Chi->bragTh;
		 D->extincL=Chi->extincL;
		 D->chi0=Chi->chi0;

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
  D->chi_delay=0.0;
  D->chi_SSON=OFF;
  D->chi_d=0;
  D->bragTh=0;
  D->extincL=0;
  D->chi0=0;
}
