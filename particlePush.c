#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_bessel.h>
#include <complex.h>

void transversePush_3D(Domain *D,int iteration);
void push_theta_gamma_3D(Domain *D,int iteration);
void push_theta_gamma_1D(Domain *D,int iteration);

void transversePush(Domain *D,int iteration)
{
  switch(D->dimension)  {
  case 1 :
	  ;
//    particlePush1D(D);
    break;

  case 2 :
//    particlePush2D(&D,iteration);
    break;
  case 3 :
    transversePush_3D(D,iteration);
    break;
    
  default :
    printf("In particlePush, what dimension(%d)?\n",D->dimension);
  }

}


void push_theta_gamma(Domain *D,int iteration)
{
  switch(D->dimension)  {
  case 1 :
    push_theta_gamma_1D(D,iteration);
    break;

  case 2 :
//    particlePush2D(&D,iteration);
    break;
  case 3 :
    push_theta_gamma_3D(D,iteration);
    break;
    
  default :
    printf("In particlePush, what dimension(%d)?\n",D->dimension);
  }

}

void transversePush_3D(Domain *D,int iteration)
{
    int i,j,s,sliceI,startI,endI;
    LoadList *LL;
    double dz,ku,ks,kx,ky,K0,g;
    double x,y,z,px,py,gamma,invGam,v0[2],v1[2],M[2][2];
    double qx,qy,sqrtQx,sqrtQy,cnt,avePx,avePy,aveGam,send[4],recv[4];
    ptclList *p;

    int myrank, nTasks;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    dz=D->dz;
    ku=D->ku;
    ks=D->ks;
    kx=D->kx;    ky=D->ky;
    K0=D->K0;
    g=D->g;
    startI=1;  endI=D->subSliceN+1;

    avePx=avePy=aveGam=cnt=0.0;
    for(sliceI=startI; sliceI<endI; sliceI++)
    {    
      for(s=0; s<D->nSpecies; s++)  {
        p=D->particle[sliceI].head[s]->pt;
        while(p) {
       	  gamma=p->gamma;  invGam=1.0/gamma;
	
          qx=(K0*K0*invGam*kx*kx+eCharge*g/eMass/velocityC)*invGam;
          qy=(K0*K0*invGam*ky*ky-eCharge*g/eMass/velocityC)*invGam;
          sqrtQx=sqrt(fabs(qx));
          sqrtQy=sqrt(fabs(qy));

	  // calculate for x direction
  	  if(qx>0) {
            M[0][0]=cos(sqrtQx*dz*0.5);
            M[0][1]=sin(sqrtQx*dz*0.5)/sqrtQx*invGam;
            M[1][0]=-gamma*sqrtQx*sin(sqrtQx*dz*0.5);
            M[1][1]=M[0][0];
  	  } else if (qx<0) {
            M[0][0]=cosh(sqrtQx*dz*0.5);
            M[0][1]=sinh(sqrtQx*dz*0.5)/sqrtQx*invGam;
            M[1][0]=gamma*sqrtQx*sinh(sqrtQx*dz*0.5);
            M[1][1]=M[0][0];
  	  } else {
            M[0][0]=1.0;
            M[0][1]=dz*0.5*invGam;
            M[1][0]=0.0;
            M[1][1]=1.0;
  	  }
          v0[0]=p->x; v0[1]=p->px;
  	  v1[0]=0.0;  v1[1]=0.0;
	  for(i=0; i<2; i++)
	    for(j=0; j<2; j++)
              v1[i]+=M[i][j]*v0[j];
          p->x=v1[0];
          p->px=v1[1];
		
          // calculate for y direction
  	  if(qy>0) {
            M[0][0]=cos(sqrtQy*dz*0.5);
            M[0][1]=sin(sqrtQy*dz*0.5)/sqrtQy*invGam;
            M[1][0]=-gamma*sqrtQy*sin(sqrtQy*dz*0.5);
            M[1][1]=M[0][0];
  	  } else if (qy<0) {
            M[0][0]=cosh(sqrtQy*dz*0.5);
            M[0][1]=sinh(sqrtQy*dz*0.5)/sqrtQy*invGam;
            M[1][0]=gamma*sqrtQy*sinh(sqrtQy*dz*0.5);
            M[1][1]=M[0][0];
  	  } else {
            M[0][0]=1.0;
            M[0][1]=dz*0.5*invGam;
            M[1][0]=0.0;
            M[1][1]=1.0;
  	  }
          v0[0]=p->y; v0[1]=p->py;
  	  v1[0]=0.0;  v1[1]=0.0;
	  for(i=0; i<2; i++)
	    for(j=0; j<2; j++)
              v1[i]+=M[i][j]*v0[j];
          p->y=v1[0];
          p->py=v1[1];

	  avePx+=p->px*p->weight;
	  avePy+=p->py*p->weight;
	  aveGam+=p->gamma*p->weight;
	  cnt+=p->weight;

	  p=p->next;
        }
      }		//End of for(s)
    }		//End of for(i)
/*
    send[0]=avePx;
    send[1]=avePy;
    send[2]=aveGam;
    send[3]=cnt;
    if(myrank==0) {
      for(i=1; i<nTasks; i++) {
        MPI_Recv(recv,4,MPI_DOUBLE,i,i, MPI_COMM_WORLD,&status);
	send[0]+=recv[0];	//avePx
	send[1]+=recv[1];	//avePy
	send[2]+=recv[2];	//aveGam
	send[3]+=recv[3];	//cnt
      }
    } else 
      MPI_Send(send,4,MPI_DOUBLE,0,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if(send[3]==0.0) { send[0]=0.0; send[1]=0.0; send[2]=D->gamma0; send[3]=1.0; }
    else             { send[0]/=send[3]; send[1]/=send[3]; send[2]/=send[3]; }
    MPI_Bcast(send,nTasks,MPI_DOUBLE,0,MPI_COMM_WORLD);
    D->avePx=send[0];
    D->avePy=send[1];
    D->aveGam=send[2];
    D->totalCnt=send[3];
    */
         
}

void push_theta_gamma_3D(Domain *D,int iteration)
{
    int i,j,N,ii,jj,s,h,H,numHarmony,nx,ny,order,ll,L,f,F;
    int startI,endI,minI,maxI,sliceI,indexJ,idx;
    LoadList *LL;
    double complex U[D->numHarmony],compVal,Em[D->SCLmode];
    double dz,dx,dy,ku,ks,kx,ky,K0,K,xi,e_mc2,r,dr,dBessel;
    double x,y,z,px,py,gamma,theta,invGam,minX,minY,wakeE,invBeta,phi;
    double coef,JJ,J1,J2,wx[2],wy[2],sumTh,sumGam,sumEzPart,w[2];
    double k1,k2,k3,k4,l1,l2,l3,l4,tmp,dPhi,amp,arg,chi;
    ptclList *p;

    dz=D->dz;    K0=D->K0;
    ku=D->ku;    ks=D->ks;
    numHarmony=D->numHarmony;
    kx=D->kx; ky=D->ky;
    dx=D->dx; dy=D->dy; dr=D->dr;
    nx=D->nx; ny=D->ny;
	 dBessel = D->dBessel;
    minX=D->minX;  minY=D->minY;
    dPhi=2*M_PI*D->numSlice;
    e_mc2 = eCharge/eMass/velocityC/velocityC;	 
    
    L = D->SCLmode; F = D->SCFmode;	 
    startI=1;       endI=D->subSliceN+1;
    minI=D->minI;   maxI=D->maxI;
    N=D->nx*D->ny;

    for(sliceI=startI; sliceI<endI; sliceI++)
    {
      if(D->wakeONOFF==ON) wakeE=D->wakeE[sliceI-startI+minI]/mc2*1e-6;
      else                wakeE=0.0;      

      for(s=0; s<D->nSpecies; s++)  {
        p=D->particle[sliceI].head[s]->pt;
        while(p) {
          x=p->x;    y=p->y;  r=sqrt(x*x+y*y);
          px=p->px;  py=p->py;
		    if(x==0) phi = 0;
          else     phi = atan2(y,x);


  	       i=(int)((x-minX)/dx);
	       j=(int)((y-minY)/dy);
	       if(i>=0 && i<nx-1 && j>=0 && j<ny-1)  {
			   indexJ = (int)(r/dr);
				wy[1]=(r/dr-indexJ); wy[0]=1.0-wy[1];
     			for(ll=0; ll<L; ll++) {
				  Em[ll]=0.0+I*0.0;
              for(f=0; f<F; f++) {
//					 amp=0.0; arg=0.0;
				    for(jj=0; jj<2; jj++) {
//				      amp+=cabs(D->Ez[sliceI][indexJ][ll][f])*wy[jj];
//				      arg+=carg(D->Ez[sliceI][indexJ][ll][f])*wy[jj];
				      Em[ll]+=D->Ez[sliceI][indexJ][ll][f]*wy[jj];
					 }
//				    Em[ll] += amp*cexp(I*(f*phi+arg));
				  }
				}


       	   wx[1]=(x-minX)/dx-i; wx[0]=1.0-wx[1];
	         wy[1]=(y-minY)/dy-j; wy[0]=1.0-wy[1];
	         for(h=0; h<numHarmony; h++)  {
				
				  U[h]=0.0+I*0.0;
	           for(ii=0; ii<2; ii++) 
                for(jj=0; jj<2; jj++)  
      			   U[h]+=D->U[h][sliceI][(j+jj)*nx+(i+ii)]*wx[ii]*wy[jj];
					
//      		  U[h]=D->U[h][sliceI][(j)*nx+(i)];
	         }

	         K=K0*(1.0+kx*kx*x*x+ky*ky*y*y);
            
				// Step 1
	         theta=p->theta;
	         gamma=p->gamma; invGam=1.0/gamma;
            invBeta = 1.0-(1.0 + K*K + px*px + py*py)*0.5*invGam*invGam;
            invBeta = 1.0/invBeta;			 
            sumTh=sumGam=0.0;
	         xi=ks/ku*0.25*K*K*invGam*invGam;
            for(h=0; h<numHarmony; h++)  {
				  H = D->harmony[h];
              if(H%2==1)  {  //odd harmony
                coef=pow(-1.0,(H-1)*0.5);
                idx=(int)(H*xi/dBessel);
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-1)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+1)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=coef*(J1-J2);
              } else {
                chi = H*sqrt(2.0)*K*ks/ku*px*invGam*invGam;
                coef=pow(-1.0,(H-2)*0.5);
                idx=(int)(H*xi/dBessel);
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-2)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+2)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=chi*0.5*coef*(J1-J2);
              }
              compVal=U[h]*cexp(I*H*theta);
              sumTh -=2*JJ*K*cimag(compVal);
              sumGam-=2*JJ*K*creal(compVal);
            }
            sumEzPart = 0.0;
            for(ll=0; ll<L; ll++)  {
              tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
              sumEzPart += 2.0*tmp;
            }
            k1=(ku-ks*(1.0+K*K+px*px+py*py+sumTh)*0.5*invGam*invGam)*invBeta;
	         l1=ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart - wakeE;

				// Step 2 2 
	         theta=p->theta+0.5*dz*k1;
	         gamma=p->gamma+0.5*dz*l1; invGam=1.0/gamma;
            invBeta = 1.0-(1.0 + K*K + px*px + py*py)*0.5*invGam*invGam;
            invBeta = 1.0/invBeta;			 
	         xi=ks/ku*0.25*K*K*invGam*invGam;
            sumTh=sumGam=0.0;
            for(h=0; h<numHarmony; h++)  {
	           H = D->harmony[h];
              if(H%2==1)  {  //odd harmony
                coef=pow(-1.0,(H-1)*0.5);
                idx=(int)(H*xi/dBessel);
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-1)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+1)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=coef*(J1-J2);
              } else {
                chi = H*sqrt(2.0)*K*ks/ku*px*invGam*invGam;
                coef=pow(-1.0,(H-2)*0.5);
                idx=(int)(H*xi/dBessel);					 
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-2)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+2)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=chi*0.5*coef*(J1-J2);
              }
              compVal=U[h]*cexp(I*H*theta);
              sumTh -=2*JJ*K*cimag(compVal);
              sumGam-=2*JJ*K*creal(compVal);
            }
            sumEzPart = 0.0;
            for(ll=0; ll<L; ll++)  {
              tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
              sumEzPart += 2.0*tmp;
            }
            k2=(ku-ks*(1.0+K*K+px*px+py*py+sumTh)*0.5*invGam*invGam)*invBeta;
	         l2=ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart - wakeE;

				// Step 3 
	         theta=p->theta+0.5*dz*k2;
	         gamma=p->gamma+0.5*dz*l2; invGam=1.0/gamma;
            invBeta = 1.0-(1.0 + K*K + px*px + py*py)*0.5*invGam*invGam;
            invBeta = 1.0/invBeta;			 
	         xi=ks/ku*0.25*K*K*invGam*invGam;
            sumTh=sumGam=0.0;
            for(h=0; h<numHarmony; h++)  {
	           H = D->harmony[h];
              if(H%2==1)  {  //odd harmony
                coef=pow(-1.0,(H-1)*0.5);
                idx=(int)(H*xi/dBessel);
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-1)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+1)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=coef*(J1-J2);
              } else {
                chi = H*sqrt(2.0)*K*ks/ku*px*invGam*invGam;
                coef=pow(-1.0,(H-2)*0.5);
                idx=(int)(H*xi/dBessel);
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-2)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+2)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=chi*0.5*coef*(J1-J2);
              }
              compVal=U[h]*cexp(I*H*theta);
              sumTh -=2*JJ*K*cimag(compVal);
              sumGam-=2*JJ*K*creal(compVal);
            }
            sumEzPart = 0.0;
            for(ll=0; ll<L; ll++)  {
              tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
              sumEzPart += 2.0*tmp;
            }
            k3=(ku-ks*(1.0+K*K+px*px+py*py+sumTh)*0.5*invGam*invGam)*invBeta;
	         l3=ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart - wakeE;

				// Step 4 
	         theta=p->theta+dz*k3;
	         gamma=p->gamma+dz*l3; invGam=1.0/gamma;
            invBeta = 1.0-(1.0 + K*K + px*px + py*py)*0.5*invGam*invGam;
            invBeta = 1.0/invBeta;			 
	         xi=ks/ku*0.25*K*K*invGam*invGam;
            sumTh=sumGam=0.0;
            for(h=0; h<numHarmony; h++)  {
	           H = D->harmony[h];
              if(H%2==1)  {  //odd harmony
                coef=pow(-1.0,(H-1)*0.5);
                idx=(int)(H*xi/dBessel);
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-1)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+1)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=coef*(J1-J2);
              } else {
                chi = H*sqrt(2.0)*K*ks/ku*px*invGam*invGam;
                coef=pow(-1.0,(H-2)*0.5);
                idx=(int)(H*xi/dBessel);
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-2)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+2)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=chi*0.5*coef*(J1-J2);
              }
              compVal=U[h]*cexp(I*H*theta);
              sumTh -=2*JJ*K*cimag(compVal);
              sumGam-=2*JJ*K*creal(compVal);
            }
            sumEzPart = 0.0;
            for(ll=0; ll<L; ll++)  {
              tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
              sumEzPart += 2.0*tmp;
            }
            k4=(ku-ks*(1.0+K*K+px*px+py*py+sumTh)*0.5*invGam*invGam)*invBeta;
	         l4=ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart - wakeE;

            tmp=dz*1.0/6.0*(k1+2*k2+2*k3+k4);
	         if(tmp>=dPhi || tmp<=-dPhi) {
              printf("iteration=%d, dTheta=%g, r=%g,sumEzPart=%g\n",iteration,tmp,r,sumEzPart); 
              exit(0);
	         } else;
            p->theta+=tmp;
            p->gamma+=dz*1.0/6.0*(l1+2*l2+2*l3+l4);

	       } else ;

	       p=p->next;
        }
      }		//End of for(s)
    }		//End of for(sliceI)

}



void push_theta_gamma_1D(Domain *D,int iteration)
{
    int i,s,h,H,numHarmony,order,startI,endI,minI,maxI,ll,L,idx,intThe,bn;
    LoadList *LL;
    double complex U[D->numHarmony],Em[D->SCLmode],compVal;
    double dz,ku,ks,K0,K,xi,e_mc2,dBessel,w[2];
    double z,gamma,theta,invGam,invBeta,tmp,dPhi,sumGam,sumTh,sumEzPart,prevThe;
    double coef,JJ,J1,J2,wakeE,absU,absU2;
    double k1,k2,k3,k4,l1,l2,l3,l4;
	 double beta1,beta2,beta3,beta4;
	 double sumGam1,sumGam2,sumGam3,sumGam4;
    ptclList *p;

    L = D->SCLmode;
    startI=1;       endI=D->subSliceN+1;
    minI=D->minI;   maxI=D->maxI;
	 dz = D->dz;     K0=D->K0;
    ku=D->ku;       ks=D->ks;
    numHarmony = D->numHarmony;
    dBessel = D->dBessel;
    dPhi=2*M_PI*D->numSlice;
    e_mc2 = eCharge/eMass/velocityC/velocityC;	 
	 bn=D->bn;

    for(i=startI; i<endI; i++)
    {
      if(D->wakeONOFF==ON) wakeE=D->wakeE[i-startI+minI]/mc2*1e-6;
      else                wakeE=0.0;      

      for(s=0; s<D->nSpecies; s++)  {
        p=D->particle[i].head[s]->pt;
        while(p) {

          for(ll=0; ll<L; ll++) 
			   Em[ll]=D->Ez[i][0][ll][0];

          for(h=0; h<numHarmony; h++)  
		      U[h]=D->U[h][i][0];
      
          K=K0;

  //--------------------------- theta --------------------------------
          prevThe=p->theta;
          theta=p->theta;
          gamma=p->gamma; invGam=1.0/gamma;
          invBeta = 1.0-(1.0 + K*K)*0.5*invGam*invGam;
          invBeta = 1.0/invBeta;
          sumTh=sumGam=0.0;
          xi=ks/ku*0.25*K*K*invGam*invGam;
          for(h=0; h<numHarmony; h++)  {
			   H = D->harmony[h];
			   if(H%2==1) {
              coef=pow(-1.0,(H-1)*0.5);
              idx=(int)(H*xi/dBessel);
				  if(idx>bn) printf("1 : idx=%d\n",idx);
				  if(idx>bn-1) idx=bn-2;
              w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
              order=(H-1)*0.5;
              J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              order=(H+1)*0.5;
              J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              JJ=coef*(J1-J2);
				}
            compVal=U[h]*cexp(I*H*theta);
            sumTh -=2*JJ*K*cimag(compVal);
          }
          sumEzPart = 0.0;
          for(ll=0; ll<L; ll++)  {
            tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
            sumEzPart += 2.0*tmp;
          }
          k1=(ku-ks*(1.0+K*K+sumTh)*0.5*invGam*invGam)*invBeta;
			 
	       theta=p->theta+0.5*dz*k1;
          invBeta = 1.0-(1.0 + K*K)*0.5*invGam*invGam;
          invBeta = 1.0/invBeta;
          sumTh=sumGam=0.0;
          xi=ks/ku*0.25*K*K*invGam*invGam;
          for(h=0; h<numHarmony; h++)  {
			   H = D->harmony[h];
			   if(H%2==1) {
              coef=pow(-1.0,(H-1)*0.5);
              idx=(int)(H*xi/dBessel);
				  if(idx>bn) printf("2 : idx=%d\n",idx);
				  if(idx>bn-1) idx=bn-2;
              w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
              order=(H-1)*0.5;
              J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              order=(H+1)*0.5;
              J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              JJ=coef*(J1-J2);
				}
            compVal=U[h]*cexp(I*H*theta);
            sumTh -=2*JJ*K*cimag(compVal);
          }
          sumEzPart = 0.0;
          for(ll=0; ll<L; ll++)  {
            tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
            sumEzPart += 2.0*tmp;
          }
          k2=(ku-ks*(1.0+K*K+sumTh)*0.5*invGam*invGam)*invBeta;

	       theta=p->theta+0.5*dz*k2;
          invBeta = 1.0-(1.0 + K*K)*0.5*invGam*invGam;
          invBeta = 1.0/invBeta;
          sumTh=sumGam=0.0;
          xi=ks/ku*0.25*K*K*invGam*invGam;
          for(h=0; h<numHarmony; h++)  {
			   H = D->harmony[h];
			   if(H%2==1) {
              coef=pow(-1.0,(H-1)*0.5);
              idx=(int)(H*xi/dBessel);
				  if(idx>bn) printf("3 : idx=%d\n",idx);
				  if(idx>bn-1) idx=bn-2;
              w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
              order=(H-1)*0.5;
              J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              order=(H+1)*0.5;
              J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              JJ=coef*(J1-J2);
				}
            compVal=U[h]*cexp(I*H*theta);
            sumTh -=2*JJ*K*cimag(compVal);
          }
          sumEzPart = 0.0;
          for(ll=0; ll<L; ll++)  {
            tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
            sumEzPart += 2.0*tmp;
          }
          k3=(ku-ks*(1.0+K*K+sumTh)*0.5*invGam*invGam)*invBeta;

	       theta=p->theta+dz*k3;
          invBeta = 1.0-(1.0 + K*K)*0.5*invGam*invGam;
          invBeta = 1.0/invBeta;
          sumTh=sumGam=0.0;
          xi=ks/ku*0.25*K*K*invGam*invGam;
          for(h=0; h<numHarmony; h++)  {
			   H = D->harmony[h];
			   if(H%2==1) {
              coef=pow(-1.0,(H-1)*0.5);
              idx=(int)(H*xi/dBessel);
				  if(idx>bn) printf("4 : idx=%d\n",idx);
				  if(idx>bn-1) idx=bn-2;
              w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
              order=(H-1)*0.5;
              J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              order=(H+1)*0.5;
              J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              JJ=coef*(J1-J2);
				}
            compVal=U[h]*cexp(I*H*theta);
            sumTh -=2*JJ*K*cimag(compVal);
          }
          sumEzPart = 0.0;
          for(ll=0; ll<L; ll++)  {
            tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
            sumEzPart += 2.0*tmp;
          }
          k4=(ku-ks*(1.0+K*K+sumTh)*0.5*invGam*invGam)*invBeta;

          tmp=dz*1.0/6.0*(k1+2*k2+2*k3+k4);
	       if(tmp>=dPhi) {
//              printf("iteration=%d, dTheta=%g, sumEzPart=%g,k1*dz=%g, k2*dz=%g, k3*dz=%g, k4*dz=%g\n",iteration,tmp,sumEzPart,k1*dz,k2*dz,k3*dz,k4*dz);
				  tmp-=dPhi;
	       } else if(tmp<=-1*dPhi) {
//              printf("iteration=%d, dTheta=%g, sumEzPart=%g,k1*dz=%g, k2*dz=%g, k3*dz=%g, k4*dz=%g\n",iteration,tmp,sumEzPart,k1*dz,k2*dz,k3*dz,k4*dz);
				  tmp+=dPhi;
			 }
          p->theta+=tmp;


  //--------------------------- gamma --------------------------------
           // Step 1
          theta=0.5*(p->theta+prevThe);
          sumEzPart = 0.0;
          for(ll=0; ll<L; ll++)  {
            tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
            sumEzPart += 2.0*tmp;
          }

          gamma=p->gamma; invGam=1.0/gamma;
          invBeta = 1.0-(1.0 + K*K)*0.5*invGam*invGam;
          invBeta = 1.0/invBeta;
          sumTh=sumGam=0.0;
          xi=ks/ku*0.25*K*K*invGam*invGam;
          for(h=0; h<numHarmony; h++)  {
			   H = D->harmony[h];
			   if(H%2==1) {
              coef=pow(-1.0,(H-1)*0.5);
              idx=(int)(H*xi/dBessel);
				  if(idx>bn) printf("gamma 1 : idx=%d\n",idx);
				  if(idx>bn-1) idx=bn-2;
              w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
              order=(H-1)*0.5;
              J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              order=(H+1)*0.5;
              J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              JJ=coef*(J1-J2);
				}
            compVal=U[h]*cexp(I*H*theta);
            sumGam-=2*JJ*K*creal(compVal);
          }
          l1=ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart - wakeE;
			 
	       gamma=p->gamma+0.5*dz*l1; invGam=1.0/gamma;
          invBeta = 1.0-(1.0 + K*K)*0.5*invGam*invGam;
          invBeta = 1.0/invBeta;
          sumTh=sumGam=0.0;
          xi=ks/ku*0.25*K*K*invGam*invGam;
          for(h=0; h<numHarmony; h++)  {
			   H = D->harmony[h];
			   if(H%2==1) {
              coef=pow(-1.0,(H-1)*0.5);
              idx=(int)(H*xi/dBessel);
				  if(idx>bn) printf("gamma 2 : idx=%d\n",idx);
				  if(idx>bn-1) idx=bn-2;
              w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
              order=(H-1)*0.5;
              J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              order=(H+1)*0.5;
              J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              JJ=coef*(J1-J2);
				}
            compVal=U[h]*cexp(I*H*theta);
            sumGam-=2*JJ*K*creal(compVal);
          }
          l2=ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart - wakeE;

	       gamma=p->gamma+0.5*dz*l2; invGam=1.0/gamma;
          invBeta = 1.0-(1.0 + K*K)*0.5*invGam*invGam;
          invBeta = 1.0/invBeta;
          sumTh=sumGam=0.0;
          xi=ks/ku*0.25*K*K*invGam*invGam;
          for(h=0; h<numHarmony; h++)  {
			   H = D->harmony[h];
			   if(H%2==1) {
              coef=pow(-1.0,(H-1)*0.5);
              idx=(int)(H*xi/dBessel);
				  if(idx>bn) printf("gamma 3 : idx=%d\n",idx);
				  if(idx>bn-1) idx=bn-2;
              w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
              order=(H-1)*0.5;
              J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              order=(H+1)*0.5;
              J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              JJ=coef*(J1-J2);
				}
            compVal=U[h]*cexp(I*H*theta);
            sumGam-=2*JJ*K*creal(compVal);
          }
          l3=ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart - wakeE;

	       gamma=p->gamma+dz*l3; invGam=1.0/gamma;
          invBeta = 1.0-(1.0 + K*K)*0.5*invGam*invGam;
          invBeta = 1.0/invBeta;
          sumTh=sumGam=0.0;
          xi=ks/ku*0.25*K*K*invGam*invGam;
          for(h=0; h<numHarmony; h++)  {
			   H = D->harmony[h];
			   if(H%2==1) {
              coef=pow(-1.0,(H-1)*0.5);
              idx=(int)(H*xi/dBessel);
				  if(idx>bn) printf("gamma 4 : idx=%d\n",idx);
				  if(idx>bn-1) idx=bn-2;
              w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
              order=(H-1)*0.5;
              J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              order=(H+1)*0.5;
              J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
              JJ=coef*(J1-J2);
				}
            compVal=U[h]*cexp(I*H*theta);
            sumGam-=2*JJ*K*creal(compVal);
          }
          l4=ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart - wakeE;

          p->gamma+=dz*1.0/6.0*(l1+2*l2+2*l3+l4);

	       p=p->next;
        }
      }		//End of for(s)
    }		//End of for(i)

}


//lala
/*	
void push_theta_gamma_1D(Domain *D,int iteration)
{
    int i,s,h,harmony,order,startI,endI,minI,maxI;
    LoadList *LL;
    double complex U[D->harmony],Ez[D->harmony],compValU,compValEz,expITh,expIhTh;
    double dz,ku,ks,K0,K,xi,ar[D->harmony],psi[D->harmony];
    double z,gamma,theta,invGam,invBeta,tmp,w,dPhi;
    double coef,JJ[D->harmony],J1,J2,sumTh,sumGam,wakeE;
    double k1,k2,k3,k4,l1,l2,l3,l4,totalEz;
    ptclList *p;

    startI=1;       endI=D->subSliceN+1;
    minI=D->minI;   maxI=D->maxI;
    dz=D->dz;
    ku=D->ku;
    ks=D->ks;
    K0=D->K0;
    harmony=D->harmony;
    dPhi=2*M_PI*D->numSlice;

    for(i=startI; i<endI; i++)
    {
      for(h=0; h<harmony; h++) { 
        ar[h]=cabs(D->U[h][i][0]);
        psi[h]=carg(I*D->U[h][i][0]);
        Ez[h]=D->Ez[h][i][0];
      }
      if(D->wakeONOFF=ON) wakeE=D->wakeE[i-startI+minI]/mc2*1e-6;
      else                wakeE=0.0;      

      for(s=0; s<D->nSpecies; s++)  {
        p=D->particle[i].head[s]->pt;
        while(p) {
          K=K0;
	  xi=K*K*0.5/(1+K*K);

          for(h=1; h<=harmony; h++)  {
            if(h%2==1)  {  //odd harmony
              coef=pow(-1.0,(h-1)/2);
              order=(h-1)/2;
              J1=gsl_sf_bessel_Jn(order,h*xi);
              order=(h+1)/2;
              J2=gsl_sf_bessel_Jn(order,h*xi);
              JJ[h-1]=coef*(J1-J2);
            } else {
              JJ[h-1]=0.0;
            }
          }

	  theta=p->theta;
	  gamma=p->gamma; invGam=1.0/gamma;
//	  invBeta=1.0+(1+K*K)*0.5*invGam*invGam;
	  invBeta=1.0+(1+K*K+K*JJ[0]*2*ar[0]*cos(theta+psi[0]))*0.5*invGam*invGam;
          compValEz=0.0+I*0.0;
          for(h=1; h<=harmony; h++)  {
	    expIhTh=cexp(I*h*theta);
            compValEz+=Ez[h-1]*expIhTh;
	  }
	  totalEz=eCharge/eMass/velocityC/velocityC*2*creal(compValEz)-wakeE;	  
          k1=ku+ks*(1.0-invBeta);
	  l1=-ks*JJ[0]*K*ar[0]*invBeta*invGam*sin(theta+psi[0])+totalEz;

	  theta=p->theta+0.5*dz*k1;
	  gamma=p->gamma+0.5*dz*l1; invGam=1.0/gamma;
//	  invBeta=1.0+(1+K*K)*0.5*invGam*invGam;
	  invBeta=1.0+(1+K*K+K*JJ[0]*ar[0]*2*cos(theta+psi[0]))*0.5*invGam*invGam;
          compValEz=0.0+I*0.0;
          for(h=1; h<=harmony; h++)  {
	    expIhTh=cexp(I*h*theta);
            compValEz+=Ez[h-1]*expIhTh;
	  }
	  totalEz=eCharge/eMass/velocityC/velocityC*2*creal(compValEz)-wakeE;
          k2=ku+ks*(1.0-invBeta);
	  l2=-ks*JJ[0]*K*ar[0]*invBeta*invGam*sin(theta+psi[0])+totalEz;

	  theta=p->theta+0.5*dz*k2;
	  gamma=p->gamma+0.5*dz*l2; invGam=1.0/gamma;
//	  invBeta=1.0+(1+K*K)*0.5*invGam*invGam;
	  invBeta=1.0+(1+K*K+K*JJ[0]*ar[0]*2*cos(theta+psi[0]))*0.5*invGam*invGam;
          compValEz=0.0+I*0.0;
          for(h=1; h<=harmony; h++)  {
	    expIhTh=cexp(I*h*theta);
            compValEz+=Ez[h-1]*expIhTh;
	  }
	  totalEz=eCharge/eMass/velocityC/velocityC*2*creal(compValEz)-wakeE;
          k3=ku+ks*(1.0-invBeta);
	  l3=-ks*JJ[0]*K*ar[0]*invBeta*invGam*sin(theta+psi[0])+totalEz;

	  theta=p->theta+dz*k3;
	  gamma=p->gamma+dz*l3; invGam=1.0/gamma;
//	  invBeta=1.0+(1+K*K)*0.5*invGam*invGam;
	  invBeta=1.0+(1+K*K+K*JJ[0]*ar[0]*2*cos(theta+psi[0]))*0.5*invGam*invGam;
          compValEz=0.0+I*0.0;
          for(h=1; h<=harmony; h++)  {
	    expIhTh=cexp(I*h*theta);
            compValEz+=Ez[h-1]*expIhTh;
	  }
	  totalEz=eCharge/eMass/velocityC/velocityC*2*creal(compValEz)-wakeE;
          k4=ku+ks*(1.0-invBeta);
	  l4=-ks*JJ[0]*K*ar[0]*invBeta*invGam*sin(theta+psi[0])+totalEz;

          tmp=dz*1.0/6.0*(k1+2*k2+2*k3+k4);
	  if(tmp>dPhi || tmp<-1.0*dPhi) {
		  printf("tmp=%g,dPhi=%g\n",tmp,dPhi);
	  } else ;
          p->theta+=tmp;
          p->gamma+=dz*1.0/6.0*(l1+2*l2+2*l3+l4);

//          // energy loss by wake field
//          delGam=D->wakeE[i-startI+minI]/mc2;	  
//          p->gamma-=delGam*dz;

          tmp=p->theta/dPhi;
	  w=tmp-(int)tmp;
	  p->theta=w*dPhi;
          if(p->theta>dPhi)   p->theta=p->theta-dPhi;
          else if(p->theta<0) p->theta=dPhi+p->theta;
	  

	  p=p->next;
        }
      }		//End of for(s)
    }		//End of for(i)
}
*/

double Runge_Kutta_gamma(double complex *U,double theta,double harmony,double ks,double K,double xi,double gamma,double invBeta0,double dz)
{
  double sum,sinX,cosX,coef,J1,J2,JJ,k1,k2,k3,k4;
  double complex compVal;
  int h,order;

  sum=0.0;
  sinX=sin(theta); cosX=cos(theta);
  for(h=1; h<=harmony; h++)  {
    if(h%2==1)  {  //odd harmony
      coef=pow(-1.0,h-1);
      order=(h-1)/2;
      J1=gsl_sf_bessel_Jn(order,h*xi);
      order=(h+1)/2;
      J2=gsl_sf_bessel_Jn(order,h*xi);
      JJ=coef*(J1-J2);
      compVal=U[h-1]*cexp(-I*theta);
      sum=JJ*2*creal(compVal)*h;
    } else {
      JJ=0.0;
      compVal=U[h-1]*cexp(-I*theta);
      sum=JJ*2*creal(compVal)*h;
    }
  }

  k1=-1.0*K*ks*invBeta0*sum/(gamma);
  
  k2=-1.0*K*ks*invBeta0*sum/(gamma+0.5*k1);
  
  k3=-1.0*K*ks*invBeta0*sum/(gamma+0.5*k2);
  
  k4=-1.0*K*ks*invBeta0*sum/(gamma+k3);

  return (k1/6.0+k2/3.0+k3/3.0+k4/6.0)*dz;
}

void phaseShift(Domain *D,int iteration)
{
	 int s,i,sliceI,startI,endI;
	 double shiftValue,theta;
    LoadList *LL;
    ptclList *p;
	 PhaseShifter *PS;
    int myrank, nTasks;
    MPI_Status status;
	 MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	 MPI_Comm_size(MPI_COMM_WORLD, &nTasks); 

    startI=1;       endI=D->subSliceN+1;

    PS=D->psList;
    while(PS->next) {
      for(i=0; i<PS->num; i++) {
        if(iteration==PS->step[i]) {
          shiftValue=PS->phase;
          if(myrank==0) printf("phase shift with %g is done at step%d.\n",shiftValue,iteration);  else ;

          for(sliceI=startI; sliceI<endI; sliceI++)  {
            for(s=0; s<D->nSpecies; s++)  {
              p=D->particle[sliceI].head[s]->pt;
              while(p) {
			       theta=p->theta;
                p->theta=theta+shiftValue;

	             p=p->next;
              }
            }		//End of for(s)
          }		//End of for(sliceI)

		  } else ;
	   }
		PS=PS->next;
	 }

}

