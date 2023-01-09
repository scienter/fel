#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

//void saveParticleHDF(Domain *D,int iteration,int s,double minPx,double density);
//void saveDensityHDF(Domain *D,int iteration);
//void saveEFieldHDF(Domain *D,int iteration);
//void saveBFieldHDF(Domain *D,int iteration);

void saveParticle(Domain *D,int iteration,int sliceI);
void saveField(Domain *D,int iteration);

void saveFile(Domain D,int iteration,int sliceI)
{
  int myrank, nTasks,s;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  LoadList *LL;

  //save field
  if(D.fieldSave==ON)   {
    if(D.saveFieldMode==TXT)  
      saveField(&D,iteration);
//    else if(D.saveFieldMode==HDF)   {
//      saveEFieldHDF(&D,iteration);
//      saveBFieldHDF(&D,iteration);
//    }  else	;
  }	else	;

/*
  //save density
  if(D.densitySave==ON)
  {
    if(D.saveDensityMode==TXT)
      saveDensity(&D,iteration);
    else if(D.saveDensityMode==HDF)
      saveDensityHDF(&D,iteration);
    if(myrank==0)
      printf("density%d is made.\n",iteration); 
  }  else	;
 
  //save current
  if(D.currentSave==ON)
  {
    if(D.saveCurrentMode==TXT)
      saveCurrent(&D,iteration);
//    else if(D.saveCurrentMode==HDF)
//      saveCurrentHDF(&D,iteration);
    if(myrank==0)
      printf("current%d is made.\n",iteration); 
  }  else	;
 */

}

void saveField(Domain *D,int iteration)
{
    int i,h,harmony,sliceN;
    double bucketX,x,UR,UI,P,coef,coef2,minX,x0,area;
    char name[100];
    FILE *out;
    int myrank, nprocs;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    minX=D->minX;
    x0=D->dx*iteration+D->minX+D->lambda0*0.5;
    bucketX=D->numSlice*D->lambda0;
    sliceN=D->sliceN;
    harmony=D->harmony;

    sprintf(name,"field%d_%d",iteration,myrank);
    out = fopen(name,"w");

    area=2*M_PI*D->spotSigR*D->spotSigR;
    coef=eMass*velocityC*velocityC*D->ks/eCharge;
    coef2=coef*coef*eps0*velocityC*area;
//    coef=eMass*D->ks*velocityC*velocityC*sqrt(2.0)/eCharge;

    switch(D->dimension) {
    case 1:
      for(i=0; i<sliceN; i++) {
        x=i*bucketX+x0;
        fprintf(out,"%.10g ",x);
        for(h=0; h<harmony; h++) {
          UR=D->UR[iteration][i][h];
          UI=D->UI[iteration][i][h];
	  P=(UR*UR+UI*UI)*coef2;
          fprintf(out,"%g %g %g",UR,UI,P);
//          fprintf(out,"%g %g ",E);
        }
        fprintf(out,"\n");
      }
      break;
    default :
      printf("what field_type? and what dimension?\n");
    }

    printf("%s is made.\n",name);
    fclose(out);
}

void savePower(Domain *D)
{
    int n,i,h,harmony;
    double x,P,coef1,coef2,coef,area,UR,UI,sum;
    char name[100];
    FILE *out;
    int myrank, nprocs;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    harmony=D->harmony;
    area=2*M_PI*D->spotSigR*D->spotSigR;
    coef1=eMass*velocityC*velocityC*D->ks/eCharge;
    coef2=coef1*coef1*eps0*velocityC*area;

    if(D->mode==Static) coef=coef2/(D->sliceN*1.0);
    else                coef=coef2*(D->lambda0/velocityC*D->numSlice);


    sprintf(name,"power");
    out = fopen(name,"w");

//    area=2*M_PI*D->spotSigR*D->spotSigR;
//    coef=eMass*velocityC*velocityC*D->ks/eCharge;
//    coef2=coef*coef*0.5*eps0*velocityC*area;

    switch(D->dimension) {
    case 1:
      for(n=0; n<=D->maxStep; n++)
      {
        x=n*D->dx;
        fprintf(out,"%g ",x);
        for(h=0; h<harmony; h++)  {
          sum=0.0;
          for(i=0; i<D->sliceN; i++)  {    
            UR=D->UR[n][i][h];
            UI=D->UI[n][i][h];
            sum+=UR*UR+UI*UI;
          }
          P=sum*coef;
          fprintf(out,"%g ",P);
        }
        fprintf(out,"\n");
      }
      break;
    default :
      printf("In savePower, what dimension?\n");
    }

    printf("%s is made.\n",name);
    fclose(out);
}

void saveParticle(Domain *D,int iteration,int sliceI)
{
  int i,s,index,totalCnt,subCnt;
  double x,theta,gamma,weight,dPhi,x0,bucketX,minX,tmp,w;
  char name[100];
  LoadList *LL;
  FILE *out;
  int myrank, nprocs;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  minX=D->minX;
  dPhi=D->numSlice*2*M_PI;
  x0=iteration*D->dx;
  bucketX=D->numSlice*D->lambda0;

  LL=D->loadList;
  s=0;
  while(LL->next) {
    sprintf(name,"%dParticle%d_%d",s,iteration,myrank);
    out = fopen(name,"a");    

    totalCnt=LL->totalCnt;
    subCnt=LL->subCnt;

    switch (D->dimension)  {
    //1D
    case 1:
      for(i=0; i<subCnt; i++) {
        theta=D->particle[s][i].theta;
        tmp=theta/dPhi;
        w=tmp-(int)tmp;
        theta=w*dPhi;
        if(theta>dPhi)   theta=theta-dPhi;
        else if(theta<0) theta=dPhi+theta;

        x=x0+(sliceI+theta/dPhi)*bucketX+minX; 
        gamma=D->particle[s][i].gamma;
        weight=D->particle[s][i].weight;
        index=D->particle[s][i].index;
        fprintf(out,"%.10g %g %g %d\n",x,gamma,weight,index);
      }	//End of while(p)
      break;
    default:
      ;
    }

//  printf("%s is made.\n",name);
    fclose(out);


    LL=LL->next;
    s++;
  }
}

/*
void saveDensity(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,s,ii,jj,kk,i1,j1,k1;
    int l,m,n;
    double x,y,z,Wx[3],Wy[3],Wz[3],weight,charge;
    char name[100];
    Particle ***particle;
    particle=D->particle;
    ptclList *p;
    LoadList *LL;
    FILE *out;
    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    double rho0[D->nSpecies];
    s=0;
    LL=D->loadList;
    while(LL->next)
    {
       rho0[s]=LL->density;
       LL=LL->next;
       s++;
    }

    switch (D->dimension)  {
    //1D
    case 1 :
      j=k=0;
      //initializing density
      for(s=0; s<D->nSpecies; s++)
      {
        for(i=0; i<iend+3; i++)   D->Rho[i][j][k]=0.0;
      
        for(i=istart; i<iend; i++)
        {
          p=particle[i][j][k].head[s]->pt;
          while(p)            {
            weight=p->weight; charge=p->charge;
            x=p->x;
            i1=(int)(i+x+0.5);
            x=i+x-i1;
            Wx[0]=0.5*(0.5-x)*(0.5-x);
            Wx[1]=0.75-x*x;
            Wx[2]=0.5*(x+0.5)*(x+0.5);
            for(ii=0; ii<3; ii++)  {
              l=i1-1+ii;
              if(istart<=l && l<iend)
                D->Rho[l][j][k]+=Wx[ii]*rho0[s]*weight*charge;
              else	;
            }
            p=p->next;
          }
        }

        sprintf(name,"%ddensity%d_%d",s,iteration,myrank);
        out = fopen(name,"w");    
        for(i=istart-1; i<=iend; i++)     {
          x=(i-istart+D->minXSub)*D->dx*D->lambda;
          fprintf(out,"%g %g\n",x,D->Rho[i][j][k]);    
        }
        fclose(out);
      }		//End for(S)
      break;

    //2D
    case 2 :
      k=0;
      //initializing density
      for(i=0; i<=iend; i++)
        for(j=jstart-1; j<=jend; j++)
          D->Rho[i][j][k]=0.0;

      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(s=0; s<D->nSpecies; s++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              weight=p->weight;
              x=p->x; y=p->y; z=p->z;
              i1=(int)(i+x+0.5);
              j1=(int)(j+y+0.5);
              x=i+x-i1;
              y=j+y-j1;
              Wx[0]=0.5*(0.5-x)*(0.5-x);
              Wx[1]=0.75-x*x;
              Wx[2]=0.5*(x+0.5)*(x+0.5);
              Wy[0]=0.5*(0.5-y)*(0.5-y);
              Wy[1]=0.75-y*y;
              Wy[2]=0.5*(y+0.5)*(y+0.5);

              for(jj=0; jj<3; jj++)
                for(ii=0; ii<3; ii++)
                {
                  l=i1-1+ii;
                  m=j1-1+jj;
                  if(istart<=l && l<iend && jstart<=m && m<jend)
                    D->Rho[l][m][k]+=Wx[ii]*Wy[jj]*rho0[s]*weight;
                }
              p=p->next;
            }
          }

      sprintf(name,"density%d_%d",iteration,myrank);
      out = fopen(name,"w");    
      for(i=istart-1; i<=iend; i++)
      {
        for(j=jstart-1; j<=jend; j++)
        {
          x=(i-istart+D->minXSub)*D->dx*D->lambda;
          y=(j-jstart+D->minYSub)*D->dy*D->lambda;
          fprintf(out,"%g %g %g\n",x,y,D->Rho[i][j][k]);    
        }           
        fprintf(out,"\n");    
      }
      fclose(out);

      break;

    case 3 :
      //initializing density
      for(i=0; i<=iend; i++)
        for(j=jstart-1; j<=jend; j++)
          for(k=kstart-1; k<=kend; k++)
            D->Rho[i][j][k]=0.0;

      for(i=istart-1; i<iend-1; i++)
        for(j=jstart-1; j<jend-1; j++)
          for(k=kstart-1; k<=kend; k++)
            for(s=0; s<D->nSpecies; s++)
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                x=p->x; y=p->y; z=p->z;
                i1=(int)(i+x+0.5);
                j1=(int)(j+y+0.5);
                k1=(int)(k+z+0.5);
                x=i+x-i1;
                y=j+y-j1;
                z=k+z-k1;
                Wx[0]=0.5*(0.5-x)*(0.5-x);
                Wx[1]=0.75-x*x;
                Wx[2]=0.5*(x+0.5)*(x+0.5);
                Wy[0]=0.5*(0.5-y)*(0.5-y);
                Wy[1]=0.75-y*y;
                Wy[2]=0.5*(y+0.5)*(y+0.5);
                Wz[0]=0.5*(0.5-z)*(0.5-z);
                Wz[1]=0.75-z*z;
                Wz[2]=0.5*(z+0.5)*(z+0.5);

                for(ii=0; ii<3; ii++)
                  for(jj=0; jj<3; jj++)
                    for(kk=0; kk<3; kk++)
                      D->Rho[i1-1+ii][j1-1+jj][k1-1+kk]
                            +=Wx[ii]*Wy[jj]*Wz[kk]*rho0[s];
                p=p->next;
              }
            }

      sprintf(name,"density%d_%d",iteration,myrank);
      out = fopen(name,"w");    
      for(i=istart-1; i<=iend; i++)
      {
        for(j=jstart-1; j<=jend; j++)
        {
          for(k=kstart-1; k<=kend; k++)
          {
            x=(i-istart+D->minXSub)*D->dx*D->lambda;
            y=(j-jstart+D->minYSub)*D->dy*D->lambda;
            z=(k-kstart+D->minZSub)*D->dz*D->lambda;
            fprintf(out,"%g %g %g %g\n",x,y,z,D->Rho[i][j][k]);    
          }           
          fprintf(out,"\n");    
        }           
        fprintf(out,"\n");    
      }
      fclose(out);
      break;

    }

}
*/

/*
void saveCurrent(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;
    char name[100];
    double x,y,z,Jx,Jy,Jz;
    FILE *out;
    int myrank, nprocs;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    sprintf(name,"current%d_%d",iteration,myrank);
    out = fopen(name,"w");

    switch(D->dimension) {

    case 1:
      j=k=0;
      for(i=istart; i<iend; i++)      {
        x=(i-2+D->minXSub)*D->dx*D->lambda;
        Jx=D->Jx[i][j][k]; Jy=D->Jy[i][j][k]; Jz=D->Jz[i][j][k];
        fprintf(out,"%g %g %g %g\n",x,Jx,Jy,Jz);
      }
      fclose(out);
      break;

    case 2:
      k=0;
      for(i=istart; i<iend; i++)    {
        for(j=jstart; j<jend; j++)        {
          x=(i-2+D->minXSub)*D->dx*D->lambda;
          y=(j-2+D->minYSub)*D->dy*D->lambda;
          Jx=D->Jx[i][j][k]; Jy=D->Jy[i][j][k]; Jz=D->Jz[i][j][k];
          fprintf(out,"%g %g %g %g %g\n",x,y,Jx,Jy,Jz);
        }
        fprintf(out,"\n");                 
      }
      fclose(out);
      break;
    default :
      printf("In saveCurrent, what dimension?\n");
    }
}
*/

