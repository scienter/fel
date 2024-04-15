#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void rearrangeParticles(Domain *D,int iteration)
{
    Particle *particle;
    particle=D->particle;

    int i,s,intZ,cnt,deleteFlag=0;
    int startI,endI,nSpecies;
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
            intZ=1;
	    theta-=dPhi;
            deleteFlag=1;
          }
          else if(theta<0) {              
            intZ=-1;
	    theta+=dPhi;
            deleteFlag=1;
          } 
          else   intZ=0;

          if(deleteFlag==1)  {
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

void periodicParticles(Domain *D,int iteration)
{
    Particle *particle;
    particle=D->particle;

    int i,s,intThe;
    int startI,endI,nSpecies;
    double dPhi,theta;
    ptclList *p;

    startI=1;  endI=D->subSliceN+1;
    nSpecies=D->nSpecies;
    dPhi=D->numSlice*2*M_PI;

    for(i=startI; i<endI; i++)
    {
      for(s=0; s<nSpecies; s++) {
        p=particle[i].head[s]->pt;
        while(p)  {
          theta=p->theta;
//			 intThe=(int)(theta/dPhi);
//          if(intThe>=1) { printf("In rearrange, intThe=%d,theta=%g,iteration=%d,dPhi=%g\n",intThe,theta,iteration,dPhi); }
          if(theta>=dPhi)  {
            intThe=(int)(theta/dPhi);
//	         theta-=dPhi*intThe;
	         theta-=dPhi;
          } else if(theta<0) {              
            intThe=(int)(theta/dPhi);
//	         theta-=dPhi*(intThe-1);
	         theta+=dPhi;
          } else ;

          if(theta>=dPhi || theta<0) { printf("In rearrange, theta=%g,iteration=%d,dPhi=%g\n",theta,iteration,dPhi);  }
          p->theta=theta;    

	       p=p->next;
	     }
      }
    }		//End of for(i)
}

