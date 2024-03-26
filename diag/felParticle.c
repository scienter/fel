#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <complex.h>


void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void restoreDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt);
void restore_attr_HDF(int *cnt,char *fileName,char *dataName,char *attrName);
void restore_Particle_HDF(double *data,char *fileName,char *dataName,int totalCnt,int subCnt,int start,int N);

void main(int argc, char *argv[])
{
   int initial,final,timeStep,division,skipCnt,step,i,j,sliceI,sliceN,s,totalCnt,n;
   int numData,nSpecies,subP,skip,*subCnt,*start;
   int ii,jj,indexI,indexJ,numXY,numGam,idxGam;
   double theta,z,x,y,gamma,px,py,weight,index,core,range,curr,charge;
   double n0,dz,minZ,w,sum,dPhi,tmp,wx[2],wy[2],minX,minY,dx,dy,den;
   double th,lambda0,bucketZ,minGam,rangeGam,dGam,gamRange,gamma0;
   double *data,*denZ,*gam,**denXY,**denGam,*gam0List;
   char fileName[100],dataName[100],attrName[100],outFile[100],outDen[100];
   FILE *out,*out2;

   if(argc<9) {
      printf("felParticle division initial final step skipCnt rangeXY numXY gamRange numGam\n");
      exit(0);
   } else ;

   division=atoi(argv[1]);
   initial=atoi(argv[2]);
   final=atoi(argv[3]);
   timeStep=atoi(argv[4]);
   skipCnt=atoi(argv[5]);
   range=atof(argv[6]);
   numXY=atoi(argv[7]);
   gamRange=atof(argv[8]);
   numGam=atoi(argv[9]);

   dGam=gamRange/(numGam*1.0);
   minX=minY=-0.5*range;
   dx=dy=range/(numXY*1.0);
   denXY=(double **)malloc((numXY+1)*sizeof(double *));
   for(i=0; i<=numXY; i++)
     denXY[i]=(double *)malloc((numXY+1)*sizeof(double ));

   for(step=initial; step<=final; step+=timeStep)
   {
      sprintf(fileName,"particle%d.h5",step);
      restoreIntMeta(fileName,"nSpecies",&nSpecies,1);
      printf("nSpecies=%d\n",nSpecies);
      gam0List=(double *)malloc(nSpecies*sizeof(double ));
      restoreIntMeta(fileName,"numData",&numData,1);
      restoreIntMeta(fileName,"sliceN",&sliceN,1);
      restoreDoubleMeta(fileName,"minZ",&minZ,1);
      restoreDoubleMeta(fileName,"gamma0",gam0List,nSpecies);
      restoreDoubleMeta(fileName,"dz",&dz,1);
      restoreDoubleMeta(fileName,"lambda0",&lambda0,1);      
      restoreDoubleMeta(fileName,"dPhi",&dPhi,1);      
      bucketZ=lambda0*dPhi*0.5/M_PI;

      denZ=(double *)malloc(sliceN*sizeof(double ));
      gam=(double *)malloc(sliceN*sizeof(double ));
      denGam=(double **)malloc(sliceN*sizeof(double *));
      for(i=0; i<sliceN; i++) { 
        denZ[i]=0.0; gam[i]=0.0; 
        denGam[i]=(double *)malloc(numGam*sizeof(double ));
      }

      //nSpecies=1;
      for(s=0; s<nSpecies; s++) {
        gamma0 = gam0List[s];
	minGam = gamma0-gamRange*0.5;

	for(i=0; i<=numXY; i++)
          for(j=0; j<=numXY; j++)
	    denXY[i][j]=0.0;
        for(i=0; i<sliceN; i++)
          for(j=0; j<numGam; j++)
            denGam[i][j]=0.0;

	sprintf(dataName,"%d",s);	
	restore_attr_HDF(&totalCnt,fileName,dataName,"totalCnt");
 
  	subCnt=(int *)malloc(division*sizeof(int ));
        start=(int *)malloc(division*sizeof(int ));
        subP=totalCnt/division;
        for(i=0; i<division-1; i++) 
   	  subCnt[i]=subP;
        subCnt[division-1]=totalCnt-subP*(division-1);
        start[0]=0;
        sum=0;
        for(n=0; n<division; n++) {
	  start[n]=sum;
          sum+=subCnt[n];
        }
        for(n=0; n<division; n++) 
          printf("particle counts : subCnt[%d]=%d, start[%d]=%d\n",n,subCnt[n],n,start[n]);

	sprintf(outFile,"%dParticle%d",s,step);	
        out=fopen(outFile,"w");

	skip=0;
	charge=0.0;
        for(n=0; n<division; n++) {
          data=(double *)malloc(subCnt[n]*numData*sizeof(double ));  

          restore_Particle_HDF(data,fileName,dataName,totalCnt,subCnt[n],start[n],numData);

	  for(i=0; i<subCnt[n]; i++) {
            // getting theta, and rearrange theta.
            theta=data[i*numData+0];
//	    tmp=theta/dPhi;
//            w=tmp-(int)tmp;
//            theta=w*dPhi;
//            if(theta>dPhi)   theta=theta-dPhi;
//            else if(theta<0) theta=dPhi+theta;
//	    if(theta>=dPhi || theta<0) { printf("theta=%g\n",theta); exit(0); }; 
	    // index of slice
            sliceI=data[i*numData+8];
            // calculating z
//            z=(sliceI+theta/dPhi)*bucketZ+step*dz+minZ;
            th=sliceI*dPhi+theta;            
            x=data[i*numData+1];
            y=data[i*numData+2];
            gamma=data[i*numData+3];
            px=data[i*numData+4];
            py=data[i*numData+5];
            weight=data[i*numData+6];
            index=data[i*numData+7];
            core=data[i*numData+9];
	    charge+=weight;

	    denZ[sliceI]+=weight;
	    gam[sliceI]+=gamma*weight;
            idxGam=(int)((gamma-minGam)/dGam);
	    if (idxGam>=0 && idxGam<numGam-1)
              denGam[sliceI][idxGam]+=1;

	    // calculation for density at beam center
	    if(sliceI==sliceN/2) {
              indexI=(int)((x-minX)/dx);
              indexJ=(int)((y-minY)/dy);
  	      if(indexI>=0 && indexI<numXY && indexJ>=0 && indexJ<numXY) {
	        wx[1]=(x-minX)/dx-indexI; wx[0]=1.0-wx[1];
	        wy[1]=(y-minY)/dy-indexJ; wy[0]=1.0-wy[1];
	        for(ii=0; ii<2; ii++)
	          for(jj=0; jj<2; jj++)
		    denXY[indexI+ii][indexJ+jj]+=wx[ii]*wy[jj]*weight;
  	      } else ;
	    } else ;

	    if(skip%skipCnt==0)
              fprintf(out,"%g %g %g %.15g %g %g %g %.9g %d %g\n",th,x,y,gamma,px,py,weight,index,sliceI,core);
	    else ;
	    skip++;
	  }

	  free(data);
	  printf("step=%d, division status is %d/%d.\n",step,n,division);
        }
       	fclose(out);
	printf("%s is made.\n",outFile);
        free(subCnt);
        free(start);

	// Den file is saved.
	sprintf(outDen,"%dDen%d",s,step);	
        out2=fopen(outDen,"w");
	sum=0.0;
        for(ii=0; ii<=numXY; ii++) {
	  x=minX+ii*dx;
          for(jj=0; jj<=numXY; jj++) {
	    y=minY+jj*dy;
	    den=denXY[ii][jj]/dx/dy/bucketZ;
	    sum+=denXY[ii][jj]*1.602e-19;
	    fprintf(out2,"%g %g %g\n",x,y,den);
	  }
	  fprintf(out2,"\n");
	}	
       	fclose(out2);
	curr=sum/bucketZ*3e8;
	printf("%s is made. peakCurrent=%g\n",outDen,curr);

        sprintf(outFile,"%dDensityGam%d",s,step);
        out=fopen(outFile,"w");
        for(i=0; i<sliceN; i++) {
          for(j=0; j<numGam; j++) {
            z=i*bucketZ+step*dz+minZ;
            gamma = minGam + dGam*j;
            fprintf(out,"%.10g %g %g\n",z,gamma,denGam[i][j]);
          }
          fprintf(out,"\n");
        }
        fclose(out);
        printf("%s is made.\n",outFile);

	
      }		//End of for(s)

      sprintf(outFile,"ptclSlice%d",step);	
      out=fopen(outFile,"w");
      for(i=0; i<sliceN; i++) {
        z=i*bucketZ+step*dz+minZ;
	n0=denZ[i]/bucketZ;
	curr=n0*1.602e-19*3e8;
	if(n0==0.0) gamma=0.0;
	else        gamma=gam[i]/n0;
	if(n0>0.0) fprintf(out,"%.10g %g %g %g\n",z,curr,n0,gamma);
	else ;
      }
      fclose(out);
      printf("%s is made. peakCurrent=%g\n",outFile,curr);

      free(denZ);
      free(gam);
      for(i=0; i<sliceN; i++) free(denGam[i]);
      free(denGam);
      free(gam0List);
      
   }

   for(i=0; i<=numXY; i++) free(denXY[i]); free(denXY);

   printf("total charge = %g [pC]\n",charge*1.602e-7);

}

void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void restoreDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void restore_Particle_HDF(double *data,char *fileName,char *dataName,int totalCnt,int subCnt,int start,int N)
{
   hid_t file_id,group_id,plist_id,dset_id,attr_id,as_id;
   hid_t dataspace,memspace;
   hsize_t dimsf[2],offset[2],count[2];
   herr_t ierr;

   //open file
   file_id=H5Fopen(fileName,H5F_ACC_RDONLY,H5P_DEFAULT);

   //set dataset
   dset_id=H5Dopen(file_id,dataName,H5P_DEFAULT);

   dataspace=H5Dget_space(dset_id);

   // define hyperslab in the dataset
   offset[0]=start;  offset[1]=0;
   count[0]=subCnt;  count[1]=N;
   H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,offset,NULL,count,NULL);

   //memory dataspace
   dimsf[0]=subCnt;
   dimsf[1]=N;
   memspace=H5Screate_simple(2,dimsf,NULL);

   // define memory hyperslab
   offset[0]=0;  offset[1]=0;
   count[0]=subCnt;  count[1]=N;
   H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset,NULL,count,NULL);

   //read the dataset
   H5Dread(dset_id,H5T_NATIVE_DOUBLE,memspace,dataspace,H5P_DEFAULT,data);

//   H5Pclose(plist_id);
   H5Dclose(dset_id);
   H5Sclose(memspace);
   H5Sclose(dataspace);
   //close the file
//   H5Gclose(group_id);
   H5Fclose(file_id);
}

void restore_attr_HDF(int *cnt,char *fileName,char *dataName,char *attrName)
{
   hid_t file_id,group_id,plist_id,dataset_id,attribute_id,dataspace_id;
   hsize_t dims;
   herr_t ierr;

   //open file
   file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
//   ierr=H5Pclose(plist_id);

   //open a group
//   group_id=H5Gopen2(file_id,groupName,H5P_DEFAULT);

   //Open an existing dataset.
   dataset_id = H5Dopen2(file_id,dataName,H5P_DEFAULT);

   // Create dataspace for attribute  
//   dims=1;
//   dataspace_id=H5Screate_simple(1,&dims,NULL);

   // open a dataset attribute
   attribute_id=H5Aopen(dataset_id,attrName,H5P_DEFAULT);

   // read the dataset
   H5Aread(attribute_id,H5T_NATIVE_INT,cnt);

   H5Aclose(attribute_id);
//   H5Sclose(dataspace_id);
   H5Dclose(dataset_id);
//   H5Gclose(group_id);
   H5Fclose(file_id);
}
