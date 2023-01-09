#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>


//void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
//void save_Slice_Particle_HDF(double *data,char *fileName,char *groupName,char *dataName,int cnt,int dataCnt);
void restore_ptclCnt_attr_HDF(int *cnt,char *fileName,char *groupName,char *dataName,char *attrName);
void restore_Slice_Particle_HDF(double *data,char *fileName,char *groupName,char *dataName,int cnt,int dataCnt);
void restoreFieldComp(double complex ****data,char *fileName,char *dataName,int h,int sliceN,int N,int subSliceN,int step,int minI);

void restore_Field_HDF(Domain *D,int iteration)
{
   int i,h,N,subSliceN,sliceN;
   double *data;
   char fileName[100],dataName[100],attrName[100];
   hid_t file_id,group_id;
   int myrank,nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   N=D->nx*D->ny; subSliceN=D->subSliceN; sliceN=D->sliceN;

   sprintf(fileName,"dump/field%d.h5",iteration);
   for(h=0; h<D->harmony; h++) {
     sprintf(dataName,"U%d",h);
     restoreFieldComp(D->U,fileName,dataName,h,sliceN,N,subSliceN,D->saveStep,D->minI);
     MPI_Barrier(MPI_COMM_WORLD);
     sprintf(dataName,"Ez%d",h);
     restoreFieldComp(D->Ez,fileName,dataName,h,sliceN,N,subSliceN,D->saveStep,D->minI);
     MPI_Barrier(MPI_COMM_WORLD);
   }
}

void restoreFieldComp(double complex ****data,char *fileName,char *dataName,int h,int sliceN,int N,int subSliceN,int step,int minI)
{
   int start,i,j;
   double *field,realV,imagV;
   double complex val;
   int myrank, nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   hid_t file_id,dset_id,plist_id,tic_id;
   herr_t status;
   hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
   hsize_t dimsf[2],count[2],offset[2];

   plist_id=H5Pcreate(H5P_FILE_ACCESS);
   H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

   file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
   H5Pclose(plist_id);

   dimsf[0]=sliceN;
   dimsf[1]=N*2;
   filespace=H5Screate_simple(2,dimsf,NULL);

   count[0]=subSliceN;
   count[1]=N*2;
   offset[0]=minI;
   offset[1]=0;
   memspace=H5Screate_simple(2,count,NULL);

   field = (double *)malloc(subSliceN*2*N*sizeof(double ));

   dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
   subfilespace=H5Dget_space(dset_id);
   H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);

   plist_id=H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
   status = H5Dread(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);

   start=0;
   for(i=0; i<subSliceN; i++) {
     for(j=0; j<N; j++)  {
       realV=field[start+2*j+0];
       imagV=field[start+2*j+1];

       data[0][h][i][j]=realV+I*imagV;
     }
     start+=2*N;
   }

   H5Pclose(plist_id);
   H5Sclose(subfilespace);
   H5Dclose(dset_id);

   H5Sclose(memspace);
   H5Sclose(filespace);
   H5Fclose(file_id);
   free(field);
}


void restore_Particle_HDF(Domain *D,int s,int sliceI,int iteration)
{
   int i,cnt,dataCnt=9,*ptclCnt;
   ptclList *p;
   double *data;
   char fileName[100],groupName[100],dataName[100],attrName[100];
   hid_t file_id,group_id;
   int myrank,nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   sprintf(fileName,"dump/particle%d_%d.h5",iteration,myrank);
   sprintf(groupName,"%d",s);
   sprintf(dataName,"%d",sliceI);
   sprintf(attrName,"ptclCnt");
   restore_ptclCnt_attr_HDF(&cnt,fileName,groupName,dataName,attrName);

   data=(double *)malloc(cnt*dataCnt*sizeof(double ));
   restore_Slice_Particle_HDF(data,fileName,groupName,dataName,cnt,dataCnt);

   for(i=0; i<cnt; i++) {
     p=(ptclList *)malloc(sizeof(ptclList));
     p->next=D->head[s]->pt;
     D->head[s]->pt=p;

     p->theta=data[i*dataCnt+0];
     p->x=data[i*dataCnt+1];
     p->y=data[i*dataCnt+2];
     p->gamma=data[i*dataCnt+3];
     p->px=data[i*dataCnt+4];
     p->py=data[i*dataCnt+5];
     p->weight=data[i*dataCnt+6];
     p->sliceI=data[i*dataCnt+7];
     p->index=data[i*dataCnt+8];
   }
   free(data);

}


void restore_Slice_Particle_HDF(double *data,char *fileName,char *groupName,char *dataName,int cnt,int dataCnt)
{
   int myrank, nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   hid_t file_id,group_id,plist_id,dset_id,attr_id,as_id;
   hid_t dataspace,memspace;
   hsize_t dimsf[2],offset[2],stride[2],block[2],count[2];
   herr_t ierr;

   //open file
   file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);

   //open a group
   group_id=H5Gopen2(file_id,groupName,H5P_DEFAULT);

   //set dataset
   dset_id=H5Dopen2(group_id,dataName,H5P_DEFAULT);

   //write the dataset
//   plist_id=H5Pcreate(H5P_DATASET_XFER);
//   H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
//   H5Dwrite(dset_id,H5T_NATIVE_DOUBLE,memspace,dataspace,plist_id,data);
   H5Dread(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);

//   H5Pclose(plist_id);
   H5Dclose(dset_id);

   //close the file
   H5Gclose(group_id);
   H5Fclose(file_id);
}

void restore_ptclCnt_attr_HDF(int *cnt,char *fileName,char *groupName,char *dataName,char *attrName)
{
   int myrank, nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   hid_t file_id,group_id,plist_id,dataset_id,attribute_id,dataspace_id;
   hsize_t dims;
   herr_t ierr;

   //open file
   file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
   ierr=H5Pclose(plist_id);

   //open a group
   group_id=H5Gopen2(file_id,groupName,H5P_DEFAULT);

   //Open an existing dataset.
   dataset_id = H5Dopen2(group_id,dataName,H5P_DEFAULT);

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
   H5Gclose(group_id);
   H5Fclose(file_id);
}

