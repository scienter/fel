#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>


void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void save_Slice_Particle_HDF(double *data,char *fileName,char *groupName,char *dataName,int cnt,int dataCnt);
void save_Slice_Particle_Cnt_HDF(int *data,char *fileName,char *groupName,char *dataName,int cnt,int subCnt,int minI);

void save_Particle_HDF(Domain *D,int s,int sliceI,int iteration)
{
   int i,cnt,dataCnt=9,*ptclCnt;
   ptclList *p;
   double *data;
   char fileName[100],groupName[100],dataName[100];
   hid_t file_id,group_id;
   int myrank,nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   // initialize dataset
   p=D->head[s]->pt;
   cnt=0;
   while(p) {
     cnt++;
     p=p->next;
   }

   data=(double *)malloc(cnt*dataCnt*sizeof(double ));
   p=D->head[s]->pt;
   i=0;
   while(p) {
     data[i*dataCnt+0]=p->theta;
     data[i*dataCnt+1]=p->x;
     data[i*dataCnt+2]=p->y;
     data[i*dataCnt+3]=p->gamma;
     data[i*dataCnt+4]=p->px;
     data[i*dataCnt+5]=p->py;
     data[i*dataCnt+6]=p->weight;
     data[i*dataCnt+7]=p->sliceI;
     data[i*dataCnt+8]=p->index;
     i++;
     p=p->next;
   }

   sprintf(fileName,"particle%d.h5",iteration);
   sprintf(groupName,"%d",s);
   sprintf(dataName,"%d",sliceI);
   printf("sliceI=%d\n",sliceI);
   save_Slice_Particle_HDF(data,fileName,groupName,dataName,cnt,dataCnt);
   free(data);

   sprintf(dataName,"ptclCnt");
   save_Slice_Particle_Cnt_HDF(&cnt,fileName,groupName,dataName,D->sliceN,1,D->minI+sliceI);

}

void save_Slice_Particle_Cnt_HDF(int *data,char *fileName,char *groupName,char *dataName,int cnt,int subCnt,int minI)
{
   int myrank, nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   hid_t file_id,group_id,plist_id,dset_id,attr_id,as_id;
   hid_t dataspace,memspace;
   hsize_t dimsf[2],count[2],offset[2],block[2],stride[2];
   herr_t ierr;

   plist_id=H5Pcreate(H5P_FILE_ACCESS);
   ierr=H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

   //open file
   file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
   ierr=H5Pclose(plist_id);

   //open a group
   group_id=H5Gopen2(file_id,groupName,H5P_DEFAULT);
    
   //file space
   dimsf[0]=cnt;	//column
   dimsf[1]=1;		
   dataspace=H5Screate_simple(2,dimsf,NULL);

   // Create a dataset in group
   dset_id=H5Dcreate2(group_id,dataName,H5T_NATIVE_INT,dataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

   //memory space
   dimsf[0]=subCnt;
   dimsf[1]=1;
   memspace=H5Screate_simple(2,dimsf,NULL);
   stride[0]=1;
   stride[1]=1;
   count[0]=1;
   count[1]=1;

   //hyperslab in file space
   block[0]=subCnt;
   block[1]=1;
   offset[0]=minI;
   offset[1]=0;
   H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,offset,stride,count,block);

   //hyperslab in memory space
   offset[0]=0;
   H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset,stride,count,block);

   //write the dataset
   plist_id=H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
   H5Dwrite(dset_id,H5T_NATIVE_INT,memspace,dataspace,plist_id,data);

   H5Pclose(plist_id);
   H5Dclose(dset_id);
   H5Sclose(dataspace);
   H5Sclose(memspace);

   //close the file
//   MPI_Barrier(MPI_COMM_WORLD);
   H5Gclose(group_id);
   H5Fclose(file_id);
}

void save_Slice_Particle_HDF(double *data,char *fileName,char *groupName,char *dataName,int cnt,int dataCnt)
{
   int myrank, nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   hid_t file_id,group_id,plist_id,dset_id,attr_id,as_id;
   hid_t dataspace;
   hsize_t dimsf[2];
   herr_t ierr;

   plist_id=H5Pcreate(H5P_FILE_ACCESS);
   ierr=H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

   //open file
   file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
   ierr=H5Pclose(plist_id);

   //open a group
   group_id=H5Gopen2(file_id,groupName,H5P_DEFAULT);
    
   //file space
   dimsf[0]=cnt;	//column
   dimsf[1]=dataCnt;	//row
   dataspace=H5Screate_simple(2,dimsf,NULL);

   // Create a dataset in group
   dset_id=H5Dcreate2(group_id,dataName,H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

   //write the dataset
   plist_id=H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
   H5Dwrite(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);

   H5Pclose(plist_id);
   H5Dclose(dset_id);
   H5Sclose(dataspace);

   //close the file
//   MPI_Barrier(MPI_COMM_WORLD);
   H5Gclose(group_id);
   H5Fclose(file_id);
}

void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void createFile(Domain *D,int iteration)
{
   int s;
   hid_t file_id,group_id;
   char fileName[100],groupName[100];

   sprintf(fileName,"particle%d.h5",iteration);
   //create a file
   file_id=H5Fcreate(fileName,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
   for(s=0; s<D->nSpecies; s++) {
     sprintf(groupName,"%d",s);
     group_id = H5Gcreate2(file_id,groupName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     H5Gclose(group_id);
   }
   H5Fclose(file_id);
}
