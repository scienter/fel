#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <hdf5.h>
#include <hdf5_hl.h>

int main(int argc, char *argv[])
{
    int s,iteration=0,rnk,i,startI,endI,subMaxStep,tmpInt,subIter;
    int progress,testProgress,maxProgress;
    double time_spent,dPhi,bucketZ,shiftZ;
    clock_t begin,end;
	 struct tm *t_now;
	 time_t timer;  //measure time
    char fileName[100],name[100];
    FILE *out;
    Domain D; 
    LoadList *LL; 
    hid_t file_id;
    int myrank, nTasks;
    MPI_Status status; 

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    begin=clock();

    if(argc < 2) 
    {  
      printf("mpirun -np N show [inputFile] [dumpNum]\n"); 
      exit(0); 
    }

    timer=time(NULL);
    t_now=localtime(&timer);

    //parameter setting
    parameterSetting(&D,argv[1]);
//printf("myrank=%d,iteration=%d, parameterSetting is done\n",myrank,iteration);

    //create mesh
    boundary(&D);
//printf("myrank=%d,iteration=%d, boundary is done\n",myrank,iteration);
 
//    removeFile(&D);

    loadSeed(&D,iteration);

    if(myrank==0) testK_quadG(&D); else ;
    MPI_Barrier(MPI_COMM_WORLD);

    //loading  beam
    iteration=0;
    LL=D.loadList;
    s=0;
    while(LL->next) {
      loadBeam(D,LL,s,iteration);
      LL=LL->next;
      s++;
    }

    // setting the wake_function
    wakeFunction(&D,iteration);
//printf("myrank=%d,iteration=%d, wakefunction is done\n",myrank,iteration);

	 shiftZ=0.0;
    while(iteration<D.maxStep) 
    {
      // updating wake_field
      if(iteration%D.wakeFieldStep==0) updateWakeField(&D,iteration); else ;
//printf("myrank=%d,iteration=%d, updateWakeField is done\n",myrank,iteration);

      // update total energy
      updateTotalEnergy(&D,iteration);
//printf("myrank=%d,iteration=%d, updateTotalEnergy is done\n",myrank,iteration);


      // save files
      if(iteration>=D.saveStart) {
        if(iteration%D.saveStep==0 || iteration==D.maxStep-1) {
	       if(D.particleSave==ON)  saveParticleHDF(&D,iteration); else ;
          if(D.fieldSave==ON)     saveFieldHDF(&D,iteration); else ;

          end=clock();
		    time_spent=(end-begin)/CLOCKS_PER_SEC/60.0;
		    if(myrank==0) {
		      printf("Time duration at %ddump:%.4gmin\n",iteration,time_spent);
		    } else ;

        }      else ;
      }  else ;


      if(D.mode==Twiss)  ;
      else  solveField(D,iteration);
//printf("myrank=%d,iteration=%d, solveField is done\n",myrank,iteration);

		//Chicane
      chicane_test(&D,iteration);
      if(D.calChicaneFlag==ON) {
       	calParticleDelay(&D);
       	rearrangeChicaneParticle(&D);
       	shiftChicaneField(&D);
			if(myrank==0) printf("Chicane is performed. at step=%d\n",iteration); else ;
      } else {
			
			set_chicane_zero(&D);		

			shiftZ+=(1.0*D.numLambdaU)/(1.0*D.numSlice);
			if(shiftZ>=0) { 
      		D.shiftFlag=ON; shiftZ-=1.0; 
        		if(D.mode==Time_Dependent)  shiftField(D,iteration); else ;
			} else ;
//printf("myrank=%d,iteration=%d, shiftField is done\n",myrank,iteration);

      	updateK_quadG(&D,iteration,0);
//printf("myrank=%d,iteration=%d, 1st updateK is done\n",myrank,iteration);

      	// twiss parameters
      	if(D.mode==Static || D.mode==Twiss) {
         	calculate_twiss(&D,iteration);
      	} else ;
//printf("myrank=%d,iteration=%d, calculate_twiss is done\n",myrank,iteration);


      	transversePush(&D,iteration);
//printf("myrank=%d,iteration=%d, 1st transversePush is done\n",myrank,iteration);

      	updateK_quadG(&D,iteration,0.5);
//printf("myrank=%d,iteration=%d, 2st updateK is done\n",myrank,iteration);

			if(D.shiftFlag==ON) push_theta_gamma(&D,iteration); else ;
//printf("myrank=%d,iteration=%d, push_theta_gamma is done\n",myrank,iteration);

      	transversePush(&D,iteration);
//printf("myrank=%d,iteration=%d, 2st transversePush is done\n",myrank,iteration);

      	//phase shifter
      	phaseShift(&D,iteration);

      	periodicParticles(&D,iteration);
//printf("myrank=%d,iteration=%d, periodicParticles is done\n",myrank,iteration);

//         rearrangeParticles(&D,iteration);

		}	//End of chicaneFlag==OFF



      if(myrank==0 && iteration%10==0) printf("iteration=%d\n",iteration); else ;
      
      iteration+=1;
    }		//End of for(iteration<maxStep)

    // save total energy for all steps
    saveTotalEnergy(&D);
    if(myrank==0) {
//      if(D.mode==Static || D.mode==Twiss) 
	      save_twiss(&D);
//      else;
    } else ;

    end=clock();
	 time_spent=(end-begin)/CLOCKS_PER_SEC/60.0;
	 if(myrank==0) {
	   //fprintf(out,"Time duration at %ddump:%.4gmin\n",iteration,time_spent);
	   printf("Time duration at %ddump:%.4gmin\n",iteration,time_spent);
    } else ;

    cleanMemory(&D);

    MPI_Finalize();

    return 0;
}
