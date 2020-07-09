#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include<math.h>
#include <mpi.h>

#define RAND_MAX       10000

void outer(int argc, char *argv[])
{
    int size, rank; //variables for size and rank
	int i,j,n=100, timeduration=20;   //Number of particles and time duration
	
    MPI_Init (&argc,&argv);      // Initialize the MPI environment

    MPI_Comm_size (MPI_COMM_WORLD, &size); //gives the number of processes in the group comm
    MPI_Comm_rank (MPI_COMM_WORLD, &rank); //gives the rank of the calling process in the group comm

    struct particleObj   //declaring the  structure of positions,velocities and mass
	{
	  double arrx[n];
	  double arry[n];
	  double arrz[n];
	  double mass[n];
	  double vx[n];
	  double vy[n];
	  double vz[n];

	};
    double arrx[n] , arry[n],arrz[n], dt=0.5;  //positions for x,y,z and step time
    double delta_x=0, delta_y=0, delta_z=0, ax=0,ay=0,az=0, x=0 , y=0, z=0, invr3=0; //initalize positions and other variables to zero
    double invr=0;
    double force, mass[n], xnew[n], vx[n],vy[n], vz[n], ynew[n],znew[n];  //force, mass and new positions
    int t;
    t=(int)(timeduration/dt);
    int sizen_particle;  
    sizen_particle = n / size;  // creating sizen_particle from n for mpi split
  
    
    double time_parallel_comp = 0;  //initalize time to zero
   struct particleObj *p1;  //creating objects of structure p1
   p1=malloc(sizeof(*p1) * sizen_particle); //allocate memory for p1
   //double startTime = MPI_Wtime();
  
   int pIndex= rank* sizen_particle;  //determine index into array structure for each process

  //initialsation
   if (rank == 0) {
	  time_parallel_comp -= MPI_Wtime();
     
           //initialsation of the structure 
		for(i=0;i<n;i++){

		   p1->arrx[i]= arrx[i]= ((double)rand())/RAND_MAX;
		   p1->arry[i]= arry[i]= ((double)rand())/RAND_MAX;
		   p1->arrz[i]= arrz[i]= ((double)rand())/RAND_MAX;
		   p1->mass[i]=  mass[i]= 0.1;
		   p1->vx[i]=  vx[i]= ((double)rand())/RAND_MAX;
		   p1->vy[i]=  vy[i]=  ((double)rand())/RAND_MAX;
		   p1->vz[i]=  vz[i]=   ((double)rand())/RAND_MAX;

	    }	   
   }
   struct particleObj *values,*splitParticles,*gatherParticles; //creating objects of structure
   splitParticles=  malloc(sizeof(*splitParticles) * sizen_particle); //allocate memory for splitParticles
   gatherParticles=  malloc(sizeof(*gatherParticles) * sizen_particle);  //allocate memory for gatherParticles

      for(i=0;i<n;i++){
	   splitParticles->arrx[i]= 0;
           splitParticles->arry[i]= 0;   //initialise split particles to zero
           splitParticles->arrz[i]= 0;
           splitParticles->mass[i]= 0;
           splitParticles->vx[i]=   0;
           splitParticles->vy[i]=   0;
           splitParticles->vz[i]=   0;
     }
    
    FILE *fp;     //for creating a file
    char* str = "string";
    int e = 10;

    fp=fopen("test.txt", "w");  //output file text created
    if(fp == NULL)
        exit(-1);

     int count=0;
	 
	// Create and commit new MPI datatype
	MPI_Datatype tpos;
    MPI_Type_contiguous(3, MPI_DOUBLE, &tpos);  //creates contiguos datatype
    MPI_Type_commit(&tpos);
	MPI_Barrier (MPI_COMM_WORLD);  //Blocks until all the process in the communicator have reached this routine
     	 
    MPI_Scatter(p1,sizen_particle,tpos,
           splitParticles,sizen_particle,tpos,0,   //sends data p1 containing positions to all other processes in the communicator
           MPI_COMM_WORLD
  	  );

    for(i=0;i<n;i++)
    {
          ax = 0.0; //Initialising accelertions
          ay = 0.0;
          az = 0.0;
          //fprintf(fp, "%f\t%f\t%f\t", arrx[i], arry[i], arrz[i]);

        for(j=0;j<n;j++)
        {

          delta_x = (arrx[j] - arrx[i]); //difference between two  x positions
          delta_y = (arry[j] - arry[i]); //difference between two  y positions
          delta_z = (arrz[j] - arrz[i]); //difference between two  z positions
          //printf("%f\n", delta_x);        prints delta_x positions
          invr= 1.0/sqrt((delta_x*delta_x) + (delta_y*delta_y)+(delta_z*delta_z)+ 0.3);
          invr3=pow(invr, 3); //raise power of  invr 
          force= mass[j]*invr3; // calculate the force
          ax= ax + force*delta_x; //acceleration for x particles
          ay= ay + force*delta_y; //acceleration for y particles
          az= az + force*delta_z; //acceleration for z particles
          //printf("ax   %f\n", invr3);  
          
          fprintf(fp, "%f\t%f\t%f\t", arrx[j], arry[j], arrz[j]);   // Writes the current positions of the particles to file
          fprintf(fp,"\n");
        }
		
        xnew[i] = (arrx[i]) + (dt*vx[i]) + (0.5*dt*dt*ax); // calculates the new x positions
        ynew[i] = (arry[i]) + (dt*vy[i]) + (0.5*dt*dt*ay); // calculates the new y positions
        znew[i] = (arrz[i]) + (dt*vz[i]) + (0.5*dt*dt*az); // calculates the new z positions
        
       // printf("%f\n", xnew[i]); prints updated positions
        vx[i] = (vx[i])+ (dt*ax); // calculates the velocity in x direction
        vy[i] = (vy[i])+ (dt*ay); // calculates the velocity in y direction
        vz[i] = (vz[i])+ (dt*az); // calculates the velocity in z direction
        //printf("%f\n", vx[i]); prints velocities
        arrx[i]=xnew[i];
        arry[i]=ynew[i]; //original positions become new positions calculated
        arrz[i]=znew[i];

       MPI_Gather(gatherParticles, sizen_particle, tpos, 
                    p1, sizen_particle, tpos,0, MPI_COMM_WORLD);  // gathers the value gatherParticles from group of processes
        if (rank == 0) {
	time_parallel_comp += MPI_Wtime();   //returns time elapsed on calling processor
    		}  
 
    }
    fclose(fp);
 
    printf("time %f\n",time_parallel_comp); //print the time for parallel completion
    MPI_Finalize ();  //terminates the MPI execution environment
}


int main(int argc, char *argv[]){
    printf("main\n");
    outer(argc,argv);
}
