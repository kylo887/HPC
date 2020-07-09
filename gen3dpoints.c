
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include <math.h>

//  int main()
// //int plot()
// {
//   printf("x coordinates\n");
//    int x[101];
//    int counterx,countery,counterz; //counter variable
//    FILE *fp;
//    fp = fopen("coord.dat", "w+");
//    if(fp==NULL)
//    {
//      printf("\nError:file cannot be opened");
//      exit(1);
//    }
//    fprintf(fp, "This is testing for fprintf...\n");
//    fprintf(fp,"x coordinates\n");
//    //fputs("This is testing for fputs...\n", fp);
//    //fclose(fp);
//    for( counterx = 1; counterx <= 100; ++counterx )
//    {
//       x[counterx]= (rand()%100);
//       printf("%5d", x[counterx]);
//       fprintf(fp,"%5d",x[counterx]); //for the output in a file
//       //for displaying value in new line
//       if (counterx % 25 == 0)
//       puts("");
//       fputs("\n",fp);
//    }
   
//    fputs("\n",fp);
//    //printf("%d\n",counterx);
//    printf("%d\n",x[2]);

//    printf("y coordinates\n");
//    fprintf(fp,"y coordinates\n");
//    int y[101];
//    for (countery=1;countery<=100;++countery)
//    {
//      y[countery]=(rand()%100);
//      printf("%5d",y[countery]);
//      fprintf(fp,"%5d",y[countery]);
//      if (countery%25==0)
//      puts("");
//      fputs("\n",fp);
//    }
//    fputs("\n",fp);
//    printf("%d\n",y[2]);
//    printf("8 days later, comment this later\n");
//    printf("blah blah blah\n");

//    printf("z coordinates\n");
//    fprintf(fp,"z coordinates\n");
//    int z[101];
//    for (counterz=1;counterz<=100;++counterz)
//    {
//      z[counterz]=(rand()%100);
//      printf("%5d",z[counterz]);
//      fprintf(fp,"%5d",z[counterz]);
//      if (counterz%25==0)
//      puts("");
//      fputs("\n",fp);
//    }
//    printf("%d\n",z[2]);
//    fclose(fp);
//    system("gnuplot -p file.gp");
//    return 0;
   
// }

int number_of_particles,number_of_particles_default=100;
float total_time,total_time_default=60;
float time_step,time_step_default=0.02;
//printf("%d\t%d\t%f",number_of_particles_default,total_time_default,time_step_default);

typedef struct box {
	float x, y , z;	// position
	float ax, ay , az;	// acceleration
	float vx, vy , vz;	// velocity
	float mass;	// mass
} box;


int lengths[3] = { 1, 1, 1 };
const MPI_Aint displacements[3] = { 0, sizeof(float), sizeof(float)*2 };
MPI_Datatype types[3] = { MPI_FLOAT, MPI_FLOAT, MPI_FLOAT }; /* to create MPI data types for position,velocity and acceleration */
MPI_Datatype POSITION;
MPI_Datatype VELOCITY;

int start,end;
void pr(int number_of_particles,int size , int rank ,int *start ,int *end,int *sendcounts,int *displs)
{
int r = (number_of_particles)%size;
int s = 0;

for (int p = 0; p < size; p++)
    {
         sendcounts[p] = ((number_of_particles)/size); /* array of number of particles assigned for current rank*/
         if (r > 0) {
            sendcounts[p]=sendcounts[p]+1;
            r--;}

         displs[p] = s; /*start boundary of particles in current rank */
         s += sendcounts[p];
    }

if (rank==(size-1))
{
    *start = displs[rank];
    *end = number_of_particles;
}
else
{
    *start = displs[rank]; /* start n end in case of parallelizing loops*/
    *end = displs[rank+1];
}
}
void integrate(box *box) {
        box->vx += box->ax * time_step;
        box->vy += box->ay * time_step;
        box->vz += box->az * time_step;

        box->x += box->vx * time_step +(1 / 2) * box->ax * time_step * time_step;
        box->y += box->vy * time_step+ (1 / 2) * box->ay * time_step * time_step;

        box->z += box->vz * time_step + (1 / 2) * box->az * time_step * time_step;

int x,y,z,delx,dely,delz,inverse,inversecube,f,mass,xnew,ynew,znew,vx,vy,vz,i,j eps;
void loop()
{
  for (i=0,i<number_of_particles,i++)
  {
    float ax=0.0;
    float ay=0.0;
    float az=0.0;
    for (j=0,j<number_of_particles,j++)
    {
      delx=x[j]-x[i];
      dely=y[j]-y[i];
      delz=z[j]-z[i];
      inverse=1/sqrt(delx**2+dely**2+delz**2+eps);
      inversecube=inverse**2
      f=mass[j]*inversecube;
      ax+=f*delx;
      ay+=f*dely;
      az+=f*delz;
    }
    xnew[i]=x[i]+time_step*vx[i]+0.5*time_step*time_step*ax;
    ynew[i]=y[i]+time_step*vy[i]+0.5*time_step*time_step*ay;
    znew[i]=z[i]+time_step*vz[i]+0.5*time_step*time_step*az;
    vx[i]+=time_step*ax;
    vy[i]+=time_step*ay;
    vz[i]+=time_step*az;
  }
  for (i=0,i<=number_of_particles,i++)
  {
    x[i]=xnew[i];
    y[i]=ynew[i];
    z[i]=znew[i];
  }
}

void write(position)
{
    FILE *dd = fopen("pt.txt", "a");
    for(int i=0; i<number_of_particles; i++)
    {
        fprintf(dd,"%f %f %f \n",position[i].x,position[i].y,position[i].z);
    }
    fprintf(dd,"\n");
    fclose(dd);
}
void sim(position,velocity)
{
    MPI_Bcast(position,number_of_particles,POSITION,0,MPI_COMM_WORLD);
    MPI_Bcast(velocity,number_of_particles,VELOCITY,0,MPI_COMM_WORLD);

    int iterations = total_time / time_step;

    for(int t=0;t<iterations;t++)
    {
        loop()

        //#MPI_Allgatherv(position+(displs[rank]),sendcounts[rank],POSITION,position,sendcounts,displs,POSITION,MPI_COMM_WORLD);
        //#MPI_Allgatherv(velocity+(displs[rank]),sendcounts[rank],VELOCITY,velocity,sendcounts,displs,VELOCITY,MPI_COMM_WORLD);

    if(rank==0)
       { write_p();}
    }
}
int main(int argc, char **argv)
{
   srand(time(0));
   MPI_Init(&argc, &argv); /* Initiate MPI*/

   MPI_Barrier(MPI_COMM_WORLD);
   double tcal = MPI_Wtime(); /* Computation start time*/

   if (argc >= 2)
      Np = atoi(argv[1]);
   else
      Np = Np_d;

   if (argc >= 3)
      total_time = atoi(argv[2]);     /* user input through command line arguments */
   else
      total_time = time_d;

   if (argc >= 4)
      dt = atof(argv[3]);
   else
      dt = dt_d;

   MPI_Comm_size(MPI_COMM_WORLD,&size); /* size and rank of MPI*/
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   MPI_Type_create_struct(3, lengths, displacements, types, &POSITION);
   MPI_Type_commit(&POSITION);                                         /*committing user defined MPI data types*/

   MPI_Type_create_struct(3, lengths, displacements, types, &VELOCITY);
   MPI_Type_commit(&VELOCITY);

   sendcounts = malloc(sizeof(int)*size); /* allocating memory for arrays*/
   displs = malloc(sizeof(int)*size);

   P_r(Np,size,rank,&start,&end,sendcounts,displs); /* calling function to calculate send counts and displacements*/

   position = (Pos *) malloc(Np * sizeof(Pos));   /* allocating memory for structs*/
   velocity = (Vel *) malloc(Np * sizeof(Vel));
   a        = (Acc *) malloc(sendcounts[rank] * sizeof(Acc));

   if (rank == 0)
        {
          generate_pos_vel(); /* initializing position and velocity in root process*/

         FILE *fp = fopen("pos.txt", "w"); /*writing initial positions*/
         for(int i=0; i<Np; i++)
           {
                fprintf(fp,"%f %f %f \n",position[i].x,position[i].y,position[i].z);
           }
         fprintf(fp,"\n");
         fclose(fp);
        }

   sim(); /* function call for particle simulation */

   MPI_Barrier(MPI_COMM_WORLD);
   double elapsedt = MPI_Wtime() - tcal; /* computational time taken to complete the process*/

   if (rank==0)
   	{
            printf("Np=%d,Total time= %f,dt=%f \n",Np,total_time,dt);
        	printf("run time = %f \n",elapsedt);
 	}

   MPI_Finalize();
}

