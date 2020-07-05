// Author: Wes Kendall
// Copyright 2011 www.mpitutorial.com
// This code is provided freely with the tutorials on mpitutorial.com. Feel
// free to modify it for your own use. Any distribution of the code must
// either provide a link to www.mpitutorial.com or keep this header intact.
//
// An intro MPI hello world program that uses MPI_Init, MPI_Comm_size,
// MPI_Comm_rank, MPI_Finalize, and MPI_Get_processor_name.
//
#include <mpi.h>
#include <stdio.h>

/*int main(int argc, char** argv) {
  // Initialize the MPI environment. The two arguments to MPI Init are not
  // currently used by MPI implementations, but are there in case future
  // implementations might need the arguments.
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  // Print off a hello world message
  printf("Hello world from processor %s, rank %d out of %d processors\n",
         processor_name, world_rank, world_size);

  // Finalize the MPI environment. No more MPI calls can be made after this
  MPI_Finalize();
}*/

int main(int argc, char **argv)
{
  int rank, size, err;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  printf("rank=%d size=%d\n",rank,size);
  if(rank==0)
  {
  int n=4711;
  MPI_Send(&n,1,MPI_INT,1,112,MPI_COMM_WORLD);
  }
  if(rank==1)
  {
  int n=0;
  MPI_Status status;
  MPI_Recv(&n,1,MPI_INT,MPI_ANY_SOURCE,112,MPI_COMM_WORLD,&status);
  printf("n=%d\n",n);
  int version, subversion;
  MPI_Get_version(&version, &subversion);
  printf("version=%d\n",version);
  printf("subversion=%d\n",subversion);
  
  printf("source=%d\n",status.MPI_SOURCE);
  printf("tag=%d\n",status.MPI_TAG);
  printf("error=%d\n",status.MPI_ERROR);
  }
  // err=MPI_Barrier(MPI_COMM_WORLD);
  // double starttime=MPI_Wtime();
  // int n=2;
  // //printf("%d\n",n);
  // int sum;
  // MPI_Reduce(&n,&sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  // //printf(" Received reduce=%d\n",sum);
  // MPI_Barrier(MPI_COMM_WORLD);
  // MPI_Allreduce(&n,&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  // printf(" Received Allreduce =%d\n",sum);

  // //printf(" Received =%d\n",n);
  // double endtime=MPI_Wtime();
  MPI_Finalize();
  //printf("%f\t,%f\n",starttime,endtime);
  return 0;
}
/*{
int rank, size, err, n;
MPI_Init(&argc, &argv);  
int buffer[4]={0,1,2,3};
int count=4;
root=0;
err=MPI_Bcast(buffer, count, MPI_INT, root, MPI_COMM_WORLD);
n=(rank+1)*4711;

  //err=MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
  printf(" Received =%d\n",n);
  MPI_Finalize();
  return 0;
}*/