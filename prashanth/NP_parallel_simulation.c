#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>


typedef struct _psystem {
	float x, y , z;	// position
	float ax, ay , az;	// acceleration
	float vx, vy , vz;	// velocity
	float mass;	// mass
} psystem;

void integrate(psystem *psystem, float deltaTime) {
        psystem->vx += psystem->ax * deltaTime;
        psystem->vy += psystem->ay * deltaTime;
        psystem->vz += psystem->az * deltaTime;

        psystem->x += psystem->vx * deltaTime +(1 / 2) * psystem->ax * deltaTime * deltaTime;
        psystem->y += psystem->vy * deltaTime+ (1 / 2) * psystem->ay * deltaTime * deltaTime;

        psystem->z += psystem->vz * deltaTime + (1 / 2) * psystem->az * deltaTime * deltaTime;

//output positions are saved in a textfile
FILE *output;
output= fopen("outputpositions.txt","a");
//fprintf(output,"\n");
//fprintf(output,"%f\t %f\t %f\n", psystem->x,psystem->y,psystem->z );
//printf("%f\t %f\t %f\n", psystem->x,psystem->y,psystem->z );
fclose(output);

}

void acceleration(psystem *a, psystem *b, float *ax, float *ay , float *az) {
	float eps = 10000; //10000;
	float distanceX = fabs(b->x - a->x);
        float distanceY =fabs(b->y - a->y);
        float distanceZ = fabs(b->z - a->z);
        float vectorDistance = sqrt(a->x * a->x + a->y * a->y+a->z * a->z+ eps);
        float vectorDistanceCubed = vectorDistance * vectorDistance * vectorDistance;
        float inverse = 1.0 / sqrt(vectorDistanceCubed);
        float scale = b->mass * inverse;
        *ax = (distanceX * scale);
        *ay = (distanceY * scale);
        *az = (distanceZ * scale);
}

void simulate(int rank, int totalparticles, int nparticles, psystem *particles, psystem *local_particles, float dt) {
        for(size_t i = 0; i < nparticles; i++) {
                float total_ax = 0, total_ay = 0, total_az = 0;
                for (size_t j = 0; j < totalparticles; j++) {
                        if (i == nparticles * rank + i) {
                                continue;
                        }
                        float ax, ay , az;
                        acceleration(&local_particles[i], &particles[j], &ax, &ay , &az);
                        total_ax += ax;
                        total_ay += ay;
                        total_az += az;

                }
                local_particles[i].ax = total_ax;
                local_particles[i].ay = total_ay;
                local_particles[i].az = total_az;

                integrate(&local_particles[i], dt);
        }
}

psystem *initparticles (int nparticles) {
        srand(time(NULL));
        psystem *particles = (psystem *)malloc(sizeof(*particles) * nparticles);
        for (int i = 0; i < nparticles; i++) {
                float initialMass = 2;
                psystem object = {
                        .x = ((double) rand() / (RAND_MAX))*2 - 1.15, 
                        .y = ((double) rand() / (RAND_MAX))*2 - 1.15,
                        .z = ((double) rand() / (RAND_MAX))*2 - 1.15,
                        .vx = (((double) rand() - 0.5) * 0.5)*1.2*0.0002,
                        .vy = (((double) rand() - 0.5) * 0.5)*2.1*0.0002,
                        .vz = (((double) rand() - 0.5) * 0.5)*0.8*0.0002,
                        .mass = (double) rand() + initialMass * 0.5
                };
        particles[i] = object;
    }
    return particles;
}

int main(int argc, char **argv) {
	float delta_time = 0.01;
      int nparticles = 128;
	float simulation_time = 10.0;
	if (argc > 3) {
		delta_time = atof(argv[1]);
		nparticles = atoi(argv[2]);
		simulation_time = atof(argv[3]);
	}
	double parallel_average_time = 0.0;
	MPI_Init(&argc, &argv);
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Datatype dt_psystem;
	MPI_Type_contiguous(7, MPI_FLOAT, &dt_psystem);
	MPI_Type_commit(&dt_psystem);

	if (rank == 0) {
		parallel_average_time -= MPI_Wtime();
	}
    	size_t items_per_process = nparticles / size;
    	psystem *particles = initparticles(nparticles);
    	if (rank == 0) {
    	    particles = initparticles(nparticles);
    	}
	psystem *local_particles = (psystem *) malloc(sizeof(*local_particles) * items_per_process);


    	MPI_Scatter(particles,items_per_process,dt_psystem,
        	local_particles,items_per_process,dt_psystem,0,
    	    	MPI_COMM_WORLD
    	);
	simulate(rank, nparticles, items_per_process, particles, local_particles, delta_time);

	psystem *gathered_particles = NULL;
	if(rank == 0){
		gathered_particles = (psystem*) malloc(sizeof(*gathered_particles) * nparticles);
	}
	MPI_Gather(local_particles,items_per_process,dt_psystem,
        	particles,items_per_process,dt_psystem,0,
        	MPI_COMM_WORLD
    	);
	if (rank == 0) {
        	parallel_average_time += MPI_Wtime();
          // printf("%f\n",parallel_average_time );
		printf("%d\n%.5f\n%.5f\n", nparticles, simulation_time, delta_time);
        	for (size_t i = 0; i < nparticles*2; ++i){
			printf("%.5f, %.5f,%.5f\n", particles[i].x, particles[i].y,particles[i].z);
			//printf("%.5f, %.5f, %.5f\n", particles[i].ax, particles[i].ay,particles[i].az );
			//printf("%.5f, %.5f\n, %.5f\n", particles[i].vx, particles[i].vy, particles[i].mass);
        	}
		for(float j = 0; j < simulation_time; j += delta_time){
			for(size_t i = 0; i < nparticles; ++i){
				//printf("%.5f, %.5f,%.5f\n", gathered_particles[i].ax, gathered_particles[i].ay, gathered_particles[i].az );
				integrate(&gathered_particles[i], delta_time);
			}
		}
    	}
	//if (particles != NULL) {
     //   	free(particles);
	//}
   // free(local_particles);
	//if(gathered_particles != NULL){
	//	free(gathered_particles);
	//}
 printf("%f\n",parallel_average_time );
    MPI_Finalize();
	return 0;
}
