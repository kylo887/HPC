#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#define MATH_PI 3.1415

typedef struct _psystem {
	float x, y, z;	// position
	float ax, ay, az;	// acceleration
	float vx, vy, vz;	// velocity
	float mass;	// mass
} psystem;

void integrate(psystem *psystem, float deltaTime) {
	  psystem->vx += psystem->ax * deltaTime;
        psystem->vy += psystem->ay * deltaTime;
        psystem->vz += psystem->az * deltaTime;
        psystem->x += psystem->vx * deltaTime + (1 / 2) * psystem->ax * deltaTime * deltaTime;
        psystem->y += psystem->vy * deltaTime + (1 / 2) * psystem->ay * deltaTime * deltaTime;
        psystem->z += psystem->vz * deltaTime + (1 / 2) * psystem->az * deltaTime * deltaTime;
//printf("%f\t %f\t %f\n", psystem->x,psystem->y,psystem->z );
//printf("%f\t %f\t %f\n", psystem->vx,psystem->vy,psystem->vz );
//output positions are saved in a textfile
FILE *output;
output= fopen("outputpositions.txt","a");
//fprintf(output,"\n");
fprintf(output,"%f\t %f\t %f\n", psystem->x,psystem->y,psystem->z );
//printf("%f\t %f\t %f\n", psystem->x,psystem->y,psystem->z );
fclose(output);
        
}

void acceleration(psystem *a, psystem *b, float *ax, float *ay , float *az) {
        float distanceX = fabsf(b->x - a->x);
        float distanceY = fabsf(b->y - a->y);
        float distanceZ = fabsf(b->z - a->z);
        float vectorDistance = sqrt(a->x * a->x + a->y * a->y +a->z * a->z  );
        float vectorDistanceCubed = vectorDistance * vectorDistance * vectorDistance;
        float inverse = 1.0 / vectorDistanceCubed;
        float scale = b->mass * inverse;
        *ax = (distanceX * scale);
        *ay = (distanceY * scale);
        *az = (distanceZ * scale);

}

void simulate(int h,int nparticles, psystem *particles, float dt) {
        double total_time = 0.0;
	for(size_t i = 0; i < nparticles; i++) {
		double start = clock();
                float total_ax = 0, total_ay = 0, total_az = 0;
                for (size_t j = 0; j < h; j++) {
                       // if (i == j) {
                            //    continue;
                        //}
                        float ax, ay, az;
//printf("%f\t %f\t %f\n", ax,ay,az );
                       acceleration(&particles[i], &particles[j], &ax, &ay, &az);
                        total_ax += ax;
                        total_ay += ay;
                        total_az += az;
                }
                particles[i].ax = total_ax;
                particles[i].ay = total_ay;
                particles[i].az = total_az;
printf("%f\n", total_ax );


                integrate(&particles[i], dt);
		double time_elapsed = ((double) clock() - start) / CLOCKS_PER_SEC;
		total_time += time_elapsed;
        }
	printf("%f\n", total_time);
}

void initparticles (int h,int nparticles, psystem *particles, float delta_time) {
        srand(time(NULL));
        
        for (int i = 0; i < nparticles; i++) {
        for(int j = 0; j < h; j++) {
                float initialMass = 2;
                psystem object = {
                        .x = ((double) rand() / (RAND_MAX))*2 - 1.15, 
                        .y = ((double) rand() / (RAND_MAX))*2 - 1.15,
                        .z = ((double) rand() / (RAND_MAX))*2 - 1.15,
                        .vx = (((double) rand() - 0.5) * 0.5)*1.2,
                        .vy = (((double) rand() - 0.5) * 0.5)*2.1,
                        .vz = (((double) rand() - 0.5) * 0.5)*0.8,
                        .mass = (double) rand() + initialMass * 0.5

                };

        particles[i] = object;


    }
//integrate(&particles[i], delta_time);
        
}
}
int main(int argc, char **argv) {
	float delta_time = 0.5;
     double duration=10;
     int h= duration/ delta_time;
	int nparticles = 100;
	//if (argc > 2) {
	//	delta_time = atof(argv[1]);
	//	nparticles = atoi(argv[2]);
	//}
	//printf("%d\n", argc);
	psystem *particles = (psystem*) malloc(nparticles * h*sizeof(*particles));
	initparticles(h,nparticles, particles,delta_time);
	simulate(h,nparticles, particles, delta_time);
	return 0;
}
