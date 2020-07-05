
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

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

void position()
typedef struct position,velocity,acceleration
{
  float x,y,z,vx,vy,vz,ax,ay,az;
}



