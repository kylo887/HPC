#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>


int particles;                                          //initializing no of particles//
float ttime,delt;                                       //initializing total time and time step//
int default_particles=20,boundary=20;                   //default values//
float default_ttime=10,default_delt=1,default_mass=10;  //default values//

struct position{                                        //structure defined for position//
    float x_coord,y_coord,z_coord;                      //coordinates//
    float x_coord_new,y_coord_new,z_coord_new;
}*po;                                                   //pointer//

struct velocity{                                        //structure defined for velocity//
    float v_x_coord,v_y_coord,v_z_coord;                //velocity values//
}*ve;

struct acceleration{                                    //structure defined for acceleration//
    float a_x_coord,a_y_coord,a_z_coord;                //acc values//
}*ac;

int *c,*d,e,f,g ;                                       //initialize variables to divide the number//
                                                        //of particles to processors//

void fn(){                                              //function to assign particles//                                           
int h=g%e,i=0;
for (int j=0;j<e;j++)
{
         c[j] = (g/e); 
         if (h>0)
         {
            c[j]+=1;
            h--;
         }
         d[j] =i; 
         i+=c[j];
}}


void fn_2(){                                           //function to randomize starting positions//
                                                       //and initial velocities// 
    for(int k=0;k<particles;k++){
        po[k].x_coord=((double)rand()/RAND_MAX*1.0-0.0)*boundary; 
        po[k].y_coord=((double)rand()/RAND_MAX*1.0-0.0)*boundary;
        po[k].z_coord=((double)rand()/RAND_MAX*1.0-0.0)*boundary;
        ve[k].v_x_coord=(double)rand()/RAND_MAX*0.0001-0.0;
        ve[k].v_y_coord=(double)rand()/RAND_MAX*0.0001-0.0;
        ve[k].v_z_coord=(double)rand()/RAND_MAX*0.0001-0.0;}}


void fn_3(){                                          //function to write coordinates to a text file//
    FILE *l = fopen("ll.txt","a");                    //appends//
    for(int m=0; m<particles;m++){
        fprintf(l,"%f %f %f\n",po[m].x_coord_new,po[m].y_coord_new,po[m].z_coord_new);}
    fprintf(l,"\n");
    fclose(l);}


void fn_4(){                                            //function where the simulation happens//
    MPI_Bcast(po,particles,MPI_FLOAT ,0,MPI_COMM_WORLD);//broacasts the structure values to all processors// 
    MPI_Bcast(ve,particles,MPI_FLOAT,0,MPI_COMM_WORLD);
    int n=ttime/delt;
    for(int o=0;o<n;o++){
        float d_x_coord,d_y_coord,d_z_coord,inv,inv_3,forc,x_coord_new,y_coord_new,z_coord_new; //acceleration//
        for(int p=0;p<particles;p++){
            ac[p].a_x_coord=0.0;
            ac[p].a_y_coord=0.0;
            ac[p].a_z_coord=0.0;
            for(int q=0;q<particles;q++){
                if(p+d[f]!=q){
                    d_x_coord=po[p+d[f]].x_coord-po[q].x_coord;
                    d_y_coord=po[p+d[f]].y_coord-po[q].y_coord;
                    d_z_coord=po[p+d[f]].z_coord-po[q].z_coord;
                    inv=1.0/sqrt((d_x_coord*d_x_coord)+(d_y_coord*d_y_coord)+(d_z_coord*d_z_coord)+__DBL_EPSILON__);
                    inv_3=pow(inv,3);
                    forc=default_mass*inv_3;
                    ac[p].a_x_coord+=forc*d_x_coord;
                    ac[p].a_y_coord+=forc*d_y_coord;
                    ac[p].a_z_coord+=forc*d_z_coord;
            
            po[p+(d[f])].x_coord_new=(delt*ve[p+d[f]].v_x_coord)+(0.5*delt*delt*ac[p].a_x_coord); //updating position from vel//
            po[p+(d[f])].y_coord_new=(delt*ve[p+d[f]].v_y_coord)+(0.5*delt*delt*ac[p].a_x_coord); //and acceleration          //
            po[p+(d[f])].z_coord_new=(delt*ve[p+d[f]].v_z_coord)+(0.5*delt*delt*ac[p].a_x_coord);
            
            ve[p+d[f]].v_x_coord+=abs((delt*ac[p].a_x_coord));                                    //function to calculate resultant//
            ve[p+d[f]].v_y_coord+=abs((delt*ac[p].a_y_coord));                                    //velocity//  
            ve[p+d[f]].v_z_coord+=abs((delt*ac[p].a_z_coord));}}}
                
         MPI_Allgatherv( po+d[f],c[f],MPI_FLOAT,po,c,d,MPI_FLOAT, MPI_COMM_WORLD);                
         MPI_Allgatherv( ve+d[f],c[f],MPI_FLOAT,ve,c,d,MPI_FLOAT, MPI_COMM_WORLD);
            if(f==0){
                fn_3();}}}



int main(int argc,char **argv){                               //main function takes arguments//
    srand(time(0)); 
    MPI_Init(&argc,&argv);                                    //starting MPI//  
    MPI_Barrier(MPI_COMM_WORLD);
    double tt=MPI_Wtime();                                    //time starts//  
    if (argc>=2)                                              //assigning values acc to indexes//              
        particles=atoi(argv[1]);
    else
        particles=default_particles;
    if (argc>=3)
        ttime=atoi(argv[2]);
    else
        ttime=default_ttime;   
    if  (argc>=4)
        delt=atof(argv[3]);
    else
        delt=default_delt;
    MPI_Comm_size(MPI_COMM_WORLD,&e);                       //size and rank//
    MPI_Comm_rank(MPI_COMM_WORLD,&f);
   int array_0[3]={1,1,1};                                  //arrays initialized for mpi create struct//
   MPI_Aint const array_1[3]={0,sizeof(float),sizeof(float)*2};
   MPI_Datatype array_2[3]={MPI_FLOAT,MPI_FLOAT,MPI_FLOAT};
   MPI_Datatype P;
   MPI_Datatype V;
   MPI_Type_create_struct(3,array_0,array_1,array_2,&P);
   MPI_Type_commit(&P);
   MPI_Type_create_struct(3,array_0,array_1,array_2,&V);
   MPI_Type_commit(&V);
   c=malloc(sizeof(int)*e);                                //memory allocation for variables//
   d=malloc(sizeof(int)*e);
   fn();                                                   //calling first fn to divide// 
   po=malloc(particles*sizeof(po));
   ve=malloc(particles*sizeof(ve));
   ac=malloc(c[f]*sizeof(ac));
   if(f==0){                                              //writing particle coordinates to the file//
       fn_2();                                            //second function called to initialize position// 
       FILE *l=fopen("ll.txt","w");                       //and velocity//
       for(int t=0;t<particles;t++){
           fprintf(l,"%f %f %f \n",po[t].x_coord,po[t].y_coord,po[t].z_coord);}
        fprintf(l,"\n");
        fclose(l);}
    fn_4();                                               //perform simulation//
    MPI_Barrier( MPI_COMM_WORLD);                         //to sync//
    double ttt=MPI_Wtime()-tt;                            //time ends//
    if(f==0){                                             //final output//
        printf("particles=%i,time=%f,delt=%f\n",particles,ttime,delt);
        printf("time=%f\n",ttt);}

 MPI_Finalize();}                                         //MPI ends//