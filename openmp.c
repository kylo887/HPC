#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main()
{
    int i,n=8;
#pragma omp parallel for
    
    for (i=0; i<n; i++)
    printf("Thread %d executes loop iteration %d\n", omp_get_thread_num(),i); 
    
}

/*#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
int main ()
{
int nthreads, tid;
/* Fork a team of threads giving them their own copies of variables */
//#pragma omp parallel private (tid, nthreads)
//{
/* Obtain thread number */
/*printf("%d\n",tid);
tid = omp_get_thread_num();

printf("Hello World from thread = %d\n", tid);
/* Only master thread does this */
/*if (tid == 0)
{
nthreads = omp_get_num_threads();
printf("Number of threads = %d\n", nthreads);
}
} /* All threads join master thread and disband */
//}
//int main()

//{
/*#pragma omp for schedule(static)
for(int n=0;n<10;n++)
    {
        int nn;
        nn=omp_get_thread_num();
        printf("%d\t%d\n", nn,n);
    }
    printf("dyn\n");
#pragma omp for schedule(dynamic)
for(int n=0;n<10;n++)
    {
        int nn;
        nn=omp_get_thread_num();
        printf("%d\t%d\n", nn,n);
  
    }
    printf("chunk\n");
int chunk =2;
#pragma omp for schedule(dynamic,chunk)
for(int n=0;n<10;n++)
    {
        int nn;
        nn=omp_get_thread_num();
        printf("%d\t%d\n", nn,n);
  
    }*/
// #ifdef _OPENMP
// omp_set_num_threads(2);
// #endif
// int i=1,j;
// printf("i=%d,j=%d\n",i,j);
// #pragma omp parallel for num_threads(3)  firstprivate(i) lastprivate(j)
// for(j=i;j<10;j++)
//     {
//         printf("%d\t",omp_get_thread_num());
//         printf("i=%d,j=%d\n",i,j);
//     }
//     printf("ii,i=%d,j=%d\n",i,j);
// }


