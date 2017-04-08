#include <iostream>
#include <omp.h>
#include <cassert>
#include <stdlib.h>


#define OMP 4

using namespace std;
//static long num_steps = 1000000;

double step;

int main (int argc, char *argv[])
{

    assert (argc >= 3 );

    long num_steps = atol(argv[1]);
    int nThreads = atoi(argv[2]);


    double start = omp_get_wtime();

    omp_set_num_threads(nThreads);

    int i;
    double x, pi;
    step = 1.0/(double) num_steps;

    //seriale
#if OMP==0
    double sum = 0.0;
    for (int i = 0; i < num_steps; ++i) {
        x= (i-0.5)*step;
        sum= sum + 4.0/(1.0+x*x);
    }
    pi= step * sum;

/*[parcuri@localhost piGreco seriale]$ g++ main.cpp -fopenmp && ./a.out 1000000000 1
Time: 4.377777163004794
PI: 3.141592655589971
*/
#elif OMP==1 //parallelo
    double sum [nThreads] = {0.0};
#pragma omp parallel
    {
        int i, id, nthrds;
        id = omp_get_thread_num ();
        nthrds = omp_get_num_threads();

        if (id ==0)
            nThreads= nthrds;
        for (i=id, sum[id]=0.0; i< num_steps; i=i+nthrds)
        {
            x=(i+0.5)*step;
            sum[id]+= 4.0 /(1.0+x*x);
            //static long num_steps = 1000000;
        }

    }

    for(i=0, pi=0.0; i<nThreads; i++)
    {
        pi+= sum[i]*step;
    }

#elif OMP==2
    double sum [nThreads] = {0.0};
#pragma omp parallel
    {
        int id= omp_get_thread_num();
#pragma omp for
        for (i=id;i< num_steps; i++){
            x = (i+0.5)*step;
            sum[id] += 4.0/(1.0+x*x);
        }
    }
    for(i=0, pi=0.0;i<nThreads;i++)
        pi += sum[i] * step;


    /*
[parcuri@localhost piGreco seriale]$ g++ main.cpp -fopenmp && ./a.out 1000000000 1
Time: 4.982200301004923
PI: 3.141592653589971
[parcuri@localhost piGreco seriale]$ g++ main.cpp -fopenmp && ./a.out 1000000000 2
Time: 16.06234886799939
PI: 3.141593816114751
[parcuri@localhost piGreco seriale]$ g++ main.cpp -fopenmp && ./a.out 1000000000 3
Time: 19.50799470100173
PI: 3.141365535760906
[parcuri@localhost piGreco seriale]$ g++ main.cpp -fopenmp && ./a.out 1000000000 4
Time: 19.52501809000387
PI: 3.141594577335071

     */

#elif OMP==3
    double sum;
#pragma omp parallel private (x, sum)
    {
        int id = omp_get_thread_num();
        for (i=id,sum=0.0;i< num_steps;i=i+nThreads){
            x = (i+0.5)*step;
            sum += 4.0/(1.0+x*x);
        }
#pragma omp critical
        pi += sum;
    }


    /*
[parcuri@localhost piGreco seriale]$ g++ main.cpp -fopenmp && ./a.out 1000000000 1
Time: 4.979461818998971
PI: 3141592653.589971
[parcuri@localhost piGreco seriale]$ g++ main.cpp -fopenmp && ./a.out 1000000000 2
Time: 6.803113484995265
PI: 3031481410.987821
[parcuri@localhost piGreco seriale]$ g++ main.cpp -fopenmp && ./a.out 1000000000 3
Time: 6.79366057599691
PI: 3302926718.182592
[parcuri@localhost piGreco seriale]$ g++ main.cpp -fopenmp && ./a.out 1000000000 4
Time: 7.625502839000546
PI: 3031995515.266843
*/

#elif OMP==4
    long long numIn = 0;
    double y;
   #pragma omp parallel private(x, y) reduction(+:numIn)
   {
       #pragma omp for
       for (i = 0; i <= num_steps; i++) {
           x = (double)rand()/RAND_MAX;
           y = (double)rand()/RAND_MAX;

           if (x*x + y*y <= 1)
               numIn++;
       }
   }
   pi = 4.0 *((double)numIn/(double)num_steps);
   printf("pi %f\n", pi);



#elif OMP==5
    double sum=0.0;
    #pragma omp parallel for private(x) reduction(+:sum)
    for (i = 0; i < num_steps; i++)
    {
      x = (i + 0.5)* step;
      sum += 4.0/(1.0 + x*x);
    }
    pi = sum / num_steps;

//    [parcuri@localhost piGreco seriale]$ g++ main.cpp -fopenmp && ./a.out 1000000000 4
//    Time: 1.210655486000178
//    PI: 3.141592653589821
//    [parcuri@localhost piGreco seriale]$ g++ main.cpp -fopenmp && ./a.out 1000000000 3
//    Time: 1.542759682999531
//    PI: 3.141592653589961
//    [parcuri@localhost piGreco seriale]$ g++ main.cpp -fopenmp && ./a.out 1000000000 2
//    Time: 2.201832018999994
//    PI: 3.141592653589901
//    [parcuri@localhost piGreco seriale]$ g++ main.cpp -fopenmp && ./a.out 1000000000 1
//    Time: 4.231976058001237
//    PI: 3.141592653589971

#endif





    double end = omp_get_wtime();
    printf("Time: %.16g\n", end-start);

    printf("PI: %.16g \n", pi);
}
