#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include<cassert>

#include<iostream>

//0 atomic | 1 critical | 2 reduction | 3 altra | 4 dymamic | 5 static |
#define TYPE 3

using namespace std;

int main(int argc, char *argv[])
{

    assert(argc>=4);
    int dimension = atoi(argv[1]);
    int nThreads = atoi(argv[2]);
    int chunkSize = atoi(argv[3]);


    int* a = new int [dimension];
    for (int n = 0; n < dimension; ++n) {
        a[n]=n;

    }


    printf("inizia l'esecuzione...\n");

    double start = omp_get_wtime();
    int i;
    omp_set_num_threads(nThreads);

    long long sum=0;

#if TYPE==0
#pragma omp parallel for shared(sum)
    for (i = 0; i < dimension; i++) {
#pragma omp atomic
        sum = sum + a[i];

    }


    //    [parcuri@localhost SommaElementiVettore]$ g++ -fopenmp main.cpp && ./a.out 10000000000 4
    //    Parallel time execution  51.87878787299996
    //    Serial time execution  3.772340397000335

    //    The Serial And Parallel Sums Are Equal


#endif

#if TYPE==1
#pragma omp parallel for shared(sum)
    for (i = 0; i < dimension; i++) {
#pragma omp critical
        sum = sum + a[i];

    }

    //        [parcuri@localhost SommaElementiVettore]$ g++ -fopenmp main.cpp && ./a.out 10000000000 4

    //        Parallel time execution  192.433580289
    //        Serial time execution  3.759798951999983

    //        The Serial And Parallel Sums Are Equal

#endif



#if TYPE==2

#pragma omp parallel for reduction (+:sum)
    for (i=0;i<dimension;i++)
        sum=sum+a[i];


    //    [parcuri@localhost SommaElementiVettore]$ g++ -fopenmp main.cpp && ./a.out 10000000000 4
    //    Parallel time execution  1.040887769000165
    //    Serial time execution  3.762261832000149

    //    The Serial And Parallel Sums Are Equal



    //    [parcuri@localhost SommaElementiVettore]$ g++ -fopenmp main.cpp && ./a.out 10000000000 3
    //    Parallel time execution  1.338144080999882
    //    Serial time execution  3.69525283899975

    //    The Serial And Parallel Sums Are Equal

    //    [parcuri@localhost SommaElementiVettore]$ g++ -fopenmp main.cpp && ./a.out 10000000000 2
    //    Parallel time execution  1.973348238000199
    //    Serial time execution  3.781234048000442

    //    The Serial And Parallel Sums Are Equal


    //    [parcuri@localhost SommaElementiVettore]$ g++ -fopenmp main.cpp && ./a.out 10000000000 1
    //    Parallel time execution  3.843764377999833
    //    Serial time execution  3.680074147000141

    //    The Serial And Parallel Sums Are Equal



#endif


#if TYPE==3

    long long sumN [nThreads] = {0};
#pragma omp parallel shared (sumN)
    {
        int id = omp_get_thread_num ();
        int nT = omp_get_num_threads();


#pragma omp for
        for (i = 0; i < dimension; i++) {
            sumN[id] = sumN[id]+ a[i];

        }
    }


    for (i = 0; i < nThreads; ++i) {
        sum += sumN[i];
    }
    printf ("%llu\n", sum);
#endif


#if TYPE==4
#pragma omp parallel for schedule (dynamic, chunkSize)

    for(i=0; i<dimension; i++)
    {
#pragma omp critical
        sum = sum + a[i];
    }

#endif

#if TYPE==5
#pragma omp parallel for schedule (static, chunkSize)

    for(i=0; i<dimension; i++)
    {
#pragma omp critical
        sum = sum + a[i];
    }

#endif


    double end = omp_get_wtime();
    printf("Parallel time execution  %.16g\n", end-start);

    start = omp_get_wtime();

    long long serial_sum = 0;

    /* Serial Calculation */
    for (i = 0; i < dimension; i++)
        serial_sum = serial_sum + a[i];


    end = omp_get_wtime();
    printf("Serial time execution  %.16g\n", end-start);

    if (serial_sum == sum)
        printf("\nThe Serial And Parallel Sums Are Equal %llu\n", sum);
    else {
        printf("\nThe Serial And Parallel Sums Are UnEqual\n");
    }




    delete [] a;




}



