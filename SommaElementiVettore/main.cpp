#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include<cassert>

#include<iostream>
#include<functional>
#include "Time.h"

using namespace std;

const int nTypes = 7;
string types [nTypes] ={"Serial", "Atomic","Critical", "Reduction", "Sum[N_THREADS]", "Dynamic", "Static"};


//0 atomic | 1 critical | 2 reduction | 3 altra | 4 dymamic | 5 static |

auto checkResult = [](int a, int b)->bool  {
    return a == b;
};

using namespace std;

int main(int argc, char *argv[])
{

    assert(argc>=5);
    int dimension = atoi(argv[1]);
    int chunkSize = atoi(argv[2]);
    int nThreadsMin = atoi(argv[3]);
    int nThreadsMax = atoi(argv[4]);


    Time<nTypes> time (types,nThreadsMax-nThreadsMin+1);

    int* a = new int [dimension];
    for (int n = 0; n < dimension; ++n) {
        a[n]=n;

    }

    long long serial_sum = 0;

    int i;
    /* Serial Calculation */
    for (i = 0; i < dimension; i++)
        serial_sum = serial_sum + a[i];

    time.stop(0,0,std::bind (checkResult,serial_sum, serial_sum));
    int indexTypes;


    //*****************************************************************************//
    for (int n = nThreadsMin; n<=nThreadsMax; n++)
    {

        long long sum=0;
        indexTypes = 1;
        omp_set_num_threads(n);


        time.start();


#pragma omp parallel for shared(sum)
        for (i = 0; i < dimension; i++) {
#pragma omp atomic
            sum = sum + a[i];

        }

        time.stop(indexTypes++,n-nThreadsMin,std::bind (checkResult,sum, serial_sum));
        cout<<types[indexTypes-1]<<" - ";
        printf ("SUM: %llu\n", sum);




        //    [parcuri@localhost SommaElementiVettore]$ g++ -fopenmp main.cpp && ./a.out 10000000000 4
        //    Parallel time execution  51.87878787299996
        //    Serial time execution  3.772340397000335

        //    The Serial And Parallel Sums Are Equal

        //*****************************************************************************//
        time.start();

        sum=0;
#pragma omp parallel for shared(sum)
        for (i = 0; i < dimension; i++) {
#pragma omp critical
            sum = sum + a[i];

        }

        //        [parcuri@localhost SommaElementiVettore]$ g++ -fopenmp main.cpp && ./a.out 10000000000 4

        //        Parallel time execution  192.433580289
        //        Serial time execution  3.759798951999983

        //        The Serial And Parallel Sums Are Equal

        time.stop(indexTypes++,n-nThreadsMin,std::bind (checkResult,sum, serial_sum));
        cout<<types[indexTypes-1]<<" - ";
        printf ("SUM: %llu\n", sum);

        //*****************************************************************************//
        time.start();

        sum=0;

#pragma omp parallel for reduction (+:sum)
        for (i=0;i<dimension;i++)
            sum=sum+a[i];

        time.stop(indexTypes++,n-nThreadsMin,std::bind (checkResult,sum, serial_sum));
        cout<<types[indexTypes-1]<<" - ";
        printf ("SUM: %llu\n", sum);
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

        //*****************************************************************************//
        time.start();


        int nThreads = n;
        sum = 0.0;

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


        time.stop(indexTypes++,n-nThreadsMin,std::bind (checkResult,sum, serial_sum));
        cout<<types[indexTypes-1]<<" - ";
        printf ("SUM: %llu\n", sum);

        //*****************************************************************************//
        time.start();
        sum=0;

#pragma omp parallel for schedule (dynamic, chunkSize)

        for(i=0; i<dimension; i++)
        {
#pragma omp critical
            sum = sum + a[i];
        }

        time.stop(indexTypes++,n-nThreadsMin,std::bind (checkResult,sum, serial_sum));
        cout<<types[indexTypes-1]<<" - ";
        printf ("SUM: %llu\n", sum);

        //*****************************************************************************//
        time.start();

        sum=0;
#pragma omp parallel for schedule (static, chunkSize)

        for(i=0; i<dimension; i++)
        {
#pragma omp critical
            sum = sum + a[i];
        }


        time.stop(indexTypes++,n-nThreadsMin,std::bind (checkResult,sum, serial_sum));
        cout<<types[indexTypes-1]<<" - ";
        printf ("SUM: %llu\n", sum);

    }
    time.print();
    delete [] a;
}


