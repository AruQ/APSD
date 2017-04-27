#include <cassert>
#include <functional>

#include "Time.h"
#define PAD 64

//PI GRECO
const int nTypes = 10;
string types [nTypes] ={"Seriale","Simple Parallel (SUM[N_THREADS])", "Parallel for", "Dynamic", "Static", "Critical","Atomic","Padding", "Reduction", "Montecarlo"};

using namespace std;

auto checkResult = [](double PI)->bool  {
    return PI<3.149 && PI >3.141;
};

int main (int argc, char *argv[])
{

    assert (argc >= 4 );

    long num_steps = atol(argv[1]);
    int nThreadsMin = atoi(argv[2]);
    int nThreadsMax = atoi(argv[3]);

    Time<nTypes> time (types,nThreadsMax-nThreadsMin+1);

    double step;

    int i;
    double x, pi;
    step = 1.0/(double) num_steps;


    time.start();
    //*****************************************************************************//
    //SERIALE

    double sum = 0.0;
    for (int i = 0; i < num_steps; ++i) {
        x= (i-0.5)*step;
        sum= sum + 4.0/(1.0+x*x);
    }
    pi= step * sum;

    time.stop(0,0,std::bind (checkResult,pi));
    cout<<types[0]<<" - ";
    printf("PI: %.16g \n", pi);


    int indexTypes;
    //*****************************************************************************//
    for (int n = nThreadsMin; n<=nThreadsMax; n++)
    {

        indexTypes = 1;
        omp_set_num_threads(n);


        time.start();
        //*****************************************************************************//
        //SIMPLE PARALLEL
        int nThreads=n;
        double sum_array [nThreads] = {0.0};
#pragma omp parallel
        {
            int i, id, nthrds;
            id = omp_get_thread_num ();
            nthrds = omp_get_num_threads();

            if (id ==0)
                nThreads= nthrds;
            for (i=id, sum_array[id]=0.0; i< num_steps; i=i+nthrds)
            {
                x=(i+0.5)*step;
                sum_array[id]+= 4.0 /(1.0+x*x);
            }

        }

        for(i=0, pi=0.0; i<nThreads; i++)
        {
            pi+= sum_array[i]*step;
        }
        time.stop(indexTypes++,n-nThreadsMin,std::bind (checkResult,pi));
        cout<<types[indexTypes-1]<<" - ";
        printf("PI: %.16g \n", pi);
        time.start();

        //*****************************************************************************//
        //PARALLEL FOR

        for(int k=0; k<nThreads; k++)
        {
            sum_array[k] = 0.0;
        }
#pragma omp parallel
        {
            int id= omp_get_thread_num();
#pragma omp for
            for (i=id;i< num_steps; i++){
                x = (i+0.5)*step;
                sum_array[id] += 4.0/(1.0+x*x);
            }
        }
        for(i=0, pi=0.0;i<nThreads;i++)
            pi += sum_array[i] * step;

        time.stop(indexTypes++,n-nThreadsMin,std::bind (checkResult,pi));
        cout<<types[indexTypes-1]<<" - ";
        printf("PI: %.16g \n", pi);
        time.start();



        //*****************************************************************************//
        //DYNAMIC

        for(int k=0; k<nThreads; k++)
        {
            sum_array[k] = 0.0;
        }
#pragma omp parallel
        {
            int id= omp_get_thread_num();
#pragma omp for schedule (dynamic, 1000)
            for (i=id;i< num_steps; i++){
                x = (i+0.5)*step;
                sum_array[id] += 4.0/(1.0+x*x);
            }
        }
        for(i=0, pi=0.0;i<nThreads;i++)
            pi += sum_array[i] * step;

        time.stop(indexTypes++,n-nThreadsMin,std::bind (checkResult,pi));
        cout<<types[indexTypes-1]<<" - ";
        printf("PI: %.16g \n", pi);
        time.start();



        //*****************************************************************************//
        //STATIC

        for(int k=0; k<nThreads; k++)
        {
            sum_array[k] = 0.0;
        }
#pragma omp parallel
        {
            int id= omp_get_thread_num();
#pragma omp for schedule (static, 1000)
            for (i=id;i< num_steps; i++){
                x = (i+0.5)*step;
                sum_array[id] += 4.0/(1.0+x*x);
            }
        }
        for(i=0, pi=0.0;i<nThreads;i++)
            pi += sum_array[i] * step;

        time.stop(indexTypes++,n-nThreadsMin,std::bind (checkResult,pi));
        cout<<types[indexTypes-1]<<" - ";
        printf("PI: %.16g \n", pi);
        time.start();



        //*****************************************************************************//
        //CRITICAL

        nThreads= n;
#pragma omp parallel
        {
            int i, id,nthrds; double x, sum;
            id = omp_get_thread_num();
            nthrds = omp_get_num_threads();
            if (id == 0) nThreads = nthrds;
            id = omp_get_thread_num();
            nthrds = omp_get_num_threads();
            for (i=id, sum=0.0;i< num_steps; i=i+nThreads){
                x = (i+0.5)*step;
                sum += 4.0/(1.0+x*x);
            }
#pragma omp critical
            pi += sum * step;
        }

        time.stop(indexTypes++,n-nThreadsMin,std::bind (checkResult,pi));
        cout<<types[indexTypes-1]<<" - ";
        printf("PI: %.16g \n", pi);
        time.start();


        //*****************************************************************************//
        //ATOMIC

        pi=0.0;
        nThreads= n;
#pragma omp parallel
        {
            int i, id,nthrds; double x, sum;
            id = omp_get_thread_num();
            nthrds = omp_get_num_threads();
            if (id == 0) nThreads = nthrds;
            id = omp_get_thread_num();
            nthrds = omp_get_num_threads();
            for (i=id, sum=0.0;i< num_steps; i=i+nThreads){
                x = (i+0.5)*step;
                sum += 4.0/(1.0+x*x);
            }
#pragma omp atomic
            pi += sum * step;
        }

        time.stop(indexTypes++,n-nThreadsMin,std::bind (checkResult,pi));
        cout<<types[indexTypes-1]<<" - ";
        printf("PI: %.16g \n", pi);
        time.start();

        //*****************************************************************************//
        //PADDING

       double sum_padding [nThreads][PAD] = {0.0};
#pragma omp parallel
        {
            int id= omp_get_thread_num();
#pragma omp for
            for (i=id;i< num_steps; i++){
                x = (i+0.5)*step;
                sum_padding[id][0] += 4.0/(1.0+x*x);
            }
        }
        for(i=0, pi=0.0;i<nThreads;i++)
            pi += sum_padding[i][0] * step;

        time.stop(indexTypes++,n-nThreadsMin,std::bind (checkResult,pi));
        cout<<types[indexTypes-1]<<" - ";
        printf("PI: %.16g \n", pi);
        time.start();

        //*****************************************************************************//
        //REDUCTION
        sum=0.0;
#pragma omp parallel for private(x) reduction(+:sum)
        for (i = 0; i < num_steps; i++)
        {
            x = (i + 0.5)* step;
            sum += 4.0/(1.0 + x*x);
        }
        pi = sum / num_steps;

        time.stop(indexTypes++,n-nThreadsMin,std::bind (checkResult,pi));
        cout<<types[indexTypes-1]<<" - ";
        printf("PI: %.16g \n", pi);
        time.start();

        //*****************************************************************************//
        //MONTECARLO

        long long numIn = 0;
        double y;
#pragma omp parallel private(x, y) reduction(+:numIn)
        {
            unsigned int id = omp_get_thread_num();

#pragma omp for
            for (i = 0; i <= num_steps; i++) {
                x = (double)rand_r(&id)/RAND_MAX;
                y = (double)rand_r(&id)/RAND_MAX;

                if (x*x + y*y <= 1)
                    numIn++;
            }
        }
        pi = 4.0 *((double)numIn/(double)num_steps);

        time.stop(indexTypes++,n-nThreadsMin,std::bind (checkResult,pi));
        cout<<types[indexTypes-1]<<" - ";
        printf("PI: %.16g \n", pi);

        cout<<std::endl;
    }


    time.print();


}
