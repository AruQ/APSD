#ifndef TIME_H
#define TIME_H

#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <omp.h>
using namespace std;


template<int DIMENSION>
class Time
{

public:
    Time (std::string* types, int nThreads) : types(types), nThreads(nThreads),partial_time(0)
    {
        for (int k = 0; k < DIMENSION; ++k) {
            time[k] = new double[nThreads];

            for(int i=0; i< nThreads; i++)
            {
                time[k][i] = -1;
            }
        }


    }
    ~Time ()
    {
        for (int n = 0; n <DIMENSION; ++n) {
            delete [] time [n];
        }

    }

    void start ()
    {
        this->start_time = omp_get_wtime();
    }

	void wait()
	{
		end_time = omp_get_wtime();
		partial_time += end_time-start_time;
	}
	
	void continue_()
	{
		start();
	}

    void stop (int i, int j, auto checkResult)
    {

        end_time = omp_get_wtime();
        if(checkResult())
            time [i][j] = end_time-start_time + partial_time;
        else
            time [i][j] = -1;
		partial_time = 0;
    }


    void print ()
    {
        for(int k=0; k< DIMENSION; k++)
        {
            cout<<types[k]<< "; ";
            for (int n=0; n< nThreads; n++)

                printf("%.4g;  ",time[k][n]);
            printf ("\n");
        }
    }



private:
    std::string* types;
    double* time[DIMENSION];
    double start_time;
    double end_time;
    int nThreads;
	double partial_time;
};


#endif
