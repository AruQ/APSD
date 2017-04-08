#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include<cassert>

using namespace std;
#define TYPE 0

bool checkResults (int ** a, int **b,int nRows, int nColumns);

int main(int argc, char *argv[])
{

    assert(argc>=3);
    int nRows = atoi(argv[1]);
    int nColumns = atoi(argv[2]);
    int nThreads = atoi(argv[3]);


    int** a = new int* [nRows];
    int** b = new int* [nRows];
    int** c_parallel = new int* [nRows];
    int** c_serial = new int* [nRows];
    for(int i = 0; i<nRows; i++)
    {
        a[i] = new int [nColumns];
        b[i] = new int [nColumns];
        c_parallel[i] = new int [nColumns];
        c_serial[i] = new int [nColumns];
    }

    for(int i = 0; i<nRows; i++)
    {
        for(int j = 0; j<nColumns; j++)
        {
            a[i][j] = i+j;
            b[i][j]= i-j;
        }
    }


    double start = omp_get_wtime();
    omp_set_num_threads(nThreads);

    printf("sono qui\n");

#if TYPE==0
#pragma omp parallel for

    for(int i=0; i<nRows; i++)
    {
        for(int j=0; j<nColumns; j++)
        {
            c_parallel[i][j] = a[i][j] + b[i][j];
        }
    }

    printf("sono qui\n");

#endif
#if TYPE==1
#pragma omp parallel for schedule(dynamic,1000) collapse(2)
    for(int i=0; i<nRows; i++){
        for(int j=0; j<nColumns; j++) {
            c_parallel[i][j] = a[i][j] + b[i][j];
        }

    }


#endif

#if TYPE==2
#pragma omp parallel for schedule(static,1000) collapse(2)
    for(int i=0; i<nRows; i++){
        for(int j=0; j<nColumns; j++) {
            c_parallel[i][j] = a[i][j] + b[i][j];
        }

    }


#endif

    int x, y;
#if TYPE==3
#pragma omp parallel for private (x,y)
    for(int i=0; i<nRows*nColumns; i++){
        x = i/nColumns;
        y= i%nColumns;
        c_parallel [x][y] = a[x][y] + b[x][y];


    }


#endif

#if TYPE==4
#pragma omp parallel for schedule(static,1000) private (x,y)
    for(int i=0; i<nRows*nColumns; i++){
        x = i/nColumns;
        y= i%nColumns;
        c_parallel [x][y] = a[x][y] + b[x][y];


    }


#endif

    double end = omp_get_wtime();
    printf("Parallel \nTIME: %.16g\n", end-start);



    start = omp_get_wtime();

    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nColumns; ++j) {
            c_serial[i][j] = a[i][j] +b [i][j];
        }
    }

    end = omp_get_wtime();
    printf("Serial \nTIME: %.16g\n", end-start);
    if (checkResults(c_parallel, c_serial, nRows, nColumns))

        printf("il risultato è corretto \n");

    else
        printf("il risultato non è corretto \n");


    for (int i = 0; i < nRows; ++i) {
        delete a[i];
        delete b[i];
        delete c_serial[i];
        delete c_parallel[i];

    }

    delete [] a;
    delete [] b;
    delete [] c_serial;
    delete [] c_parallel;




}

bool checkResults (int ** a, int **b, int nRows, int nColumns)
{
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nColumns; ++j) {
            if (a[i][j] != b[i][j])
                return false;
        }

    }
    return true;

}




