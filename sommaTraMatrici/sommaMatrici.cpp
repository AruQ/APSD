#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include<cassert>
#include <string>
#include<iostream>

using namespace std;


const int nTypes = 6;
string types [nTypes] ={"parallel for", "collapse (dynamic 1000)", "collapse (static 1000)", "Linearizzato", "Linearizzato (static 1000)","seriale"};


bool checkResults (int ** a, int nRows, int nColumns);

void computeSerial (int ** a, int** c_serial, const int & nRows, const int & nColumns)
{
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nColumns; ++j) {
            c_serial[i][j] = a[i][j] - a[i][j];
        }
    }

}

void printTime (double& start, int i, int j, double ** time, int ** a, const int & nRows, const int & nColumns)
{

    double end = omp_get_wtime();
    if(checkResults(a, nRows, nColumns))
        time [i][j] = end-start;
    else
        time [i][j] = -1;
    cerr<< end-start << endl;
    start = omp_get_wtime();
}

int main(int argc, char *argv[])
{

    assert(argc>=4);
    int nRows = atoi(argv[1]);
    int nColumns = atoi(argv[2]);
    int nThreadsMin = atoi(argv[3]);
    int nThreadsMax = atoi(argv[4]);

    double* time[nTypes];
    for (int k = 0; k < nTypes; ++k) {
        time[k] = new double[nThreadsMax-nThreadsMin];
        for(int n=0; n<= nThreadsMax-nThreadsMin; n++ )
        {
            time[k][n] = 0.0;
        }
    }


    int** a = new int* [nRows];
    int** result = new int* [nRows];

    for(int i = 0; i<nRows; i++)
    {
        a[i] = new int [nColumns];
        result[i] = new int [nColumns];

    }

    for(int i = 0; i<nRows; i++)
    {
        for(int j = 0; j<nColumns; j++)
        {
            a[i][j] = i+j;
        }
    }


    double start = omp_get_wtime();

    computeSerial(a,result,nRows, nColumns);
    printTime(start, 5,0, time,result, nRows, nColumns );




    for (int n = nThreadsMin; n<=nThreadsMax; n++)
    {
        omp_set_num_threads(n);


#pragma omp parallel for

        for(int i=0; i<nRows; i++)
        {
            for(int j=0; j<nColumns; j++)
            {
                result[i][j] = a[i][j] - a[i][j];
            }
        }

        printTime(start, 0,n-nThreadsMin, time, result, nRows, nColumns);

#pragma omp parallel for schedule(dynamic,1000) collapse(2)
        for(int i=0; i<nRows; i++){
            for(int j=0; j<nColumns; j++) {
                result[i][j] = a[i][j] - a[i][j];
            }

        }

        printTime(start, 1,n-nThreadsMin, time,result, nRows, nColumns);

#pragma omp parallel for schedule(static,1000) collapse(2)
        for(int i=0; i<nRows; i++){
            for(int j=0; j<nColumns; j++) {
                result[i][j] = a[i][j] - a[i][j];
            }

        }

        printTime(start, 2,n-nThreadsMin, time, result, nRows, nColumns);


        int x, y;
#pragma omp parallel for private (x,y)
        for(int i=0; i<nRows*nColumns; i++){
            x = i/nColumns;
            y= i%nColumns;
            result [x][y] = a[x][y] - a[x][y];
        }

        printTime(start, 3,n-nThreadsMin, time,result, nRows, nColumns);


#pragma omp parallel for schedule(static,1000) private (x,y)
        for(int i=0; i<nRows*nColumns; i++){
            x = i/nColumns;
            y= i%nColumns;
            result [x][y] = a[x][y] - a[x][y];
        }
        printTime(start, 4,n-nThreadsMin, time, result, nRows, nColumns);

    }



    for (int k = 0; k < nTypes; ++k) {
        cout<<types[k]<<"; ";
        for (int n = 0; n <= nThreadsMax-nThreadsMin; ++n) {
            printf("%.4g;  ",time[k][n]);
        }
        printf ("\n");
    }


    for (int i = 0; i < nRows; ++i) {
        delete [] a[i];
        delete [] result[i];

    }

    delete [] a;
    delete [] result;

}

bool checkResults (int ** a, int nRows, int nColumns)
{
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nColumns; ++j) {
            if (a[i][j] != 0)
                return false;
        }

    }
    return true;

}




