#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include<cassert>
#include <string>
#include<iostream>

using namespace std;
#define TYPE 0


const int nTypes = 6;
string types [nTypes] ={"parallel for", "collapse (dynamic 1000)", "collapse (static 1000)", "Linearizzato", "Linearizzato (static 1000)","seriale"};


bool checkResults (int ** a, int **b,int nRows, int nColumns);

void computeSerial (int ** a, int **b, int** c_serial, const int & nRows, const int & nColumns)
{
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nColumns; ++j) {
            c_serial[i][j] = a[i][j] +b [i][j];
        }
    }

}

void printTime (double& start, int i, int j, double ** time, int ** a, int **b,const int & nRows, const int & nColumns)
{

    double end = omp_get_wtime();
    if(checkResults(a,b,nRows, nColumns))
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
    }


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

    computeSerial(a,b,c_serial,nRows, nColumns);
    printTime(start, 5,0, time,c_serial,c_serial, nRows, nColumns );




    for (int n = nThreadsMin; n<=nThreadsMax; n++)
    {
        omp_set_num_threads(n);

//        printf("Numero threads: %d\n", n);



        //#if TYPE==0
#pragma omp parallel for

        for(int i=0; i<nRows; i++)
        {
            for(int j=0; j<nColumns; j++)
            {
                c_parallel[i][j] = a[i][j] + b[i][j];
            }
        }




        printTime(start, 0,n-nThreadsMin, time, c_parallel, c_serial, nRows, nColumns);
        //#endif
        //#if TYPE==1
#pragma omp parallel for schedule(dynamic,1000) collapse(2)
        for(int i=0; i<nRows; i++){
            for(int j=0; j<nColumns; j++) {
                c_parallel[i][j] = a[i][j] + b[i][j];
            }

        }

        printTime(start, 1,n-nThreadsMin, time,c_parallel, c_serial, nRows, nColumns);
        //#endif

        //#if TYPE==2
#pragma omp parallel for schedule(static,1000) collapse(2)
        for(int i=0; i<nRows; i++){
            for(int j=0; j<nColumns; j++) {
                c_parallel[i][j] = a[i][j] + b[i][j];
            }

        }

        printTime(start, 2,n-nThreadsMin, time, c_parallel, c_serial, nRows, nColumns);
        //#endif

        int x, y;
        //#if TYPE==3
#pragma omp parallel for private (x,y)
        for(int i=0; i<nRows*nColumns; i++){
            x = i/nColumns;
            y= i%nColumns;
            c_parallel [x][y] = a[x][y] + b[x][y];


        }


        //#endif
        printTime(start, 3,n-nThreadsMin, time,c_parallel, c_serial, nRows, nColumns);

        //#if TYPE==4
#pragma omp parallel for schedule(static,1000) private (x,y)
        for(int i=0; i<nRows*nColumns; i++){
            x = i/nColumns;
            y= i%nColumns;
            c_parallel [x][y] = a[x][y] + b[x][y];


        }


        //#endif

        printTime(start, 4,n-nThreadsMin, time, c_parallel, c_serial, nRows, nColumns);

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
        delete [] b[i];
        delete [] c_serial[i];
        delete [] c_parallel[i];

    }

     for (int n = 0; n <nTypes; ++n) {
//         delete time [n];
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




