#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include<cassert>
#include <string>
#include<iostream>
#include <functional>

using namespace std;
#define TYPE 0


const int nTypes = 5;
string types [nTypes] = {"parallel for occurences", "parallel for found atomic","reduction","parallel for occurences atomic", "seriale"};

int computeSerial(int * a, int nElements, int toFind) {
    int occurences = 0;
    for (int i = 0; i < nElements; ++i) {
        if (a[i] == toFind)
            occurences++;
    }
    return occurences;

}

auto checkFoundResultOccurences = [](int * a, int nElements, int toFind, int occurrences)->bool  {
    return computeSerial(a, nElements, toFind) == occurrences;
};

auto checkFoundResultFound = [](int * a, int nElements, int toFind, bool found )->bool  {
    return (computeSerial(a, nElements, toFind) == 0 && !found) || (computeSerial(a, nElements, toFind) > 0 && found);
};


void printTime(double& start, int i, int j, double ** time, auto function) {

    double end = omp_get_wtime();
    if (function())
    {
        time [i][j] = end - start;
    }
    else
    {
        time [i][j] = -1;
    }
    start = omp_get_wtime();
}

int main(int argc, char *argv[]) {

    assert(argc >= 4);
    int nElements = atoi(argv[1]);
    int toFind = atoi(argv[2]);
    int nThreadsMin = atoi(argv[3]);
    int nThreadsMax = atoi(argv[4]);

    double* time[nTypes];
    for (int k = 0; k < nTypes; ++k) {
        time[k] = new double[nThreadsMax - nThreadsMin];
    }

    int* a = new int [nElements];
    for (int i = 0; i < nElements; i++) {
        a[i] = 1;
    }

    double start = omp_get_wtime();

    int occurences;
    occurences = computeSerial(a, nElements, toFind);
    //    printTime(start, nTypes - 1, 0, time, a, nElements, toFind, occurences, false);
    printTime(start, nTypes - 1, 0, time,std::bind (checkFoundResultOccurences,a,nElements,toFind,occurences));


    int indexTypes = 0;

    for (int n = nThreadsMin; n <= nThreadsMax; n++) {
        omp_set_num_threads(n);

        int arrayOccurences[n] = {0};
#pragma omp parallel
        {
            int id = omp_get_thread_num();
#pragma omp  for
            for (int i = 0; i < nElements; i++) {
                if (a[i] == toFind)
                    arrayOccurences[id]++;
            }
            int occurences = 0;
            for (int i = 0; i < n; i++)
                occurences += arrayOccurences[i];

        }
        printTime(start, indexTypes++, n - nThreadsMin, time,std::bind (checkFoundResultOccurences,a,nElements,toFind,occurences));

        bool found = false;
#pragma omp parallel for
        for (int i = 0; i < nElements; i++) {
            if (!found && a[i] == toFind) {
#pragma omp atomic
                found++;
            }
        }
        printTime(start, indexTypes++, n - nThreadsMin, time, std::bind (checkFoundResultFound,a,nElements,toFind,found));


        occurences = 0;
#pragma omp parallel for reduction (+:occurences)
        for (int i = 0; i < nElements; i++)
            if (a[i] == toFind)
                occurences++;


        printTime(start, indexTypes++, n - nThreadsMin, time,std::bind (checkFoundResultOccurences,a,nElements,toFind,occurences));

        occurences = 0;
#pragma omp parallel
        {
            int id = omp_get_thread_num();
#pragma omp  for
            for (int i = 0; i < nElements; i++) {
                if (a[i] == toFind)
#pragma omp atomic
                    occurences++;
            }

        }
        printTime(start, indexTypes++, n - nThreadsMin, time, std::bind (checkFoundResultOccurences,a,nElements,toFind,occurences));

        indexTypes = 0;
    }

    for (int k = 0; k < nTypes; ++k) {
        cout << types[k] << "; ";
        for (int n = 0; n <= nThreadsMax - nThreadsMin; ++n) {
            printf("%.4g;  ", time[k][n]);
        }
        printf("\n");
    }

    delete [] a;
    return 0;
}





