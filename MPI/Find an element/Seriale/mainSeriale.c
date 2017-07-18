#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int computeSerial(int * a, int nElements, int toFind) {
    int occurences = 0;
    for (int i = 0; i < nElements; ++i) {
        if (a[i] == toFind)
            occurences++;
    }
    return occurences;

}

int main(int argc, char* argv[]) {

    int DIMENSION = 200000000;

    clock_t tStart = clock();

    int * array = (int*) malloc (sizeof(int)*DIMENSION);
    double start = 0, finish = 0, partialTime = 0;

    int rank,size;

    for (int i= 0; i< DIMENSION; i++)
    {
        array[i] = i%(DIMENSION/2);

    }


    computeSerial(array,DIMENSION,1);

    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    free(array);
}
