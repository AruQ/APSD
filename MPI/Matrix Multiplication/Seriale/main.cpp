#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define DIMENSION 10000

int main (int argc, char** argv)
{
    int ** a = new int*[DIMENSION];
    int ** b = new int*[DIMENSION];
    int ** c = new int*[DIMENSION];



    for (int i=0; i< DIMENSION; i++)
    {
        a[i] = new int [DIMENSION];
        b[i] = new int [DIMENSION];

        for (int j = 0; j < DIMENSION; ++j) {
            a[i][j]= 1;
            b[i][j]= 2;
        }
        c[i] = new int [DIMENSION];
    }

    clock_t tStart = clock();
    for (int i = 0; i < DIMENSION; ++i) {
        for (int j = 0; j < DIMENSION; ++j) {
            c[i][j] = 0;
            for (int k = 0; k < DIMENSION; ++k) {
                c[i][j]+= a[i][k]*b[k][j];
            }

        }
    }


     printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);


    return 0;


}
