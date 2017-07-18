#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 240000000
//#define SIZE 1000000000
int main (int argc, char** argv)
{

    clock_t tStart = clock();

    double* data = (double*) malloc (sizeof(double) * SIZE);
    double* index = data;
    for (double i = 0.0; i < SIZE; i+=1.0)
    {
        *index=i;
        index ++;
    }

    double sum=0.0;
    index=data;

    //     Compute the correct result
    for (double i = 0.0; i < SIZE; i+=1.0)
    {
	sum += *index;
	index++;
        
    }

    printf("La somma seriale è: %15.0f\n",sum);
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    free (data);


}
