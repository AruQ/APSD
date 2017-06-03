#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
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
  int * array = (int*) malloc (sizeof(int)*DIMENSION);
  double start = 0, finish = 0, partialTime = 0;

  int rank,size;

  for (int i= 0; i< DIMENSION; i++)
  {
      array[i] = i%(DIMENSION/2);

  }
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  start = MPI_Wtime();
  computeSerial(array,DIMENSION,200000000);
  finish = MPI_Wtime();
  partialTime = finish - start;

  MPI_Finalize();
  printf("Time = %f\n",partialTime );
  free(array);
}
