//Find an element in a vector / matrix (with blocking and non-blocking messages)

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#define DIMENSION 200000000
#define mpi_root 0
#define toFind 3

int main(int argc, char* argv[]) {

    double start = 0, finish = 0, partialTime = 0;

    int rank,size;
    MPI_Status status;

    int *array, *local_array;
    int globalOccurences;

    int countsSerial =0;


    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    int local_size = DIMENSION/(size);
    if (rank == size-2)
    {
        local_size += DIMENSION%(size);
    }

    if (rank==mpi_root) {
        array = (int*) malloc (sizeof(int)*DIMENSION);

        for (int i= 0; i< DIMENSION; i++)
        {

            array[i] = i%(DIMENSION/2);
//            if (array[i] == toFind)
//                countsSerial++;

        }
    }

    local_array = (int*) malloc (sizeof(int)*local_size);

    MPI_Scatter(array,local_size, MPI_INT,local_array,local_size,MPI_INT, mpi_root,MPI_COMM_WORLD);


    int occurences = 0;
    for(int i = 0; i< local_size; i++)
    {
        if (local_array[i] == toFind)
        {
            occurences++;
        }
    }


    MPI_Reduce(&occurences, &globalOccurences, 1, MPI_INT, MPI_SUM, mpi_root,MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();
    partialTime = finish - start;

    if(rank == mpi_root)
    {
      printf("Time = %f\n",partialTime );
    }
    MPI_Finalize();



    if (rank==mpi_root)
    {

        printf("Number of occurrences parallel: %d\n",globalOccurences);
//        printf("Number of occurrences serial: %d\n",countsSerial);
        free(array);

    }
    free (local_array);

    return 0;
}

