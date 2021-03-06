//Find an element in a vector / matrix (with blocking and non-blocking messages)

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#define DIMENSION 200000000
#define mpi_root 0
#define toFind 300

int main(int argc, char* argv[]) {

    int rank,size;
    MPI_Status status;
    MPI_Request request;
    int *array, *local_array;

    double start = 0, finish = 0, partialTime = 0;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    int local_size = DIMENSION/size;
    if (rank == size-1)
    {
        local_size += DIMENSION%size;
    }

    if (rank==mpi_root) {
        array = (int*) malloc (sizeof(int)*DIMENSION);

        for (int i= 0; i< DIMENSION; i++)
        {
            array[i] = i%(DIMENSION/10);

        }
    }

    local_array = (int*) malloc (sizeof(int)*local_size);


    MPI_Scatter(array,local_size, MPI_INT,local_array,local_size,MPI_INT, mpi_root,MPI_COMM_WORLD);


    int found = -1;
    int done;
    int globalFound = -1;

    MPI_Irecv(&globalFound,1,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&request);
    MPI_Test(&request,&done,&status);

    for(int i = 0; i< local_size && !done; i++)
    {
        if (local_array[i] == toFind)
        {
            found = rank* local_size+i;
            for(int n=0;n<size;++n) {
                MPI_Send(&found,1,MPI_INT,n,0,MPI_COMM_WORLD);
            }
        }
        MPI_Test(&request,&done,&status);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();
    partialTime = finish - start;

    if(rank == mpi_root)
    {
      printf("Time = %f\n",partialTime );
    }

    if (found != -1) {
        printf("P:%d stopped at global index %d \n",rank, found);
    }

    MPI_Finalize();

    if (rank==0)
        free(array);
    free (local_array);

    return 0;
}
