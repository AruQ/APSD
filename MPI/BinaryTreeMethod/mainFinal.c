#include <mpi.h>
#include<cmath>

#include <stdio.h>
#include <stdlib.h>

const unsigned long long int SIZE= 400000000;


MPI_Status status;

#define mpi_root 0
inline int log_2(unsigned long long int value)
{
    return log(value) / log(2);
}

long double binary_tree_method (unsigned long long int value, int rank, int numproc)
{

    int treeDepth = log_2(SIZE);
    int increment = 1;

    unsigned long long int receivedValue;

    for (int i = 0; i< treeDepth; i++)
    {
        for (int j = 0;  j+increment< numproc; j+=increment)
        {
            int receiver = j;
            int sender = j+increment;

            if (rank == receiver)
            {
                MPI_Recv(&receivedValue, 1,
                         MPI_UNSIGNED_LONG_LONG, sender, 0, MPI_COMM_WORLD, &status);
                value += receivedValue;
            }
            else if (rank == sender)
            {
                MPI_Send(&value, 1, MPI_UNSIGNED_LONG_LONG,
                         receiver, 0, MPI_COMM_WORLD);
            }

        }


        increment*=2;
    }
    return value;
}



int main(int argc, char* argv[]) {


    int numproc;
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);

    double start = 0, finish = 0, partialTime = 0;


    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();


    unsigned long long int * data;

    int local_size = (SIZE/numproc);

    if (rank == mpi_root)
    {
        data = (unsigned long long int*) malloc (sizeof(unsigned long long int) * SIZE);
        unsigned long long int * index = data;
        for (long long i = 0; i < SIZE; i++)
        {
            *index = i;
            index++;
        }
        for (int i=1; i<numproc; i++)
        {
            if (i == numproc-1)
                MPI_Send(&data[i*local_size], local_size+SIZE%numproc, MPI_UNSIGNED_LONG_LONG,
                        i, mpi_root, MPI_COMM_WORLD);
            else
                MPI_Send(&data[i*local_size], local_size, MPI_UNSIGNED_LONG_LONG,
                        i, mpi_root, MPI_COMM_WORLD);

        }
    }

    if (rank == numproc-1)
    {
        local_size += SIZE%numproc;
    }

    if (rank != mpi_root)
    {
        data = (unsigned long long int*) malloc (sizeof(unsigned long long int) * local_size);
        MPI_Recv(data, local_size,
                 MPI_UNSIGNED_LONG_LONG, mpi_root, 0, MPI_COMM_WORLD, &status);
    }

    unsigned long long int localSum =0;
    unsigned long long int * index = data;
    for (unsigned long long int i = 0; i < local_size; i++) {
        localSum+= *index;
        index++;
    }


    unsigned long long int result = binary_tree_method(localSum, rank, numproc);
        MPI_Barrier(MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();
    partialTime = finish - start;


    if(rank == mpi_root)
    {
        unsigned long long int sum =0.0;
        printf("Time = %f\n",partialTime );

        //        start = MPI_Wtime();

        index = data;
        for (unsigned long long int i = 0; i < SIZE; i++)
        {

            sum += *index;
            index++;
        }

        //        finish = MPI_Wtime();
        //        partialTime = finish - start;

        //        printf("Time Serial = %f\n",partialTime );

        printf("MPI Result     = %llu\n", result);
        printf("Correct Result = %llu\n", sum);
    }




    MPI::Finalize();


    free (data);

}
