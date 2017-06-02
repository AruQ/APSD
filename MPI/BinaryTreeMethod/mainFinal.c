#include <mpi.h>
#include<cmath>

#include <stdio.h>
#include <stdlib.h>

#define SIZE 8

MPI_Status status;

#define mpi_root 0
inline int log2(int value)
{
    return log(value) / log(2);
}

int binary_tree_method (int value, int rank, int numproc)
{

    int treeDepth = log2(SIZE);
    int increment = 1;

    int receivedValue;

    for (int i = 0; i< treeDepth; i++)
    {
        for (int j = 0;  j+increment< numproc; j+=increment)
        {
            int receiver = j;
            int sender = j+increment;

            if (rank == receiver)
            {
                MPI_Recv(&receivedValue, 1,
                         MPI_INT, sender, 0, MPI_COMM_WORLD, &status);

                //                printf("sono il processo %d ho ricevuto da %d il numero %d e lo sommo a %d\n", rank, sender, receivedValue, value);
                value += receivedValue;
            }
            else if (rank == sender)
            {
                MPI_Send(&value, 1, MPI_INT,
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


    int * data;
    //    int * local_data;
    int local_size = (SIZE/numproc);
    //    if (rank == numproc-1)
    //    {
    //        local_size += SIZE%numproc;
    //    }

//    printf ("local size %d\n", local_size);
    data = (int*) malloc (sizeof(int) * SIZE);
    for (int i = 0; i < SIZE; i++)
        data[i] = i;
    if (rank == mpi_root)
    {

        for (int i=1; i<numproc; i++)
        {
            //            if (i == numproc-1)
            //                MPI_Send(&data[i*local_size], local_size+SIZE%numproc, MPI_INT,
            //                        i, mpi_root, MPI_COMM_WORLD);

            //            else
            MPI_Send(&data[i*local_size], local_size, MPI_INT,
                    i, mpi_root, MPI_COMM_WORLD);

        }
    }

    if (rank != mpi_root)
    {
//        data = (int*) malloc (sizeof(int) * local_size);
//        MPI_Recv(data, local_size,
//                 MPI_INT, mpi_root, 0, MPI_COMM_WORLD, &status);

    }

//    int localSum =0;
//    for (int i = 0; i < local_size; ++i) {
//        //        if (rank==0)
//        localSum+= data[i];


//    }

    //    printf ("rank %d local size %d local sum %d \n",rank, local_size, localSum);

    // each rank only gets one entry, and
    // they need to sum them by sending messages
    int result = binary_tree_method(data[rank], rank, numproc);
    MPI_Barrier(MPI_COMM_WORLD);

    // Compute the correct result
    int sum = 0;
    for (int i = 0; i < SIZE; i++)
        sum += data[i];


    if (rank == mpi_root) {
        printf("MPI Result     = %d\n", result);
        printf("Correct Result = %d\n", sum);
    }
    MPI::Finalize();

    //    if (rank==0)
    free (data);

    //    else
    //        free(local_data);
}
