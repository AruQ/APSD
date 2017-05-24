#include <stdio.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char* argv[]) {
    int         my_rank;
    int         numprocs;
    int         source;
    int         dest;
    int         data=0;
    int         tag = 0;
    MPI_Status  status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    if (my_rank == 0) {
        for (dest = 1; dest < numprocs; dest++) {

            MPI_Send(&dest, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
        }

        for (dest = 1; dest < numprocs; dest++) {

            MPI_Recv(&data, 1, MPI_INT, dest, tag, MPI_COMM_WORLD, &status);
            printf("Greetings from process %d! Value: %d\n",
                   dest, data);

        }

    } else { /* my_rank != 0 */
        MPI_Recv(&data, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        data = data*2;
        MPI_Send(&data, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
    }
    MPI_Finalize();
} /* main */
