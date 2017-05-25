#include <stdio.h>
#include "mpi.h"
#define DIMENSION               5        /* number of rows and columns in matrix */

MPI_Status status;


main(int argc, char **argv)
{
  int ** a = new int*[DIMENSION];
  int ** b = new int*[DIMENSION];
  int ** c = new int*[DIMENSION];
  for (int i=0; i< DIMENSION; i++)
  {
    a[i] = new int [DIMENSION];
    b[i] = new int [DIMENSION];
    c[i] = new int [DIMENSION];
  }



  //    printf("sono rows = %d e sono rank= %d\n", a[0][0], a[1]);
  double start = 0, finish = 0, partialTime = 0;
  int rank, numprocs, numworkers, rows, source, offset;


  //    int numtasks,taskid,numworkers,source,dest,rows,offset,i,j,k;



  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  numworkers = numprocs-1;



  /*---------------------------- master ----------------------------*/
  if (rank == 0)
  {
    for (int i=0; i< DIMENSION; i++)
    {
      for (int j = 0; j < DIMENSION; ++j) {
        a[i][j]= 1;
        b[i][j]= 2;
      }

    }


    /* send matrix data to the worker tasks */
    rows = DIMENSION/numworkers;


    offset = 0;

    // MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    for (int dest=1; dest<=numworkers; dest++)
    {
      MPI_Send(&offset, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
      MPI_Send(&rows, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
      for(int row = 0;row<rows;row++)
        MPI_Send(&a[offset+row][0], DIMENSION, MPI_INT,dest,1, MPI_COMM_WORLD);
      for(int row = 0;row<DIMENSION;row++)
        MPI_Send(&b[row][0], DIMENSION, MPI_INT, dest, 1, MPI_COMM_WORLD);
      offset += rows;

    }


    /* wait for results from all worker tasks */
    for (int i=1; i<=numworkers; i++)
    {

      source = i;
      MPI_Recv(&offset, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
      MPI_Recv(&rows, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
      for(int row = 0;row<rows;row++)
        MPI_Recv(&c[offset+row][0], rows*DIMENSION, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
    }


    // MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();
    partialTime = finish - start;


    printf("Here is the result matrix:\n");
    for (int i=0; i<DIMENSION; i++) {
      for (int j=0; j<DIMENSION; j++)
      printf("%d   ", c[i][j]);
      printf ("\n");
    }

    printf("Time elapsed: %f\n", partialTime);



  }

  /*---------------------------- worker----------------------------*/
  if (rank > 0) {
    source = 0;
    MPI_Recv(&offset, 1, MPI_INT, source, 1, MPI_COMM_WORLD, &status);
    MPI_Recv(&rows, 1, MPI_INT, source, 1, MPI_COMM_WORLD, &status);
    for(int row = 0;row<rows;row++)
    MPI_Recv(&a[row][0], DIMENSION, MPI_INT, source, 1, MPI_COMM_WORLD, &status);
    for(int row = 0;row<DIMENSION;row++)
    MPI_Recv(&b[row][0],  DIMENSION, MPI_INT, source, 1, MPI_COMM_WORLD, &status);

    /* Matrix multiplication */
    for (int i=0; i<rows; i++) {
      for (int j=0; j<DIMENSION; j++)
      {
        c[i][j] = 0.0;
        for (int k=0; k<DIMENSION; k++)
        c[i][j] = c[i][j] + a[i][k] * b[k][j];
      }
    }


    MPI_Send(&offset, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
    MPI_Send(&rows, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
    for(int row = 0;row<rows;row++)
      MPI_Send(&c[row][0], DIMENSION, MPI_INT, 0, 2, MPI_COMM_WORLD);
  }



  MPI_Finalize();
  for (int i = 0; i < DIMENSION; ++i) {
    delete [] a[i];
    delete [] b[i];
    delete [] c[i];

  }

  delete [] a;
  delete [] b;
  delete [] c;
}
