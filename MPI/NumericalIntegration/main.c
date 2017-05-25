#include <stdio.h>

#include "mpi.h"

float calculateIntegral(float a, float b, int n, float h);
double get_max_time(double partialTime, int my_rank, int p);

int main (int argc, char** argv)
{
    int my_rank, numprocs;

    float a=0.0, b=9.0;
    int n= 999999999; //numeri di trapezoidi
    float h; //base singoli trapeziodi
    float local_a, local_b, local_n;

    float integral; //
    float total;
    int source;
    int dest=0;
    int tag=0;
    double start = 0, finish = 0, partialTime = 0;

    MPI_Status status;

    MPI_Init (&argc,&argv);

    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);


    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    h= (b-a)/n;
    local_n = n/numprocs;

    local_a = a + h*local_n* my_rank;
    local_b = local_a + local_n * h;

    integral =calculateIntegral(local_a, local_b, local_n, h);

    if(my_rank == 0)
    {
      total = integral;
      for(source = 1; source < numprocs ;source++)
      {
        MPI_Recv(&integral,1,MPI_FLOAT,source,tag,MPI_COMM_WORLD,&status);
        total += integral;
      }
    }
    else
    {
      MPI_Send(&integral,1,MPI_FLOAT,dest,tag,MPI_COMM_WORLD);
    }

    finish = MPI_Wtime();

    partialTime = finish - start;

    partialTime = get_max_time(partialTime, my_rank, numprocs);

    if(my_rank == 0)
    {
      printf("With n = %d number of trapezioids, our estimate of integral from %f to %f = %f\n",n,a,b,total);
      printf("Time = %f\n",partialTime );
    }



    MPI_Finalize();


}

double get_max_time(double partialTime, int my_rank, int p) {
   MPI_Status status;
   double temp;

   if (my_rank == 0) {
      for (int source = 1; source < p; source++) {
         MPI_Recv(&temp, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
         if (temp > partialTime) partialTime = temp;
      }
   } else {
      MPI_Send(&partialTime, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
   }
   return partialTime;
}

float f(float x)
{
    return x*x;


}

float calculateIntegral(float a, float b, int n, float h)
{

    float integral =(f(a) + f(b))/2.0;;
    float x = a;

    for (int i = 1; i < n; ++i) {
        x+=h;
        integral = integral + f(x);
    }
    integral*=h;
    return integral;





}
