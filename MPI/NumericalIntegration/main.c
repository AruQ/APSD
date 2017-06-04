#include <stdio.h>

#include "mpi.h"

double calculateIntegral(double a, double b, int n, double h);
double get_max_time(double partialTime, int my_rank, int p);

#define N 1000000000

int main (int argc, char** argv)
{
    int my_rank, numprocs;

    double a=0.0, b=10000000.0;
    double h; //base singoli trapeziodi

    double local_a, local_b, local_n;

    double integral; //
    double total;
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


    h= (b-a)/N;
    local_n = N/numprocs;

    local_a = a + h*local_n* my_rank;

    if (my_rank == numprocs-1)
    {
        local_n += N%numprocs;
    }
    local_b = local_a + local_n * h;


    integral =calculateIntegral(local_a, local_b, local_n, h);

    if(my_rank == 0)
    {
        total = integral;
        for(source = 1; source < numprocs ;source++)
        {
            MPI_Recv(&integral,1,MPI_DOUBLE,source,tag,MPI_COMM_WORLD,&status);
            total += integral;
        }
    }
    else
    {
        MPI_Send(&integral,1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();
    partialTime = finish - start;

    if(my_rank == 0)
    {
        printf("With n = %d number of trapezioids, our estimate of integral from %f to %f = %f\n",N,a,b,total);
        printf("Time = %f\n",partialTime );
    }



    MPI_Finalize();


}


double f(double x)
{
    return x*x;


}

double calculateIntegral(double a, double b, int n, double h)
{

    double integral =(f(a) + f(b))/2.0;
    double x = a;

    for (int i = 1; i < n; ++i) {
        x+=h;
        integral = integral + f(x);
    }
    integral*=h;
    return integral;





}
