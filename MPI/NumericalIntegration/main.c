#include <stdio.h>

#include "mpi.h"

float calculateIntegral(float& a, float& b, int& n, float& h);

int main (int argc, char** argv)
{
    int my_rank, numprocs;

    float a=0.0, b=1.0;
    int n= 1024; //numeri di trapezoidi
    float h; //base singoli trapeziodi
    float local_a, local_b, local_n;

    float integral; //
    float total;
    int source;
    int dest=0;
    int tag=0;

    MPI_Status status;

    MPI_Init (&argc,&argv);

    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);


    h= (b-a)/n;
    local_n = n/numprocs;

    local_a = a + h*local_n* my_rank;
    local_b = local_a + local_n * h;

    integral =calculateIntegral(local_a, local_b, local_n, h);






}
float f(float x)
{
    return x*x;


}

float calculateIntegral(float& a, float& b, int& n, float& h)
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
