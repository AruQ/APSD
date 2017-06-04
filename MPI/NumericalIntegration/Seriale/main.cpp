#include <stdio.h>
#include <time.h>
double calculateIntegral(double a, double b, int n, double h);

#define N 1000000000
int main (int argc, char** argv)
{

    clock_t tStart = clock();

    double a=0.0, b=10000000.0;

    double h; //base singoli trapeziodi

    double integral; //

    double start = 0, finish = 0, partialTime = 0;


    h= (b-a)/N;
    integral =calculateIntegral(a, b, N, h);

    printf("With n = %d number of trapezioids, our estimate of integral from %f to %f = %f\n",N,a,b,integral);


    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);


}


double f(double x)
{
    return x*x;


}

double calculateIntegral(double a, double b, int n, double h)
{

    double integral =(f(a) + f(b))/2.0;;
    double x = a;

    for (int i = 1; i < n; ++i) {
        x+=h;
        integral = integral + f(x);
    }
    integral*=h;
    return integral;





}
