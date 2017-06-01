#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


MPI_Status status;
inline int get (int i, int j, int nCol)
{
    return i*nCol+j;
}


//COUNT E DISPLACEMENT SERVE A TUTTI I PROCESSI?

int numnodes,myid,mpi_err;
#define mpi_root 0

#define DIMENSION 13


void init_it(int  *argc, char ***argv);

void init_it(int  *argc, char ***argv) {
    mpi_err = MPI_Init(argc,argv);
    mpi_err = MPI_Comm_size( MPI_COMM_WORLD, &numnodes );
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
}

int main(int argc,char *argv[]){
    int * local_a;
    int * local_c;
    int * b, *a, *c;
    int* counts;

    init_it(&argc,&argv);

    b = (int*)malloc(sizeof(int)*DIMENSION * DIMENSION);

    int* displs =(int*) malloc( numnodes * sizeof(int) );

    if(myid == mpi_root){

        a = (int*)malloc(sizeof(int)*DIMENSION * DIMENSION);

        for (int i= 0; i< DIMENSION*DIMENSION; i++)
        {
            a[i] = 1;
            b[i] = 2;
        }
        c = (int*)malloc(sizeof(int)*DIMENSION * DIMENSION);


        for (int dest = 1; dest<numnodes; ++dest) {
            MPI_Send(b, DIMENSION*DIMENSION, MPI_INT, dest, mpi_root, MPI_COMM_WORLD);

        }

    }

    counts=(int*)malloc(sizeof(int)*numnodes);
    for (int i = 0; i < numnodes; ++i) {
        counts[i] = DIMENSION/(numnodes);
        if (i == numnodes-1)
            counts[i] += DIMENSION%numnodes;

        counts[i] *=DIMENSION;


    }

    displs[0]=0;
    for(int i=1;i<numnodes;i++){
        displs[i]=counts[i-1]+displs[i-1];
    }



    local_a = (int*)malloc(sizeof(int)*counts[myid]);
    local_c = (int*)malloc(sizeof(int)*counts[myid]);


//    printf("sono count di myid  %d %d\n", myid,counts[myid]);

    MPI_Scatterv(a,counts,displs, MPI_INT,local_a,counts[myid],MPI_INT, mpi_root,MPI_COMM_WORLD);

    if (myid != mpi_root)
    {

        MPI_Recv(b, DIMENSION*DIMENSION, MPI_INT, mpi_root, 0, MPI_COMM_WORLD, &status);


    }
    for (int i = 0; i < counts[myid]; ++i) {

        local_c[i] = 0;
    }

    for (int i = 0; i < counts[myid]/DIMENSION; ++i) {
        for (int j = 0; j < DIMENSION; ++j) {
            for (int k = 0; k < DIMENSION; ++k) {
                local_c[get(i,j,DIMENSION)] += local_a[get(i,k,DIMENSION)] * b[get(k,j,DIMENSION)];

            }
        }

    }

    mpi_err = MPI_Gatherv(local_c,counts[myid], MPI_INT, c,counts ,displs, MPI_INT, mpi_root,MPI_COMM_WORLD);

    if (myid == mpi_root)
    {
        for (int i = 0; i< DIMENSION*DIMENSION; i++)
        {
            printf ("%d ", c[i]);
            if ((i+1)%DIMENSION ==0)
                printf ("\n");
        }
    }

    mpi_err = MPI_Finalize();

    if (myid == mpi_root)
    {
        delete [] a;
        delete [] b;
        delete [] c;
    }
    delete [] local_a;
    delete [] local_c;
    delete [] counts;
    delete [] displs;
}
