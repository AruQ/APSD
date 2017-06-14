#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>


MPI_Status status;
inline int get (int i, int j, int nCol)
{
    return i*nCol+j;
}



int numnodes,myid,mpi_err;
#define mpi_root 0

#define DIMENSION 8


void init_it(int  *argc, char ***argv);

void init_it(int  *argc, char ***argv) {
    mpi_err = MPI_Init(argc,argv);
    mpi_err = MPI_Comm_size( MPI_COMM_WORLD, &numnodes );
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
}

int main(int argc,char *argv[]){
    int * local_a;
    int * local_c;
    int * local_b;
    int * b, *a, *c;
    int* counts;
    double start = 0, finish = 0, partialTime = 0;

    init_it(&argc,&argv);

    MPI_Datatype columnsType;
    MPI_Datatype subMatrixType;
    MPI_Status status;


    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();


    b = (int*)malloc(sizeof(int)*DIMENSION * DIMENSION);

    //    int* displs =(int*) malloc( numnodes * sizeof(int) );





    //    counts=(int*)malloc(sizeof(int)*numnodes);

    //    for (int i = 0; i < numnodes; ++i) {
    //        counts[i] = DIMENSION/(numnodes);
    //        if (i == numnodes-1)
    //            counts[i] += DIMENSION%numnodes;

    //        counts[i] *=DIMENSION;


    //    }

    //    displs[0]=0;
    //    for(int i=1;i<numnodes;i++){
    //        displs[i]=counts[i-1]+displs[i-1];
    //    }



    //    printf("sono count di myid  %d count: %d e il mio displ è %d\n", myid,counts[myid], displs[myid]);


    int nColumns = DIMENSION / sqrt(numnodes-1);
    int stride = DIMENSION-nColumns;
    if(myid == mpi_root){

        a = (int*)malloc(sizeof(int)*DIMENSION * DIMENSION);

        for (int i= 0; i< DIMENSION*DIMENSION; i++)
        {
            a[i] = 1;
            b[i] = 2;
        }
        c = (int*)malloc(sizeof(int)*DIMENSION * DIMENSION);

        printf("sono qui\n");

        MPI_Type_vector(DIMENSION, nColumns, stride,MPI_INT,&columnsType);
        MPI_Type_commit(&columnsType);

        for (int dest =1; dest<numnodes; dest++ )
        {

            int starterIndexB = (dest-1)*nColumns %DIMENSION;
            int starterIndexA = (dest-1)*nColumns *DIMENSION;
            printf ("starter index %d nCol %d stride %d \n", starterIndexB, nColumns, stride);
            //            MPI_Send (&a[starterIndexA], nColumns*DIMENSION, MPI_INT, dest, 0, MPI_COMM_WORLD);
            MPI_Send (&b[starterIndexB], nColumns, columnsType, dest, 0, MPI_COMM_WORLD);
        }
        //        for (int dest = 1; dest<numnodes; ++dest) {
        //            MPI_Send(b, DIMENSION*DIMENSION, MPI_INT, dest, mpi_root, MPI_COMM_WORLD);

        //        }

    }


    local_a = (int*)malloc(sizeof(int)*nColumns* DIMENSION);
    local_b = (int*)malloc(sizeof(int)*nColumns* DIMENSION);
    local_c = (int*)malloc(sizeof(int)*nColumns* nColumns);




    printf ("sizeof(int)*nColumns* DIMENSION %d \n", sizeof(int)*nColumns* DIMENSION);

    if (myid != mpi_root)
    {

        printf ("sono qui %d\n", nColumns * DIMENSION);
//        MPI_Recv(local_a, nColumns*DIMENSION, MPI_INT, mpi_root, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(local_b, DIMENSION, MPI_INT, mpi_root, 0, MPI_COMM_WORLD, &status);

        //        for (int i = 0; i < nColumns; ++i) {
        //            for (int j = 0; j < nColumns; ++j) {
        //                for (int k = 0; k < DIMENSION; ++k) {
        //                    local_c[get(i,j,DIMENSION)] += local_a[get(i,k,DIMENSION)] * b[get(i,k,DIMENSION)];

        //                }
        //            }

        MPI_Type_vector(nColumns, nColumns, stride, MPI_INT, &subMatrixType);
        MPI_Type_commit(&subMatrixType);

        MPI_Send (local_c, nColumns, subMatrixType, mpi_root, 1, MPI_COMM_WORLD);



    }


    if (myid == mpi_root)
    {
        int starterIndex = 0;
        for (int source = 1; source < numnodes; ++source) {


            MPI_Recv(&c[starterIndex], DIMENSION, MPI_INT, source, 1, MPI_COMM_WORLD, &status);

            if ((starterIndex+nColumns) % DIMENSION == 0)
            {
                starterIndex += DIMENSION * (nColumns-1);
            }

            starterIndex+=nColumns;
        }
    }

    //    MPI_Scatterv(a,counts,displs, MPI_INT,local_a,counts[myid],MPI_INT, mpi_root,MPI_COMM_WORLD);

    //    for (int i = 0; i < counts[myid]; ++i) {

    //        local_c[i] = 0;
    //    }

    //    for (int i = 0; i < counts[myid]/DIMENSION; ++i) {
    //        for (int j = 0; j < DIMENSION; ++j) {
    //            for (int k = 0; k < DIMENSION; ++k) {
    //                local_c[get(i,j,DIMENSION)] += local_a[get(i,k,DIMENSION)] * b[get(k,j,DIMENSION)];

    //            }
    //        }

    //    }

    //    mpi_err = MPI_Gatherv(local_c,counts[myid], MPI_INT, c,counts ,displs, MPI_INT, mpi_root,MPI_COMM_WORLD);

    if (myid == mpi_root)
    {
        //        for (int i = 0; i< DIMENSION*DIMENSION; i++)
        //        {
        //            printf ("%d ", c[i]);
        //            if ((i+1)%DIMENSION ==0)
        //                printf ("\n");
        //        }
    }

    //    MPI_Barrier(MPI_COMM_WORLD);
    //    finish = MPI_Wtime();
    //    partialTime = finish - start;

    //    if(myid == mpi_root)
    //    {
    //      printf("Time = %f\n",partialTime );
    //    }

    mpi_err = MPI_Finalize();

    if (myid == mpi_root)
    {
        delete [] a;
        delete [] b;
        delete [] c;
    }
    //    delete [] local_a;
    //    delete [] local_c;
    //    delete [] local_b;
    //    delete [] counts;
    //    delete [] displs;
}
