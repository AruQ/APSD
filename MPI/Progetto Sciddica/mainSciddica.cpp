#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include "Reader.h"

#define NUMBER_OF_DIMENSIONS 2

#define VERTICAL 0
#define HORIZONTAL 1


enum HALO_TYPE{UP=0,DOWN,LEFT,RIGHT};


struct InfoBlock
{
    int first_x = 0;
    int first_y = 0;
    int size_x = 0;
    int size_y = 0;
    int last_x = 0;
    int last_y = 0;
    bool halo[4] = {false, false, false,false};
    int* cart_coordinates;
    int* cart_dimensions;
    int rank;
};


struct InfoHalo
{
    int   	neighbor[4];
    MPI_Request requests[4];
};



class Substate
{
private:

    static int ID;
    const int neighborhoodPattern[4][2]=
    {{-1,0},{0,-1}, {0,1},{1,0}}; //CONTROLARE ORDINE
    float * current;
    float * next;
    int size;
    InfoBlock* infoBlock;
    InfoHalo infoSendHalo;
    InfoHalo infoRecvHalo;
    unsigned int myid;

    MPI_Datatype MPI_VERTICAL_BORDER;
    MPI_Datatype MPI_HORIZONTAL_BORDER;





    float * halos[4];

public:

    Substate(InfoBlock* infoBlock)
    {
        setInfoBlock(infoBlock);



    }

    Substate()
    {

    }

    ~Substate()
    {

        if(current != NULL)
            delete [] current;
        if(next != NULL)
            delete [] next;

        for (int i = 0; i < 4; ++i) {
            if (halos[i] != NULL)
            {
                delete halos[i];
            }
        }
    }

    void setMatrix (float * matrix) //non andare in segfault su submatrix non allocato(se è vero)
    {

        for (int i = 0; i < size; ++i) {
            current[i] = matrix[i];
            next[i] = matrix[i];
        }
    }

    void initRoot(float * matrix, int size)
    {
        int starterIndex = infoBlock->first_y*size + infoBlock->first_x;


        for (int i = 0; i < infoBlock->size_y; ++i) {
            for (int j = 0; j < infoBlock->size_x; ++j) {
                current[i*infoBlock->size_x+j] = matrix[starterIndex];
                next[i*infoBlock->size_x+j] = matrix[starterIndex];
                starterIndex++;

            }
            starterIndex=infoBlock->first_y*size + infoBlock->first_x +(size*(i+1));
        }
    }


    void parallelForSwap() {
#pragma omp parallel for
        for (int i = 0; i < size; i++) {

            current[i] = next[i];
        }

    }

    void setInfoBlock(InfoBlock * _infoBlock)
    {

        myid = ID++;
        infoBlock = _infoBlock;
        size = infoBlock->size_x * infoBlock->size_y;
        current= new float[size];
        next = new float [size];

        for (int i = 0; i < 4; ++i)
        {
            if (infoBlock->halo[i] && ((i == HALO_TYPE::UP)|| (i == HALO_TYPE::DOWN)))
            {
                halos[i] = new float [infoBlock->size_x];

                for (int j = 0; j < infoBlock->size_x; ++j) {
                    halos[i][j] = 0.0f;
                }
            }

            else if (infoBlock->halo[i] && ((i == HALO_TYPE::LEFT)|| (i == HALO_TYPE::RIGHT)))
            {
                halos[i] = new float [infoBlock->size_y];
                for (int j = 0; j < infoBlock->size_y; ++j) {
                    halos[i][j] = 0.0f;
                }
            }
            else/* if (!infoBlock->halo[i])*/
            {
                halos[i]= NULL;
            }
        }

        MPI_Type_vector(infoBlock->size_y, 1, infoBlock->size_x, MPI_FLOAT, &MPI_VERTICAL_BORDER);
        MPI_Type_commit(&MPI_VERTICAL_BORDER);

        MPI_Type_contiguous(infoBlock->size_x, MPI_FLOAT, &MPI_HORIZONTAL_BORDER);
        MPI_Type_commit(&MPI_HORIZONTAL_BORDER);
    }

    bool get(int i,int j,  float & val)
    {
        if (i>=0 && i<infoBlock->size_y && j>=0 && infoBlock->size_x)
        {
            val= current [i*infoBlock->size_x + j];
            return true;
        }
        return false;

    }
    bool get(int i, float & val)
    {
        if (i>=0 && i<infoBlock->size_y *infoBlock->size_x)
        {
            val= current [i];
            return true;
        }
        return false;

    }

    bool getX(int i, int j, int n, float & val)
    {

        if (n>4)
            return false;
        int newI = i+neighborhoodPattern[n][0];
        int newJ = j+neighborhoodPattern[n][1];
        if (get(newI, newJ, val)) //se il vicino si trova nella sottomatrice
        {
            return true;
        }
        if(halos[n]== NULL) //se il vicino è fuori dai bordi
        {
            return false;
        }
        val= halos[n][j]; //assegna il valore
        return true;

    }

    void init(float val)
    {
        for (int i = 0; i < size; ++i) {
            current[i] = val;
            next[i] = val;
        }
    }
    bool set (int i, int j, float & val)
    {
        if (i>=0 && i<infoBlock->size_y && j>=0 && infoBlock->size_x)
        {
            next [i*infoBlock->size_x + j]= val;
            return true;
        }
        return false;
    }

    bool set (int i,float & val)
    {
        if (i>=0 && i<infoBlock->size_y *infoBlock->size_x)
        {
            next [i]= val;
            return true;
        }
        return false;
    }



    void block_receiving(int & MPI_root, MPI_Comm & MPI_COMM_CUBE, int tag)
    {
        MPI_Status status;
        MPI_Recv(current, infoBlock->size_x*infoBlock->size_y, MPI_FLOAT, MPI_root, tag, MPI_COMM_CUBE, &status);
        for (int i = 0; i < infoBlock->size_x* infoBlock->size_y; ++i) {
            next[i] = current[i];
        }
    }

    void send_halos(MPI_Comm & MPI_COMM_CUBE)
    {

        MPI_Cart_shift(MPI_COMM_CUBE, VERTICAL, 1, &infoSendHalo.neighbor[0], &infoSendHalo.neighbor[1]);
        MPI_Cart_shift(MPI_COMM_CUBE, HORIZONTAL, 1, &infoSendHalo.neighbor[2], &infoSendHalo.neighbor[3]);


        int starterIndex [4] = {0, (infoBlock->size_y-1)*infoBlock->size_x, 0,infoBlock->size_x-1};
        for (int i = 0; i < 4; ++i) {
            //            cout<<infoSendHalo.neighbor[i]<<endl;

            if (infoSendHalo.neighbor[i] != -1)
            {
                int tag = myid+myid*i+infoSendHalo.neighbor[i];

                if (i ==UP || i == DOWN)
                    MPI_Isend(&current[starterIndex[i]],1,MPI_HORIZONTAL_BORDER,infoSendHalo.neighbor[i],tag,MPI_COMM_CUBE,&infoSendHalo.requests[i]);
                else
                    MPI_Isend(&current[starterIndex[i]],1,MPI_VERTICAL_BORDER,infoSendHalo.neighbor[i],tag,MPI_COMM_CUBE,&infoSendHalo.requests[i]);


            }
        }
    }



    void recv_halos(MPI_Comm & MPI_COMM_CUBE)
    {
        MPI_Cart_shift(MPI_COMM_CUBE, VERTICAL, 1, &infoRecvHalo.neighbor[0], &infoRecvHalo.neighbor[1]);
        MPI_Cart_shift(MPI_COMM_CUBE, HORIZONTAL, 1, &infoRecvHalo.neighbor[2], &infoRecvHalo.neighbor[3]);




        for (int i = 0; i < 4; ++i) {


            if (infoSendHalo.neighbor[i] != -1)
            {
                int sender = (i%2==0? i+1:i-1);
                int tag = myid+myid*sender+infoBlock->rank;

                MPI_Irecv(halos[i],infoBlock->size_x,MPI_FLOAT,infoRecvHalo.neighbor[i],tag,MPI_COMM_CUBE,&infoRecvHalo.requests[i]);

            }
        }


    }

    bool sendCompleted()
    {
        for (int i = 0; i < 4; ++i) {
            MPI_Status status;
            if (infoBlock->halo[i])
                MPI_Wait(&infoSendHalo.requests[i], &status);
        }
    }

    bool recvCompleted()
    {

        for (int i = 0; i < 4; ++i) {
            MPI_Status status;
            if (infoBlock->halo[i])
            {
                MPI_Wait(&infoRecvHalo.requests[i], &status);
            }
        }

    }

    void stampaHalos ()
    {
        for (int i = 0; i < 4; ++i) {
            if (infoBlock->halo[i])
            {
                printf("halo i=%d del rank (%d,%d)\n", i, infoBlock->cart_coordinates[0], infoBlock->cart_coordinates[1]);
                if (i == UP ||i == DOWN)
                {

                    for(int j=0; j< infoBlock->size_x; j++)
                    {

                        printf("%f ", halos[i][j]);

                    }
                }
                else
                    for(int j=0; j< infoBlock->size_y; j++)
                    {

                        printf("%f ", halos[i][j]);

                    }


            }
            printf("\n");


        }
    }

    friend std::ostream & operator <<( std::ostream &os, const Substate &substate )
    {

        os<<"Matrix di size: "<<substate.infoBlock->size_y<<" X "<<substate.infoBlock->size_x<<"\n";

        for (int i = 0; i < substate.infoBlock->size_y; ++i) {
            for (int j = 0; j < substate.infoBlock->size_x; ++j) {
                os<<substate.current[i*substate.infoBlock->size_x+ j]<<" ";


            }
            os<<"\n";
        }
        return os;
    }




};

int Substate::ID =10;

#define NUMBER_OF_OUTFLOWS 4

class CellularAutomata
{
private:
    Substate altitude;
    Substate debrids;
    Substate f[NUMBER_OF_OUTFLOWS];

    double epsilon;
    double r;

    InfoBlock* infoBlock;

public:
    CellularAutomata (InfoBlock* infoBlock) : infoBlock(infoBlock),altitude(infoBlock),
        debrids(infoBlock)
    {
        f[0].setInfoBlock(infoBlock);
        f[1].setInfoBlock(infoBlock);
        f[2].setInfoBlock(infoBlock);
        f[3].setInfoBlock(infoBlock);
    }

    void transitionFunction (MPI_Comm & MPI_COMM_CUBE)
    {

        for (int i = 0; i < infoBlock->size_x * infoBlock->size_y; ++i) {
            float val= 0.0f;
            debrids.get(i,val);
            val++;
            debrids.set(i,val);

        }

    }


    void init (int & MPI_root, MPI_Comm & MPI_COMM_CUBE)
    {
        altitude.block_receiving(MPI_root, MPI_COMM_CUBE, 0);
        debrids.block_receiving(MPI_root, MPI_COMM_CUBE, 1);
    }


    void init (float * z, float * h)
    {
        this->debrids.setMatrix(h);
        this->altitude.setMatrix(z);
    }

    void initRoot (float * z, float * h, int size)
    {
        this->debrids.initRoot(h,size);
        this->altitude.initRoot(z,size);
    }

    void run (unsigned int STEPS, MPI_Comm & MPI_COMM_CUBE)
    {
        int step = 0;


        altitude.send_halos(MPI_COMM_CUBE);
        altitude.recv_halos(MPI_COMM_CUBE);

        recvCompleted(&altitude);

        debrids.send_halos(MPI_COMM_CUBE);
        debrids.recv_halos(MPI_COMM_CUBE);


        //        //        sendCompleted(&debrids);
        recvCompleted(&debrids);

        while(step < STEPS)
        {
            transitionFunction(MPI_COMM_CUBE);


            debrids.parallelForSwap();
            //            cout<<debrids<<endl;


            debrids.send_halos(MPI_COMM_CUBE);
            debrids.recv_halos(MPI_COMM_CUBE);

            sendCompleted(&debrids);
            recvCompleted(&debrids);

            steering();




            step++;
        }
        debrids.stampaHalos();
    }

    void steering()
    {

    }

    void sendCompleted(Substate* substate)
    {
        substate->sendCompleted();
    }

    void recvCompleted(Substate* substate)
    {
        substate->recvCompleted();
    }




    friend std::ostream & operator <<( std::ostream &os, const CellularAutomata &CA )
    {

        os<<CA.debrids;
        return os;
    }





};

void stampa (InfoBlock infoBlock,int rank)
{
    printf("\n rank: %d \n first_col %d first_row %d \n", rank,infoBlock.first_x, infoBlock.first_y);
    printf("last_col %d last_row %d \n", infoBlock.last_x, infoBlock.last_y);
    printf("size_col %d size_row %d \n\n", infoBlock.size_x, infoBlock.size_y);
}

void block_distribution (InfoBlock & infoBlock,unsigned int size_x, unsigned int size_y);

void block_sending(float* data, InfoBlock * infoBlock, int num_procs,int size_x, MPI_Datatype & MPI_BLOCK_TYPE, int & MPI_root, MPI_Comm & MPI_COMM_CUBE, int tag);

int main(int argc, char *argv[])
{


    int num_procs = 0;
    int  rank = 0;
    //    int size = 0;



    int   	neighbor_down   	= 0;
    int   	neighbor_up   	= 0;
    int   	neighbor_left   	= 0;
    int   	neighbor_right   	= 0;


    int cart_coordinates[NUMBER_OF_DIMENSIONS] = {0, 0};
    int  cart_dimensions[NUMBER_OF_DIMENSIONS] = {0, 0};
    int  cart_periodicity[NUMBER_OF_DIMENSIONS]= {0, 0};

    MPI_Comm MPI_COMM_CUBE;
    //    MPI_Datatype MPI_BORDER;
    //    MPI_Datatype MPI_COORDINATES;
    MPI_Datatype MPI_BLOCK_TYPE;
    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    //    MPI_Barrier(MPI_COMM_WORLD);

    const char* pathZ = argv[1];
    const char* pathH = argv[2];

    int totalSteps = atoi(argv[3]);
    int stepOffset = atoi(argv[4]);

    /* Get the number of processes created by MPI and their rank */
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int MPI_root = num_procs -1;

    /* Create a 2D cartesian topology of the processes */
    MPI_Dims_create(num_procs, NUMBER_OF_DIMENSIONS, cart_dimensions);

    MPI_Cart_create(MPI_COMM_WORLD, NUMBER_OF_DIMENSIONS, cart_dimensions, cart_periodicity, 0, &MPI_COMM_CUBE);



    /* Get relevant data from the created topology */
    MPI_Cart_coords(MPI_COMM_CUBE, rank, NUMBER_OF_DIMENSIONS, cart_coordinates);
    //    printf("rank %d : %d cciao %d \n", rank, cart_coordinates[0], cart_coordinates[1]);
    MPI_Cart_rank(MPI_COMM_CUBE, cart_coordinates, &rank);
    MPI_Cart_shift(MPI_COMM_CUBE, VERTICAL, 1, &neighbor_up, &neighbor_down);
    //    printf("rank %d : %d cciao %d \n", rank, cart_coordinates[0], cart_coordinates[1]);
    //    printf("rank %d : up %d down %d \n", rank, neighbor_up, neighbor_down);

    MPI_Cart_shift(MPI_COMM_CUBE, HORIZONTAL, 1, &neighbor_left, &neighbor_right);



    InfoBlock infoBlock;
    infoBlock.cart_coordinates= cart_coordinates;
    infoBlock.cart_dimensions = cart_dimensions;
    infoBlock.rank = rank;
    //    unsigned int size = ;

    Reader readerZ(pathZ);
    Reader readerH(pathH);
    if (rank == MPI_root)
    {
        readerZ.loadFromFile();
        readerH.loadFromFile();
        readerZ.fillMatrix();
        readerH.fillMatrix();

                std::cout<<readerZ<<std::endl;
                std::cout<<readerH<<std::endl;


    }
    MPI_Bcast(&readerZ.nCols, 1, MPI_INT,MPI_root, MPI_COMM_CUBE);
    MPI_Bcast(&readerZ.nRows, 1, MPI_INT,MPI_root, MPI_COMM_CUBE);

    if (rank != MPI_root)
    {
        //        MPI_Status status;
        //        MPI_Recv(&reader.getCols(), 1, MPI_INT, MPI_root, 0, MPI_COMM_CUBE, &status);
        //        MPI_Recv(&reader.getRows(), 1, MPI_INT, MPI_root, 0, MPI_COMM_CUBE, &status);
    }

    block_distribution(infoBlock, readerZ.nCols, readerZ.nRows);
    //        stampa(infoBlock, rank);
    CellularAutomata sciddica (&infoBlock);


    if (rank == MPI_root)
    {
        block_sending(readerZ.getDataLinear(),&infoBlock,num_procs,readerZ.nCols,MPI_BLOCK_TYPE,MPI_root,MPI_COMM_CUBE,0);
        block_sending(readerH.getDataLinear(),&infoBlock,num_procs,readerZ.nCols,MPI_BLOCK_TYPE,MPI_root,MPI_COMM_CUBE,1);
        sciddica.initRoot(readerZ.getDataLinear(),readerH.getDataLinear(), readerZ.nCols);

    }
    if (rank != MPI_root)
    {
        //TODO LA ROOT DEVE RIEMPIRE A MANO LA SUA

        //fare metodo che ricolleziona i dati nella matrice globale


        sciddica.init(MPI_root, MPI_COMM_CUBE);
        //            substate.block_receiving(MPI_root,MPI_COMM_CUBE,1);



    }

    //    cout<<"sono rank : "<< rank<<" \n"<<sciddica<<std::endl;

    MPI_Barrier (MPI_COMM_CUBE);
    sciddica.run(10,MPI_COMM_CUBE);


    //    sciddica.transitionFunction(MPI_COMM_CUBE);


    //    block_sending(reader, infoBlock, MPI_BLOCK_TYPE, size);

    //TODO CREARE BLOCCHI PER OGNI PROCESSO

    MPI_Finalize();



    return 0;
}







void block_distribution (InfoBlock & infoBlock,unsigned int size_x, unsigned int size_y)
{

    infoBlock.size_x = size_x / infoBlock.cart_dimensions[1];
    infoBlock.size_y = size_y / infoBlock.cart_dimensions[0];

    infoBlock.first_x = infoBlock.size_x * infoBlock.cart_coordinates[1];
    infoBlock.first_y = infoBlock.size_y * infoBlock.cart_coordinates[0];

    infoBlock.last_x = infoBlock.first_x+ infoBlock.size_x -1;
    infoBlock.last_y = infoBlock.first_y+ infoBlock.size_y -1;


    if(infoBlock.first_x !=0)
    {

        infoBlock.halo[LEFT] = true;
        //        infoBlock.first_x --;
        //        infoBlock.size_x++;
        //        infoBlock.halo_size_x ++;
    }

    if(infoBlock.first_y !=0)
    {
        infoBlock.halo[UP] = true;
        //        infoBlock.first_y --;
        //        infoBlock.size_y++;
        //        infoBlock.halo_size_y ++;
    }

    if(infoBlock.last_x !=size_x -1)
    {
        infoBlock.halo[RIGHT] = true;
        //        infoBlock.last_x ++;
        //        infoBlock.size_x++;
        //        infoBlock.halo_size_x ++;
    }

    if(infoBlock.last_y !=size_y -1)
    {
        infoBlock.halo[DOWN] = true;
        //        infoBlock.last_y ++;
        //        infoBlock.size_y++;
        //        infoBlock.halo_size_y ++;
    }
}


void block_sending(float* data,InfoBlock * infoBlock, int num_procs, int size_x, MPI_Datatype & MPI_BLOCK_TYPE, int & MPI_root, MPI_Comm & MPI_COMM_CUBE, int tag)
{


    //DEVE FARLO IL ROOT

    MPI_Type_vector(infoBlock->size_y, infoBlock->size_x, size_x, MPI_FLOAT, &MPI_BLOCK_TYPE);
    MPI_Type_commit(&MPI_BLOCK_TYPE);


    int starterIndex = 0;
    for (int dest =0; dest<num_procs; dest++ )
    {
        if(dest== MPI_root)
            continue;


        //        cout<<"sono rank="<<rank <<" e inizio da "<<starterIndex<<endl;
        MPI_Send (&data[starterIndex], 1, MPI_BLOCK_TYPE, dest, tag, MPI_COMM_CUBE);
        if ((dest+1) % infoBlock->cart_dimensions[0] ==0)
        {
            starterIndex = (dest+1) / infoBlock->cart_dimensions[1] * size_x * infoBlock->size_y;

        }
        else
            starterIndex += infoBlock->size_x;
    }







}


