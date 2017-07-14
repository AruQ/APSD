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
    bool halo[4] = {false};
    int* cart_coordinates;
    int* cart_dimensions;
};


class Substate
{
private:
    const int neighborhoodPattern[4][2]=
    {{-1,0},{0,-1}, {0,1},{1,0}}; //CONTROLARE ORDINE
    float * subMatrix;
    int size;
    InfoBlock* infoBlock;

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

        if(subMatrix != NULL)
            delete [] subMatrix;

        for (int i = 0; i < 4; ++i) {
            if (halos[i] != NULL)
            {
                delete halos[i];
            }
        }
    }

    void setInfoBlock(InfoBlock * _infoBlock)
    {
        infoBlock = _infoBlock;
        size = infoBlock->size_x * infoBlock->size_y;
        subMatrix= new float[size];

        for (int i = 0; i < 4; ++i)
        {
            if (infoBlock->halo[i] && ((i == HALO_TYPE::UP)|| (i == HALO_TYPE::DOWN)))
            {
                halos[i] = new float [infoBlock->size_x];
            }

            else if (infoBlock->halo[i] && ((i == HALO_TYPE::LEFT)|| (i == HALO_TYPE::RIGHT)))
            {
                halos[i] = new float [infoBlock->size_y];
            }
            else if (!infoBlock->halo[i])
            {
                halos[i]= NULL;
            }
        }
    }

    bool get(int i,int j,  float & val)
    {
        if (i>=0 && i<infoBlock->size_y && j>=0 && infoBlock->size_x)
        {
            val= subMatrix [i*infoBlock->size_x + j];
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
            subMatrix[i] = val;
        }
    }
    bool set (int i, int j, float & val)
    {
        if (i>=0 && i<infoBlock->size_y && j>=0 && infoBlock->size_x)
        {
            subMatrix [i*infoBlock->size_x + j]= val;
            return true;
        }
        return false;
    }



    void block_receiving(int & MPI_root, MPI_Comm & MPI_COMM_CUBE, int tag)
    {
        MPI_Status status;
        MPI_Recv(subMatrix, infoBlock->size_x*infoBlock->size_y, MPI_FLOAT, MPI_root, tag, MPI_COMM_CUBE, &status);
    }

    friend std::ostream & operator <<( std::ostream &os, const Substate &substate )
    {

        os<<"Matrix di size: "<<substate.infoBlock->size_y<<" X "<<substate.infoBlock->size_x<<"\n";

        for (int i = 0; i < substate.infoBlock->size_y; ++i) {
            for (int j = 0; j < substate.infoBlock->size_x; ++j) {
                os<<substate.subMatrix[i*substate.infoBlock->size_x+ j]<<" ";


            }
            os<<"\n";
        }
        return os;
    }




};

#define NUMBER_OF_OUTFLOWS 4

class CellularAutomata
{
private:
    Substate z;
    Substate h;
    Substate f[NUMBER_OF_OUTFLOWS];

    double epsilon;
    double r;

    InfoBlock* infoBlock;

public:
    CellularAutomata (InfoBlock* infoBlock) : infoBlock(infoBlock),z(infoBlock),
        h(infoBlock)
    {
        f[0].setInfoBlock(infoBlock);
        f[1].setInfoBlock(infoBlock);
        f[2].setInfoBlock(infoBlock);
        f[3].setInfoBlock(infoBlock);
    }

    void transitionFunction ()
    {

    }


    void init (int & MPI_root, MPI_Comm & MPI_COMM_CUBE)
    {
        z.block_receiving(MPI_root, MPI_COMM_CUBE, 0);
        h.block_receiving(MPI_root, MPI_COMM_CUBE, 1);
    }

    friend std::ostream & operator <<( std::ostream &os, const CellularAutomata &CA )
    {

        os<<CA.h;
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

void block_sending(Reader* reader, InfoBlock * infoBlock, int num_procs,int size_x, MPI_Datatype & MPI_BLOCK_TYPE, int & MPI_root, MPI_Comm & MPI_COMM_CUBE, int tag);

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

    const char* path = argv[1];
    int totalSteps = atoi(argv[2]);
    int stepOffset = atoi(argv[3]);

    /* Get the number of processes created by MPI and their rank */
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int MPI_root = num_procs -1;

    /* Create a 2D cartesian topology of the processes */
    MPI_Dims_create(num_procs, NUMBER_OF_DIMENSIONS, cart_dimensions);

    MPI_Cart_create(MPI_COMM_WORLD, NUMBER_OF_DIMENSIONS, cart_dimensions, cart_periodicity, 0, &MPI_COMM_CUBE);



    /* Get relevant data from the created topology */
    MPI_Cart_coords(MPI_COMM_CUBE, rank, NUMBER_OF_DIMENSIONS, cart_coordinates);
    printf("rank %d : %d cciao %d \n", rank, cart_coordinates[0], cart_coordinates[1]);
    MPI_Cart_rank(MPI_COMM_CUBE, cart_coordinates, &rank);
    MPI_Cart_shift(MPI_COMM_CUBE, VERTICAL, 1, &neighbor_up, &neighbor_down);
    printf("rank %d : %d cciao %d \n", rank, cart_coordinates[0], cart_coordinates[1]);
    MPI_Cart_shift(MPI_COMM_CUBE, HORIZONTAL, 1, &neighbor_left, &neighbor_right);



    InfoBlock infoBlock;
    infoBlock.cart_coordinates= cart_coordinates;
    infoBlock.cart_dimensions = cart_dimensions;
    //    unsigned int size = ;

    Reader reader(path);
    if (rank == MPI_root)
    {
        reader.loadFromFile();

        std::cout<<reader<<std::endl;



    }
    MPI_Bcast(&reader.nCols, 1, MPI_INT,MPI_root, MPI_COMM_CUBE);
    MPI_Bcast(&reader.nRows, 1, MPI_INT,MPI_root, MPI_COMM_CUBE);

    if (rank != MPI_root)
    {
        //        MPI_Status status;
        //        MPI_Recv(&reader.getCols(), 1, MPI_INT, MPI_root, 0, MPI_COMM_CUBE, &status);
        //        MPI_Recv(&reader.getRows(), 1, MPI_INT, MPI_root, 0, MPI_COMM_CUBE, &status);
    }

    block_distribution(infoBlock, reader.nCols, reader.nRows);
    stampa(infoBlock, rank);
    CellularAutomata sciddica (&infoBlock);


        if (rank == MPI_root)
        {
            block_sending(&reader,&infoBlock,num_procs,reader.nCols,MPI_BLOCK_TYPE,MPI_root,MPI_COMM_CUBE,0);
            block_sending(&reader,&infoBlock,num_procs,reader.nCols,MPI_BLOCK_TYPE,MPI_root,MPI_COMM_CUBE,1);

        }
        if (rank != MPI_root)
        {
//TODO LA ROOT DEVE RIEMPIRE A MANO LA SUA

            //fare metodo che ricolleziona i dati nella matrice globale


            sciddica.init(MPI_root, MPI_COMM_CUBE);
//            substate.block_receiving(MPI_root,MPI_COMM_CUBE,1);

            cout<<"sono rank : "<< rank<<" \n"<<sciddica<<std::endl;


        }


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


void block_sending(Reader* reader, InfoBlock * infoBlock, int num_procs, int size_x, MPI_Datatype & MPI_BLOCK_TYPE, int & MPI_root, MPI_Comm & MPI_COMM_CUBE, int tag)
{


    //DEVE FARLO IL ROOT

    MPI_Type_vector(infoBlock->size_y, infoBlock->size_x, size_x, MPI_FLOAT, &MPI_BLOCK_TYPE);
    MPI_Type_commit(&MPI_BLOCK_TYPE);


    int starterIndex = 0;
    for (int dest =0; dest<num_procs; dest++ )
    {
        if(dest== MPI_root)
            continue;

        MPI_Send (&reader->getData()[starterIndex], 1, MPI_BLOCK_TYPE, dest, tag, MPI_COMM_CUBE);
        starterIndex+= infoBlock->size_y* infoBlock->size_x;
    }







}


