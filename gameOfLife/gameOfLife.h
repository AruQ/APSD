#include <iostream>
#include <unistd.h>
using namespace std;

#define getToroidalX(index, size) (   (index)<0?((size)+(index)):( (index)>((size)-1)?((index)-(size)):(index) )   )

class GameOfLife {
private:
    bool ** matrix [2];

    int nRows;
    int nCol;

public:

    GameOfLife(int nRows, int nCol) : nRows(nRows), nCol(nCol) {
        matrix[0] = new bool*[nRows];
        matrix[1] = new bool*[nRows];

        for (int i = 0; i < nRows; i++) {
            matrix[0][i] = new bool[nCol];
            matrix[1][i] = new bool[nCol];
            for (int j = 0; j < nCol; j++) {
                matrix[0][i][j] = 0;
                matrix[1][i][j] = 0;

            }

        }
    }

    void print() {
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCol; j++)
                cout << matrix[0][i][j] << "  ";
            cout << endl;
        }
    }

    void init(int i, int j, bool value) {
        assert(i < nRows && j < nCol);
        matrix[0][i][j] = value;
    }

    void set(int i, int j, bool value) {
        assert(i < nRows && j < nCol);
        matrix[1][i][j] = value;
    }

    bool get(int i, int j) {
        assert(i < nRows && j < nCol);
        return matrix[0][i][j];
    }

    int sumNeighborhood(int i, int j) {
        int sum = 0;
        for (int k = -1; k <= 1; k++)
            for (int z = -1; z <= 1; z++) {
                if (k == 0 && z == 0)
                    continue;
                sum += matrix[0][getToroidalX(i+k,nRows)][getToroidalX(j+z,nCol)];
            }
        return sum;
    }

    void swap() {
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCol; j++) {
                matrix[0][i][j] = matrix[1][i][j];
            }
        }
    }

    ~GameOfLife() {
        for (int i = 0; i < nRows; i++) {
            delete matrix[0][i];
            delete matrix[1][i];
        }
        delete [] matrix[0];
        delete [] matrix[1];
    }

    void run(int steps) {
        for (int s = 0; s < steps; s++) {
#pragma omp for
            for (int i = 0; i < nRows; i++) {
                for (int j = 0; j < nCol; j++) {

                    transitionFunction(i, j);
                }


            }
            swap();
            print();
            //usleep(30000);
            //cout<<"-------------------------"<<endl;
        }


    }
private:
    
    

    void transitionFunction(int i, int j) {
        int sum = sumNeighborhood(i, j);
        if (sum == 3 || (sum == 2 && get(i, j))) {
            set(i, j, 1);

        } else
            set(i, j, 0);
    }

};
