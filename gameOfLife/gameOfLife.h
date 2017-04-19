#include <iostream>
#include <unistd.h>
using namespace std;

#define getToroidalX(index, size) (   (index)<0?((size)+(index)):( (index)>((size)-1)?((index)-(size)):(index) )   )

class GameOfLife {
private:
    bool** next;
    bool** curr;
    int nRows;
    int nCol;

public:

    GameOfLife(int nRows, int nCol) : nRows(nRows), nCol(nCol) {
        next = new bool*[nRows];
        curr = new bool*[nRows];
        for (int i = 0; i < nRows; i++) {
            curr[i] = new bool[nCol];
            next[i] = new bool[nCol];
            for (int j = 0; j < nCol; j++) {
                curr[i][j] = 0;
                next[i][j] = 0;

            }

        }
    }

    void print() {
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCol; j++)
                cout << curr[i][j] << "  ";
            cout << endl;
        }
    }

    void init(int i, int j, bool value) {
        assert(i < nRows && j < nCol);
        curr[i][j] = value;
    }

    void set(int i, int j, bool value) {
        assert(i < nRows && j < nCol);
        next[i][j] = value;
    }

    bool get(int i, int j) {
        assert(i < nRows && j < nCol);
        return curr[i][j];
    }

    int sumNeighborhood(int i, int j) {
        int sum = 0;
        for (int k = -1; k <= 1; k++)
            for (int z = -1; z <= 1; z++) {
                if (k == 0 && z == 0)
                    continue;
                sum += curr[getToroidalX(i+k,nRows)][getToroidalX(j+z,nCol)];
            }
        return sum;
    }

    void swap() {
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCol; j++) {
                curr[i][j] = next[i][j];
            }
        }
    }

    ~GameOfLife() {
        for (int i = 0; i < nRows; i++) {
            delete next[i];
            delete curr[i];
        }
        delete [] next;
        delete [] curr;
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