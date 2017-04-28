#include <iostream>
#include <unistd.h>
#include "../Time.h"
#include <functional>
using namespace std;

#define getToroidalX(index, size) (   (index)<0?((size)+(index)):( (index)>((size)-1)?((index)-(size)):(index) )   )
const int nTypes = 4;
string types[nTypes] = {"seriale", "parallel for", "static", "dynamic"};
auto checkResult = []()->bool {
    return true;
};

class GameOfLife {
public:

    GameOfLife(int nRows, int nCol, int nThreadsMin, int nThreadsMax) : nRows(nRows), nCol(nCol), nThreadsMin(nThreadsMin), nThreadsMax(nThreadsMax), time(types, nThreadsMax - nThreadsMin + 1) {
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
                sum += matrix[0][getToroidalX(i + k, nRows)][getToroidalX(j + z, nCol)];
            }
        return sum;
    }

    ~GameOfLife() {
        for (int i = 0; i < nRows; i++) {
            delete matrix[0][i];
            delete matrix[1][i];
        }
        delete [] matrix[0];
        delete [] matrix[1];
    }

    void run(int steps, int printSteps) {

        serialRun(steps, printSteps);
        for (int i = nThreadsMin; i <= nThreadsMax; i++) {
            parallelForRun(steps, i, printSteps);
            parallelForStaticRun(steps, i, printSteps);
            parallelForDynamicRun(steps, i, printSteps);
        }
        time.print();
    }


private:

    bool ** matrix [2];

    int nRows;
    int nCol;


    int nThreadsMin;
    int nThreadsMax;
    Time<nTypes> time;

    void swap() {
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCol; j++) {
                matrix[0][i][j] = matrix[1][i][j];
            }
        }
    }

    void parallelForSwap() {
#pragma omp for
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCol; j++) {
                matrix[0][i][j] = matrix[1][i][j];
            }
        }
    }

    void parallelForStaticSwap() {
#pragma omp for schedule(static, 10)
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCol; j++) {
                matrix[0][i][j] = matrix[1][i][j];
            }
        }
    }

    void parallelForDynamicSwap() {
#pragma omp for schedule(dynamic, 10)
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCol; j++) {
                matrix[0][i][j] = matrix[1][i][j];
            }
        }
    }

    void transitionFunction(int i, int j) {
        int sum = sumNeighborhood(i, j);
        if (sum == 3 || (sum == 2 && get(i, j))) {
            set(i, j, 1);

        } else
            set(i, j, 0);
    }

    void serialRun(int steps, int printStep) {
        time.start();
        for (int s = 0; s < steps; s++) {
            for (int i = 0; i < nRows; i++) {
                for (int j = 0; j < nCol; j++) {

                    transitionFunction(i, j);
                }


            }
            swap();
            time.wait();
            if ((s + 1) % printStep == 0)
                print();
            time.continue_();
        }
        time.stop(0, 0, std::bind(checkResult));
    }

    void parallelForRun(int steps, int nThreads, int printStep) {
        time.start();
        for (int s = 0; s < steps; s++) {
#pragma omp parallel 
            {
                omp_set_num_threads(nThreads);
#pragma omp for
                for (int i = 0; i < nRows; i++) {
                    for (int j = 0; j < nCol; j++) {

                        transitionFunction(i, j);
                    }


                }
                parallelForSwap();
            }
            time.wait();
            if ((s + 1) % printStep == 0)
                print();
            time.continue_();
        }
        time.stop(1, nThreads - nThreadsMin, std::bind(checkResult));
    }

    void parallelForStaticRun(int steps, int nThreads, int printStep) {
        time.start();
        for (int s = 0; s < steps; s++) {
#pragma omp parallel 
            {
                omp_set_num_threads(nThreads);
#pragma omp for schedule(static, 10)
                for (int i = 0; i < nRows; i++) {
                    for (int j = 0; j < nCol; j++) {

                        transitionFunction(i, j);
                    }


                }
                parallelForStaticSwap();
            }
            time.wait();
            if ((s + 1) % printStep == 0)
                print();
            time.continue_();
        }
        time.stop(2, nThreads - nThreadsMin, std::bind(checkResult));
    }

    void parallelForDynamicRun(int steps, int nThreads, int printStep) {
        time.start();
        for (int s = 0; s < steps; s++) {
#pragma omp parallel 
            {
                omp_set_num_threads(nThreads);
#pragma omp for schedule(dynamic, 10)
                for (int i = 0; i < nRows; i++) {
                    for (int j = 0; j < nCol; j++) {

                        transitionFunction(i, j);
                    }


                }
                parallelForDynamicSwap();
            }
            time.wait();
            if ((s + 1) % printStep == 0)
                print();
            time.continue_();
        }
        time.stop(3, nThreads - nThreadsMin, std::bind(checkResult));
    }

};
