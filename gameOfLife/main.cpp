#include <iostream>
#include <cstdlib>
#include <cassert>
#include "gameOfLife.h"
using namespace std;


void run();

int main(int argc, char ** argv) {


    assert(argc >= 7);
    int nRow = atoi(argv[1]);
    int nCol = atoi(argv[2]);
    int steps = atoi(argv[3]);
	int nThreadMin = atoi(argv[4]);
	int nThreadMax = atoi(argv[5]);
	int printSteps = atoi(argv[6]);
    GameOfLife g(nRow, nCol,nThreadMin,nThreadMax);
    g.init(0, 2, 1);
    g.init(1, 0, 1);
    g.init(1, 2, 1);
    g.init(2, 1, 1);
    g.init(2, 2, 1);

    g.run(steps,printSteps);

    //    g.print();
    return 0;
}



