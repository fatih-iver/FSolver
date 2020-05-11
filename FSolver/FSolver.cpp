#include <mpi.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <thread>
#include <chrono>

using namespace std;

struct Offset {
    int I, J, K;
};

int calculateFlatIndex(int dim, int i, int j, int k) {
    return i * dim * dim + j * dim + k;
}

void fillWithZeros(double* cube, int dim) {
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                cube[calculateFlatIndex(dim, i, j, k)] = 0;
            }
        }
    }
}

void print(double* cube, int dim) {
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                cout << cube[calculateFlatIndex(dim, i, j, k)] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}

bool isWholeNumber(double num)
{
    return num == static_cast<int>(num);
}

void calculateMyOffset(int myRank, int cubeDimension, int subcubeDimension, Offset &myOffset) {

    myOffset.I = 0;
    myOffset.J = 0;
    myOffset.K = 0;

    for (int i = 0; i < myRank; i++) {

        myOffset.K += subcubeDimension;

        if (myOffset.K == cubeDimension) {

            myOffset.K = 0;

            myOffset.J += subcubeDimension;

            if (myOffset.J == cubeDimension) {

                myOffset.J = 0;

                myOffset.I += subcubeDimension;
            }

        }

    }



}

int main(int argc, char* argv[]) {

    int numberOfTasks, myRank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfTasks);

    double subcubePerDimension = cbrt(numberOfTasks);

    if (!isWholeNumber(subcubePerDimension)) {
        cout << "number of tasks is not cube root of a whole number!" << endl;
        return 0;
    }

    int cubeDimension = stoi(argv[1]);

    double subcubeDimension = cubeDimension / subcubePerDimension;

    if (!isWholeNumber(subcubeDimension)) {
        cout << "subcube dimension is not a whole number!" << endl;
        return 0;
    }

    const int pointsPerEdge = subcubeDimension + 1;

    double *subcube = new double[pointsPerEdge * pointsPerEdge * pointsPerEdge];

    fillWithZeros(subcube, pointsPerEdge);

    std::this_thread::sleep_for(std::chrono::milliseconds(myRank * 250));

    struct Offset myOffset;

    calculateMyOffset(myRank, cubeDimension, subcubeDimension, myOffset);

    cout << "My Rank: " << myRank << endl;

    cout << "My Offset: " << myOffset.I << " " << myOffset.J << " " << myOffset.K << endl;

    MPI_Finalize();

    return 0;
}