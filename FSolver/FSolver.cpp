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

int calculateFlatIndex(int pointsPerEdge, int i, int j, int k) {
    return i * pointsPerEdge * pointsPerEdge + j * pointsPerEdge + k;
}

int calculateFlatOffset(int pointsPerEdge, Offset myOffset) {
    return calculateFlatIndex(pointsPerEdge, myOffset.I, myOffset.J, myOffset.K);
}

void fillWithZeros(double* cube, int pointsPerEdge) {
    for (int i = 0; i < pointsPerEdge; i++) {
        for (int j = 0; j < pointsPerEdge; j++) {
            for (int k = 0; k < pointsPerEdge; k++) {
                cube[calculateFlatIndex(pointsPerEdge, i, j, k)] = 0;
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

    // MPI-Initialization -----------------------------------------------------------------------------

    int numberOfTasks, myRank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfTasks);

    // Pre-Checks -------------------------------------------------------------------------------------

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

    // Point - Calculations -------------------------------------------------------------

    const int pointsPerEdge = cubeDimension + 1;

    const int pointsPerSubEdge = subcubeDimension + 1;

    // Subcube - Initialization -------------------------------------------------------------

    double *subcube = new double[pointsPerSubEdge * pointsPerSubEdge * pointsPerSubEdge];

    fillWithZeros(subcube, pointsPerSubEdge);

    // Offsets -------------------------------------------------------------------------------------

    struct Offset myOffset;

    calculateMyOffset(myRank, cubeDimension, subcubeDimension, myOffset);

    Offset myNormalizedOffset;

    myNormalizedOffset.I = myOffset.I / subcubeDimension;
    myNormalizedOffset.J = myOffset.J / subcubeDimension;
    myNormalizedOffset.K = myOffset.K / subcubeDimension;

    int offsetFlatIndex = calculateFlatOffset(pointsPerEdge, myOffset);

    // Neighbours -----------------------------------------------------------------------------------

    int*** neighbours  = new int** [subcubePerDimension];

    for (int i = 0; i < subcubePerDimension; i++) {
        
        neighbours[i] = new int* [subcubePerDimension];

        for (int j = 0; j < subcubePerDimension; j++) {
           
            neighbours[i][j] = new int[subcubePerDimension];
       
        }

    }

    for (int i = 0; i < subcubePerDimension; i++) {
        for (int j = 0; j < subcubePerDimension; j++) {
            for (int k = 0; k < subcubePerDimension; k++) {
                neighbours[i][j][k] = calculateFlatIndex(subcubePerDimension, i, j, k);
            }
        }
    }

    /*

    if (myRank == 0) {
        for (int i = 0; i < subcubePerDimension; i++) {
            for (int j = 0; j < subcubePerDimension; j++) {
                for (int k = 0; k < subcubePerDimension; k++) {
                    cout << " " << i << " " << j << " " << k << " : " << neighbours[i][j][k] << endl;
                }
            }
        }
    }

    */

    int ME = myRank;

    int UP = -1;
    int DOWN = -1;
    
    int LEFT = -1;
    int RIGHT = -1;

    int FRONT = - 1;
    int BACK = -1;

    if (myNormalizedOffset.I - 1 >= 0) {
        UP = neighbours[myNormalizedOffset.I - 1][myNormalizedOffset.J][myNormalizedOffset.K];
    }

    if (myNormalizedOffset.I + 1 < subcubePerDimension) {
        DOWN = neighbours[myNormalizedOffset.I + 1][myNormalizedOffset.J][myNormalizedOffset.K];
    }

    if (myNormalizedOffset.J - 1 >= 0) {
        LEFT = neighbours[myNormalizedOffset.I][myNormalizedOffset.J - 1][myNormalizedOffset.K];
    }

    if (myNormalizedOffset.J + 1 < subcubePerDimension) {
        RIGHT = neighbours[myNormalizedOffset.I][myNormalizedOffset.J + 1][myNormalizedOffset.K];
    }

    if (myNormalizedOffset.K - 1 >= 0) {
        BACK = neighbours[myNormalizedOffset.I][myNormalizedOffset.J][myNormalizedOffset.K - 1];
    }

    if (myNormalizedOffset.K + 1 < subcubePerDimension) {
        FRONT = neighbours[myNormalizedOffset.I][myNormalizedOffset.J][myNormalizedOffset.K + 1];
    }
    
    // Sleep ---------------------------------------------------------------------------------

    std::this_thread::sleep_for(std::chrono::milliseconds(myRank * 250));

    if (myRank == 0) {
        
        cout << "Number of tasks: " << numberOfTasks << endl;
        cout << "Cube Dimension: " << cubeDimension << endl;
        cout << "Subcube Dimension: " << subcubeDimension << endl;

        cout << "Cube Volume: " << cubeDimension * cubeDimension * cubeDimension << endl;
        cout << "Subcube Volume: " << subcubeDimension * subcubeDimension * subcubeDimension << endl;
    
    }

    cout << endl << "My Rank: " << myRank << endl;

    cout << "My Offset: " << myOffset.I << " " << myOffset.J << " " << myOffset.K << endl;

    cout << "Flat : " << offsetFlatIndex << endl;

    cout << "Neighbours:" << " U " << UP << " D " << DOWN << " L " << LEFT << " R " << RIGHT << " F " << FRONT << " B " << BACK << endl;
   
    // Finalize -----------------------------------------------------------------------------------

    MPI_Finalize();

    return 0;
}