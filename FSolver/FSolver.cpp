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

void print(double* cube, int pointsPerEdge) {
    for (int i = 0; i < pointsPerEdge; i++) {
        for (int j = 0; j < pointsPerEdge; j++) {
            for (int k = 0; k < pointsPerEdge; k++) {
                cout << cube[calculateFlatIndex(pointsPerEdge, i, j, k)] << " ";
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

double u(double x, double y, double z) {
    return x * y * z;
}

double f(double x, double y, double z) {
    return 1 + 1 + 1;
}

int main(int argc, char* argv[]) {

    // MPI-Initialization -----------------------------------------------------------------------------

    int numberOfTasks, myRank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfTasks);

    // Points per Edge Calculations -------------------------------------------------------------------------------------

    double subcubePerEdge = cbrt(numberOfTasks);

    if (!isWholeNumber(subcubePerEdge)) {
        cout << "cube root of number of tasks is not a whole number!" << endl;
        return 0;
    }

    int n = stoi(argv[1]);

    int pointsPerEdge = n + 1;

    double _pointsPerSubEdge = pointsPerEdge / subcubePerEdge;

    if (!isWholeNumber(_pointsPerSubEdge)) {
        cout << "points per subcube is not a whole number!" << endl;
        return 0;
    }

    int pointsPerSubEdge = (int) _pointsPerSubEdge;

    int pointsPerGhostEdge = 1 + pointsPerSubEdge + 1;

    // Subcube - Initialization -------------------------------------------------------------

    double *subcube = new double[pointsPerSubEdge * pointsPerSubEdge * pointsPerSubEdge];

    fillWithZeros(subcube, pointsPerSubEdge);

    double* ghostCube = new double[pointsPerGhostEdge * pointsPerGhostEdge * pointsPerGhostEdge];

    fillWithZeros(ghostCube, pointsPerGhostEdge);

    

    // Offsets -------------------------------------------------------------------------------------

    struct Offset myOffset;

    calculateMyOffset(myRank, pointsPerEdge, pointsPerSubEdge, myOffset);

    int offsetFlatIndex = calculateFlatOffset(pointsPerEdge, myOffset);

    Offset myNormalizedOffset = {myOffset.I / pointsPerSubEdge, myOffset.J / pointsPerSubEdge, myOffset.K / pointsPerSubEdge };

    // Neighbours -----------------------------------------------------------------------------------

    int*** neighbours  = new int** [subcubePerEdge];

    for (int i = 0; i < subcubePerEdge; i++) {
        
        neighbours[i] = new int* [subcubePerEdge];

        for (int j = 0; j < subcubePerEdge; j++) {
           
            neighbours[i][j] = new int[subcubePerEdge];
       
        }

    }

    for (int i = 0; i < subcubePerEdge; i++) {
        for (int j = 0; j < subcubePerEdge; j++) {
            for (int k = 0; k < subcubePerEdge; k++) {
                neighbours[i][j][k] = calculateFlatIndex(subcubePerEdge, i, j, k);
            }
        }
    }

    // Find Neighbours ---------------------------------------------------------------------------------------

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

    if (myNormalizedOffset.I + 1 < subcubePerEdge) {
        DOWN = neighbours[myNormalizedOffset.I + 1][myNormalizedOffset.J][myNormalizedOffset.K];
    }

    if (myNormalizedOffset.J - 1 >= 0) {
        LEFT = neighbours[myNormalizedOffset.I][myNormalizedOffset.J - 1][myNormalizedOffset.K];
    }

    if (myNormalizedOffset.J + 1 < subcubePerEdge) {
        RIGHT = neighbours[myNormalizedOffset.I][myNormalizedOffset.J + 1][myNormalizedOffset.K];
    }

    if (myNormalizedOffset.K - 1 >= 0) {
        BACK = neighbours[myNormalizedOffset.I][myNormalizedOffset.J][myNormalizedOffset.K - 1];
    }

    if (myNormalizedOffset.K + 1 < subcubePerEdge) {
        FRONT = neighbours[myNormalizedOffset.I][myNormalizedOffset.J][myNormalizedOffset.K + 1];
    }
    
    // DEBUG INFO -------------------------------------------------------------------------------------------

    std::this_thread::sleep_for(std::chrono::milliseconds(myRank * 250));

    if (myRank == 0) {
    
        cout << "-------- DEBUG INFO --------" << endl;

        cout << endl;
        cout << "Number of Tasks: " << numberOfTasks << endl;
        cout << "Subcube Per Edge: " << subcubePerEdge << endl;
        cout << "n: " << n << endl;
        cout << "Points per Edge: " << pointsPerEdge << endl;
        cout << "Points per Sub Edge: " << pointsPerSubEdge << endl;
        cout << "Points per Ghost Edge: " << pointsPerGhostEdge << endl;
        cout << endl;
    
    }

    cout << "My Rank: " << myRank << endl;

    cout << "My Offset: " << myOffset.I << " " << myOffset.J << " " << myOffset.K << endl;

    cout << "My Flat Offset : " << offsetFlatIndex << endl;

    cout << "My Normalized Offset: " << myNormalizedOffset.I << " " << myNormalizedOffset.J << " " << myNormalizedOffset.K << endl;

    cout << "Neighbours:" << " U " << UP << " D " << DOWN << " L " << LEFT << " R " << RIGHT << " F " << FRONT << " B " << BACK << endl;

    cout << endl;

    print(subcube, pointsPerSubEdge);

    //print(ghostCube, pointsPerGhostEdge);

    // Finalize -------------------------------------------------------------------------------------------

    MPI_Finalize();

    return 0;
}